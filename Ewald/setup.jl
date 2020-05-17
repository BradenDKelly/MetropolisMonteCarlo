#module Setup
include("structs.jl")
function FindMolType(type::String,moleculeList::Array)

    for i in eachindex(moleculeList)
        item = moleculeList[i]
        if lowercase(item.resnm[1]) == lowercase(type)
            return i
        end
    end
end
function FindAtomInMol(atomType::String,molType::Int64)
    for (i,item) in enumerate(systemTop.molParams[molType].atoms)
        if lowercase(item.atomnm) == lowercase(atomType)
            return i
        end
    end
end

function FindNumericAtomType(molType::Int64,atomType::Int64,systemTop::FFParameters)
    name = systemTop.molParams[molType].atoms[atomType].type
    for (i,item) in enumerate(systemTop.atomTypes)
        if lowercase(item.name) == lowercase(name)
            return i
        end
    end
end
export ReadPDB, ReadTopFile

function ReadPDB(pdbName)

    molName = split(pdbName,".")[1]
    box = Vector{Float64}(undef,3)
    coords = []
    atomnm = []
    resnm = []
    resnr = []
    elem = []
    df = 1.0  # factor from converting between Angstrom and nm
    #conect = [] # currently disabled
    open(pdbName) do file
        for line in eachline(file)
            if occursin("ATOM",line) || occursin("HETATM",line)
                push!(coords,
                    [
                    df*parse(Float64,line[31:38]),
                    df*parse(Float64,line[40:46]),
                    df*parse(Float64,line[48:55])
                    ])
                push!(atomnm,strip(line[12:15]))
                push!(resnm,strip(line[17:21]))
                push!(resnr,parse(Int64,(line[22:27])) )
                if length(line) >= 77
                    push!(elem,strip(line[77:78]))
                else
                    push!(elem,"")
                end
            elseif occursin("CRYST1",line)
                box = df*[parse(Float64,split(line)[2]),parse(Float64,split(line)[3]),
                    parse(Float64,split(line)[4]) ]
            #elseif occursin("CONECT",line)
            #    push!(conect,[parse(Int64,line[7:12]), parse(Int64,line[13:18]) ])
            end # if ATOM/HETATM/CRYST/CONECT
        end # for line in eachline(file)
    end # close open do file

    allsame(x) = all(y -> y == first(x), x) # checks if all elements in an array are the same
    if allsame(resnm)
        println(resnm[1]," is alone in ", pdbName)
    else
        println(resnm[1]," is not alone in ", pdbName)
    end

    # Create object to store all of this info (i.e., it is called Topology)
    # Initialize Topology object
    #x = length(elem)
    topology=Topology(
        molName,
        [box[i] for i in 1:3],
        [SVector(coord...) for coord in coords],
        [atomnm[i] for i in eachindex(atomnm)],
        [resnm[i] for i in eachindex(resnm)],
        [resnr[i] for i in eachindex(resnr)],
        [elem[i] for i in eachindex(elem)]
        )
    return topology
end # function ReadPDB

function ReadTopFile(top_file::AbstractString)
#=
This function is passed a topology file name. The topology is scanned and all
topology information is saved into a "topology" object. This object contains
all sections of a gromacs topology file. For every molecule in the topology file
a molecule object is created to store information specific to a molecule (MolParam).
All of these molecule objects are stored in a Vector inside the topology object.
First pass through all objects are stored in temporary lists since lists are
easily appended to (push!), then transferred to permanent object structures once
vector lengths are known.
This will read either a flat topology file (no #includes) or a topology file
which stores specific molecule information in "#include" a.k.a. itp files.
It will not read #includes for information about FF's.
=#
    # Read forcefield and topology file
    #ffParameters = FFParameters()
    local system, defaults # these are temporary objects (structs) for this function
    atomtypes = []
    moleculetypes = []                  # 1 object for each molecule type
    #moleculetypesExcl = []
    #moleculetypesName = []
    topName::String = "?"
    topnrexcl::Int64 = 0
    molparams = MolParam[]
    topatoms = Atoms[]
    topbonds = Bonds[]
    toppairs = Pairs[]
    topangles = Angles[]
    topdihedrals = Dihedrals[]
    zone = "?"
    #ffParameters.molecules = Dict{String,Int64}()
    molecules=Dict{String,Int64}()
    #molName2Num=Dict{String,Int64}()
    molNumber = 0
    molNumberCheck = 0
    molType = 0 # dummy counter
    # First, scan to the bottom and get the number of molecules in the system
    # double check is implemented, we count number of molecules defined along the way
    # Actually, we do the reverse of these two things :)
    open(top_file) do file
        for lin in eachline(file)
            line = split(strip(lin),";")[1]
            if occursin("[ moleculetype ]",line) molNumber += 1 end
            if occursin("#include",line)
                filename = split(line,'"')[2]
                open(filename) do file
                    for lins in eachline(file)
                        lines = split(strip(lins),";")[1]
                        if occursin("[ moleculetype ]",lines) molNumber += 1 end
                    end
                end
            end
            #=Implement a double check on the number of molecule types in
            topology file by checking the [ molecule ] section =#
            if occursin("[ molecules ]",line) zone = strip(line[2:end-1]) end
            if zone == "molecules" && length(split(line)) == 2
                molNumberCheck +=1
            end
        end
    end
    if molNumber != molNumberCheck
        error("The number of moleculetypes in topology file ",top_file,
        " does not match the number in section [ molecules ] of said file.")
     end

################################################################################
#                  Read parameters from Topology File
################################################################################

    open(top_file) do file

        for lin in eachline(file)
            line = strip(split(strip(lin),";")[1])   # I only care about what is on the LHS
            if length(line) == 0              # of the first semi-colon
                continue
            end
            if startswith(line, '[') && endswith(line, ']')
                if ( (zone == "moleculetype") || (zone == "system") ) && (length(topatoms) > 0)
                    push!(molparams,MolParam(topName,topnrexcl,
                        [topatoms[i] for i in 1:length(topatoms)],
                        [topbonds[i] for i= 1:length(topbonds)],
                        [topangles[i] for i in 1:length(topangles)],
                        [topdihedrals[i] for i in 1:length(topdihedrals)]
                        ))
                    topatoms = Atoms[]
                    topbonds = Bonds[]
                    toppairs = Pairs[]
                    topangles = Angles[]
                    topdihedrals = Dihedrals[]
                    molType += 1
                elseif ( (zone == "moleculetype") || (zone == "system") ) && length(topatoms) == 0
                    molType += 1
                end
                zone = strip(line[2:end-1])
                continue
            end
            if zone == "defaults" && length(split(line) ) > 2 #&& defaults == false
                s = split(line)

                #=; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
                      1             2               yes             0.5       0.8333 =#
                defaults = Defaults( parse(Int64,s[1]),
                    parse(Int64,s[2]),s[3],parse(Float64,s[4]),parse(Float64,s[5]) )
            elseif zone == "atomtypes" && length(split(line) ) > 2
                inline = split(line)
            #= ; name   atomicnr    mass      charge     ptype  sigma(nm)     epsilon(kJ/mol)
                       O1        O1    15.99940   -0.8340      A   0.315061     0.6364000 =#
                push!( atomtypes, Atomtypes(inline[1],inline[2],
                    parse(Float64,inline[3]), parse(Float64,inline[4]),inline[5],
                    parse(Float64,inline[6]),parse(Float64,inline[7]) ) )
            elseif zone == "molecules" # last item in topology
                inline = split(line)
                molecules[inline[1]] = parse(Int64,inline[2])
            ####################################################################
            #
            #        Get the tricker stuff from either itp or top file
            #
            ####################################################################
            elseif occursin("#include",line)
                zone = "itp"
                # open itp file and get atoms/bonds/angles/pairs/dihedrals
                filename = split(line,'"')[2] # name is in quotes, which is the 2nd segment ...     #include "name"
                itpzone = "?"                 # using double quotation mark as delimiter -> segment:   (1)     (2)
                println("Reading Topology information from itp file for ", filename)
                open(filename) do file
                    # First store everything in lists, then add to appropriate
                    # objects as vectors once lengths are known
                    itpName::String = "?"
                    itpnrexcl::Int64 = 0
                    itpatoms = Atoms[]
                    itpbonds = Bonds[]
                    itppairs = Pairs[]
                    itpangles =Angles[]
                    itpdihedrals = Dihedrals[]
                    for itplin in eachline(file)                # read each line of file

                        itpline = strip(split(itplin,";")[1])   # I only care about what is on the LHS
                        if length(itpline) == 0                 # of the first semi-colon
                            continue
                        end
                        if startswith(itpline, '[') && endswith(itpline, ']')
                            itpzone = strip(itpline[2:end-1])  # there is a word between the brackets... get that word...
                            continue
                        end
                        if itpzone == "moleculetype"  && length(split(itpzone)) >= 1#Make new molecule object
                            inline = split( itpline )
                            itpName = inline[1]
                            itpnrexcl = parse( Int64,inline[2] )

                        elseif itpzone == "atoms" && length(split(itpline) ) > 2
                #= ;   nr type  resnr residue  atom   cgnr     charge    mass
                        1  O1     1    SOL      O1      1      -0.8340  15.99940
                =#
                            s = split( itpline )
                            push!(itpatoms,Atoms( parse(Int64,s[1]),s[2],
                                parse(Int64,s[3]),s[4],s[5],parse(Int64,s[6]),
                                parse(Float64,s[7]), parse(Float64,s[8]) ) )

                        elseif itpzone == "bonds" && length(split(itpline) ) > 2
                            #=;   ai     aj funct   r             k
                                    1      2   1    1.0860e-01  2.8937e+05  =#
                                s = split( itpline )
                                push!(itpbonds,Bonds( parse(Int64,s[1]),
                                    parse(Int64,s[2]),parse(Int64,s[3]),
                                    parse(Float64,s[4]), parse(Float64,s[5]) ) )

                        elseif itpzone == "pairs" && length(split(itpline) ) > 2
                            #=;   ai     aj    funct
                                    1      6      1  =#
                                s = split( itpline )
                                push!(itppairs,Pairs( parse(Int64,s[1]),
                                    parse(Int64,s[2]), parse(Int64,s[3]) ) )

                        elseif itpzone == "angles" && length(split(itpline) ) > 2
                            #=
                            ai  aj ak  funct   theta         cth
                            1    3  4     1    1.1988e+02  4.0317e+02
                            =#
                            s = split( itpline )
                            push!(itpangles,Angles( parse(Int64,s[1]),
                                parse(Int64,s[2]),parse(Int64,s[3]),
                                parse(Int64,s[4]),deg2rad(parse(Float64,s[5])),
                                parse(Float64,s[6]) ) )

                        elseif itpzone == "dihedrals" && length(split(itpline) ) > 2

                            s = split( itpline )
                            if length(s) > 9   # improper dihedral
                            #=
                            i j k l func   C0   C1     C2     C3   C4    C5
                            1 3 5 6  3  30.334 0.00 -30.334  0.00 0.00 0.000
                            =#
                                push!(itpdihedrals,Dihedrals( parse(Int64,s[1]),
                                    parse(Int64,s[2]),parse(Int64,s[3]),
                                    parse(Int64,s[4]),parse(Int64,s[5]),
                                    [ parse(Float64,s[6]),parse(Float64,s[7]),
                                    parse(Float64,s[8]),parse(Float64,s[9]),
                                    parse(Float64,s[10]),parse(Float64,s[11]) ] ) )
                            else # proper dihedral
                            #=
                            i  j  k  l  func phase     kd   pn
                            1  5  3  4   1   180.0  4.6024   2
                            =#
                                push!(itpdihedrals,Dihedrals( parse(Int64,s[1]),
                                    parse(Int64,s[2]),parse(Int64,s[3]),
                                    parse(Int64,s[4]),parse(Int64,s[5]),
                                    parse(Float64,s[6]),parse(Float64,s[7]),
                                    parse(Int64,s[8]) ) )
                            end  # end proper vs imprper dihedral check
                        end # ends molecule properties checks

                    end # ends looping through each line of file
                    push!(molparams,MolParam(itpName,itpnrexcl,
                        [itpatoms[i] for i in 1:length(itpatoms)],
                        [itpbonds[i] for i in 1:length(itpbonds)],
                        [itpangles[i] for i in 1:length(itpangles)],
                        [itpdihedrals[i] for i in 1:length(itpdihedrals)]
                        ))
                end # closes itp file

            elseif zone == "moleculetype"

                inline = split( line )
                topName = inline[1]       #instantiate a molecule object and at to the type list
                topnrexcl = parse( Int64,inline[2] )

            elseif zone == "atoms" && length(split(line) ) > 3
            #= ;   nr type  resnr residue  atom   cgnr     charge    mass
            1  O1     1    SOL      O1      1      -0.8340  15.99940
            =#
                s = split( line )
                push!(topatoms,Atoms( parse(Int64,s[1]),s[2],
                    parse(Int64,s[3]),s[4],s[5],parse(Int64,s[6]),
                    parse(Float64,s[7]), parse(Float64,s[8]) ) )

            elseif zone == "bonds" && length(split(line) ) > 2
                #=;   ai     aj funct   r             k
                        1      2   1    1.0860e-01  2.8937e+05  =#
                    s = split( line )

                    push!(topbonds,Bonds( parse(Int64,s[1]),
                        parse(Int64,s[2]),parse(Int64,s[3]),
                        parse(Float64,s[4]), parse(Float64,s[5]) ) )

            elseif zone == "pairs" && length(split(line) ) > 2
                #=;   ai     aj    funct
                        1      6      1  =#
                    s = split( line )
                    push!(toppairs,Pairs( parse(Int64,s[1]),
                        parse(Int64,s[2]), parse(Int64,s[3]) ) )

            elseif zone == "angles" && length(split(line) ) > 2
                #=
                ai  aj ak  funct   theta         cth
                1    3  4     1    1.1988e+02  4.0317e+02
                =#
                s = split( line )
                push!(topangles,Angles( parse(Int64,s[1]),
                    parse(Int64,s[2]),parse(Int64,s[3]),
                    parse(Int64,s[4]),deg2rad(parse(Float64,s[5])),
                    parse(Float64,s[6]) ) )

            elseif zone == "dihedrals" && length(split(line) ) > 2

                s = split( line )
                if length(s) > 9   # improper dihedral
                #=
                i j k l func   C0   C1     C2     C3   C4    C5
                1 3 5 6  3  30.334 0.00 -30.334  0.00 0.00 0.000
                =#
                    push!(topdihedrals,Dihedrals( parse(Int64,s[1]),
                        parse(Int64,s[2]),parse(Int64,s[3]),
                        parse(Int64,s[4]),parse(Int64,s[5]),
                        [ parse(Float64,s[6]),parse(Float64,s[7]),
                        parse(Float64,s[8]),parse(Float64,s[9]),
                        parse(Float64,s[10]),parse(Float64,s[11]) ] ) )
                else # proper dihedral
                #=
                i  j  k  l  func phase     kd   pn
                1  5  3  4   1   180.0  4.6024   2
                =#
                    push!(topdihedrals,Dihedrals( parse(Int64,s[1]),
                        parse(Int64,s[2]),parse(Int64,s[3]),
                        parse(Int64,s[4]),parse(Int64,s[5]),
                        parse(Float64,s[6]),parse(Float64,s[7]),
                        parse(Int64,s[8]) ) )
                end  # end proper vs imprper dihedral check
            elseif zone == "system" # second last item in topology
                system = line
                println("System name is: ", system)
            end # ends molecule properties checks
                    # search in this file for atoms/bonds/angles/pairs/dihedrals
        end
    end # closes top_file

    systemTopology = FFParameters(defaults,
                       [atomtypes[i] for i in 1:length(atomtypes)],
                       [molparams[i] for i in 1:length(molparams)],
                       system,
                       molecules)
    return systemTopology
end #ReadTopFile

" The molecular Dynamics Version"
function MakeAtomArrays(systemTop::FFParameters,atomsPDB::Topology)

    atomArray = Vector{PerAtomStruct}(undef, length(atomsPDB.resnr))
    pArray = Vector{ParticleAtom{Float64,Int64}}(undef, length(atomsPDB.resnr))
    # atomArray
    if atomsPDB.resnr[1] != 1                # just a double check on molecule numbering
        diff = Int(1 - atomsPDB.resnr[1])
        for i in eachindex(atomsPDB.resnr)
            atomsPDB.resnr[i] += diff
        end
    end
    for i in 1:length(atomsPDB.r)
        #assign coords,mass,atomtype,molecule atom belongs to, neighborlist
        mol = FindMolType(String(atomsPDB.resnm[i]),moleculeList )
        atom = FindAtomInMol( String(atomsPDB.atomnm[i]),mol )
        atomTypeNumber = FindNumericAtomType(mol,atom,systemTop)   # used for tables
        mass = systemTop.molParams[mol].atoms[atom].mass
        charge = systemTop.molParams[mol].atoms[atom].charge
        atomArray[i] = PerAtomStruct(atomsPDB.atomnm[i],
            atomTypeNumber,
            atomsPDB.resnm[i],
            atomsPDB.resnr[i],
            atomsPDB.r[i],
            XYZ( Velocity(mass,temperature)... ),   # initialize velocity
            XYZ(0.0, 0.0, 0.0),                # initialize force
            mass,
            charge)
    end
    # pArray
    for i in 1:length(atomsPDB.r)
        mol = FindMolType(String(atomsPDB.resnm[i]),moleculeList ) # returns molecule #
        atom = FindAtomInMol( String(atomsPDB.atomnm[i]),mol )     # returns atom # in that molecule
        atomTypeNumber = FindNumericAtomType(mol,atom,systemTop)   # returns atom # in all atoms list # used for tables
        mass = systemTop.molParams[mol].atoms[atom].mass
        charge = systemTop.molParams[mol].atoms[atom].charge
        pArray[i] = ParticleAtom(SVector{3}([atomsPDB.r[i]...] ),
                                 SVector{3}(Velocity(mass,298.15)...),
                                 zeros(SVector{3}),
                                 atomTypeNumber,
                                 atomsPDB.resnr[i],
                                 mass,
                                 charge
                                 )
    end

    push!(warnings,"Temperature is hardcoded to 298.15 for maxwell-boltzmann velocity distribution on line 704")
    for i in eachindex(atomArray)
        if atomArray[i].atomType > 13 println(atomArray[i].atomType) end
    end

    return atomArray, pArray
end

" The KMC Version"
function MakeAtomArrays(systemTop::FFParameters,atomsPDB::Topology, mode::String)

    pArray = Vector{ParticleAtomKMC{Float64,Int64}}(undef, length(atomsPDB.resnr))
    # pArray
    if atomsPDB.resnr[1] != 1
        diff = Int(1 - atomsPDB.resnr[1])
        for i in eachindex(atomsPDB.resnr)
            atomsPDB.resnr[i] += diff
        end
    end

    # Check that all atoms are in a positive box. Usually the box is centered at 0,
    # so half the molecules have negative coordinates... I use a box that starts at 0
    xshift = 0.0
    yshift = 0.0
    zshift = 0.0

    xshift = minimum([atomsPDB.r[i][1] for i=1:length(atomsPDB.r)])
    yshift = minimum([atomsPDB.r[i][2] for i=1:length(atomsPDB.r)])
    zshift = minimum([atomsPDB.r[i][3] for i=1:length(atomsPDB.r)])

    if xshift < 0.0 || yshift < 0.0 || zshift < 0.0
        println("X,Y,Z shifts to coords are: ", xshift, " ", yshift, " ", zshift)

        for i in 1:length(atomsPDB.r)
            #atomsPDB.r[i] = XYZ(atomsPDB.r[i].x + abs(xshift),
            #                atomsPDB.r[i].y + abs(yshift), atomsPDB.r[i].z + abs(zshift) )
            atomsPDB.r[i] =SVector(atomsPDB.r[i][1] + abs(xshift),
                                   atomsPDB.r[i][2] + abs(yshift),
                                   atomsPDB.r[i][3] + abs(zshift))
        end
    end

    # right, not carry on.
    for i in 1:length(atomsPDB.r)
        mol = FindMolType(String(atomsPDB.resnm[i]),moleculeList )
        atom = FindAtomInMol( String(atomsPDB.atomnm[i]),mol )
        atomTypeNumber = FindNumericAtomType(mol,atom,systemTop)   # used for tables
        charge = systemTop.molParams[mol].atoms[atom].charge
        mass = systemTop.molParams[mol].atoms[atom].mass

        pArray[i] = ParticleAtomKMC( atomsPDB.resnr[i],
                                mol,
                                atomTypeNumber,
                                mass,
                                SVector{3}([atomsPDB.r[i]...] ),
                                charge
                                )
    end
    ##############################
    #   Make molecule list
    ##############################
    nmol = pArray[end].molNum
    molTracker = Vector{Molecule}(undef,nmol)
    mcount = 0
    acount = 0
    iter = 2
    "an ugly while loop for making molecule list. I do not like it, but it does work"
    while true && mcount < nmol && pArray[iter-1] != pArray[end]
        mcount += 1
        COM = []
        mass=[]
        start = iter - 1
        while pArray[iter].molNum == pArray[iter-1].molNum
            acount += 1
            push!(COM,pArray[iter-1].coords)
            push!(mass,pArray[iter-1].mass)
            iter += 1
            if  pArray[iter-1] == pArray[end] break end # break if we are at the end of the array
        end
        push!(COM,pArray[iter-1].coords)
        push!(mass,pArray[iter-1].mass)

        centerOfMass = Center_of_Mass(COM,mass)

        ed = iter - 1
        mol = FindMolType(String(atomsPDB.resnm[acount]),moleculeList )
        molTracker[mcount] = Molecule( start, ed, centerOfMass, 0.0, mol)
        iter += 1

    end
    soa = StructArray(pArray)
    moa = StructArray(molTracker)
    return  soa, moa
end # MakeAtomArrays
################################################################################
#
#                 Make parameter tables
#                   i.e., epsilon, sigma, (alpha if EXP-6)
#
################################################################################


function MakeTables(systemTop::FFParameters,atomsPDB::Topology)
    x = length(systemTop.atomTypes)
    ϵ = [systemTop.atomTypes[i].ϵ for i=1:x]
    σ = [systemTop.atomTypes[i].σ for i=1:x]
    #=
    ϵ_ij = Array{Float64,2}(undef,x,x)
    σ_ij = Array{Float64,2}(undef,x,x)
    α_ij = Array{Float64,2}(undef,x,x)
    qqTable = Array{Float64,2}(undef,x,x)
    qqholder = []
    for i = 1:length(systemTop.molParams)
        for j = 1:length(systemTop.molParams[i].atoms)
            push!(qqholder,systemTop.molParams[i].atoms[j].charge)
        end
    end

    for i=1:x
        for j=1:x
            ϵ_ij[i,j] = sqrt(systemTop.atomTypes[i].ϵ * systemTop.atomTypes[j].ϵ)
            #ϵ_ij[i,j] = (systemTop.atomTypes[i].ϵ + systemTop.atomTypes[j].ϵ) / 2
            σ_ij[i,j] = (systemTop.atomTypes[i].σ + systemTop.atomTypes[j].σ) / 2
            α_ij[i,j] = 1.1
            qqTable[i,j]= qqholder[i] * qqholder[j]
        end
    end

    vdwTable = Tables([ϵ_ij[i,j] for i=1:size(ϵ_ij,1),j=1:size(ϵ_ij,1)],
                [σ_ij[i,j] for i=1:size(σ_ij,1),j=1:size(σ_ij,1)]) #,
                #[α_ij[i,j] for i=1:size(α_ij,1),j=1:size(α_ij,1)])
    =#
    vdwTable = Tables(ϵ, σ)
    push!(warnings,"Only epsilon and sigma have been defined, not alpha.")

    IntraBond = Bonds[]
    IntraAngle = Angles[]
    IntraDihedral = Dihedrals[]

    atomcount = 0
    skip = 0

    for i=1:(length(atomsPDB.resnr)-1) # holds last molecules index
        if skip > 0
            skip -= 1
            continue
        end


        molNum = FindMolType(String(atomsPDB.resnm[i]),moleculeList)   # Int64
        for j=1:length(systemTop.molParams[molNum].bonds)
            deepCopyBond = deepcopy(systemTop.molParams[molNum].bonds[j])
            push!(IntraBond,deepCopyBond)
            index = length(IntraBond)
            ai = systemTop.molParams[molNum].bonds[j].ai
            aj = systemTop.molParams[molNum].bonds[j].aj
            IntraBond[index].ai = atomcount + ai
            IntraBond[index].aj = atomcount + aj
        end
        for j=1:length(systemTop.molParams[molNum].angles)
            deepCopyAngle = deepcopy(systemTop.molParams[molNum].angles[j])
            push!(IntraAngle,deepCopyAngle)
            index = length(IntraAngle)
            ai = systemTop.molParams[molNum].angles[j].ai
            aj = systemTop.molParams[molNum].angles[j].aj
            ak = systemTop.molParams[molNum].angles[j].ak
            IntraAngle[index].ai = atomcount + ai
            IntraAngle[index].aj = atomcount + aj
            IntraAngle[index].ak = atomcount + ak
        end
        for j=1:length(systemTop.molParams[molNum].dihedrals)
            deepCopyDihedral = deepcopy(systemTop.molParams[molNum].dihedrals[j])
            push!(IntraDihedral,deepCopyDihedral)
            index = length(IntraDihedral)
            ai = systemTop.molParams[molNum].dihedrals[j].ai
            aj = systemTop.molParams[molNum].dihedrals[j].aj
            ak = systemTop.molParams[molNum].dihedrals[j].ak
            al = systemTop.molParams[molNum].dihedrals[j].al
            IntraDihedral[index].ai = atomcount + ai
            IntraDihedral[index].aj = atomcount + aj
            IntraDihedral[index].ak = atomcount + ak
            IntraDihedral[index].al = atomcount + al
        end
        atomcount += length(systemTop.molParams[molNum].atoms)
        skip = length(systemTop.molParams[molNum].atoms) -1 # skip all the atoms we just did for this molecule
    end

    intraFF = IntraForceField( [ IntraBond[i] for i = 1:length(IntraBond) ],
                [ IntraAngle[i] for i = 1:length(IntraAngle) ],
                [ IntraDihedral[i] for i = 1:length(IntraDihedral) ] )

    " make non_bonded list"
    n = length(atomsPDB.resnr)
    nonbonded_matrix = [ true for i=1:n, j=1:n ]# where i==j is false]
    for i=1:n
        nonbonded_matrix[i,i] = false
    end

    "exclude bonded atoms from non-bonded interactions"
    for i=1:length(intraFF.bonds)
        ai = intraFF.bonds[i].ai
        aj = intraFF.bonds[i].aj
        nonbonded_matrix[ai,aj] = false
        nonbonded_matrix[aj,ai] = false
    end

    "exclude atoms sharing an angle from non-bonded interactions"
    for i=1:length(intraFF.angles)
        ai = intraFF.angles[i].ai
        ak = intraFF.angles[i].ak
        nonbonded_matrix[ai,ak] = false
        nonbonded_matrix[ak,ai] = false
    end

    "exclude atoms sharing a 1-4 dihedral from non-bonded interactions"
    scaled_pairListTemp = []
    for i=1:length(intraFF.dihedrals)
        ai = intraFF.dihedrals[i].ai
        al = intraFF.dihedrals[i].al
        nonbonded_matrix[ai,al] = false
        nonbonded_matrix[al,ai] = false
        push!(scaled_pairListTemp,[ai,al])   # calculate 1-4 interactions seperately
    end
    scaled_pairs = Vector{SVector{2,Int64}}(undef, length(scaled_pairListTemp))
    for i=1:length(scaled_pairListTemp)
        scaled_pairs[i] = SVector(scaled_pairListTemp[i]...)
    end

    return intraFF, vdwTable, nonbonded_matrix, scaled_pairs # qqTable
end # end MakeTables

function readNIST(filename)

    LJcoords = []
    QQcoords = []
    charge = []
    atomType = []
    atomName = []
    molName = []
    num = 7
    println(filename)
    open(filename) do file
        for line in eachline(file)

            if length(split(strip(line) ) ) > 2
                line = strip(line)

                push!(QQcoords,
                    [
                    parse(Float64,split(line)[2]),
                    parse(Float64,split(line)[3]),
                    parse(Float64,split(line)[4])
                    ])
                push!(molName, "WAT")
                if split(line)[5] == "O"
                    num = 7
                    push!(atomName, "O1")
                    push!(atomType,1)
                    push!(charge, -2 * 0.42380)
                    push!(LJcoords,
                        [
                        parse(Float64,split(line)[2]),
                        parse(Float64,split(line)[3]),
                        parse(Float64,split(line)[4])
                        ])
                elseif split(line)[5] == "H"
                    num += 1
                    push!(atomName,"H" * string(num) )
                    push!(atomType,2)
                    push!(charge, 0.42380)

                end
            end # if line has 5 columns
        end # loop over all lines in files
    end # open files
    # Make Vector of Structs
    indx = 1
    ind = 1
    L = Int(floor(length(atomType)/3) )
    SoA = Vector{Molecule}(undef,L)
    SoArray = Vector{Particle}(undef,L*3+1)
    atomTracker = Vector{MVector{2,Int64}}(undef,L)

    for m = 1:L
        SoA[m] = Molecule(1,
                [atomType[i] for i=ind:(ind + 2)],
                [ SVector(LJcoords[m][1],LJcoords[m][2],LJcoords[m][3]) ],
                [charge[i] for i=ind:(ind + 2) ],
                [SVector(QQcoords[i][1],QQcoords[i][2],QQcoords[i][3]) for i=indx:(indx+2)]
                )
        atomTracker[m] = MVector(indx,indx+2)

        indx += 3
        ind += 3

    end
    molnum=1
    for m = 1:length(charge)
        SoArray[m] = Particle(molnum,
                atomType[m],
                SVector(QQcoords[m][1],QQcoords[m][2],QQcoords[m][3]),
                charge[m],
                )
        if m % 3 == 0 molnum += 1 end
    end
    SoArray[length(charge) + 1] = Particle(1,
            2,
            SVector(0.0,0.0,0.0),0.0 )
    qq_q = [charge[i] for i in eachindex(charge)]
    qq_r = [SVector{3,Float64}(QQcoords[i]...) for i=1:length(QQcoords)] #Vector{SVector{3,Float64}}(undef,length(QQcoords))
    return SoA, qq_r, qq_q, SoArray, atomTracker
end # function

"Tables to prevent numerical overflow from atoms/charges overlapping"
function OverflowTable(vdwTable,qqTable, kappa)
    push!(logFile,"Start of vdw interaction overflow cutoffs")
    num_atom_types = size(vdwTable.ϵij, 1)
    LJCutoff = zeros(num_atom_types,num_atom_types)
    LJCutoff² = zeros(num_atom_types,num_atom_types)
    @inbounds for i = 1:num_atom_types
        @inbounds for j = 1:num_atom_types
            ϵ = vdwTable.ϵij[i,j]
            σ = vdwTable.σij[i,j]
            if abs(ϵ) > 0.00 && abs(σ) > 0.00 # the abs() is overkill and unphysical
                LJCutoff[i,j] = Newtonf(ϵ,σ,overflow )
                LJCutoff²[i,j] = LJCutoff[i,j]^2
                push!(logFile,"Newton result for interaction: $i - $j $(LJCutoff²[i,j]) ")
            end
        end
    end
    qqCutoff = zeros(size(qqTable,1),size(qqTable,1))
    qqCutoff² = zeros(size(qqTable,1),size(qqTable,1))
    for i=1:size(qqTable,1)
        for j=1:size(qqTable,1)
            qqCutoff[i,j] = Newtonf(qqTable[i,j], kappa, qqOverflow ,"qq") # located in auxilaryMaths.jl
            qqCutoff²[i,j] = qqCutoff[i,j]^2
        end
    end
    len = size(LJCutoff²,1)
    overflowCutoff² = OverFlowTable([LJCutoff²[i,j] for i=1:len , j=1:len], # located in structs.jl
                                    [qqCutoff²[i,j] for i=1:len , j=1:len])
    push!(logFile,"End of vdw interaction overflow cutoffs")
    return overflowCutoff², num_atom_types
end
#end # module

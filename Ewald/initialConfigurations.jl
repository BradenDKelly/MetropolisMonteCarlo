#####
#
#          Starting Configurations
#
############
include("auxillary.jl")

export InitCubicGrid, PrintPDB, ReadCNF

function InitCubicGrid(n::Int, rho::Real)

    #!------------------------------------------------------------------------
    #! Created by Braden Kelly
    #!------------------------------------------------------------------------
    #! Creates an initial configuration
    #!------------------------------------------------------------------------
    #! input:  n       number of particles
    #!         rho     density
    #! output: coords  coordinates of n particles
    #!------------------------------------------------------------------------

    #! Calculate box length (L)
    L = (n / rho)^(1.0 / 3.0)

    #! Calculate the lowest perfect cube that will contain all of the particles
    nCube = 2

    while (nCube^3 < n)
        nCube = nCube + 1
    end
    coords = Array{Float64,2}(undef, 3, n)
    # initial position of particle 1
    posit = [0, 0, 0]

    #!begin assigning particle positions
    for i = 1:n
        coords[:, i] = (posit + [0.01, 0.01, 0.01]) * (L / nCube)
        #! Advancing the index (posit)
        posit[1] = posit[1] + 1
        if (posit[1] == nCube)
            posit[1] = 0
            posit[2] = posit[2] + 1
            if (posit[2] == nCube)
                posit[2] = 0
                posit[3] = posit[3] + 1
            end
        end
    end
    return [
        SVector{3}(coords[1, i], coords[2, i], coords[3, i])
        for i = 1:size(coords, 2)
    ]
end #function

""" Generates .pdb files """
function PrintPDB(r, box, step = 1, filename = "pdbOutput")

    open(filename * "_" * string(step) * ".pdb", "w") do file

        line = @sprintf(
            "%-7s %7.3f %7.3f %7.3f %30s \n",
            "CRYST1",
            10.0 * box,
            10.0 * box,
            10.0 * box,
            "90.00  90.00  90.00 P 1           1"
        )
        write(file, line)
        mol = 1
        atomName = "O"
        f = 1.0  # scaling factor for coordinates
        for (i, atom) in enumerate(r)
            if i == 2
                atomName = "H"
            end

            molName = "Sol"
            #atomName = systemTop.molParams[soa.mt[i]].atoms[soa.at[i]].atomnm
            #molName = systemTop.molParams[soa.mt[i]].atoms[soa.at[i]].resnm

            line = @sprintf(
                "%-6s %4d %3s %4s %5d %3s %7.3f %7.3f %7.3f %5.2f %5.2f \n",
                "ATOM",
                mol, #i,
                atomName,
                molName,
                mol,
                " ",
                f * r[i][1],
                f * r[i][2],
                f * r[i][3],
                1.00,
                0.00
            )
            write(file, line)
            if i % 3 == 0
                mol += 1
                atomName = "O"
            else
                atomName = "H"
            end
        end
    end
end

function PrintPDB(soa::StructArray, boxSize, step=1, filename="pdbOutput")
    df = 1.0  # disctance conversion if needed to go from nm to Angstrom
    open(filename * "_" * string(step) * ".pdb", "w") do file

        line = @sprintf("%-7s %7.3f %7.3f %7.3f %30s \n","CRYST1", df * boxSize[1],
               df * boxSize[2], df * boxSize[3], "90.00  90.00  90.00 P 1           1")
        write(file,line)

        for (i,atom) in enumerate(soa)

            atomName = systemTop.atomTypes[soa.atype[i]].name
            #atomName = systemTop.molParams[soa.molType[i]].atoms[soa.atype[i]].atomnm
            #molName = systemTop.molParams[soa.molType[i]].atoms[soa.atype[i]].resnm
            molName = systemTop.molParams[soa.molType[i]].name

            line = @sprintf("%-6s %4d %3s %4s %5d %3s %7.3f %7.3f %7.3f %5.2f %5.2f \n",
            "ATOM",i, atomName, molName, soa.molNum[i], " ", df * soa.coords[i][1],
            df * soa.coords[i][2],df * soa.coords[i][3], 1.00, 0.00   )
            write(file,line)
        end
    end
end

function PrintOutput(
    system::Requirements,
    totProps::Properties2,
    atomType,
    atomName,
    qq_r,
    qq_q,
    box,
    step = 1,
    filename = "xyz_quat",
)

    open(filename * "_" * string(step) * ".pdb", "w") do file

        line = @sprintf("%-7s %7.3f %7.3f %7.3f\n", "Output", box, box, box,)
        write(file, line)
        write(file, " Molecular coordinates and quaternions\n")
        write(
            file,
            " #, mol name, atom Start, atom End, x, y, z, q0, q1, q2, q3\n",
        )
        for i = 1:length(system.rm)
            line = @sprintf(
                "%4d %-7s %4d %4d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
                i,
                system.molNames[i],
                system.thisMol_theseAtoms[i][1],
                system.thisMol_theseAtoms[i][2],
                system.rm[i][1],
                system.rm[i][2],
                system.rm[i][3],
                totProps.quat[i][1],
                totProps.quat[i][2],
                totProps.quat[i][3],
                totProps.quat[i][4]
            )
            write(file, line)
        end
        write(file, "Atom coordinates\n")
        write(file, "#, name, type, charge, x, y, z\n")
        for i = 1:length(system.ra)
            line = @sprintf(
                "%4d %-7s %-7s %7.3f %7.3f %7.3f %7.3f \n",
                i,
                system.atomNames[i],
                system.atomTypes[i],
                qq_q[i],
                system.ra[i][1],
                system.ra[i][2],
                system.ra[i][3]
            )
            write(file, line)
        end
    end
end

function ReadCNF(input = "cnf_input.inp")
    r = []
    e = []
    i = 0
    box1 = 0.0
    #box=9.42953251
    open(input) do file
        for line in eachline(file)
            i += 1

            if i == 2 #length(split(line)) == 1  && split(line) != typeof(Flo)
                box1 = parse(Float64, strip(line)) # 9.42953251 #parse(Float64, split(line)) # this should be on the 2nd Line
                #println("hardcoded box at line 257 in ReadCNF: ", box1)
            end

            if i >= 3 #length(split(line)) > 4

                lin = split(line)

                push!(
                    r,
                    [
                        parse(Float64, strip(lin[1])),
                        parse(Float64, strip(lin[2])),
                        parse(Float64, strip(lin[3])),
                    ],
                )

                push!(
                    e,
                    [
                        parse(Float64, strip(lin[4])),
                        parse(Float64, strip(lin[5])),
                        parse(Float64, strip(lin[6])),
                        parse(Float64, strip(lin[7])),
                    ],
                )
            end
        end
    end
    return r, e, box1
end

function ReadNIST(filename)

    LJcoords = []
    QQcoords = []
    charge = []
    atomType = []
    atomName = []
    molName = []
    num = 7
    println(filename)
    i = 0
    box1 = 0
    open(filename) do file
        for line in eachline(file)
            i += 1
            if i == 1 #length(split(line)) == 1  && split(line) != typeof(Flo)
                box1 = parse(Float64, split(line)[1]) # 9.42953251 #parse(Float64, split(line)) # this should be on the 2nd Line
            end

            if length(split(strip(line))) > 2 && i > 2
                line = strip(line)
                push!(
                    QQcoords,
                    [
                        parse(Float64, split(line)[2]),
                        parse(Float64, split(line)[3]),
                        parse(Float64, split(line)[4]),
                    ],
                )
                push!(molName, "WAT")
                if split(line)[5] == "O"
                    num = 7
                    push!(atomName, "O1")
                    push!(atomType, 1)
                    push!(charge, -2 * 0.42380)
                    push!(
                        LJcoords,
                        [
                            parse(Float64, split(line)[2]),
                            parse(Float64, split(line)[3]),
                            parse(Float64, split(line)[4]),
                        ],
                    )
                elseif split(line)[5] == "H"
                    num += 1
                    push!(atomName, "H" * string(num))
                    push!(atomType, 2)
                    push!(charge, 0.42380)

                end
            end # if line has 5 columns
        end # loop over all lines in files
    end # open files
    L = Int(floor(length(atomType) / 3))
    indx = 1
    atomTracker = Vector{MVector{2,Int64}}(undef, L)
    rm = []   # center-of-mass of the water molecule
    for m = 1:L
        a = []
        push!(a, SVector(QQcoords[indx]...))     # oxygen
        push!(a, SVector(QQcoords[indx+1]...))   # hydrogen
        push!(a, SVector(QQcoords[indx+2]...))   # hydrogen
        mass = [15.99, 1.009, 1.009]   # [O, H, H]
        push!(rm, COM(a, mass))

        atomTracker[m] = MVector(indx, indx + 2)

        indx += 3
    end

    qq_q = [charge[i] for i in eachindex(charge)]
    qq_r = [SVector{3,Float64}(QQcoords[i]...) for i = 1:length(QQcoords)] #Vector{SVector{3,Float64}}(undef,length(QQcoords))
    return qq_r, qq_q, rm, qq_r, atomTracker, box1, atomName, atomType
end # function

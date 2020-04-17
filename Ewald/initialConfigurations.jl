#####
#
#          Starting Configurations
#
############

export InitCubicGrid, PrintPDB, ReadCNF

function InitCubicGrid(n::Int,rho::Real)

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
        L = (n/rho)^(1.0/3.0)

        #! Calculate the lowest perfect cube that will contain all of the particles
        nCube = 2

        while (nCube^3 < n)
           nCube = nCube+1
        end
        coords = Array{Float64,2}(undef, 3, n)
        # initial position of particle 1
        posit=[0, 0, 0]

        #!begin assigning particle positions
        for i=1:n
           coords[:,i] = (posit+[0.01,0.01,0.01])*(L/nCube)
           #! Advancing the index (posit)
                posit[1]=posit[1]+1
                if (posit[1] == nCube)
                   posit[1] = 0
                   posit[2] = posit[2]+1
                   if (posit[2] == nCube)
                      posit[2] = 0
                      posit[3] = posit[3] + 1
                   end
                end
        end
    return [SVector{3}(coords[1,i],coords[2,i],coords[3,i]) for i = 1:size(coords,2)]
end #function

""" Generates .pdb files """
function PrintPDB(r, box, step=1, filename="pdbOutput")

    open(filename * "_" * string(step) * ".pdb", "w") do file

        line = @sprintf("%-7s %7.3f %7.3f %7.3f %30s \n","CRYST1", 10.0 * box,
               10.0 * box, 10.0 * box, "90.00  90.00  90.00 P 1           1")
        write(file,line)
        mol = 1
        atomName = "O"
        f = 1.0  # scaling factor for coordinates
        for (i,atom) in enumerate(r)
            if i==2
                atomName = "H"
            end

            molName = "Sol"
            #atomName = systemTop.molParams[soa.mt[i]].atoms[soa.at[i]].atomnm
            #molName = systemTop.molParams[soa.mt[i]].atoms[soa.at[i]].resnm

            line = @sprintf("%-6s %4d %3s %4s %5d %3s %7.3f %7.3f %7.3f %5.2f %5.2f \n",
            "ATOM",i, atomName, molName, mol, " ", f*r[i][1],
            f * r[i][2],f * r[i][3], 1.00, 0.00   )
            write(file,line)
            if i % 3 == 0
                mol += 1
                atomName = "O"
            else
                atomName = "H"
            end
        end
    end
end

function ReadCNF(input="cnf_input.inp")
    r = []
    e = []
    i = 0
    box1=0.0
    #box=9.42953251
    open(input) do file
        for line in eachline(file)
            i += 1

            if i == 2 #length(split(line)) == 1  && split(line) != typeof(Flo)
                box1 = parse(Float64, strip(line) ) # 9.42953251 #parse(Float64, split(line)) # this should be on the 2nd Line
                println("hardcoded box at line 257 in ReadCNF: ", box1)
            end

            if i >= 3 #length(split(line)) > 4

                lin = split(line)

                push!(r, [ parse(Float64, strip(lin[1])),
                    parse(Float64, strip(lin[2])),
                    parse(Float64, strip(lin[3])) ] )

                push!(e, [ parse(Float64, strip(lin[4])),
                    parse(Float64, strip(lin[5])),
                    parse(Float64, strip(lin[6])),
                    parse(Float64, strip(lin[7])) ] )
            end
        end
    end
    return r, e, box1
end

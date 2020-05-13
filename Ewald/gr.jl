function makeRDF(pdbdir, side, numbins, cm=False, numAT=0):

    maxbin = numbins # number of bins
    sideh = side/2.0
    dr = float(sideh/maxbin) # bin width
    hist = [0]*(maxbin+1) # initialize list
    rdf = {}
    count = 0
    pdbs = listdir(pdbdir) # list of files
    nstep = len(pdbs)

    println( "Directory "+pdbdir+" has "+str(len(pdbs))+" files.")

    # loop over pdb files
    for pdb in pdbs
        if not pdb.endswith(".pdb")
            nstep -= 1
            continue # skip other files
        count += 1
        print "Reading file ... "+pdb+" ("+ str(count) +"/"+ str(len(pdbs)) +")"

        # read atom coordinates from PDB
        atoms = []
        cmm = [0,0,0] # center of mass of molecule
        atc = 1 # atom counter
        lines = open(pdbdir+"/"+pdb)
        mass = [15.999, 0.000, 0.000]
        for line in lines:
            if line[0:4] != "ATOM":
                continue # skip other lines
            coords = map(float, (line[31:54]).split())

            if cm == true # calculate center of mass
                # we assume masses of all particles are 1.0 (not true, altered BDK May 10, 2020)
                cmm[0] += coords[0] * mass[atc-1]
                cmm[1] += coords[1] * mass[atc-1]
                cmm[2] += coords[2] * mass[atc-1]
                if atc < numAT
                    atc += 1
                else:
                    atc = 1
                    cmm[0] /= sum(mass) #numAT
                    cmm[1] /= sum(mass) #numAT
                    cmm[2] /= sum(mass) #numAT

                    # fold coordinates
                    for i in range(3)
                        tmp = floor(cmm[i] * (1.0/side))
                        cmm[i] -= tmp * side
                    end

                    atoms.append((cmm[0],cmm[1],cmm[2]))
                    cmm = [0,0,0]
                end
            else # no cm calculation
                atoms.append((coords[0], coords[1], coords[2]))
            end
        end

        # loop over particle pairs
        npart = len(atoms)
        print " looping over particle pairs (" +str(npart)+ "^2) ... "
        for i in range(npart)

            xi = (atoms[i])[0]
            yi = (atoms[i])[1]
            zi = (atoms[i])[2]

            for j in range(i+1, npart)
                xx = xi - (atoms[j])[0]
                yy = yi - (atoms[j])[1]
                zz = zi - (atoms[j])[2]

                # minimum image
                if (xx < -sideh)   xx = xx + side end
                if (xx > sideh)    xx = xx - side end
                if (yy < -sideh)  yy = yy + side end
                if (yy > sideh)   yy = yy - side end
                if (zz < -sideh)   zz = zz + side end
                if (zz > sideh)   zz = zz - side end

                # distance between i and j
                rd  = xx * xx + yy * yy + zz * zz
                rij = sqrt(rd)

                bin = int(ceil(rij/dr)) # determine in which bin the distance falls
                if (bin <= maxbin)
                    hist[bin] += 1
                end
            end
        end
    end

    # normalize
    print "Normalizing ... "
    phi = npart/pow(side, 3.0) # number density (N*V)
    norm = 2.0 * pi * dr * phi * nstep * npart

    for i in range(1, maxbin+1)
        rrr = (i - 0.5) * dr
        val = hist[i]/ norm / ((rrr * rrr) + (dr * dr) / 12.0)
        rdf.update({rrr:val})
    end

    return rdf
end

#-------------------------------------------------------------------#

# write RDF into file
boxsize = 31.1448
numbins = 384 # number of bins
cm = True # calculate RDF from center of mass of molecule
numAT = 3 # number of atoms in molecule
pdbsdir = os.getcwd() #"./pdbs/" # directory with PDB files
outfile = "./rdf.out"

rdf = makeRDF(pdbsdir, boxsize, numbins, cm, numAT)
print "Writing output file ... " +outfile
outfile = open(outfile, "w")
for r in sorted(rdf.iterkeys()): # sort before writing into file
    outfile.write("%15.8g %15.8g\n"%(r, rdf[r]))
outfile.close()

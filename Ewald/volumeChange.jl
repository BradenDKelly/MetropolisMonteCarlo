######
#
#                       Volume Change Routine(s)
#
#########


"""
function VolumeChange()
    #Subroutine MC_vol(vol_attempt,vol_accept,vmax,L)
    #!coords,atom_coords,,
    # !use MT19937 ! Uses the Mersenne Twister random number generator grnd()
    #use random_generators
    #use Global_Variables
    #use energy_routines

    #implicit none

    #integer, intent(inout)                   :: vol_accept, vol_attempt
    #real, intent(inout)                      :: L
    #! real, dimension(num_mol)			     :: Xb,Xnew,Yb,Ynew,Zb,Znew
    #real, dimension(3,num_mol)		         :: coords_new
    #real, intent(in)                         :: vmax
    #!integer, intent(in)                      :: Em
    #real                                     :: vol_old,vol_new
    #!real                                     :: lnvol_new
    #real                                     :: L_new,L_old,energy_old,vir_old
    #real                                     :: test,f,energy_new,vir_new,rho_new,coru_new
    #!integer                                  :: Nb
    #!real                                     :: BoxSize
    #integer	                                 :: i,j, start,endd,begin,finish

    # ! cori is a vector containing the energy tail correction for each molecule type
    #! TypeAtom is a matrix containing the types of each atom in each molecule
    #! maxlenthmol is the number of atoms in the largest molecule
    #! atoms_per_molecule a vector containing the number of atoms in each type of molecule
    #! i.e. if CO2 is molecule type #1 then the first entry in atoms_per_molecule is 3 (C+O+O = 3 atoms)
    #! real, dimension(num_atom_types)        :: cori
    #! integer, dimension(maxlengthmol,nm)     :: TypeAtom
    #! integer	                             :: maxlengthmol
    #! integer, dimension(num_mol_types)                  :: atoms_per_molecule
    #! integer, intent(in)									:: FF_Flag ! 1 = Lennard-Jones, 2 = EXP-6
    #! ALPHA is the matrix containing the cross parameters of the EXP-6 ALPHA parameters
    #!real, dimension (num_atom_types,num_atom_types), intent(in)							:: ALPHA
    #real	                        :: enn, virn
    #real, dimension(3,num_mol)                               :: change_in_coords
    #real, dimension(3,num_atoms)                             ::atom_XYZ !,atom_YT,atom_ZT
    #!==============================================================================================

      vol_attempt = vol_attempt + 1

      L_old = L
      vir_old = vir
      energy_old = energy

      vol_old = L_old^3

    #!  lnvol_new = LOG(vol_old) + (rranf()-0.5)*LOG(vmax) #! if using lnV modify acceptance to (num_mol + 1)
      vol_new = vol_old + (rand()-0.5)*vmax #!EXP(lnvol_new)
      L_new = vol_new^(1.0/3.0)

      f = L_new/L_old #! Scaling factor

      coords_new(:,:) = f*coords(:,:)

     change_in_coords(1,:) = coords_new(1,:) - coords(1,:)
     change_in_coords(2,:) = coords_new(2,:) - coords(2,:)
     change_in_coords(3,:) = coords_new(3,:) - coords(3,:)

      for i = 1:num_mol

        start = ThisMol_startAtom_endAtom(i,1) #! from the array of atoms, this is the atom that starts molecule 1
        finish = ThisMol_startAtom_endAtom(i,2)

         for j = start:finish

             atom_XYZ(:,j) = atom_coords(:,j) + change_in_coords(:,i)

         end
     end
      #!Need to insert IF Statement for SHIFT
      start = 0;finish=0;endd=0;begin=0
      energy_new = 0.0
      vir_new = 0.0
      enn = 0.0
      virn = 0.0
     #!=======================================================================================================
     #!				Calculate total energy of new configuration
     #!=======================================================================================================

       for i = 1:num_mol-1

        start = ThisMol_startAtom_endAtom(i,1) #! from the array of atoms, this is the atom that starts molecule 1
        finish = ThisMol_startAtom_endAtom(i,2) #! from the array of atoms, this is the atom that ends molecule 1

        for j = i+1:num_mol

    #!        if(i == j) cycle ! That's right, skidaddle if you can't be different from i!

            begin = ThisMol_startAtom_endAtom(j,1) #! from the array of atoms, this is the atom that starts molecule 2
            endd = ThisMol_startAtom_endAtom(j,2) #! from the array of atoms, this is the atom that ends molecule 2

            call ener_single(finish-start+1,atom_XYZ(:,start:finish),coords_new(:,i),endd-begin+1,
                atom_XYZ(:,begin:endd),coords_new(:,j),L_new,enn,virn,
    	        atom_type_list(start:finish),atom_type_list(begin:endd)) #!energy of new coord

            energy_new = energy_new + enn
            vir_new = vir_new + virn

        end
     end


     #!===============================================================================================
     rho_new = num_mol/vol_new
     #!===============================================================================================
     #!				Add tail correction to potential energy
     #!===============================================================================================

      if( cut_type == 'tail_corr')
         call ener_corr(rho_new,vol_new,coru_new)
         energy_new = energy_new + coru_new
      end


     #!=====================================================================================================
     #! 					Acceptance criteria
     #!=====================================================================================================
      test = exp(-beta*(P*(vol_new-vol_old)-(num_mol*log(vol_new/vol_old)/beta)
           + (energy_new-energy_old) ) )

      if (rand() < test)
          #accept volume move
          vol_accept = vol_accept + 1

          if( cut_type == 'tail_corr')
             energy_new = energy_new - coru_new
          elseif( cut_type(1:9) == 'cut_shift') then
              energy = energy_new
          end
          vir = vir_new
          #Rescale box
     #!     f = L_new/L_old ! Scaling factor
          coords(:,:) = f*coords(:,:)
          atom_coords(:,:) = atom_XYZ(:,:)
          L = L_new
      end

end #Subroutine MC_vol
"""

      implicit none
      integer i,natoms,ios,nstep,iatom
      character(len = 4) atom_name
      integer nmol,natoms_in_mol
      real(8) :: x,y,z,fx,fy,fz
      real(8), dimension(200,3) :: rr,ff
      real(8), dimension(3) :: com
      real(8), dimension(200) :: mass
      real(8) :: fmolx,fmoly,fmolz,tmolx,tmoly,tmolz
      real(8) :: aufac
      integer imol
      open(10,file = 'molecular_geometry.mannitol')
      open(20,file = 'pf.xyz')
      open(30,file = 'target_forces.dat')

      aufac = 2625.0d0 / 0.52918d0 

      call read_molecules(natoms_in_mol,mass)

      nstep = 0
      nmol = 4

      do while(.true.)
         read(20,*,iostat = ios) natoms
         if(ios.lt.0) exit
         nstep = nstep + 1
         write(*,*) nstep
         write(30,*) nstep
         read(20,*)

         do imol = 1,nmol
            fmolx = 0.0d0 
            fmoly = 0.0d0 
            fmolz = 0.0d0 
            do iatom = 1,natoms_in_mol
               read(20,*) atom_name,x,y,z,fx,fy,fz
               rr(iatom,1) = x
               rr(iatom,2) = y
               rr(iatom,3) = z

               ff(iatom,1) = fx
               ff(iatom,2) = fy
               ff(iatom,3) = fz

               fmolx = fmolx + ff(iatom,1)
               fmoly = fmoly + ff(iatom,2)
               fmolz = fmolz + ff(iatom,3)
            end do 
!     find the center of mass of this molecule
            call calc_com(rr,com,natoms_in_mol,mass)

            do iatom = 1,natoms_in_mol
               rr(iatom,1) = rr(iatom,1) - com(1)
               rr(iatom,2) = rr(iatom,2) - com(2)
               rr(iatom,3) = rr(iatom,3) - com(3)
            end do 

            tmolx = 0.0d0 
            tmoly = 0.0d0 
            tmolz = 0.0d0 
            
            do iatom = 1,natoms_in_mol
               tmolx = tmolx + rr(iatom,2) * ff(iatom,3) - rr(iatom,3) * ff(iatom,2)
               tmoly = tmoly + rr(iatom,3) * ff(iatom,1) - rr(iatom,1) * ff(iatom,3)
               tmolz = tmolz + rr(iatom,1) * ff(iatom,2) - rr(iatom,2) * ff(iatom,1)
            end do 
            write(30,*) imol,fmolx*aufac,fmoly*aufac,fmolz*aufac,tmolx*aufac,tmoly*aufac,tmolz*aufac

         end do 
      end do 
      end

      subroutine calc_com(rr,com,natoms_in_mol,mass)
      implicit none
      real(8), dimension(200,3) :: rr
      real(8), dimension(3) :: com
      real(8), dimension(200) :: mass
      real(8) :: totmass
      integer natoms_in_mol,iatom

      com = 0.0d0 
      totmass = 0.0d0 
      do iatom = 1,natoms_in_mol
         com(1) = com(1) + rr(iatom,1) * mass(iatom)
         com(2) = com(2) + rr(iatom,2) * mass(iatom)
         com(3) = com(3) + rr(iatom,3) * mass(iatom)
         totmass = totmass + mass(iatom)
      end do
      com(1) = com(1) / totmass
      com(2) = com(2) / totmass
      com(3) = com(3) / totmass

      end subroutine calc_com

      subroutine read_molecules(natoms_in_mol,mass)
      implicit none
      integer iatom,natoms_in_mol
      character(len = 4) atom_name
      character(len = 8) molname
      real(8), dimension(200) :: mass
      real(8) :: x,y,z
      integer ios

      read(10,*,iostat = ios) molname,natoms_in_mol
      do iatom = 1,natoms_in_mol
         read(10,*) atom_name,mass(iatom),x,y,z
      end do 
   
      end subroutine read_molecules

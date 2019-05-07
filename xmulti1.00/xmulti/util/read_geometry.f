      implicit none
      integer natoms,iatom,ios
      character(len=3) atname
      character(len = 8) :: molname,aa
      real(8) :: x,y,z,energy,emin
      molname = 'mannitol'

      emin = -688.351560538d0

      open(10,file ='20-lowest-conformer-coordinates.xyz')
      do while(.true.)
         read(10,*,iostat = ios) natoms
         if(ios.lt.0) exit
         read(10,*) aa,energy
         write(*,*) molname,natoms,(energy - emin)* 2625.50d0
         do iatom = 1,natoms
            read(10,*) atname,x,y,z
            write(*,*) atname,x,y,z
         end do 
      end do 
      end


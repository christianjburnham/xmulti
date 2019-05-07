      implicit none
      real(8),dimension(32) :: vec
      character(len = 3) :: atname
      character(len = 6) :: aaa
      character(len = 10) :: molname
      integer natoms,n,ios,k,i,nmax
      real(8) :: ee,emin
      real(8), dimension(2048) :: energy

      open(10,file = '20-lowest-conformer-multipoles.dat')
      molname = 'mannitol'
      natoms = 26
      read(10,*) 
      n = 0 
      emin = 1.d+23
      do while(.true.)
         read(10,*,iostat = ios) aaa,ee
         if(ee.lt.emin) emin = ee
         if(ios.lt.0) exit
         n = n + 1
         write(*,*) molname,natoms,n

         energy(n) = ee
         do i = 1,natoms
            read(10,*) atname,(vec(k),k=1,16)
            write(*,17) atname,(vec(k),k=1,16)
 17         format(a3,16f8.4)
         end do 
      end do 

      nmax = n

      do n = 1,nmax
         write(20,*) n,(energy(n) - emin) * 2625.50d0
      end do 

      end

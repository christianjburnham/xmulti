!********************************************************
!**     FIND CYCLES
!**    ring recognition code
!**    Christian J Burnham, UCD, 2017
!********************************************************

      implicit none
      real(8), dimension(3,3) :: cvec,cveci
      real(8) :: x,y,z,uu
      real(8) :: xdif,ydif,zdif,rcut2,rr2
      real(8), dimension(100,3) :: ringxyz
      real(8), dimension(80000,3) :: rox,rhy,rr
      integer natoms,nox,nhy,i,j,k,mm,ncell
      integer, dimension(100) :: ringcount
      real(8), dimension(80000,3,3) :: rwater
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer, dimension(10,80000,0:10) :: ringlist
      integer :: maxring
      integer, dimension(100) :: ring1,ring2,ring12
      integer n1,n2,n12,fused,ipdc
      integer nstep,ios
      
      ipdc = 1
      maxring = 10

      open(10,file = 'super.xyz') 

      open(23,file = 'ring3.dat')
      open(24,file = 'ring4.dat')
      open(25,file = 'ring5.dat')
      open(26,file = 'ring6.dat')
      open(27,file = 'ring7.dat')
      open(28,file = 'ring8.dat')
      open(29,file = 'ring9.dat')
      open(30,file = 'ring10.dat')

      open(43,file = 'nring3.dat')
      open(44,file = 'nring4.dat')
      open(45,file = 'nring5.dat')
      open(46,file = 'nring6.dat')
      open(47,file = 'nring7.dat')
      open(48,file = 'nring8.dat')
      open(49,file = 'nring9.dat')
      open(50,file = 'nring10.dat')

      open(60,file = 'ring_results.dat')
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     READ IN GRAPH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      nstep = 0 
      do while(.true.)
         nstep = nstep + 1
!     comment as appropriate.

!     EITHER get the graph from atomic coordinates
      call readcoordinates(nox,rwater,cvec,cveci,ipdc,ios,rox,uu,ncell)
      if(ios < 0) exit
      call get_hbondlist(nox,rwater,cvec,cveci,ipdc,neighbors,num_neighbors)

!     OR get the graph directly from graph.dat
!      call readgraph(nox,neighbors,num_neighbors)


      call findrings(nox,neighbors,num_neighbors,ringcount,ringlist,maxring)

!     clean up results - remove rings which are formed from pairs of smaller rings.

!      write(*,*) 'cleanup'
      call cleanup(ringcount,ringlist,neighbors,num_neighbors)

!     print results
      write(*,*) nstep,uu
      write(60,*) nstep,uu
      call flush(60)
      do i = 3,maxring
         write(*,*) i,ringcount(i)/ncell
         write(60,*) i,ringcount(i)/ncell
         write(40+i,*) ringcount(i)/ncell
         call flush(60)
         do j = 1,ringcount(i)
            do k = 1,i
               mm = ringlist(i,j,k)
               ringxyz(k,1) = rox(mm,1)
               ringxyz(k,2) = rox(mm,2)
               ringxyz(k,3) = rox(mm,3)
            end do 
            write(20+i,*) nstep,(ringlist(i,j,k),k=1,i),(ringxyz(k,1),ringxyz(k,2),ringxyz(k,3),k=1,i)
         end do 
      end do 
      end do 
      

      end


      subroutine readgraph(nox,neighbors,num_neighbors)
      implicit none
      integer i,j,i1,i2,ios
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer:: nox

      open(50,file = 'graph.dat')

      nox = 0 
      do while(.true.)
         read(50,*,iostat = ios) i,j
         if(ios < 0) exit
         if(i > nox) nox = i 
         if(j > nox) nox = j 
         num_neighbors(i) = num_neighbors(i) + 1
         num_neighbors(j) = num_neighbors(j) + 1
         neighbors(i,num_neighbors(i)) = j 
         neighbors(j,num_neighbors(j)) = i
      end do 

      end subroutine readgraph


      subroutine readcoordinates(nox,rwater,cvec,cveci,ipdc,ios,rox,uu,ncell)
      real(8), dimension(3,3) :: cvec,cveci
      real(8), dimension(80000,3,3) :: rwater
      integer, dimension(80000,2) :: water_index
      real(8), dimension(80000,3) :: rox,rhy,rr
      real(8) :: xdif,ydif,zdif,rr2,uu
      character(len = 2048) :: cvec_buffer
      character(len = 6), dimension(80000) :: atomicname
      integer nox,ipdc,ncell
      character(len = 1) :: atname

      ios = 0 
      
      rcut2 = 1.2d0**2 

!     read in the coordinates and pack O atom coordinates into rox vector.

      nox = 0 
      nhy = 0 
      read(10,*,iostat = ios) natoms
      read(10,'(A)',iostat = ios) cvec_buffer
      call get_cvec(cvec_buffer,cvec,cveci,uu,ncell)
!      write(*,*) cvec(1,1),cvec(1,2),cvec(1,3)
!      write(*,*) cvec(2,1),cvec(2,2),cvec(2,3)
!      write(*,*) cvec(3,1),cvec(3,2),cvec(3,3)
!     write(*,*) 'natoms = ',natoms

      if(ios .ne. 0) stop
      do i = 1,natoms
         read(10,*,iostat = ios) atomicname(i),x,y,z
         if(ios < 0) exit
         atname = atomicname(i)
         rr(i,1) = x
         rr(i,2) = y
         rr(i,3) = z
         if(atname .eq. 'O') then 
            nox = nox + 1
            rox(nox,1) = x
            rox(nox,2) = y
            rox(nox,3) = z
!            write(*,*) nox,x,y,z
         else if(atname .eq. 'H') then 
            nhy = nhy + 1
            rhy(nhy,1) = x
            rhy(nhy,2) = y 
            rhy(nhy,3) = z
         endif
      end do 

!     now find all the waters molecule indices
!     NOTE, this only needs to be on the first step, assuming the bonding doesn't change.

      do i = 1,nox
         k = 0 
         do j = 1,nhy
            xdif = rox(i,1) - rhy(j,1)
            ydif = rox(i,2) - rhy(j,2)
            zdif = rox(i,3) - rhy(j,3)

            if(ipdc.eq.1) then 
               call nearest(xdif,ydif,zdif,cvec,cveci)
            endif

            rr2 = xdif**2 + ydif**2 + zdif**2

            if(rr2 < rcut2) then 
               k = k + 1
               water_index(i,k) = j
            endif

         end do 
      end do 

      do i = 1,nox

         rwater(i,1,1) = rox(i,1)
         rwater(i,1,2) = rox(i,2)
         rwater(i,1,3) = rox(i,3)

         j = water_index(i,1)

         rwater(i,2,1) = rhy(j,1)
         rwater(i,2,2) = rhy(j,2)
         rwater(i,2,3) = rhy(j,3)

         j = water_index(i,2)

         rwater(i,3,1) = rhy(j,1)
         rwater(i,3,2) = rhy(j,2)
         rwater(i,3,3) = rhy(j,3)

      end do 


      end subroutine readcoordinates


      subroutine get_hbondlist(nox,rwater,cvec,cveci,ipdc,neighbors,num_neighbors)
      implicit none
      integer nox,i,j,ipdc
      integer, dimension(100) :: ringcount
      real(8), dimension(80000,3,3) :: rwater
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      real(8) :: xdif,ydif,zdif,rcut2,x,y,z,rr2
      real(8), dimension(3,3) :: cvec,cveci
      real(8) :: x1dif,y1dif,z1dif,x2dif,y2dif,z2dif
      real(8) :: r1dis,r2dis,costheta,ff
      real(8) :: rohdis, roodis,theta,dot
      real(8) :: theta1,theta2,theta3,theta4
      real(8) :: angfac,pi
      logical :: hbonded

      rcut2 = 3.5d0**2
      neighbors = 0 
      num_neighbors = 0 

      ringcount = 0 

      pi = 4.0 * atan(1.0d0) 

      angfac = 180.0d0 / pi 


!     build up an array of all the nearest neighbors of each O site.
      do i = 1,nox
         do j = i+1,nox

            xdif = rwater(i,1,1) - rwater(j,1,1)
            ydif = rwater(i,1,2) - rwater(j,1,2)
            zdif = rwater(i,1,3) - rwater(j,1,3)

            if(ipdc.eq.1) then 
               call nearest(xdif,ydif,zdif,cvec,cveci)
            endif
            rr2 = xdif**2 + ydif**2 + zdif**2
               
            if(rr2 < rcut2) then 
!     nearest neighbor found

!     now check the OHO angle, and only accept it as a H-bond if r_OH ^ r_OO < 30 degrees.

               roodis = dsqrt(rr2)
               hbonded = .false.

               x1dif = rwater(i,2,1) - rwater(i,1,1)
               y1dif = rwater(i,2,2) - rwater(i,1,2)
               z1dif = rwater(i,2,3) - rwater(i,1,3)
               
               if(ipdc.eq.1) then 
                  call nearest(x1dif,y1dif,z1dif,cvec,cveci)
               endif

               rohdis = dsqrt(x1dif**2 + y1dif**2 + z1dif**2)

               dot = xdif * x1dif + ydif * y1dif + zdif * z1dif
               dot = dot / (roodis * rohdis)
               theta1 = 180.0d0 - acos(dot) * angfac

               if(abs(theta1) < 30.0d0) then 
                  
                  hbonded = .true.
               else
                  
                  
                  x2dif = rwater(i,3,1) - rwater(i,1,1)
                  y2dif = rwater(i,3,2) - rwater(i,1,2)
                  z2dif = rwater(i,3,3) - rwater(i,1,3)
                  
                  if(ipdc.eq.1) then 
                     call nearest(x2dif,y2dif,z2dif,cvec,cveci)
                  endif

                  rohdis = dsqrt(x2dif**2 + y2dif**2 + z2dif**2)
                  
                  dot = xdif * x2dif + ydif * y2dif + zdif * z2dif
                  dot = dot / (roodis * rohdis)
                  theta2 = 180.0d0 - acos(dot) * angfac
                  
               if(abs(theta2) < 30.0d0) then 
                  
                  hbonded = .true.
               else
                  
                  x1dif = rwater(j,2,1) - rwater(j,1,1)
                  y1dif = rwater(j,2,2) - rwater(j,1,2)
                  z1dif = rwater(j,2,3) - rwater(j,1,3)
               
                  if(ipdc.eq.1) then 
                     call nearest(x1dif,y1dif,z1dif,cvec,cveci)
                  endif

                  rohdis = dsqrt(x1dif**2 + y1dif**2 + z1dif**2)
                  
                  dot = xdif * x1dif + ydif * y1dif + zdif * z1dif
                  dot = -dot / (roodis * rohdis)
                  theta3 = 180.0d0 - acos(dot) * angfac

                  
                  if(abs(theta3) < 30.0d0) then 
                     
                     hbonded = .true.
                  else
                     
                     x2dif = rwater(j,3,1) - rwater(j,1,1)
                     y2dif = rwater(j,3,2) - rwater(j,1,2)
                     z2dif = rwater(j,3,3) - rwater(j,1,3)
                     
                     if(ipdc.eq.1) then 
                        call nearest(x2dif,y2dif,z2dif,cvec,cveci)
                     endif

                     rohdis = dsqrt(x2dif**2 + y2dif**2 + z2dif**2)
                     
                     dot = xdif * x2dif + ydif * y2dif + zdif * z2dif
                     dot = -dot / (roodis * rohdis)
                     theta4 = 180.0d0 -  acos(dot) * angfac
                     
                     
                     if(abs(theta4) < 30.0d0) then
                        hbonded = .true.
                     endif
                  endif
               endif
            endif
            
            if(hbonded) then
               num_neighbors(i) = num_neighbors(i) + 1
               num_neighbors(j) = num_neighbors(j) + 1
               neighbors(i,num_neighbors(i)) = j 
               neighbors(j,num_neighbors(j)) = i
            endif
            
         endif
            
      end do


      end do 


      end subroutine get_hbondlist


      subroutine findrings(nox,neighbors,num_neighbors,ringcount,ringlist,maxring)

      implicit none
      integer i,j,k,l,m1,n1,m2,n2,m3,n3,m4,n4,m5,n5,m6
      integer i1,i2,i3,i4,i5,i6,i7,natoms,nox,count
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer, dimension(100) :: ringcount
      integer, dimension(10,80000,0:10) :: ringlist
      integer depth,ii,iox,npath,maxring
      integer, dimension(1000) :: pathlist

      do i = 1,10
         ringcount(i) = 0 
      end do 


      do iox = 1,nox
         depth = 0 
         pathlist = 0 
         npath = 0 
         call tree(num_neighbors,neighbors,ringcount,nox,depth,iox,pathlist,npath,ringlist,maxring)
      end do 


      end subroutine findrings


      recursive subroutine tree(num_neighbors,neighbors,ringcount,nox,depth,iox,pathlist,npath,ringlist,maxring)
      integer, dimension(80000) :: num_neighbors
      integer, dimension(80000,10) :: neighbors
      integer, dimension(100) :: ringcount
      integer nox,depth,iox,iox2,neighbor,noxlist,npath,dd,count
      integer, dimension(1000) :: oxlist,pathlist
      integer, dimension(10,80000,0:10) :: ringlist
      integer :: maxring,sum
      logical visited

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  RECURSIVE SUBROUTINE TO FIND ALL CYCLES IN GRAPH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(depth .ge. maxring) return

      pathlist(depth+1) = iox

!      write(*,*) 'path',(pathlist(k),k=1,depth+1)

      do neighbor = 1, num_neighbors(iox)
         iox2 = neighbors(iox,neighbor)


!     check if a ring has been found
         if(iox2 .eq. pathlist(1).and.depth.gt.1) then 

!     the following line ensures that not both directions around the ring are counted
            if(pathlist(depth+1).gt.pathlist(2)) then 

               ringcount(depth + 1) = ringcount(depth + 1) + 1
               count = ringcount(depth + 1)
               sum = 0 
               do k = 1,depth + 1
                  ringlist(depth+1,count,k) = pathlist(k)
                  sum = sum + pathlist(k)
               end do 
               ringlist(depth+1,count,0) = sum

            endif
         endif

!     check if this oxygen is on the list
         visited = .false.
         do j = 1,depth+1
            if(pathlist(j).eq.iox2) then 
               visited = .true.
               exit
            endif
         end do 

!     the iox2 > pathlist(1) ensures that each ring is only counted once

         if(.not.visited.and.iox2 > pathlist(1)) then
!            write(*,*) "iox = ",iox,iox2,depth
            dd = depth + 1
            call tree(num_neighbors,neighbors,ringcount,nox,dd,iox2,pathlist,npath2,ringlist,maxring)
        endif

      end do 

      return

      end subroutine tree



      subroutine cleanup(ringcount,ringlist,neighbors,num_neighbors)
      integer, dimension(100) :: ringcount,ringcount2
      integer, dimension(10,80000,0:10) :: ringlist,ringlist2
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer, dimension(100) :: ring1,ring2
      integer, dimension(1000) :: pathlist
      integer, dimension(100,100) :: shortestpath,shortestpath2
      integer :: ringsize
      integer, dimension(10,80000) :: exists
      logical shortcutfound

      exists = 1

      do i1 = 1,10
         do n1 = 1,ringcount(i1)

               ring1 = 0 
               do k1 = 1,i1
                  ring1(k1) = ringlist(i1,n1,k1)
               end do 

               do k = 1,i1
                  do l = 1,i1
                     shortestpath(k,l) = 1000
                  end do 
               end do 

               do k1 = 1,i1
                  pathlist = 0 

!     a maxdepth of only half the ringsize seems to work.
                  maxdepth = i1/2
                  ringsize = i1
                  iox = ring1(k1)
                  nn = k1
                  call tree2(maxdepth,ring1,ringsize,iox,nn,neighbors,num_neighbors,0,pathlist,shortestpath)
               end do
               
!               do k = 1,i1
!                  write(*,*) 'short',(shortestpath(k,l),l=1,i1)
!               end do 
!               write(*,*)

               shortcutfound = .false.
               do k = 1,i1
                  do l = 1,i1
                     diff = k - l
                     ii = int(abs(diff - i1 * anint(diff / real(i1))))
                     if(ii.ne.shortestpath(k,l)) shortcutfound = .true.
                     shortestpath2(k,l) = ii
                  end do 
               end do 

               if(shortcutfound) then 
!                  write(*,*) '****** SHORTCUT FOUND!! ********'
                  exists(i1,n1) = 0
               endif

         end do 
      end do 


      ringlist2 = 0 

      do i = 1,10
         nn = 0 
         do j = 1,ringcount(i)
            if(exists(i,j).eq.1) then 
               nn = nn + 1
               do k = 0,i
                  ringlist2(i,nn,k) = ringlist(i,j,k)
               end do 
            endif
         end do 
         ringcount2(i) = nn
      end do 
      
      ringlist = ringlist2
      ringcount = ringcount2



!      write(*,*) 'stopping'

!      stop
      end subroutine cleanup



      subroutine remove(n12,ring12,ringcount,ringlist,exists,imatch)
      integer n12,sum
      integer, dimension(100) :: ring12,ring0,ring1
      integer, dimension(100) :: ringcount
      integer, dimension(10,80000,0:10) :: ringlist
      integer, dimension(10,80000) :: exists
      integer imatch

      if(n12 .gt.10) return


      imatch = 0 

      sum = 0 
      do j = 1,n12
         ring1(j) = ring12(j)
         sum = sum + ring12(j)
      end do 

!     sort this ring

      call piksrt(n12,ring1)

      do i = 1,ringcount(n12)

!     check if the sum of their elements matches
         if(sum.eq.ringlist(n12,i,0)) then 

            do j = 1,n12
               ring0(j) = ringlist(n12,i,j)
            end do 
            
!     now sort second ring
            
            call piksrt(n12,ring0)
            nmatch = 0 
            do j = 1,n12
!     write(*,*) "CC",j,ring1(j),ring0(j)
               if(ring1(j) .eq. ring0(j)) nmatch = nmatch + 1
            end do 
            if(nmatch .eq. n12) then 
!     write(*,*) "MATCH!!"
               imatch = 1
               exists(n12,i) = 0
               return
            endif
!     write(*,*) 
            
         endif
      end do 
      end subroutine remove
      


      SUBROUTINE piksrt(n,arr)
      INTEGER n
      INTEGER arr(n)
      INTEGER i,j
      REAL a
      do 12 j=2,n
        a=arr(j)
        do 11 i=j-1,1,-1
          if(arr(i).le.a)goto 10
          arr(i+1)=arr(i)
 11             continue
        i=0
 10           arr(i+1)=a
 12               continue
      return
      END


      recursive subroutine tree2(maxdepth,ring,ringsize,iox,nn,neighbors,num_neighbors,depth,pathlist,shortestpath)
      integer, dimension(100) :: ring
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer, dimension(1000) :: pathlist
      integer depth
      integer, dimension(100,100) :: shortestpath
      logical visited
      integer :: ringsize

      do k = 1,ringsize
         if(iox.eq.ring(k)) then 
            if(depth.lt.shortestpath(k,nn)) then 
               shortestpath(k,nn) = depth
               shortestpath(nn,k) = depth
            endif
         endif
      end do 


      if(depth.ge.maxdepth) return

      pathlist(depth+1) = iox

      do neighbor = 1,num_neighbors(iox)
         iox2 = neighbors(iox,neighbor)


!     check if this oxygen is on the list
         visited = .false.
         do j = 1,depth+1
!            write(*,*) 'pathlist',pathlist(j),j
            if(pathlist(j).eq.iox2) then 
               visited = .true.
               exit
            endif
         end do 

         if(.not.visited) then
            newdepth = depth + 1
!            write(*,*) depth,'tree2',iox,iox2
            call tree2(maxdepth,ring,ringsize,iox2,nn,neighbors,num_neighbors,newdepth,pathlist,shortestpath)
         endif

      end do 

      end subroutine tree2


      subroutine get_cvec(buffer,cvec,cveci,uu,ncell)
      implicit none
      character(len = 2048) :: buffer
      real(8), dimension(10) :: vec
      real(8) :: amag,bmag,cmag,alpha,beta,gamma,uu
      real(8), dimension(3,3) :: cvec,cveci
      integer jcell1max,jcell2max,jcell3max,ncell
      character(len = 8) :: type
      logical okflag

      if(buffer.ne.'    ') then 
         read(buffer,*) type,vec
         if(type.eq.'XYZ') then
            cvec(1,1) = vec(1)
            cvec(1,2) = vec(2)
            cvec(2,2) = vec(3)
            cvec(1,3) = vec(4)
            cvec(2,3) = vec(5)
            cvec(3,3) = vec(6)
         else if(type.eq.'ANGLE') then 
            amag = vec(1)
            bmag = vec(2)
            cmag = vec(3)
            alpha = vec(4)
            beta = vec(5)
            gamma = vec(6)
            uu = vec(7)
            jcell1max = vec(8)
            jcell2max = vec(9) 
            jcell3max = vec(10) 
            amag = amag * jcell1max
            bmag = bmag * jcell2max
            cmag = cmag * jcell3max
            ncell = jcell1max*jcell2max*jcell3max
            call get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec)
         else
            write(*,*) 'ERROR: MUST SPECIFY EITHER XYZ&
     & OR ANGLE CVEC TYPE IN INPUT_FILE'
            stop
         endif

         call mat3inverse(cvec,cveci,okflag)
      else 
      endif

      end subroutine get_cvec


      subroutine get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec0)
      implicit none
      real(8) :: alpha,beta,gamma,amag,bmag,cmag
      real(8) :: alpha_rad,beta_rad,gamma_rad
      real(8), dimension(6) :: param
      real(8) :: angfac,pi
      real(8), dimension(3,3) :: cvec0

      pi = 4.0d0*datan(1.d0) 
      angfac = 180.0d0 / pi

      alpha_rad = alpha / angfac
      beta_rad = beta / angfac
      gamma_rad = gamma / angfac

      cvec0 = 0.0d0

      cvec0(1,1) = amag
      cvec0(1,2) = bmag * dcos(gamma_rad)
      cvec0(2,2) = dsqrt(bmag**2 - cvec0(1,2)**2)
      cvec0(1,3) = cmag * dcos(beta_rad)
      cvec0(2,3) = (bmag * cmag * dcos(alpha_rad) - cvec0(1,2) * cvec0(1,3))/cvec0(2,2)
      cvec0(3,3) = dsqrt(cmag**2 - cvec0(1,3)**2 - cvec0(2,3)**2)

      end subroutine get_cartesian_cvec

      SUBROUTINE mat3inverse (A, AINV, OK_FLAG)

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************


      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE mat3inverse

      subroutine nearest(xdif,ydif,zdif,cmat,cmati)
      implicit none
      real(8),dimension(3,3) :: cmat,cmati
      real(8) :: xdif,ydif,zdif
      real(8),dimension(3) :: vec1,vec2

!     applies the minimum image convention for triclinic cells
!     note, assumes that both cmat and cmati are upper triangular

      vec1(1) = anint(cmati(1,1) * xdif + cmati(1,2) * ydif + cmati(1,3) * zdif)
      vec1(2) = anint(cmati(2,2) * ydif + cmati(2,3) * zdif)
      vec1(3) = anint(cmati(3,3) * zdif)

      xdif = xdif - (cmat(1,1) * vec1(1) + cmat(1,2) * vec1(2) + cmat(1,3) * vec1(3))
      ydif = ydif - (cmat(2,2) * vec1(2) + cmat(2,3) * vec1(3))
      zdif = zdif - (cmat(3,3) * vec1(3))

      end subroutine nearest

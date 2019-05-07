      implicit none
      integer imol,i,n,ios,nlist,nfound,index,nf,nn
      real(8), dimension(262144) :: amag,bmag,cmag,alpha,beta,gamma,uu,dens,press
      integer, dimension(262144) :: nmol
      real(8) :: u,d
      real(8),dimension(262144,20) :: x_mol,y_mol,z_mol,phi_mol,theta_mol,psi_mol
      integer,dimension(262144,20) :: conformer
      real(8) :: cos_a,cos_b,cos_c
      real(8) :: pi,angfac,vol
      real(8) :: denstol,utol
      real(8), dimension(10) :: coord
      real(8), dimension(262144) :: ulist,denslist,u_list
      integer, dimension(262144) :: index_list,order_list
      real(8) :: d_u,d_dens,pmax
      character(len = 8) :: string1,string2,string3,string4
      character(len = 8), dimension(262144,20) :: molstring
      logical accept,print_flag,per_mol

      per_mol = .false.

      pi = 4.0d0 * datan(1.0d0)
      angfac = pi / 180.0d0 

      denstol = 0.01
      utol = 0.5d0 
      pmax = 200.0d0 

      open(10,file = 'coord.rigid')
      open(77,file = 'best.rigid')

      n = 0 
      do while(.true.)
         read(10,*,iostat = ios) nn
         if(ios.lt.0) exit
         n = n + 1
         
         nmol(n) = nn
         read(10,*) string1,amag(n),bmag(n),cmag(n),alpha(n),beta(n),gamma(n),uu(n),press(n),dens(n)

         if(per_mol) uu(n) = uu(n) / dble(nmol(n))

         write(88,*) uu(n),dens(n)
         write(89,*) uu(n),press(n)
         write(90,*) press(n),dens(n)

         ulist(n) = uu(n)
         denslist(n) = dens(n)

         read(10,*) string2,string3
         do imol = 1,nmol(n)
            read(10,*) molstring(n,imol),x_mol(n,imol),y_mol(n,imol),z_mol(n,imol) &
     &           ,phi_mol(n,imol),theta_mol(n,imol),psi_mol(n,imol),conformer(n,imol)
         end do 
      end do 

      nlist = n

      nfound = 0 
      do n = 1,nlist
         u = ulist(n)

         accept = .true.
         do nf = 1,nfound
            index = index_list(nf)
            d_u = ulist(n) - ulist(index)
            d_dens = denslist(n) - denslist(index)
            if((dabs(d_u).lt.utol.and.dabs(d_dens).lt.denstol).or.abs(press(n)).gt.pmax) then 
               accept = .false.
            exit
            endif
         end do 

         if(accept) then
            nfound = nfound + 1
            index_list(nfound) = n
            u_list(nfound) = u
         endif
      end do 


      do nf = 1,nfound
         index = index_list(nf)
         write(55,*) ulist(index),denslist(index)
      end do 
      
      
      call piksrt(nfound,u_list,index_list)
      
      do nf = 1,nfound
         index = index_list(nf)
         write(*,*) u_list(nf),index_list(nf),dens(index)
      end do 
      
      do nf = 1,nfound
         index = index_list(nf)
         write(77,*) nmol(index)
         write(77,*) string1,amag(index),bmag(index),cmag(index) &
     &,alpha(index),beta(index),gamma(index),uu(index),press(index),dens(index)
         write(77,*) string2,string3
         do imol = 1,nmol(index)
         write(77,*) molstring(index,imol),x_mol(index,imol),y_mol(index,imol),z_mol(index,imol) &
     &           ,phi_mol(index,imol),theta_mol(index,imol),psi_mol(index,imol),conformer(index,imol)
         end do 
      end do 


      write(*,*) 'nfound = ',nfound

      stop


      end


      SUBROUTINE piksrt(n,arr,index)
      INTEGER n
      REAL(8) arr(n)
      integer index(n)
      INTEGER i,j,ind
      REAL(8) a
      integer a2
      do 12 j=2,n
        a=arr(j)
        ind = index(j)
        do 11 i=j-1,1,-1
          if(arr(i).le.a)goto 10
          arr(i+1)=arr(i)
          index(i+1) = index(i)
 11             continue
        i=0
 10           arr(i+1)=a
              index(i+1) = ind
 12               continue
      return
      END

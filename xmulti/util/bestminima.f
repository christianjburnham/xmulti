      implicit none
      character(len=2048) buffer
      character(len = 2048), dimension(64) :: textlist
      real(8) :: u,u2,p,d,u1,d1,utol,dtol
      real(8), dimension(32768):: ulist,denslist
      real(8) :: g1,g2,g3,g4,g5,g6
      real(8) :: h1,h2,h3,h4,h5,h6
      integer, dimension(32768) :: index
      integer nlist,i,ios,ndata,natoms,n,nmax,naccept,j,k,nmol,nitems
      character (len = 6), dimension(1024) :: atname,molname
      character (len = 6) :: type,type2,type3,type4
      real(8), dimension(512,1024,3) :: rr
      real(8), dimension(512,1024,16) :: rr_rigid
      real(8), dimension(1024,3) :: rr0
      real(8), dimension(512,6) :: cvec
      real(8), dimension(512) :: uuu
      real(8), dimension(32768) :: uu

      logical accept

      nmax = 64
      
      open(10,file = "energy.dat")
      open(12,file = 'coord.xyz')
      open(13,file = 'coord.rigid')
      open(30,file = 'best.xyz')
      open(31,file = 'best.rigid')
      open(15,file = 'best.energy')

      utol = 0.05d0 
      dtol = 0.02d0

      nlist = 0 
      ndata = 0 
      accept = .false.
      do while(.true.)
         read(10,*,iostat = ios) u,p,d
         if(ios.lt.0) exit
         ndata = ndata + 1
         uu(ndata) = u
!     check to see if this is on the list of minima
         accept = .true.
         do i = 1,nlist
            u1 = ulist(i)
            d1 = denslist(i)

            if(abs(u1-u).lt.utol.and.(abs(d1-d).lt.dtol)) then 
               accept = .false.
               exit
            endif
         end do 
         if(accept) then 
            nlist = nlist+1
            ulist(nlist) = u
            denslist(nlist) = d
            index(nlist) = ndata
         endif
      end do 
      
      nmax = min(nmax,nlist)
      call piksrt(nlist,ulist,denslist,index)

      do i = nmax,1,-1
         write(*,*) i,ulist(i),denslist(i)
         write(15,*) i,ulist(i),denslist(i)
      end do 

      ndata = 0 
      naccept = 0
      do while(.true.)
         read(12,*,iostat = ios) natoms
         if(ios.lt.0) exit
         ndata = ndata + 1
         read(12,*) type,g1,g2,g3,g4,g5,g6,u

         read(13,*) nmol
         if(ios.lt.0) exit
         read(13,*) type2,h1,h2,h3,h4,h5,h6,u2
         read(13,*) type3,type4
!     check if it's on the list
         accept = .false.
         do n = 1,nlist
            if(index(n).eq.ndata.and.n.le.nmax) then 
               accept = .true.
               exit
            endif
         end do 
         if(accept) then 
            naccept = naccept + 1
            cvec(n,1) = g1
            cvec(n,2) = g2
            cvec(n,3) = g3
            cvec(n,4) = g4
            cvec(n,5) = g5
            cvec(n,6) = g6
            uuu(n) = u
!            write(*,*) 'xxx',naccept,uuu(n),uu(ndata)
            if(abs(uuu(n)-uu(ndata)).gt.1.e-5) write(*,*) &
     & 'warning! Possible mismatch in energy.dat and coord.xyz'

            do i = 1,natoms
               read(12,*) atname(i),rr(n,i,1),rr(n,i,2),rr(n,i,3)
            end do 
            do i = 1,nmol
               read(13,'(A)') buffer
               call split_text_to_list(buffer,textlist,nitems)
               molname(i) = textlist(1)
               do j = 2,nitems
                  read(textlist(j),*) rr_rigid(n,i,j-1)
               end do 
            end do 
         else
            do i = 1,natoms
               read(12,*) 
            end do 
            do i = 1,nmol
               read(13,*) 
            end do 
         endif


      end do 

      do n = nmax,1,-1
         write(30,*) natoms
         write(30,*) type,(cvec(n,j),j=1,6),uuu(n)
         do i = 1,natoms
            write(30,*) atname(i),rr(n,i,1),rr(n,i,2),rr(n,i,3)
         end do 
         write(31,*) nmol
         write(31,*) type2,(cvec(n,j),j=1,6),uuu(n)
         write(31,*) type3,type4
         do i = 1,nmol
            write(31,*) molname(i),(rr_rigid(n,i,j),j=1,nitems-1)
         end do 
      end do 


      end


      SUBROUTINE piksrt(n,arr,arr2,index)
      INTEGER n
      REAL(8) arr(n),arr2(n)
      integer index(n)
      INTEGER i,j,ind
      REAL(8) a,a2
      do 12 j=2,n
        a=arr(j)
        a2 = arr2(j)
        ind = index(j)
        do 11 i=j-1,1,-1
          if(arr(i).le.a)goto 10
          arr(i+1)=arr(i)
          arr2(i+1) = arr2(i)
          index(i+1) = index(i)
 11             continue
        i=0
 10           arr(i+1)=a
              arr2(i+1) = a2
              index(i+1) = ind
 12               continue
      return
      END



      subroutine split_text_to_list(buffer,textlist,n)
      character(len=2048) buffer,text1,text2,string
      character(len = 2048), dimension(64) :: textlist
      character(len=2048) :: justifyl
      integer i,n

      string = justifyl(buffer)
      i = 0 
      do while(.true.)
         call split_text(string,text1,text2)
         i = i + 1
         textlist(i) = text1
         if(text2.eq.'  ') exit
         string = text2
      end do 
      n = i 
      end

      subroutine split_text(buffer,text1,text2)
      implicit none
      character(len=2048) :: buffer,text2,text1
      integer pos,first_tab,first_space
      character(len=2048) :: justifyl

! Find the first instance of whitespace.  Split label and data.

      first_space = scan(buffer,' ')
      first_tab = scan(buffer,char(9))
      if(first_tab.eq.0) first_tab = 512
      if(first_space.eq.0) first_space = 512

      pos = min(first_tab,first_space)
      text1 = buffer(1:pos)
      text2 = buffer(pos+1:)

      text1 = justifyl(text1)
      text2 = justifyl(text2)

      end subroutine split_text


     FUNCTION JUSTIFYL(STRING) 
      IMPLICIT NONE 
      CHARACTER(LEN = *), INTENT(IN) :: STRING 
      CHARACTER(LEN = LEN(STRING)) JUSTIFYL 
      INTEGER I, L 
      L=LEN(STRING) 
      JUSTIFYL = REPEAT(' ',L)  ! Clear out any rubbish (try leaving this out) 
      I=1 
      DO WHILE ((STRING(I:I) == ' ' .OR. STRING(I:I) == CHAR(9)) &
     &          .AND. I < L)              ! Look for first non-(space or tab) character 
        I=I+1 
      END DO 
      JUSTIFYL(1:L-I+1) = STRING(I:L)       ! Shift left 
      JUSTIFYL(L-I+2:L) = REPEAT(' ',I-1)   ! Replace end with spaces 
      END FUNCTION JUSTIFYL

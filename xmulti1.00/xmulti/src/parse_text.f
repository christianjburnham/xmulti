      module parse_text

      contains

      subroutine split_text_to_list(buffer,textlist,n)
      character(len=2048) buffer,text1,text2,string
      character(len = 2048), dimension(64) :: textlist
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
      CHARACTER(LEN = 2048) JUSTIFYL 
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

      function upcase(string) result(upper)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: upper
      integer :: j
      do j = 1,len(string)
         if(string(j:j) >= "a" .and. string(j:j) <= "z") then
            upper(j:j) = achar(iachar(string(j:j)) - 32)
         else
            upper(j:j) = string(j:j)
         end if
      end do
      end function upcase

      end module parse_text

      subroutine file_doesnt_exist(filename)
      implicit none
      character(len = *) :: filename
      write(*,*) 'ERROR: FILE ',filename,' DOES NOT EXIST'
      stop
      end 

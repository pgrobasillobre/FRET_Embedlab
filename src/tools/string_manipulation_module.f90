module string_manipulation_module
!
  implicit none
  public
!
contains
!----------------------------------------------------------------------
  subroutine lower(str)
!
!   lower-case
!
    character :: str*(*)
    integer :: i
!   
    do i = 1, len(str)
      select case(str(i:i))
        case("A":"Z")
          str(i:i) = achar(iachar(str(i:i))+32)
      end select
    end do  
!    
  end subroutine lower
!-----------------------------------------------------------------------
  character(200) function sweep_blanks(in_str)
!
!   sweep_blanks : https://www.tek-tips.com/viewthread.cfm?qid=1629864
!
    character(*), intent(in) :: in_str
    character(200) :: out_str
    character :: ch
    integer :: j
!    
    out_str = " "
    do j=1, len_trim(in_str)
!    
!   get j-th char
!
      ch = in_str(j:j)
      if (ch .ne. " ".and.ch.ne.char(9)) then
        out_str = trim(out_str) // ch
      endif
      sweep_blanks = out_str
    end do
!
  end function sweep_blanks
!-----------------------------------------------------------------------
  integer function get_n_in_line(line) result(number_)
!
!   Function for getting how many elements are in a line
!
    implicit none
!
    character(len=200), intent(inout) :: line
    character(len=200)                :: line_nospaces
!
    line_nospaces = sweep_blanks(line)
    number_ = len_trim(line) - len_trim(line_nospaces) + 1
!
  end function get_n_in_line
!-----------------------------------------------------------------------
  logical function commented_line(line) result(commented)
!
!   Function to get a commented line
!
    implicit none
!
    character(len=200), intent(in) :: line
    character(len=200) :: line_adj
!
    line_adj = trim(adjustl(line))
    if ( line_adj(1:1).eq.'!'.or.line_adj(1:1).eq.'#') then
        commented = .true.
    else
        commented = .false.
    endif
!
  end function commented_line
!-----------------------------------------------------------------------
  integer function move_to_character(string,string_2) result(point)
! 
!   move to specific character 
! 
    implicit none
! 
    character(len=200), intent(inout) :: string
    character(len=*), intent(in)      :: string_2
!    
    string = adjustl(string)
    point = 1
    do while (point .lt. 200)
       if (string(point:point) .eq. string_2) then
          exit
       else
          point = point + 1
          cycle
       endif
    enddo
! 
  end function move_to_character
!----------------------------------------------------------------------
  subroutine go_to_string(IIn,string,found_string,num_lines,num_string)
!
!   go to string
!
    implicit none
!
    integer                           :: IIn
    integer                           :: num_lines
    integer, optional                 :: num_string
    character(len=*), intent(in)      :: string
    logical                           :: found_string
!
!   internal 
!
    integer            :: i
    character(len=200) :: line
!
    found_string = .false.
    num_string   = 1
!
    do i = 1, num_lines
!
       Read(IIn,'(a)') line
       call lower(line)
!
       if (trim(adjustl(line)).eq.string) then
!
          found_string = .true.
          exit
!
       endif
!     
       num_string = num_string + 1
!     
    enddo
!   
  end subroutine go_to_string
!-----------------------------------------------------------------------
  subroutine get_number_lines(IIn,nlines)
!
!   get number lines
!
    implicit none
!
    integer :: IIn
    integer :: nlines
    integer :: io
!     
    nlines = 0
    do
      read(IIn,*,iostat=io)
      if (io.ne.0) exit
      nlines = nlines + 1
    end do
! 
  end subroutine get_number_lines
!-----------------------------------------------------------------------
  subroutine check_float(IWrt,line_what,line)
!
!   check if a string contains a float
!
    implicit none
!
    integer, intent(in) :: IWrt
    character(len=200),intent(in) :: line_what
    character(len=200),intent(in) :: line
!     
    If(index(line,".").le.0) then
!
       write(IWrt,'(1x,a)') trim(line_what)//": "//trim(line)//" is not a float"
       stop
!
    endif
! 
  end subroutine check_float
!-----------------------------------------------------------------------
  pure recursive function substitute_string(string,search,substitute) result(modified_string)
!  
!     Recursive function
!     Adapted from https://stackoverflow.com/questions/58938347/
!     how-do-i-replace-a-character-in-the-string-with-another-charater-in-fortran
!
      implicit none
!      
!     input variables
!
      character(len=*), intent(in)  :: string, search, substitute
!
!     output variables
!
      character(len=:), allocatable :: modified_string
!
!     integer variables
!
      integer                       :: i, stringLen, searchLen
!      
      stringLen = len(string)
      searchLen = len(search)
!      
      if (stringLen.eq.0 .or. searchLen.eq.0) then
!      
         modified_string = ""
         return
!      
      elseif (stringLen.lt.searchLen) then
!      
         modified_string = string
         return
!      
      else 
!      
         i = 1
         do
!         
            if (string(i:i+searchLen-1).eq.search) then
                modified_string = string(1:i-1) // substitute // &
                substitute_string(string(i+searchLen:stringLen),search,substitute)
                exit
            end if
!         
            if (i+searchLen.gt.stringLen) then
                modified_string = string
                exit
            end if
!         
            i = i + 1
            cycle
!             
         end do
!      
      end if
!      
!      
  end function substitute_string
!-----------------------------------------------------------------------
end module string_manipulation_module

!----------------------------------------------------------------------
module solvent_module
!      
!   Module solvent
!
    use output_module
!!    use target_module
    use string_manipulation_module
    use parameters_module
!
    Implicit None
!
    public solvent
!
!   solvent type
!
    type solvent_type
! 
      integer                                   :: natoms
!
      real(dp), dimension(:),   allocatable  :: q
      real(dp), dimension(:,:), allocatable  :: mu
      real(dp), dimension(:,:), allocatable    :: xyz
!
      logical :: charges, dipoles
!
    end type solvent_type
!
    type (solvent_type), Save :: solvent
!
   contains
!----------------------------------------------------------------------
   subroutine read_solvent(infile, solv)
!
!    Read input solv file
!
     implicit none
!
     character(len=*),    intent(in)        :: infile
     type (solvent_type), intent(out)  :: solv
!
!
!    internal variables
!
     character(len=200) :: line
!
     integer :: nlines
     integer :: num_string_initial, num_string_end, dum
     integer :: IIn
     integer :: i,j
     integer :: iost
!
     logical :: found_string
!
     solv%charges = .false.
     solv%dipoles = .false.
!
     IIn = 21
!
     open (unit=IIn, file=trim(infile), status="old", action="read", iostat=iost,err=01)
!
!      Check number of lines
       Rewind(IIn) 
       call get_number_lines(IIn,nlines)
!
!      Check number of entries for FRET quantities
       Rewind(IIn) 
       call go_to_string(IIn,fret_start,found_string,nlines,num_string_initial)
!
       if(.not.found_string) then
          call out_%error('No FRET charges/dipoles found in "'//trim(infile)//'" input')
       endif
!
       Rewind(IIn) 
       call go_to_string(IIn,fret_end,found_string,nlines,num_string_end)
!
       if(.not.found_string) then
          call out_%error('No end for FRET charges/dipoles found in "'//trim(infile)//'" input')
       endif
!
       solv%natoms = num_string_end - num_string_initial - 2
!
!      Check if we have charges or charges + dipoles
       Rewind(IIn) 
       call go_to_string(IIn,charges_header,solv%charges,nlines,dum)

       Rewind(IIn) 
       call go_to_string(IIn,dipoles_header,solv%dipoles,nlines,dum)

       if (.not.solv%charges .and. .not.solv%dipoles) call out_%error("Charges and/or dipoles not found ")
!
       call allocate_solvent(solv)
!
       Rewind(IIn) 
       if(solv%charges) then !FQ case
          call go_to_string(IIn,charges_header,found_string,nlines,dum)
       elseif(solv%dipoles) then ! FQFMu case
          call go_to_string(IIn,dipoles_header,found_string,nlines,dum)
       endif
!
       do i = 1, solv%natoms
!
         read(IIn,'(a)') line
!       
         if(len_trim(line).eq.0) then
            call out_%error("Blank line found in corrupt FRET charges/dipoles file: "//trim(infile))
         endif
!
         !print *, line
         if(solv%charges) then
            read(line,'(5(f25.16, 5x))') solv%q(i), solv%q(i), solv%xyz(1,i), solv%xyz(2,i), solv%xyz(3,i)
         elseif(solv%dipoles) then
            read(line,'(11(f25.16, 5x))') solv%q(i),    solv%q(i),                 &
                                          solv%mu(1,i), solv%mu(2,i), solv%mu(3,i),&
                                          solv%mu(1,i), solv%mu(2,i), solv%mu(3,i),&
                                          solv%xyz(1,i),  solv%xyz(2,i),  solv%xyz(3,i)
         endif
!
       enddo
!
     01 continue
     close(IIn)
!
   end subroutine read_solvent
!----------------------------------------------------------------------
   subroutine allocate_solvent(solv)
!
!    Read input solv file
!
     implicit none
!
     type (solvent_type), intent(inout)  :: solv
!
     allocate(solv%xyz(3,solv%natoms))
!
     if (solv%charges .and. .not. solv%dipoles) then 
        allocate(solv%q(solv%natoms))

     elseif (solv%dipoles .and. .not. solv%charges) then
        allocate(solv%mu(3,solv%natoms))

     else
        call out_%error("Not recognised whether a charges or charges + dipoles is requested")
     endif
!
   end subroutine allocate_solvent
!----------------------------------------------------------------------
   subroutine delete_solvent(solv)
!
!    Read input solv file
!
     implicit none
!
     type (solvent_type), intent(inout)  :: solv
!
     deallocate(solv%xyz)
!
     if (solv%charges .and. .not. solv%dipoles) then 
        deallocate(solv%q)

     elseif (solv%dipoles .and. .not. solv%charges) then
        deallocate(solv%mu)

     endif
!
   end subroutine delete_solvent
!----------------------------------------------------------------------
end module solvent_module


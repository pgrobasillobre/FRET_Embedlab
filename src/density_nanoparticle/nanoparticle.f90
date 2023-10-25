!----------------------------------------------------------------------
module nanoparticle_module
!      
!   Module nanoparticle
!
    use output_module
!!    use target_module
    use string_manipulation_module
    use parameters_module
!
    Implicit None
!
    public nanoparticle
!
!   nanoparticle type
!
    type nanoparticle_type
! 
      integer                                   :: natoms
!
      real(dp), dimension(:,:),   allocatable  :: q
      real(dp), dimension(:,:,:), allocatable  :: mu
      real(dp), dimension(:,:), allocatable    :: xyz
!
      logical :: charges, dipoles
!
    end type nanoparticle_type
!
    type (nanoparticle_type), Save :: nanoparticle
!
   contains
!----------------------------------------------------------------------
   subroutine read_nanoparticle(infile, np)
!
!    Read input np file
!
     implicit none
!
     character(len=*),    intent(in)        :: infile
     type (nanoparticle_type), intent(out)  :: np
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
     np%charges = .false.
     np%dipoles = .false.
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
       np%natoms = num_string_end - num_string_initial - 2
!
!      Check if we have charges or charges + dipoles
       Rewind(IIn) 
       call go_to_string(IIn,charges_header,np%charges,nlines,dum)
       if (.not.np%charges) call out_%error("ONLY CHARGES ARE SUPPORTED || NOT DIPOLES YET ")
!
       call allocate_nanoparticle(np)
!
       Rewind(IIn) 
       call go_to_string(IIn,charges_header,found_string,nlines,dum)
!
       do i = 1, np%natoms
!
         read(IIn,'(a)') line
!       
         if(len_trim(line).eq.0) then
            call out_%error("Blank line found in corrupt FRET charges/dipoles file: "//trim(infile))
         endif
!
         !print *, line
         if(np%charges) then
            read(line,'(5(f25.16, 5x))') np%q(i,1), np%q(i,2), np%xyz(1,i), np%xyz(2,i), np%xyz(3,i)
         elseif(np%dipoles) then
            call out_%error("charges+dipoles not yet supported")
         endif
         !write(*,'(a,5(e20.10, 5x))'), ' ', np%q(i,1), np%q(i,2), np%xyz(1,i), np%xyz(2,i), np%xyz(3,i)
         !print *, ''
         !print *, ''
!
       enddo
!
     01 continue
     close(IIn)
!
   end subroutine read_nanoparticle
!----------------------------------------------------------------------
   subroutine allocate_nanoparticle(np)
!
!    Read input np file
!
     implicit none
!
     type (nanoparticle_type), intent(inout)  :: np
!
     allocate(np%xyz(3,np%natoms))
!
     if (np%charges .and. .not. np%dipoles) then 
        allocate(np%q(np%natoms,2))

     elseif (np%dipoles .and. .not. np%charges) then
        allocate(np%mu(3,np%natoms,2))

     else
        call out_%error("Not recognised whether a charges or charges + dipoles is requested")
     endif
!
   end subroutine allocate_nanoparticle
!----------------------------------------------------------------------
   subroutine delete_nanoparticle(np)
!
!    Read input np file
!
     implicit none
!
     type (nanoparticle_type), intent(inout)  :: np
!
     deallocate(np%xyz)
!
     if (np%charges .and. .not. np%dipoles) then 
        deallocate(np%q)

     elseif (np%dipoles .and. .not. np%charges) then
        deallocate(np%mu)

     endif
!
   end subroutine delete_nanoparticle
!----------------------------------------------------------------------
end module nanoparticle_module


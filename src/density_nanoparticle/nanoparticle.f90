!----------------------------------------------------------------------
module nanoparticle_module
!      
!   Module nanoparticle
!
!!    use output_module
!!    use target_module
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
      integer                                :: natoms
!
      real(dp), dimension(:),   allocatable  :: q
      real(dp), dimension(:,:), allocatable  :: mu
      real(dp), dimension(:,:), allocatable  :: xyz
!
      logical :: charges, dipoles
!
    end type nanoparticle_type
!
    type (nanoparticle_type), Save :: nanoparticle
!
   contains
!----------------------------------------------------------------------
!!!   subroutine read_nanoparticle(infile, np)
!!!!
!!!!    Read input np file
!!!!
!!!     implicit none
!!!!
!!!     character(len=*),    intent(in)   :: infile
!!!     type (nanoparticle_type), intent(out)  :: np
!!!!
!!!!
!!!!    internal variables
!!!!
!!!     integer                         :: IIn
!!!     integer                         :: i,j
!!!     integer                         :: iost
!!!!
!!!     IIn = 21
!!!!
!!!!    Read variables
!!!     open (unit=IIn, file=trim(infile), status="old", action="read", iostat=iost,err=01)
!!!     read(IIn,'(A)') np%str1
!!!     read(IIn,'(A)') np%str2
!!!     read(IIn,*)     np%natoms, np%xmin,       np%ymin,       np%zmin
!!!     read(IIn,*)     np%nx,     np%dx
!!!     read(IIn,*)     np%ny,     np%dummy_real, np%dy
!!!     read(IIn,*)     np%nz,     np%dummy_real, np%dummy_real, np%dz
!!!!
!!!     allocate(np%atomic_number(np%natoms),np%atomic_charge(np%natoms),&
!!!              np%x(np%natoms),np%y(np%natoms),np%z(np%natoms),np%rho(np%nx,np%ny,np%nz))
!!!!
!!!     np%atomic_number = zero
!!!     np%atomic_charge = zero
!!!     np%x = zero 
!!!     np%y = zero
!!!     np%z = zero
!!!     np%rho = zero
!!!     np%nelectrons = 0
!!!!
!!!     do i=1,np%natoms
!!!        read(IIn,*) np%atomic_number(i),np%atomic_charge(i),np%x(i),np%y(i),np%z(i) 
!!!        np%nelectrons = np%nelectrons + np%atomic_number(i)
!!!     enddo
!!!!
!!!     do i=1,np%nx
!!!        do j=1,np%ny
!!!           read(IIn,*) np%rho(i,j,:)
!!!        enddo
!!!     enddo  
!!!!
!!!     np%rho = np%rho*np%dx*np%dy*np%dz ! weight nanoparticle by np volume
!!!!
!!!     01 continue
!!!     close(IIn)
!!!!
!!!   end subroutine read_nanoparticle
!----------------------------------------------------------------------
!!   subroutine delete_nanoparticle(np)
!!!
!!!    Read input np file
!!!
!!     implicit none
!!!
!!     type (nanoparticle_type), intent(inout)  :: np
!!!
!!     deallocate(np%atomic_number)
!!     deallocate(np%atomic_charge)
!!     deallocate(np%x)
!!     deallocate(np%y)
!!     deallocate(np%z)
!!     deallocate(np%rho)
!!!
!!   end subroutine delete_nanoparticle
!----------------------------------------------------------------------
end module nanoparticle_module


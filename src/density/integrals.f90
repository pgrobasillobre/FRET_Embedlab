!----------------------------------------------------------------------
module integrals_module
!      
!   Module integrals
!
!!    use output_module
!!    use target_module
    use parameters_module
!
    Implicit None
!
    public integrals
!
!   integrals type
!
    type integrals_type
! 
      private
!
      character (len=200)                     :: str1,str2
!
      integer                                 :: natoms,nx,ny,nz
      integer                                 :: nelectrons
      integer, dimension(:), allocatable      :: atomic_number
!
      real(dp)                                :: xmin,ymin,zmin,dx,dy,dz,dummy_real
      real(dp), dimension(:), allocatable     :: atomic_charge,x,y,z
      real(dp), dimension(:,:,:), allocatable :: rho 
!
    end type integrals_type
!
    type (integrals_type), Save :: integrals 
!
   contains
!----------------------------------------------------------------------
!!   subroutine int_density(cube)
!!!
!!!    Integrate density from cube file
!!!
!!     implicit none
!!!
!!     type (density_type), intent(in)  :: cube
!!!
!!!
!!!    internal variables
!!!
!!     real (dp) :: integral
!!     integer   :: i,j,k
!!!
!!     integral = zero
!!!
!!     do i = 1,cube%nx
!!        do j = 1,cube%ny
!!           do k = 1,cube%nz
!!              integral = integral + cube%rho(i,j,k)*cube%dx*cube%dy*cube%dz
!!           enddo
!!        enddo
!!     enddo
!!!
!!    call out_%print_density_integral(cube%natoms,                   &
!!                                     cube%nx,cube%ny,cube%nz,       &
!!                                     cube%dx,cube%dy,cube%dz,       &
!!                                     cube%xmin,cube%ymin,cube%zmin, &
!!                                     cube%nelectrons,integral)
!!!
!!  end subroutine int_density 
!----------------------------------------------------------------------
end module integrals_module


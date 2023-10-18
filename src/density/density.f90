!----------------------------------------------------------------------
module density_module
!      
!   Module density
!
    use output_module
    use target_module
    use parameters_module
!
    Implicit None
!
    public density
!
!   density type
!
    type density_type
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
    end type density_type
!
    type (density_type), Save :: density
!
   contains
!----------------------------------------------------------------------
   subroutine read_density (infile, cube)
!
!    Read input cube file
!
     implicit none
!
     character(len=*),    intent(in)   :: infile
     type (density_type), intent(out)  :: cube
!
!
!    internal variables
!
     integer                         :: IIn
     integer                         :: i,j
     integer                         :: iost
!
     IIn = 21
!
!    Read variables
     open (unit=IIn, file=trim(infile), status="old", action="read", iostat=iost,err=01)
     read(IIn,'(A)') cube%str1
     read(IIn,'(A)') cube%str2
     read(IIn,*)     cube%natoms, cube%xmin,       cube%ymin,       cube%zmin
     read(IIn,*)     cube%nx,     cube%dx
     read(IIn,*)     cube%ny,     cube%dummy_real, cube%dy
     read(IIn,*)     cube%nz,     cube%dummy_real, cube%dummy_real, cube%dz
!
     allocate(cube%atomic_number(cube%natoms),cube%atomic_charge(cube%natoms),&
              cube%x(cube%natoms),cube%y(cube%natoms),cube%z(cube%natoms),cube%rho(cube%nx,cube%ny,cube%nz))
!
     cube%atomic_number = zero
     cube%atomic_charge = zero
     cube%x = zero 
     cube%y = zero
     cube%z = zero
     cube%rho = zero
     cube%nelectrons = 0
!
     do i=1,cube%natoms
        read(IIn,*) cube%atomic_number(i),cube%atomic_charge(i),cube%x(i),cube%y(i),cube%z(i) 
        cube%nelectrons = cube%nelectrons + cube%atomic_number(i)
     enddo
!
     do i=1,cube%nx
        do j=1,cube%ny
           read(IIn,*) cube%rho(i,j,:)
        enddo
     enddo  
!
     01 continue
     close(IIn)
!
   end subroutine read_density
!----------------------------------------------------------------------
   subroutine int_density(cube)
!
!    Integrate density from cube file
!
     implicit none
!
     type (density_type), intent(in)  :: cube
!
!
!    internal variables
!
     real (dp) :: integral
     integer   :: i,j,k
!
     integral = zero
!
     do i = 1,cube%nx
        do j = 1,cube%ny
           do k = 1,cube%nz
              integral = integral + cube%rho(i,j,k)*cube%dx*cube%dy*cube%dz
           enddo
        enddo
     enddo
!
    call out_%print_density_integral(cube%natoms,                   &
                                     cube%nx,cube%ny,cube%nz,       &
                                     cube%dx,cube%dy,cube%dz,       &
                                     cube%xmin,cube%ymin,cube%zmin, &
                                     cube%nelectrons,integral)
!
  end subroutine int_density 
!----------------------------------------------------------------------
end module density_module


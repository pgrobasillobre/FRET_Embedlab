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
      character (len=200)                       :: str1,str2
!
      integer                                   :: natoms,nx,ny,nz
      integer                                   :: n_points_reduced
      integer                                   :: nelectrons
      integer,   dimension(:), allocatable      :: atomic_number
      character, dimension(:), allocatable      :: atomic_label
!
      real(dp)                                  :: xmin,ymin,zmin,dx,dy,dz,dummy_real
      real(dp)                                  :: maxdens
      real(dp), dimension(:), allocatable       :: atomic_charge,x,y,z
      real(dp), dimension(:,:,:), allocatable   :: rho
      real(dp), dimension(:), allocatable       :: rho_reduced
      real(dp), dimension(:,:), allocatable     :: xyz 
!
    end type density_type
!
    type (density_type), Save :: density
!
   contains
!----------------------------------------------------------------------
   subroutine read_density(infile, cube)
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
     integer                         :: i,j,k
     integer                         :: iost
!
     real(dp)                        :: x_tmp, y_tmp, z_tmp
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
     allocate(cube%atomic_number(cube%natoms),cube%atomic_label(cube%natoms),cube%atomic_charge(cube%natoms),&
              cube%x(cube%natoms),cube%y(cube%natoms),cube%z(cube%natoms),cube%rho(cube%nx,cube%ny,cube%nz), &
              cube%xyz(3,ncellmax), cube%rho_reduced(ncellmax))
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
        cube%atomic_label(i) = map_atomic_number_to_label(cube%atomic_number(i))
     enddo
!
     do i=1,cube%nx
        do j=1,cube%ny
           read(IIn,*) cube%rho(i,j,:)
        enddo
     enddo  
!
     cube%rho = cube%rho*cube%dx*cube%dy*cube%dz ! weight density by cube volume
!
!    Find the maximum value of the density
! 
     cube%maxdens = zero
     do i=1,cube%nx
        do j=1,cube%ny
           do k=1,cube%nz
              if (abs(cube%rho(i,j,k)) .gt. cube%maxdens) cube%maxdens = abs(cube%rho(i,j,k))
           enddo
        enddo
     enddo  
!
!    Save reduce density of the cube, and calculate associated coordinates
!
     cube%xyz = zero
!
     cube%n_points_reduced = 0
     do i = 1, cube%nx
        x_tmp = cube%xmin + cube%dx*(i-1)
        do j = 1, cube%ny
           y_tmp = cube%ymin + cube%dy*(j-1)
           do k = 1, cube%nz
              z_tmp = cube%zmin + cube%dz*(k-1)
!
!             If we calculate the overlap integral we should not reduce the density cube
              if (.not. target_%calc_overlap_int .and.  cube%n_points_reduced .gt. ncellmax) &
                                     call out_%error("Increase cutoff, huge density to be managed")
!
              if (.not. target_%calc_overlap_int .and. abs(cube%rho(i,j,k)) .gt. cube%maxdens * target_%cutoff .or.&
                  target_%calc_overlap_int) then
!
                 cube%n_points_reduced = cube%n_points_reduced + 1
!
                 cube%rho_reduced(cube%n_points_reduced) = cube%rho(i,j,k)
!
                 cube%xyz(1,cube%n_points_reduced) = x_tmp 
                 cube%xyz(2,cube%n_points_reduced) = y_tmp
                 cube%xyz(3,cube%n_points_reduced) = z_tmp
!
              endif
!
           enddo
        enddo
     enddo
!
     01 continue
     close(IIn)
!
   end subroutine read_density
!----------------------------------------------------------------------
   subroutine delete_density(cube)
!
!    Read input cube file
!
     implicit none
!
     type (density_type), intent(inout)  :: cube
!
     deallocate(cube%atomic_number)
     deallocate(cube%atomic_charge)
     deallocate(cube%x)
     deallocate(cube%y)
     deallocate(cube%z)
     deallocate(cube%rho)
!
   end subroutine delete_density
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
              integral = integral + cube%rho(i,j,k)!*cube%dx*cube%dy*cube%dz
           enddo
        enddo
     enddo
!
    call out_%print_density(target_%density_file,cube%natoms,cube%n_points_reduced, &
                            cube%nx,cube%ny,cube%nz,                                &
                            cube%dx,cube%dy,cube%dz,                                &
                            cube%xmin,cube%ymin,cube%zmin,cube%nelectrons,          &
                            cube%atomic_label,cube%x,cube%y,cube%z,                 &
                            integral=integral)
!
!!    call out_%print_density(cube)
!
  end subroutine int_density 
!----------------------------------------------------------------------
  function map_atomic_number_to_label(atnum) result(label)
!
!   Map atomic number to atomic label
!
    implicit none
!
    integer          :: atnum
    character(len=2) :: label
!
    if (atnum.eq.1) then
       label = 'H'
    else if (atnum.eq.2) then 
       label = 'He'
    else if (atnum.eq.3) then
       label = 'Li'
    else if (atnum.eq.4) then
       label = 'Be'
    else if (atnum.eq.5) then
       label = 'B'
    else if (atnum.eq.6) then
       label = 'C'
    else if (atnum.eq.7) then
       label = 'N'
    else if (atnum.eq.8) then
       label = 'O'
    else if (atnum.eq.9) then
       label = 'F'
    else if (atnum.eq.10) then
       label = 'Ne'
    else if (atnum.eq.11) then
       label = 'Na'
    else if (atnum.eq.12) then
       label = 'Mg'
    else if (atnum.eq.13) then
       label = 'Al'
    else if (atnum.eq.14) then
       label = 'Si'
    else if (atnum.eq.15) then
       label = 'P'
    else if (atnum.eq.16) then
       label = 'S'
    else if (atnum.eq.17) then
       label = 'Cl'
    else if (atnum.eq.18) then
       label = 'Ar'
    else
       call out_%error('  We only consider the mapping of atomic numbers until Ar')
    endif
!
  end function
!----------------------------------------------------------------------
end module density_module


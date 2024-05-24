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
   subroutine read_density(infile, cube, rotation)
!
!    Read input cube file
!
     implicit none
!
     character(len=*),    intent(in)  :: infile
!
     logical,             intent(in)  :: rotation
!
     type (density_type), intent(out) :: cube
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
     if (rotation) call rotate_cube_coordinates(infile,cube)
     call print_cube_coordinates('aceptor_cube_points_ROTATED',cube%n_points_reduced,cube%xyz(1:3,1:cube%n_points_reduced)) 
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
   subroutine rotate_cube_coordinates(infile,cube)
!
!    Rotate cube coordinates
!
     implicit none
!
     character(len=*), intent(in) :: infile
!
     type (density_type), intent(inout) :: cube
!
!    internal variables
!
     integer :: i
!
     real(dp) :: cos_theta, sin_theta
!
     real(dp), dimension(3) :: geom_center
!
     real(dp), dimension(3,cube%n_points_reduced) :: xyz_000, xyz_rot
!
print *, ' Rotation for cube file ', trim(adjustl(infile))
print *, ' Rotation angle ' , target_%aceptor_density_rotation_angle
print *, ' Transdip ', target_%aceptor_transdip(1:3) 
!
     if (target_%debug) call print_cube_coordinates('aceptor_cube_points',cube%n_points_reduced,cube%xyz)
!
     geom_center(1) = sum(cube%xyz(1,1:cube%n_points_reduced))/cube%n_points_reduced 
     geom_center(2) = sum(cube%xyz(2,1:cube%n_points_reduced))/cube%n_points_reduced
     geom_center(3) = sum(cube%xyz(3,1:cube%n_points_reduced))/cube%n_points_reduced
!
!    Translate density to the origin of coordinates
!
     do i = 1,cube%n_points_reduced
        xyz_000(1,i) = cube%xyz(1,i) - geom_center(1)
        xyz_000(2,i) = cube%xyz(2,i) - geom_center(2)
        xyz_000(3,i) = cube%xyz(3,i) - geom_center(3)
     enddo
!
!    Rotate density
!
     cos_theta = dcos(target_%aceptor_density_rotation_angle)
     sin_theta = dsin(target_%aceptor_density_rotation_angle) 
!
     if (target_%rotation_axys.eq.'x') then
!
        do i = 1,cube%n_points_reduced
           xyz_rot(1,i) = xyz_000(1,i)
           xyz_rot(2,i) = xyz_000(2,i) * cos_theta - xyz_000(3,i) * sin_theta 
           xyz_rot(3,i) = xyz_000(2,i) * sin_theta + xyz_000(3,i) * cos_theta
        enddo

     else if (target_%rotation_axys.eq.'y') then
!
        do i = 1,cube%n_points_reduced
           xyz_rot(1,i) =  xyz_000(1,i) * cos_theta + xyz_000(3,i) * sin_theta
           xyz_rot(2,i) =  xyz_000(2,i)  
           xyz_rot(3,i) = -xyz_000(1,i) * sin_theta + xyz_000(3,i) * cos_theta
        enddo
!
     else if (target_%rotation_axys.eq.'z') then
!
        do i = 1,cube%n_points_reduced
           xyz_rot(1,i) = xyz_000(1,i) * cos_theta - xyz_000(2,i) * sin_theta
           xyz_rot(2,i) = xyz_000(1,i) * sin_theta + xyz_000(2,i) * cos_theta  
           xyz_rot(3,i) = xyz_000(3,i)
        enddo
!
     endif
!
!    Translate density to initial position
!
     do i = 1,cube%n_points_reduced
        cube%xyz(1,i) = xyz_rot(1,i) + geom_center(1)
        cube%xyz(2,i) = xyz_rot(2,i) + geom_center(2)
        cube%xyz(3,i) = xyz_rot(3,i) + geom_center(3)
     enddo
!
!     call print_cube_coordinates('aceptor_cube_points_ROTATED',cube%n_points_reduced,cube%xyz(1:3,1:cube%n_points_reduced)) 
!
! por algunha razon as coordeadas rotadas non saen da subrutina
print *, 'por algunha razon as coordeadas rotadas non saen da subrutina'
stop
!
!
  end subroutine rotate_cube_coordinates
!----------------------------------------------------------------------
   subroutine print_cube_coordinates(infile,n_points,xyz)
!
!    Print  cube coordinates
!
     implicit none
!
     character(len=*), intent(in) :: infile
!
     integer, intent(in)          :: n_points
!
     real(dp), dimension(3,n_points), intent(inout) :: xyz
!
!    internal variables
!
     integer :: i
     integer :: iin = 21
!
     open(unit=iin,file=trim(infile)//'.xyz',status="unknown")
!
        write(iin,*) n_points
        write(iin,*) ' Original cube coordinates'
!
        do i = 1, n_points
           write(iin, '(a,f25.16,2x,f25.16,2x,f25.16)' ) 'H' ,xyz(1,i), xyz(2,i), xyz(3,i) 
        enddo
!
     close(iin)
!
  end subroutine print_cube_coordinates
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
    else if (atnum.eq.19) then
       label = 'K'
    else if (atnum.eq.20) then
       label = 'Ca'
    else if (atnum.eq.21) then
       label = 'Sc'
    else if (atnum.eq.22) then
       label = 'Ti'
    else if (atnum.eq.23) then
       label = 'V'
    else if (atnum.eq.24) then
       label = 'Cr'
    else if (atnum.eq.25) then
       label = 'Mn'
    else if (atnum.eq.26) then
       label = 'Fe'
    else if (atnum.eq.27) then
       label = 'Co'
    else if (atnum.eq.28) then
       label = 'Ni'
    else if (atnum.eq.29) then
       label = 'Cu'
    else if (atnum.eq.30) then
       label = 'Zn'
    else if (atnum.eq.31) then
       label = 'Ga'
    else if (atnum.eq.32) then
       label = 'Ge'
    else if (atnum.eq.33) then
       label = 'As'
    else if (atnum.eq.34) then
       label = 'Se'
    else if (atnum.eq.35) then
       label = 'Br'
    else if (atnum.eq.36) then
       label = 'Kr'
    else if (atnum.eq.37) then
       label = 'Rb'
    else if (atnum.eq.38) then
       label = 'Sr'
    else if (atnum.eq.39) then
       label = 'Y'
    else if (atnum.eq.40) then
       label = 'Zr'
    else if (atnum.eq.41) then
       label = 'Nb'
    else if (atnum.eq.42) then
       label = 'Mo'
    else if (atnum.eq.43) then
       label = 'Tc'
    else if (atnum.eq.44) then
       label = 'Ru'
    else if (atnum.eq.45) then
       label = 'Rh'
    else if (atnum.eq.46) then
       label = 'Pd'
    else if (atnum.eq.47) then
       label = 'Ag'
    else if (atnum.eq.48) then
       label = 'Cd'
    else if (atnum.eq.49) then
       label = 'In'
    else if (atnum.eq.50) then
       label = 'Sn'
    else if (atnum.eq.51) then
       label = 'Sb'
    else if (atnum.eq.52) then
       label = 'Te'
    else if (atnum.eq.53) then
       label = 'I'
    else if (atnum.eq.54) then
       label = 'Xe'
    else if (atnum.eq.55) then
       label = 'Cs'
    else if (atnum.eq.56) then
       label = 'Ba'
    else if (atnum.eq.57) then
       label = 'La'
    else if (atnum.eq.58) then
       label = 'Ce'
    else if (atnum.eq.59) then
       label = 'Pr'
    else if (atnum.eq.60) then
       label = 'Nd'
    else if (atnum.eq.61) then
       label = 'Pm'
    else if (atnum.eq.62) then
       label = 'Sm'
    else if (atnum.eq.63) then
       label = 'Eu'
    else if (atnum.eq.64) then
       label = 'Gd'
    else if (atnum.eq.65) then
       label = 'Tb'
    else if (atnum.eq.66) then
       label = 'Dy'
    else if (atnum.eq.67) then
       label = 'Ho'
    else if (atnum.eq.68) then
       label = 'Er'
    else if (atnum.eq.69) then
       label = 'Tm'
    else if (atnum.eq.70) then
       label = 'Yb'
    else if (atnum.eq.71) then
       label = 'Lu'
    else if (atnum.eq.72) then
       label = 'Hf'
    else if (atnum.eq.73) then
       label = 'Ta'
    else if (atnum.eq.74) then
       label = 'W'
    else if (atnum.eq.75) then
       label = 'Re'
    else if (atnum.eq.76) then
       label = 'Os'
    else if (atnum.eq.77) then
       label = 'Ir'
    else if (atnum.eq.78) then
       label = 'Pt'
    else if (atnum.eq.79) then
       label = 'Au'
    else if (atnum.eq.80) then
       label = 'Hg'
    else
       call out_%error('  We only consider the mapping of atomic numbers until Ar')
    endif
!
  end function
!----------------------------------------------------------------------
end module density_module


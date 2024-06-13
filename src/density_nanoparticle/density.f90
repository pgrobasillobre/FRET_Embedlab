!----------------------------------------------------------------------
module density_module
!      
!   Module density
!
    !$ use omp_lib
    use target_module
    use parameters_module
!
    Implicit None
!
    public rotate_cube
    public rotate_transdip
    public density
!
!   density type
!
    type density_type
! 
      character (len=200)                         :: str1,str2
!
      integer                                     :: natoms,nx,ny,nz
      integer                                     :: n_points_reduced
      integer                                     :: nelectrons
      integer,   dimension(:), allocatable        :: atomic_number
      character(len=2), dimension(:), allocatable :: atomic_label
!
      real(dp)                                    :: xmin,ymin,zmin
      real(dp)                                    :: maxdens, volume
      real(dp), dimension(3)                      :: geom_center, geom_center_mol, dx, dy, dz
      real(dp), dimension(:), allocatable         :: atomic_charge,x,y,z
      real(dp), dimension(:,:,:), allocatable     :: rho
      real(dp), dimension(:), allocatable         :: rho_reduced
      real(dp), dimension(:,:), allocatable       :: xyz 
!
    end type density_type
!
    type (density_type), Save :: density
!
    Interface rotate_cube
       Module Procedure rotate_cube_coordinates
    End Interface
!
    Interface rotate_transdip
       Module Procedure rotate_transdip_coordinates
    End Interface
!
    Interface delete_cube_density
       Module Procedure delete_density
    End Interface
!
   contains
!----------------------------------------------------------------------
   subroutine read_density(infile, cube, rotation, what_dens)
!
!    Read input cube file
!
     implicit none
!
     character(len=*),    intent(in)        :: infile
!
     logical,             intent(in)        :: rotation
!
     character(len=*), optional, intent(in) :: what_dens
!
     type (density_type), intent(out) :: cube
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
     read(IIn,*)     cube%nx,     cube%dx(1),      cube%dx(2),      cube%dx(3)
     read(IIn,*)     cube%ny,     cube%dy(1),      cube%dy(2),      cube%dy(3)
     read(IIn,*)     cube%nz,     cube%dz(1),      cube%dz(2),      cube%dz(3)
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
!    Check we have an appropriate cube file to weight the density by cube volume
!
     if (cube%dx(2) .gt. 0.0 .or. cube%dx(3) .gt. 0.0 .or. &
         cube%dy(1) .gt. 0.0 .or. cube%dy(3) .gt. 0.0 .or. &
         cube%dz(1) .gt. 0.0 .or. cube%dz(2) .gt. 0.0) call error("Cube file conflict: dx,dy,dz matrix is not diagonal")
!
     cube%volume = cube%dx(1)*cube%dy(2)*cube%dz(3)
     cube%rho = cube%rho*cube%dx(1)*cube%dy(2)*cube%dz(3) ! weight density by cube volume
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
        x_tmp = cube%xmin + cube%dx(1)*(i-1)
        do j = 1, cube%ny
           y_tmp = cube%ymin + cube%dy(2)*(j-1)
           do k = 1, cube%nz
              z_tmp = cube%zmin + cube%dz(3)*(k-1)
!
!             If we calculate the overlap integral we should not reduce the density cube
              if (.not. target_%calc_overlap_int .and.  cube%n_points_reduced .gt. ncellmax) &
                                     call error("Increase cutoff, huge density to be managed")
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
!    Calculate geometrical center of cube and molecule
!
     cube%geom_center(1) = sum(cube%xyz(1,1:cube%n_points_reduced))/cube%n_points_reduced 
     cube%geom_center(2) = sum(cube%xyz(2,1:cube%n_points_reduced))/cube%n_points_reduced
     cube%geom_center(3) = sum(cube%xyz(3,1:cube%n_points_reduced))/cube%n_points_reduced
!
     cube%geom_center_mol(1) = sum(cube%x(1:cube%natoms))/cube%natoms 
     cube%geom_center_mol(2) = sum(cube%y(1:cube%natoms))/cube%natoms
     cube%geom_center_mol(3) = sum(cube%z(1:cube%natoms))/cube%natoms
!
     if (abs(cube%geom_center_mol(1) - cube%geom_center(1)) > 0.3 .or. &
         abs(cube%geom_center_mol(2) - cube%geom_center(2)) > 0.3 .or. &
         abs(cube%geom_center_mol(3) - cube%geom_center(3)) > 0.3) &
         call error("Cube file corrupt: cube center and molecular center do not coincide")
!
!    Rotate aceptor and/or donor cube if requested 
!
     if (rotation .and. what_dens.eq.'aceptor') then 
!
        if (target_%aceptor_transdip_rotate) call rotate_transdip_coordinates(what_dens,                                 &
                                                                              cube%geom_center_mol,                      &
                                                                              target_%aceptor_ref_vector,                &
                                                                              target_%aceptor_transdip_rotate_align_with,&
                                                                              target_%aceptor_transdip,                  &
                                                                              target_%aceptor_transdip_rot,              &
                                                                              target_%aceptor_density_rotation_angle) 
!
        call rotate_cube_coordinates(what_dens,target_%aceptor_density_rotation_angle,cube)
!
     else if (rotation .and. what_dens.eq.'donor') then 
!
        if (target_%donor_transdip_rotate) call rotate_transdip_coordinates(what_dens,                               &
                                                                            cube%geom_center_mol,                    &
                                                                            target_%donor_ref_vector,                &
                                                                            target_%donor_transdip_rotate_align_with,&
                                                                            target_%donor_transdip,                  &
                                                                            target_%donor_transdip_rot,              &
                                                                            target_%donor_density_rotation_angle) 
!
        call rotate_cube_coordinates(what_dens,target_%donor_density_rotation_angle,cube)
!
     endif
!
     if (target_%debug.ge.2) call print_cube_density(outfile='debug/'//what_dens//'_fret.cube',scale_volume=.true.,cube=cube)
!
     01 continue
     close(IIn)
!
   end subroutine read_density
!----------------------------------------------------------------------
   subroutine print_cube_density(outfile, scale_volume, cube)
!
!    Read input cube file
!
     implicit none
!
     character(len=*),    intent(in) :: outfile
!
     logical, intent(in)             :: scale_volume
!
     type (density_type), intent(in) :: cube
!
!
!    internal variables
!
     integer                         :: IIn = 21
     integer                         :: i,j,k
     integer                         :: iost
!
     1000 Format(I5,1x,F11.6,1x,F11.6,1x,F11.6)
     1001 Format(I5,1x,F11.6,1x,F11.6,1x,F11.6,1x,F11.6)
!
!    Read variables
     open (unit=IIn, file=trim(outfile), status="new", action="write", iostat=iost,err=01)
     write(IIn,'(A)') '*** adf ***'
     write(IIn,'(A)') 'FRET_Embedlab printed density'
     write(IIn,1000)     cube%natoms, cube%xmin,       cube%ymin,       cube%zmin
     write(IIn,1000)     cube%nx,     cube%dx(1),      cube%dx(2),      cube%dx(3)
     write(IIn,1000)     cube%ny,     cube%dy(1),      cube%dy(2),      cube%dy(3)
     write(IIn,1000)     cube%nz,     cube%dz(1),      cube%dz(2),      cube%dz(3)
!
!
     do i=1,cube%natoms
        write(IIn,1001) cube%atomic_number(i),cube%atomic_charge(i),cube%x(i),cube%y(i),cube%z(i) 
     enddo
!
     do i=1,cube%nx
        do j=1,cube%ny
          do k = 1,cube%nz
              if (scale_volume) write(IIn,'(6E13.5)', advance='no') cube%rho(i,j,k)/cube%volume
              if (.not.scale_volume) write(IIn,'(6E13.5)', advance='no') cube%rho(i,j,k)
              if (mod(k,6) == 0) write(IIn,*) ! New line every 6 numbers              
           enddo
           if (mod(k,6).ne.0) write(IIn,*) ! New line if not already done in previous iteration
        enddo
     enddo  
!
     01 continue
     close(IIn)
!
   end subroutine print_cube_density
!----------------------------------------------------------------------
   subroutine delete_density(cube)
!
!    Read input cube file
!
     implicit none
!
     type (density_type), intent(inout)  :: cube
!
     if (allocated(cube%atomic_number)) deallocate(cube%atomic_number)
     if (allocated(cube%atomic_label) ) deallocate(cube%atomic_label)
     if (allocated(cube%atomic_charge)) deallocate(cube%atomic_charge)
     if (allocated(cube%x)            ) deallocate(cube%x)
     if (allocated(cube%y)            ) deallocate(cube%y)
     if (allocated(cube%z)            ) deallocate(cube%z)
     if (allocated(cube%rho)          ) deallocate(cube%rho)
     if (allocated(cube%xyz)          ) deallocate(cube%xyz)
     if (allocated(cube%rho_reduced)  ) deallocate(cube%rho_reduced)
!
   end subroutine delete_density
!----------------------------------------------------------------------
   subroutine int_density(cube,integral)
!
!    Integrate density from cube file
!
     implicit none
!
     type (density_type), intent(in)  :: cube
!
     real(dp), intent(out)            :: integral
!
!
!    internal variables
!
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
  end subroutine int_density 
!----------------------------------------------------------------------
   subroutine rotate_cube_coordinates(what_dens,rot_angle,cube)
!
!    Rotate cube coordinates
!
     implicit none
!
     character(len=*), intent(in) :: what_dens
!
     real(dp), intent(in)         :: rot_angle
!
     type (density_type), intent(inout) :: cube
!
!    internal variables
!
     integer :: i, j, k, l
!
     real(dp) :: cos_theta, sin_theta
!
     real(dp)                                     :: x_tmp, y_tmp, z_tmp, oversize
     real(dp), dimension(3,cube%n_points_reduced) :: xyz_rot
!
     type (density_type) :: cube_rot
!
     if (target_%debug.ge.1) call print_cube_coordinates('debug/'//what_dens//'_cube_points',cube%n_points_reduced,cube%xyz)
!
!    Calculate and save angle cosine and sine
!
     cos_theta = dcos(rot_angle)
     sin_theta = dsin(rot_angle) 
!
!
!    ================ MODIFY DENSITY ================= 
!
!    Translate density center to the origin of coordinates
!
     do i = 1,cube%n_points_reduced
        cube%xyz(1,i) = cube%xyz(1,i) - cube%geom_center(1)
        cube%xyz(2,i) = cube%xyz(2,i) - cube%geom_center(2)
        cube%xyz(3,i) = cube%xyz(3,i) - cube%geom_center(3)
     enddo
!
!    Rotate
!
     if (target_%rotation_axys.eq.'x') then
!
        do i = 1,cube%n_points_reduced
           xyz_rot(1,i) = cube%xyz(1,i)
           xyz_rot(2,i) = cube%xyz(2,i) * cos_theta - cube%xyz(3,i) * sin_theta 
           xyz_rot(3,i) = cube%xyz(2,i) * sin_theta + cube%xyz(3,i) * cos_theta
        enddo

     else if (target_%rotation_axys.eq.'y') then
!
        do i = 1,cube%n_points_reduced
           xyz_rot(1,i) =  cube%xyz(1,i) * cos_theta + cube%xyz(3,i) * sin_theta
           xyz_rot(2,i) =  cube%xyz(2,i)  
           xyz_rot(3,i) = -cube%xyz(1,i) * sin_theta + cube%xyz(3,i) * cos_theta
        enddo
!
     else if (target_%rotation_axys.eq.'z') then
!
        do i = 1,cube%n_points_reduced
           xyz_rot(1,i) = cube%xyz(1,i) * cos_theta - cube%xyz(2,i) * sin_theta
           xyz_rot(2,i) = cube%xyz(1,i) * sin_theta + cube%xyz(2,i) * cos_theta  
           xyz_rot(3,i) = cube%xyz(3,i)
        enddo
!
     endif
!
!    Translate density to initial position
!
     do i = 1,cube%n_points_reduced
        cube%xyz(1,i) = xyz_rot(1,i) + cube%geom_center(1)
        cube%xyz(2,i) = xyz_rot(2,i) + cube%geom_center(2)
        cube%xyz(3,i) = xyz_rot(3,i) + cube%geom_center(3)
     enddo
!
     if (target_%debug.ge.1) call print_cube_coordinates('debug/'//what_dens//'_cube_points_ROTATED', &
                                                     cube%n_points_reduced,              &
                                                     cube%xyz(1:3,1:cube%n_points_reduced)) 
!
!
!    =============== ROTATED CUBE TO PRINT AS CHECK ================ 
!
     if (target_%debug.ge.2) then
!
!       Create new cube
        cube_rot%geom_center = zero
        cube_rot%geom_center_mol = zero
        cube_rot%dx = zero 
        cube_rot%dy = zero
        cube_rot%dz = zero
!
        cube_rot%natoms = cube%natoms
!
        oversize = 2.0 ! Oversize in angstroms for rotated cube reconstruction
!
!       Establish limits of new cube
        cube_rot%xmin = (minval(cube%xyz(1,1:cube%n_points_reduced)) * ToAng - oversize)*ToBohr 
        cube_rot%ymin = (minval(cube%xyz(2,1:cube%n_points_reduced)) * ToAng - oversize)*ToBohr
        cube_rot%zmin = (minval(cube%xyz(3,1:cube%n_points_reduced)) * ToAng - oversize)*ToBohr
!
        cube_rot%dx(1) = cube%dx(1)*2.0
        cube_rot%dy(2) = cube%dy(2)*2.0
        cube_rot%dz(3) = cube%dz(3)*2.0
!
        cube_rot%volume = cube_rot%dx(1) * cube_rot%dy(2) * cube_rot%dz(3)
!
        cube_rot%nx = int((oversize*two*ToBohr + maxval(cube%xyz(1,1:cube%n_points_reduced)) - &
                                                 minval(cube%xyz(1,1:cube%n_points_reduced)))/cube_rot%dx(1))
        cube_rot%ny = int((oversize*two*ToBohr + maxval(cube%xyz(2,1:cube%n_points_reduced)) - &
                                                 minval(cube%xyz(2,1:cube%n_points_reduced)))/cube_rot%dy(2))
        cube_rot%nz = int((oversize*two*ToBohr + maxval(cube%xyz(3,1:cube%n_points_reduced)) - &
                                                 minval(cube%xyz(3,1:cube%n_points_reduced)))/cube_rot%dz(3))
!
!
!       Fill more cube_rot attributes
        allocate(cube_rot%atomic_number(cube_rot%natoms),             &
                 cube_rot%atomic_label(cube_rot%natoms),              &
                 cube_rot%atomic_charge(cube_rot%natoms),             &
                 cube_rot%x(cube_rot%natoms),                         &
                 cube_rot%y(cube_rot%natoms),                         &
                 cube_rot%z(cube_rot%natoms),                         &
                 cube_rot%rho(cube_rot%nx,cube_rot%ny,cube_rot%nz),   &
                 cube_rot%xyz(3,cube_rot%nx*cube_rot%ny*cube_rot%nz))    
!
        cube_rot%atomic_number(1:cube_rot%natoms) = cube%atomic_number(1:cube%natoms)
        cube_rot%atomic_label(1:cube_rot%natoms)  = cube%atomic_label(1:cube%natoms)
        cube_rot%atomic_charge(1:cube_rot%natoms) = cube%atomic_charge(1:cube%natoms)
        cube_rot%x(1:cube_rot%natoms)             = cube%x(1:cube%natoms) 
        cube_rot%y(1:cube_rot%natoms)             = cube%y(1:cube%natoms)        
        cube_rot%z(1:cube_rot%natoms)             = cube%z(1:cube%natoms)        
        cube_rot%rho = zero
        cube_rot%xyz = zero
!
!       
!       Map rotated reduced cube%rho to new cube_rot%rho
!
        do i = 1, cube_rot%nx
           x_tmp = cube_rot%xmin + cube_rot%dx(1)*(i-1)
           do j = 1, cube_rot%ny
              y_tmp = cube_rot%ymin + cube_rot%dy(2)*(j-1)
              do k = 1, cube_rot%nz
                 z_tmp = cube_rot%zmin + cube_rot%dz(3)*(k-1)
!
                 cube_rot%xyz(1,i) = x_tmp
                 cube_rot%xyz(2,j) = y_tmp
                 cube_rot%xyz(3,k) = z_tmp
!
!                Map new coordinates to rotated coordinates and assign density values.
!                Threshold for mapping is the smallest cube volume in original cube
!
                 do l = 1, cube%n_points_reduced
!
                    if (abs(x_tmp - cube%xyz(1,l)) .lt. cube%dx(1) .and. &
                        abs(y_tmp - cube%xyz(2,l)) .lt. cube%dy(2) .and. &
                        abs(z_tmp - cube%xyz(3,l)) .lt. cube%dz(3)) then
!
!                      Rotated density spanning over new volume
                       cube_rot%rho(i,j,k) = cube_rot%volume * (cube%rho_reduced(l)/cube%volume)
                       go to 10

                    endif
!
                 enddo
!
                 10 continue
!
              enddo
           enddo
        enddo
!
        call print_cube_density(outfile='debug/'//what_dens//'_fret_ROTATED_density_only.cube', &
                                scale_volume=.false.,cube=cube_rot)
!
        call delete_density(cube_rot)

     endif
!
  end subroutine rotate_cube_coordinates
!----------------------------------------------------------------------
   subroutine rotate_transdip_coordinates(what_dens,geom_center_mol,ref_vector,align_dips,transdip,transdip_rot,theta)
!
!    Rotate cube coordinates
!
     implicit none
!
     character(len=*), intent(in)          :: what_dens
     real(dp), dimension(3), intent(in)    :: geom_center_mol,ref_vector
     logical, intent(in)                   :: align_dips
!
     real(dp), dimension(3), intent(inout) :: transdip,transdip_rot
     real(dp), intent(inout)               :: theta
!
!    internal variables
!
     real(dp), dimension(3) :: ref_vector_trans

     real(dp) :: cos_theta, sin_theta
     real(dp) :: theta_check
!
     logical  :: rotation_changed
!
     theta_check  = 0.0
!
     rotation_changed = .false.
!
!    Move reference vector to the origin of coordinates and calculate angle if needed
!
     if (align_dips) then 
        ref_vector_trans(1) = ref_vector(1) - geom_center_mol(1)*ToAng 
        ref_vector_trans(2) = ref_vector(2) - geom_center_mol(2)*ToAng
        ref_vector_trans(3) = ref_vector(3) - geom_center_mol(3)*ToAng
!
        theta = vectors_angle(transdip,ref_vector_trans)
     endif
!
!    Rotate transdip 
!
10   cos_theta = dcos(theta)
     sin_theta = dsin(theta) 
!
     if (target_%rotation_axys.eq.'x') then
!
        transdip_rot(1) = transdip(1)
        transdip_rot(2) = transdip(2) * cos_theta - transdip(3) * sin_theta 
        transdip_rot(3) = transdip(2) * sin_theta + transdip(3) * cos_theta

     else if (target_%rotation_axys.eq.'y') then
!
        transdip_rot(1) =  transdip(1) * cos_theta + transdip(3) * sin_theta
        transdip_rot(2) =  transdip(2)  
        transdip_rot(3) = -transdip(1) * sin_theta + transdip(3) * cos_theta
!
     else if (target_%rotation_axys.eq.'z') then
!
        transdip_rot(1) = transdip(1) * cos_theta - transdip(2) * sin_theta
        transdip_rot(2) = transdip(1) * sin_theta + transdip(2) * cos_theta  
        transdip_rot(3) = transdip(3)
!
     endif
!
!    Check angle between rotated transdip and reference vector
!    If they are not coincident, change rotation direction.
!
!    angle_thres: parameter threshold for angle mismatch
!
     if (align_dips) then
!
        theta_check = vectors_angle(transdip_rot,ref_vector_trans)
!
        if (abs(theta_check) .gt. angle_thresh .and. .not. rotation_changed) then 
           theta = -theta
           rotation_changed = .true.
           GOTO 10
!
        else if (abs(theta_check) .gt. angle_thresh .and. rotation_changed) then 
           write(iuout,'(/1x,a,f8.3)') 'Calculated angle: ', theta_check*to_degrees 
           write(iuout,'(a)') ' '
           write(iuout,'(a,1x,f8.3,2x,f8.3,2x,f8.3)') ' Geom-centered reference vector:', &
               ref_vector_trans(1), ref_vector_trans(2), ref_vector_trans(3)
           call error('Alignment with reference vector was not possible')
        endif
!
        if (what_dens.eq.'aceptor') target_%aceptor_angle_check = theta_check*to_degrees 
        if (what_dens.eq.'donor')   target_%donor_angle_check   = theta_check*to_degrees
!
     endif

!
!    Print rotated transdip
!
     if (target_%debug.ge.1) call print_transdip_nmd('debug/transition_dipole_'//what_dens,transdip,geom_center_mol)
     if (target_%debug.ge.1) call print_transdip_nmd('debug/transition_dipole_'//what_dens//'_ROTATED',transdip_rot,geom_center_mol)
!
  end subroutine rotate_transdip_coordinates
!----------------------------------------------------------------------
   subroutine print_cube_coordinates(infile,n_points,xyz)
!
!    print  cube coordinates
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
        write(iin,*) 'cube coordinates'
!
        do i = 1, n_points
           write(iin, '(a,f25.16,2x,f25.16,2x,f25.16)' ) 'h' ,xyz(1,i)*ToAng, &
                                                              xyz(2,i)*ToAng, &
                                                              xyz(3,i)*ToAng
        enddo
!
     close(iin)
!
  end subroutine print_cube_coordinates
!----------------------------------------------------------------------
   subroutine print_transdip_nmd(infile,transdip_vector,center)
!
!    print  cube coordinates
!
     implicit none
!
     character(len=*), intent(in) :: infile
!
     real(dp), dimension(3), intent(in) :: transdip_vector, center
!
!    internal variables
!
     integer :: iin = 21
!
     open(unit=iin,file=trim(infile)//'.nmd',status="unknown")
!
        write(iin,'(a,1x,f10.5,2x,f10.5,2x,f10.5)') 'coordinates ', center(1)*ToAng, center(2)*ToAng, center(3)*ToAng 
        write(iin, '(a,f25.16,2x,f25.16,2x,f25.16)' ) 'mode 1' , transdip_vector(1), transdip_vector(2), transdip_vector(3)
!
     close(iin)
!
  end subroutine print_transdip_nmd
!----------------------------------------------------------------------
   subroutine error(string)
!
!    ERROR subroutine needed to not use output_module (avoid circular dependency)
!
     implicit none
!     
     character(len=*),intent(in) :: string
!
!    check if the file is opened
!
     write(iuout,'(/1x,a)') "Error during the execution of FRET_Embedlab"
     write(iuout,'(1x,a/)') trim(string)
     flush(iuout)
     stop
!       
   end subroutine error
!-----------------------------------------------------------------------
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
       call error('  We only consider the mapping of atomic numbers until Ar')
    endif
!
  end function
!----------------------------------------------------------------------
function vectors_angle(vec_1, vec_2) result(angle)
!
!    Calculate angle between vectors
!
     implicit none
!
     real(dp), dimension(3), intent(in) :: vec_1
     real(dp), dimension(3), intent(in) :: vec_2
!
     real(dp) :: division, angle
!
     real(dp) :: vec_1_X_vec_2, vec_1_mod, vec_2_mod
!
     vec_1_X_vec_2 = dot_product(vec_1, vec_2)
!
     vec_1_mod = dsqrt(dot_product(vec_1, vec_1))
     vec_2_mod = dsqrt(dot_product(vec_2, vec_2))
!
     division = vec_1_X_vec_2 / (vec_1_mod * vec_2_mod)
!
!    Handle decimal points to avoid floating point errors
     if (division .gt.  1.00000000000000d0) division =  1.0d0
     if (division .lt. -1.00000000000000d0) division = -1.0d0
!
     angle = dacos(division) 
!
end function vectors_angle
!----------------------------------------------------------------------
end module density_module


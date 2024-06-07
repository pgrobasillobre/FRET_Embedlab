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
      real(dp), dimension(3)                   :: geom_center
!
      logical :: charges, dipoles
!
    end type nanoparticle_type
!
    type (nanoparticle_type), Save :: nanoparticle
!
    Interface rotate_np
       Module Procedure rotate_np_coords_dips
    End Interface rotate_np
!
   contains
!----------------------------------------------------------------------
   subroutine read_nanoparticle(infile, np, rotation)
!
!    Read input np file
!
     implicit none
!
     character(len=*),    intent(in)        :: infile
!
     logical, intent(in)                    :: rotation
!
     type (nanoparticle_type), intent(out)  :: np
!
!
!    internal variables
!
     character(len=400) :: line
!
     integer :: nlines
     integer :: num_string_initial, num_string_end, dum
     integer :: IIn
     integer :: i
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
!
       Rewind(IIn)
       call go_to_string(IIn,dipoles_header,np%dipoles,nlines,dum)
!
       if (.not.np%charges .and. .not.np%dipoles) call out_%error("Charges and/or dipoles not found")
!
       call allocate_nanoparticle(np)
!
       Rewind(IIn) 
!
       if(np%charges) then !BEM of wFQ case
          call go_to_string(IIn,charges_header,found_string,nlines,dum)
       elseif(np%dipoles) then ! wFQFMu case
          call go_to_string(IIn,dipoles_header,found_string,nlines,dum)
       endif
!
       do i = 1, np%natoms
!
         read(IIn,'(a)') line
!       
         if(len_trim(line).eq.0) then
            call out_%error("Blank line found in corrupt FRET charges/dipoles file: "//trim(infile))
         endif
!
!        Read XYZ coordinates from input file (in Bohr)
!
         if(np%charges) then
            read(line,'(5(f25.16, 5x))') np%q(i,1), np%q(i,2), np%xyz(1,i), np%xyz(2,i), np%xyz(3,i)
!
         elseif(np%dipoles) then
            read(line,'(11(f25.16, 5x))') np%q(i,1),    np%q(i,2),                 &
                                          np%mu(1,i,1), np%mu(2,i,1), np%mu(3,i,1),&
                                          np%mu(1,i,2), np%mu(2,i,2), np%mu(3,i,2),&
                                          np%xyz(1,i),  np%xyz(2,i),  np%xyz(3,i)  
!
         endif
!
       enddo
!
     01 continue
     close(IIn)
!
!    Rotate np geometry and (eventually) dipoles if requested
!
     if (rotation) then
!
!       Calculate np's geometrical center
!
        np%geom_center(1) = sum(np%xyz(1,1:np%natoms))/np%natoms 
        np%geom_center(2) = sum(np%xyz(2,1:np%natoms))/np%natoms
        np%geom_center(3) = sum(np%xyz(3,1:np%natoms))/np%natoms
!
        call rotate_np_coords_dips(target_%donor_density_rotation_angle,np)
!
     endif
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
     if (np%charges) then 
        allocate(np%q(np%natoms,2))
     elseif (np%dipoles) then
        allocate(np%q(np%natoms,2))
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
   subroutine rotate_np_coords_dips(rot_angle,np)
!
!    Rotate cube coordinates
!
     implicit none
!
     real(dp), intent(in) :: rot_angle
!
     type (nanoparticle_type), intent(inout) :: np
!
!    internal variables
!
     integer :: i, j, k, l
!
     real(dp) :: cos_theta, sin_theta
!
     real(dp)                              :: x_tmp, y_tmp, z_tmp, oversize
     real(dp), allocatable, dimension(:,:) :: xyz_rot
!
     type (nanoparticle_type) :: np_rot
!
     if (target_%debug) call print_np_coords_dips('debug/np',np)
!
!    Calculate and save angle cosine and sine
!
     cos_theta = dcos(rot_angle)
     sin_theta = dsin(rot_angle)
!
!
!    ================ MODIFY NP GEOMETRY ================= 
!
     allocate(xyz_rot(3,np%natoms))
!
!    Translate density center to the origin of coordinates
!
     do i = 1,np%natoms
        np%xyz(1,i) = np%xyz(1,i) - np%geom_center(1)
        np%xyz(2,i) = np%xyz(2,i) - np%geom_center(2)
        np%xyz(3,i) = np%xyz(3,i) - np%geom_center(3)
     enddo
!
!    Rotate
!
     if (target_%rotation_axys.eq.'x') then
!
        do i = 1,np%natoms
           xyz_rot(1,i) = np%xyz(1,i)
           xyz_rot(2,i) = np%xyz(2,i) * cos_theta - np%xyz(3,i) * sin_theta
           xyz_rot(3,i) = np%xyz(2,i) * sin_theta + np%xyz(3,i) * cos_theta
        enddo

     else if (target_%rotation_axys.eq.'y') then
!
        do i = 1,np%natoms
           xyz_rot(1,i) =  np%xyz(1,i) * cos_theta + np%xyz(3,i) * sin_theta
           xyz_rot(2,i) =  np%xyz(2,i)
           xyz_rot(3,i) = -np%xyz(1,i) * sin_theta + np%xyz(3,i) * cos_theta
        enddo
!
     else if (target_%rotation_axys.eq.'z') then
!
        do i = 1,np%natoms
           xyz_rot(1,i) = np%xyz(1,i) * cos_theta - np%xyz(2,i) * sin_theta
           xyz_rot(2,i) = np%xyz(1,i) * sin_theta + np%xyz(2,i) * cos_theta
           xyz_rot(3,i) = np%xyz(3,i)
        enddo
!
     endif
!
!    Translate density to initial position
!
     do i = 1,np%natoms
        np%xyz(1,i) = xyz_rot(1,i) + np%geom_center(1)
        np%xyz(2,i) = xyz_rot(2,i) + np%geom_center(2)
        np%xyz(3,i) = xyz_rot(3,i) + np%geom_center(3)
     enddo
!
!!!!
!!!!
!!!!    =============== ROTATED CUBE TO PRINT AS CHECK ================ 
!!!!
!!!     if (target_%debug) then
!!!!
!!!!       Create new cube
!!!        cube_rot%geom_center = zero
!!!        cube_rot%geom_center_mol = zero
!!!        cube_rot%dx = zero
!!!        cube_rot%dy = zero
!!!        cube_rot%dz = zero
!!!!
!!!        cube_rot%natoms = cube%natoms
!!!!
!!!        oversize = 2.0 ! Oversize in angstroms for rotated cube reconstruction
!!!!
!!!!       Establish limits of new cube
!!!        cube_rot%xmin = (minval(cube%xyz(1,1:cube%n_points_reduced)) * ToAng - oversize)*ToBohr
!!!        cube_rot%ymin = (minval(cube%xyz(2,1:cube%n_points_reduced)) * ToAng - oversize)*ToBohr
!!!        cube_rot%zmin = (minval(cube%xyz(3,1:cube%n_points_reduced)) * ToAng - oversize)*ToBohr
!!!!
!!!        cube_rot%dx(1) = cube%dx(1)*2.0
!!!        cube_rot%dy(2) = cube%dy(2)*2.0
!!!        cube_rot%dz(3) = cube%dz(3)*2.0
!!!!
!!!        cube_rot%volume = cube_rot%dx(1) * cube_rot%dy(2) * cube_rot%dz(3)
!!!!
!!!        cube_rot%nx = int((oversize*two*ToBohr + maxval(cube%xyz(1,1:cube%n_points_reduced)) - &
!!!                                                 minval(cube%xyz(1,1:cube%n_points_reduced)))/cube_rot%dx(1))
!!!        cube_rot%ny = int((oversize*two*ToBohr + maxval(cube%xyz(2,1:cube%n_points_reduced)) - &
!!!                                                 minval(cube%xyz(2,1:cube%n_points_reduced)))/cube_rot%dy(2))
!!!        cube_rot%nz = int((oversize*two*ToBohr + maxval(cube%xyz(3,1:cube%n_points_reduced)) - &
!!!                                                 minval(cube%xyz(3,1:cube%n_points_reduced)))/cube_rot%dz(3))
!!!!
!!!!
!!!!       Fill more cube_rot attributes
!!!        allocate(cube_rot%atomic_number(cube_rot%natoms),             &
!!!                 cube_rot%atomic_label(cube_rot%natoms),              &
!!!                 cube_rot%atomic_charge(cube_rot%natoms),             &
!!!                 cube_rot%x(cube_rot%natoms),                         &
!!!                 cube_rot%y(cube_rot%natoms),                         &
!!!                 cube_rot%z(cube_rot%natoms),                         &
!!!                 cube_rot%rho(cube_rot%nx,cube_rot%ny,cube_rot%nz),   &
!!!                 cube_rot%xyz(3,cube_rot%nx*cube_rot%ny*cube_rot%nz))
!!!!
!!!        cube_rot%atomic_number(1:cube_rot%natoms) = cube%atomic_number(1:cube%natoms)
!!!        cube_rot%atomic_label(1:cube_rot%natoms)  = cube%atomic_label(1:cube%natoms)
!!!        cube_rot%atomic_charge(1:cube_rot%natoms) = cube%atomic_charge(1:cube%natoms)
!!!        cube_rot%x(1:cube_rot%natoms)             = cube%x(1:cube%natoms)
!!!        cube_rot%y(1:cube_rot%natoms)             = cube%y(1:cube%natoms)
!!!        cube_rot%z(1:cube_rot%natoms)             = cube%z(1:cube%natoms)
!!!        cube_rot%rho = zero
!!!        cube_rot%xyz = zero
!!!!
!!!!       
!!!!       Map rotated reduced cube%rho to new cube_rot%rho
!!!!
!!!        do i = 1, cube_rot%nx
!!!           x_tmp = cube_rot%xmin + cube_rot%dx(1)*(i-1)
!!!           do j = 1, cube_rot%ny
!!!              y_tmp = cube_rot%ymin + cube_rot%dy(2)*(j-1)
!!!              do k = 1, cube_rot%nz
!!!                 z_tmp = cube_rot%zmin + cube_rot%dz(3)*(k-1)
!!!!
!!!                 cube_rot%xyz(1,i) = x_tmp
!!!                 cube_rot%xyz(2,j) = y_tmp
!!!                 cube_rot%xyz(3,k) = z_tmp
!!!!
!!!!                Map new coordinates to rotated coordinates and assign density values.
!!!!                Threshold for mapping is the smallest cube volume in original cube
!!!!
!!!                 do l = 1, cube%n_points_reduced
!!!!
!!!                    if (abs(x_tmp - cube%xyz(1,l)) .lt. cube%dx(1) .and. &
!!!                        abs(y_tmp - cube%xyz(2,l)) .lt. cube%dy(2) .and. &
!!!                        abs(z_tmp - cube%xyz(3,l)) .lt. cube%dz(3)) then
!!!!
!!!!                      Rotated density spanning over new volume
!!!                       cube_rot%rho(i,j,k) = cube_rot%volume * (cube%rho_reduced(l)/cube%volume)
!!!                       go to 10
!!!
!!!                    endif
!!!!
!!!                 enddo
!!!!
!!!                 10 continue
!!!!
!!!              enddo
!!!           enddo
!!!        enddo
!!!!
!!!        call print_cube_density(outfile='debug/'//what_dens//'_fret_ROTATED_density_only.cube', &
!!!                                scale_volume=.false.,cube=cube_rot)
!!!!
!!!        call delete_density(cube_rot)
!!!
!!!     endif
!!!!
     if (target_%debug) call print_np_coords_dips('debug/np_rot',np)
!
stop
  end subroutine rotate_np_coords_dips
!----------------------------------------------------------------------
   subroutine print_np_coords_dips(infile,np)
!
!    print  cube coordinates
!
     implicit none
!
     character(len=*), intent(in) :: infile
!
     type (nanoparticle_type), intent(inout) :: np
!
!    internal variables
!
     integer :: i
     integer :: iin = 21
!
     open(unit=iin,file=trim(infile)//'.xyz',status="unknown")
!
        write(iin,*) np%natoms
        write(iin,*) 'NP coordinates'
!
        do i = 1, np%natoms
           write(iin, '(a,f25.16,2x,f25.16,2x,f25.16)' ) 'Xx' , np%xyz(1,i)*ToAng, &
                                                                np%xyz(2,i)*ToAng, &
                                                                np%xyz(3,i)*ToAng
        enddo
!
     close(iin)
!
!
     if (np%dipoles) then
!
!       Real part
!
        open(unit=iin,file=trim(infile)//'_re.nmd',status="unknown")
!
           write(iin, '(a)', advance='no') 'coordinates'
           do i = 1, np%natoms
               write(iin, '(3(f10.5, 2x))', advance='no') np%xyz(1,i)*ToAng, np%xyz(2,i)*ToAng, np%xyz(3,i)*ToAng
           end do
           write(iin,*)
           write(iin, '(a)', advance='no') 'mode 1'
           do i = 1, np%natoms
               write(iin, '(3(f10.5, 2x))', advance='no') np%mu(1:3,i,1)
           end do
!
        close(iin)
!
!       Imaginary part
!
        open(unit=iin,file=trim(infile)//'_im.nmd',status="unknown")
!
           write(iin, '(a)', advance='no') 'coordinates'
           do i = 1, np%natoms
               write(iin, '(3(f10.5, 2x))', advance='no') np%xyz(1,i)*ToAng, np%xyz(2,i)*ToAng, np%xyz(3,i)*ToAng
           end do
           write(iin,*)
           write(iin, '(a)', advance='no') 'mode 1'
           do i = 1, np%natoms
               write(iin, '(3(f10.5, 2x))', advance='no') np%mu(1:3,i,2)
           end do
!
        close(iin)
!
     endif
!
  end subroutine print_np_coords_dips
!----------------------------------------------------------------------
end module nanoparticle_module


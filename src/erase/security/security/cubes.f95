MODULE cubes
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE

  TYPE, PUBLIC :: cube
    PRIVATE
    CHARACTER (LEN=72)                           :: str1,str2
    INTEGER                                      :: natoms,nx,ny,nz
    INTEGER, DIMENSION(:), ALLOCATABLE           :: atomic_number
    REAL(KIND=wp)                                :: xmin,ymin,zmin,dx,dy,dz,dummy_real
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE     :: atomic_charge,x,y,z
    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: rho 
  END TYPE cube

  PUBLIC :: cube_get,   &
            cube_write, &
            cube_add,   &
            cube_sub,   &
            cube_int,   &
            cube_cdz,   &
            cube_del
  !--------------------------------------------------
  !--------------------------------------------------
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cube_add
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cube_sub
  END INTERFACE
  !--------------------------------------------------
  !--------------------------------------------------
  CONTAINS
  !*************************************************
  SUBROUTINE cube_get (mycube, infile)
    !
    ! In/out variables
    !
    CHARACTER(LEN=*), INTENT(IN)               :: infile
    TYPE (cube), INTENT(OUT)                   :: mycube
    !
    ! Attributes of mycube (mycube%str1 = ...)
    !
    CHARACTER (LEN=72)                           :: str1,str2
    INTEGER                                      :: natoms,nx,ny,nz
    INTEGER, DIMENSION(:), ALLOCATABLE           :: atomic_number
    REAL(KIND=wp)                                :: xmin,ymin,zmin,dx,dy,dz,dummy_real
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE     :: atomic_charge,x,y,z
    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: rho
    !
    ! Other internal variables
    !
    INTEGER                                    :: IIn
    INTEGER                                    :: i,j
    !
    IIn = 10
    !
    ! Read variables
    OPEN (UNIT=IIn, FILE=infile, STATUS="old", ACTION="read")
    read(IIn,'(A29)') str1
    read(IIn,'(A40)') str2
    read(IIn,*) natoms,xmin,ymin,zmin
    read(IIn,*) nx, dx
    read(IIn,*) ny, dummy_real, dy
    read(IIn,*) nz, dummy_real, dummy_real, dz
    !
    allocate(atomic_number(natoms),atomic_charge(natoms),&
             x(natoms),y(natoms),z(natoms),rho(nx,ny,nz))
    !
    do i=1,natoms
       read(IIn,*) atomic_number(i),atomic_charge(i),x(i),y(i),z(i) 
    enddo
    !
    do i=1,nx
       do j=1,ny
          read(IIn,*) rho(i,j,:)
       enddo
    enddo  
    !
    close(IIn)
    !
    ! Assign attributes to mycube
    mycube%str1          = str1
    mycube%str2          = str2
    mycube%natoms        = natoms
    mycube%nx            = nx
    mycube%ny            = ny
    mycube%nz            = nz
    mycube%atomic_number = atomic_number
    mycube%xmin          = xmin
    mycube%ymin          = ymin
    mycube%zmin          = zmin
    mycube%dx            = dx
    mycube%dy            = dy
    mycube%dz            = dz
    mycube%atomic_charge = atomic_charge
    mycube%x             = x
    mycube%y             = y
    mycube%z             = z
    mycube%rho           = rho 
  END SUBROUTINE cube_get
  !--------------------------------------------------
  !--------------------------------------------------
  SUBROUTINE cube_write (mycube, outfile)
    !
    ! In/out variables
    !
    CHARACTER(LEN=*), INTENT(IN)               :: outfile
    TYPE (cube), INTENT(IN)                    :: mycube
    !
    ! Other internal variables
    !
    INTEGER                                    :: OOut
    INTEGER                                    :: i,j
    REAL(KIND=wp)                              :: dummy_real
    
    !
    OOut = 11
    dummy_real = 0.0d0
    !
    OPEN (UNIT=OOut, FILE=outfile, STATUS="unknown", ACTION="write")
    write(Oout,'(1x,a)') mycube%str1          
    write(Oout,'(1x,a)') mycube%str2          
    write(Oout,'(I5,1x,E15.7,1x,E15.7,1x,E15.7)') mycube%natoms,mycube%xmin,mycube%ymin,mycube%zmin            
    write(Oout,'(I5,1x,E15.7,1x,E15.7,1x,E15.7)') mycube%nx, mycube%dx, dummy_real, dummy_real
    write(Oout,'(I5,1x,E15.7,1x,E15.7,1x,E15.7)') mycube%ny, dummy_real, mycube%dy, dummy_real
    write(Oout,'(I5,1x,E15.7,1x,E15.7,1x,E15.7)') mycube%nz, dummy_real, dummy_real, mycube%dz
    !
    do i =1,mycube%natoms
       write(OOut,'(I5,1x,E15.7,1x,E15.7,1x,E15.7,1x,E15.7)') mycube%atomic_number(i),mycube%atomic_charge(i),&
                                                              mycube%x(i),mycube%y(i),mycube%z(i)
    enddo
    !
    do i=1,mycube%nx
       do j=1,mycube%ny
          write(OOut,'(E15.7)') mycube%rho(i,j,:)
       enddo
    enddo
    close(11)
  END SUBROUTINE cube_write
  !--------------------------------------------------
  !--------------------------------------------------
  FUNCTION cube_add (mycube1, mycube2)
    TYPE(cube) :: cube_add
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    !
    ! Internal variables
    INTEGER                :: i,j,k,counter
    INTEGER                :: total_atoms 
    !
    LOGICAL                :: unique_atom
    !
    cube_add%rho  = mycube1%rho + mycube2%rho
    !
    ! Assign remaining attributes to cube_add
    cube_add%str1          = mycube1%str1
    cube_add%str2          = mycube1%str2
    !
    cube_add%nx            = mycube1%nx
    cube_add%ny            = mycube1%ny
    cube_add%nz            = mycube1%nz
    !
    cube_add%xmin          = mycube1%xmin
    cube_add%ymin          = mycube1%ymin
    cube_add%zmin          = mycube1%zmin
    cube_add%dx            = mycube1%dx
    cube_add%dy            = mycube1%dy
    cube_add%dz            = mycube1%dz
    !
    ! Fill other quantities
    !
    !How many non-repeated atoms do we have?
    total_atoms = mycube1%natoms + mycube2%natoms
    do i = 1,mycube1%natoms
       do j = 1,mycube2%natoms
          if ((mycube1%x(i).eq. mycube2%x(j)) .and. (mycube1%y(i) .eq. mycube2%y(j)) &
                                              .and. (mycube1%z(i) .eq. mycube2%z(j))) then 
             ! If all the positions are the same one atom is repeated
             total_atoms = total_atoms - 1   
          endif
       enddo
    enddo
    !
    cube_add%natoms = total_atoms 
    allocate(cube_add%x(cube_add%natoms),cube_add%y(cube_add%natoms), &
             cube_add%z(cube_add%natoms),cube_add%atomic_number(cube_add%natoms), &
             cube_add%atomic_charge(cube_add%natoms))
   !
   ! Assign non-repeated atomic numbers, charges and positions
    ! From mycube1: assign all atoms
    k = 0
    do i = 1,mycube1%natoms
       k = k + 1
       cube_add%atomic_number(k) = mycube1%atomic_number(i)
       cube_add%atomic_charge(k) = mycube1%atomic_charge(i)
       cube_add%x(k)             = mycube1%x(i)
       cube_add%y(k)             = mycube1%y(i)
       cube_add%z(k)             = mycube1%z(i)            
    enddo
    !
    ! From mycube2: assign non-repeated atoms
    do i = 1,mycube2%natoms
       unique_atom = .false.
       counter = 0
       do j = 1,mycube1%natoms
          if ((mycube2%x(i).ne. mycube1%x(j)) .or. (mycube2%y(i) .ne. mycube1%y(j)) &
                                              .or. (mycube2%z(i) .ne. mycube1%z(j))) then
             ! If one of the positions are different, the atom is unique
             counter = counter + 1
                if (counter .eq. mycube1%natoms) unique_atom = .true.
          endif
       enddo 
       ! If the atom is unique we take all its descriptors
       if (unique_atom) then
          k = k + 1
          cube_add%atomic_number(k) = mycube2%atomic_number(i)
          cube_add%atomic_charge(k) = mycube2%atomic_charge(i)
          cube_add%x(k)             = mycube2%x(i)
          cube_add%y(k)             = mycube2%y(i)
          cube_add%z(k)             = mycube2%z(i)
       endif
    enddo

  END FUNCTION cube_add
  !--------------------------------------------------
  !--------------------------------------------------
  FUNCTION cube_sub (mycube1, mycube2)
    TYPE(cube) :: cube_sub
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    !
    ! Internal variables
    INTEGER                :: i,j,k,counter
    INTEGER                :: total_atoms
    !
    LOGICAL                :: unique_atom
    !
    cube_sub%rho  = mycube1%rho - mycube2%rho
    !
    !
    ! Assign remaining attributes to cube_sub
    cube_sub%str1          = mycube1%str1
    cube_sub%str2          = mycube1%str2
    !
    cube_sub%nx            = mycube1%nx
    cube_sub%ny            = mycube1%ny
    cube_sub%nz            = mycube1%nz
    !
    cube_sub%xmin          = mycube1%xmin
    cube_sub%ymin          = mycube1%ymin
    cube_sub%zmin          = mycube1%zmin
    cube_sub%dx            = mycube1%dx
    cube_sub%dy            = mycube1%dy
    cube_sub%dz            = mycube1%dz
    !
    ! Fill other quantities
    !
    !How many non-repeated atoms do we have?
    total_atoms = mycube1%natoms + mycube2%natoms
    do i = 1,mycube1%natoms
       do j = 1,mycube2%natoms
          if ((mycube1%x(i).eq. mycube2%x(j)) .and. (mycube1%y(i) .eq. mycube2%y(j)) &
                                              .and. (mycube1%z(i) .eq. mycube2%z(j))) then
             ! If all the positions are the same one atom is repeated
             total_atoms = total_atoms - 1
          endif
       enddo
    enddo
    !
    cube_sub%natoms = total_atoms
    allocate(cube_sub%x(cube_sub%natoms),cube_sub%y(cube_sub%natoms), &
             cube_sub%z(cube_sub%natoms),cube_sub%atomic_number(cube_sub%natoms),&
             cube_sub%atomic_charge(cube_sub%natoms))
   !
   ! Assign non-repeated atomic numbers, charges and positions
    ! From mycube1: assign all atoms
    k = 0
    do i = 1,mycube1%natoms
       k = k + 1
       cube_sub%atomic_number(k) = mycube1%atomic_number(i)
       cube_sub%atomic_charge(k) = mycube1%atomic_charge(i)
       cube_sub%x(k)             = mycube1%x(i)
       cube_sub%y(k)             = mycube1%y(i)
       cube_sub%z(k)             = mycube1%z(i)
    enddo
    !
    ! From mycube2: assign non-repeated atoms
    do i = 1,mycube2%natoms
       unique_atom = .false.
       counter = 0
       do j = 1,mycube1%natoms
          if ((mycube2%x(i).ne. mycube1%x(j)) .or. (mycube2%y(i) .ne. mycube1%y(j)) &
                                              .or. (mycube2%z(i) .ne. mycube1%z(j))) then
             ! If one of the positions are different, the atom is unique
             counter = counter + 1
                if (counter .eq. mycube1%natoms) unique_atom = .true.
          endif
       enddo
       ! If the atom is unique we take all its descriptors
       if (unique_atom) then
          k = k + 1
          cube_sub%atomic_number(k) = mycube2%atomic_number(i)
          cube_sub%atomic_charge(k) = mycube2%atomic_charge(i)
          cube_sub%x(k)             = mycube2%x(i)
          cube_sub%y(k)             = mycube2%y(i)
          cube_sub%z(k)             = mycube2%z(i)
       endif
    enddo

  END FUNCTION cube_sub
  !--------------------------------------------------
  !--------------------------------------------------
  FUNCTION cube_int (mycube)
    REAL (KIND=wp) :: cube_int
    TYPE (cube), INTENT(IN) :: mycube
    !
    ! Attributes of mycube (mycube%str1 = ...)
    !
    CHARACTER (LEN=72)                           :: str1,str2
    INTEGER                                      :: natoms,nx,ny,nz
    INTEGER, DIMENSION(:), ALLOCATABLE           :: atomic_number
    REAL(KIND=wp)                                :: xmin,ymin,zmin,dx,dy,dz,dummy_real
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE     :: atomic_charge,x,y,z
    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: rho
    !
    ! Other internal variables
    !
    INTEGER                                      :: i,j,k
    REAL(KIND=wp)                                :: x1,x2,y1,y2,z1,z2 
    !
    cube_int = 0d0
    do i = 1,mycube%nx
       do j = 1,mycube%ny
          do k = 1,mycube%nz
             cube_int = cube_int + mycube%rho(i,j,k)*mycube%dx*mycube%dy*mycube%dz
          enddo
       enddo
    enddo
    !
  END FUNCTION cube_int
  !--------------------------------------------------
  !--------------------------------------------------
  SUBROUTINE cube_del (mycube)
    TYPE (cube), INTENT(IN) :: mycube
    ! ...
  END SUBROUTINE cube_del
  !--------------------------------------------------
  SUBROUTINE cube_cdz (mycube,infile)
    !
    ! In/out variables
    !
    CHARACTER(LEN=*), INTENT(IN)               :: infile
    TYPE (cube),      INTENT(IN)               :: mycube
    !
    ! Attributes of mycube (mycube%str1 = ...)
    !
    CHARACTER (LEN=72)                           :: str1,str2
    INTEGER                                      :: natoms,nx,ny,nz
    INTEGER, DIMENSION(:), ALLOCATABLE           :: atomic_number
    REAL(KIND=wp)                                :: xmin,ymin,zmin,dx,dy,dz,dummy_real
    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE     :: atomic_charge,x,y,z
    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: rho
    !
    ! Other internal variables
    !
    INTEGER                                      :: OOut
    INTEGER                                      :: i,j,k
    !
    REAL (KIND=wp),ALLOCATABLE, DIMENSION(:)     :: cube_int_xy, cube_int_xyz
    !
    !
    ALLOCATE(cube_int_xy(mycube%nz),cube_int_xyz(mycube%nz))

    ! One dimensional function, integration ONLY in XY plane
    cube_int_xy(:) = 0d0
    do i = 1,mycube%nx
       do j = 1,mycube%ny
          do k = 1,mycube%nz
             cube_int_xy(k) = cube_int_xy(k) + mycube%rho(i,j,k)*mycube%dx*mycube%dy
          enddo
       enddo
    enddo
    !
    !
    ! One dimensional function, integration in XYZ
    ! BEFORE
    !!cube_int_xyz(:) = 0d0
    !!!do i = 1,mycube%nx
    !!!   do j = 1,mycube%ny
    !!!      do k = 1,mycube%nz
    !!!         cube_int_xyz(k) = cube_int_xyz(k) + mycube%rho(i,j,k)*mycube%dx*mycube%dy*mycube%dz
    !!!      enddo
    !!!   enddo
    !!!enddo



    ! NOW TRYING
    cube_int_xyz(:) = 0d0
    do i = 1,mycube%nx
       do j = 1,mycube%ny
          do k = 1,mycube%nz
             cube_int_xyz(k) = cube_int_xyz(k) + mycube%rho(i,j,k)*mycube%dx*mycube%dy*mycube%dz
          enddo
       enddo
    enddo






    !
    ! Print final result for the two 1-D arrays
    OOut = 11
    OPEN (UNIT=OOut, FILE=infile, STATUS="unknown", ACTION="write")   
    do i = 1,mycube%nz
       write(OOut,'(I3,1x,E15.8,1x,E15.8)') i, cube_int_xy(i), cube_int_xyz(i)
    enddo
    print *, cube_int_xy(1), cube_int_xyz(1)
    close(11)










  END SUBROUTINE cube_cdz
  !--------------------------------------------------
  !*************************************************
END MODULE cubes

!----------------------------------------------------------------------
module output_module
!      
!   Module output
!
    use target_module
    use parameters_module
!!!    use string_manipulation_module
!!!    use array_manipulation_module
!
    Implicit None
!
!   Public Variables
!
    public out_
!
!   output type
!
    type out_type
!
      integer  :: iunit = 12 ! unit of the file
!!      integer  :: ivrb  = 0  ! Verbose modality 
!
      character (len=12)  :: sticks = "(1x,80(1h-))"

      character (len=200) :: filename
!
      contains
!
      procedure :: out_file_fill
      procedure :: print_banner
      procedure :: print_density
      procedure :: print_results_integrals
!!      procedure :: print_matrix
      procedure :: warning
      procedure :: error
!
    end type out_type
!    
    type (out_type), Save :: out_
!
   contains
!----------------------------------------------------------------------
   subroutine out_file_fill(out_,in_file)
!
     implicit none
!     
!    define the name of output file  
!
     class(out_type)  :: out_
     character(len=*) :: in_file
!
     integer :: nlen
!    
     nlen = len_trim(in_file)
     write(out_%filename,'(a)') in_file(1:nlen-4) // '.log'
!
     out_%filename = trim(out_%filename)
!       
   end subroutine out_file_fill
!-----------------------------------------------------------------------
   subroutine print_banner(out_)
!
     implicit none
!     
!    print banner FRET Embedlab 
!
     class(out_type)  :: out_
!
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(11x,a)') "       ______ _____  ______ _______                    " 
     Write(out_%iunit,'(11x,a)') "      |  ____|  __ \|  ____|__   __|                   "                           
     Write(out_%iunit,'(11x,a)') "      | |__  | |__) | |__     | |                      "
     Write(out_%iunit,'(11x,a)') "      |  __| |  _  /|  __|    | |                      "
     Write(out_%iunit,'(11x,a)') "      | |    | | \ \| |____   | |                      "
     Write(out_%iunit,'(11x,a)') "      |_|    |_|  \_\______|  |_|                      "
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,'(11x,a)') "       ______           _              _ _       _     "
     Write(out_%iunit,'(11x,a)') "      |  ____|         | |            | | |     | |    "
     Write(out_%iunit,'(11x,a)') "      | |__   _ __ ___ | |__   ___  __| | | __ _| |__  "
     Write(out_%iunit,'(11x,a)') "      |  __| | '_ ` _ \| '_ \ / _ \/ _` | |/ _` | '_ \ "
     Write(out_%iunit,'(11x,a)') "      | |____| | | | | | |_) |  __/ (_| | | (_| | |_) |"
     Write(out_%iunit,'(11x,a)') "      |______|_| |_| |_|_.__/ \___|\__,_|_|\__,_|_.__/ "
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,'(25x,a)') "Program by Pablo Grobas Illobre"
     Write(out_%iunit,'(23x,a)') "with contributions by Sveva Sodomaco"
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " "
     Flush(out_%iunit)
!
!
   end subroutine print_banner
!----------------------------------------------------------------------
!!!   subroutine print_matrix(out_,string,matrix,idim1,idim2)
!!!!
!!!     implicit none
!!!!     
!!!!    print matrix 
!!!!
!!!     class(out_type)  :: out_
!!!!
!!!     integer, intent(in)                             :: idim1
!!!     integer, intent(in)                             :: idim2
!!!     real(dp), dimension(idim1,idim2), intent(in)    :: matrix
!!!     character(len=*)                                :: string
!!!!     
!!!!    internal
!!!!
!!!     integer                                 :: i, j, k
!!!     integer                                 :: ierr
!!!     integer                                 :: NVec
!!!     integer                                 :: NElements
!!!     integer                                 :: len_string
!!!     character(len=7)                        :: format_string
!!!!     
!!!     real(dp), dimension(:,:,:), allocatable :: dummy
!!!     real(dp), dimension(:,:,:), allocatable :: dummy2
!!!!     
!!! 2000 Format(3x,5(8x,i5))
!!!
!!!!
!!!!    Count how many vectors are there
!!!!
!!!     NVec      = ceiling(float(idim2)/five)
!!!     NElements = mod(idim2,5)
!!!     if(nelements.eq.0) nelements = 5
!!!!     
!!!     allocate(dummy(idim1,NVec-1,5),stat=ierr)
!!!     if(ierr.gt.0) call out_%error('not enough space in print_matrix')
!!!!     
!!!     allocate(dummy2(idim1,1,NElements),stat=ierr)
!!!     if(ierr.gt.0) call out_%error('not enough space in print_matrix')
!!!!     
!!!     if(NVec.ne.1) then
!!!        !$omp parallel do collapse(3) 
!!!        do i = 1, NVec-1
!!!           do j = 1, idim1
!!!              do k = 1, 5
!!!                 dummy(j,i,k) = matrix(j,k+(i-1)*5)
!!!              enddo
!!!           enddo
!!!        enddo
!!!       !$omp end parallel do
!!!     endif
!!!!     
!!!     do j = 1, idim1
!!!        do k = 1, nelements
!!!           dummy2(j,1,k) = matrix(j,k+(NVec-1)*5)
!!!        enddo
!!!     enddo    
!!!!
!!!     len_string = len(string)
!!!     write(format_string,'(a,i2,a)') "(",40-len_string/2,"X,a)"
!!!!
!!!     write(out_%iunit,out_%sticks)
!!!     write(out_%iunit,format_string) string
!!!     write(out_%iunit,out_%sticks)
!!!!     
!!!     do i = 1, nvec
!!!!     
!!!        write(out_%iunit,'(a)') ' '
!!!!     
!!!        if(i.ne.nvec) write(out_%iunit,2000) (k+(i-1)*5,k=1,5)
!!!        if(i.eq.nvec) write(out_%iunit,2000) (k+(i-1)*5,k=1,NElements)
!!!!     
!!!        do j = 1, idim1
!!!!     
!!!           if(i.ne.NVec) write(out_%iunit,'(i4,4x,5(e11.4,2x))') j,(dummy(j,i,k), k = 1, 5 )
!!!!     
!!!           if(i.eq.nvec) write(out_%iunit,'(i4,4x,5(e11.4,2x))') j,(dummy2(j,1,k), k = 1,NElements)
!!!!     
!!!        enddo
!!!!     
!!!     enddo
!!!!     
!!!     write(out_%iunit,out_%sticks)
!!!     flush(out_%iunit)
!!!!    
!!!   end subroutine print_matrix
!----------------------------------------------------------------------
   subroutine warning(out_,string)
!
     implicit none
!     
!    define the name of output file  
!
     class(out_type)  :: out_
     character(len=*) :: string
!
     integer          :: unit_
!    
!    check if the file is opened
!
     inquire(file=out_%filename,number=unit_)
!     
     if (unit_.eq.out_%iunit) then
        write(out_%iunit,'(/1x,a/)') "Warning! "//trim(string)
     else
        write(*,'(/1x,a/)') "Warning! "//trim(string)
     endif
!       
   end subroutine warning
!----------------------------------------------------------------------
   subroutine error(out_,string)
!
     implicit none
!     
!    define the name of output file  
!
     class(out_type)  :: out_
     character(len=*) :: string
!
     integer          :: unit_
!    
!    check if the file is opened
!
     inquire(file=out_%filename,number=unit_)
!     
     if (unit_.eq.out_%iunit) then
        write(out_%iunit,'(/1x,a)') "Error during the execution of FRET_Embedlab"
        write(out_%iunit,'(1x,a/)') trim(string)
        flush(out_%iunit)
        stop
     else
        write(*,'(/1x,a)') "Error during the execution of FRET_Embedlab"
        write(*,'(1x,a/)') trim(string)
        flush(out_%iunit)
        stop
     endif
!       
   end subroutine error
!-----------------------------------------------------------------------
   subroutine print_density(out_,cube_file,natoms,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,nelectrons,integral,header)
!   subroutine print_density(out_,cube)
!
!
     implicit none
!     
!    print density integral
!
!
!!    type (density_type), intent(in)  :: cube 
!!
     character(len=*), intent(in)  :: cube_file
     integer,  intent(in)          :: natoms,nelectrons,nx,ny,nz
     real(dp), intent(in)          :: dx,dy,dz
     real(dp), intent(in)          :: xmin,ymin,zmin

     real(dp),         intent(in), optional  :: integral
     character(len=*), intent(in), optional  :: header

     class(out_type)  :: out_
!
!    formats
!
     1000 Format(3x,I5,1x,E15.7,1x,E15.7,1x,E15.7)
!
     Write(out_%iunit,'(a)') " "
     if(PRESENT(header)) then
        Write(out_%iunit,'(22x,a)') header
     else
        Write(out_%iunit,'(22x,a)') '       Density Information                    ' 
     endif
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,'(3x,a)') "Density File: "//trim(cube_file)
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,'(3x,a)') "Density Grid (CUBE format): "
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit, 1000) natoms, xmin, ymin, zmin
     Write(out_%iunit, 1000) nx,   dx, zero, zero
     Write(out_%iunit, 1000) ny, zero,   dy, zero
     Write(out_%iunit, 1000) nz, zero, zero,   dz
     Write(out_%iunit,'(a)') " " 
     Write(out_%iunit,*)     "    Total number of grid points: ", nx*ny*nz
     Write(out_%iunit,'(a)') " " 
     if(PRESENT(integral)) then
        Write(out_%iunit,'(a)') "    ======================================================="
        Write(out_%iunit,*) "    Integrated electron density --> ", integral
        Write(out_%iunit,*) "    Total electrons in molecule --> ", nelectrons
        Write(out_%iunit,'(a)') "    ======================================================="
        Write(out_%iunit,'(a)') " " 
     endif
     Write(out_%iunit,out_%sticks) 

     Flush(out_%iunit)
!
!
   end subroutine print_density
!-----------------------------------------------------------------------
   subroutine print_results_integrals(out_,aceptor_donor_coulomb,aceptor_donor_overlap)
!
!    print integral results
!
     implicit none
!
!    input variables
!     
     class(out_type)  :: out_
!
     real(dp), optional  :: aceptor_donor_coulomb
     real(dp), optional  :: aceptor_donor_overlap
!
!
!    internal variables
! 
     real(dp), dimension(2)  :: v_tot ! 1: real, 2: imaginary
!
     real(dp)                :: v_mod
!
     v_tot = zero
!
     Write(out_%iunit,'(a)') " " 
     Write(out_%iunit,'(18x,a)') '                  RESULTS                    ' 
     Write(out_%iunit,'(a)') " " 
     Write(out_%iunit,out_%sticks) 
     if(PRESENT(aceptor_donor_coulomb)) then
        Write(out_%iunit,'(a)') " "
        Write(out_%iunit,'(5x,a,e15.6,a)') "Aceptor-Donor Coulomb:  ", aceptor_donor_coulomb, '  a.u.'
        Write(out_%iunit,'(5x,a,e15.6,a)') "Aceptor-Donor Overlap:  ", aceptor_donor_overlap, '  a.u.'
        v_tot(1) = aceptor_donor_coulomb + aceptor_donor_overlap
     endif
!
     v_mod = dsqrt(DOT_PRODUCT(v_tot,v_tot))
     Write(out_%iunit,'(30x,a)') " ------------------------------------"
     Write(out_%iunit,'(5x,a,e15.6,e15.6,a)') "Total Potential      :  ", v_tot(1), v_tot(2), ' i  a.u.'
     !Write(out_%iunit,'(5x,a,e15.6,a)') "Total Potential Modulus:", v_mod, ' a.u.'
     Write(out_%iunit,'(a)') " " 
     if(target_%name_.ne.'aceptor_np') then
        Write(out_%iunit,'(5x,a,e15.6,a)') "Keet:", two * pi * (v_mod**two) * target_%spectral_overlap, '  a.u.'
     endif
     Write(out_%iunit,'(a)') " " 
     Write(out_%iunit,out_%sticks) 
     Flush(out_%iunit)
!
!
   end subroutine print_results_integrals
!----------------------------------------------------------------------
end module output_module

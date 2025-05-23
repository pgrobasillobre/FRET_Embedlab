!----------------------------------------------------------------------
module output_module
!      
!   Module output
!
    use target_module
    use parameters_module
!
    Implicit None
!
!   Public Variables
!
    public print_density_file
    public print_error
    public out_
!
!   output type
!
    type out_type
!
      integer  :: iunit = iuout ! unit of the file
!
      character (len=12)  :: sticks = "(1x,80(1h-))"

      character (len=200) :: filename
!
      contains
!
      procedure :: out_file_fill
      procedure :: print_banner
      procedure :: print_density
      procedure :: print_nanoparticle
      procedure :: print_results_integrals
      procedure :: warning
      procedure :: error
!
    end type out_type
!    
    type (out_type), Save :: out_
!
    Interface print_density_file
       Module Procedure print_density 
    End Interface print_density_file
!
    Interface print_error
       Module Procedure error
    End Interface print_error
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
!    print banner FretLab 
!
     class(out_type)  :: out_
!
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,'(20x,a)') "    ______          __  __          __  " 
     Write(out_%iunit,'(20x,a)') "   / ____/_______  / /_/ /   ____ _/ /_ "                           
     Write(out_%iunit,'(20x,a)') "  / /_  / ___/ _ \/ __/ /   / __ `/ __ \"
     Write(out_%iunit,'(20x,a)') " / __/ / /  /  __/ /_/ /___/ /_/ / /_/ /"
     Write(out_%iunit,'(20x,a)') "/_/   /_/   \___/\__/_____/\__,_/_.___/ "
     Write(out_%iunit,'(20x,a)') "                                         "
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,'(25x,a)') "Program by Pablo Grobas Illobre"
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " "
     Flush(out_%iunit)
!
   end subroutine print_banner
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
        write(out_%iunit,'(/1x,a)') "Error during the execution of FretLab"
        write(out_%iunit,'(1x,a/)') trim(string)
        flush(out_%iunit)
        stop
     else
        write(*,'(/1x,a)') "Error during the execution of FretLab"
        write(*,'(1x,a/)') trim(string)
        flush(out_%iunit)
        stop
     endif
!       
   end subroutine error
!-----------------------------------------------------------------------
   subroutine print_density(out_,cube_file,cube,integral,header)
!
     use density_module, ONLY : density_type
!
     implicit none
!     
!    print density information
!
     class (density_type), intent(in)  :: cube 
!
     character(len=*), intent(in)  :: cube_file
!
     real(dp),         intent(in), optional  :: integral
     character(len=*), intent(in), optional  :: header

     integer :: i

     class(out_type)  :: out_
!
!    formats
!
     1000 Format(3x,I5,1x,E15.7,1x,E15.7,1x,E15.7)
     1001 Format(7x,a,2x,f12.6,2x,f12.6,2x,f12.6)
!
     write(out_%iunit,out_%sticks) 
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
     Write(out_%iunit, 1000) cube%natoms, cube%xmin ,   cube%ymin ,   cube%zmin
     Write(out_%iunit, 1000) cube%nx,     cube%dx(1),   cube%dx(2),   cube%dx(3)
     Write(out_%iunit, 1000) cube%ny,     cube%dy(1),   cube%dy(2),   cube%dy(3)
     Write(out_%iunit, 1000) cube%nz,     cube%dz(1),   cube%dz(2),   cube%dz(3)
     Write(out_%iunit,'(a)') " " 
     Write(out_%iunit,*)     "    Total number of grid points: ", cube%nx*cube%ny*cube%nz
     If (target_%name_ .ne. "integrate_density") then
        Write(out_%iunit,*)     "    ---> Reduced density points:", cube%n_points_reduced
     Endif
     Write(out_%iunit,'(a)') " " 
     Write(out_%iunit,'(3x,a)') "Associated molecular coordinates (Å): "
     Write(out_%iunit,'(a)') " "
     Do i=1,cube%natoms
        Write(out_%iunit,1001) cube%atomic_label(i), cube%x(i)*ToAng, cube%y(i)*ToAng, cube%z(i)*ToAng
     Enddo
     Write(out_%iunit,'(a)') " "

     If (header.eq.aceptor_header .and. target_%aceptor_transdip_rotate) then
        Write(out_%iunit,'(3x,a,1x,f10.5,1x,f10.5,1x,f10.5)') "INPUT   Transition density dipole (x,y,z): ", &
                                                                target_%aceptor_transdip(1),                   &
                                                                target_%aceptor_transdip(2),                   &
                                                                target_%aceptor_transdip(3)
        Write(out_%iunit,'(a)') " "
        Write(out_%iunit,'(3x,a,1x,f10.5,1x,f10.5,1x,f10.5)') "ROTATED Transition density dipole (x,y,z): ", &
                                                                target_%aceptor_transdip_rot(1),               &
                                                                target_%aceptor_transdip_rot(2),               &
                                                                target_%aceptor_transdip_rot(3)
        Write(out_%iunit,'(a)') " "
        Write(out_%iunit,'(3x,a,1x,f8.3,1x,a)')  "Alignment angle:", target_%aceptor_angle_check, "°"
        Write(out_%iunit,'(a)') " "

     Else If (header.eq.donor_header .and. target_%donor_transdip_rotate) then
        Write(out_%iunit,'(3x,a,1x,f10.5,1x,f10.5,1x,f10.5)') "INPUT   Transition density dipole (x,y,z): ", &
                                                                target_%donor_transdip(1),                   &
                                                                target_%donor_transdip(2),                   &
                                                                target_%donor_transdip(3)
        Write(out_%iunit,'(a)') " "
        Write(out_%iunit,'(3x,a,1x,f10.5,1x,f10.5,1x,f10.5)') "ROTATED Transition density dipole (x,y,z): ", &
                                                                target_%donor_transdip_rot(1),               &
                                                                target_%donor_transdip_rot(2),               &
                                                                target_%donor_transdip_rot(3)
        Write(out_%iunit,'(a)') " "
        Write(out_%iunit,'(3x,a,1x,f8.3,1x,a)')  "Alignment angle:", target_%donor_angle_check, "°"
        Write(out_%iunit,'(a)') " "
    EndIf

     if(PRESENT(integral)) then
        Write(out_%iunit,'(a)') "    ============================================================"
        Write(out_%iunit,*) "    Integrated electron density --> ", integral
        Write(out_%iunit,'(a)') "    ============================================================"
        Write(out_%iunit,'(a)') " " 
     endif

     Flush(out_%iunit)
!
   end subroutine print_density
!-----------------------------------------------------------------------
   subroutine print_nanoparticle(out_,natoms,xyz,charges,dipoles)
!
     implicit none
!     
!    print nanoparticle information
!
!
     integer,  intent(in)                       :: natoms
!
     real(dp), dimension(3,natoms), intent(in)  :: xyz
!
     logical :: charges, dipoles
!
     class(out_type)  :: out_
!
!    internal variables
!
     integer :: i
!
!    formats
!
     character(len=76) :: format_1 = "(13x,'Atom',15x,'X',19x,'Y',19x,'Z')"
     character(len=50) :: format_2 = "(12x,a4,1x,3(f20.6))"
!

     if (charges) write(out_%iunit,'(23x,a)') "NP model = charges"
     if (dipoles) write(out_%iunit,'(23x,a)') "NP model = charges + dipoles"
     Write(out_%iunit,'(a)') " "

     write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " "
     Write(out_%iunit,'(21x,a)') '       Nanoparticle Geometry (Å)                    ' 
     Write(out_%iunit,'(a)') " "

     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,format_1)
     Write(out_%iunit,out_%sticks) 

     Write(out_%iunit,'(a)') " "
     do i = 1, natoms
        write(out_%iunit,format_2) 'Xx', xyz(1,i)*ToAng, &
                                         xyz(2,i)*ToAng, &
                                         xyz(3,i)*ToAng
     enddo
     Write(out_%iunit,'(a)') " "

     Flush(out_%iunit)
!
!
   end subroutine print_nanoparticle
!-----------------------------------------------------------------------
 
   subroutine print_results_integrals(out_,aceptor_donor_coulomb,aceptor_donor_overlap,&
                                      aceptor_np_int)
!
!    print integral results
!
     implicit none
!
!    input variables
!     
     class(out_type)  :: out_
!
     real(dp), optional                :: aceptor_donor_coulomb
     real(dp), optional                :: aceptor_donor_overlap
     real(dp), dimension(2), optional  :: aceptor_np_int
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

     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " " 
     Write(out_%iunit,'(18x,a)') '                  RESULTS                    ' 
     Write(out_%iunit,'(a)') " " 
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(a)') " "
!
     if(PRESENT(aceptor_donor_coulomb)) then
        Write(out_%iunit,'(5x,a,f25.16,a)') "Aceptor-Donor Coulomb :  ", aceptor_donor_coulomb, '  a.u.'
        if (target_%calc_overlap_int) then
            Write(out_%iunit,'(5x,a,f25.16,a)') "Aceptor-Donor Overlap :  ", aceptor_donor_overlap, '  a.u.'
        endif
        v_tot(1) = aceptor_donor_coulomb + aceptor_donor_overlap
     endif
!
!    IMPORTANT: ADF prints the transition densities with opposite sign,
!               so we correct the sign of the potential integral
     if(PRESENT(aceptor_np_int)) then
        Write(out_%iunit,'(5x,a,f25.16,a,f25.16,a)') "Aceptor-NP Interaction:  ", -aceptor_np_int(1), '    + ',&
                                                                                  -aceptor_np_int(2) , ' i  a.u.'
        v_tot(1) = v_tot(1) - aceptor_np_int(1)
        v_tot(2) = v_tot(2) - aceptor_np_int(2)
     endif
!
     v_mod = dsqrt(DOT_PRODUCT(v_tot,v_tot))
!
     if (target_%name_.ne.'aceptor_np') then
         if (target_%name_.eq.'aceptor_np_donor') then 
            Write(out_%iunit,'(34x,a)') " -------------------------------------------------------------"
            Write(out_%iunit,'(5x,a,f25.16,a,f25.16,a)') "Total Potential       :  ", v_tot(1), '    + ', v_tot(2), ' i  a.u.'
         elseif (target_%name_.eq.'aceptor_donor') then 
            Write(out_%iunit,'(35x,a)') " --------------------------"
            Write(out_%iunit,'(5x,a,f25.16,a)') "Total Potential       :  ", v_tot(1), '  a.u.'
         endif
         Write(out_%iunit,'(a)') " "
         Write(out_%iunit,'(5x,a,f25.16,a)') "Total Potential Modulus:", v_mod, ' a.u.'
     endif
!
     Write(out_%iunit,'(a)') " " 
     if(target_%name_.ne.'aceptor_np') then
        Write(out_%iunit,'(5x,a,f25.16,a)') "Keet:", two * pi * (v_mod**two) * target_%spectral_overlap, '  a.u.'
        Write(out_%iunit,'(a)') " " 
     endif
     Write(out_%iunit,out_%sticks) 
     Flush(out_%iunit)
!
!
   end subroutine print_results_integrals
!----------------------------------------------------------------------
end module output_module

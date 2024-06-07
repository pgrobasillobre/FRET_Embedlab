!----------------------------------------------------------------------
module input_module
!      
!   Module input
!
    !$ use omp_lib
    use output_module
    use parallel_module
    use string_manipulation_module
    use target_module
!
    Implicit None
!
!   Public Variables
!
    public inp_
!
!   bem type
!
    type inp_type
!
      integer  :: iunit = 10 ! unit of the input file
!
      integer  :: NArg
!
      character (len=200) :: filename
!
      contains
!
      procedure :: get_arguments
      procedure :: check_input_file
      procedure :: read_
      procedure :: get_target
      procedure :: print_input_info
!
    end type inp_type
!    
    type (inp_type), Save :: inp_
!
   contains
!----------------------------------------------------------------------
   subroutine get_arguments(inp_)
!   
!    Grep the input file
!   
     implicit none
!
     class(inp_type)  :: inp_
!
!     
     character(len=100)  :: junk
!     
     !$ parallel%n_threads_OMP = omp_get_max_threads()
!     
     inp_%NArg=command_argument_count()
!
     If (inp_%NArg .eq. 1) then 
!
        Call get_command_argument(inp_%narg,inp_%filename)
        call out_%out_file_fill(inp_%filename)
!
     else If (inp_%NArg .gt. 1) then
!
        Call get_command_argument(1,inp_%filename)
        call out_%out_file_fill(inp_%filename)
        Call get_command_argument(2,junk)

!
!       Give number of procs for parallel
!
        If(trim(junk).eq.'-omp') then 
!
           If(inp_%NArg .ge. 3) then
!
              Call get_command_argument(3,junk)
              read(junk,'(i3)') parallel%n_threads_OMP
              If(inp_%NArg .eq. 4) Call get_command_argument(4,out_%filename)
!
           else
!
              Write(*,*) '-omp but no number given'
              Stop
!
           endif
!
!       Give input_name output_name [-omp #]
!
        else
!
           Call get_command_argument(2,out_%filename)
!
           If(inp_%NArg.gt.2) then
!
              Call get_command_argument(3,junk)
!
              If(trim(junk).ne.'-omp') then
!
                 Write(*,*) 'I do not know what option is: ', trim(junk)
!
              else
!
                 If(inp_%NArg .eq. 4) then
!
                    Call get_command_argument(4,junk)
                    read(junk,'(i3)') parallel%n_threads_OMP
!
                 else
!
                    Write(*,*) '-omp but no number given'
                    Stop
!
                 endif
!
              endif
!
           endif
!
        endif
!
     else If (inp_%NArg .eq. 0) then 
!
        Write(*,'(/a/)') "Type the inp filename (e.g. filename.inp)"
        Read(*,*) inp_%filename
        Write(*,'(/a/)') "Type the log filename (e.g. filename.log)"
        Read(*,*) out_%filename
        inp_%NArg = 2
!
     endIf
!
     inp_%filename = trim(inp_%filename)
!  
   end subroutine get_arguments
!-----------------------------------------------------------------------
   subroutine check_input_file(inp_)
!   
!    Check if file correctly exists and has the right extension
!   
     implicit none
!
     class(inp_type)  :: inp_
!
     integer          :: NLen
     character(len=3) :: check
     logical          :: exists
!    
     Inquire(file=inp_%filename,exist=exists)
!    
     If(.not.exists) call out_%error("File "//trim(inp_%filename)//" does not exist. Check the name")
!    
     NLen = len_trim(inp_%filename)
!    
     Write(check,'(a)') inp_%filename(Nlen-2:Nlen)
!    
     call lower(check)
!    
     If(check.ne.'inp') call out_%warning("The extension of your file is wrong. Use .inp")
!  
   end subroutine check_input_file
!---------------------------------------------------------------------------
   subroutine read_(inp_)
!  
!    Reading input_file
!     
     implicit none
     class(inp_type)  :: inp_
!
!    stuff for new input
!
     integer :: nlines
     integer :: num_string_input
     logical :: found_string
!
     character(len=200), dimension(:), allocatable             :: AtomName_geometry
!
     integer :: iost
!
!
     inp_%iunit                = 10! file
     
!    paramters for Iterative defaults
!
     open(unit=inp_%iunit,file=trim(inp_%filename),status="old",iostat=iost,err=01) 
!
          Rewind(inp_%iunit) 
          call get_number_lines(inp_%iunit,nlines)
!
          call inp_%get_target(nlines)
!
     01 continue 
     close(10) 
!
   end subroutine read_
!-----------------------------------------------------------------------
   subroutine get_target(inp_,nlines)
!
!    get target 
!
     implicit none
!
!    input variables
!
     class(inp_type)        :: inp_
     integer, intent(in)    :: nlines
!
!
!    internal variables
!
     integer :: num_string_initial, num_string_end, num_keywords
     integer :: i
!
     logical :: found_string
     logical :: aceptor_density_rotate
     logical :: donor_density_rotate
!
     character(len=200) :: line
     character(len=200) :: line_what
     character(len=200) :: line_keyword
     character(len=200) :: target_name = ""
!
     logical :: cutoff                 = .false.
     logical :: integrate_cube         = .false.
     logical :: fret_donor_aceptor     = .false.
     logical :: fret_donor_aceptor_NP  = .false.
     logical :: fret_aceptor_NP        = .false.
     logical :: aceptor_density        = .false.
     logical :: donor_density          = .false.
     logical :: nanoparticle           = .false.
     logical :: omega_0                = .false.
     logical :: spectral_overlap       = .false.
     logical :: rotation_axys          = .false.
!
     target_%cutoff  = zero
     target_%omega_0 = zero
     target_%spectral_overlap = zero
     target_%aceptor_density_rotation_angle = zero
     target_%donor_density_rotation_angle = zero
!
     target_%aceptor_transdip     = zero
     target_%aceptor_transdip_rot = zero 
     target_%aceptor_ref_vector   = zero
!
     target_%donor_transdip     = zero
     target_%donor_transdip_rot = zero 
     target_%donor_ref_vector   = zero
!
     target_%calc_overlap_int   = .false.
!
     aceptor_density_rotate                     = .false.
     target_%rotate_aceptor                     = .false.
     target_%aceptor_transdip_rotate            = .false.
     target_%aceptor_transdip_rotate_align_with = .false.
!
     donor_density_rotate                       = .false.
     target_%rotate_donor                       = .false.
     target_%donor_transdip_rotate              = .false.
     target_%donor_transdip_rotate_align_with   = .false.
!
     target_%rotate_np = .false.
!
     target_%debug = .false.
!
     rewind(inp_%iunit)
!
        do i = 1, nlines
!
          read(inp_%iunit,'(a)') line
!        
          if(.not.commented_line(line)) then
!        
            call lower(line)
!          
            if(index(line,":").gt.0) then ! : in input
!        
              line_what = trim(adjustl(line(1:index(line,":")-1) ))
!            
              line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!        
              if(line_what.eq.'integrate cube file')  then 
!
                 integrate_cube  = .true.
                 target_%density_file = line_keyword
                 if (.not. file_exists(line_keyword)) call out_%error('File "'//trim(line_keyword)//'" not found')
!
              elseif(line_what.eq.'aceptor density') then
!
                 aceptor_density = .true.
                 target_%aceptor_density = line_keyword
                 if (.not. file_exists(line_keyword)) call out_%error('File "'//trim(line_keyword)//'" not found')
!
              elseif(line_what.eq.'aceptor density rotate') then
!
                 aceptor_density_rotate = .true.
                 read(line_keyword,'(f25.16)') target_%aceptor_density_rotation_angle
                 target_%aceptor_density_rotation_angle = target_%aceptor_density_rotation_angle * to_radians ! Assumed to be in degrees
!
              elseif(line_what.eq.'aceptor transition dipole') then
!
                 target_%aceptor_transdip_rotate = .true.
                 read(line_keyword, *) target_%aceptor_transdip(1), &
                                       target_%aceptor_transdip(2), &
                                       target_%aceptor_transdip(3)
!
              elseif(line_what.eq.'aceptor transition dipole align with') then
!
                 target_%aceptor_transdip_rotate_align_with = .true.
                 read(line_keyword, *) target_%aceptor_ref_vector(1), &
                                       target_%aceptor_ref_vector(2), &
                                       target_%aceptor_ref_vector(3)
!
              elseif(line_what.eq.'donor density') then
!
                 donor_density = .true.
                 target_%donor_density = line_keyword
                 if (.not. file_exists(line_keyword)) call out_%error('File "'//trim(line_keyword)//'" not found')
!
              elseif(line_what.eq.'donor density rotate') then
!
                 donor_density_rotate = .true.
                 read(line_keyword,'(f25.16)') target_%donor_density_rotation_angle
                 target_%donor_density_rotation_angle = target_%donor_density_rotation_angle * to_radians ! Assumed to be in degrees
!
              elseif(line_what.eq.'donor transition dipole') then
!
                 target_%donor_transdip_rotate = .true.
                 read(line_keyword, *) target_%donor_transdip(1), &
                                       target_%donor_transdip(2), &
                                       target_%donor_transdip(3)
!
              elseif(line_what.eq.'donor transition dipole align with') then
!
                 target_%donor_transdip_rotate_align_with = .true.
                 read(line_keyword, *) target_%donor_ref_vector(1), &
                                       target_%donor_ref_vector(2), &
                                       target_%donor_ref_vector(3)
!
              elseif(line_what.eq.'nanoparticle') then
!
                 nanoparticle = .true.
                 target_%nanoparticle = line_keyword
                 if (.not. file_exists(line_keyword)) call out_%error('File "'//trim(line_keyword)//'" not found')
!
              elseif(line_what.eq.'cutoff') then
                 cutoff = .true.
!
                 call check_float(out_%iunit,line_what,line_keyword)
                 read(line_keyword,'(f25.16)') target_%cutoff
!
                 if (target_%cutoff .lt. zero) call out_%error("Cutoff cannot be negative")
!
              elseif(line_what.eq.'omega_0') then
                 omega_0 = .true.
!
                 call check_float(out_%iunit,line_what,line_keyword)
                 read(line_keyword,'(f25.16)') target_%omega_0
!
              elseif(line_what.eq.'spectral overlap') then
                 spectral_overlap = .true.
!
                 call check_float(out_%iunit,line_what,line_keyword)
                 read(line_keyword,'(f25.16)') target_%spectral_overlap
!
              elseif(line_what.eq.'debug') then
                 if (line_keyword.eq.'true') target_%debug = .true.
!
              elseif(line_what.eq.'rotation axys') then
                 rotation_axys = .true.
!
                 target_%rotation_axys = line_keyword
                 if (target_%rotation_axys .ne. 'x' .and. &
                     target_%rotation_axys .ne. 'y' .and. &
                     target_%rotation_axys .ne. 'z') call out_%error( 'Rotation axys "'//trim(line_keyword)//'" not recognized')
!
              else  
!
                 call out_%error( 'Input entry "'//trim(line)//'" not recognized')

              endif
!
            else ! : in input
!          
              call out_%error( 'I am not able to find : in entry "'//trim(line)//'"')
!            
            endif
!
          endif
!        
        end do
!
!    Checks and rising errors aceptor
!
     if (aceptor_density_rotate .or. target_%aceptor_transdip_rotate) target_%rotate_aceptor = .true. 
!
     if(target_%aceptor_transdip_rotate .and. .not. aceptor_density_rotate .and. &
        target_%aceptor_transdip_rotate .and. .not. target_%aceptor_transdip_rotate_align_with) then
        call out_%error("Provide angle/reference alignment vector for transition dipole (and density) rotation")
!
     elseif(aceptor_density_rotate             .and. .not. rotation_axys .or. &
            target_%aceptor_transdip_rotate_align_with .and. .not. rotation_axys) then 
        call out_%error("Rotation of aceptor/donor density but no rotation axys given in input")
!
     elseif(target_%aceptor_transdip_rotate_align_with .and. .not. target_%aceptor_transdip_rotate) then
        call out_%error("No acceptor transition density provided in input")

     elseif(target_%aceptor_transdip_rotate_align_with .and. target_%aceptor_density_rotation_angle > zero) then
        call out_%error("Acceptor rotation angle provided but conflicts with 'align with' option")
!
     endif
!
!    Checks and rising errors donor
!
     if (donor_density_rotate .or. target_%donor_transdip_rotate) target_%rotate_donor = .true. 
!
     if(target_%donor_transdip_rotate .and. .not. donor_density_rotate .and. &
        target_%donor_transdip_rotate .and. .not. target_%donor_transdip_rotate_align_with) then
        call out_%error("Provide angle/reference alignment vector for transition dipole (and density) rotation")
!
     elseif(donor_density_rotate             .and. .not. rotation_axys .or. &
            target_%donor_transdip_rotate_align_with .and. .not. rotation_axys) then 
        call out_%error("Rotation of aceptor/donor density but no rotation axys given in input")
!
     elseif(target_%donor_transdip_rotate_align_with .and. .not. target_%donor_transdip_rotate) then
        call out_%error("No donor transition density provided in input")

     elseif(target_%donor_transdip_rotate_align_with .and. target_%donor_density_rotation_angle > zero) then
        call out_%error("donor rotation angle provided but conflicts with 'align with' option")
!
     endif
!
!    Check and rising errors np
!
     if (donor_density_rotate) target_%rotate_np = .true.
!
!    Check rotation axys
!
     if (rotation_axys .and. .not. target_%rotate_aceptor .and. .not. target_%rotate_donor) &
          call out_%error("Rotation axys given but no rotate aceptor/donor selected")
!
!    Create folder to store debug quantities
!
     if (target_%debug) then
        if (file_exists("debug")) call execute_command_line("rm -rf debug")
        call execute_command_line("mkdir debug")
     endif
!
!    assign the different targets
!
     if (.not. cutoff .and. .not. integrate_cube .and. .not. cutoff .and. omega_0) then
        call out_%error("Cutoff needed in input")
!
     elseif(integrate_cube .and. aceptor_density .or. &
            integrate_cube .and. donor_density   .or. &
            integrate_cube .and. nanoparticle) then         
! 
        call out_%error( "You are requesting a cube integration with other type of calculation")
!
     elseif(integrate_cube) then
!
        target_%name_ = "integrate_density"
!
     elseif(donor_density .and. aceptor_density .and. .not. nanoparticle) then
!
        target_%name_ = "aceptor_donor"
        if (target_%omega_0 < zero) call out_%error("Omega_0 cannot be negative")
        if (target_%omega_0 > 1.0E-15) target_%calc_overlap_int = .true.
        !if (target_%omega_0 < 1E-14) call out_%error("Aceptor-donor calculation requested but no Omega_0 in input")
!
        if (target_%spectral_overlap < zero) call out_%error("Spectral overlap cannot be negative")
        if (target_%spectral_overlap < 1E-14)&
           call out_%error("Aceptor-donor calculation requested but no spectral overlap in input")
!
     elseif(aceptor_density .and. nanoparticle .and. .not. donor_density) then
!
        target_%name_ = "aceptor_np"
!
     elseif(aceptor_density .and. nanoparticle .and. donor_density) then
!
        target_%name_ = "aceptor_np_donor"
        if (target_%omega_0 < zero) call out_%error("Omega_0 cannot be negative")
        if (target_%omega_0 > 1.0E-15) target_%calc_overlap_int = .true.
        !if (target_%omega_0 < 1E-14) call out_%error("Aceptor-NP-donor calculation requested but no Omega_0 in input")
!
        if (target_%spectral_overlap < zero) call out_%error("Spectral overlap cannot be negative")
        if (target_%spectral_overlap < 1E-14)&
           call out_%error("Aceptor-NP-donor calculation requested but no spectral overlap in input")
!
     else
!         
        call out_%error( "You are requesting an unkown calculation")
!
     endif
!
   end subroutine get_target
!-----------------------------------------------------------------------
  subroutine print_input_info(inp_)
!
!   Printing input information on output
!
    implicit none
!
    class(inp_type) :: inp_
!
    integer :: i
!    
     write(out_%iunit,'(23x,a)') "Input  File: "//trim(inp_%filename)
     write(out_%iunit,'(23x,a)') "Output File: "//trim(out_%filename)
     write(out_%iunit,'(/,23x,a,i5)') "OMP Threads: ", parallel%n_threads_OMP
     write(out_%iunit,'(a)') ""
     write(out_%iunit,out_%sticks) 
!
!    general informations
!
     write(out_%iunit,'(a)') ""
     write(out_%iunit,'(23x,a)') "Calculation --> "//trim(target_%name_)
     write(out_%iunit,'(a)') ""
     !write(out_%iunit,'(23x,a,i1)') "Verbose       : ", out_%ivrb

     if(target_%name_.eq."integrate_density") then
        write(out_%iunit,'(23x,a)') "Density File: "//trim(target_%density_file)

     elseif(target_%name_.eq."aceptor_donor") then
        write(out_%iunit,'(23x,a)') "Aceptor Density : "//trim(target_%aceptor_density)
        write(out_%iunit,'(23x,a)') "Donor   Density : "//trim(target_%donor_density)
        write(out_%iunit,'(a)')
        write(out_%iunit,'(23x,a,e11.4)') "Cutoff            = ", target_%cutoff
        if (target_%calc_overlap_int) write(out_%iunit,'(23x,a,e11.4)') "Omega_0          = ", target_%omega_0
        write(out_%iunit,'(23x,a,e11.4)') "Spectral Overlap = ", target_%spectral_overlap

     elseif(target_%name_.eq."aceptor_np") then
        write(out_%iunit,'(23x,a)') "Aceptor Density  : "//trim(target_%aceptor_density)
        write(out_%iunit,'(23x,a)') "Nanoparticle File: "//trim(target_%nanoparticle)
        write(out_%iunit,'(a)')
        write(out_%iunit,'(23x,a,e11.4)') "Cutoff            = ", target_%cutoff

     elseif(target_%name_.eq."aceptor_np_donor") then
        write(out_%iunit,'(23x,a)') "Aceptor Density   : "//trim(target_%aceptor_density)
        write(out_%iunit,'(23x,a)') "Donor   Density   : "//trim(target_%donor_density)
        write(out_%iunit,'(23x,a)') "Nanoparticle File : "//trim(target_%nanoparticle)
        write(out_%iunit,'(a)')
        write(out_%iunit,'(23x,a,e11.4)') "Cutoff            = ", target_%cutoff
        if (target_%calc_overlap_int) write(out_%iunit,'(23x,a,e11.4)') "Omega_0           = ", target_%omega_0
        write(out_%iunit,'(23x,a,e11.4)') "Spectral Overlap  = ", target_%spectral_overlap

     endif

     write(out_%iunit,'(a)') ""
     !write(out_%iunit,out_%sticks) 
!     
     flush(out_%iunit)
! 
  end subroutine print_input_info
!----------------------------------------------------------------------
  function file_exists(filename) result(res)
!
!   Return true if file exist, false otherwise
!
    implicit none
!
    character(len=*),intent(in) :: filename
    logical                     :: res
!
    inquire( file=trim(filename), exist=res )
!
  end function
!----------------------------------------------------------------------
end module input_module

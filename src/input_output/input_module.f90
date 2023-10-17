!----------------------------------------------------------------------
module input_module
!      
!   Module input
!
    !$ use omp_lib
    !!!use output_module
    !!!use parameters_module
    !!!use string_manipulation_module
    !!!use array_manipulation_module
    !!!use target_module
    !!!use bem_module
    !!!use fq_module
    !!!use fqfmu_module
    !!!use wfq_module
    !!!use wfqfmu_module
    !!!use wfq_bem_module
    !!!use wfqfmu_bem_module
    !!!use algorithm_module
    !!!use general_module
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
      !!procedure :: check_input_file
      !!procedure :: read_
      !!procedure :: get_section_what
      !!procedure :: get_target
      !!procedure :: get_section_forcefield
      !!procedure :: get_section_algorithm
      !!procedure :: get_section_field
      !!procedure :: get_section_output
      !!procedure :: get_section_parameters
      !!procedure :: get_section_atomtypes
      !!procedure :: get_section_bem
      !!procedure :: get_section_general
      !!procedure :: get_section_geometry
      !!procedure :: get_number_of_atoms
      !!procedure :: print_input_info
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
print *, 'oli'
!!!!!     
!!!!     character(len=100)  :: junk
!!!!!     
!!!!     !$ algorithm%n_threads_OMP = omp_get_max_threads()
!!!!!     
!!!!     inp_%NArg=command_argument_count()
!!!!!
!!!!     If (inp_%NArg .eq. 1) then 
!!!!!
!!!!        Call get_command_argument(inp_%narg,inp_%filename)
!!!!        call out_%out_file_fill(inp_%filename)
!!!!!
!!!!     else If (inp_%NArg .gt. 1) then
!!!!!
!!!!        Call get_command_argument(1,inp_%filename)
!!!!        call out_%out_file_fill(inp_%filename)
!!!!        Call get_command_argument(2,junk)
!!!!!
!!!!        If(trim(junk).eq.'-omp') then 
!!!!!
!!!!           If(inp_%NArg .ge. 3) then
!!!!!
!!!!              Call get_command_argument(3,junk)
!!!!              read(junk,'(i3)') algorithm%n_threads_OMP
!!!!              If(inp_%NArg .eq. 4) Call get_command_argument(4,out_%filename)
!!!!!
!!!!           else
!!!!!
!!!!              Write(*,*) '-omp but no number given'
!!!!!
!!!!           endif
!!!!!
!!!!        else
!!!!!
!!!!           Call get_command_argument(2,out_%filename)
!!!!!
!!!!           If(inp_%NArg.gt.2) then
!!!!!
!!!!              Call get_command_argument(3,junk)
!!!!!
!!!!              If(trim(junk).ne.'-omp') then
!!!!!
!!!!                 Write(*,*) 'I do not know what option is: ', trim(junk)
!!!!!
!!!!              else
!!!!!
!!!!                 If(inp_%NArg .eq. 4) then
!!!!!
!!!!                    Call get_command_argument(4,junk)
!!!!                    read(junk,'(i3)') algorithm%n_threads_OMP
!!!!!
!!!!                 else
!!!!!
!!!!                    Write(*,*) '-omp but no number given'
!!!!!
!!!!                 endif
!!!!!
!!!!              endif
!!!!!
!!!!           endif
!!!!!
!!!!        endif
!!!!!
!!!!     else If (inp_%NArg .eq. 0) then 
!!!!!
!!!!        Write(*,'(/a/)') "Type the mfq filename (e.g. filename.mfq)"
!!!!        Read(*,*) inp_%filename(1:99)
!!!!        Write(*,'(/a/)') "Type the log filename (e.g. filename.log)"
!!!!        Read(*,*) out_%filename(1:99)
!!!!        inp_%NArg = 2
!!!!!
!!!!     endIf
!!!!!
!!!!     inp_%filename = trim(inp_%filename)
!!!!!  
   end subroutine get_arguments
!!!!!-----------------------------------------------------------------------
!!!!   subroutine check_input_file(inp_)
!!!!!   
!!!!!    Check if file correctly exists and has the right exstension
!!!!!   
!!!!     implicit none
!!!!!
!!!!     class(inp_type)  :: inp_
!!!!!
!!!!     integer          :: NLen
!!!!     character(len=3) :: check
!!!!     logical          :: exists
!!!!!    
!!!!     Inquire(file=inp_%filename,exist=exists)
!!!!!    
!!!!     If(.not.exists) call out_%error("File "//trim(inp_%filename)//" Does Not Exist. Check the name")
!!!!!    
!!!!     NLen = len_trim(inp_%filename)
!!!!!    
!!!!     Write(check,'(a)') inp_%filename(Nlen-2:Nlen)
!!!!!    
!!!!     call lower(check)
!!!!!    
!!!!     If(check.ne.'mfq') call out_%warning("The extension of your file is wrong. Use .mfq")
!!!!!  
!!!!   end subroutine check_input_file
!!!!!---------------------------------------------------------------------------
!!!!   subroutine read_(inp_)
!!!!!  
!!!!!    Reading input_file
!!!!!     
!!!!     implicit none
!!!!     class(inp_type)  :: inp_
!!!!!
!!!!!    stuff for new input
!!!!!
!!!!     integer :: nlines
!!!!     integer :: num_string_input
!!!!     logical :: found_string
!!!!!
!!!!     character(len=200), dimension(:), allocatable             :: AtomName_geometry
!!!!!
!!!!     integer :: iost
!!!!!
!!!!!
!!!!     inp_%iunit                = 10! file
!!!!     
!!!!!    paramters for Iterative defaults
!!!!!
!!!!     open(unit=inp_%iunit,file=trim(inp_%filename),status="old",iostat=iost,err=01) 
!!!!!
!!!!          Rewind(inp_%iunit) 
!!!!          call get_number_lines(inp_%iunit,nlines)
!!!!!
!!!!          call inp_%get_target(nlines)
!!!!!
!!!!          if(target_%name_.ne."bem") then
!!!!!        
!!!!             Rewind(inp_%iunit) 
!!!!             call go_to_string(inp_%iunit,"input geometry",found_string,nlines,num_string_input)
!!!!!
!!!!             if(.not.found_string) then 
!!!!                call out_%error('No input geometry section in input')
!!!!             endif
!!!!!
!!!!          else 
!!!!!
!!!!              num_string_input = nlines
!!!!!
!!!!          endif
!!!!!
!!!!          rewind(inp_%iunit)
!!!!!
!!!!          call inp_%get_section_what(num_string_input)
!!!!!
!!!!          if(target_%name_.ne."bem") &
!!!!             call inp_%get_section_forcefield(num_string_input)
!!!!!
!!!!          call inp_%get_section_algorithm(num_string_input)
!!!!!
!!!!          call inp_%get_section_field(num_string_input)
!!!!!
!!!!          call inp_%get_section_output(num_string_input)
!!!!!
!!!!          call inp_%get_section_general(num_string_input)
!!!!!
!!!!          call inp_%get_section_bem(num_string_input)         
!!!!!
!!!!          if(target_%name_.ne."bem") then
!!!!!
!!!!              call inp_%get_section_atomtypes(num_string_input)
!!!!!
!!!!              call inp_%get_section_parameters(num_string_input)
!!!!!
!!!!              call inp_%get_section_geometry(nlines,AtomName_geometry)
!!!!!
!!!!              call assign_parameters(AtomName_geometry)
!!!!!
!!!!             if(target_%name_.eq.'wfq'          .or. &
!!!!                target_%name_.eq.'wfqfmu'       .or. &
!!!!                target_%name_.eq.'wfqib'        .or. &
!!!!                target_%name_.eq.'wfq_bem'      .or. &
!!!!                target_%name_.eq.'wfqfmu_bem'   .or. &
!!!!                target_%name_.eq.'wfqib_bem') then
!!!!!                
!!!!                call target_%assign_model_parameters()
!!!!!
!!!!!lookhere
!!!!                if(target_%name_.eq.'wfqfmu' .or.    &
!!!!                   target_%name_.eq.'wfqfmu_bem') call wfqfmu%read_alpha()
!!!!!end here
!!!!!          
!!!!              endif
!!!!!          
!!!!              allocate(fq%n_atoms_per_molecule(fq%n_mol))
!!!!!
!!!!              call check_inconsistencies_input_file()
!!!!!
!!!!              call array_scale(tobohr,3*target_%n_atoms,target_%coord)
!!!!!                                   
!!!!              deallocate(AtomName_geometry)
!!!!!              
!!!!          endif
!!!!!
!!!!          If(bem%exists) call bem%create_tesserae()
!!!!!
!!!!          if(target_%name_.eq.'wfq_bem'   .or. &
!!!!             target_%name_.eq.'wfqfmu_bem'.or. &
!!!!             target_%name_.eq.'wfqib_bem') then
!!!!             target_%n_var = target_%n_var + bem%n_var
!!!!             call wfq_bem%check_boundaries()
!!!!          endif             
!!!!!     
!!!!     01 continue 
!!!!     close(10) 
!!!!!
!!!!   end subroutine read_
!!!!!-----------------------------------------------------------------------
!!!!   subroutine get_section_what(inp_,nlines)
!!!!!
!!!!!    get section what
!!!!!
!!!!     implicit none
!!!!!
!!!!!    input variables
!!!!!
!!!!     class(inp_type)  :: inp_
!!!!     integer, intent(in)             :: nlines
!!!!!
!!!!!    input/output variables
!!!!!
!!!!!
!!!!!    internal variables
!!!!!
!!!!     integer :: num_string_initial, num_string_end, num_keywords
!!!!     logical :: found_string
!!!!     integer :: i
!!!!     character(len=200) :: line
!!!!     character(len=200) :: line_what
!!!!!
!!!!     out_%what = ""
!!!!     call go_to_string(inp_%iunit,"what"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!     if(.not.found_string) call out_%error('No what section in input')
!!!!!
!!!!     rewind(inp_%iunit)
!!!!     call go_to_string(inp_%iunit,"end what",found_string,nlines,num_string_end)
!!!!!
!!!!     if(.not.found_string) then
!!!!        call out_%error('No end for what section in input')
!!!!     endif
!!!!!
!!!!     num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!     rewind(inp_%iunit)
!!!!     call go_to_string(inp_%iunit,"what"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!     do i = 1, num_keywords
!!!!!
!!!!       read(inp_%iunit,'(a)') line
!!!!!     
!!!!       if(len_trim(line).eq.0) then 
!!!!!          
!!!!          call out_%error("There is a blank line (what)")
!!!!!          
!!!!       else if(.not.commented_line(line)) then
!!!!!     
!!!!         call LoWeR(line)
!!!!!       
!!!!         line_what = trim(adjustl(line))
!!!!!     
!!!!         if(line_what.eq.'restart') then 
!!!!            out_%what = line_what
!!!!            general%restart = .true.
!!!!         else
!!!!            out_%what = line_what
!!!!!            
!!!!            if(out_%what.ne.'energy'.and.&
!!!!               out_%what.ne.'static response'.and.&
!!!!               out_%what.ne.'cross section'.and.&
!!!!               out_%what.ne.'restart'.and.&
!!!!               out_%what.ne.'current') then
!!!!               write(out_%iunit,'(1x,a/)') 'Keyword: "'//trim(out_%what)//'" in what section not recognized'
!!!!               write(out_%iunit,'(1x,a)' ) 'Possible keywords in what section:'
!!!!               write(out_%iunit,'(1x,a)' ) '  - energy'
!!!!               write(out_%iunit,'(1x,a)' ) '  - static response'
!!!!               write(out_%iunit,'(1x,a)' ) '  - cross section'
!!!!               write(out_%iunit,'(1x,a)' ) '  - current'
!!!!               write(out_%iunit,'(1x,a)' ) '  - restart'
!!!!               stop
!!!!            endif
!!!!!            
!!!!         endif
!!!!!
!!!!       endif
!!!!!     
!!!!     end do
!!!!!
!!!!!
!!!!   end subroutine get_section_what
!!!!!-----------------------------------------------------------------------
!!!!   subroutine get_target(inp_,nlines)
!!!!!
!!!!!    get section forcefield
!!!!!
!!!!     implicit none
!!!!!
!!!!!    input variables
!!!!!
!!!!     class(inp_type)  :: inp_
!!!!     integer, intent(in)             :: nlines
!!!!!
!!!!!
!!!!!    internal variables
!!!!!
!!!!     integer :: num_string_initial, num_string_end, num_keywords
!!!!     logical :: found_string
!!!!     integer :: i
!!!!     character(len=200) :: line
!!!!     character(len=200) :: line_what
!!!!     character(len=200) :: line_keyword
!!!!     character(len=200) :: target_name = ""
!!!!!
!!!!     logical :: forcefield_present = .false.
!!!!     logical :: dynamic_ff = .true.
!!!!!
!!!!     rewind(inp_%iunit)
!!!!     call go_to_string(inp_%iunit,"forcefield"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!     if(.not.found_string) then ! let's just see if it is a bem calculation
!!!!        rewind(inp_%iunit)
!!!!        call go_to_string(inp_%iunit,"bem",found_string,nlines,num_string_initial)
!!!!        if(found_string) then 
!!!!           bem%exists  = .true.
!!!!        else
!!!!           call out_%error( 'No forcefield or bem section in input')
!!!!        endif
!!!!     else 
!!!!        forcefield_present = .true.
!!!!     endif
!!!!!
!!!!     if(forcefield_present) then
!!!!!
!!!!        rewind(inp_%iunit)
!!!!        call go_to_string(inp_%iunit,"end forcefield",found_string,nlines,num_string_end)
!!!!!
!!!!        if(.not.found_string) &
!!!!           call out_%error( 'No end for forcefield section in input')
!!!!!
!!!!        num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!        rewind(inp_%iunit)
!!!!        call go_to_string(inp_%iunit,"forcefield"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!        do i = 1, num_keywords
!!!!!
!!!!          read(inp_%iunit,'(a)') line
!!!!!        
!!!!          if(len_trim(line).eq.0) then 
!!!!!             
!!!!             call out_%error( "There is a blank line (forcefield)")
!!!!!             
!!!!          else if(.not.commented_line(line)) then
!!!!!        
!!!!            call lower(line)
!!!!!          
!!!!            if(index(line,":").gt.0) then ! : in input
!!!!!        
!!!!              line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!            
!!!!              line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!        
!!!!              if(line_what.eq.'static')  then 
!!!!                 target_name  = trim(adjustl(line_keyword))
!!!!                 if(target_name.ne.'fq'.and. &
!!!!                    target_name.ne.'fq_pqeq'.and. &
!!!!                    target_name.ne.'fqfmu'.and.   &
!!!!                    target_name.ne.'fqfmu_pqeq') &
!!!!                    call out_%error( "ForceField "//trim(target_name)//" not recognised")
!!!!              else if(line_what.eq.'dynamic') then 
!!!!                 dynamic_ff = .true.
!!!!                 target_name  = trim(adjustl(line_keyword))
!!!!                 if(line_keyword.ne.'wfq'   .and.&
!!!!                    line_keyword.ne.'wfqfmu'.and.&
!!!!                    line_keyword.ne.'wfqib') &
!!!!                    call out_%error( "wForceField = "//trim(line_keyword)//" not recognized")
!!!!              endif
!!!!!
!!!!              if(line_what.ne.'static'.and.&
!!!!                 line_what.ne.'dynamic'.and.&
!!!!                 line_what.ne.'kernel') then
!!!!                 write(out_%iunit,'(1x,a/)') 'Keyword : "'//trim(line_what)//'" in ForceField section not recognized'
!!!!                 write(out_%iunit,'(1x,a)' ) 'Possible Keywords in ForceField section:'
!!!!                 write(out_%iunit,'(1x,a)' ) '  - static:  [fq,fqfmu,fq_pqeq,fqfmu_pqeq]'
!!!!                 write(out_%iunit,'(1x,a)' ) '  - dynamic: [wfq,wfqfmu,wfqib]'
!!!!                 write(out_%iunit,'(1x,a)' ) '  - kernel:  [coulomb,ohno,gaussian]'
!!!!                 stop
!!!!              endif
!!!!!        
!!!!            else ! : in input
!!!!!          
!!!!              call out_%error( "I am not able to find : in keyword (forcefield)")
!!!!!            
!!!!            endif
!!!!!
!!!!          endif
!!!!!        
!!!!        end do
!!!!!     
!!!!        rewind(inp_%iunit)
!!!!        call go_to_string(inp_%iunit,"bem",found_string,nlines,num_string_initial)
!!!!        if(found_string) bem%exists = .true.
!!!!!
!!!!     endif
!!!!!     
!!!!!    assign the different targets
!!!!!
!!!!     if(target_name.eq.'fq'.or.target_name.eq.'fq_pqeq') then         
!!!!        target_ => fq
!!!!        target_%name_ = "fq"
!!!!     else if(target_name.eq.'fqfmu'.or.target_name.eq.'fqfmu_pqeq') then         
!!!!        target_ => fqfmu
!!!!        target_%name_ = "fqfmu"
!!!!     else if(target_name.eq.'wfq') then 
!!!!        if(bem%exists) then 
!!!!           target_ => wfq_bem
!!!!           target_%name_ = "wfq_bem"
!!!!        else
!!!!           target_ => wfq
!!!!           target_%name_ = "wfq"
!!!!        endif
!!!!     else if(target_name.eq.'wfqfmu') then 
!!!!        if(bem%exists) then 
!!!!           target_ => wfqfmu_bem
!!!!           target_%name_ = "wfqfmu_bem"
!!!!        else
!!!!           target_ => wfqfmu
!!!!           target_%name_ = "wfqfmu"
!!!!        endif
!!!!     else if(target_name.eq.'wfqib') then 
!!!!        if(bem%exists) then 
!!!!           target_ => wfq_bem
!!!!           target_%name_ = "wfqib_bem"
!!!!        else
!!!!           target_ => wfq
!!!!           target_%name_ = "wfqib"
!!!!        endif
!!!!     else if (target_name.eq."") then
!!!!        if(bem%exists) then
!!!!           target_ => bem
!!!!           target_%name_ = "bem"
!!!!        else 
!!!!           call out_%error( "No target model defined")
!!!!        endif
!!!!     endif
!!!!!
!!!!!
!!!!   end subroutine get_target
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_section_forcefield(inp_,nlines)
!!!!!
!!!!!   get section forcefield
!!!!!
!!!!    implicit none
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    integer :: i
!!!!    character(len=200) :: line
!!!!    character(len=200) :: line_what
!!!!    character(len=200) :: line_keyword
!!!!!
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"forcefield"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) then 
!!!!       if(target_%name_ .eq. "bem") then
!!!!          return
!!!!       else
!!!!          call out_%error( 'No forcefield section in input')
!!!!       endif
!!!!    endif
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end forcefield",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) call out_%error( 'No end for forcefield section in input')
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"forcefield"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!    do i = 1, num_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!!    
!!!!      if(len_trim(line).eq.0) then 
!!!!!         
!!!!         call out_%error( "There is a blank line (forcefield)")
!!!!!         
!!!!      else if(.not.commented_line(line)) then
!!!!!    
!!!!        call LoWeR(line)
!!!!!      
!!!!        if(index(line,":").gt.0) then ! : in input
!!!!!    
!!!!          line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!        
!!!!          line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!    
!!!!          if(line_what.eq.'static') then 
!!!!             target_%ForceField  = trim(adjustl(line_keyword))
!!!!             if(target_%ForceField.ne.'fq'.and.      &
!!!!                target_%ForceField.ne.'fq_pqeq'.and. &
!!!!                target_%ForceField.ne.'fqfmu'.and.   &
!!!!                target_%ForceField.ne.'fqfmu_pqeq')  &
!!!!                call out_%error( "ForceField not recognised")
!!!!          endif
!!!!!          
!!!!          if(line_what.eq.'dynamic') then 
!!!!             if(line_keyword.ne.'wfq'.and.    & 
!!!!                line_keyword.eq.'wfqfmu'.and. &
!!!!                line_keyword.eq.'wfqib') &
!!!!                call out_%error( "wForceField = "//trim(line_keyword)//" not recognized")
!!!!          endif
!!!!!          
!!!!          if(line_what.eq.'kernel') then
!!!!             target_%kernel = line_keyword
!!!!             if(target_%kernel.ne.'coulomb'.and. &
!!!!                target_%kernel.ne.'ohno'.and.    &
!!!!                target_%kernel.ne.'gaussian')    &
!!!!                call out_%error( "kernel = "//trim(line_keyword)//" not recognized")
!!!!          endif
!!!!!
!!!!          if(line_what.ne.'static'.and.&
!!!!             line_what.ne.'dynamic'.and.&
!!!!             line_what.ne.'kernel') then
!!!!             write(out_%iunit,'(1x,a/)') 'Keyword : "'//trim(line_what)//'" in ForceField section not recognized'
!!!!             write(out_%iunit,'(1x,a)' ) 'Possible Keywords in ForceField section:'
!!!!             write(out_%iunit,'(1x,a)' ) '  - static:  [fq,fqfmu,fq_pqeq,fqfmu_pqeq]'
!!!!             write(out_%iunit,'(1x,a)' ) '  - dynamic: [wfq,wfqfmu,wfqib]'
!!!!             write(out_%iunit,'(1x,a)' ) '  - kernel:  [coulomb,ohno,gaussian]'
!!!!             stop
!!!!          endif
!!!!!    
!!!!        else ! : in input
!!!!!      
!!!!          call out_%error( "I am not able to find : in keyword (forcefield)")
!!!!!        
!!!!        endif
!!!!!
!!!!      endif
!!!!!    
!!!!    end do
!!!!!    
!!!!!
!!!!!
!!!!  end subroutine get_section_forcefield
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_section_algorithm(inp_,nlines)
!!!!!
!!!!!   get section algorithm
!!!!!
!!!!    implicit none
!!!!    integer, parameter :: dp=kind(1.0d0)
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: i
!!!!    integer :: ierr
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    character(len=200) :: what
!!!!    character(len=200) :: line
!!!!    character(len=200) :: line_what
!!!!    character(len=200) :: line_keyword
!!!!!
!!!!    what = 'inv'
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"algorithm",found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) return
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end algorithm",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) then
!!!!       call out_%error('No end for algorithm section in input')
!!!!    endif
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"algorithm"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!!
!!!!    do i = 1, num_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!!    
!!!!      if(len_trim(line).eq.0) then 
!!!!!         
!!!!         call out_%error("There is a blank line (algorithm)")
!!!!!         
!!!!      else if(.not.commented_line(line)) then
!!!!!    
!!!!         call LoWeR(line)
!!!!!       
!!!!         if(index(line,":").gt.0) then ! : in input
!!!!!
!!!!           line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!         
!!!!           line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!         
!!!!           if(line_what.eq.'n_iter') then
!!!!              read(line_keyword,'(i5)',iostat=ierr) algorithm%n_iter
!!!!              if(ierr.ne.0) call out_%error("n_iter: "//trim(line_keyword)//&
!!!!                                            " is not an integer")
!!!!           else if(line_what.eq.'n_diis') then
!!!!              read(line_keyword,'(i5)',iostat=ierr) algorithm%n_diis
!!!!              if(ierr.ne.0) call out_%error("n_diis: "//trim(line_keyword)//&
!!!!                                            " is not an integer")
!!!!           else if(line_what.eq.'iter_norm') then
!!!!              read(line_keyword,'(i5)',iostat=ierr) algorithm%norm_for_convergence
!!!!              if(ierr.ne.0) call out_%error("iter_norm: "//trim(line_keyword)//&
!!!!                                            " is not an integer")
!!!!           else if(line_what.eq.'tolerance') then
!!!!              call check_float(out_%iunit,line_what,line_keyword)
!!!!              read(line_keyword,'(f25.16)') algorithm%threshold
!!!!           else 
!!!!              what = line_what
!!!!           endif
!!!!!
!!!!         else ! no : so only line_what
!!!!!
!!!!           line_what = trim(adjustl(line))
!!!!!
!!!!           if(line_what.eq.'inversion'.or.line_what.eq.'inv') then 
!!!!               algorithm%inversion = .true.
!!!!           else if(line_what.eq.'iterative'.or.line_what.eq.'iter') then 
!!!!               algorithm%iterative = .true.
!!!!           else if(line_what.eq.'onthefly'.or.line_what.eq.'otf') then 
!!!!               algorithm%iterative = .true.
!!!!               algorithm%on_the_fly = .true.
!!!!           else if(line_what.eq.'no change param'.or. &
!!!!               line_what.eq.'no change parameters') then
!!!!               algorithm%do_not_change_thresholds = .True.
!!!!           else if(line_what.eq.'matrix in parallel') then 
!!!!               algorithm%matrix_in_parallel  = .True.
!!!!           else if(line_what.eq.'freq in parallel') then
!!!!               algorithm%freq_in_parallel = .True.
!!!!           else if(line_what.eq.'rmse') then
!!!!               algorithm%rmse = .True.
!!!!           else 
!!!!               what = line_what
!!!!           endif
!!!!!
!!!!         endif
!!!!!
!!!!      endif
!!!!!    
!!!!    end do
!!!!!
!!!!!   check keywords
!!!!!
!!!!    if(what.ne.'inversion'.and.&
!!!!       what.ne.'inv'.and.&
!!!!       what.ne.'iterative'.and.&
!!!!       what.ne.'iter'.and.&
!!!!       what.ne.'onthefly'.and.&
!!!!       what.ne.'otf'.and.&
!!!!       what.ne.'no change parameters'.and.&
!!!!       what.ne.'matrix in parallel'.and.&
!!!!       what.ne.'freq in parallel'.and.&
!!!!       what.ne.'n_iter'.and.&
!!!!       what.ne.'n_diis'.and.&
!!!!       what.ne.'iter_norm'.and.&
!!!!       what.ne.'tolerance') then
!!!!       write(out_%iunit,'(1x,a/)') 'Keyword : "'//trim(what)//'" in Algorithm section not recognized'
!!!!       write(out_%iunit,'(1x,a)' ) 'Possible Keywords in Algorithm section:'
!!!!       write(out_%iunit,'(1x,a)' ) '  - inversion(inv)  '
!!!!       write(out_%iunit,'(1x,a)' ) '  - iterative(iter) '
!!!!       write(out_%iunit,'(1x,a)' ) '  - onthefly (otf)  '
!!!!       write(out_%iunit,'(1x,a)' ) '  - no change parameters '
!!!!       write(out_%iunit,'(1x,a)' ) '  - matrix in parallel '
!!!!       write(out_%iunit,'(1x,a)' ) '  - freq in parallel  '
!!!!       write(out_%iunit,'(1x,a)' ) '  - n_iter:               [integer]'
!!!!       write(out_%iunit,'(1x,a)' ) '  - n_diis:               [integer]'
!!!!       write(out_%iunit,'(1x,a)' ) '  - iter_norm:            [integer]'
!!!!       write(out_%iunit,'(1x,a)' ) '  - tolerance:            [ float ]'
!!!!       stop
!!!!    endif
!!!!!
!!!!!
!!!!  end subroutine get_section_algorithm
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_section_field(inp_,nlines)
!!!!!
!!!!!   get section field
!!!!!
!!!!    implicit none
!!!!    integer, parameter :: dp=kind(1.0d0)
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: nostep = 0
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    integer :: i, j, ierr
!!!!    integer :: n_read_freq
!!!!    character(len=200) :: what
!!!!    character(len=200) :: line
!!!!    character(len=200) :: line_what
!!!!    character(len=200) :: line_keyword
!!!!    real(dp) :: freqeV
!!!!    real(dp), parameter:: one = 1.0d0
!!!!!
!!!!    what = 'efield'
!!!!!
!!!!    n_read_freq = 0
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"field",found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) return
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end field",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) call out_%error('No end for field section in input')
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"field"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!    do i = 1, num_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!!    
!!!!      if(len_trim(line).eq.0) then 
!!!!!         
!!!!         call out_%error("There is a blank line (field)")
!!!!!         
!!!!      else if (.not.commented_line(line)) then
!!!!!    
!!!!         call lower(line)
!!!!!       
!!!!         if(index(line,":").gt.0) then ! : in input
!!!!!
!!!!           line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!         
!!!!           line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!         
!!!!           if(line_what.eq.'efield') then 
!!!!              if(line_keyword.eq.'static')  field_%static  = .true.
!!!!              if(line_keyword.eq.'dynamic') field_%dynamic = .true.
!!!!           else if(line_what.eq.'rhs') then 
!!!!              field_%rhs_form = trim(line_keyword)
!!!!              if(field_%rhs_form.ne.'potential'.and.field_%rhs_form.ne.'field') &
!!!!                 call out_%error('rhs : '//trim(field_%rhs_form)//' not recognised in field section')
!!!!           else if(line_what.eq.'polarization') then 
!!!!              field_%polarization = trim(line_keyword)
!!!!              if(field_%polarization.ne.'x'.and. &
!!!!                 field_%polarization.ne.'y'.and. &
!!!!                 field_%polarization.ne.'z'.and. &
!!!!                 field_%polarization.ne.'all') &
!!!!                 call out_%error('polarization : '//trim(field_%polarization)//' not recognised in field section')
!!!!           else if(line_what.eq.'field intensity') then
!!!!              call check_float(out_%iunit,line_what,line_keyword)
!!!!              read(line_keyword,'(f25.16)') field_%e_0
!!!!           else if(line_what.eq.'nfreq') then
!!!!              read(line_keyword,'(i5)',iostat=ierr) field_%n_freq
!!!!              if(ierr.ne.0) call out_%error("nfreq : "//trim(line_keyword)//&
!!!!                                            " is not an integer")
!!!!           else if(line_what.eq.'min freq') then
!!!!              call check_float(out_%iunit,line_what,line_keyword)
!!!!              read(line_keyword,'(f25.16)') field_%min_freq
!!!!              if(field_%min_freq.lt.0.0d0) call out_%error("min freq is less than 0.0d0")
!!!!           else if(line_what.eq.'max freq') then
!!!!              call check_float(out_%iunit,line_what,line_keyword)
!!!!              read(line_keyword,'(f25.16)') field_%max_freq
!!!!              if(field_%max_freq.lt.0.0d0) call out_%error("max freq is less or equal to 0.0d0")
!!!!           else if(line_what.eq.'step') then
!!!!              call check_float(out_%iunit,line_what,line_keyword)
!!!!              read(line_keyword,'(f25.16)') field_%step_freq
!!!!              if(field_%step_freq.lt.0.0d0) call out_%error("step is less or equal to 0.0d0")
!!!!           else if(line_what.eq.'external freq') then
!!!!              n_read_freq = get_n_in_line(line_keyword)
!!!!           else
!!!!              what = line_what
!!!!           endif
!!!!!
!!!!         endif
!!!!!
!!!!      endif
!!!!!    
!!!!    end do
!!!!!
!!!!!   check keywords
!!!!!
!!!!    if(what.ne.'efield'.and.&
!!!!       what.ne.'rhs'.and.&
!!!!       what.ne.'field intensity'.and.&
!!!!       what.ne.'nfreq'.and.&
!!!!       what.ne.'min freq'.and.&
!!!!       what.ne.'max freq'.and.&
!!!!       what.ne.'step'.and.&
!!!!       what.ne.'external freq') then
!!!!       write(out_%iunit,'(1x,a/)') 'Keyword : "'//trim(what)//'" in Field section not recognized'
!!!!       write(out_%iunit,'(1x,a)' ) 'Possible Keywords in Field section:'
!!!!       write(out_%iunit,'(1x,a)' ) '  - efield:      [static,dynamic] '
!!!!       write(out_%iunit,'(1x,a)' ) '  - rhs:         [potential,field]'
!!!!       write(out_%iunit,'(1x,a)' ) '  - field intensity:     [float]'
!!!!       write(out_%iunit,'(1x,a)' ) '  - nfreq:              [integer] '
!!!!       write(out_%iunit,'(1x,a)' ) '  - min freq:            [float]'
!!!!       write(out_%iunit,'(1x,a)' ) '  - max freq:            [float]  '
!!!!       write(out_%iunit,'(1x,a)' ) '  - step:                [float]'
!!!!       write(out_%iunit,'(1x,a)' ) '  - external freq:      [float(s)]'
!!!!       stop
!!!!    endif
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"field"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!!
!!!!    If(field_%dynamic) then
!!!!!
!!!!!      check input frequencies
!!!!!
!!!!       call field_%check_frequencies(n_read_freq)
!!!!!
!!!!       allocate(field_%freq(field_%n_freq),stat=ierr)
!!!!       If(ierr.gt.0) call out_%error('Insufficient space for Freq allocation')
!!!!       Call array_clear(field_%n_freq,field_%freq)
!!!!!
!!!!       do i = 1, num_keywords
!!!!!
!!!!         read(inp_%iunit,'(a)') line
!!!!!    
!!!!         if(len_trim(line).eq.0) then 
!!!!!         
!!!!            call out_%error("There is a blank line (field)")
!!!!!         
!!!!         else if (.not.commented_line(line)) then
!!!!!    
!!!!            call LoWeR(line)
!!!!!       
!!!!            if(index(line,":").gt.0) then ! : in input
!!!!!
!!!!              line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!         
!!!!              line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!      
!!!!              if(line_what.eq.'external freq') then
!!!!!      
!!!!                 read(line_keyword,*) field_%freq(1:field_%n_freq)
!!!!!                 
!!!!                 do j = 1, field_%n_freq
!!!!                    if(field_%freq(j).lt.0.0d0) then
!!!!                       call out_%error("One or more external freq are less than zero")
!!!!                    else if(field_%freq(j).eq.0.0d0) then 
!!!!                       field_%freq(j) = 1.0d+10
!!!!                    endif
!!!!                 enddo
!!!!                 nostep = 1
!!!!!
!!!!              endif
!!!!!
!!!!            endif
!!!!!
!!!!         endif
!!!!!    
!!!!       end do
!!!!!
!!!!       call field_%set_frequencies(nostep)
!!!!!
!!!!    else if (field_%static.and. &
!!!!            (target_%name_.eq.'wfq'.or. &
!!!!             target_%name_.eq.'wfqfmu'.or.&
!!!!             target_%name_.eq.'wfqib')) then
!!!!       field_%n_freq = 1
!!!!       field_%dynamic = .true.
!!!!       Allocate(field_%freq(field_%n_freq))
!!!!       FreqeV  = 1.0d-6
!!!!       call freqeVtonm(FreqeV,field_%freq(1))
!!!!!
!!!!    endif
!!!!!
!!!!!
!!!!  end subroutine get_section_field
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_section_output(inp_,nlines)
!!!!!
!!!!!   get section output
!!!!!
!!!!    implicit none
!!!!    integer, parameter :: dp=kind(1.0d0)
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: ierr
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    integer :: i
!!!!    character(len=200) :: what
!!!!    character(len=200) :: line
!!!!    character(len=200) :: line_what
!!!!    character(len=200) :: line_keyword
!!!!    real(dp), parameter:: one = 1.0d0
!!!!!
!!!!    what = 'verbose'
!!!!!
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"output",found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) return
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end output",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) call out_%error('No end for output section in input')
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"output"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!!
!!!!    do i = 1, num_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!!    
!!!!      if(len_trim(line).eq.0) then 
!!!!!         
!!!!         call out_%error("There is a blank line (output)")
!!!!!         
!!!!      else if(.not.commented_line(line)) then
!!!!!    
!!!!         call LoWeR(line)
!!!!!       
!!!!         if(index(line,":").gt.0) then ! : in input
!!!!!
!!!!           line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!         
!!!!           line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!         
!!!!           if(line_what.eq.'maxima') then 
!!!!              if(line_keyword.eq.'iso alpha')   general%i_solution = 1
!!!!              if(line_keyword.eq.'long alpha')  general%i_solution = 2
!!!!              if(line_keyword.eq.'absorption')  general%i_solution = 3
!!!!              if(line_keyword.eq.'scattering')  general%i_solution = 4
!!!!              if(line_keyword.eq.'exctinction') general%i_solution = 5
!!!!           else if(line_what.eq.'verbose') then
!!!!              read(line_keyword,'(i5)',iostat=ierr) out_%ivrb
!!!!              if(ierr.ne.0) call out_%error("verbose : "//trim(line_keyword)//&
!!!!                                            " is not an integer")
!!!!           else 
!!!!              what = line_what
!!!!           endif
!!!!!
!!!!         else
!!!!!
!!!!           line_what = trim(adjustl(line))
!!!!!
!!!!           what = line_what
!!!!!
!!!!         endif
!!!!!
!!!!      endif
!!!!!    
!!!!    end do
!!!!!    
!!!!    if(what.ne.'verbose'.and.&
!!!!       what.ne.'lorentzian') then
!!!!       write(out_%iunit,'(1x,a/)') 'Keyword : "'//trim(what)//'" in Output section not recognized'
!!!!       write(out_%iunit,'(1x,a)' ) 'Possible Keywords in Output section:'
!!!!       write(out_%iunit,'(1x,a)' ) '  - maxima :   [iso alpha,long alpha,absorption,scattering,extinction]'
!!!!       write(out_%iunit,'(1x,a)' ) '  - verbose:      [integer]'
!!!!       stop
!!!!    endif
!!!!!
!!!!!
!!!!  end subroutine get_section_output
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_section_parameters(inp_,nlines)
!!!!!
!!!!!   get section parameters
!!!!!
!!!!    implicit none
!!!!    integer, parameter :: dp=kind(1.0d0)
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!    
!!!!!   internal variables
!!!!!
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    integer :: i, j
!!!!    integer :: ierr
!!!!    character(len=200) :: line
!!!!    character(len=200) :: line_what
!!!!    character(len=200) :: line_keyword
!!!!    real(dp), parameter:: one = 1.0d0
!!!!!
!!!!    integer, dimension(:), allocatable :: num_keywords_atomtype
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"parameters",found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) return
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end parameters",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) call out_%error('No end for parameters section in input')
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!!   get number of keywords per n_atomtypes
!!!!!
!!!!    allocate(num_keywords_atomtype(target_%n_atomtypes))
!!!!!
!!!!!   allocation of vectors wfq/wfqfmu
!!!!!
!!!!    if(target_%name_.eq.'wfq'          .or. &
!!!!       target_%name_.eq.'wfqfmu'       .or. &
!!!!       target_%name_.eq.'wfqib'        .or. &
!!!!       target_%name_.eq.'wfq_bem'      .or. &
!!!!       target_%name_.eq.'wfqfmu_bem'   .or. &
!!!!       target_%name_.eq.'wfqib_bem') then
!!!!!    
!!!!       Allocate(parameters%tau           (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for tau')
!!!!       parameters%tau = zero
!!!!!
!!!!       Allocate(parameters%sigma_0       (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for sigma_0')
!!!!       parameters%sigma_0 = zero
!!!!!
!!!!       Allocate(parameters%scaling       (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for scaling')
!!!!       parameters%scaling = one ! scaling default
!!!!!
!!!!       Allocate(parameters%A_ij          (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for A_ij')
!!!!       parameters%A_ij = zero
!!!!!
!!!!       Allocate(parameters%fermi_d       (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for fermi_d')
!!!!       parameters%fermi_d = zero
!!!!!
!!!!       Allocate(parameters%fermi_s       (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for fermi_s')
!!!!       parameters%fermi_s = zero
!!!!!
!!!!       Allocate(parameters%first_length  (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for first_length')
!!!!       parameters%first_length = zero
!!!!!
!!!!       Allocate(parameters%second_length (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for second_length')
!!!!       parameters%second_length = zero
!!!!!
!!!!       Allocate(parameters%third_length  (target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for third length')
!!!!       parameters%third_length = zero
!!!!!
!!!!       Allocate(parameters%factor_density(target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for factor_density')
!!!!       parameters%factor_density = zero
!!!!!
!!!!       Allocate(parameters%n_structures(target_%n_atomtypes), stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for n_structures')
!!!!       parameters%n_structures = 1 !n_structures default
!!!!!
!!!!       Allocate(parameters%graphene_geometry(target_%n_atomtypes),stat=ierr)
!!!!       if(ierr.gt.0) call out_%error('not enough space for n_structures')
!!!!!
!!!!       if(target_%name_.eq.'wfqfmu'.or.target_%name_.eq.'wfqfmu_bem') then
!!!!!
!!!!          allocate(wfqfmu%epsilon_w(target_%n_atomtypes))
!!!!          wfqfmu%epsilon_w = ""
!!!!!
!!!!          allocate(wfqfmu%file_in(target_%n_atomtypes))
!!!!          wfqfmu%file_in = ""
!!!!!
!!!!       endif
!!!!!
!!!!    endif
!!!!!
!!!!    do i = 1 , target_%n_atomtypes
!!!!!
!!!!       rewind(inp_%iunit)
!!!!       call go_to_string(inp_%iunit,"atomtype "//trim(target_%atom_type(i)),found_string,nlines,num_string_initial)
!!!!!
!!!!       if(.not.found_string) call out_%error("Not found atomtype "//trim(target_%atom_type(i))//" in parameters section")
!!!!!
!!!!       rewind(inp_%iunit)
!!!!       call go_to_string(inp_%iunit,"end atomtype "//trim(target_%atom_type(i)),found_string,nlines,num_string_end)
!!!!!
!!!!       if(.not.found_string) call out_%error('No end for atomtype '//trim(target_%atom_type(i))//'in parameters section')
!!!!!
!!!!       num_keywords_atomtype(i) = num_string_end - num_string_initial - 1
!!!!!
!!!!    enddo
!!!!!
!!!!    do i = 1, target_%n_atomtypes
!!!!!
!!!!       rewind(inp_%iunit)
!!!!       call go_to_string(inp_%iunit,"atomtype "//trim(target_%atom_type(i)),found_string,nlines,num_string_initial)
!!!!!
!!!!       do j = 1, num_keywords_atomtype(i)
!!!!!
!!!!         read(inp_%iunit,'(a)') line
!!!!!       
!!!!         if(len_trim(line).eq.0) then 
!!!!!            
!!!!            call out_%error("There is a blank line (parameters)")
!!!!!            
!!!!         else if(.not.commented_line(line)) then
!!!!!       
!!!!            call lower(line)
!!!!!          
!!!!            if(index(line,":").gt.0) then ! : in input
!!!!!
!!!!              line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!            
!!!!              line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!            
!!!!              if(line_what.eq.'graphene geometry') then 
!!!!!              
!!!!                 parameters%graphene_geometry(i) = line_keyword
!!!!!                 
!!!!                 if(line_keyword.ne.'disk'.and. & 
!!!!                    line_keyword.ne.'triang'.and. &
!!!!                    line_keyword.ne.'hole'.and. &
!!!!                    line_keyword.ne.'ribbon'.and. &
!!!!                    line_keyword.ne.'nanotube'.and. &
!!!!                    line_keyword.ne.'cnt'.and. &
!!!!                    line_keyword.ne.'ring'.and. &
!!!!                    line_keyword.ne.'bowtie-cnt') then
!!!!                    Write(out_%iunit,'(1x,a)')  'graphene geometry: '//trim(line_keyword)//' not recognised'
!!!!                    write(out_%iunit,'(1x,a)' ) 'Possible geometries:'
!!!!                    write(out_%iunit,'(1x,a)' ) '  - disk'
!!!!                    write(out_%iunit,'(1x,a)' ) '  - triang'
!!!!                    write(out_%iunit,'(1x,a)' ) '  - hole'
!!!!                    write(out_%iunit,'(1x,a)' ) '  - ribbon'
!!!!                    write(out_%iunit,'(1x,a)' ) '  - nanotube (cnt)'
!!!!                    write(out_%iunit,'(1x,a)' ) '  - ring'
!!!!                    write(out_%iunit,'(1x,a)' ) '  - bowtie-cnt'
!!!!                    stop
!!!!!                    
!!!!                 endif
!!!!!            
!!!!              else if(line_what.eq.'epsilon') then 
!!!!                 wfqfmu%epsilon_w(i) = line_keyword
!!!!                 wfqfmu%alpha_inside = .true.
!!!!!            
!!!!              else if(line_what.eq.'wfqfmu file') then 
!!!!                 wfqfmu%file_in(i) = line_keyword
!!!!                 wfqfmu%alpha_inside = .false.
!!!!!                 
!!!!              else if(line_what.eq.'tau')         then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%tau(i)
!!!!!            
!!!!              else if(line_what.eq.'sigma0')           then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%sigma_0(i)
!!!!                  parameters%sigma_0(i) = parameters%sigma_0(i)*FSmAu
!!!!!            
!!!!              else if(line_what.eq.'scaling')         then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%scaling(i)
!!!!!            
!!!!              else if(line_what.eq.'ri')               then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%A_ij(i)
!!!!                  parameters%A_ij(i) = parameters%A_ij(i)**2
!!!!!            
!!!!              else if(line_what.eq.'fermi function d') then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%fermi_d(i)
!!!!!            
!!!!              else if(line_what.eq.'fermi function s') then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%fermi_s(i)
!!!!!            
!!!!              else if(line_what.eq.'primary length')   then 
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%first_length(i)
!!!!                  parameters%first_length(i) = parameters%first_length(i)*ToBohr
!!!!!            
!!!!              else if(line_what.eq.'secondary length') then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%second_length(i)
!!!!                  parameters%second_length(i) = parameters%second_length(i)*ToBohr
!!!!!            
!!!!              else if(line_what.eq.'third length') then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%third_length(i)
!!!!                  parameters%third_length(i) = parameters%third_length(i)*ToBohr
!!!!!            
!!!!              else if(line_what.eq.'factor density')   then
!!!!                  call check_float(out_%iunit,line_what,line_keyword)
!!!!                  read(line_keyword,'(f25.16)') parameters%factor_density(i)
!!!!!            
!!!!              else if(line_what.eq.'n structures') then
!!!!                  read(line_keyword,'(i10)',iostat=ierr) parameters%n_structures(i)
!!!!                  if(ierr.ne.0) then
!!!!                     call out_%error("n structures: "//trim(line_keyword)//&
!!!!                                     " is not an integer")
!!!!                  else 
!!!!                     if(parameters%n_structures(i).le.0) call out_%error("n structures must be positive")
!!!!                  endif
!!!!!            
!!!!              else 
!!!!!            
!!!!                  write(out_%iunit,'(1x,a/)') 'Keyword : "'//trim(line_what)//'" in Parameters section not recognized'
!!!!                  write(out_%iunit,'(1x,a)' ) 'Possible Keywords in Parameters section:'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - factor redshift:   [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - scaling:           [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - tau:               [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - sigma0:            [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - RI:                [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - Fermi function d:  [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - Fermi function s:  [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - Factor density:    [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - Primary Length:    [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - Secondary Length:  [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - Third Length:      [float]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - n structures:      [integer]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - graphene geometry: [disk,triang,ring,ribbon,nanotube(cnt),ring,bowtie cnt]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - wfqfmu file:       [string]'
!!!!                  write(out_%iunit,'(1x,a)' ) '  - epsilon:           [string]'
!!!!                  stop
!!!!!            
!!!!              endif
!!!!!
!!!!            endif
!!!!!
!!!!         endif
!!!!!       
!!!!       end do !keyword per atomtype
!!!!!       
!!!!    end do !target_%n_atomtypes
!!!!!
!!!!!
!!!!  end subroutine get_section_parameters
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_section_atomtypes(inp_,nlines)
!!!!!
!!!!!   get section atomtypes
!!!!!
!!!!    implicit none
!!!!    integer, parameter :: dp=kind(1.0d0)
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: ierr
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    integer :: i, j, k
!!!!    character(len=200) :: line
!!!!    character(len=200) :: line_what
!!!!    character(len=200) :: line_keyword
!!!!!
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"atom types" ,found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) call out_%error('No atom types section in input')
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end atom types",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) call out_%error('No end for atom type section in input')
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"atom types",found_string,nlines,num_string_initial)
!!!!!
!!!!!
!!!!    j = 0
!!!!    k = 0
!!!!    do i = 1, num_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!!    
!!!!      if(len_trim(line).eq.0) then 
!!!!!         
!!!!         call out_%error("There is a blank line (atom types)")
!!!!!         
!!!!      else if(.not.commented_line(line)) then
!!!!!    
!!!!         call LoWeR(line)
!!!!         j = j + 1 
!!!!!       
!!!!         if(index(line,":").gt.0) then ! : in input
!!!!!
!!!!           line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!         
!!!!           line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!         
!!!!           if(j.eq.1.and.line_what.ne.'number') &
!!!!              call out_%error("First keyword in Atomtype should be number: [integer]")
!!!!!         
!!!!           if(line_what.eq.'number') then
!!!!              read(line_keyword,'(i5)',iostat=ierr) target_%n_atomtypes
!!!!              if(ierr.ne.0) call out_%error("number : "//trim(line_keyword)//&
!!!!                                            " is not an integer")
!!!!              if(target_%n_atomtypes.gt.num_keywords-1) &
!!!!                 call out_%error("more atom types than keywords")
!!!!!
!!!!              allocate(target_%atom_type(target_%n_atomtypes),stat=ierr)
!!!!              if(ierr.ne.0) call out_%error('Not Enough Space for atomtype')
!!!!!
!!!!              allocate(fq%chi(target_%n_atomtypes),stat=ierr)
!!!!              if(ierr.ne.0) call out_%error('Not Enough Space for electronegativity')
!!!!!
!!!!              allocate(fq%eta(target_%n_atomtypes),stat=ierr)
!!!!              if(ierr.ne.0) call out_%error('Not Enough Space for chemical hardnesses')
!!!!!
!!!!              if(target_%forcefield.eq.'fq_pqeq') then
!!!!!
!!!!                 allocate(fq%r_q(target_%n_atomtypes),Stat=IErr)
!!!!                 if(ierr.ne.0) call out_%error('Not Enough Space for R_q (FQ)')
!!!!!
!!!!              else if(target_%forcefield.eq.'fqfmu'.or.target_%forcefield.eq.'fqfmu_pqeq') then 
!!!!!
!!!!                 allocate(fqfmu%alpha(target_%n_atomtypes),Stat=IErr)
!!!!                 if(ierr.ne.0) call out_%error('Not Enough Space for polarizabilities (FQFMu)')
!!!!                 allocate(fq%r_q(target_%n_atomtypes),Stat=IErr)
!!!!                 if(ierr.ne.0) call out_%error('Not Enough Space for R_q (FQFMu)')
!!!!                 allocate(fqfmu%r_mu(target_%n_atomtypes),Stat=IErr)
!!!!                 if(ierr.ne.0) call out_%error('Not Enough Space for R_mu (FQFMu)')
!!!!!
!!!!              else if(target_%forcefield.eq.'fq'.and.target_%kernel.eq.'gaussian') then
!!!!!
!!!!                 allocate(fq%r_q(target_%n_atomtypes),Stat=IErr)
!!!!                 if(ierr.ne.0) call out_%error('Not Enough Space for R_q (FQ)')
!!!!!
!!!!              endif
!!!!!
!!!!           else
!!!!              k = k + 1
!!!!              call read_atomtype_parameters(line_what,    &
!!!!                                            line_keyword, &
!!!!                                            k)
!!!!!
!!!!           endif
!!!!           
!!!!!
!!!!         endif
!!!!!
!!!!      endif
!!!!!    
!!!!    end do
!!!!!
!!!!!
!!!!  end subroutine get_section_atomtypes
!!!!!-----------------------------------------------------------------------
!!!!  subroutine read_atomtype_parameters(line_what,    &
!!!!                                      line_keyword, &
!!!!                                      k)
!!!!!
!!!!!   read atomtype parameters
!!!!!
!!!!    implicit none
!!!!!
!!!!!   input variables
!!!!!
!!!!    character(len=200), intent(in)  :: line_what
!!!!    character(len=200), intent(in)  :: line_keyword
!!!!    integer :: k
!!!!!
!!!!    integer :: init,end_, i_name
!!!!    character(len=200), dimension(:), allocatable :: param_type, param
!!!!    character(len=200) :: string, substring 
!!!!!
!!!!    integer :: i,count_parameters
!!!!!
!!!!!
!!!!    target_%atom_type(k) = line_what
!!!!!
!!!!    count_parameters = 0
!!!!    do i = 1, 200
!!!!       if(line_keyword(i:i).eq.'=') then
!!!!          count_parameters=count_parameters+1
!!!!       endif
!!!!    enddo
!!!!    if(count_parameters.eq.0) call out_%error('atom type section, but no parameters')
!!!!!
!!!!    Allocate(param_type(count_parameters))
!!!!    Allocate(param(count_parameters))
!!!!!
!!!!    string = line_keyword
!!!!!
!!!!    init = move_to_character(string,'[')+1
!!!!    end_ = move_to_character(string,']')-1
!!!!!
!!!!    string = string(init:end_)
!!!!    string = sweep_blanks(string)
!!!!!
!!!!    do i = 1, count_parameters
!!!!!
!!!!       if(i.eq.1) then
!!!!          init = 1
!!!!       else 
!!!!          init = move_to_character(string,',')+1
!!!!       endif
!!!!!
!!!!       string = string(init:200)
!!!!!
!!!!       if(i.eq.count_parameters) then
!!!!          end_ = 200
!!!!       else
!!!!          end_ = move_to_character(string,',')-1
!!!!       endif
!!!!!
!!!!       substring =  string(1:end_)
!!!!!
!!!!       i_name = move_to_character(string,'=')
!!!!       param_type(i)   = trim(substring(1:i_name-1))
!!!!       param     (i)   = trim(substring(i_name+1:) )
!!!!!
!!!!       string    =  string(end_+1:200)
!!!!!
!!!!    enddo
!!!!!
!!!!    do i = 1, count_parameters
!!!!!
!!!!       if(param_type(i).eq.'chi') then
!!!!          call check_float(out_%iunit,param_type(i),param(i))
!!!!          read(param(i),'(f25.16)') fq%chi(k)
!!!!       else if(param_type(i).eq.'eta') then
!!!!          call check_float(out_%iunit,param_type(i),param(i))
!!!!          read(param(i),'(f25.16)') fq%eta(k)
!!!!          if(fq%eta(k).lt.1.0d-14) &
!!!!             call out_%error("Eta = 0.0d0 for atomtype: "// trim(target_%atom_type(k)))
!!!!       else if(param_type(i).eq.'alpha') then
!!!!          call check_float(out_%iunit,param_type(i),param(i))
!!!!          read(param(i),'(f25.16)') fqfmu%alpha(k)
!!!!          if(fqfmu%alpha(k).lt.1.0d-14) &
!!!!             call out_%error("alpha = 0.0d0 for atomtype: "// trim(target_%atom_type(k)))
!!!!       else if(param_type(i).eq.'rq') then
!!!!          call check_float(out_%iunit,param_type(i),param(i))
!!!!          read(param(i),'(f25.16)') fq%r_q(k)
!!!!          if(fq%r_q(k).lt.1.0d-14) &
!!!!             call out_%error("Rq = 0.0d0 for atomtype: "// trim(target_%atom_type(k)))
!!!!       else if(param_type(i).eq.'rmu') then
!!!!          call check_float(out_%iunit,param_type(i),param(i))
!!!!          read(param(i),'(f25.16)') fqfmu%r_mu(k)
!!!!          if(fqfmu%r_mu(k).lt.1.0d-14) &
!!!!             call out_%error("Rmu = 0.0d0 for atomtype: "// trim(target_%atom_type(k)))
!!!!       else
!!!!          call out_%error('Parameter '//trim(param_type(i))//' not recognized for atom name: '//trim(target_%atom_type(k)))
!!!!       endif
!!!!!
!!!!    enddo
!!!!!
!!!!    if(target_%forcefield.eq.'fq'.and.target_%kernel.eq.'gaussian') then
!!!!       fq%r_q(k)  = dsqrt(two/pi)/fq%eta(k)
!!!!    else if(target_%forcefield.eq.'fqfmu') then
!!!!       if(fqfmu%alpha(k).lt.1.0d-14) &
!!!!          call out_%error("FQFMu ForceField but not Alpha in atomtypes")
!!!!       fq%r_q(k)     = dsqrt(two/pi)/fq%eta(k)
!!!!       fqfmu%r_mu(k) = (fqfmu%alpha(k)/(dsqrt(pi/two)*three))**(one/three)
!!!!    else if(target_%forcefield.eq.'fqfmu_pqeq') then
!!!!       fqfmu%pqeq = .true.
!!!!       if(fqfmu%alpha(k).lt.1.0d-14) &
!!!!          call out_%error("FQFMu_PQEq ForceField but no Alpha in atomtypes")
!!!!       if(fq%r_q(k).lt.1.0d-14) &
!!!!          call out_%error("FQFMu_PQEq ForceField but no Rq in atomtypes")
!!!!       if(fqfmu%r_mu(k).lt.1.0d-14) &
!!!!          call out_%error("FQFMu_PQEq ForceField but no Rmu in atomtypes")
!!!!       fq%r_q(k)     = dsqrt(fq%r_q(k)**2/0.2341d0)
!!!!       fqfmu%r_mu(k) = dsqrt(fqfmu%r_mu(k)**2/0.2341d0)
!!!!    else if(target_%forcefield.eq.'fq_pqeq') then
!!!!       fq%pqeq = .true.
!!!!       if(fq%r_q(k).lt.1.0d-14) &
!!!!          call out_%error("FQ_PQEq ForceField but no Rq in atomtypes")
!!!!       fq%r_q(k) = dsqrt(fq%r_q(k)**2/0.2341d0)
!!!!    endif
!!!!!
!!!!!
!!!!  end subroutine read_atomtype_parameters
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_section_bem(inp_,nlines)
!!!!!
!!!!!   get section general
!!!!!
!!!!    implicit none
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    integer :: i
!!!!    character(len=200) :: line
!!!!    character(len=200) :: line_what
!!!!    character(len=200) :: line_keyword
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"bem",found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) return
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end bem",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) call out_%error('No end for bem section in input')
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"bem"    ,found_string,nlines,num_string_initial)
!!!!    do i = 1, num_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!      if(len_trim(line).eq.0) then
!!!!!
!!!!         call out_%error("There is a blank line (BEM)")
!!!!!
!!!!      else if(.not.commented_line(line)) then
!!!!!
!!!!        call LoWeR(line)
!!!!!
!!!!        if(index(line,":").gt.0) then ! : in input
!!!!!
!!!!           line_what = trim(adjustl(line(1:index(line,":")-1) ))
!!!!!
!!!!           line_keyword = trim(adjustl( line(index(line,":")+1:200) ))
!!!!!
!!!!           if(line_what.eq.'normal scalar factor') then
!!!!!
!!!!              call check_float(out_%iunit,line_what,line_keyword)
!!!!              read(line_keyword,'(f25.16)') bem%normal_factor
!!!!!
!!!!           else if(line_what.eq.'gmsh file') then ! gmsh file with tesserae
!!!!!
!!!!              bem%gmsh_file = line_keyword
!!!!!
!!!!           else if(line_what.eq.'gmsh file units') then ! gmsh file with tesserae
!!!!!
!!!!              if (trim(line_keyword) .eq. 'angstroms' .or. trim(line_keyword) .eq. 'aa') then 
!!!!!
!!!!                 bem%tess_scale = 1.0d0
!!!!!
!!!!              else if (trim(line_keyword) .eq. 'nm') then 
!!!!!
!!!!                 bem%tess_scale = 10.0d0
!!!!!
!!!!              else
!!!!!
!!!!                 call out_%error('Keyword: "'//trim(line_keyword)//'" in gmsh file units (bem) not recognized')
!!!!!
!!!!              endif
!!!!!
!!!!           else if(line_what.eq.'solvent') then !external solvent
!!!!!
!!!!              bem%solvent = line_keyword
!!!!!
!!!!           else if(line_what.eq.'permittivity') then !permittivity function tabulated inside nanoFQ
!!!!!
!!!!              bem%permittivity_inside = .true.
!!!!              bem%permittivity_type   = line_keyword
!!!!!
!!!!           else if(line_what.eq.'permittivity file') then !permittivity function tabulated inside nanoFQ
!!!!!
!!!!              bem%permittivity_inside = .false.
!!!!              bem%permittivity_file   = line_keyword
!!!!!
!!!!           else if(line_what.eq.'green function') then ! -- approx or exact -- mennucci or trugler
!!!!!
!!!!              bem%green_function = line_keyword
!!!!              if(bem%green_function.ne.'approx'.and.&
!!!!                 bem%green_function.ne.'exact') &
!!!!                 call out_%error('Keyword: "'//trim(line_keyword)//'" in green function (bem) not recognized')
!!!!!
!!!!           else if(line_what.eq.'sphere radius') then !radius of eventual sphere if green_function = approx
!!!!!
!!!!              call check_float(out_%iunit,line_what,line_keyword)
!!!!              read(line_keyword,'(f25.16)') bem%sphere_r
!!!!!
!!!!              if(bem%sphere_r.le.0.0d0) call out_%error("Sphere radius: bad value in input")
!!!!              bem%sphere_r = bem%sphere_r * ToBohr
!!!!!
!!!!           endif
!!!!!        
!!!!           if(line_what.ne.'normal scalar factor'  .and. &
!!!!              line_what.ne.'gmsh file'             .and. &
!!!!              line_what.ne.'gmsh file units'       .and. &
!!!!              line_what.ne.'permittivity file'     .and. &
!!!!              line_what.ne.'permittivity'          .and. &
!!!!              line_what.ne.'solvent'               .and. &
!!!!              line_what.ne.'green function'        .and. &
!!!!              line_what.ne.'sphere radius') then
!!!!              write(out_%iunit,'(1x,a/)') 'Keyword: "'//trim(line_what)//'" in bem section not recognized'
!!!!              write(out_%iunit,'(1x,a)' ) 'Possible keywords in bem section:'
!!!!              write(out_%iunit,'(1x,a)' ) '  - normal scalar factor' 
!!!!              write(out_%iunit,'(1x,a)' ) '  - gmsh file'            
!!!!              write(out_%iunit,'(1x,a)' ) '  - permittivity'
!!!!              write(out_%iunit,'(1x,a)' ) '  - permittivity file'
!!!!              write(out_%iunit,'(1x,a)' ) '  - solvent'              
!!!!              write(out_%iunit,'(1x,a)' ) '  - green function'       
!!!!              write(out_%iunit,'(1x,a)' ) '  - sphere radius (in Ang)'
!!!!              stop
!!!!           endif
!!!!!           
!!!!        else ! no : so only line_what
!!!!!
!!!!           line_what = trim(adjustl(line))
!!!!!
!!!!           if(line_what.eq.'charge constraint') then !charge constraint 
!!!!!
!!!!              bem%charge_constraint = .true.
!!!!!
!!!!!lookhere
!!!!           else if(line_what.eq.'ief') then ! IEF for rhs
!!!!!
!!!!              bem%rhs_ief = .true.
!!!!!
!!!!           else if(line_what.eq.'reduced path') then ! Correction to epsilon
!!!!!
!!!!              bem%reduced_path = .true.
!!!!!
!!!!           else if(line_what.eq.'split polar') then
!!!!!
!!!!              bem%split_polar = .true.
!!!!!endhere
!!!!           else
!!!!! 
!!!!              call out_%error("Line what "//trim(line_what)//" not recognised in BEM section")
!!!!!
!!!!           endif
!!!!!        
!!!!        endif
!!!!      endif
!!!!    enddo
!!!!!todo          
!!!!!   subroutine che faccia check delle keyword
!!!!!
!!!!    call bem%read_gmsh_file()
!!!!!
!!!!    call bem%read_permittivity_bem()
!!!!!
!!!!  end subroutine get_section_bem
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_section_general(inp_,nlines)
!!!!!
!!!!!   get section general
!!!!!
!!!!    implicit none
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!
!!!!!   input/output variables
!!!!!
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    integer :: i
!!!!    character(len=200) :: line
!!!!    character(len=200) :: line_what
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"general"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) return
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end general",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) call out_%error('No end for general section in input')
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"general"    ,found_string,nlines,num_string_initial)
!!!!!
!!!!    do i = 1, num_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!!    
!!!!      if(len_trim(line).eq.0) then 
!!!!!         
!!!!         call out_%error("There is a blank line (general)")
!!!!!         
!!!!      else if(.not.commented_line(line)) then
!!!!!    
!!!!        call LoWeR(line)
!!!!!      
!!!!        line_what = trim(adjustl(line))
!!!!!    
!!!!        if(line_what.eq.'principal axes') then 
!!!!           general%principal_axis = .true.
!!!!!    
!!!!        else if(line_what.eq.'save info file') then 
!!!!           general%save_info = .true.
!!!!        else
!!!!           write(out_%iunit,'(1x,a/)') 'Keyword : "'//trim(line_what)//'" in General section not recognized'
!!!!           write(out_%iunit,'(1x,a)' ) 'Possible Keywords in General section:'
!!!!           write(out_%iunit,'(1x,a)' ) '  - principal axes'
!!!!           write(out_%iunit,'(1x,a)' ) '  - save info file'
!!!!           stop
!!!!        endif
!!!!!
!!!!      endif
!!!!!    
!!!!    end do
!!!!!
!!!!!
!!!!  end subroutine get_section_general
!!!!!----------------------------------------------------------------------
!!!!  subroutine get_section_geometry(inp_,nlines,AtomName)
!!!!!
!!!!!   get section geometry
!!!!!
!!!!    implicit none
!!!!    integer, parameter :: dp=kind(1.0d0)
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!!
!!!!!   input/output variables
!!!!!
!!!!    character(len=200), dimension(:), allocatable, intent(inout) :: AtomName
!!!!!
!!!!!    internal variables
!!!!!
!!!!    integer :: num_string_initial, num_string_end, num_keywords
!!!!    logical :: found_string
!!!!    integer :: i, k
!!!!    character(len=200) :: line
!!!!!
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"input geometry" ,found_string,nlines,num_string_initial)
!!!!!
!!!!    if(.not.found_string) call out_%error('No input geometry section in input')
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"end input geometry",found_string,nlines,num_string_end)
!!!!!
!!!!    if(.not.found_string) call out_%error('No end for input geometry section in input')
!!!!!
!!!!    num_keywords = num_string_end - num_string_initial - 1
!!!!!
!!!!    call inp_%get_number_of_atoms(nlines,num_keywords)
!!!!!
!!!!    Allocate(AtomName(target_%n_atoms))
!!!!    Allocate(fq%i_mol(target_%n_atoms))
!!!!    Allocate(target_%coord(3,target_%n_atoms))
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"input geometry",found_string,nlines,num_string_initial)
!!!!!
!!!!    k = 0
!!!!    do i = 1, num_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!!    
!!!!      if(len_trim(line).eq.0) then 
!!!!!         
!!!!         call out_%error("There is a blank line (input geometry)")
!!!!!         
!!!!      else if(.not.commented_line(line)) then
!!!!!    
!!!!         call lower(line)
!!!!         line = substitute_string(line,char(9)," ") !eliminate any tab from string
!!!!!       
!!!!         k = k + 1
!!!!         call read_atoms_imols_xyz(line,atomname(k),fq%i_mol(k), &
!!!!                                   target_%coord(1,k),target_%coord(2,k),target_%coord(3,k))
!!!!!
!!!!!
!!!!      endif
!!!!!    
!!!!    end do
!!!!!
!!!!!
!!!!  end subroutine get_section_geometry
!!!!!-----------------------------------------------------------------------
!!!!  subroutine get_number_of_atoms(inp_,nlines,n_keywords)
!!!!!
!!!!!   get number of atoms
!!!!!
!!!!    implicit none
!!!!!
!!!!!   input variables
!!!!!
!!!!    class(inp_type)  :: inp_
!!!!    integer, intent(in)             :: nlines
!!!!    integer, intent(in)             :: n_keywords
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: num_string
!!!!    integer :: i
!!!!    logical :: found_string
!!!!    character(len=200) :: line
!!!!!
!!!!    rewind(inp_%iunit)
!!!!    call go_to_string(inp_%iunit,"input geometry",found_string,nlines,num_string)
!!!!!
!!!!    target_%n_atoms = 0 
!!!!    do i = 1, n_keywords
!!!!!
!!!!      read(inp_%iunit,'(a)') line
!!!!!    
!!!!      if(len_trim(line).eq.0) then 
!!!!!         
!!!!         call out_%error("There is a blank line (input geometry)")
!!!!!         
!!!!      else if(.not.commented_line(line)) then
!!!!!         
!!!!          target_%n_atoms = target_%n_atoms + 1
!!!!!         
!!!!      endif
!!!!!
!!!!    enddo
!!!!!
!!!!!
!!!!  end subroutine get_number_of_atoms
!!!!!-----------------------------------------------------------------------
!!!!  subroutine read_atoms_imols_xyz(line,atomname,IMol,x,y,z)
!!!!!
!!!!!   read atoms imols xyz
!!!!!
!!!!    implicit none
!!!!    integer, parameter :: dp=kind(1.0d0)
!!!!!
!!!!!   input variables
!!!!!
!!!!    character(len=200), intent(inout) :: line
!!!!!
!!!!!   input/output variables
!!!!!
!!!!    character(len=200), intent(inout)      :: AtomName
!!!!    integer, intent(inout)                 :: IMol
!!!!    real(dp), intent(inout)                :: x,y,z
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: ierr
!!!!    integer :: istart, iend
!!!!    real(dp), parameter :: zero = 0.0d0
!!!!    character(len=200)  :: line_x,line_y,line_z
!!!!    character(len=200)  :: line_imol, keyword_imol, value_imol
!!!!!
!!!!    x = zero
!!!!    y = zero
!!!!    z = zero
!!!!!
!!!!    if (index(line,'[').gt.0) then !lettura [IMol]
!!!!!
!!!!       istart    = move_to_character(line,'[')
!!!!       atomname  = line(1:istart-1)
!!!!       iend      = move_to_character(line,']')-1
!!!!       line_imol = line(istart+1:iend)
!!!!!
!!!!       istart       = move_to_character(line_imol,'=')
!!!!       keyword_imol = trim(adjustl(line_imol(1:istart-1)))
!!!!!
!!!!       if(keyword_imol.ne.'imol') call out_%error('Keyword '//keyword_imol//' not recognized')
!!!!!
!!!!       value_imol = trim(adjustl(line_imol(istart+1:)))
!!!!       value_imol = sweep_blanks(value_imol)
!!!!!
!!!!       read(value_imol,'(i5)',iostat=ierr) IMol
!!!!       if(ierr.ne.0) call out_%error("IMol= "//trim(value_imol)//&
!!!!                                     " is not an integer")
!!!!!    
!!!!       iend = move_to_character(line,']')
!!!!       line = line(iend+1:)
!!!!!     
!!!!       iend      = move_to_character(line,' ')
!!!!       line_X    = line(1:iend)
!!!!       line_X    = sweep_blanks(line_X)
!!!!       call check_float(out_%iunit,atomname(1:2)//" X (input geometry)",line_X)
!!!!       read(line_X, '(f25.16)') x
!!!!!
!!!!       line      = line(iend+1:200)
!!!!       iend      = move_to_character(line,' ')
!!!!       line_Y    = line(1:iend)
!!!!       line_Y    = sweep_blanks(line_Y)
!!!!       call check_float(out_%iunit,atomname(1:2)//" Y (input geometry)",line_Y)
!!!!       read(line_Y, '(f25.16)') y 
!!!!!
!!!!       line      = line(iend+1:200)
!!!!       line_Z    = adjustl(line)
!!!!       call check_float(out_%iunit,atomname(1:2)//" Z (input geometry)",line_Z)
!!!!       read(line_Z, '(f25.16)') z 
!!!!!
!!!!    else !no IMol
!!!!!
!!!!       IMol      = 1
!!!!       iend      = move_to_character(line,' ')
!!!!       atomname  = line(1:iend-1)
!!!!       line      = line(iend+1:)
!!!!!     
!!!!       iend      = move_to_character(line,' ')
!!!!       line_X    = line(1:iend)
!!!!       line_X    = sweep_blanks(line_X)
!!!!       call check_float(out_%iunit,atomname(1:2)//" X (input geometry)",line_X)
!!!!       read(line_X, '(f25.16)') x
!!!!!
!!!!       line      = line(iend+1:200)
!!!!       iend      = move_to_character(line,' ')
!!!!       line_Y    = line(1:iend)
!!!!       line_Y    = sweep_blanks(line_Y)
!!!!       call check_float(out_%iunit,atomname(1:2)//" Y (input geometry)",line_Y)
!!!!       read(line_Y, '(f25.16)') y 
!!!!!
!!!!       line      = line(iend+1:200)
!!!!       line_Z    = adjustl(line)
!!!!       call check_float(out_%iunit,atomname(1:2)//" Z (input geometry)",line_Z)
!!!!       read(line_Z, '(f25.16)') z 
!!!!!
!!!!    endif
!!!!!
!!!!  end subroutine read_atoms_imols_xyz
!!!!!-----------------------------------------------------------------------
!!!!  subroutine assign_parameters(AtomName_geometry)
!!!!!
!!!!!   assign parameters
!!!!!
!!!!    implicit none
!!!!!
!!!!!   input/output variables
!!!!!
!!!!    character(len=200), dimension(target_%n_atoms), intent(in)       :: AtomName_geometry
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer :: i, j
!!!!    integer :: IErr
!!!!    logical :: found_atomname, found_atomtype
!!!!!
!!!!!   first allocation of the different vectors
!!!!!
!!!!    allocate(target_%atomic_number(target_%n_atoms),Stat=IErr)
!!!!    if(ierr.gt.0) call out_%error('Not Enough Space for atomic numbers')
!!!!!
!!!!    allocate(target_%map_atomtypes(target_%n_atoms),stat=ierr)
!!!!    if(ierr.ne.0) call out_%error('not enough space for map_atomtypes')
!!!!!
!!!!!   assign atomic variables according to atomic parameters
!!!!!
!!!!    do i = 1, target_%n_atoms
!!!!!
!!!!!    assign atomic numbers
!!!!!
!!!!       target_%atomic_number(i) = assign_atomic_number(AtomName_geometry(i))
!!!!!     
!!!!       found_atomname = .false.
!!!!!
!!!!       do j = 1, target_%n_atomtypes
!!!!!     
!!!!          if(AtomName_geometry(i).eq.target_%atom_type(j)) then
!!!!!         
!!!!            target_%map_atomtypes(i) = j
!!!!!
!!!!            found_atomname = .true.
!!!!            exit
!!!!!
!!!!          endif
!!!!!
!!!!          if(j.eq.target_%n_atomtypes.and..not.found_atomname) &
!!!!             call out_%error("No parameters given for atomtype: "//trim(AtomName_geometry(i)))
!!!!!
!!!!       enddo
!!!!!     
!!!!    enddo
!!!!!
!!!!!   Assign NMol and target_%n_var
!!!!!
!!!!    fq%n_mol = maxval(fq%i_mol)
!!!!!
!!!!    If(trim(target_%forcefield).eq.'fq'.or.trim(target_%forcefield).eq.'fq_pqeq') then
!!!!!
!!!!       target_%n_var = target_%n_atoms + fq%n_mol
!!!!!
!!!!    else
!!!!!
!!!!       target_%n_var = 4*target_%n_atoms + fq%n_mol
!!!!!
!!!!    endif
!!!!!
!!!!!   check that you did not provide parameters for different atoms
!!!!!
!!!!    do i = 1, target_%n_atomtypes
!!!!!     
!!!!       found_atomtype = .false.
!!!!!     
!!!!       do j = 1, target_%n_atoms
!!!!!     
!!!!          if(target_%atom_type(i).eq.AtomName_geometry(j)) found_atomtype = .true.
!!!!!
!!!!          if(j.eq.target_%n_atoms.and..not.found_atomtype) &
!!!!             call out_%error("No Atomtype "//trim(target_%atom_type(i))//" in input geometry.")
!!!!!
!!!!       enddo
!!!!!     
!!!!    enddo 
!!!!!
!!!!!
!!!!  end subroutine assign_parameters
!!!!!-----------------------------------------------------------------------
!!!!  real*8 function assign_atomic_number(AtomName) result(AtmNumber)
!!!!!
!!!!!   Function to assign atomic number
!!!!!
!!!!    implicit none
!!!!!
!!!!    character(len=*), intent(in)      :: AtomName
!!!!!
!!!!    AtmNumber = 0.0d0
!!!!    If(trim(AtomName).eq.'h') then
!!!!       AtmNumber = 1.0d0
!!!!    else If(trim(AtomName).eq.'hw') then
!!!!       AtmNumber = 1.0d0
!!!!    else If(trim(AtomName).eq.'h-hw') then
!!!!       AtmNumber = 1.0d0
!!!!    else If(trim(AtomName).eq.'he') then
!!!!       AtmNumber = 2.0d0
!!!!    else If(trim(AtomName).eq.'li') then
!!!!       AtmNumber = 3.0d0
!!!!    else If(trim(AtomName).eq.'be') then
!!!!       AtmNumber = 4.0d0
!!!!    else If(trim(AtomName).eq.'b') then
!!!!       AtmNumber = 5.0d0
!!!!    else If(trim(AtomName).eq.'c') then
!!!!       AtmNumber = 6.0d0
!!!!    else If(trim(AtomName).eq.'n') then
!!!!       AtmNumber = 7.0d0
!!!!    else If(trim(AtomName).eq.'o') then
!!!!       AtmNumber = 8.0d0
!!!!    else If(trim(AtomName).eq.'ow') then
!!!!       AtmNumber = 8.0d0
!!!!    else If(trim(AtomName).eq.'o-ow') then
!!!!       AtmNumber = 8.0d0
!!!!    else If(trim(AtomName).eq.'f') then
!!!!       AtmNumber = 9.0d0
!!!!    else If(trim(AtomName).eq.'ne') then
!!!!       AtmNumber = 10.0d0
!!!!    else If(trim(AtomName).eq.'na') then
!!!!       AtmNumber = 11.0d0
!!!!    else If(trim(AtomName).eq.'mg') then
!!!!       AtmNumber = 12.0d0
!!!!    else If(trim(AtomName).eq.'al') then
!!!!       AtmNumber = 13.0d0
!!!!    else If(trim(AtomName).eq.'si') then
!!!!       AtmNumber = 14.0d0
!!!!    else If(trim(AtomName).eq.'p') then
!!!!       AtmNumber = 15.0d0
!!!!    else If(trim(AtomName).eq.'s') then
!!!!       AtmNumber = 16.0d0
!!!!    else If(trim(AtomName).eq.'cl') then
!!!!       AtmNumber = 17.0d0
!!!!    else If(trim(AtomName).eq.'ar') then
!!!!       AtmNumber = 18.0d0
!!!!    else If(trim(AtomName).eq.'ag') then
!!!!       AtmNumber = 47.0d0
!!!!    else If(trim(AtomName).eq.'cu') then
!!!!       AtmNumber = 29.0d0
!!!!    else If(trim(AtomName).eq.'au') then
!!!!       AtmNumber = 79.0d0
!!!!    else 
!!!!       call out_%error('AtomName: '//trim(AtomName)//' not recognized for atomic symbol')
!!!!    endif
!!!!!
!!!!  end function assign_atomic_number
!!!!!-----------------------------------------------------------------------
!!!!  subroutine check_inconsistencies_input_file()
!!!!!
!!!!!   check input file for inconsistencies
!!!!!
!!!!    implicit none
!!!!!
!!!!!   internal variables
!!!!!
!!!!    integer  :: i,j,k
!!!!    real(dp) :: distIJ
!!!!    real(dp), parameter :: zero = 0.0d0
!!!!!
!!!!    If(out_%what.eq.'current'.and..not.general%principal_axis) &
!!!!       call out_%error('Current calculation but no principal axes')
!!!!    If(out_%what.eq.'energy'.and.(field_%static.or.field_%dynamic)) &
!!!!       call out_%error('Energy calculation but external field')
!!!!    If(out_%what.eq.'static response'.and..not.field_%static)&
!!!!       call out_%error('Static Response calculation but no static field')
!!!!    If(out_%what.eq.'dynamic response'.and..not.field_%dynamic)&
!!!!       call out_%error('Cross section calculation but no dynamic field')
!!!!    If(general%restart.and..not.field_%dynamic) &
!!!!       call out_%error('Requested restart but no dynamic field')
!!!!    If(algorithm%on_the_fly.and..not.algorithm%iterative)&
!!!!       call out_%error('Onthefly but not Iterative' )
!!!!    If(field_%dynamic.and.field_%n_freq.eq.0) &
!!!!       call out_%error('Dynamic Field but no external frequency' )
!!!!    If(field_%static.or.field_%dynamic) then 
!!!!       If(field_%rhs_form.ne.'potential'.and.field_%rhs_form.ne.'field') &
!!!!          call out_%error('External Field but no RHS selected (field or potential)')
!!!!    endif
!!!!    If(target_%forcefield.eq.'fqfmu'.or.target_%forcefield.eq.'fqfmu_pqeq') then
!!!!       If(target_%kernel .eq. 'coulomb') then
!!!!          call out_%error('FQFMu or FQFMu_PQEq force field but Kernel = coulomb')
!!!!       else if(target_%kernel .eq. 'ohno') then 
!!!!          call out_%error('FQFMu or FQFMu_PQEq force field but Kernel = ohno')
!!!!       endif
!!!!    endif
!!!!          
!!!!    If(target_%name_.eq.'wfq'.or.target_%name_.eq.'wfqib') then 
!!!!       If(.not.field_%dynamic) &
!!!!          call out_%error('wFQ or wFQIB but no dynamic field')
!!!!       If(target_%forcefield.ne.'fq'.and.target_%forcefield.ne.'fq_pqeq') &
!!!!          call out_%error('wFQ or wFQIB but forcefield is not set to FQ or FQ_PQEq')
!!!!       If(target_%kernel .eq. 'coulomb') then
!!!!          call out_%error('wFQ or wFQIB force field but Kernel = coulomb')
!!!!       else if(target_%kernel .eq. 'ohno') then 
!!!!          call out_%error('wFQ or wFQIB force field but Kernel = ohno')
!!!!       endif
!!!!    else If(target_%name_.eq.'wfqfmu') then
!!!!       If(.not.field_%dynamic) &
!!!!          call out_%error('wFQFMu but no dynamic field')
!!!!       If(target_%forcefield.ne.'fqfmu'.and.target_%forcefield.ne.'fqfmu_pqeq') &
!!!!          call out_%error('wFQFMu but forcefield is not set to FQFMu')
!!!!       If(target_%kernel .eq. 'coulomb') then
!!!!          call out_%error('wFQFMu force field but Kernel = coulomb')
!!!!       else if(target_%kernel .eq. 'ohno') then 
!!!!          call out_%error('wFQFMu force field but Kernel = ohno')
!!!!       endif
!!!!    endif
!!!!!
!!!!    If(target_%name_.eq.'wfqib') then
!!!!!
!!!!       If(len_trim(wfqfmu%epsilon_w(1)).eq.0) &
!!!!          call out_%error('wFQIB but not epsilon specification')
!!!!       If(trim(wfqfmu%epsilon_w(1)).ne.'etchegoin') then
!!!!          Write(out_%iunit,'(1x,a)') 'wFQIB: epsilon '//trim(wfqfmu%epsilon_w(1))//' not recognised'
!!!!          call out_%error('wFQIB: only Etchegoin epsilon implemented (only for Ag)')
!!!!       else
!!!!          if(trim(wfqfmu%epsilon_w(1)).eq.'etchegoin'.and.any(target_%atomic_number.ne.47.0d0)) &
!!!!              call out_%error('wFQIB: Etchegoin epsilon implemented only for Ag')
!!!!       endif
!!!!!
!!!!    endif
!!!!!    
!!!!    If(field_%dynamic) target_%n_var = target_%n_var - fq%n_mol
!!!!!
!!!!!   check no repetition of atoms
!!!!!
!!!!    do i = 1, target_%n_atoms
!!!!!
!!!!       do j = 1, i-1
!!!!!
!!!!         distIJ = dsqrt((target_%coord(1,I)-target_%coord(1,J))**2 + &
!!!!                        (target_%coord(2,I)-target_%coord(2,J))**2 + &
!!!!                        (target_%coord(3,I)-target_%coord(3,J))**2 )
!!!!!                      
!!!!         if(distIJ.lt.1.0d-7) call out_%error("Two atoms are one over the other")
!!!!!
!!!!       enddo
!!!!!
!!!!    enddo
!!!!!
!!!!!   check imol
!!!!!
!!!!    fq%n_atoms_per_molecule = 1
!!!!    k = 1
!!!!    do i = 1, target_%n_atoms - 1
!!!!!
!!!!       if(fq%i_mol(i).eq.fq%i_mol(i+1)) then
!!!!!     
!!!!          if(fq%i_mol(i).eq.fq%i_mol(i+1)) then
!!!!!
!!!!            fq%n_atoms_per_molecule(k) = fq%n_atoms_per_molecule(k) + 1
!!!!!
!!!!          endif
!!!!!
!!!!       else
!!!!!
!!!!          if(fq%i_mol(i+1).ne.(fq%i_mol(i)+1)) then
!!!!!
!!!!             Write(out_%iunit,'(1x,a,i5,a)') "Molecule ",fq%i_mol(i)," is not present"
!!!!             Stop
!!!!!
!!!!          endif
!!!!!
!!!!          k = k + 1
!!!!!
!!!!       endif
!!!!!
!!!!    enddo
!!!!!
!!!!    If(target_%name_.eq.'wfq'.or. & 
!!!!       target_%name_.eq.'wfqib'.or. & 
!!!!       target_%name_.eq.'wfqfmu') then
!!!!!
!!!!       if(any(parameters%scaling.eq.zero)) &
!!!!         call out_%error(trim(target_%name_)//': scaling = 0.0d0 ')
!!!!       if(any(parameters%tau.eq.zero).and.target_%name_.ne.'wfqib') &
!!!!         call out_%error(trim(target_%name_)//': Tau = 0.0d0 ')
!!!!       if(any(parameters%sigma_0.eq.zero).and.any(target_%atomic_number.ne.6.0d0).and.target_%name_.ne.'wfqib') &
!!!!         call out_%error(trim(target_%name_)//': Sigma0 = 0.0d0 ')
!!!!       if(any(parameters%first_length.eq.zero).and.any(target_%atomic_number.eq.6.0d0)) &
!!!!         call out_%error(trim(target_%name_)//': Primary Length = 0.0d0 ')
!!!!       if(any(parameters%A_ij.eq.zero)) &
!!!!         call out_%error(trim(target_%name_)//': RI = 0.0d0')
!!!!       if(any(parameters%fermi_d.eq.zero)) &
!!!!         call out_%error(trim(target_%name_)//': fermi function d = 0.0d0')
!!!!       if(any(parameters%fermi_s.eq.zero)) &
!!!!         call out_%error(trim(target_%name_)//': fermi function s = 0.0d0')
!!!!       if(any(parameters%factor_density.eq.zero).and.any(target_%atomic_number.eq.6.0d0)) &
!!!!         call out_%error(trim(target_%name_)//' for graphene: factor density = 0.0d0')
!!!!       if(any(parameters%second_length.eq.zero).and. &
!!!!          any(parameters%graphene_geometry.eq.'hole').and. &
!!!!          any(parameters%graphene_geometry.eq.'ribbon').and. &
!!!!          any(parameters%graphene_geometry.eq.'nanotube').and. &
!!!!          any(parameters%graphene_geometry.eq.'cnt').and. &
!!!!          any(parameters%graphene_geometry.eq.'ring').and. &
!!!!          any(parameters%graphene_geometry.eq.'bowtie-cnt')) &
!!!!          call out_%error(trim(target_%name_)//' for graphene: secondary length = 0.0d0')
!!!!       if(any(parameters%third_length.eq.zero).and. &
!!!!          any(parameters%graphene_geometry.eq.'bowtie-cnt')) &
!!!!          call out_%error(trim(target_%name_)//' for graphene: third length = 0.0d0')
!!!!!
!!!!    endif
!!!!!
!!!!!
!!!!  end subroutine check_inconsistencies_input_file
!!!!!---------------------------------------------------------------------------
!!!!  subroutine print_input_info(inp_)
!!!!!
!!!!!   Printing input information on output
!!!!!
!!!!    implicit none
!!!!
!!!!    class(inp_type) :: inp_
!!!!!
!!!!    integer :: i
!!!!!    
!!!!     ! 
!!!!     ! !$omp parallel do private(i)
!!!!     ! do i = 1, field_%n_freq
!!!!     !    call  freqnmtoev(field_%freq(i),field_%freq(i))
!!!!     ! enddo
!!!!     ! !$omp end parallel do
!!!!     write(out_%iunit,'(23x,a)') "Input File:  "//trim(inp_%filename)
!!!!     write(out_%iunit,'(23x,a)') "Output File: "//trim(out_%filename)
!!!!     if(general%save_info) write(out_%iunit,'(23x,a)') &
!!!!       "Info File:   "//out_%filename(1:len_trim(out_%filename))//'tar.gz'
!!!!     write(out_%iunit,'(/,23x,a,i5)') "OMP Threads: ", algorithm%n_threads_OMP
!!!!!
!!!!!    general informations
!!!!!
!!!!     write(out_%iunit,out_%sticks) 
!!!!     write(out_%iunit,'(23x,a)') "What          : "//trim(out_%what)
!!!!     write(out_%iunit,'(23x,a,i1)') "Verbose       : ", out_%ivrb
!!!!     write(out_%iunit,'(23x,a,l1)') "Principal Axis: ",general%principal_axis
!!!!     write(out_%iunit,'(23x,a,l1)') "Restart       : ",general%restart
!!!!     write(out_%iunit,out_%sticks) 
!!!!!     
!!!!!    field informations
!!!!!
!!!!     if (field_%static) then
!!!!!     
!!!!        write(out_%iunit,'(23x,a)') "Static External Field"
!!!!        write(out_%iunit,'(23x,a,e11.4,a)') "Field Intensity :", field_%e_0, " au"
!!!!!     
!!!!     else if (field_%dynamic) then
!!!!!     
!!!!        write(out_%iunit,'(23x,a)') "Dynamic External Field"
!!!!        write(out_%iunit,'(23x,a,e11.4,a)') "Field Intensity :", field_%e_0, " au"
!!!!!     
!!!!     else
!!!!!     
!!!!        write(out_%iunit,'(23x,a)') "Ground State"
!!!!!     
!!!!     endif
!!!!!     
!!!!     if (field_%n_freq.ge.1) write(out_%iunit,'(23x,a,i5)') "Num. Freq.      :", field_%n_freq
!!!!!     
!!!!     if (field_%n_freq.gt.1) then 
!!!!!     
!!!!        write(out_%iunit,'(23x,a,e11.4,a)') "Min. Freq.      :", field_%min_freq," eV"
!!!!        write(out_%iunit,'(23x,a,e11.4,a)') "Max. Freq.      :", field_%max_freq," eV"
!!!!        if (out_%ivrb.ge.2) then
!!!!           do i = 1, field_%n_freq
!!!!              write(out_%iunit,'(23x,a,i5,a,e10.3,a)') "Freq (",i,")    :",field_%freq(i), " eV"
!!!!           enddo
!!!!        endif
!!!!!
!!!!     else if(field_%n_freq.eq.1) then
!!!!!
!!!!        write(out_%iunit,'(23x,a,e11.4,a)') "Freq.           :", field_%freq(1), " nm"
!!!!!
!!!!     endif
!!!!!
!!!!     write(out_%iunit,out_%sticks)
!!!!!     
!!!!!    algorithm informations
!!!!!
!!!!     write(out_%iunit,'(23x,a,l1)') "Matrix in parallel : ", algorithm%matrix_in_parallel
!!!!     write(out_%iunit,'(23x,a,l1)') "Freq in parallel   : ", algorithm%freq_in_parallel
!!!!     if (algorithm%inversion) then
!!!!        write(out_%iunit,'(23x,a)') "Inversion Solution "
!!!!     else if (algorithm%iterative) then
!!!!        write(out_%iunit,'(23x,a,/)') "Iterative Solution "
!!!!        if(algorithm%on_the_fly) write(out_%iunit,'(23x,a/)') "On the Fly algorithm"
!!!!        write(out_%iunit,'(23x,a,l1)') "Do Not Change Thresholds : ", algorithm%do_not_change_thresholds
!!!!        write(out_%iunit,'(23x,a,l1)') "Changed Thresholds       : ", algorithm%changed_thresholds
!!!!        write(out_%iunit,'(23x,a,i5)') "Number of Iterations     :", algorithm%n_iter
!!!!        write(out_%iunit,'(23x,a,i5)') "Number of DIIS cycles    :", algorithm%n_diis
!!!!        write(out_%iunit,'(23x,a,i5)') "Norm to check converge   :", algorithm%norm_for_convergence
!!!!        write(out_%iunit,'(23x,a,e10.3)') "Threshold to converge    :", algorithm%threshold
!!!!     endif
!!!!     write(out_%iunit,out_%sticks) 
!!!!!
!!!!!    Force Field informations
!!!!!
!!!!     if(target_%name_.eq.'fq'        .or. &
!!!!        target_%name_.eq.'fqfmu'     .or. &
!!!!        target_%name_.eq.'wfq'       .or. &
!!!!        target_%name_.eq.'wfqfmu'    .or. &
!!!!        target_%name_.eq.'wfqib'     .or. &
!!!!        target_%name_.eq.'wfq_bem'   .or. &
!!!!        target_%name_.eq.'wfqfmu_bem'.or. &
!!!!        target_%name_.eq.'wfqib_bem') then
!!!!!        
!!!!        write(out_%iunit,'(23x,a,/)') "Target     : "// trim(target_%name_)
!!!!!        
!!!!        write(out_%iunit,'(23x,a)') "ForceField : "// trim(target_%forcefield)
!!!!        write(out_%iunit,'(23x,a)') "Kernel     : "// trim(target_%kernel)
!!!!!        
!!!!        if(target_%forcefield.eq.'fq_pqeq') then
!!!!           write(out_%iunit,'(23x,a,l1)') "PQEq       : ",fq%pqeq
!!!!        else if(target_%forcefield.eq.'fqfmu_pqeq') then
!!!!           write(out_%iunit,'(23x,a,l1)') "PQEq       : ",fqfmu%pqeq
!!!!        endif
!!!!!
!!!!        if(field_%static.or.field_%dynamic) &
!!!!           write(out_%iunit,'(23x,a)') "RHS        : "// trim(field_%rhs_form)
!!!!!
!!!!!       atomtypes informations
!!!!!
!!!!        write(out_%iunit,'(/,23x,a,i3)') "Num. AtomTypes : ", target_%n_atomtypes
!!!!!        
!!!!        call target_%print_atomtypes()
!!!!!
!!!!        write(out_%iunit,out_%sticks)
!!!!!
!!!!     endif
!!!!!
!!!!     if(bem%exists) call bem%print_parameters()
!!!!!
!!!!     if(target_%name_.ne.'bem') call target_%print_coord("Input")
!!!!!
!!!!     if(bem%exists) call bem%print_coord_bem("Input")
!!!!!
!!!!     flush(out_%iunit)
!!!!! 
!!!!  end subroutine print_input_info
!!!!!----------------------------------------------------------------------
!!!!  subroutine initialize_save_info()
!!!!!
!!!!!
!!!!    implicit none
!!!!!
!!!!!   internal variables 
!!!!!
!!!!    integer              :: iost
!!!!    integer              :: i, j !lookhere: PGI I have added j
!!!!    real(dp)             :: freqeV
!!!!    real(dp), parameter  :: ToAngs = 1.0d0/1.8897261254578281d0
!!!!!
!!!!    write(out_%info_file,'(a)') out_%filename(1:len_trim(out_%filename)-4) // '.info'
!!!!    If(out_%ivrb.ge.1) then
!!!!       write(out_%iunit,'(1x,a)') "Created file info : " // trim(out_%info_file)
!!!!       write(out_%iunit,out_%sticks) 
!!!!    endif
!!!!!    
!!!!    open(unit=out_%unit_info,file=trim(out_%info_file),status="unknown",iostat=iost,err=03)
!!!!!
!!!!!      write info in the info file
!!!!!      first we write all the infos in order
!!!!!
!!!!       write(out_%unit_info,'(a,i10)')   'NAtoms     : ', target_%n_atoms
!!!!       write(out_%unit_info,'(a,i10)')   'NMolec     : ', fq%n_mol
!!!!       write(out_%unit_info,'(a,i10)')   'NVar       : ', target_%n_var
!!!!       if(field_%static) then
!!!!          write(out_%unit_info,'(a,i10)')   'IEFld      : ', 1
!!!!       else if(field_%dynamic) then
!!!!          write(out_%unit_info,'(a,i10)')   'IEFld      : ', 2
!!!!       else
!!!!          write(out_%unit_info,'(a,i10)')   'IEFld      : ', 0
!!!!       endif
!!!!       if(target_%kernel.eq.'coulomb') then
!!!!          write(out_%unit_info,'(a,i10)')   'IKern      : ', 0
!!!!       else if(target_%kernel.eq.'ohno') then
!!!!          write(out_%unit_info,'(a,i10)')   'IKern      : ', 1
!!!!       else if(target_%kernel.eq.'gaussian') then
!!!!          write(out_%unit_info,'(a,i10)')   'IKern      : ', 2
!!!!       endif
!!!!       write(out_%unit_info,'(a,a)')  'Model      : ',trim(target_%name_)
!!!!       write(out_%unit_info,'(a,i10)')   'NumExFreq  : ', field_%n_freq
!!!!       write(out_%unit_info,'(a,E13.6)') 'Input |E|  : ', field_%e_0
!!!!       if(target_%name_.ne.'bem') then
!!!!          write(out_%unit_info,'(a,l10)')   'PrincAxes  : ', general%principal_axis
!!!!       elseif(target_%name_.eq.'bem') then
!!!!          write(out_%unit_info,'(a,l10)')   'BEM-IEF    : ', bem%rhs_ief   
!!!!       endif
!!!!       write(out_%unit_info,'(a)')       'ForceField : '//target_%forcefield
!!!!       write(out_%unit_info,'(a)') 'Input Geometry (angstrom)'
!!!!!   
!!!!       if(target_%name_.ne.'bem') then
!!!!!lookhere: pgi campo
!!!!          do i = 1, target_%n_atoms
!!!!             write(out_%unit_info,'(f4.1,2x,i8,3(1x,f25.16))')  &
!!!!                              target_%atomic_number(i), &
!!!!                              fq%i_mol(i),              &
!!!!                              target_%coord(1,i)*ToAngs,       &
!!!!                              target_%coord(2,i)*ToAngs,       &
!!!!                              target_%coord(3,i)*ToAngs
!!!!          enddo 
!!!!!
!!!!          if(target_%name_.eq.'wfq_bem'.or.target_%name_.eq.'wfqfmu_bem') then
!!!!             write(out_%unit_info,'(a)') 'Input BEM Geometry (angstrom)'
!!!!             do i = 1, bem%n_var 
!!!!                write(out_%unit_info,'(4(1x,f25.16))')  &
!!!!                       bem%area(i)*(ToAngs**two),   &
!!!!                       bem%coord(1,i)*ToAngs,       &
!!!!                       bem%coord(2,i)*ToAngs,       &
!!!!                       bem%coord(3,i)*ToAngs
!!!!             enddo
!!!!          endif
!!!!!end here
!!!!!   
!!!!          write(out_%unit_info,'(a)') "Parameters: "//target_%forcefield
!!!!          !write(out_%unit_info,'(a,i3)') "Num. AtomTypes : ", target_%n_atomtypes
!!!!          if(target_%forcefield.eq.'fq') then
!!!!!lookhere: PGI I am doing this fast
!!!!             do j = 1, target_%n_atoms
!!!!                do i = 1, target_%n_atomtypes
!!!!                   write(out_%unit_info,'(a2,2(1x,f25.16))')     &
!!!!                                     target_%atom_type(i), &
!!!!                                     fq%chi,               &
!!!!                                     fq%eta(i)
!!!!                enddo 
!!!!             enddo
!!!!!end here
!!!!          else if(target_%forcefield.eq.'fq_pqeq') then
!!!!!lookhere: PGI I am doing this fast
!!!!             do j = 1, target_%n_atoms
!!!!                do i = 1, target_%n_atomtypes
!!!!                   write(out_%unit_info,'(a2,3(1x,f25.16))')     &
!!!!                                     target_%atom_type(i), &
!!!!                                     fq%chi(i),            &
!!!!                                     fq%eta(i),            &
!!!!                                     fq%r_q(i)
!!!!                enddo 
!!!!             enddo
!!!!!end here
!!!!          else if(target_%forcefield.eq.'fqfmu') then
!!!!!lookhere: PGI I am doing this fast
!!!!             do j = 1, target_%n_atoms
!!!!                do i = 1, target_%n_atomtypes
!!!!                   write(out_%unit_info,'(a2,3(1x,f25.16))')     &
!!!!                                     target_%atom_type(i), &
!!!!                                     fq%chi(i),            &
!!!!                                     fq%eta(i),            &
!!!!                                     fqfmu%alpha(i)
!!!!                enddo 
!!!!             enddo
!!!!!end here
!!!!          else if(target_%forcefield.eq.'fqfmu_pqeq') then
!!!!!lookhere: PGI I am doing this fast
!!!!             do j = 1, target_%n_atoms
!!!!                do i = 1, target_%n_atomtypes
!!!!                   write(out_%unit_info,'(a2,5(1x,f25.16))')     &
!!!!                                     target_%atom_type(i), &
!!!!                                     fq%chi(i),            &
!!!!                                     fq%eta(i),            &
!!!!                                     fqfmu%alpha(i),       &
!!!!                                     fq%r_q(i),            &
!!!!                                     fqfmu%r_mu(i)
!!!!                enddo 
!!!!             enddo
!!!!!end here
!!!!          endif
!!!!!
!!!!!         NAtoms_per_IMol
!!!!!
!!!!          write(out_%unit_info,'(a)') 'Number of atoms per IMol'
!!!!          do i = 1, fq%n_mol
!!!!             write(out_%unit_info,'(i8,2x,i8)') i, fq%n_atoms_per_molecule(i)
!!!!          enddo
!!!!!
!!!!       elseif(target_%name_.eq.'bem') then
!!!!!
!!!!          do i = 1, bem%n_var
!!!!             write(out_%unit_info,'(4(1x,f25.16))')  &
!!!!                    bem%area(i)*(ToAngs**two),   &
!!!!                    bem%coord(1,i)*ToAngs,       &
!!!!                    bem%coord(2,i)*ToAngs,       &
!!!!                    bem%coord(3,i)*ToAngs
!!!!!
!!!!          enddo
!!!!!
!!!!       endif
!!!!!
!!!!       if(target_%name_.eq.'wfq'.or. &
!!!!          target_%name_.eq.'wfqib'.or. &
!!!!          target_%name_.eq.'wfqfmu'.or. &
!!!!          target_%name_.eq.'wfq_bem'.or. &
!!!!          target_%name_.eq.'wfqib_bem'.or. &
!!!!          target_%name_.eq.'wfqfmu_bem') then
!!!!!
!!!!          write(out_%unit_info,'(a)') "Dynamic parameters: 9"
!!!!          do i = 1, target_%n_atomtypes
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%tau(i)
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%sigma_0(i)      
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%A_ij(i)         
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%fermi_d(i)      
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%fermi_s(i)      
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%first_length(i) 
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%second_length(i)
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%third_length(i) 
!!!!             write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%factor_density(i)
!!!!          enddo
!!!!!
!!!!       endif
!!!!!                                   
!!!!!      Frequenze
!!!!!
!!!!       if(field_%n_freq.gt.0) then
!!!!          write(out_%unit_info,'(a)') 'Frequencies (eV)'
!!!!          do i = 1, field_%n_freq
!!!!             call FreqnmtoeV(field_%freq(i),FreqeV)
!!!!             write(out_%unit_info,'(f25.16)') FreqeV
!!!!          enddo
!!!!       endif
!!!!!
!!!!!
!!!!    03   Continue      
!!!!    Close(out_%unit_info)
!!!!!
!!!!!
!!!!  end subroutine initialize_save_info
!----------------------------------------------------------------------
end module input_module

!----------------------------------------------------------------------
module output_module
!      
!   Module output
!
!!!    use parameters_module
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
!!      integer  :: unit_info = 8 ! unit of the info file
      integer  :: iunit = 12 ! unit of the file
!!      integer  :: ivrb  = 0  ! Verbose modality 
!!!
!!      logical  :: exists = .false.
!!!
      character (len=12)  :: sticks = "(1x,80(1h-))"
      character (len=200) :: filename
!!      character (len=200) :: info_file
!!      character (len=200) :: what
!!!
      contains
!!!
      procedure :: out_file_fill
      procedure :: print_banner
!!      procedure :: print_matrix
      procedure :: warning
      procedure :: error
!!      procedure :: clean_up_scratch
!!      procedure :: make_targz
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
!    print banner nanoFQ  
!
     class(out_type)  :: out_
!
     Write(out_%iunit,out_%sticks) 
     Write(out_%iunit,'(22x,a)') "                                      "
     Write(out_%iunit,'(22x,a)') "           =============              "
     Write(out_%iunit,'(22x,a)') "           FRET EMBEDLAB              "
     Write(out_%iunit,'(22x,a)') "           =============              "
     Write(out_%iunit,'(22x,a)') "                                      "
     Write(out_%iunit,'(22x,a)') "                                      "
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
!----------------------------------------------------------------------
!!!   subroutine clean_up_scratch(out_)
!!!!   
!!!     implicit none
!!!!
!!!!    cleaning up scratch
!!!!
!!!     class(out_type)  :: out_
!!!!
!!!     write(out_%iunit,out_%sticks)
!!!     write(out_%iunit,'(/ a /)') ' I am cleaning up the backup (*.nanofq.bk)'
!!!     write(out_%iunit,out_%sticks)
!!!     call system('rm *.nanofq.bk')
!!!!    
!!!   end subroutine clean_up_scratch
!!!!-----------------------------------------------------------------------
!!!   subroutine make_targz(out_)
!!!!   
!!!     implicit none 
!!!!
!!!!    create tar.gz archive 
!!!!
!!!     class(out_type)  :: out_
!!!!     
!!!!    internal variables
!!!!     
!!!     integer :: iost
!!!     integer :: unit_freq
!!!     integer :: unit_info = 8
!!!!     
!!!     character(len=100)  :: line_file_freq
!!!     character(len=100)  :: file_csv
!!!     character(len=100)  :: file_info
!!!     character(len=100)  :: file_tar
!!!     character(len=318)  :: command
!!!!     
!!!     logical             :: exist_csv
!!!     logical             :: exist_freq =.false.
!!!     logical             :: info_freq_exist =.false.
!!!!
!!!     write(file_info,'(a)') out_%filename(1:len_trim(out_%filename)-4) // '.info'
!!!     write(file_csv,'(a)') out_%filename(1:len_trim(out_%filename)-4) // '.csv'
!!!     Inquire(file=file_csv,exist=exist_csv)     
!!!!     
!!!!    check if freq file exists
!!!!     
!!!     call execute_command_line('for i in *.freq ; do test -f "$i" '// &
!!!                               '&& echo "exists one or more files" > info_freq.txt && break; done')
!!!     inquire(file="info_freq.txt",exist=info_freq_exist)
!!!     if(info_freq_exist) then
!!!        unit_freq = 13
!!!        Open(unit=unit_freq,file="info_freq.txt",status="OLD",IOSTAT=IOST,ERR=02)
!!!!
!!!             read(unit_freq,'(a)') line_file_freq
!!!             if(index(line_file_freq,"exists one or more files").gt.0) exist_freq = .true.
!!!!    
!!!        02   Continue      
!!!        Close(unit_freq)
!!!        call execute_command_line("rm info_freq.txt")
!!!     endif
!!!!
!!!!      
!!!     Open(unit=unit_info,file=file_info,status="OLD",IOSTAT=IOST, &
!!!          ACCESS='APPEND',ERR=03)
!!!!
!!!          write(unit_info,out_%sticks)
!!!          write(unit_info,'(24x,a)') 'Normal termination of nanoFQ'
!!!          write(unit_info,out_%sticks)
!!!!    
!!!     03   Continue      
!!!     Close(unit_info)
!!!!
!!!     write(file_tar,'(a)') out_%filename(1:len_trim(out_%filename)-4) // '.tar.gz'
!!!!
!!!     if(exist_csv) then 
!!!        if(exist_freq) then
!!!           command = "tar -czf "//file_tar//" "//file_info//" "//file_csv//" *.freq"
!!!        else
!!!           command = "tar -czf "//file_tar//" "//file_info//" "//file_csv
!!!        endif
!!!     else
!!!        if(exist_freq) then
!!!           command = "tar -czf "//file_tar//" "//file_info//" *.freq"
!!!        else
!!!           command = "tar -czf "//file_tar//" "//file_info
!!!        endif
!!!     endif
!!!!
!!!     call execute_command_line(trim(command))
!!!     call execute_command_line("rm "//file_info)
!!!     if(exist_freq) call execute_command_line("rm *freq")
!!!!
!!!     If(out_%ivrb.ge.1) then
!!!        write(out_%iunit,'(1x,a)') "Created file tar  : "//file_tar
!!!        write(out_%iunit,out_%sticks) 
!!!     endif
!!!!
!!!!
!!!   end subroutine make_targz
!-----------------------------------------------------------------------
end module output_module

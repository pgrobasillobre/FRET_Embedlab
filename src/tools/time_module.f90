!----------------------------------------------------------------------
module time_module
!      
!   Module time
!
    !$ use omp_lib
    use parameters_module
    !!use string_manipulation_module
    use output_module
!
    Implicit None
!
!   Public Variables
!
    public time
!
!   time type
!
    type time_type
!
      character(len=100), dimension(1) :: strings_timer
      real(dp), dimension(1)           :: sec_i !elapsed seconds initial
      real(dp), dimension(1)           :: sec_f !elapsed seconds final
      real(dp), dimension(1)           :: cpu_i !cpu seconds initial
      real(dp), dimension(1)           :: cpu_f !cpu seconds final
!
      contains
!
      procedure :: initialize  => initialize_time
      procedure :: start       => start_time
!!!!      procedure :: finish      => finish_time
!!!!      procedure :: conclude    => conclude_time
!
    end type time_type
!    
    type (time_type), Save :: time
!
   contains
!----------------------------------------------------------------------
! debugpgi: do we need this?
   subroutine initialize_time(time)
!   
!    initialize timer (create strings_ and data_timer)
!
!    1      : total runtime
!
     class(time_type) :: time
!     
     time%strings_timer(1) = "total"
!  
   end subroutine initialize_time
!----------------------------------------------------------------------
   subroutine start_time(time,string)
!    
!    start timer for specific task
!
     class(time_type)  :: time
!
     character(len=*) :: string
!
     integer :: element
!     
     if(any(string.eq.time%strings_timer)) then
!
       element = findloc(time%strings_timer,string,1)
       call get_time(time%sec_i(element),time%cpu_i(element))
!
     else
!
       call out_%error("Timer "//trim(string)//" is not initialized")
!
     endif
!  
   end subroutine start_time
!----------------------------------------------------------------------
!!!   subroutine finish_time(time,string)
!!!!   
!!!!    finish timer for specific task
!!!!
!!!     class(time_type)  :: time
!!!!
!!!     character(len=*) :: string
!!!!
!!!     integer :: element
!!!!     
!!!     if(any(string.eq.time%strings_timer)) then
!!!!
!!!       element = findloc(time%strings_timer,string,1)
!!!       call get_time(time%sec_f(element),time%cpu_f(element))
!!!!
!!!     else
!!!!
!!!       call out_%error("Timer "//trim(string)//" is not initialized")
!!!!
!!!     endif
!!!!  
!!!   end subroutine finish_time
!----------------------------------------------------------------------
   subroutine get_time(sec,cpu)
!   
!    get time
!
     real(dp) :: sec
     real(dp) :: cpu
!
     character(len=10) :: TimeINFO
     character(len=8)  :: DateINFO ! ccyymmdd
     character(len=6)  :: Second
     character(len=4)  :: Year
     character(len=2)  :: Hour, Minute, Month, Day
!
     integer           :: IYear
     integer           :: IMonth
     integer           :: IDay
     integer           :: IHour
     integer           :: IMin
     real(dp)          :: Seconds
!
     CALL DATE_AND_TIME(DateINFO,TimeINFO)
!
     Year  = DateINFO(1:4)
     Read( Year,'(i4)')  IYear 
!
     Month = DateINFO(5:6)
     Read( Month,'(i2)')  IMonth 
!
     Day   = DateINFO(7:8)
     Read( Day,'(i2)')  IDay
!
     Hour   = TimeINFO(1:2)
     Read( Hour,'(i2)')  IHour 
!
     Minute = TimeINFO(3:4)
     Read( Minute,'(i2)')  IMin 
!
     Second = TimeINFO(5:10)
     Read(Second,'(f6.3)')  Seconds
!
     If(Mod(IYear,4).eq.0) then
        sec = 366*86400
     else
        sec = 365*86400
     endif
!
     If(IMonth.eq.4.or.IMonth.eq.6.or.IMonth.eq.9.or.IMonth.eq.11) then
!
        sec = sec + IMonth*30*86400
!
     else If(IMonth.eq.2) then
!
        If(Mod(IYear,4).eq.0) then
           sec = sec + IMonth*29*86400
        else
           sec = sec + IMonth*28*86400
        endif
!
     else
!
        sec = sec + IMonth*31*86400
!
     endif
!
     sec = sec + IDay*86400 + IHour*3600 + IMin*60 + Seconds
!
     call cpu_time(cpu)
!  
   end subroutine get_time
!----------------------------------------------------------------------
!!!   subroutine conclude_time(time)
!!!!
!!!     implicit none
!!!!
!!!     class(time_type) :: time
!!!!
!!!     real(dp)         :: SecTot
!!!     real(dp)         :: ElapHours   
!!!     real(dp)         :: ElapMinutes 
!!!     real(dp)         :: ElapSeconds 
!!!!     
!!!     real(dp)         :: CPUTot      
!!!     real(dp)         :: CPUHours    
!!!     real(dp)         :: CPUMinutes  
!!!     real(dp)         :: CPUSeconds  
!!!!
!!!     character(len=10) :: TimeINFO, final_date, final_time
!!!     character(len=8)  :: DateINFO ! ccyymmdd
!!!     character(len=6)  :: Second
!!!     character(len=4)  :: Year
!!!     character(len=2)  :: Hour, Minute, Month, Day
!!!!     
!!!     1343 Format(31X,'Il programma fa sempre quello che gli dici di fare',&
!!!         /,60X,'-- C.CAPPELLI, SEMPRE')
!!!     1004 Format(43X,'     CPU Time: ',i5,' h ',i2,' min ',i2,' sec')
!!!     1005 Format(43X,' Elapsed Time: ',i5,' h ',i2,' min ',i2,' sec')
!!!     1010 Format(14x,'Normal Termination of nanoFQ program in date ',10a,2x,&
!!!          2a,2x,8a)
!!!!    
!!!     CALL DATE_AND_TIME(DateINFO,TimeINFO)
!!!!     
!!!     Year   = DateINFO(1:4)
!!!     Month  = DateINFO(5:6)
!!!     Day    = DateINFO(7:8)
!!!     Hour   = TimeINFO(1:2)
!!!     Minute = TimeINFO(3:4)
!!!     Second = TimeINFO(5:10)
!!!!     
!!!     Write(final_date,'(a)') Day // '/' // Month // '/' // Year
!!!     Write(final_time,'(a)') Hour // ':' // Minute // ':' // Second(1:2)
!!!!    
!!!     SecTot = time%sec_f(1) - time%sec_i(1)
!!!     ElapHours    = SecTot / 3600.0d0
!!!     ElapMinutes  = mod(SecTot/60.0d0,60.0d0)
!!!     ElapSeconds  = mod(SecTot,60.0d0)
!!!!     
!!!     CPUTot       = time%cpu_f(1) - time%cpu_i(1)
!!!     CPUHours     = CPUTot / 3600.0d0
!!!     CPUMinutes   = mod(CPUTot/60.0d0,60.0d0)
!!!     CPUSeconds   = mod(CPUTot,60.0d0)
!!!!     
!!!     Write(out_%iunit,'(a)')
!!!     Write(out_%iunit,1343)
!!!     Write(out_%iunit,'(a)')
!!!     Write(out_%iunit,out_%sticks)
!!!     Write(out_%iunit,'(a)')
!!!     Write(out_%iunit,1004) Int(CPUHours), Int(CPUMinutes), NInt(CPUSeconds)
!!!     Write(out_%iunit,1005) Int(ElapHours), Int(ElapMinutes), NInt(ElapSeconds)
!!!     Write(out_%iunit,'(a)')
!!!     Write(out_%iunit,out_%sticks)
!!!     Write(out_%iunit,'(a)')
!!!     Write(out_%iunit,1010) final_date,' at ',final_time
!!!     Write(out_%iunit,'(a)')
!!!     Write(out_%iunit,out_%sticks)
!!!!     
!!!   end subroutine conclude_time
!----------------------------------------------------------------------
end module time_module

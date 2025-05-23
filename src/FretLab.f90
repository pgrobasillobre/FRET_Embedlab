Program FretLab
!
use input_module
use output_module
use time_module
use target_module
use algorithm_module
!
!$ use omp_lib
!
!---------------------------------------------------------------------------
!                     ______          __  __          __          
!                    / ____/_______  / /_/ /   ____ _/ /_         
!                   / /_  / ___/ _ \/ __/ /   / __ `/ __ \        
!                  / __/ / /  /  __/ /_/ /___/ /_/ / /_/ /        
!                 /_/   /_/   \___/\__/_____/\__,_/_.___/         
!
!---------------------------------------------------------------------------
!
!                       Program by Pablo Grobas Illobre
!                             
!                         For any problem write to:
!                         pgrobasillobre@gmail.com
!   
!---------------------------------------------------------------------------
!
implicit none
!
call inp_%get_arguments()
call time%initialize()
call time%start("total")

!$ call omp_set_num_threads(parallel%n_threads_OMP) 

open(unit=out_%iunit,file=out_%filename,status="unknown")
!
    call out_%print_banner()
!
    call inp_%check_input_file()
!
    call inp_%read_() 
!
    call inp_%print_input_info()
!
!
    if(target_%name_.eq.'integrate_density') then 
!
       call algorithm%integrate_density()
!
    elseif(target_%name_.eq.'aceptor_donor') then
!
       call algorithm%eet_aceptor_donor()
!
    elseif(target_%name_.eq.'aceptor_np') then
!
       call algorithm%aceptor_np_interaction()
!
    elseif(target_%name_.eq.'aceptor_np_donor') then
!
       call algorithm%aceptor_np_donor()
!
    endif
!
!
    call time%finish("total")
    call time%conclude()
!
close(out_%iunit)
!
End Program FretLab
!---------------------------------------------------------------------------

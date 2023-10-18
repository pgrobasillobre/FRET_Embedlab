Program FRET_Embedlab
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
!                ______ _____  ______ _______                    
!               |  ____|  __ \|  ____|__   __|                   
!               | |__  | |__) | |__     | |                      
!               |  __| |  _  /|  __|    | |                      
!               | |    | | \ \| |____   | |                      
!               |_|    |_|  \_\______|  |_|                      
!
!                ______           _              _ _       _
!               |  ____|         | |            | | |     | |    
!               | |__   _ __ ___ | |__   ___  __| | | __ _| |__  
!               |  __| | '_ ` _ \| '_ \ / _ \/ _` | |/ _` | '_ \ 
!               | |____| | | | | | |_) |  __/ (_| | | (_| | |_) |
!               |______|_| |_| |_|_.__/ \___|\__,_|_|\__,_|_.__/ 
!                                                    
!---------------------------------------------------------------------------
!
!                       Program by Pablo Grobas Illobre
!                    with contributions by Sveva Sodomaco
!                             
!                         For any problem write to:
!                            pgrobasillobre@sns.it
!   
!---------------------------------------------------------------------------
!
implicit none
!
call inp_%get_arguments()
call time%initialize()
call time%start("total")

!$ call omp_set_num_threads(algorithm%n_threads_OMP) 

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
       !call algorithm%FRET_aceptor_donor()
       call out_%error("aceptor_np not supported")
!
    elseif(target_%name_.eq.'aceptor_np_donor') then
!
       !call algorithm%FRET_aceptor_donor()
       call out_%error("aceptor_np_donor not supported")
!
    endif
!
!
    call time%finish("total")
    call time%conclude()
!
close(out_%iunit)
!
End Program FRET_Embedlab
!---------------------------------------------------------------------------

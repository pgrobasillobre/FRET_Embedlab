Program FRET_Embedlab
!
use input_module
use output_module
use time_module
use target_module
!!use general_module
!$ use omp_lib
!---------------------------------------------------------------------------
!                                                        
!                             FRET EMBEDLAB                                                      
!                                                        
!                                                        
!                                                        
!                                                        
!
!---------------------------------------------------------------------------
!
!                       Program by Pablo Grobas Illobre
!                    with contributions by Sveva Sodomaco
!                             
!                       For any problem write to:
!                          pgrobasillobre@sns.it
! 
!---------------------------------------------------------------------------
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
!
    call inp_%print_input_info()
!
! CASES FRET
!!!    if(field_%static) then 
!!!!
!!!       call algorithm%solve_static_field()
!!!!
!!!    else if(field_%dynamic) then
!!!!
!!!       call algorithm%solve_dynamic_field()
!!!!
!!!    else 
!!!!
!!!       call algorithm%solve_ground_state()
!!!!
!!!    endif
!
!
    call time%finish("total")
    call time%conclude()
!
close(out_%iunit)
!
End Program FRET_Embedlab
!---------------------------------------------------------------------------

!----------------------------------------------------------------------
module target_module
!      
!   Module target
!
    !$ use omp_lib
    !use output_module
    !use parameters_module
    !use string_manipulation_module
!
    Implicit None
!
!   Public Variables
!
    public target_
!
!   target type
!
    type :: target_type
!
      character(len=200)                    :: name_      
!      
      contains
!
!!      procedure :: assign_model_parameters
!
    end type target_type
!    
!    class (target_type), pointer, save :: target_
    type (target_type), save :: target_
!
!----------------------------------------------------------------------
!contains
!----------------------------------------------------------------------
!!!   subroutine assign_model_parameters(target_)
!!!!
!!!     implicit none
!!!!     
!!!!    subroutine for assigning model parameters
!!!!
!!!     class(target_type)    :: target_
!!!!     
!!!     if(target_%name_.ne.'wfq'       .and. &
!!!        target_%name_.ne.'wfqfmu'    .and. &
!!!        target_%name_.ne.'wfqib'     .and. &
!!!        target_%name_.ne.'bem'       .and. &
!!!        target_%name_.ne.'wfq_bem'   .and. &
!!!        target_%name_.ne.'wfqfmu_bem'.and. &
!!!        target_%name_.ne.'wfqib_bem') &
!!!        call out_%error("Target name: "//trim(target_%name_)//" not recognised")
!!!!
!!!   end subroutine assign_model_parameters
!-----------------------------------------------------------------------
end module target_module

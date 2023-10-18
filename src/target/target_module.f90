!----------------------------------------------------------------------
module target_module
!      
!   Module target
!
    !$ use omp_lib
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
      character(len=200)                    :: name_, density_file, aceptor_density, donor_density, nanoparticle      
!      
      contains
!
    end type target_type
!    
    type (target_type), save :: target_
!
!-----------------------------------------------------------------------
end module target_module

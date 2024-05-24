!----------------------------------------------------------------------
module target_module
!      
!   Module target
!
    use parameters_module
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
      character(len=200) :: name_, density_file, aceptor_density, donor_density, nanoparticle      
!
      real(dp)               :: cutoff 
      real(dp)               :: omega_0                  
      real(dp)               :: spectral_overlap
      real(dp)               :: aceptor_density_rotation_angle
!
      real(dp), dimension(3) :: aceptor_transdip
!
      logical                :: calc_overlap_int
      logical                :: aceptor_density_rotate
      logical                :: aceptor_transdip_rotate
!      
      contains
!
    end type target_type
!    
    type (target_type), save :: target_
!
!-----------------------------------------------------------------------
end module target_module

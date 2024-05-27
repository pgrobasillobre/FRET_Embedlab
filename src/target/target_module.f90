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
      character(len=200) :: name_
      character(len=200) :: density_file 
      character(len=200) :: aceptor_density 
      character(len=200) :: donor_density 
      character(len=200) :: nanoparticle      
      character(len=200) :: rotation_axys
!
      real(dp)               :: cutoff 
      real(dp)               :: omega_0                  
      real(dp)               :: spectral_overlap
      real(dp)               :: aceptor_density_rotation_angle
!
      real(dp), dimension(3) :: aceptor_transdip, aceptor_transdip_rot, ref_vector
!
      logical                :: calc_overlap_int
      logical                :: aceptor_density_rotate
      logical                :: aceptor_transdip_rotate
      logical                :: aceptor_transdip_rotate_align_with
      logical                :: rotate_aceptor
!
      logical                :: debug
!      
      contains
!
    end type target_type
!    
    type (target_type), save :: target_
!
!-----------------------------------------------------------------------
end module target_module

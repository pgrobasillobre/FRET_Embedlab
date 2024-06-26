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
      real(dp)               :: donor_density_rotation_angle
      real(dp)               :: transdip_refvec_theta
      real(dp)               :: aceptor_angle_check, donor_angle_check
!
      real(dp), dimension(3) :: aceptor_transdip, aceptor_transdip_rot, aceptor_ref_vector
      real(dp), dimension(3) :: donor_transdip, donor_transdip_rot, donor_ref_vector
!
      integer                :: debug
!
      logical                :: calc_overlap_int
      logical                :: aceptor_density_rotate, donor_density_rotate
      logical                :: aceptor_transdip_rotate, donor_transdip_rotate
      logical                :: aceptor_transdip_rotate_align_with, donor_transdip_rotate_align_with
      logical                :: rotate_aceptor, rotate_donor, rotate_np
!
      contains
!
    end type target_type
!    
    type (target_type), save :: target_
!
!-----------------------------------------------------------------------
end module target_module

!----------------------------------------------------------------------
module integrals_module
!      
!   Module integrals
!
    !$ use omp_lib
    use density_module
    use nanoparticle_module
    use target_module
    use parameters_module
!
    Implicit None
!
    public integrals
!
!   integrals type
!
    type integrals_type
! 
      real(dp)                  :: aceptor_donor_coulomb
      real(dp)                  :: aceptor_donor_overlap
      real(dp), dimension(2)    :: aceptor_np_int !1: real, 2: imaginary
!
    end type integrals_type
!
    type (integrals_type), Save :: integrals 
!
   contains
!----------------------------------------------------------------------
   subroutine eet_aceptor_donor_integral(integrals,aceptor,donor)
!
!    Compute donor-aceptor potential for EET rate calculation
!
     implicit none
!
     type (density_type), intent(in)  :: aceptor, donor 
!
     type (integrals_type), intent(out) :: integrals
!
!
!    internal variables
!
     real(dp) :: x_a,y_a,z_a !position of aceptor
     real(dp) :: x_d,y_d,z_d !position of donor
     real(dp) :: r(3)        !aceptor-donor position     
     real(dp) :: dist        
     real(dp) :: invdst 
     real(dp) :: int_coulomb, int_overlap 
     real(dp) :: sf,sf0,screen_pot !for screening
!
     real(dp), dimension(:), allocatable   :: rho_aceptor, rho_donor
     real(dp), dimension(:,:), allocatable :: xyz_aceptor, xyz_donor
!
     integer  :: i,j,k,l,m,n
!
!    To allocate these quantities internally requires more memory but makes the calculation
!    in parallel faster. FORTRAN has problems accessing object-types in parallell. 
     allocate(rho_aceptor(aceptor%n_points_reduced),   rho_donor(donor%n_points_reduced),&
              xyz_aceptor(3,aceptor%n_points_reduced), xyz_donor(3,donor%n_points_reduced))
!
     xyz_aceptor(1:3,1:aceptor%n_points_reduced) = aceptor%xyz(1:3,1:aceptor%n_points_reduced)
     xyz_donor(1:3,1:donor%n_points_reduced) = donor%xyz(1:3,1:donor%n_points_reduced)
!
     rho_aceptor(1:aceptor%n_points_reduced) = aceptor%rho_reduced(1:aceptor%n_points_reduced)
     rho_donor(1:donor%n_points_reduced) = donor%rho_reduced(1:donor%n_points_reduced)
!
     integrals%aceptor_donor_coulomb = zero
     integrals%aceptor_donor_overlap = zero
!
     int_coulomb = zero
     int_overlap = zero
!
     r   = zero
!
     !$omp parallel do private(i,j,r,dist,invdst,sf,sf0,screen_pot) &
     !$omp collapse(1) & 
     !$omp reduction(+:int_coulomb,int_overlap) 
!    aceptor
     do i = 1, aceptor%n_points_reduced
!       donor
        do j = 1, donor%n_points_reduced
!
!          Aceptor-donor overlap 
           if (i.eq.j .and. target_%calc_overlap_int) then
              int_overlap = int_overlap +&
                            rho_aceptor(i) * rho_donor(j)
           endif 
!
           r(1) = (xyz_aceptor(1,i)-xyz_donor(1,j))
           r(2) = (xyz_aceptor(2,i)-xyz_donor(2,j))
           r(3) = (xyz_aceptor(3,i)-xyz_donor(3,j))
!
           dist = dsqrt(DOT_PRODUCT(r,r))
!
!          Skip when grid points are coincident to avoid instabilities
           if (dist.le.1.0e-14) go to 01
!
           invdst = one/dist
!
!          Screening function 
           sf  = dist / QMscrnFact

           sf0        = erf(sf)
           screen_pot = sf0
!
!          Integrate charges as rho_aceptor * rho_donor * (1/dist) * screening
!            --> the density has been already weigthed by the cube volume
           int_coulomb = int_coulomb +&
                         rho_aceptor(i) * rho_donor(j) * invdst * screen_pot
!   
   
           01 continue                  
!
        enddo
     enddo
     !$omp end parallel do
!
!    Sum all the contributions to the integrals
     integrals%aceptor_donor_coulomb = int_coulomb
!
!    Weigth overlap by -omega_0
     if (target_%calc_overlap_int) integrals%aceptor_donor_overlap = &
                                           -target_%omega_0 * int_overlap
!
!    Deallocate and exit
!
     deallocate(rho_aceptor)
     deallocate(rho_donor)
     deallocate(xyz_aceptor)
     deallocate(xyz_donor)
!
  end subroutine eet_aceptor_donor_integral 
!----------------------------------------------------------------------
   subroutine aceptor_nanoparticle_interaction_integral(integrals,aceptor,np)
!
!    Compute donor-aceptor potential for EET rate calculation
!
     implicit none
!
     type (density_type), intent(in)      :: aceptor
     type (nanoparticle_type), intent(in) :: np
!
     type (integrals_type), intent(out)   :: integrals
!
!
!    internal variables
!
     real(dp) :: x_a,y_a,z_a !position of aceptor
     real(dp) :: r(3)        !aceptor-np position     
     real(dp) :: dist        
     real(dp) :: invdst 
     real(dp) :: sf,sf0,sf1,screen_pot, screen_fld !for screening
     real(dp) :: aceptor_np_int_re, aceptor_np_int_im 
     real(dp) :: aceptor_np_int_re_mu, aceptor_np_int_im_mu 
!
     real(dp), dimension(:), allocatable     :: rho_aceptor
     real(dp), dimension(:,:), allocatable   :: xyz_aceptor, xyz_np
     real(dp), dimension(:,:), allocatable   :: mm_q
     real(dp), dimension(:,:,:), allocatable :: mm_mu
!
     integer  :: i,j
!
!    To allocate these quantities internally requires more memory but makes the calculation
!    in parallel faster. FORTRAN has problems accessing object-types in parallell. 
     if (np%charges) then
        allocate(rho_aceptor(aceptor%nx*aceptor%ny*aceptor%nz),   mm_q(np%natoms,2),&
                 xyz_aceptor(3,aceptor%n_points_reduced), xyz_np(3,np%natoms))
     else if (np%dipoles) then
        allocate(rho_aceptor(aceptor%nx*aceptor%ny*aceptor%nz), mm_q(np%natoms,2),  mm_mu(3,np%natoms,2),&
                 xyz_aceptor(3,aceptor%n_points_reduced), xyz_np(3,np%natoms))
     endif
!
     xyz_aceptor(1:3,1:aceptor%n_points_reduced) = aceptor%xyz(1:3,1:aceptor%n_points_reduced)
     xyz_np(1:3,1:np%natoms) = np%xyz(1:3,1:np%natoms)
!
     rho_aceptor(1:aceptor%n_points_reduced) = aceptor%rho_reduced(1:aceptor%n_points_reduced)
!
     if (np%charges) then
        mm_q(1:np%natoms,1:2) = np%q(1:np%natoms,1:2)
     else if(np%dipoles) then
        mm_q(1:np%natoms,1:2) = np%q(1:np%natoms,1:2)
        mm_mu(1:3,1:np%natoms,1:2) = np%mu(1:3,1:np%natoms,1:2)
     endif
!
     aceptor_np_int_re = zero
     aceptor_np_int_im = zero
     aceptor_np_int_re_mu = zero
     aceptor_np_int_im_mu = zero
!
     integrals%aceptor_np_int = zero
!
     r   = zero
     sf  = zero
     sf0 = zero
     sf1 = zero 
     screen_pot = zero
     screen_fld = zero
!
     !$omp parallel do private(i,j,r,dist,invdst,sf,sf0,sf1,screen_pot,screen_fld) &
     !$omp collapse(1) & 
     !$omp reduction(+:aceptor_np_int_re,aceptor_np_int_im) & 
     !$omp reduction(+:aceptor_np_int_re_mu,aceptor_np_int_im_mu)
!
!    aceptor
     do i = 1, aceptor%n_points_reduced
!       nanoparticle
        do j = 1, np%natoms
!
           r(1) = (xyz_np(1,j)-xyz_aceptor(1,i))
           r(2) = (xyz_np(2,j)-xyz_aceptor(2,i))
           r(3) = (xyz_np(3,j)-xyz_aceptor(3,i))
!
           dist = dsqrt(DOT_PRODUCT(r,r))
!
!          Skip when grid points are coincident to avoid instabilities
           if (dist.le.1.0e-14) go to 10
!
           invdst = one/dist
!
!          Screening function 
           sf  = dist / QMscrnFact
!
           sf0        = erf(sf)
           screen_pot = sf0
           if (np%dipoles) then
              sf1 = ( two * sf / sqrtpi ) * EXP(-(sf**2))
              screen_fld = sf0 - sf1
           endif
!
!          Integrate rho * q (imaginary charges)
!            --> the density has been already weigthed by the cube volume
!
           if (np%charges) then
              aceptor_np_int_re = aceptor_np_int_re +&
                                  rho_aceptor(i) * mm_q(j,1) * invdst * screen_pot
!
              aceptor_np_int_im = aceptor_np_int_im +&
                                  rho_aceptor(i) * mm_q(j,2) * invdst * screen_pot
           else if (np%dipoles) then 
!   
!             Charges
              aceptor_np_int_re = aceptor_np_int_re +&
                                  rho_aceptor(i) * mm_q(j,1) * invdst * screen_pot
!
              aceptor_np_int_im = aceptor_np_int_im +&
                                  rho_aceptor(i) * mm_q(j,2) * invdst * screen_pot
!
!             Dipoles
              aceptor_np_int_re_mu = aceptor_np_int_re_mu +&
                                  (-rho_aceptor(i) * mm_mu(1,j,1) * r(1) * (invdst**3) * screen_fld -&
                                    rho_aceptor(i) * mm_mu(2,j,1) * r(2) * (invdst**3) * screen_fld -&
                                    rho_aceptor(i) * mm_mu(3,j,1) * r(3) * (invdst**3) * screen_fld) 
!
              aceptor_np_int_im_mu = aceptor_np_int_im_mu +&
                                  (-rho_aceptor(i) * mm_mu(1,j,2) * r(1) * (invdst**3) * screen_fld -&
                                    rho_aceptor(i) * mm_mu(2,j,2) * r(2) * (invdst**3) * screen_fld -&
                                    rho_aceptor(i) * mm_mu(3,j,2) * r(3) * (invdst**3) * screen_fld)
!
           end if
!
           10 continue                  
!
        enddo
     enddo
     !$omp end parallel do
!
     integrals%aceptor_np_int(1) = aceptor_np_int_re + aceptor_np_int_re_mu
     integrals%aceptor_np_int(2) = aceptor_np_int_im + aceptor_np_int_im_mu
!
!
!    Deallocate and exit
!
     if (allocated(rho_aceptor)) deallocate(rho_aceptor)
     if (allocated(mm_q)       ) deallocate(mm_q)
     if (allocated(mm_mu)      ) deallocate(mm_mu)
     if (allocated(xyz_aceptor)) deallocate(xyz_aceptor)
     if (allocated(xyz_np)     ) deallocate(xyz_np)       
!
  end subroutine aceptor_nanoparticle_interaction_integral 
!----------------------------------------------------------------------
end module integrals_module


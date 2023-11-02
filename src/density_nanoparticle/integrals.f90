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
     integer :: vector_index
     integer  :: i,j,k,l,m,n
!
!    To allocate these quantities internally requires more memory but makes the calculation
!    in parallel faster. FORTRAN is slower accessing object-types. 
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
        do j = 1, donor%n_points_reduced
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
!          Aceptor-donor overlap 
           int_overlap = int_overlap +&
                         rho_aceptor(i) * rho_donor(j)
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
     integrals%aceptor_donor_overlap = -target_%omega_0 * int_overlap
!
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
     real(dp) :: sf,sf0,screen_pot !for screening
!
     integer  :: i,j,k,l
!
     integrals%aceptor_np_int = zero
!
     r   = zero
!
!    aceptor
     do i = 1, aceptor%nx
        !Write(*,'(4x,a6,i6,a7,i6)') ' Cycle ', i, ' out of ', aceptor%nx
        x_a = aceptor%xmin + aceptor%dx*(i-1)
        do j = 1, aceptor%ny
           y_a = aceptor%ymin + aceptor%dy*(j-1)
           do k = 1, aceptor%nz
              z_a = aceptor%zmin + aceptor%dz*(k-1)
!
!             nanoparticle
              do l = 1, np%natoms
!
                 r(1) = (x_a-np%xyz(1,l))
                 r(2) = (y_a-np%xyz(2,l))
                 r(3) = (z_a-np%xyz(3,l))
!
                 dist = dsqrt(DOT_PRODUCT(r,r))
!
!                Skip when grid points are coincident to avoid instabilities
                 if (dist.le.1.0e-14) then
                    go to 10
                 else
                    invdst = one/dist
                 endif
!
!                Screening function 
                 sf  = dist / QMscrnFact
!
                 sf0        = erf(sf)
                 screen_pot = sf0
!
!                Integrate rho * q (imaginary charges)
!                  --> the density has been already weigthed by the cube volume
!
!                WE HAVE TO UNDERSTAND IF THE DENSITY HAS THE PROPER SIGN
!
                 integrals%aceptor_np_int(1) = integrals%aceptor_np_int(1) +&
                                               aceptor%rho(i,j,k) * np%q(l,1) * invdst * screen_pot
!
                 integrals%aceptor_np_int(2) = integrals%aceptor_np_int(2) +&
                                               aceptor%rho(i,j,k) * np%q(l,2) * invdst * screen_pot
!
                 10 continue                  
!
              enddo
           enddo
        enddo
     enddo
!
  end subroutine aceptor_nanoparticle_interaction_integral 
!----------------------------------------------------------------------
end module integrals_module


!----------------------------------------------------------------------
module integrals_module
!      
!   Module integrals
!
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
     real(dp) :: sf,sf0,screen_pot !for screening
!
     real(dp), dimension(:), allocatable :: vector_results_coulomb, vector_results_overlap
!
     integer :: vector_index
     integer  :: i,j,k,l,m,n
!
     allocate(vector_results_coulomb(aceptor%nx*aceptor%ny*aceptor%nz*donor%nx*donor%ny*donor%nz),&
              vector_results_overlap(aceptor%nx*aceptor%ny*aceptor%nz*donor%nx*donor%ny*donor%nz))
!
     integrals%aceptor_donor_coulomb = zero
     integrals%aceptor_donor_overlap = zero
!
     vector_results_coulomb = zero
     vector_results_overlap = zero
!
     r   = zero
!
!!!$omp parallel do private(i,j,k,l,m,n,vector_index) &
!!!$omp collapse(6)
vector_index = 0
!    aceptor
     do i = 1, aceptor%nx
        !x_a = aceptor%xmin + aceptor%dx*(i-1)
        do j = 1, aceptor%ny
           !y_a = aceptor%ymin + aceptor%dy*(j-1)
           do k = 1, aceptor%nz
              !z_a = aceptor%zmin + aceptor%dz*(k-1)
!
!             donor
              do l = 1, donor%nx
                 !x_d = donor%xmin + donor%dx*(l-1)
                 do m = 1, donor%ny
                    !y_d = donor%ymin + donor%dy*(m-1)
                    do n = 1, donor%nz
                       !z_d = donor%zmin + donor%dz*(n-1)
!
                       !vector_index = vector_index + 1
!
                       !vector_index = (i-1)*aceptor%ny*aceptor%nz*donor%nx*donor%ny*donor%nz + &
                       !               (j-1)*aceptor%nz*donor%nx*donor%ny*donor%nz + &
                       !               (k-1)*donor%nx*donor%ny*donor%nz + &
                       !               (l-1)*donor%ny*donor%nz + &
                       !               (m-1)*donor%nz + &
                       !               n
!

                       r(1) = (aceptor%xyz(1,i,j,k)-donor%xyz(1,l,m,n))
                       r(2) = (aceptor%xyz(2,i,j,k)-donor%xyz(2,l,m,n))
                       r(3) = (aceptor%xyz(3,i,j,k)-donor%xyz(3,l,m,n))
!
                       dist = dsqrt(DOT_PRODUCT(r,r))
!
!                      Skip when grid points are coincident to avoid instabilities
                       if (dist.le.1.0e-14) then
                          !vector_results_coulomb(vector_index) = zero
                          !vector_results_overlap(vector_index) = zero
                          go to 10
!                      Skip when any density is zero
                       !elseif(abs(aceptor%rho(i,j,k)) .lt. 1.0E-15 .or. abs(donor%rho(l,m,n)) .lt. 1.0E-15) then
                       !   !vector_results_coulomb(vector_index) = zero
                       !   !vector_results_overlap(vector_index) = zero
                       !   go to 10
                       else
                          invdst = one/dist
                       endif
!
!                      Screening function 
                       sf  = dist / QMscrnFact

                       sf0        = erf(sf)
                       screen_pot = sf0
!
!                      Integrate charges as rho_aceptor * rho_donor * (1/dist) * screening
!                        --> the density has been already weigthed by the cube volume

                       !print *,vector_index
!
                       !vector_results_coulomb(vector_index) = aceptor%rho(i,j,k) * donor%rho(l,m,n) * invdst * screen_pot
!
                       !print *, vector_results_coulomb(vector_index)
!
                       !vector_results_overlap(vector_index) = aceptor%rho(i,j,k) * donor%rho(l,m,n)
!
                       !vector_results_overlap(vector_index) = 0.0d0
!
                       integrals%aceptor_donor_coulomb = integrals%aceptor_donor_coulomb +&
                                                         aceptor%rho(i,j,k) * donor%rho(l,m,n) * invdst * screen_pot
!
!                      Aceptor-donor overlap 
                       integrals%aceptor_donor_overlap = integrals%aceptor_donor_overlap +&
                                                         aceptor%rho(i,j,k) * donor%rho(l,m,n)
!
                       10 continue                  
!
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
!!!$omp end parallel do
!print *, 'sizes', size(vector_results_coulomb), size(vector_results_overlap)
!print *, vector_index
!stop
!
!    Sum all the contributions to the integrals
     !integrals%aceptor_donor_coulomb = SUM(vector_results_coulomb)
     !integrals%aceptor_donor_coulomb = zero
!
!    Weigth overlap by -omega_0
     !integrals%aceptor_donor_overlap = -target_%omega_0 * SUM(vector_results_overlap)
     !integrals%aceptor_donor_overlap = zero


!    Weight overlap by -omega_0
     integrals%aceptor_donor_overlap = -target_%omega_0 * integrals%aceptor_donor_overlap
!
     deallocate(vector_results_coulomb)
     deallocate(vector_results_overlap)
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


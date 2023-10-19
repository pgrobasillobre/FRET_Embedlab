!----------------------------------------------------------------------
module integrals_module
!      
!   Module integrals
!
    use density_module
!!    use output_module
!!    use target_module
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
     integer  :: i,j,k,l,m,n
!
     integrals%aceptor_donor_coulomb = zero
     integrals%aceptor_donor_overlap = zero
!
     r   = zero
!
!    aceptor
     do i = 1, aceptor%nx
        Write(*,'(4x,a6,i6,a7,i6)') ' Cycle ', i, ' out of ', aceptor%nx
        x_a = aceptor%xmin + aceptor%dx*(i-1)
        do j = 1, aceptor%ny
           y_a = aceptor%ymin + aceptor%dy*(j-1)
           do k = 1, aceptor%nz
              z_a = aceptor%zmin + aceptor%dz*(k-1)
!
!             donor
              do l = 1, donor%nx
                 x_d = donor%xmin + donor%dx*(l-1)
                 do m = 1, donor%ny
                    y_d = donor%ymin + donor%dy*(m-1)
                    do n = 1, donor%nz
                       z_d = donor%zmin + donor%dz*(n-1)
!
                       r(1) = (x_a-x_d)
                       r(2) = (y_a-y_d)
                       r(3) = (z_a-z_d)
!
                       dist = dsqrt(DOT_PRODUCT(r,r))
!
!                      Skip when grid points are coincident to avoid instabilities
                       if (dist.le.1.0e-14) then
                          go to 10
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
!
!    Weight overlap by -omega_0
     integrals%aceptor_donor_overlap = -target_%omega_0 * integrals%aceptor_donor_overlap
!
  end subroutine eet_aceptor_donor_integral 
!----------------------------------------------------------------------
end module integrals_module


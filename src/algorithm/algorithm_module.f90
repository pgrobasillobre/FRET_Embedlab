!----------------------------------------------------------------------
module algorithm_module
!      
!   Module for algorithm class
!
    use output_module
    use density_module
    use target_module
    use parameters_module
    use integrals_module
!
    Implicit None
!
!   Public Variables
!
    public algorithm
!
!   algorithm type
!
    type :: algorithm_type
!
       contains
!
       procedure :: integrate_density
       procedure :: eet_aceptor_donor
!
    end type algorithm_type
!    
    type (algorithm_type), target, save :: algorithm
!
   contains
!----------------------------------------------------------------------
   subroutine integrate_density(algorithm)
!
     implicit none
!     
!    subroutine for integrate density
!
     class(algorithm_type) :: algorithm
!
     type(density_type)    :: density_cube
!     
     call read_density(target_%density_file, density_cube)
!
     call int_density(density_cube)
!
     call delete_density(density_cube)
!
   end subroutine integrate_density
!----------------------------------------------------------------------
   subroutine eet_aceptor_donor(algorithm)
!
     implicit none
!     
!    subroutine for integrate density
!
     class(algorithm_type) :: algorithm
!
     type(density_type)    :: aceptor_density, donor_density
     type(integrals_type)  :: integrals_calc
!    
     call read_density(target_%aceptor_density, aceptor_density)
     call read_density(target_%donor_density, donor_density)
!
     call out_%print_density(target_%aceptor_density,aceptor_density%natoms,                 &
                             aceptor_density%nx,aceptor_density%ny,aceptor_density%nz,       &
                             aceptor_density%dx,aceptor_density%dy,aceptor_density%dz,       &
                             aceptor_density%xmin,aceptor_density%ymin,aceptor_density%zmin, &
                             aceptor_density%nelectrons,header='Aceptor Density Information')
!
     call out_%print_density(target_%donor_density,donor_density%natoms,                 &
                             donor_density%nx,donor_density%ny,donor_density%nz,       &
                             donor_density%dx,donor_density%dy,donor_density%dz,       &
                             donor_density%xmin,donor_density%ymin,donor_density%zmin, &
                             donor_density%nelectrons, header=' Donor Density Information')
!
     call eet_aceptor_donor_integral(integrals_calc,aceptor=aceptor_density,donor=donor_density)

!
!

     call delete_density(aceptor_density)
     call delete_density(donor_density)
!
   end subroutine eet_aceptor_donor
!----------------------------------------------------------------------
end module algorithm_module

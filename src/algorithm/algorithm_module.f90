!----------------------------------------------------------------------
module algorithm_module
!      
!   Module for algorithm class
!
    use output_module
    use density_module
    use nanoparticle_module
    use solvent_module
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
       integer :: n_threads_OMP
!
       contains
!
       procedure :: integrate_density
       procedure :: eet_aceptor_donor
       procedure :: aceptor_np_interaction
       procedure :: aceptor_np_donor
       procedure :: aceptor_solv_interaction
       procedure :: aceptor_solv_donor
       procedure :: aceptor_solv_np_donor
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
!    subroutine for aceptor-donor EET
!
     class(algorithm_type) :: algorithm
!
     type(density_type)    :: aceptor_density, donor_density
     type(integrals_type)  :: integrals_calc
!
!    Read input files
!
     call read_density(target_%aceptor_density, aceptor_density)
     call read_density(target_%donor_density, donor_density)
!
!    Print aceptor / donor density characteristics
!
     call out_%print_density(target_%aceptor_density,aceptor_density%natoms,aceptor_density%n_points_reduced, &
                             aceptor_density%nx,aceptor_density%ny,aceptor_density%nz,                        &
                             aceptor_density%dx,aceptor_density%dy,aceptor_density%dz,                        &
                             aceptor_density%xmin,aceptor_density%ymin,aceptor_density%zmin,                  &
                             aceptor_density%nelectrons,header='   Aceptor Density Information')
!
     call out_%print_density(target_%donor_density,donor_density%natoms,donor_density%n_points_reduced,       &
                             donor_density%nx,donor_density%ny,donor_density%nz,                              &
                             donor_density%dx,donor_density%dy,donor_density%dz,                              &
                             donor_density%xmin,donor_density%ymin,donor_density%zmin,                        &
                             donor_density%nelectrons, header='    Donor Density Information')
!
!    Calculate integrals
!
     call eet_aceptor_donor_integral(integrals_calc,aceptor=aceptor_density,donor=donor_density)
!
!    Print results and deallocate
!
     call out_%print_results_integrals(aceptor_donor_coulomb=integrals_calc%aceptor_donor_coulomb,&
                                       aceptor_donor_overlap=integrals_calc%aceptor_donor_overlap)
!
     call delete_density(aceptor_density)
     call delete_density(donor_density)
!
   end subroutine eet_aceptor_donor
!----------------------------------------------------------------------
   subroutine aceptor_np_interaction(algorithm)
!
     implicit none
!     
!    subroutine for aceptor-np interaction
!
     class(algorithm_type)    :: algorithm
!
     type(density_type)       :: aceptor_density
     type(nanoparticle_type)  :: nanoparticle
     type(integrals_type)     :: integrals_calc
!
!    Read input files
!
     call read_nanoparticle(target_%nanoparticle, nanoparticle)
     call read_density(target_%aceptor_density, aceptor_density)
!
!    Print NP / aceptor density characteristics
!
     call out_%print_nanoparticle(nanoparticle%natoms, nanoparticle%xyz,       &
                                  nanoparticle%charges,nanoparticle%dipoles)
!
     call out_%print_density(target_%aceptor_density,aceptor_density%natoms,aceptor_density%n_points_reduced, &
                             aceptor_density%nx,aceptor_density%ny,aceptor_density%nz,                        &
                             aceptor_density%dx,aceptor_density%dy,aceptor_density%dz,                        &
                             aceptor_density%xmin,aceptor_density%ymin,aceptor_density%zmin,                  &
                             aceptor_density%nelectrons,header='   Aceptor Density Information')
!
!    Calculate integrals
!
     call aceptor_nanoparticle_interaction_integral(integrals_calc,aceptor=aceptor_density,np=nanoparticle)
!
!    Print results and deallocate
!
     call out_%print_results_integrals(aceptor_np_int=integrals_calc%aceptor_np_int)
!
     call delete_nanoparticle(nanoparticle)
     call delete_density(aceptor_density)
!
   end subroutine aceptor_np_interaction
!----------------------------------------------------------------------
   subroutine aceptor_solv_interaction(algorithm)
!
     implicit none
!
!    subroutine for aceptor-solv interaction
!
     class(algorithm_type)    :: algorithm
!
     type(density_type)       :: aceptor_density
     type(solvent_type)       :: solvent
     type(integrals_type)     :: integrals_calc
!
!    Read input files
!
     call read_solvent(target_%solvent, solvent)
     call read_density(target_%aceptor_density, aceptor_density)
!
!    Print Solvent / aceptor density characteristics
!
     call out_%print_solvent(solvent%natoms,solvent%xyz,solvent%charges,solvent%dipoles)
!
     call out_%print_density(target_%aceptor_density,aceptor_density%natoms,aceptor_density%n_points_reduced, &
                             aceptor_density%nx,aceptor_density%ny,aceptor_density%nz,                        &
                             aceptor_density%dx,aceptor_density%dy,aceptor_density%dz,                        &
                             aceptor_density%xmin,aceptor_density%ymin,aceptor_density%zmin,                  &
                             aceptor_density%nelectrons,header='   Aceptor Density Information')
!
!    Calculate integrals
!
     call aceptor_solvent_interaction_integral(integrals_calc,aceptor=aceptor_density,solv=solvent)
!
     call out_%print_results_integrals(aceptor_solv_int=integrals_calc%aceptor_solv_int)
!
!    Print results and deallocate
!
     call delete_solvent(solvent)
     call delete_density(aceptor_density)
!
   end subroutine aceptor_solv_interaction
!----------------------------------------------------------------------
   subroutine aceptor_np_donor(algorithm)
!
     implicit none
!     
!    subroutine for aceptor-np-donor EET
!
     class(algorithm_type)    :: algorithm
!
     type(density_type)       :: aceptor_density, donor_density
     type(nanoparticle_type)  :: nanoparticle
     type(integrals_type)     :: integrals_calc
!
!    Read input files
!
     call read_nanoparticle(target_%nanoparticle, nanoparticle)
     call read_density(target_%aceptor_density, aceptor_density)
     call read_density(target_%donor_density, donor_density)
!
!    Print NP / densities characteristics
!
     call out_%print_nanoparticle(nanoparticle%natoms, nanoparticle%xyz,       &
                                  nanoparticle%charges,nanoparticle%dipoles)
!
     call out_%print_density(target_%aceptor_density,aceptor_density%natoms,aceptor_density%n_points_reduced, &
                             aceptor_density%nx,aceptor_density%ny,aceptor_density%nz,                        &
                             aceptor_density%dx,aceptor_density%dy,aceptor_density%dz,                        &
                             aceptor_density%xmin,aceptor_density%ymin,aceptor_density%zmin,                  &
                             aceptor_density%nelectrons,header='   Aceptor Density Information')
!
     call out_%print_density(target_%donor_density,donor_density%natoms,donor_density%n_points_reduced,       &
                             donor_density%nx,donor_density%ny,donor_density%nz,                              &
                             donor_density%dx,donor_density%dy,donor_density%dz,                              &
                             donor_density%xmin,donor_density%ymin,donor_density%zmin,                        &
                             donor_density%nelectrons, header='    Donor Density Information')
!
!    Calculate integrals
!
     call aceptor_nanoparticle_interaction_integral(integrals_calc,aceptor=aceptor_density,np=nanoparticle)
!
     call eet_aceptor_donor_integral(integrals_calc,aceptor=aceptor_density,donor=donor_density)
!
!    Print results and deallocate
!
     call out_%print_results_integrals(aceptor_donor_coulomb=integrals_calc%aceptor_donor_coulomb,&
                                       aceptor_donor_overlap=integrals_calc%aceptor_donor_overlap,&
                                       aceptor_np_int=integrals_calc%aceptor_np_int)
!
     call delete_nanoparticle(nanoparticle)
     call delete_density(aceptor_density)
     call delete_density(donor_density)
!
   end subroutine aceptor_np_donor
!----------------------------------------------------------------------
   subroutine aceptor_solv_donor(algorithm)
!
     implicit none
!
!    subroutine for aceptor-solvent-donor EET
!
     class(algorithm_type)    :: algorithm
!
     type(density_type)       :: aceptor_density, donor_density
     type(solvent_type)       :: solvent
     type(integrals_type)     :: integrals_calc
!
!    Read input files
!
     call read_solvent(target_%solvent, solvent)
     call read_density(target_%aceptor_density, aceptor_density)
     call read_density(target_%donor_density, donor_density)
!
!    Print solvent / densities characteristics
!
     call out_%print_solvent(solvent%natoms,solvent%xyz,solvent%charges,solvent%dipoles)
!
     call out_%print_density(target_%aceptor_density,aceptor_density%natoms,aceptor_density%n_points_reduced, &
                             aceptor_density%nx,aceptor_density%ny,aceptor_density%nz,                        &
                             aceptor_density%dx,aceptor_density%dy,aceptor_density%dz,                        &
                             aceptor_density%xmin,aceptor_density%ymin,aceptor_density%zmin,                  &
                             aceptor_density%nelectrons,header='   Aceptor Density Information')
!
     call out_%print_density(target_%donor_density,donor_density%natoms,donor_density%n_points_reduced, &
                             donor_density%nx,donor_density%ny,donor_density%nz,                        &
                             donor_density%dx,donor_density%dy,donor_density%dz,                        &
                             donor_density%xmin,donor_density%ymin,donor_density%zmin,                  &
                             donor_density%nelectrons, header='    Donor Density Information')
!
!    Calculate integrals
!
     !call aceptor_nanoparticle_interaction_integral(integrals_calc,aceptor=aceptor_density,np=nanoparticle)
!
     call eet_aceptor_donor_integral(integrals_calc,aceptor=aceptor_density,donor=donor_density)
!
!    Print results and deallocate
!
     call out_%print_results_integrals(aceptor_donor_coulomb=integrals_calc%aceptor_donor_coulomb,&
                                       aceptor_donor_overlap=integrals_calc%aceptor_donor_overlap,&
                                               aceptor_solv_int=integrals_calc%aceptor_solv_int)
!
     call delete_solvent(solvent)
     call delete_density(aceptor_density)
     call delete_density(donor_density)
!
   end subroutine aceptor_solv_donor
!----------------------------------------------------------------------
   subroutine aceptor_solv_np_donor(algorithm)
!
     implicit none
!
!    subroutine for aceptor-solvent-np-donor EET
!
     class(algorithm_type)    :: algorithm
!
     type(density_type)       :: aceptor_density, donor_density
     type(solvent_type)       :: solvent
     type(nanoparticle_type)  :: nanoparticle
     type(integrals_type)     :: integrals_calc
!
!    Read input files
!
     call read_solvent(target_%solvent, solvent)
     call read_nanoparticle(target_%nanoparticle, nanoparticle)
     call read_density(target_%aceptor_density, aceptor_density)
     call read_density(target_%donor_density, donor_density)
!
!    Print solvent / np / densities characteristics
!
     call out_%print_solvent(solvent%natoms,solvent%xyz,solvent%charges,solvent%dipoles)
!
     call out_%print_nanoparticle(nanoparticle%natoms,nanoparticle%xyz,nanoparticle%charges,nanoparticle%dipoles)
!
     call out_%print_density(target_%aceptor_density,aceptor_density%natoms,aceptor_density%n_points_reduced, &
                             aceptor_density%nx,aceptor_density%ny,aceptor_density%nz,                        &
                             aceptor_density%dx,aceptor_density%dy,aceptor_density%dz,                        &
                             aceptor_density%xmin,aceptor_density%ymin,aceptor_density%zmin,                  &
                             aceptor_density%nelectrons,header='   Aceptor Density Information')
!
     call out_%print_density(target_%donor_density,donor_density%natoms,donor_density%n_points_reduced, &
                             donor_density%nx,donor_density%ny,donor_density%nz,                        &
                             donor_density%dx,donor_density%dy,donor_density%dz,                        &
                             donor_density%xmin,donor_density%ymin,donor_density%zmin,                  &
                             donor_density%nelectrons, header='    Donor Density Information')
!
!    Calculate integrals
!
     !call aceptor_nanoparticle_interaction_integral(integrals_calc,aceptor=aceptor_density,np=nanoparticle)
!
     call eet_aceptor_donor_integral(integrals_calc,aceptor=aceptor_density,donor=donor_density)
!
!    Print results and deallocate
!
     call out_%print_results_integrals(aceptor_donor_coulomb=integrals_calc%aceptor_donor_coulomb,&
                                       aceptor_donor_overlap=integrals_calc%aceptor_donor_overlap,&
                                               aceptor_solv_int=integrals_calc%aceptor_solv_int,&
                                               aceptor_np_int=integrals_calc%aceptor_np_int)
!
     call delete_solvent(solvent)
     call delete_nanoparticle(nanoparticle)
     call delete_density(aceptor_density)
     call delete_density(donor_density)
!
   end subroutine aceptor_solv_np_donor
!----------------------------------------------------------------------
end module algorithm_module

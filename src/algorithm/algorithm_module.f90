!----------------------------------------------------------------------
module algorithm_module
!      
!   Module for algorithm class
!
    use output_module
    use density_module
    use nanoparticle_module
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
     real(dp)              :: integral
!     
     call read_density(target_%density_file, density_cube, rotation=.false.)
!
     call int_density(density_cube,integral)
!
     call out_%print_density(cube_file=target_%density_file,cube=density_cube,integral=integral)
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
     call read_density(target_%aceptor_density, aceptor_density, rotation=target_%rotate_aceptor, what_dens='aceptor')
     call read_density(target_%donor_density, donor_density, rotation=target_%rotate_donor, what_dens='donor')
!
!    Print aceptor / donor density characteristics
!
     call out_%print_density(target_%aceptor_density,aceptor_density,header=aceptor_header)
     call out_%print_density(target_%donor_density,donor_density,header=donor_header)
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
     call read_density(target_%aceptor_density, aceptor_density, rotation=target_%rotate_aceptor, what_dens='aceptor')
     call read_nanoparticle(target_%nanoparticle, nanoparticle, rotation=target_%rotate_np)
!
!    Print NP / aceptor density characteristics
!
     call out_%print_nanoparticle(nanoparticle%natoms, nanoparticle%xyz,       &
                                  nanoparticle%charges,nanoparticle%dipoles)
!
     call out_%print_density(target_%aceptor_density,aceptor_density,header=aceptor_header)
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
     call read_density(target_%aceptor_density, aceptor_density, rotation=target_%rotate_aceptor, what_dens='aceptor')
     call read_density(target_%donor_density, donor_density, rotation=target_%rotate_donor, what_dens='donor')
     call read_nanoparticle(target_%nanoparticle, nanoparticle, rotation=target_%rotate_np)
!
!    Print NP / densities characteristics
!
     call out_%print_nanoparticle(nanoparticle%natoms, nanoparticle%xyz,       &
                                  nanoparticle%charges,nanoparticle%dipoles)
!
     call out_%print_density(target_%aceptor_density,aceptor_density,header=aceptor_header)
     call out_%print_density(target_%donor_density,donor_density,header=donor_header)
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
end module algorithm_module

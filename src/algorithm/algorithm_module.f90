!----------------------------------------------------------------------
module algorithm_module
!      
!   Module for algorithm class
!
    !!use output_module
    use target_module
    use parameters_module
    use density_module
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
     call read_density(target_%density_file, density)
!
     call int_density(density)
!
!     call print_density_integral(density)
!




   end subroutine integrate_density
!----------------------------------------------------------------------
end module algorithm_module

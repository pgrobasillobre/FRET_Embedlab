!----------------------------------------------------------------------
module parallel_module
!      
!   Module for parallel class
!
    Implicit None
!
!   Public Variables
!
    public parallel
!
!   parallel type
!
    type :: parallel_type
!
       integer  :: n_threads_OMP
!        
       contains
!
    end type parallel_type
!    
    type (parallel_type), target, save :: parallel
!
end module parallel_module

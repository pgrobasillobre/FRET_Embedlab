!----------------------------------------------------------------------
module parameters_module
!      
!   Module Parameters
!
    Implicit None
    public 
!
    integer, parameter                       :: dp=kind(1.0d0)
!    
    real(dp), parameter                      :: zero   = 0.0d0
    real(dp), parameter                      :: one    = 1.0d0
    real(dp), parameter                      :: two    = 2.0d0
    real(dp), parameter                      :: three  = 3.0d0
    real(dp), parameter                      :: four   = 4.0d0
    real(dp), parameter                      :: five   = 5.0d0
    real(dp), parameter                      :: six    = 6.0d0
    real(dp), parameter                      :: seven  = 7.0d0
    real(dp), parameter                      :: eight  = 8.0d0
!    
    real(dp), parameter                      :: pi     = four*atan(one)
!    
    real(dp), parameter                      :: ToBohr = 1.8897261254578281d0
    real(dp), parameter                      :: ToAng  = 1.0d0/ToBohr
    real(dp), parameter                      :: Light  = 137.035227d0
!
    real(dp), parameter                      :: Half   = 0.5d0
    real(dp), parameter                      :: FSmAu  = 0.00000021739d0
!
!   nearest neighbour distances
!
    real(dp), parameter                      :: RAg0 = 2.885d0 * ToBohr
    real(dp), parameter                      :: RAu0 = 2.885d0 * ToBohr
    real(dp), parameter                      :: RNa0 = 6.9226100144600462d0 ! 3.66328d0 * ToBohr
    real(dp), parameter                      :: RC0  = 1.418d0 * ToBohr
    real(dp), parameter                      :: ROH0 = 1.65d0
    real(dp), parameter                      :: RAl0 = 2.861d0 * ToBohr
    real(dp), parameter                      :: RCu0 = 2.543463d0 * ToBohr
!
    real(dp), parameter                      :: fermi_velocity = 0.4573138778d0
!
!-----------------------------------------------------------------------
end module parameters_module

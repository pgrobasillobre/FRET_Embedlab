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
    real(dp), parameter                      :: QMscrnFact = 0.2d0
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
!
    real(dp), parameter                      :: ncellmax = 1000000 
!
    character(len=200), parameter            :: fret_start = '# fret quantities ------------------------#'
    character(len=200), parameter            :: fret_end   = '# end fret quantities ------------------------'
    character(len=200), parameter            :: charges_header ='#        q_re                     q_im                    coords_x&
                                                                           &                 coords_y                 coords_z'
!
!-----------------------------------------------------------------------
end module parameters_module

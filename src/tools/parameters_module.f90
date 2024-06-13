!----------------------------------------------------------------------
module parameters_module
!      
!   Module Parameters
!
    Implicit None
    public 
!
    integer, parameter                     :: dp=kind(1.0d0)
!    
    real(dp), parameter                    :: zero   = 0.0d0
    real(dp), parameter                    :: one    = 1.0d0
    real(dp), parameter                    :: two    = 2.0d0
    real(dp), parameter                    :: three  = 3.0d0
    real(dp), parameter                    :: four   = 4.0d0
    real(dp), parameter                    :: five   = 5.0d0
    real(dp), parameter                    :: six    = 6.0d0
    real(dp), parameter                    :: seven  = 7.0d0
    real(dp), parameter                    :: eight  = 8.0d0
!
    real(dp), parameter                    :: QMscrnFact = 0.2d0
!    
    real(dp), parameter                    :: pi         = four*atan(one)
    real(dp), parameter                    :: sqrtpi     = 1.7724538509055159d0
    real(dp), parameter                    :: to_radians = pi/180.0d0
    real(dp), parameter                    :: to_degrees = one/to_radians
!    
    real(dp), parameter                    :: angle_thresh = 0.1d0*to_radians ! = 6 degrees
!
    real(dp), parameter                    :: ToBohr = 1.8897261254578281d0
    real(dp), parameter                    :: ToAng  = 1.0d0/ToBohr
    real(dp), parameter                    :: Light  = 137.035227d0
!
    real(dp), parameter                    :: Half   = 0.5d0
    real(dp), parameter                    :: FSmAu  = 0.00000021739d0
!
    real(dp), parameter                    :: ncellmax = 10000000 ! Fretty has 1.0E06, we have 1.0E07 
!
    integer, parameter                     :: iuout = 12
!
    character(len=30),  parameter          :: aceptor_header = '   Aceptor Density Information' 
    character(len=29),  parameter          :: donor_header   = '    Donor Density Information'
!
    character(len=200), parameter          :: fret_start = '# fret quantities ------------------------#'
    character(len=200), parameter          :: fret_end   = '# end fret quantities ------------------------'
    character(len=200), parameter          :: charges_header ='#        q_re                     q_im                coords_x&
                                                                         &                 coords_y                 coords_z'
    character(len=200), parameter          :: dipoles_header ='# q_re,   q_im,   mu_re_x,   mu_re_y,   mu_re_z,  mu_im_x,   &
                                                                         &mu_im_y,   mu_im_z, coords_x, coords_y, coords_z'
!
!-----------------------------------------------------------------------
end module parameters_module

module constants

  implicit none

  public

  ! stefan: Stefan-Boltzmann constant (usually denoted by sigma)
  ! R0: Estimate of the radius of the earth from Google.
  ! T: Mean Solar Day (seconds)
  ! OMEGA: rotation speed of earth of earth (rad/sec)
  double precision :: PI, TWOPI, PIBY2, SECDAY, DAYSYR, SECSYR, STEFAN, SIGOV2
  double precision :: R0, T, OMEGA

  parameter (PI = 3.14159265358979324D0, &
       TWOPI = 2*PI, &
       PIBY2 = 0.5*PI, &
       SECDAY = 86400.0d0, &
       DAYSYR = 365.0d0, &
       SECSYR = SECDAY*DAYSYR, &
       STEFAN = 5.67040D-8, &
       SIGOV2 = 0.5d0*STEFAN, &
       R0 = 6378.1e3, &
       T = (23*60 + 56)*60 + 4.098903691d0,  &
       OMEGA = 2*PI/T &    
       )

contains

end module constants


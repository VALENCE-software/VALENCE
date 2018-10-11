
!>    \mainpage Variational Subspace Valence Bond (VSVB) method
!!     "The variational subspace valence bond method",
!!     G. D. Fletcher, J. Chem. Phys., 142, 134112 (2015).
!!    \n
!!    'Orbital Basis Set' (OBS) version recomputes the
!!    super-contracted integrals as needed or stores them
!!    in aggregate memory

program      valence_driver
  use valence_init
  use valence_finit
  implicit none
! the values are dummy variables, since the output is printed once to the screen
  integer d
  real(kind(0.d0)) e
  call valence_initialize
  call calculate_vsvb_energy( e )
  call valence_finalize
!  call valence( d, e, .false., .true.)

end program valence_driver

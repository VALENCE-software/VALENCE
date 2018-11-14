
!> initializes VSVB info
!! \param info [out] : set to 0 if initialization occurs successfully
subroutine init(info)
  integer info
  call valence_api_initialize(info)
end subroutine init

!> returns the number of atoms
!! \param n [out] : number of atoms
subroutine getn(n)
  use tools, only: dp
  use integrals
  implicit none
  integer n
  ! get number of atoms from integrals module
  n = natom
end subroutine getn

!> returns the VSVB energy for given coordinates
!! \param x [in] : cartesian coordinates to get the VSVB energy of (angstroms)
!! \param v [out] : VSVB energy for input coordinates ( cm^{-1} )
subroutine calcsurface(x,v)
  use tools, only: dp
  real(dp), dimension(*) :: x
  real(dp) v
  call valence_api_calculate_energy( x, v )
  ! convert from hartrees to cm^{-1}
  v = v*219474.631_dp

end subroutine calcsurface

subroutine finalize
  call valence_api_finalize()
end subroutine finalize

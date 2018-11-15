
!> initializes VSVB info
!! \param info [out] : set to 0 if initialization occurs successfully
!! \param call_mpi_init [in] : flag to call mpi_init or not
!!                if it's 0, do not call mpi_init.
!!                if it's 1, call mpi_init.
!! \param comm [in] : MPI communicator to be used if call_mpi_init == 0
!!                and MPI is being used
subroutine valence_api_initialize( info, call_mpi_init, comm )
  use tools, only: dp
  use integrals
  use valence_init
  implicit none
#ifdef VALENCE_MPI
  include "mpif.h"
  integer ierr, status
#endif
  integer info, call_mpi_init, comm
  ! call initialization of VSVB code

#ifdef VALENCE_MPI
  if( call_mpi_init .eq. 0 ) then
    call valence_initialize( comm )
  else
! passing no argument means that MPI_Init will be called internally
     call valence_initialize( )
  endif
#else
    call valence_initialize( )
#endif
  info = 0
end subroutine valence_api_initialize

!> returns the VSVB energy for given coordinates
!! \param x [in] : cartesian coordinates to get the VSVB energy of (angstroms)
!! \param v [out] : VSVB energy for input coordinates ( Hartrees )
subroutine valence_api_calculate_energy(x,v)
  use tools, only: dp, angs2bohr
  use integrals
  use xm
  implicit none
#ifdef VALENCE_MPI
  include "mpif.h"
#endif
  real(dp), dimension(*) :: x
  real(dp) v
  integer i,k,j
#ifdef VALENCE_NITROGEN_PARALLEL
  real(dp) energy, etol
  real(dp)      zero,         ten,            rln10
  parameter  ( zero = 0.0_dp, ten = 10.0_dp, rln10=2.30258_dp )
  integer :: myiostat, ns, np, con_length, mnshi, mxshi
  character(100) script
  logical, save :: first_time=.true.
#endif

  k = 0
  do    i  =  1,  natom
     do    j  =  1,  3
        k = k + 1
        coords( j, i )  =  x( k )
     end   do
  end   do

  !  update coordinates
#ifdef VALENCE_NITROGEN_PARALLEL

  call xm_output_dimensions_tolerances( 'input', .false. )

  write(107,*)
  write(107,*)

  call xm_output_coords_basis( 'input', .true. )

  ! add orbitals to the new_input file
  if( first_time ) then
     first_time = .false.
     call xm_output( '', energy, etol, 'input',.true. )
     call system( "mv input new_input" )
  else
     call system( "cat input orbitals > new_input" )
  endif

  ! call the script at location $VALENCE_SCRIPT to run VALENCE with MPI
  call getenv( "VALENCE_SCRIPT", script)
  call system( script )

  ! pull out the energy
  open(unit=108, file='energy_output', status='UNKNOWN', action='readwrite', &
       position='rewind', iostat=myiostat)

  read( 108, "(f19.15)") v
  close(unit=108)
#else

  call angs2bohr(natom,coords)
  call calculate_vsvb_energy( v )

#endif

#ifdef VALENCE_MPI
! passing an argument means that MPI_Finalize will not be called internally
  call xm_end( 1 )
#else
  call xm_end( )
#endif

end subroutine valence_api_calculate_energy

!> finalizes VSVB run
!! \param call_mpi_finalize [in] : flag to call mpi_finalize or not
!!                if it's 0, do not call mpi_finalize.
!!                if it's 1, call mpi_finalize.
subroutine valence_api_finalize( call_mpi_finalize )
  use valence_finit
  implicit none
  integer call_mpi_finalize
#ifdef VALENCE_MPI
  include "mpif.h"
  integer status
  if( call_mpi_finalize .eq. 0 ) then
! passing an argument means that MPI_Finalize will not be called internally
      call valence_finalize( 1 )
  else
! passing no argument means that MPI_Finalize will be called internally
      call valence_finalize(  )
  endif
#else

  call valence_finalize( )

#endif

end subroutine valence_api_finalize

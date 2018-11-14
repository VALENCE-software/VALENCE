
!> initializes VSVB info
!! \param info [out] : set to 0 if initialization occurs successfully
subroutine valence_api_initialize(info)
  use tools, only: dp
  use integrals
  use valence_init
  implicit none
#ifdef VALENCE_MPI
  include "mpif.h"
  integer ierr, mpi_new_comm,status
#endif
  integer info
  ! call initialization of VSVB code
#ifdef VALENCE_MPI
  call mpi_init(status)
!  call mpi_comm_dup( mpi_comm_world, mpi_new_comm, ierr)
!  call valence_initialize( mpi_new_comm )
  call valence_initialize( mpi_comm_world )
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
  call xm_end( mpi_comm_world )
#else
  call xm_end( )
#endif

end subroutine valence_api_calculate_energy

!> finalizes VSVB run
subroutine valence_api_finalize
  use valence_finit
  implicit none
#ifdef VALENCE_MPI
  include "mpif.h"
  integer status
  call valence_finalize( mpi_comm_world )
  call mpi_finalize( status )
#else
  call valence_finalize( )
#endif
end subroutine valence_api_finalize

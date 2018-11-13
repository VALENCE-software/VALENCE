
!> initializes VSVB info
!! \param info [out] : set to 0 if initialization occurs successfully
subroutine init(info)
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
#ifdef VALENCE_NITROGEN_READ
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
#ifdef VALENCE_NITROGEN_READ
  open(unit=107, file='input', status='UNKNOWN', action='readwrite', &
       position='rewind', iostat=myiostat)


  write(107,"(15i4)") natom,natom_t, npair,nunpd,ndocc,  &
         totlen,xpmax, nspinc, num_sh,num_pr,nang, ndf,nset,nxorb, mxctr

  write(107,*)

  write(107,"(6i4, f16.12, f16.12 )",advance="no")  int(ctol/rln10), int(-log(dtol)/log(10.0_dp)), int(-log(itol)/log(10.0_dp)), &
       ntol_e_min, ntol_e_max, max_iter, ptbnmax,feather

  do j=1,nset
     write(107,"(2i4)",advance="no") orbset( 1, j ), orbset( 2, j )
  enddo

     write(107,*)
     write(107,*)

  do j=1,natom
        write(107,"(i4,' ')",advance="no") atom_t(j)
     do i=1,3
        write(107,"(f16.12,' ')",advance="no")  coords(i,j)
        write(*,"(f16.12,' ')",advance="no")  coords(i,j)
     end do
     write(107,*)
     write(*,*)
  end do

  write(107,*)
! basis set

     ns = 1
     np = 1
       do    i  =  1,  natom_t
          write(107,"(f16.12,i4)")  nuc_charge( i ),num_shell_atom( i )
          do    j  =  1,  num_shell_atom( i )
             con_length = map_shell2prim( ns + 1 ) - map_shell2prim( ns )
             write(107,"(2i4)")  ang_mom( ns ), con_length

             !     avoid redundant input of unit weight for uncontracted GTO

             if ( con_length .eq. 1 ) then
                write(107,"(f18.12)")  exponent( np )
                np = np + 1
             else
                do k  =  1,  con_length
                   write(107,"(2f18.12)")  exponent( np ), con_coeff( np )
                   np = np + 1
                end   do
             end if

             ns = ns + 1
          end   do
       end   do

     write(107,*)

  close(unit=107)

! add orbitals to the new_input file
  if( first_time ) then
     first_time = .false.
     call xm_output( '', energy, etol, 'input',.true. )
     call system( "mv input new_input" )
  else
     call system( "cat input orbitals > new_input" )
  endif


! call the script at location $VALENCE_SCRIPT to run VALENCE
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


! convert from hartrees to cm^{-1}
  v = v*219474.631_dp
  
end subroutine calcsurface

subroutine finalize
  use valence_finit
  implicit none
#ifdef VALENCE_MPI
  include "mpif.h"
  integer status
  call mpi_finalize( status )
  call valence_finalize( mpi_comm_world )
#else
  call valence_finalize( )
#endif
end subroutine finalize


!> initializes VSVB info
!! \param info [out] : set to 0 if initialization occurs successfully
subroutine init(info)
  use tools, only: dp
  use integrals
  use valence_init
  implicit none
  include "mpif.h"
  integer info,i, status, ierr, mpi_new_comm
  integer d
  real(dp) e
  ! call initialization of VSVB code
!  call valence( d, e, .false., .false. )
  call mpi_init(status)
  call mpi_comm_dup( mpi_comm_world, mpi_new_comm, ierr)
  call valence_initialize( mpi_new_comm )
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
  include "mpif.h"
  real(dp), dimension(*) :: x
  real(dp) v
  integer i,k
  integer ::  j, myiostat

  k = 0
  do    i  =  1,  natom
     do    j  =  1,  3
        k = k + 1
        coords( j, i )  =  x( k )
     end   do
  end   do
     
!  update coordinates somehow
  ! open(unit=107, file='coords1', status='UNKNOWN', action='readwrite', &
  !      position='rewind', iostat=myiostat)

  ! do j=1,natom
  !       write(107,"(i4,' ')",advance="no") atom_t(j)
  !    do i=1,3
  !       write(107,"(f16.12,' ')",advance="no")  coords(i,j)
  !       write(*,"(f16.12,' ')",advance="no")  coords(i,j)
  !    end do
  !    write(107,*)
  !    write(*,*)
  ! end do

  ! close(unit=107)
  ! call system( "../examples/script" )

  ! open(unit=108, file='fff', status='UNKNOWN', action='readwrite', &
  !      position='rewind', iostat=myiostat)

  ! read( 108, "(f19.15)") v
  ! close(unit=108)


! if orbital file exists, read from it. otherwise, don't., 


  call angs2bohr(natom,coords)
  call calculate_vsvb_energy( v )
  call xm_end( mpi_comm_world )

! convert from hartrees to cm^{-1}
  v = v*219474.631_dp
  
end subroutine calcsurface

subroutine finalize
  use valence_finit
  implicit none
  include "mpif.h"
  integer status
  call mpi_finalize( status )
  call valence_finalize( mpi_comm_world )
end subroutine finalize

module valence_init

contains
subroutine valence_initialize( comm )
  use tools, only: angs2bohr
  use          densitywork
  use          integrals
  use xm
  implicit     none
  integer, intent(in), optional :: comm

  !     local variables
! this handles all reading of input, allocating for input arrays, and MPI initialization
! these five calls should only be done once
  if( present( comm ) ) then
     call xm_propagate( comm )
  else
     call xm_propagate
  endif
  call xm_print(   'title', 'variational subspace valence bond;' )
  call xm_print( 'comment', 'written by G. Fletcher;' )
  call read_allocate_input
  !     scale the geometry to bohrs for the energy in au
  call angs2bohr(natom,coords)

#ifndef SIMINT_INT
  call rys_init
#endif

end subroutine valence_initialize


subroutine read_allocate_input
  use          densitywork
  use          integrals
  use tools, only: dp
  use xm
  implicit     none

  integer mxctr_in, nset_in,max_iter_in
  integer      ntol_c,ntol_d,ntol_i, ntol_e_min_in,ntol_e_max_in
  integer i,k,j,ierr,mnshi,mxshi
  real(dp)      zero,         ten,            rln10
  parameter  ( zero = 0.0_dp, ten = 10.0_dp, rln10=2.30258_dp )
  character(len=1000) :: input_file

! read in the input file
#ifdef FILE_IN 
 call get_command_argument(1, input_file)
  if( len_trim(input_file) == 0 ) call xm_abort('must have one input file')
  open( unit=100, file=input_file )
! needs to be closed at some point
!  close( 100 )
#endif
  call xm_getdims( natom,natom_t, npair,nunpd,ndocc, totlen,  &
       xpmax,nspinc, num_sh,num_pr,nang, ndf,nset_in,nxorb, mxctr_in )
  if ( npair .gt. 0 .and. nspinc .lt. 1 ) call xm_abort('no spin couplings;')

  mxctr = mxctr_in
  nset = nset_in

  allocate( atom_t(    natom + xpmax ), stat = ierr )
  allocate( coords( 3, natom + xpmax ), stat = ierr )
  allocate( map_atom2shell( natom_t + xpmax ),  stat = ierr )
  allocate( num_shell_atom( natom_t + xpmax ),  stat = ierr )
  allocate( map_shell2prim( num_sh  + 1 ),  stat = ierr )
  allocate(        ang_mom( num_sh      ),  stat = ierr )
  allocate(     nuc_charge( natom_t     ),  stat = ierr )
  allocate(       exponent( num_pr      ),  stat = ierr )
  allocate(      con_coeff( num_pr      ),  stat = ierr )
  allocate( unnormalized_con_coeff( num_pr      ),  stat = ierr )

  spinopt  =  .false.
  !!  uncomment this for Blue Gene /Q
! if  (  nspinc .eq. 0  )  then
!    allocate(  pair_sc( npair, 2, 1 ),  stat = ierr )
!    allocate( coeff_sc(           1 ),  stat = ierr )
! else
     allocate(  pair_sc( npair, 2, nspinc ),  stat = ierr )
     allocate( coeff_sc(           nspinc ),  stat = ierr )
! end if

  allocate( xorb( nxorb ),  stat = ierr )
  allocate( root( nxorb ),  stat = ierr )
  allocate( orbset( 2, nset ),  stat = ierr )


  nelec  = 2*npair + 2*ndocc + nunpd
  norbs  = 2*npair + ndocc + nunpd + ndf
  nalpha =   npair + nunpd + ndocc
  nbeta  =   npair + ndocc

  !     strictly, +xpmax only needed if doing 1st-order optimization

  allocate(    map_orbs( norbs + xpmax + 1 ),  stat = ierr )
  allocate( orbas_atnum( norbs + xpmax ),  stat = ierr )
  allocate( orbas_atset( mxctr, norbs + xpmax ),  stat = ierr )
  allocate( xpset( totlen + xpmax ),  stat = ierr )
  allocate( coeff( totlen + xpmax ),  stat = ierr )


  call xm_input( ntol_c, ntol_e_min_in, ntol_e_max_in, ntol_d, ntol_i, &
       max_iter_in )
  ntol_e_min = ntol_e_min_in
  ntol_e_max = ntol_e_max_in
  max_iter = max_iter_in

  !     set tolerances from input

  ctol = rln10*dble( ntol_c  )
  dtol =     ten**( -ntol_d  )
  itol =     ten**( -ntol_i )

end subroutine read_allocate_input
end module valence_init

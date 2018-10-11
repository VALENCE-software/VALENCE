
!>    \mainpage Variational Subspace Valence Bond (VSVB) method
!!     "The variational subspace valence bond method",
!!     G. D. Fletcher, J. Chem. Phys., 142, 134112 (2015).
!!    \n
!!    'Orbital Basis Set' (OBS) version recomputes the
!!    super-contracted integrals as needed or stores them
!!    in aggregate memory

subroutine      valence( energy )
  use tools, only: dp
  use valence_init
  use valence_finit
  implicit     none
  real(dp) energy

! this handles all reading of input, allocating for input arrays, and MPI initialization
  call valence_initialize

! yes, this is overly broad, but I wanted to separate the initialization/input reading,
! finalization, and the rest of the code so that it's easier to call this from a library
  call calculate_vsvb_energy( energy )

! handles deallocation of what was allocated in valence_initialize
  call valence_finalize

end subroutine valence


subroutine      calculate_vsvb_energy( energy )
  use          densitywork
  use          integrals
  use xm
  use timing_flops
  use tools, only: dp,get_nuclear_repulsion_energy, angs2bohr
#ifdef SIMINT_INT
  use valence_simint, only: valence_initialize_simint, valence_finalize_simint
#endif
  implicit     none
#ifdef PRINT_TIMING
  include "mpif.h"
#endif

  !     local variables

  real(dp),  allocatable ::   w( : , : )
  real(dp),  allocatable :: eig( : )
  real(dp),  allocatable ::  v1( : )
  real(dp),  allocatable ::  v2( : )
  logical,  allocatable ::  coefflock( : )

  integer      ierr, i,j,k,it,n, mnshi,mxshi
  integer      nao_type, mxcf2
  real(dp) zij,rsq
  real(dp)      tokcal
  parameter  ( tokcal = 627.509469_dp )

  !     optimization control
  !     etol = total energy convergence tolerance

  integer      norbas,norbz
  integer      ntol_c,ntol_d,ntol_i
  real(dp)      energy,eprev,eprv_sc,eprv_orb
  real(dp)      etol,  relaxn
  real(dp)      zero,         ten,            rln10
  parameter  ( zero = 0.0_dp, ten = 10.0_dp, rln10=2.30258_dp )
  integer     nproc,myrank,master      ! for Nstore calc.

  !     generic print buffers

  integer   int_out(2)    !  lengths not
  real(dp)   dbl_out(2)    !   protected

! ! input parameters
!   real(8) external_coords(*)
!   logical use_external_coords

  ! if( use_external_coords ) then
  !    k = 0
  !    do    i  =  1,  natom
  !       do    j  =  1,  3
  !          k = k + 1
  !          coords( j, i )  =  external_coords( k )
  !       end   do
  !    end   do
  !    call angs2bohr(natom,coords)
  ! endif

#ifdef SIMINT_INT
     ! initialize simint, based on basis set that was read in
     call valence_initialize_simint
#endif

  !     compute the nuclear repulsion energy

  call get_nuclear_repulsion_energy(natom, atom_t, nuc_charge, coords, enucrep)

  dbl_out( 1 ) = enucrep
  if ( natom .gt. 1 ) call xm_print(  &
       'parameter', 'nuclear repulsion;', int_out, dbl_out )


  !     initialize the ao integrals

  nao_type  =  0
  do   i  =  0,   nang
     nao_type  =  nao_type  +  ( ( i + 1 )*( i + 2 ) )/2
  end  do

  allocate(    nxyz( 3, nao_type ),  stat = ierr )
  allocate(    angn(    nao_type ),  stat = ierr )
  allocate(     ashl( 0:nang ),  stat = ierr )
  allocate(     ashi( 0:nang ),  stat = ierr )

  call  cartesian( nang, nxyz )
  call  setangn()

  mxcf2  = ( ( ( nang + 1 )*( nang + 2 ) )/2 )**2

  allocate(    dij( mxcf2 ),  stat = ierr )
  allocate(    dkl( mxcf2 ),  stat = ierr )

  !     get the maximum OBS size

  max_obs  =  0
  do  i  =  1,  norbs
     n  =  0
     do  j  =  1,   orbas_atnum( i )
        it  =  atom_t( orbas_atset( j, i ) )
        mnshi =         map_atom2shell( it )
        mxshi = mnshi + num_shell_atom( it ) - 1
        do  k  =  mnshi, mxshi
           n  =  n  +  shell_size( ang_mom( k ) )
        end  do
     end  do
     max_obs  =  max( max_obs, n )
  end  do

  allocate( coeffi( max_obs ),  stat = ierr )
  allocate( coeffj( max_obs ),  stat = ierr )
  allocate( coeffk( max_obs ),  stat = ierr )
  allocate( coeffl( max_obs ),  stat = ierr )

  allocate( atom_ndf( 2, mxctr ),  stat = ierr )
  allocate(  ndf2orb( mxctr ),  stat = ierr )
  allocate(    xpnew( xpmax ),  stat = ierr )


  !     normalize the input orbitals, NDF's first

  call  normal( 2*npair + nunpd + ndocc + 1, norbs )
  call  normal( 1, norbs - ndf )

  allocate(    bra_a( nalpha ),  stat = ierr )
  allocate(    ket_a( nalpha ),  stat = ierr )
  allocate(    bra_b( nbeta  ),  stat = ierr )
  allocate(    ket_b( nbeta  ),  stat = ierr )
  allocate(    bexch( npair  ),  stat = ierr )
  allocate(    kexch( npair  ),  stat = ierr )
  allocate(       bra( nelec ),  stat = ierr )
  allocate(       ket( nelec ),  stat = ierr )
  allocate(  wdet(  nelec, nelec  ),  stat = ierr )
#ifdef GIVENS_BGQ
  if( mod(nalpha, 4) .gt. 0 ) then
     padded_size = nalpha+(4-mod(nalpha, 4))
  else
#elif GIVENS_KNL
  if( mod(nalpha, 8) .gt. 0 ) then
     padded_size = nalpha+(8-mod(nalpha, 8))
  else
#endif
     padded_size = nalpha
#if defined(GIVENS_BGQ) || defined(GIVENS_KNL)
  endif
#endif
  allocate(  aket( padded_size, nalpha ),  stat = ierr )
  allocate(  bket( padded_size, nalpha ),  stat = ierr )

  allocate(  abra_npair( npair, nelec  ),  stat = ierr )
  allocate(  bbra_npair( npair, nelec  ),  stat = ierr )
  allocate(  abra_docc_un( ndocc+nunpd, npair*2  ),  stat = ierr )
  allocate(  bbra_docc_un( ndocc+nunpd, npair*2  ),  stat = ierr )
  allocate(  aket_docc_un( ndocc+nunpd, ndocc+nunpd ),  stat = ierr )
  allocate(  bket_docc_un( ndocc+nunpd, ndocc+nunpd ),  stat = ierr )

  norbz = 2*npair + ndocc + nunpd + 1
  allocate( schwarz(  (norbz**2 + norbz)/2  ), stat = ierr )

  allocate(  ipvt( padded_size ),  stat = ierr )

  !     compute the guess energy

  call  guess_energy( energy )
#ifdef PRINT_TIMING
  guess_time = MPI_Wtime()
#endif
  dbl_out( 1 ) = energy
  call xm_print( 'parameter', 'guess energy;', int_out, dbl_out )
  call xm_output( 'save', energy,etol )


  !     max_iter=1 with nset=0 can be used to do just the spin optimization

  if  (  max_iter .gt. 0  ) then
     call xm_print( 'header', 'orbital optimization;' )
     etol  =  ten**( -ntol_e_max )/tokcal
     dem_gs  =  ptbnmax .gt. zero  .and.  nxorb .eq. 0


     !     allocate memory for storing integrals

     if  (  dem_gs  )  then
        norbz = 2*npair + ndocc + nunpd
     else
        norbz = 2*npair + ndocc + nunpd + 1
     end if

!    nstore = min( 2 + (((norbz**2+norbz)/2)*norbz**2)/nrank,2000000000/nrank )
     nstore = min( 2 + (((norbz**2+norbz)/2)*norbz**2)/nrank,60000000 )
     allocate(  eribuf( nstore ),  stat = ierr  )


     !     allocate memory for solvers

     if  (  dem_gs  )  then
        allocate( coefflock( xpmax ), stat = ierr )
        call xm_print( 'header', 'direct energy minimization;' )
        call xm_print( 'comment',   &
             'cycle micro orb coeff    relaxn.    rel./tol    ptbn.    C(i);' )
     end  if
     if  (  .not. dem_gs  .or.  nspinc .gt. 1  )  then
        hdim  =  max( xpmax, nspinc )
        allocate( ham( hdim, hdim ),  stat = ierr )
        allocate( ovl( hdim, hdim ),  stat = ierr )
        allocate(   w( hdim, hdim ),  stat = ierr )
        allocate( eig( hdim ),  stat = ierr )
        allocate(  v1( hdim ),  stat = ierr )
        allocate(  v2( hdim ),  stat = ierr )
        call xm_print( 'header', '(full) first-order method;' )
        call xm_print( 'comment',   &
             'cycle  orbital   relaxation(kCal)   (..per orb.)/tol;' )
     end  if

! minimizes energy wrt coefficients
     call minimize_energy( energy, &
          w, eig, v1, v2, coefflock,int_out, dbl_out )

     deallocate( eribuf, stat = ierr )
     if  (  dem_gs  )  then
        deallocate( coefflock, stat = ierr )
     end  if
     if  (  .not. dem_gs  .or.  nspinc .gt. 1  )  then
        deallocate(  v2,  stat = ierr )
        deallocate(  v1,  stat = ierr )
        deallocate( eig,  stat = ierr )
        deallocate(   w,  stat = ierr )
        deallocate( ovl,   stat = ierr )
        deallocate( ham,   stat = ierr )
     end if

  end if    !  single energy or optimization

  deallocate(  ipvt,  stat = ierr )

  deallocate( schwarz,   stat = ierr )
  deallocate( abra_npair,  stat = ierr )
  deallocate( bbra_npair,  stat = ierr )
  deallocate( abra_docc_un,  stat = ierr )
  deallocate( bbra_docc_un,  stat = ierr )
  deallocate( aket_docc_un,  stat = ierr )
  deallocate( bket_docc_un,  stat = ierr )

  deallocate( bket,  stat = ierr )
  deallocate( aket,  stat = ierr )
  deallocate( wdet,  stat = ierr )
  deallocate(  ket,  stat = ierr )
  deallocate(  bra,  stat = ierr )
  deallocate(    kexch,  stat = ierr )
  deallocate(    bexch,  stat = ierr )
  deallocate(    ket_b,  stat = ierr )
  deallocate(    ket_a,  stat = ierr )
  deallocate(    bra_b,  stat = ierr )
  deallocate(    bra_a,  stat = ierr )

  deallocate(    xpnew,  stat = ierr )
  deallocate(  ndf2orb,  stat = ierr )
  deallocate( atom_ndf,  stat = ierr )

  deallocate( coeffl,  stat = ierr )
  deallocate( coeffk,  stat = ierr )
  deallocate( coeffj,  stat = ierr )
  deallocate( coeffi,  stat = ierr )

  deallocate(    dkl,  stat = ierr )

  deallocate(   nxyz,  stat = ierr )
  deallocate(   angn,  stat = ierr )
  deallocate(   ashl,  stat = ierr )
  deallocate(   ashi,  stat = ierr )

  deallocate(    dij,  stat = ierr )


#ifdef SIMINT_INT
     ! deallocates simint-related things
     call valence_finalize_simint
#endif

   end subroutine calculate_vsvb_energy




!>    compute energy for system, without optimization
!!    \param energy [out]:  VSVB energy
subroutine  guess_energy ( energy )
  use         densitywork
  use         integrals
  use tools, only: dp
  use xm
  implicit    none

  real(dp)     energy(1)

  integer     num_non_docc,idocc,num_spatial_orbs,i
  real(dp)     wfnorm(1)

  !     compute a single energy with the current orbitals
  !     setup the wave function: bra and ket

  num_non_docc = 2*npair + nunpd   !  'non-DOCC' orbitals,electrons
  do i = 1, num_non_docc
     bra( i ) = i
     ket( i ) = i
  end do
  i = num_non_docc + 1
  do idocc = num_non_docc + 1, num_non_docc + ndocc
     bra( i   ) = idocc
     ket( i   ) = idocc
     bra( i+1 ) = idocc
     ket( i+1 ) = idocc
     i = i + 2
  end do
  call  wfndet
  num_spatial_orbs  =  num_non_docc + ndocc
  call  schwarz_ints( num_spatial_orbs, num_non_docc )
  store_eri  =  .false.
  call  vsvb_energy( 0 ,num_non_docc,num_spatial_orbs, energy(1), wfnorm(1), .false., .true. )
  call  xm_equalize( energy, 1 )
  call  xm_equalize( wfnorm, 1 )
  energy  =  energy/wfnorm  +  enucrep
end subroutine  guess_energy





subroutine  demgs_opt ( iorb,num_iter,cumulx,energy,etol,coefflock )
  use         densitywork
  use         integrals
  use tools, only: dp
  use xm
  implicit    none

  integer     iorb,num_iter
  real(dp)     cumulx,energy(1),etol

  !     generic print buffers

  integer     int_out(4)    !  lengths not
  real(dp)     dbl_out(4)    !   protected

  integer     i,icoeff,ncoeff,num_non_docc,idocc,eorb,k,num_spatial_orbs,micro_iter,nnrgy
  real(dp)     tokcal, relaxn, perturb,eprev,origc,wfnorm(1)
  parameter ( tokcal = 627.509469_dp )
  logical     updated, coefflock(*)

  !     Direct Energy Minimization for Ground States (DEMGS)
  !     This 'zeroth order' approach is more fine-grained than the
  !     conventional first-order method, and could be more efficient
  !     when the guess is close to the solution


  num_non_docc = 2*npair + nunpd
  do i = 1, num_non_docc
     bra( i ) = i
     ket( i ) = i
  end do
  i = num_non_docc + 1
  do idocc = num_non_docc + 1, num_non_docc + ndocc
     bra( i   ) = idocc
     ket( i   ) = idocc
     bra( i+1 ) = idocc
     ket( i+1 ) = idocc
     i = i + 2
  end do
  eorb  =  iorb

  !     map a DOCC label to the beta-spin 'slot' in the wfn

  if  (  iorb .gt. num_non_docc  )  eorb  =  2*iorb - num_non_docc
  call  wfndet
  num_spatial_orbs  =  num_non_docc + ndocc
  call  schwarz_ints( num_spatial_orbs, num_non_docc )


  nnrgy  =  0
  micro_iter  =  0
  ncoeff  =  map_orbs( iorb + 1 ) - map_orbs( iorb )
  store_eri  =  .true.
  eri_stored  =  .false.
  perturb  =  ptbnmax
  do  while  (  tokcal*(perturb**2)  .ge.  etol  )
     do  i  =  1, ncoeff
        coefflock( i ) = .false.
     end do
     updated  =  .true.
     micro_iter  =  micro_iter  +  1
     do  while  (  updated  )
        updated  =  .false.
        do  i  =  1, ncoeff
           if  (  .not. coefflock( i )  )  then
              icoeff  =  i  +  map_orbs( iorb ) - 1
              eprev  =  energy(1)
              origc  =  coeff( icoeff )

              !     add the perturbation, update wfn, compute its energy

              coeff( icoeff ) = origc + perturb
              call  normal( iorb, iorb )
              do    k = 1, nelec
                 call  ovint( iorb, ket( k ) )
                 wdet( eorb, k ) = sint
                 wdet( k, eorb ) = sint
                 if  (  iorb .gt. num_non_docc  )  then
                    wdet( eorb-1, k ) = sint
                    wdet( k, eorb-1 ) = sint
                 end  if
              end  do
              call  vsvb_energy( iorb,num_non_docc,num_spatial_orbs, energy(1), wfnorm(1), .false.,.false. )
              call  xm_equalize( energy, 1 )
              call  xm_equalize( wfnorm, 1 )
              energy(1)  =  energy(1)/wfnorm(1)  +  enucrep
              nnrgy = nnrgy + 1
              eri_stored  =  .true.

              if  (  energy(1) .lt. eprev  )  then
                 updated  =  .true.
                 relaxn  =  ( energy(1) - eprev )*tokcal
                 cumulx  =  cumulx  +  relaxn
                 int_out( 1 ) = num_iter
                 int_out( 2 ) = micro_iter
                 int_out( 3 ) = iorb
                 int_out( 4 ) = i
                 dbl_out( 1 ) = cumulx
                 dbl_out( 2 ) = relaxn/etol
                 dbl_out( 3 ) = perturb
                 dbl_out( 4 ) = coeff( icoeff )
                 call  xm_print( 'demgs', ' ', int_out, dbl_out )
                 call xm_output( 'save', energy(1),etol )
                 eprev  =  energy(1)
              else

                 !     subtract the perturbation

                 coeff( icoeff ) = origc - perturb
                 call  normal( iorb, iorb )
                 do    k = 1, nelec
                    call  ovint( iorb, ket( k ) )
                    wdet( eorb, k ) = sint
                    wdet( k, eorb ) = sint
                    if  (  iorb .gt. num_non_docc  )  then
                       wdet( eorb-1, k ) = sint
                       wdet( k, eorb-1 ) = sint
                    end  if
                 end  do
                 call  vsvb_energy( iorb,num_non_docc,num_spatial_orbs, energy(1), wfnorm(1), .false.,.false. )
                 call  xm_equalize( energy, 1 )
                 call  xm_equalize( wfnorm, 1 )
                 energy(1)  =  energy(1)/wfnorm(1)  +  enucrep
                 nnrgy = nnrgy + 1

                 if  (  energy(1) .lt. eprev  )  then
                    updated  =  .true.
                    relaxn  =  ( energy(1) - eprev )*tokcal
                    cumulx  =  cumulx  +  relaxn
                    int_out( 1 ) = num_iter
                    int_out( 2 ) = micro_iter
                    int_out( 3 ) = iorb
                    int_out( 4 ) = i
                    dbl_out( 1 ) = cumulx
                    dbl_out( 2 ) = relaxn/etol
                    dbl_out( 3 ) = -perturb
                    dbl_out( 4 ) = coeff( icoeff )
                    call  xm_print( 'demgs', ' ', int_out, dbl_out )
                    call xm_output( 'save', energy(1),etol )
                    eprev  =  energy(1)
                 else

                    !     restore the original weight and lock it for this perturbation

                    coeff( icoeff ) = origc
                    coefflock( i )  =  .true.
                    call  normal( iorb, iorb )
                    do    k = 1, nelec
                       call  ovint( iorb, ket( k ) )
                       wdet( eorb, k ) = sint
                       wdet( k, eorb ) = sint
                       if  (  iorb .gt. num_non_docc  )  then
                          wdet( eorb-1, k ) = sint
                          wdet( k, eorb-1 ) = sint
                       end  if
                    end  do
                    energy  =  eprev

                 end  if    !  subtract perturbation
              end  if    !  add perturbation
           end  if    !  weight lock
        end  do    !  orbital weights
     end  do    !  updates
     perturb  =  perturb * feather
  end  do    !  perturbations

  int_out(1) = nnrgy
  call xm_print( 'dimension', 'number of energy calculations;',  &
       int_out, dbl_out )
  store_eri  =  .false.
end subroutine  demgs_opt





subroutine  first_order_opt ( iorb, w,eig,v1,v2, energy )
  use         densitywork
  use         integrals
#ifdef SIMINT_INT
  use SimintFortran
  use valence_simint, only: shell_map
#endif
  use tools, only: dp
  use xm
  implicit    none

  integer     iorb
  real(dp)     w( hdim, *), eig(*),v1(*),v2(*), energy

  integer     num_non_docc,idocc,norbas,num_spatial_orbs,eorb
  integer     iorbas,iat,it,ish,iao,i,j,k,n,mnshi,mxshi
  integer     new_orb,new_atm,new_typ, ib,jb,idf
  integer     dovecs,ifail
  parameter ( dovecs = 1 )    !  =0 for eigenvalues only
  real(dp)     zero,           one, wfnorm
  parameter ( zero = 0.0d+00, one = 1.0d+00 )

  num_non_docc = 2*npair + nunpd   !  'non-DOCC' orbitals,electrons
  num_spatial_orbs  =  num_non_docc + ndocc

  !     construct the wave function orbital lists so as to
  !     avoid optimizing DOCC orbitals twice in succession
  !     convert the input orbital address (iorb) to an electron
  !     address (eorb) in the wfn

  eorb = 0
  if (  iorb .le. num_non_docc  )  then

     eorb  =  iorb
     do i = 1, num_non_docc
        bra( i ) = i
        ket( i ) = i
     end do

     i = num_non_docc + 1
     do idocc = num_non_docc + 1, num_non_docc + ndocc
        bra( i   ) = idocc
        ket( i   ) = idocc
        bra( i+1 ) = idocc
        ket( i+1 ) = idocc
        i = i + 2
     end do

  else

     !     to optimize a DOCC, append it to the SC list in the wfn, skip
     !     over it in the DOCC list, increment num_non_docc, decrement ndocc

     do i = 1, num_non_docc
        bra( i ) = i
        ket( i ) = i
     end do
     bra( num_non_docc + 1 ) = iorb
     ket( num_non_docc + 1 ) = iorb
     bra( num_non_docc + 2 ) = iorb
     ket( num_non_docc + 2 ) = iorb

     i = num_non_docc + 3
     do idocc = num_non_docc + 1, num_non_docc + ndocc
        if  (  idocc .ne. iorb  ) then
           bra( i   ) = idocc
           ket( i   ) = idocc
           bra( i+1 ) = idocc
           ket( i+1 ) = idocc
           i = i + 2
        end if
     end do

     !     here, the alpha-spin orbital ('1') is substituted

     eorb  =  2*npair + nunpd  +  1
     num_non_docc  =  num_non_docc  +  2
     num_spatial_orbs  =  num_spatial_orbs + 1
  end if


  if  (  iorb  .gt.  0  )   then
     new_orb  =  norbs
     new_atm  =  natom
     new_typ  =  natom_t

     !     to map from the variational subspaces to the unique
     !     basis set data, the least overhead is to make a 'dummy'
     !     NDF for each regular AO, including 1-shell atom type
     !     for simplicity, input NDF's just leave a null 'gap'
     !     in this orbital list

     do  ib = map_orbs( iorb ), map_orbs( iorb + 1 ) - 1
        new_orb  =  new_orb  +  1

        iorbas = 0
        do  iat = 1,  orbas_atnum(      iorb )
           it =  atom_t( orbas_atset( iat, iorb ) )
           mnshi =         map_atom2shell( it )
           mxshi = mnshi + num_shell_atom( it ) - 1
           do  ish = mnshi, mxshi
              do  iao = 1,  shell_size( ang_mom( ish ) )
                 iorbas = iorbas + 1
                 if (  iorbas .eq. xpset( ib )  ) then

                    new_atm  =  new_atm  +  1
                    new_typ  =  new_typ  +  1

                    orbas_atnum(    new_orb ) = 1
                    orbas_atset( 1, new_orb ) = new_atm
                    xpset( map_orbs( new_orb ) ) = iao
                    coeff( map_orbs( new_orb ) ) = one
                    coords( 1, new_atm ) = coords( 1, orbas_atset( iat, iorb ) )
                    coords( 2, new_atm ) = coords( 2, orbas_atset( iat, iorb ) )
                    coords( 3, new_atm ) = coords( 3, orbas_atset( iat, iorb ) )
                    atom_t(    new_atm ) = new_typ
                    num_shell_atom( new_typ ) = 1
                    map_atom2shell( new_typ ) = ish


#ifdef SIMINT_INT
                    ! make a new simint_shell here
                    call simint_initialize_shell( shell_map(1,new_atm) )

                    call simint_create_shell( map_shell2prim( ish + 1) - map_shell2prim( ish ), ang_mom( ish ), &
                         coords( 1, new_atm),coords( 2, new_atm),coords( 3, new_atm), &
                         exponent(map_shell2prim( ish )), con_coeff(map_shell2prim( ish )), &
                         shell_map(1,new_atm) )
#endif

                 end  if
              end  do    !  shell function
           end  do    !  shell
        end  do    !  atom

        map_orbs( new_orb + 1 ) = map_orbs( new_orb ) + 1
     end  do    !  orbital AO
  end  if    !  orbital

  call  wfndet
  call  schwarz_ints( num_spatial_orbs, num_non_docc )

  !     loop over the hamiltonian and overlap matrices

  store_eri  =  .true.
  eri_stored  =  .false.
  norbas  = map_orbs( iorb + 1 ) - map_orbs( iorb )
  do  ib  =  1,  norbas
     idf = xpset( map_orbs( iorb ) + ib - 1 )
     if  (  idf .lt. 1  )  then
        bra( eorb )  =  norbs + idf
     else
        bra( eorb )  =  norbs + ib
     end  if
     do    k = 1, nelec
        call  ovint( bra( eorb ), ket( k ) )
        wdet( eorb, k ) = sint
     end  do

     do  jb  =  1,  ib
        idf = xpset( map_orbs( iorb ) + jb - 1 )
        if  (  idf .lt. 1  )  then
           ket( eorb )  =  norbs + idf
        else
           ket( eorb )  =  norbs + jb
        end  if
        do    k = 1, nelec
           call  ovint( bra( k ), ket( eorb ) )
           wdet( k, eorb ) = sint
        end  do

        call  vsvb_energy( iorb,num_non_docc,num_spatial_orbs, energy, wfnorm, .false., .false. )
        eri_stored  =  .true.

        ham( ib, jb ) = energy
        ovl( ib, jb ) = wfnorm
     end  do   !  jb
  end  do   !  ib


  !     spin average to optimize DOCC with unpaired electrons

  if  (  nunpd .gt. 0  .and.  iorb .gt. 2*npair + nunpd  )  then

     bra( eorb ) = iorb
     ket( eorb ) = iorb
     call  wfndet

     eorb  =  2*npair + nunpd  +  2        ! beta-spin position
     store_eri  =  .true.
     eri_stored  =  .false.
     do  ib  =  1,  norbas
        idf = xpset( map_orbs( iorb ) + ib - 1 )
        if  (  idf .lt. 1  )  then
           bra( eorb )  =  norbs + idf
        else
           bra( eorb )  =  norbs + ib
        end  if
        do    k = 1, nelec
           call  ovint( bra( eorb ), ket( k ) )
           wdet( eorb, k ) = sint
        end  do

        do  jb  =  1,  ib
           idf = xpset( map_orbs( iorb ) + jb - 1 )
           if  (  idf .lt. 1  )  then
              ket( eorb )  =  norbs + idf
           else
              ket( eorb )  =  norbs + jb
           end  if
           do    k = 1, nelec
              call  ovint( bra( k ), ket( eorb ) )
              wdet( k, eorb ) = sint
           end  do

           call  vsvb_energy( iorb,num_non_docc,num_spatial_orbs, energy,wfnorm, .true., .false. )
           eri_stored  =  .true.

           ham( ib, jb ) = ham( ib, jb ) + energy
           ovl( ib, jb ) = ovl( ib, jb ) + wfnorm
        end  do   !  jb
     end  do   !  ib
  end  if   !  spin averaging


  !     symmetrize ham and ovl

  do i = 1, norbas
     do j = 1, i - 1
        ham( j, i ) = ham( i, j )
        ovl( j, i ) = ovl( i, j )
     end do
  end do

  !     gather hamiltonian and overlap matrices across all ranks

  call xm_equalize( ham, hdim*norbas )
  call xm_equalize( ovl, hdim*norbas )
#ifdef PRINT_MATRIX
  call write_matrix(ham,norbas,norbas,'ham')
  call write_matrix(ovl,norbas,norbas,'ovl')
#endif
  !     solve the generalized eigenvalue problem

  call  rsg( hdim, norbas, ham, ovl, eig,   &
       dovecs, w, v1, v2, ifail )

 if ( ifail .ne. 0 ) then
    call xm_print('error','EISPACK rsg error code; ', [ifail])
    if (ifail < hdim) then
      call xm_print('comment','Error code is the index of the unconverged eigenval;')
    else if (ifail == 7*hdim+1) then
      call xm_print('comment','Overlap matrix is not positive definite;')
      call xm_print('comment','Try tightening integral screening threshold;')
      call write_matrix(ovl,norbas,norbas,'overlap')
      call xm_print('comment','Written to disk: overlap.mm;')
    else
      call xm_print('comment','Unknown error code;')
    end if
    call xm_abort('solver failed;')
  end if 
  !     select desired root of the eigenproblem and update the wfn
  !     (rsg sorts the eigenvalues)

  energy = eig( 1 )
  do  i  =  1,  nxorb
     if ( iorb .eq. xorb( i ) ) energy = eig( root( i ) + 1 )
  end do
  energy = energy + enucrep

  do  i  =  1, norbas
     coeff( i + map_orbs( iorb ) - 1 ) = w( i, 1 )
     do  j  =  1,  nxorb
        if ( iorb .eq. xorb( j ) )  &
             coeff( i + map_orbs( iorb ) - 1 ) = w( i, root( j ) + 1 )
     end do
  end do

  !     normalize new orbital

  call  normal( iorb, iorb )
  store_eri  =  .false.

#ifdef SIMINT_INT
! free simint memory
  if  (  iorb  .gt.  0  )   then
     new_orb  =  norbs
     new_atm  =  natom

     do  ib = map_orbs( iorb ), map_orbs( iorb + 1 ) - 1
        new_orb  =  new_orb  +  1

        iorbas = 0
        do  iat = 1,  orbas_atnum(      iorb )
           it =  atom_t( orbas_atset( iat, iorb ) )
           mnshi =         map_atom2shell( it )
           mxshi = mnshi + num_shell_atom( it ) - 1
           do  ish = mnshi, mxshi
              do  iao = 1,  shell_size( ang_mom( ish ) )
                 iorbas = iorbas + 1
                 if (  iorbas .eq. xpset( ib )  ) then

                    new_atm  =  new_atm  +  1
                    ! make a new simint_shell here
                    call simint_free_shell( shell_map(1,new_atm) )


                 end  if
              end  do    !  shell function
           end  do    !  shell
        end  do    !  atom

!        map_orbs( new_orb + 1 ) = map_orbs( new_orb ) + 1
     end  do    !  orbital AO
  end  if    !  orbital
#endif

end subroutine  first_order_opt





subroutine  spin_opt ( w,eig,v1,v2, energy )
  use         densitywork
  use         integrals
  use tools, only: dp
  use xm
  implicit    none

  real(dp)     w( hdim, *), eig(*),v1(*),v2(*), energy

  integer     i,j,num_non_docc,idocc,num_spatial_orbs
  integer     dovecs,ifail
  parameter ( dovecs = 1 )
  real(dp)     wfnorm


  !     optimize the spin coupling weights

  num_non_docc = 2*npair + nunpd
  do i = 1, num_non_docc
     bra( i ) = i
     ket( i ) = i
  end do
  i = num_non_docc + 1
  do idocc = num_non_docc + 1, num_non_docc + ndocc
     bra( i   ) = idocc
     ket( i   ) = idocc
     bra( i+1 ) = idocc
     ket( i+1 ) = idocc
     i = i + 2
  end do
  call  wfndet
  num_spatial_orbs  =  num_non_docc + ndocc
  call  schwarz_ints( num_spatial_orbs, num_non_docc )

  do i = 1, nspinc
     do j = 1, nspinc
        ovl( i, j ) = 0.0d+00  
        ham( i, j ) = 0.0d+00
     end do
  end do

  store_eri  =  .false.
  spinopt  =  .true.
  call  vsvb_energy( 0,num_non_docc,num_spatial_orbs, energy, wfnorm, .false.,.false. )
  spinopt  =  .false.

  call xm_equalize( ham, nspinc*hdim )
  call xm_equalize( ovl, nspinc*hdim )

  do i = 1, nspinc
     do j = 1, i
        ovl( i, j ) = ovl( i, j )/dble( nelec )
     end do
  end do

  do i = 1, nspinc
     do j = 1, i
        ham( j, i ) = ham( i, j )
        ovl( j, i ) = ovl( i, j )
     end do
  end do

  call  rsg( hdim, nspinc, ham, ovl, eig,   &
       dovecs, w, v1, v2, ifail )
 if ( ifail .ne. 0 ) then
    call xm_print('error','EISPACK rsg error code; ', [ifail])
    if (ifail < hdim) then
      call xm_print('comment','Error code is the index of the unconverged eigenval;')
    else if (ifail == 7*hdim+1) then
      call xm_print('comment','Overlap matrix is not positive definite;')
      call xm_print('comment','Try tightening integral screening threshold;')
      call write_matrix(ovl,nspinc,nspinc,'overlap')
      call xm_print('comment','Written to disk: overlap.mm;')
    else
      call xm_print('comment','Unknown error code;')
    end if
    call xm_abort('solver failed;')
  end if 

  !     select lowest root 
  !     note, SC weights are not normalized

  energy  =  eig( 1 )  +  enucrep
  do i = 1, nspinc
     coeff_sc( i ) = w( i, 1)
  end do
end subroutine  spin_opt




!>    compute energy integral and normalization integral for system
!!
!!  For nelec electrons, and a single determinant wavefunction (no spin
!!  couplings), the energy integral is
!!  \f$ E = \sum_{ij}^{nelec} [ d^1_{ij} w_{ij} h_{s,ij} ]
!!    + \sum_i^{nelec}\sum_{j<i}\sum_k^{nelec}\sum_{l<k}
!!    [ d^2_{ikjl} w_{ikjl} (< i(1) j(2) | k(1) l(2) >_s-< i(1) j(2) | l(1) k(2) >_s) ]  \f$
!!
!!  For N_p spin coupled pairs, and N_{sc} spin couplings, there are
!!  M_{sc} = N_{sc}^2 2^{2N_p} terms like the one above, which differ in spin
!!  functions:
!!
!!  \f$ E = \sum_{ij}^{nelec} [ (\sum_q^{M_{sc}} d^1_{q,ij} w_{q,ij} ) h_{s,ij} ]
!!    + \sum_i^{nelec}\sum_{j<i}\sum_k^{nelec}\sum_{l<k}
!!        [ (\sum_q^{M_{sc}} d^2_{q,ikjl}
!!          (  w_{q,ikjl} < i(1) j(2) | k(1) l(2) >_s
!!             - w_{q,iljk} < i(1) j(2) | l(1) k(2) >_s ) ]  \f$
!!
!!  The first term in the sum is the one electron term.
!!
!!  \f$ h_{s,ij} = < i(1) | h_s | j(1) >  \f$ where h is the standard one-electron
!!  electron kinetic and nuclei-electron attraction operator.
!!  the s subscript is meant to denote that this integral does not include
!!  spin integration (this is in the "w" variable)
!!  i,j are one-electron spin orbitals.
!!
!!  \f$ d^1_{q,ij} \f$ is the first-order cofactor of the matrix of overlap integrals
!!  between the spin orbitals. That is, it is the determinant of the overlap
!!  matrix with row i and column j removed, multiplied by (-1)^{i+j}
!!  q denotes the term in the spin coupling/pairs expansion
!!
!!  \f$ w_{q,ij} \f$ is the spin function integration, where the spin functions
!!  are alpha or beta, whichever are associated with spin orbital i(1) and j(1)
!!  q denotes the term in the spin coupling/pairs expansion
!!
!!  The second term in the sum is the two electron term.
!!
!!  \f$ < i(1) j(2) | k(1) l(2) >_s \f$ is the electron-electron repulsion integral.
!!  i,j,k,l are one-electron spin orbitals.
!!  the s subscript is meant to denote that this integral does not include
!!  spin integration (this is in the "w" variable)
!!
!!  \f$ d^2_{q,ikjl} \f$ is the second-order cofactor of the matrix of overlap integrals
!!  between the spin orbitals. That is, it is the determinant of the overlap
!!  matrix with row i, column k, row j, and column l removed,
!!  multiplied by (-1)^{i+j+k+l}
!!  q denotes the term in the spin coupling/pairs expansion
!!
!!  \f$ w_{q,ijkl} \f$ is the spin function integration, where the spin functions
!!  are alpha or beta, whichever are associated with spin orbitals
!!  i(1),k(1),j(2),l(2). q denotes the term in the spin coupling/pairs expansion
!!
!!  The normalization integral is the determinant of the matrix of overlap
!!  integrals between spin orbitals.
!!  N = Determinant of spin orbital overlap matrix
!! \f$  = \sum_{j}^{nelec} [(\sum_q^{M_{sc}} d^1_{q,1j} w_{q,1j} ) < 1(1) | j(1) >_s ]
!! \f$
!!  \param iorb [in]: perturbed orbital, if called from an optimization
!!                      routine
!!  \param num_non_docc [in]: number of non-docc orbitals (2*npair+nunpd)
!!  \param num_spatial_orbs [in]: number of "spatial" orbitals, not spin orbitals
!!                        (2*npair+nunpd+ndocc)
!!  \param energy [out]: VSVB energy
!!  \param wfnorm [out]: normalization term for the wavefunction
!!  \param spinav [in]: logical flag controlling whether or not to
!!                        spin-average
!!  \param ijkl_symmetry_is_enabled [in]: logical flag controlling whether or not to
!!                        to use kl < ij symmetry. this should be false
!!                        for optimizations where bra != ket
subroutine  vsvb_energy ( iorb,num_non_docc,num_spatial_orbs, energy, wfnorm, spinav,&
     ijkl_symmetry_is_enabled )
  use         densitywork
  use         integrals
  use tools, only: dp
  use xm
  implicit    none

  integer     iorb,num_non_docc,num_spatial_orbs
  real(dp)     energy,wfnorm

  integer     i,j,k,l, io,jo,ko,lo, ib,jb,kb,lb, is,js,ks,ls
  integer     ie,je,ke,le, iso,jso,kso,lso
  real(dp)     zero, d1,d2,erep_int, exchanged_erep_int,&
       erep_sum, exchanged_erep_sum, dummy
  parameter ( zero = 0.0d+00 )

  real(dp)  erep_density,exchanged_erep_density

  integer     ijorb,klorb
  integer     nproc,myrank,master, indx,ericount,jorb
  logical     nonsub,signif, spinav

  logical erep_is_signif, exchange_is_signif, calculate_integrals, &
       ijkl_symmetry_is_enabled

  logical i_is_docc, j_is_docc, k_is_docc, l_is_docc, ik_are_docc_and_spin_is_zero,&
       ij_are_docc_and_spin_is_zero, jl_are_docc_and_spin_is_zero,&
       jk_are_docc_and_spin_is_zero,il_are_docc_and_spin_is_zero,&
       compute_term, compute_erep_term,compute_exchanged_erep_term


  !     task lists could be long 
  integer(8)  ntasks,task,nproc8,ij8,ntasks_ijo,ntasks_klo,&
       mytasks, loctask

  nproc8  =  nrank                           !  compiler
  task  =  0

  !     compute 1e- energy and wfn norm
  !     energy = \sum_{ij}^{nelec} [ (\sum_q^{M_{sc}} d^1_{q,ij} w_{q,ij} )
  !                                   h_{s,ij} ]
  !     wfnorm = \sum_{j}^{nelec} [ (\sum_q^{M_{sc}} d^1_{q,1j} w_{q,1j} )
  !                                < 1(1) | j(1) >_s ]
  !         or, as it is done here:
  !     wfnorm = \sum_{ij}^{nelec} [ (\sum_q^{M_{sc}} d^1_{q,ij} w_{q,ij} )
  !                                < i(1) | j(1) >_s ] / nelec
  !
  ! set up an array to store arrays for the unpaired and DOCC e-
  ! used in computing the density
  call set_up_unpaired_docc

  energy = zero
  wfnorm = zero

  ! i_is_docc (j_is_docc) is a flag that is true if spin orbital i (j) is part of a docc orbital.
  ! if both i and j are part of docc orbitals, that means that we know that we have
  ! alpha and beta spins and can do the spin integration easily before the integral
  ! and cofactor calculation.
  ! the spin integration can be done easily because we know that for docc orbitals
  ! alpha spins will always have the same parity of position and beta spins will
  ! also match parity. So all we need to do is check that both are the same parity.
  do  i  =  1,  nelec
     i_is_docc = .false.; if ( i .gt. num_non_docc ) i_is_docc = .true.
     do  j  =  1,  nelec
        j_is_docc = .false.; if ( j .gt. num_non_docc ) j_is_docc = .true.

        ij_are_docc_and_spin_is_zero = i_is_docc .and. j_is_docc &
             .and. mod( i, 2 ) .ne. mod( j, 2 )

        if( ij_are_docc_and_spin_is_zero ) then
           ! the spin integration is 0, so don't compute the term
           compute_term = .false.
        else
           compute_term = .true.
        endif

        if( compute_term ) then
           task   =  task  +  1
           if  (  mod( task, nproc8 ) .eq. irank  )  then

              call   ovint( bra(i), ket(j) )
              call   int1e( bra(i), ket(j) )
              dme_b(1) = i ;  dme_k(1) = j
              call density( 1 , d1, dummy,dummy,dummy, .true.,.false. )
              wfnorm  =  wfnorm  +  sint*d1
              energy  =  energy  +  hint*d1

           end  if   !  round-robin
        end  if   !  spin allowed
     end  do   !  j electron
  end  do   !  i electron


  !     scale wfn norm to unit operator
  !     this is divided by nelec since it only needs to be done for one "i"
  wfnorm  =  wfnorm/dble( nelec )

  !     compute 2e- energy
  !
  !
  !  2e- energy = \sum_i^{nelec}\sum_{j<i}\sum_k^{nelec}\sum_{l<k}
  !                 [ ( \sum_q^{M_{sc}} d^2_{q,ikjl}
  !                 ( w_{q,ikjl} < i(1) j(2) | k(1) l(2) >_s
  !                    - w_{q,iljk} < i(1) j(2) | l(1) k(2) >_s ) ]
  !
  !  Practically, to take advantage of the fact that spin integration removes
  !  a lot of terms, the code evaluates an expansion:
  !
  ! = \sum_{io}^{orbitals}\sum_{jo<=io}\sum_{ko}^{orbitals}\sum_{lo<=ko}
  !   \sum_i^{spins in io}\sum_j^{spins in jo<io}\sum_k^{spins in ko}
  !   \sum_l^{spins in lo<ko} [ ( \sum_q^{M_{sc}} d^2_{q,io+i,ko+k,jo+j,lo+l}
  !    ( w_{q,ikjl}  < io(1) jo(2) | ko(1) lo(2) >
  !    - w_{q,iljk} < io(1) jo(2) | lo(1) ko(2) > ) ) ]
  !
  !   where "orbitals" are "spatial" orbitals that have no spin coordinates
  !
  ! = \sum_{io}^{orbitals}\sum_{jo<=io}\sum_{ko}^{orbitals}\sum_{lo<=ko}
  !   \sum_i^{spins in io}\sum_j^{spins in jo<io}\sum_k^{spins in ko}
  !   \sum_l^{spins in lo<ko}
  !    [ erep_density(io,i,ko,k,jo,j,lo,l) * erep_int(io,jo,ko,lo)
  !   - exchanged_erep_density(io,i,ko,k,jo,j,lo,l) * exchanged_erep_int(io,jo,lo,ko) ]
  !
  !    where erep_density(io,i,ko,k,jo,j,lo,l) =
  !     \sum_q^{M_{sc}} d^2_{q,io+i,ko+k,jo+j,lo+l} w_{q,ikjl}
  !
  !    exchanged_erep_density(io,i,ko,k,jo,j,lo,l) =
  !     \sum_q^{M_{sc}} d^2_{q,io+i,ko+k,jo+j,lo+l} w_{q,iljk}
  !
  !    erep_int(io,jo,ko,lo) = < io(1) jo(2) | ko(1) lo(2) >
  !
  !    exchanged_erep_int(io,jo,lo,ko) = < io(1) jo(2) | lo(1) ko(2) >
  !
  !    where the spin sums are inside the terms.
  !
  !     task distribution is in contiguous 'blocks'
  !     but still 'round-robin' for even load-balancing

  !     note that this is a loop over "spatial" (that is, not including spin)
  !     two-electron docc orbitals and "spatial" one-electron unpaired and spin
  !     coupled orbitals, not spin orbitals. For each "spatial" orbital, the spin is
  !     considered in the density() call. (if the orbitals are docc,
  !     the spins and spin integration is considered before the density() call)
  ntasks_ijo  =  (num_spatial_orbs**2 + num_spatial_orbs)/2
  ntasks_klo  =  (num_spatial_orbs**2 + num_spatial_orbs)/2
  ij8     =  ntasks_ijo
  ntasks  =  ntasks_ijo * ntasks_klo
  if( ijkl_symmetry_is_enabled )  ntasks = (ntasks_klo**2 + ntasks_klo)/2
  mytasks =  ntasks/nproc8
  if ( irank .lt. mod( ntasks, nproc8 ) ) mytasks = mytasks + 1

  ericount  =  0
  task  =  irank
  do loctask  =  1, mytasks

     !     back-out orbital labels from the task index: kl; then ij

     if( ijkl_symmetry_is_enabled ) then
        call xm_dtriang8( 1+task, ijorb, klorb )
     else
        klorb  =  1  +       task/ij8
        ijorb  =  1  +  mod( task,ij8 )
     endif

     call xm_dtriang( ijorb, io, jo )
     call xm_dtriang( klorb, ko, lo )

     ! if io == jo (or ko == lo) and they are not docc
     ! (they are spin coupled or unpaired), don't calculate the integral
     ! and density, since two different electrons cannot be in the same
     ! spin orbital (docc orbitals have alpha and beta spin so this is
     ! possible for them)
     if( (io .eq. jo .and. io .le. num_non_docc) &
          .or. (ko .eq. lo .and. ko .le. num_non_docc) ) then
        task  =  task  +  nrank
        cycle
     endif

     !     screen the electron repulsion integrals  (ERI)
     erep_is_signif = schwarz(indx( io, ko ))*schwarz(indx( jo, lo )) .gt. itol
     exchange_is_signif = schwarz(indx( io, lo ))*schwarz(indx( jo, ko )) .gt. itol
     exchanged_erep_int = zero
     erep_int = zero
     gint  =  zero
     if  (  erep_is_signif .or. exchange_is_signif  )  then       !     compute the ERI

        !     in order to track the subject orbital for storage purposes
        !     and account for the case where the subject DOCC appends the
        !     non-DOCC set (for efficiency), adjust 'iorb' here

        jorb  =  iorb
        if  (  iorb .gt. 2*npair+nunpd  .and. .not.dem_gs  )  then
           jorb  =  2*npair + nunpd +  1       !  alpha spin
           if  (  spinav  )   &
                jorb  =  2*npair + nunpd +  2       !  beta  spin
        end  if
        nonsub  =                 io .ne. jorb
        nonsub  =  nonsub  .and.  jo .ne. jorb
        nonsub  =  nonsub  .and.  ko .ne. jorb
        nonsub  =  nonsub  .and.  lo .ne. jorb

        ! check if the integral is one of the schwarz integrals, which were
        ! previously calculated
        if( nonsub .and. io .eq. jo .and. ko .eq. lo ) then
           ! ko == lo, so erep_int and exchanged_erep_int are the same
           erep_int = schwarz( indx( io, ko ) )*schwarz( indx( io, ko ) )
           exchanged_erep_int = erep_int
        else

           !     map the orbital labels to an electron 'slot' in the wfn
           ie  = io  ;  if  ( io .gt. num_non_docc )  ie  =  2*io - num_non_docc -1
           je  = jo  ;  if  ( jo .gt. num_non_docc )  je  =  2*jo - num_non_docc -1
           ke  = ko  ;  if  ( ko .gt. num_non_docc )  ke  =  2*ko - num_non_docc -1
           le  = lo  ;  if  ( lo .gt. num_non_docc )  le  =  2*lo - num_non_docc -1

           calculate_integrals = .false.

           if  (  store_eri  )  then
              if  (  nonsub  )  then        !  store significant integrals NOT over
                 !     the subject orbital, up to 'nstore' storage
                 !     note that, for simplicity, storage is based on the schwarz inequality
                 !     not the actual value of the ERI
                 ericount  =  ericount  +  2
                 if  (  ericount  .le.  nstore  )  then   ! ERI can be or has been stored
                    if  (  eri_stored  )  then    !  retrieve previously stored integral
                       calculate_integrals = .false.
                       exchanged_erep_int = eribuf( ericount )
                       erep_int = eribuf( ericount - 1 )
                    else                          !  compute and store
                       calculate_integrals = .true.
                       ! stored after computing below
                    end  if
                 else          !   insufficient storage, compute it anyway
                    calculate_integrals = .true.
                 end  if       !   enough storage?
              else          !   ERI over the subject orbital
                 calculate_integrals = .true.
              end  if       !   subject/non-subject ERI
           else          !   ERI storage is deactivated
              calculate_integrals = .true.
           end  if       !   using storage? yes/no

           ! calculate integrals if necessary
           if( calculate_integrals ) then
              if( erep_is_signif ) then
                 ! calculate < io(1) jo(2) | ko(1) lo(2) >
                 call  int2e( bra(ie), ket(ke), bra(je), ket(le) )
                 erep_int = gint
              endif
              if( exchange_is_signif ) then
                 ! calculate < io(1) jo(2) | lo(1) ko(2) >
                 if( erep_is_signif .and. (ket(le) .eq. ket(ke)) ) then
                    exchanged_erep_int=gint
                 else
                    call  int2e( bra(ie), ket(le), bra(je), ket(ke) )
                    exchanged_erep_int=gint
                 endif
              endif
              ! store the integral if needed
              if( store_eri .and. nonsub .and. (ericount .le. nstore) .and.&
                   (.not. eri_stored) ) then
                 eribuf( ericount ) = exchanged_erep_int
                 eribuf( ericount - 1 ) = erep_int
              endif
           endif

        endif   ! ending check if the integral was a schwarz integral

     end  if       !   ERI is likely to be significant


     !     assess significance of actual (computed or fetched) ERI
     !     choice here is to screen the ERI for the density, could be vice-versa
     !     but the density becomes more expensive in the large limit even though
     !     sometimes the ERI is more expensive

     erep_is_signif = abs( erep_int ) .gt. itol
     exchange_is_signif = abs( exchanged_erep_int ) .gt. itol
     erep_sum = zero
     exchanged_erep_sum = zero
     erep_density = zero
     exchanged_erep_density = zero
     if  (  erep_is_signif .or. exchange_is_signif  )  then           !     compute the density

        !  loop over spins associated with the "spatial" orbitals io,jo,ko,lo and
        !  sum up the density cofactors
        !  \sum_i^{spins in io}\sum_j^{spins in jo<io}\sum_k^{spins in ko}
        !  \sum_l^{spins in lo<ko}
        !   [ erep_density(io,i,ko,k,jo,j,lo,l) * erep_int(io,jo,ko,lo)
        !   - exchanged_erep_density(io,i,ko,k,jo,j,lo,l) * exchanged_erep_int(io,jo,lo,ko) ]

        !     setup for orbital labels for density

        ! if i_is_docc = false, it means that io is not a docc orbital.
        ! if i_is_docc = true, it means that io is a docc orbital.
        ! (if io and jo (or ko and lo) are docc orbitals, the parity of the associated
        ! spin orbitals can easily tell if the spin integration is 0 or not, so it's
        ! convenient to do it here.)
        !
        ! iso is the number of spin orbitals associated with "spatial" orbital io.
        ! for docc orbitals, this is 2 (alpha and beta). for orbitals which hold
        ! unpaired electrons or spin coupled electrons, this is 1 (unpaired are always
        ! alpha, spin coupled can be alpha or beta.
        i_is_docc = .false. ; iso   =  1
        j_is_docc = .false. ; jso   =  1
        k_is_docc = .false. ; kso   =  1
        l_is_docc = .false. ; lso   =  1

        if ( io .gt. num_non_docc ) then
           i_is_docc = .true. ; iso = 2
        end if
        if ( jo .gt. num_non_docc ) then 
           j_is_docc = .true. ; jso = 2
        end if
        if ( ko .gt. num_non_docc ) then
           k_is_docc = .true. ; kso = 2
        end if
        if ( lo .gt. num_non_docc ) then
           l_is_docc = .true. ; lso = 2
        end if

        !     for the given "spatial" orbitals io,jo,ko,lo,
        !     loop over spin orbitals to get the density
        do is =  1, iso
           i = io  ;  if ( io .gt. num_non_docc )  i = 2*io - num_non_docc + is -2
           do js =  1, jso
              j = jo  ;  if ( jo .gt. num_non_docc )  j = 2*jo - num_non_docc + js -2

              ! when io > jo, "is" is always > "js", so this is always true.
              ! if io == jo, the following condition this avoids spin orbitals j < i,
              ! since we're strictly doing spin orbitals j < i
              if  (  j  .lt.  i  )  then
                 do ks =  1, kso
                    k = ko  ;  if ( ko .gt. num_non_docc )  k = 2*ko - num_non_docc + ks -2
                    do ls =  1, lso
                       l = lo  ;  if ( lo .gt. num_non_docc )  l = 2*lo - num_non_docc + ls -2

                       if  (  l  .lt.  k  )  then

                          !     integrate over spin if both are docc

                          ! docc spin integration for erep
                          ik_are_docc_and_spin_is_zero = &
                               i_is_docc .and. k_is_docc .and. mod( i, 2 ) .ne. mod( k, 2 )
                          jl_are_docc_and_spin_is_zero = &
                               j_is_docc .and. l_is_docc .and. mod( j, 2 ) .ne. mod( l, 2 )

                          if( ik_are_docc_and_spin_is_zero .or. jl_are_docc_and_spin_is_zero ) then
                             ! the spin integration is 0, so don't compute the term
                             compute_erep_term = .false.
                          else
                             compute_erep_term = .true.
                          endif

                          ! docc spin integration for the exchanged erep
                          il_are_docc_and_spin_is_zero = &
                               i_is_docc .and. l_is_docc .and. mod( i, 2 ) .ne. mod( l, 2 )
                          jk_are_docc_and_spin_is_zero = &
                               j_is_docc .and. k_is_docc .and. mod( j, 2 ) .ne. mod( k, 2 )

                          if( il_are_docc_and_spin_is_zero .or. jk_are_docc_and_spin_is_zero ) then
                             ! the spin integration is 0, so don't compute the term
                             compute_exchanged_erep_term = .false.
                          else
                             compute_exchanged_erep_term = .true.
                          endif

                          if( (compute_erep_term .and. erep_is_signif) &
                               .or. (compute_exchanged_erep_term .and. exchange_is_signif ) ) then

                             ! i,k,j,l are indexes for the spin orbitals. in the density cofactor compution,
                             ! dme_b(1) is the first row to strike, dme_k(1) is the first column to strike,
                             ! dme_b(2) is the second row to strike, and dme_k(2) is the second column to
                             ! strike

                             ! the density() routine call includes spin integration. outputs are
                             ! erep_density = \sum_q^{M_{sc}} d^2_{q,io+i,ko+k,jo+j,lo+l} w_{q,ikjl}
                             ! and
                             ! exchanged_erep_density = \sum_q^{M_{sc}} [
                             !                               d^2_{q,io+i,ko+k,jo+j,lo+l} w_{q,iljk} ]
                             ! if orbitals are docc, then spin was integrated with the variables above.
                             dme_b(1) = i ;  dme_k(1) = k ;  dme_b(2) = j ;  dme_k(2) = l

                             call density( 2, erep_density,exchanged_erep_density,&
                                  erep_int, exchanged_erep_int, erep_is_signif, exchange_is_signif )

                             ! erep_sum = \sum_i^{spins in io}\sum_j^{spins in jo<io}\sum_k^{spins in ko}
                             !           \sum_l^{spins in lo<ko} [ erep_density(io,i,ko,k,jo,j,lo,l)
                             !              * erep_int(io,jo,ko,lo) ]

                             if(erep_is_signif) erep_sum = erep_sum + erep_density * erep_int

                             if(exchange_is_signif) then
                                exchanged_erep_sum = exchanged_erep_sum &
                                     + exchanged_erep_density * exchanged_erep_int
                             endif
                          endif

                       end  if    !  k.ne.l (permutations)
                    end  do    !  l SO
                 end  do    !  k SO

              end  if    !  j<i (erep operator)
           end  do    !  j SO
        end  do    !  i SO

     else
     end  if    !  density for significant ERI


     !     accumulate the total energy
     if( ijkl_symmetry_is_enabled ) then
        if( ko .eq. io .and. lo .eq. jo ) then
           ! nothing
        else
           ! multiply by two to account for symmetry
           erep_sum = erep_sum*2.0_dp
           exchanged_erep_sum = exchanged_erep_sum*2.0_dp
        endif
     endif

     energy  =  energy  +  erep_sum - exchanged_erep_sum
     !     set up next task index
     task  =  task  +  nrank
  end  do    !  my tasks
end subroutine  vsvb_energy




!>    computes the overlap integrals and stores in densitywork::wdet
subroutine  wfndet
  use         densitywork
  use         integrals
  use tools, only: dp
  use xm
  implicit    none

  integer  :: task, i,j

  !     make the unsubstituted overlap square

  do   i = 1, nelec
     do   j = 1, nelec
        wdet( j, i ) = 0.0_dp
     end  do
  end  do


  !!  this is moderately inefficient for DOCC orbitals, but since
  !!  it is 1e- integrals and parallel, it is low priority

  task = 0
  do   i = 1, nelec
     do   j = 1, i
        task = task + 1
        if ( mod( task, nrank ) .eq. irank ) then
           call  ovint( bra( i ), ket( j ) )
           wdet( i, j ) = sint
           wdet( j, i ) = sint
        end  if
     end  do
  end  do

  call xm_equalize( wdet, nelec**2 )

#ifdef PRINT_MATRIX
  ! Writing matrix to a file
  call write_matrix(wdet,nelec,nelec,'det')
#endif

end subroutine  wfndet




!> compute integrals used in schwarz screening
!!  \param num_spatial_orbs [in]: number of "spatial" orbitals, not spin orbitals
!!                      (2*npair+nunpd+ndocc)
!!  \param num_non_docc [in]: number of non-docc orbitals (2*npair+nunpd)
subroutine  schwarz_ints ( num_spatial_orbs, num_non_docc )
  use tools, only: dp
  use         densitywork
  use         integrals
  use xm
  implicit    none

  integer     num_spatial_orbs, num_non_docc
  integer     nproc,myrank,master,task, i,j,indx,ij,ie,je


  ij  =  0
  do  i  =  1,  num_spatial_orbs
     do  j  =  1,  i
        ij  =  ij  +  1
        schwarz( ij )  =  0.0_dp
     end  do
  end  do


  task = 0
  do  i  =  1,  num_spatial_orbs
     do  j  =  1,  i
        task = task + 1
        if ( mod( task, nrank ) .eq. irank ) then
           ie  =  i  ;  if  (  i .gt. num_non_docc )  ie  =  2*i - num_non_docc
           je  =  j  ;  if  (  j .gt. num_non_docc )  je  =  2*j - num_non_docc
           call  int2e( bra( ie ), ket( je ), bra( ie ), ket( je ) )
           schwarz( indx( i, j ) )  =  sqrt( gint )
        end  if
     end  do
  end  do

  call xm_equalize( schwarz, (num_spatial_orbs**2 + num_spatial_orbs)/2 )
end  subroutine  schwarz_ints




!>    computes the density for the given spin orbitals in dme_[bk]
!!    \param nord [in]: 1 if it's a 1e density, 2 if it's a 2e density
!!    \param erep_density,exchanged_erep_density [out]: density cofactors.
!!           if nord == 1, only erep_density has output
!!    \param calc_dens, calc_exchange_dens [in] : logical flags controlling
!!            whether or not to calculate the density or exchanged density.
!!            (so for nord == 1, calc_exchange_dens should be false)
subroutine  density ( nord, erep_density,exchanged_erep_density,&
     erep_int, exchanged_erep_int, calc_dens, calc_exchange_dens )
  use tools, only: dp
  use         densitywork
  use         integrals
  implicit    none
  integer     nord,k,l,dima,dimb
  real(dp)    erep_density,exchanged_erep_density, erep_int, &
       exchanged_erep_int

  integer     isc,jsc
  real(dp)    erep_density_sc,exchanged_erep_density_sc
  logical     calc_dens, calc_exchange_dens

  if  ( spinopt ) then
     if  ( nord .eq. 1 ) then
        do isc = 1, nspinc
           do jsc = 1, isc
              call density_sc( nord, isc, jsc, erep_density_sc,&
                   exchanged_erep_density_sc, calc_dens, calc_exchange_dens )
              ham( isc, jsc ) = ham( isc, jsc ) + (erep_density_sc) * hint
              ovl( isc, jsc ) = ovl( isc, jsc ) + (erep_density_sc) * sint
           end do
        end do
     else if  ( nord .eq. 2 ) then
        do isc = 1, nspinc
           do jsc = 1, isc
              call density_sc( nord, isc, jsc, erep_density_sc,&
                   exchanged_erep_density_sc, calc_dens, calc_exchange_dens )
              ham( isc, jsc ) = ham( isc, jsc ) &
                   + erep_density_sc * erep_int &
                   - exchanged_erep_density_sc * exchanged_erep_int
           end do
        end do
     end if
  else
     if  (  nspinc  .gt.  0  )  then
        erep_density = 0.0_dp
        exchanged_erep_density = 0.0_dp
        !! densities are sums over M_{sc} ( = N_{sc}^2 2^{2N_p} ) terms for each spatial
        !! orbital integral
        do isc = 1, nspinc
           do jsc = 1, nspinc
              call density_sc( nord, isc, jsc, erep_density_sc,&
                   exchanged_erep_density_sc, calc_dens, calc_exchange_dens )

              erep_density = erep_density &
                   + erep_density_sc * coeff_sc( isc ) * coeff_sc( jsc )
              exchanged_erep_density = exchanged_erep_density &
                   + exchanged_erep_density_sc * coeff_sc( isc ) &
                   * coeff_sc( jsc )

           end do
        end do
     else

        dima=npair+ndocc+nunpd
        dimb=npair+ndocc
        erep_density = 0.0_dp
        exchanged_erep_density = 0.0_dp
        call build_abket(dima,dimb)
        call det( dima, dimb, nord,erep_density,exchanged_erep_density,&
             calc_dens, calc_exchange_dens )
     end if
  end if
end subroutine  density


!>    computes the density for the given spin orbitals in dme_[bk] and
!!    between the isc'th and jsc'th spin coupling
!!    \param nord [in]: 1 if it's a 1e density, 2 if it's a 2e density
!!    \param isc, jsc [in]: spin couplings, only meaningful if nspinc>0
!!    \param erep_density,exchanged_erep_density [out]: density cofactors.
!!           if nord == 1, only erep_density has output
!!    \param calc_dens, calc_exchange_dens [in] : logical flags controlling
!!            whether or not to calculate the density or exchanged density.
!!            (so for nord == 1, calc_exchange_dens should be false)
subroutine  density_sc ( nord, isc, jsc, erep_density,exchanged_erep_density,&
     calc_dens, calc_exchange_dens )
  use tools, only: dp
  use         densitywork
  use         integrals
  implicit    none
  integer     nord, isc,jsc, previously_computed
  real(dp)    erep_density,exchanged_erep_density

  integer     i,j,nexb,nexk, dima,dimb
  logical     calc_dens, calc_exchange_dens

  erep_density = 0.0_dp
  exchanged_erep_density = 0.0_dp
  nexb=npair
  nexk=npair
  dima=npair+ndocc+nunpd
  dimb=npair+ndocc

  call  dbra( nord, nexb,nexk, dima,dimb, erep_density,&
       exchanged_erep_density,isc,jsc,calc_dens, calc_exchange_dens )
end subroutine density_sc


!>    for the isc, jsc spin coupling pair, sets up mapping matrices to generate
!!    the 2^{N_p} density terms in the bra, and for each, calls dket (which
!!    the generates 2^{N_p} density terms for permutations in the ket)
!!    \param nord [in]: 1 if it's a 1e density, 2 if it's a 2e density
!!    \param nexb, nexk [in]: number of spin coupled pairs in the bra and ket
!!    \param dima, dimb [in]: number of alpha and beta electrons
!!    \param erep_density,exchanged_erep_density [in/out] density cofactors.
!!           if nord == 1, only erep_density has output
!!    \param isc, jsc [in]: spin couplings, only meaningful if nspinc>0
!!    \param calc_dens, calc_exchange_dens [in] : logical flags controlling
!!            whether or not to calculate the density or exchanged density.
!!            (so for nord == 1, calc_exchange_dens should be false)
subroutine  dbra ( nord, nexb,nexk, dima,dimb, erep_density,&
     exchanged_erep_density,isc,jsc,calc_dens, calc_exchange_dens )
  use tools, only: dp
  use         densitywork
  use         integrals
  implicit    none
  integer     nord, nexb,nexk, dima,dimb,isc,jsc
  real(dp)    erep_density,exchanged_erep_density

  integer    io,i,j,k,l,tmp
  logical    calc_dens, calc_exchange_dens


  !     generate the unique determinent lists with the ket list nested
  !     inside the bra list - this avoids storing the 2**Npair list
  !     begin with the identity det

  ! get the mapping of alpha and beta electrons to OLCAO for the given
  ! spin coupling
  do k = 1, npair
     bra_a(k) = pair_sc( k, 1, isc )
     bra_b(k) = pair_sc( k, 2, isc )
  end do

  ! split overlap matrix into alpha and beta parts, and store
  ! as the OLCAO associated with each spin coupled alpha (beta) e-
  ! with every other spin orbital
  do l = 1, nelec
     do k = 1, npair
        abra_npair( k, l ) = wdet( bra_a( k ), l )
        bbra_npair( k, l ) = wdet( bra_b( k ), l )
     end do
  end do

  call  dket( nord, nexk,dima,dimb, erep_density,exchanged_erep_density,&
       isc,jsc,calc_dens, calc_exchange_dens )

  !     for the remaining bra terms - mimics the action of the
  !     spin function - offset from DME labels

  do io = 1, nexb

     !     note the combinatorial algorithm can proceed from an
     !     arbitrary starting point within the list
     !     this is the first combination of 'io'th order

     do k = 1, io
        bexch(k) = k
     end do

     !     remake the bra and apply exchange(s)

     do k = 1, npair
        bra_a(k) = pair_sc( k, 1, isc )
        bra_b(k) = pair_sc( k, 2, isc )
     end do

     do k = 1, io
        tmp = bra_a(k)
        bra_a(k) = bra_b(k)
        bra_b(k) = tmp
     end do

     do l = 1, nelec
        do k = 1, npair
           abra_npair( k, l ) = wdet( bra_a( k ), l )
           bbra_npair( k, l ) = wdet( bra_b( k ), l )
        end do
     end do

     call  dket( nord, nexk,dima,dimb, erep_density,exchanged_erep_density,&
          isc,jsc,calc_dens, calc_exchange_dens )

     !     loop over remaining combinations of 'io'th order

     i = io
     do while (bexch(1).lt.1+nexb-io)
        if (bexch(i).lt.i+nexb-io) then
           bexch(i) = bexch(i) + 1
           do k = i+1, io
              bexch(k) = bexch(k-1) + 1
           end do

           !     remake the bra and apply exchange(s)

           do k = 1, npair
              bra_a(k) = pair_sc( k, 1, isc )
              bra_b(k) = pair_sc( k, 2, isc )
           end do

           do k = 1, io
              l = bexch(k)
              tmp = bra_a(l)
              bra_a(l) = bra_b(l)
              bra_b(l) = tmp
           end do

           do l = 1, nelec
              do k = 1, npair
                 abra_npair( k, l ) = wdet( bra_a( k ), l )
                 bbra_npair( k, l ) = wdet( bra_b( k ), l )
              end do
           end do

           call  dket( nord, nexk,dima,dimb, erep_density,exchanged_erep_density,&
                isc, jsc,calc_dens, calc_exchange_dens )

           i = io
        else
           i = i - 1
        end if
     end do   !  while i
  end do   !  io
end subroutine dbra


!>    for the isc, jsc spin coupling pair, and the mapping already set up in
!!    bra_[ab], this routine sets up mapping matrices to generate
!!    the 2^{N_p} density terms in the ket, and calls det() to calculate the
!!    density.
!!    \param nord [in]: 1 if it's a 1e density, 2 if it's a 2e density
!!    \param nexk [in]: number of spin coupled pairs in the bra and ket
!!    \param dima, dimb [in]: number of alpha and beta electrons
!!    \param erep_density,exchanged_erep_density [in/out] density cofactors.
!!           if nord == 1, only erep_density has output
!!    \param isc, jsc [in]: spin couplings, only meaningful if nspinc>0
!!    \param calc_dens, calc_exchange_dens [in] : logical flags controlling
!!            whether or not to calculate the density or exchanged density.
!!            (so for nord == 1, calc_exchange_dens should be false)
subroutine  dket ( nord, nexk,dima,dimb, erep_density,&
     exchanged_erep_density,isc,jsc,calc_dens, calc_exchange_dens )
  use tools, only: dp
  use         densitywork
  use         integrals
  implicit    none
  integer     nord, nexk,dima,dimb,isc,jsc
  real(dp)    erep_density,exchanged_erep_density

  integer     jo,i,j,k,l,tmp
  logical     calc_dens, calc_exchange_dens

  !     'ket' list - these operations 'mirror' the bra ones

  !     begin with the identity det

  do k = 1, npair
     ket_a(k) = pair_sc( k, 1, jsc )
     ket_b( k ) = pair_sc( k, 2, jsc )
  end do

  call build_abket( dima, dimb )

  call  det( dima,dimb,nord, erep_density,exchanged_erep_density,&
       calc_dens, calc_exchange_dens )

  !     apply 'exchange' operators between alpha,beta subdets
  !     for the remaining ket terms - mimics the action of the
  !     spin function - offset from the DME labels
  !     loop over numbers of exchanges to apply

  do jo = 1, nexk

     !     this is the first combination of 'jo'th order

     do k = 1, jo
        kexch(k) = k
     end do

     !     remake the ket and apply exchange(s)

     do k = 1, npair
        ket_a(k) = pair_sc( k, 1, jsc )
        ket_b( k ) = pair_sc( k, 2, jsc )
     end do

     do k = 1, jo
        tmp = ket_a(k)
        ket_a(k) = ket_b(k)
        ket_b(k) = tmp
     end do

     call build_abket( dima, dimb )

     call  det( dima,dimb,nord, erep_density,exchanged_erep_density,&
          calc_dens, calc_exchange_dens )

     !     loop over remaining combinations of 'jo'th order

     j = jo
     do while (kexch(1).lt.1+nexk-jo)
        if (kexch(j).lt.j+nexk-jo) then
           kexch(j) = kexch(j) + 1
           do k = j+1, jo
              kexch(k) = kexch(k-1) + 1
           end do

           !     remake the ket and apply exchange(s)

           do k = 1, npair
              ket_a(k) = pair_sc( k, 1, jsc )
              ket_b( k ) = pair_sc( k, 2, jsc )
           end do

           do k = 1, jo
              l = kexch(k)
              tmp = ket_a(l)
              ket_a(l) = ket_b(l)
              ket_b(l) = tmp
           end do


           call build_abket( dima, dimb )

           call  det( dima,dimb,nord, erep_density,exchanged_erep_density,&
                calc_dens, calc_exchange_dens )

           j = jo
        else
           j = j - 1
        end if
     end do   !  while j
  end do   !  jo
end subroutine dket



!>    perform spin integration and calculate density determinant
!!
!!    if nord == 1,  returns first order cofactor
!!    d^1_{dme_b(1),dme_k(1)} w_{dme_b(1),dme_k(1)} in density
!!
!!    If nord == 2,  returns second order cofactor
!!    d^2_{dme_b(1),dme_k(1),dme_b(2),dme_k(2)}
!!    w_{dme_b(1),dme_b(2),dme_k(1),dme_k(2)} in density and
!!    returns second order cofactor d^2_{dme_b(1),dme_k(1),dme_b(2),dme_k(2)}
!!    w_{dme_b(1),dme_b(2),dme_k(2),dme_k(1)} in exchanged_density.
!!
!!    (the densities are equivalent to cofactors, since dme_b(2) < dme_b (1)
!!    and dme_k(2) < dme_k(1) )
!!    \param dima, dimb [in] : number of alpha and beta electrons
!!    \param nord [in] : 1 if it's a 1e density, 2 if it's
!!           a 2e density
!!    \param density,exchanged_density [in/out] : density cofactors
!!           if nord == 1, only density has output
!!    \param calc_dens, calc_exchange_dens [in] : logical flags controlling
!!            whether or not to calculate the density or exchanged density.
!!            (so for nord == 1, calc_exchange_dens should be false)
subroutine  det ( dima,dimb,nord, density,exchanged_density,&
     calc_dens, calc_exchange_dens )
  use tools, only: dp
  use         densitywork
  use         integrals
  implicit    none
  integer     dima,dimb,nord
  real(dp)    density,exchanged_density

  logical     erep_spin_is_nonzero, exchange_spin_is_nonzero, &
       calc_exchange_dens,calc_dens

  logical both_are_alpha_spin, both_are_beta_spin, &
       density_setup_with_erep, density_setup_with_exchange
  integer     loc_alpha_bra,loc_alpha_ket,loc_beta_bra,loc_beta_ket

  integer     iord,i

  real(dp)    zero, one, ad,bd
  parameter  ( zero = 0.0d+00, one = 1.0d+00 )

  erep_spin_is_nonzero = .false.
  exchange_spin_is_nonzero = .false.
  ad = 0.0_dp
  bd = 0.0_dp
  !     perform spin integration

  !  dme_b(1) is the row to strike, dme_k(1) is the column to strike
  ! for second order cofactors, dme_b(2) is the second row to strike,
  ! and dme_k(2) is the second column to strike
  if( calc_dens ) then
     do iord = 1, nord

        ! returns true/false into both_are_alpha_spin, both_are_beta_spin depending on
        ! if the spin orbitals dme_b( iord ) and dme_k( iord ) are both alpha,
        ! both beta, or neither, and
        ! returns the location of the spin orbital dme_b(iord) in the
        ! alpha or beta submatrix of the overlap matrix (whichever is appropriate).
        call check_spin_and_locate( dme_b( iord ), dme_k( iord ), bra_a,&
             bra_b, ket_a, ket_b, nalpha,nbeta, &
             loc_alpha_bra, loc_alpha_ket, loc_beta_bra,&
             loc_beta_ket,both_are_alpha_spin, both_are_beta_spin )

        if( both_are_alpha_spin ) then
           !     electron spins both alpha: 'unintegrate' the spatial variables
           do i = 1, nalpha
              aket( loc_alpha_bra, i ) = zero
              aket( i, loc_alpha_ket ) = zero
           end do
           aket( loc_alpha_bra, loc_alpha_ket ) = one

        else if( both_are_beta_spin ) then
           !     electron spins both beta ...
           do i = 1, nbeta
              bket( loc_beta_bra, i ) = zero
              bket( i, loc_beta_ket ) = zero
           end do
           bket( loc_beta_bra, loc_beta_ket ) = one
        else
           !            spin integral is zero
           exit
        end if

     end do     !  iord

     erep_spin_is_nonzero = both_are_beta_spin .or. both_are_alpha_spin

  endif

  density_setup_with_erep = calc_dens .and. erep_spin_is_nonzero

  ! if it's a density for a 2-electron integral, we need the exchange term,too
  if( calc_exchange_dens ) then

     do iord = 1, nord

        call check_spin_and_locate( dme_b( iord ), dme_k( nord-iord+1 ), &
             bra_a, bra_b, ket_a,ket_b, nalpha,nbeta,&
             loc_alpha_bra, loc_alpha_ket,&
             loc_beta_bra,loc_beta_ket,both_are_alpha_spin, both_are_beta_spin )

        if( both_are_alpha_spin ) then
           !     electron spins both alpha: 'unintegrate' the spatial variables
           ! only do this if the first term didn't
           if( .not. density_setup_with_erep ) then
              do i = 1, nalpha
                 aket( loc_alpha_bra, i ) = zero
                 aket( i, loc_alpha_ket ) = zero
              end do
              aket( loc_alpha_bra, loc_alpha_ket ) = one
           endif
        else if( both_are_beta_spin ) then
           !     electron spins both beta ...
           if( .not. density_setup_with_erep ) then
              do i = 1, nbeta
                 bket( loc_beta_bra, i ) = zero
                 bket( i, loc_beta_ket ) = zero
              end do
              bket( loc_beta_bra, loc_beta_ket ) = one
           endif
        else
           !     spin integral is zero
           exit
        end if
     end do     !  iord

     exchange_spin_is_nonzero = both_are_alpha_spin .or. both_are_beta_spin

  endif

  ! now you know whether erep or exchange is necessary.
  density_setup_with_exchange = calc_exchange_dens .and. &
       exchange_spin_is_nonzero

  if( density_setup_with_erep .or. density_setup_with_exchange ) then
     !     determinant survives the spin integration

     if (dima.eq.0) then    !  possible with spin-free method
        ad = one
     else
        if (dima.eq.1) then
           ad = aket(1,1)
        else
           call givdr( padded_size, nalpha, aket, dtol, ad, ipvt )
        end if
     end if
     if (dimb.eq.0) then
        bd = one
     else
        if (dimb.eq.1) then
           bd = bket(1,1)
        else
           call givdr( padded_size, nbeta, bket, dtol, bd, ipvt )
        end if
     end if

  endif

  if( density_setup_with_erep ) density = density + ad*bd

  ! if density_setup_with_erep is true, then [ab]ket was filled using
  ! row i, col l, row j, col k, so we need to flip the sign (since the columns
  ! were flipped from what we want, and the determinant of a matrix when you flip
  ! two columns is -1 * the determinant of the original matrix.
  if( density_setup_with_exchange &
       .and. .not. density_setup_with_erep ) then
     exchanged_density = exchanged_density - ad*bd
  endif

  if( density_setup_with_exchange .and. density_setup_with_erep ) then
     exchanged_density = exchanged_density + ad*bd
  endif

end subroutine det




!> driver to calculate a determinant via givens rotations
!! \param max_n [in]: sizing, in case memory doesn't match
!! \param n [in]: size of adet
!! \param adet [in]: matrix to calculate determinant of
!! \param tol [in]: tolerance
!! \param d [out]: determinant of adet
subroutine  givdr ( max_n, n, adet, tol, d, ipvt )
  use timing_flops
  use tools, only: dp
  use xm, only: write_determinant
  implicit    none
#ifdef PRINT_TIMING
  include "mpif.h"
#endif
  integer     max_n, n, ipvt(*)
  real(dp)      adet( max_n, *), d

  integer     i
  real(dp)      tol
  real(dp) t1, t2


#ifdef PRINT_TIMING
  t1 = MPI_Wtime()
#endif

#ifdef PRINT_COUNTERS
  count_determinants = count_determinants+1
#endif

#ifdef USE_C_VERSION
  call givensc( adet, max_n, n, tol )
  d = 1.0_dp
  do i = 1, n
     d = d * adet(i,i)
  end do

#elif USE_FORTRAN_VERSION
  call givens_single( adet, max_n, n, tol )
  !     the product of the diagonal elements is det(adet)

  d = 1.0_dp
  do i = 1, n
     d = d * adet(i,i)
  end do
#elif GIVENS_KNL
  call givens( adet, max_n, n, tol )
  !     the product of the diagonal elements is det(adet)

  d = 1.0_dp
  do i = 1, n
     d = d * adet(i,i)
  end do

#elif USE_LAPACK
  do i = 1, n
     ipvt(i) = 0
  end do
  call determinant_lapack( adet, max_n, n, d, ipvt )

#elif GIVENS_BGQ
  call givens_bgq( adet, max_n, n, tol )
  !     the product of the diagonal elements is det(adet)

  d = 1.0_dp
  do i = 1, n
     d = d * adet(i,i)
  end do

#else
  print *, "should never happen"
#endif

#ifdef PRINT_MATRIX
  call write_determinant(d)
#endif

#ifdef PRINT_TIMING
  t2 = MPI_Wtime()
#endif

  kernel_time = kernel_time + (t2-t1)

end subroutine givdr


!     normalisation of LCAO coefficients
!!
!!     this routine normalizes orbitals
!!     between the two inputs (so different subgroups of
!!     orbitals can be normalized separately).
!!     Normalizes orbital "i" by multiplying all the
!!     basis function coefficients in the orbital expansion by
!!     1/ sqrt(< i | i >).
!!
!!     \param ist,ind the beginning and ending of the orbitals to normalize
subroutine  normal ( ist, ind )
  use         integrals
  implicit    none
  integer     ist, ind

  integer     i,j


  do i = ist, ind

     !     1. compute the self-overlap integral

     call ovint( i, i )

     !     2. go back and scale the coefficients

     ! calculate 1 / sqrt( S_{ii} ) for given i
     ! loop through all basis functions coefficients associated with the
     ! orbital and weight by 1 / sqrt( S_{ii} )
     sint = sint**( -0.5_dp )
     do j = map_orbs( i ), map_orbs( i + 1 ) - 1
        coeff( j ) = coeff( j )*sint
     end do

  end do
end subroutine normal





!     assume the NDF OBS is a subset of the target orbital OBS
!     - this must be set in the input

subroutine  ndf2obs ( iorb, indf )
  use         integrals
  implicit    none

  integer     indf,iorb,mnshi,mxshi,it,jt,ib,jb,ia,iat,jat,ish,i,j


  !     NDF atom AO ranges avoid strict ordering of input AOs

  ib = 0
  do  iat = 1,  orbas_atnum( indf )
     atom_ndf( 1, iat )  =  ib  +  1
     it  =  atom_t( orbas_atset( iat, indf ) )
     mnshi  =         map_atom2shell( it )
     mxshi  = mnshi + num_shell_atom( it ) - 1
     do  ish = mnshi, mxshi
        ib = ib + shell_size( ang_mom( ish ) )
     end do
     atom_ndf( 2, iat )  =  ib
  end do


  !     map the OBS of the NDF to that of the orbital using it

  ib = 1
  do  iat = 1,  orbas_atnum( iorb )
     do  jat = 1,  orbas_atnum( indf  )
        if  (  orbas_atset( jat, indf  ) .eq.    &
             orbas_atset( iat, iorb )  )  ndf2orb( jat )  =  ib
     end do
     it =  atom_t( orbas_atset( iat, iorb ) )
     mnshi  =         map_atom2shell( it )
     mxshi  = mnshi + num_shell_atom( it ) - 1
     do  ish = mnshi, mxshi
        ib = ib + shell_size( ang_mom( ish ) )
     end do
  end do


  !     map the NDF expansion

  i  =  0
  do ib  =  map_orbs( indf ), map_orbs( indf + 1 ) - 1

     !     find which atom this AO belongs to

     do ia = 1, orbas_atnum( indf )
        if  (  xpset( ib ) .ge. atom_ndf( 1, ia )  .and.   &
             xpset( ib ) .le. atom_ndf( 2, ia )  )   &
             j  =  ndf2orb( ia )  -  atom_ndf( 1, ia )
     end  do
     i  =  i  +  1
     xpnew( i )  =  xpset( ib )  +  j

  end do      !   AOs
end subroutine  ndf2obs






!> normalize a shell of primitive GTO functions
!!
!! \param ang_mom [in] : angular momentum of the shell (and primitive GTOs)
!! \param con_length [in] : number of primitives in the shell--dimension of
!!  exponent and con_coeff
!! \param exponent [in] : array of exponents of primitive GTOs in the shell
!! \param con_coeff [in/out] : array of coefficients of primitive GTOs in the
!!                             shell
subroutine norm_prim (  ang_mom, con_length,   &
     exponent, con_coeff )
  use tools, only: dp
  implicit   none
  integer    ang_mom,  con_length
  real(dp)     exponent(*),  con_coeff(*)

  integer    ig, jg
  real(dp)     two, dblfac, pi32, sovl, fac, fax
  parameter  ( two  =  2.0_dp )

  !     multiply the coefficients of the primitive functions by their
  !     normalization factors (= 1/sqrt( self-overlap integral ))
  !     the 'PT...' factor is (2L-1)!!/(2**L), 'pi32' = (Pi)**3/2

  pi32 = ( acos(-1.0_dp) )**( 1.5_dp )
  fac = pi32*dblfac( 2*ang_mom - 1 )*( two**( -ang_mom ) )
  fax = -1.5_dp -ang_mom
  do ig = 1,  con_length
     sovl  =  fac*( ( two * exponent( ig ) )**fax )
     con_coeff( ig )  =  con_coeff( ig )*( sovl**( -0.5_dp ) )
  end do

  !     compute the normalization factor for the whole CGTO
  !     using our freshly normalized coefficients
  !     not bothering to exploit shell-pair label symmetry

  sovl  =  0.0_dp
  do ig =  1,  con_length
     do jg =  1,  con_length
        sovl  =  sovl + fac*con_coeff( ig )*con_coeff( jg ) *  &
             ( ( exponent( ig ) + exponent( jg ) )**fax )
     end do
  end do

  !     multiply the same coefficients by the CGTO normalization factor

  sovl  = sovl**( -0.5_dp )
  do ig =  1,  con_length
     con_coeff( ig )  =  con_coeff( ig )*sovl
  end do
end subroutine norm_prim




!     this function returns the 'double factorial' number,
!     n!! = 1.3.5.7.9...n, for odd n
!     it assumes the input n is the odd integer it should be

real(dp)    function  dblfac( n )
  use tools, only : dp
  implicit  none
  integer   n, i
  dblfac  =  1.0_dp
  do    i  =  3,  n, 2
     dblfac  =  dblfac * dble( i )
  end   do
end function dblfac




!>    fills array nxyz with powers of x,y,z coordinates of primitive functions
!!    in order:
!!
!!    nxyz(1-3, 1) = 0 (since x, y, z powers are 0 for the first
!!    primitive s shell.
!!
!!    nxyz(1-3, 2-4) = 1,0,0 then 0,1,0, then 0,0,1 since the power
!!    of x,y,z are each one for the 2-4 primitive functions
!!    \param nang [in]: highest angular momentum in the basis set
!!    \param nxyz [in]: output array, filled
subroutine   cartesian ( nang, nxyz )
  implicit     none
  integer      nang, nxyz(3,*)

  integer      ip,ix,iy,iz,l,il,i
  integer n,j
  nxyz(1,1) = 0
  nxyz(2,1) = 0
  nxyz(3,1) = 0
  ip = 1
  do l = 1, nang
     iz = 0
     do il = l + 1, 1, -1
        ix = il
        iy = 0
        do i = 1, il
           ip = ip + 1
           ix = ix - 1
           nxyz(1,ip) = ix
           nxyz(2,ip) = iy
           nxyz(3,ip) = iz
           iy = iy + 1
        end do
        iz = iz + 1
     end do
  end do

! For CCA ordering:
#ifdef CCA_ORDER
  n = 2
  do l=1,nang
  do i = 0,l
     do j = 0,i
        nxyz(1,n) = l - i
        nxyz(2,n) = i - j
        nxyz(3,n) = j
!        print *, nxyz(1,n),nxyz(2,n),nxyz(3,n)
        n = n + 1
     end do
  enddo
  enddo
#endif

end subroutine cartesian




!>    fills arrays with coefficients of primitive cartesian GTOs
!!
!!    fills ashl(i) with sqrt*( (2*i-1)!! ) for each angular momentum i=0,nang
!!
!!    fills ashi(i) with 1/( ashl(i) )
!!
!!    fills angn(ij) with ashl(i) * 1/(ashl(power of x)) * 1/(ashl(power of y)) *
!!    1/(ashl((power of z)) for each primitive
!!    in each angular momentum i,where (power of x) + (power of y) + (power of z) = i)
!!    and ij walks over all primitives--1 for s, 3 for p, etc.
subroutine  setangn
  use tools, only: dp
  use         integrals
  implicit    none

  integer    i,j,ij
  real(dp)   one, fii
  parameter (one=1.0d+00)

  ij = 0
  do i = 0, nang
     fii = one
     do j = 1, 2*i-1, 2
        fii = fii*dble(j)
     end do
     fii = sqrt(fii)
     ashl(i) = fii
     ashi(i) = one/fii
     do j = 1, shell_size( i )
        ij = ij + 1
        angn(ij) = ashl(i)  &
             * ashi( nxyz( 1, ij ) )  &
             * ashi( nxyz( 2, ij ) )  &
             * ashi( nxyz( 3, ij ) )
     end do
  end do
end subroutine setangn





integer function indx ( i, j )
  implicit none
  integer  i,j,ix
  ix   = max( i, j )
  indx = (ix*ix-ix)/2 + min( i, j )
end function indx




!> Sets up the bra_[ab], ket_[ab], [ab]bra arrays for the
!! docc and unpaired electrons
!! since these arrays never change during a run.
!! the goal is to set up the
!! bra_[ab] and ket_[ab] contain indexes mapping alpha (beta) electrons to
!! their spin functions indices in the bra and ket part of the
!! wavefunction
subroutine set_up_unpaired_docc
  use         densitywork
  use         integrals
  implicit none
  integer i,j,dima, dimb,num_sc_elec


  dima=npair+ndocc+nunpd
  dimb=npair+ndocc
  num_sc_elec = npair*2

  ! for a given alpha electron, store the spin orbital associated with it for
  ! the bra and ket in bra_[ab] and ket_[ab]
  ! they are used to make the appropriate the density matrix for a given integral

  ! start from the unpaired electrons (after the 2*npair spin coupled electrons)
  ! (where we treat one spin coupling at a time).
  j = 2*npair + 1
  do i =   npair + 1, npair + nunpd
     bra_a( i ) = j
     ket_a( i ) = j
     j = j + 1
  end do
  !     lastly, the DOCC electrons
  j = 2*npair + nunpd + 1
  do i =   npair + nunpd + 1, npair + nunpd + ndocc
     bra_a( i ) = j
     ket_a( i ) = j
     j = j + 2
  end do

  j = 2*npair + nunpd + 2
  do i =   npair + 1, npair + ndocc
     bra_b( i ) = j
     ket_b( i ) = j
     j = j + 2
  end do

  ! this is the first step in forming the matrix which we want the determinant of
  ! (that is, [ab]ket). the first step involves splitting the overlap matrix into
  ! alpha and beta blocks and stores the overlap matrix in terms of the spin orbital
  ! associated with each non-spin coupled alpha (beta) e- and every spin coupled
  ! electron.
  ! uses bra_[ab] to map from full overlap matrix to alpha and beta blocks.
  do j = 1, num_sc_elec
     do i = 1, dima-npair
        abra_docc_un( i, j ) = wdet( bra_a( i+npair ), j )
     end do
     do i = 1, dimb-npair
        bbra_docc_un( i, j ) = wdet( bra_b( i +npair ), j )
     end do
  end do


  ! this is the second step, which forms the non-spin coupled (docc and unpaired)
  ! part of [ab]ket, the matrix whose determinant is the density. the non-spin
  ! coupled part is prepared here, since it will never change during the
  ! run. later, this will be copied into [ab]ket each time before the determinant
  ! is calculated. the second step for the sc e- is calculated elsewhere.
  ! uses bra_[ab] and ket_[ab]to map from full overlap matrix to alpha and beta
  ! blocks.
  do i = 1, dima-npair
     do j = 1, dima-npair
        aket_docc_un( j, i ) = wdet( bra_a(j+npair), ket_a( i+npair ) )
     end do
  end do
  do i = 1, dimb-npair
     do j = 1, dimb-npair
        bket_docc_un( j, i ) = wdet( bra_b(j+npair), ket_b( i+npair ) )
     end do
  end do


end subroutine set_up_unpaired_docc

!> builds the [ab]ket matrix
!! using the overlaps in [ab]bra and the configurarion in ket_[ab],
!! [ab]ket is formed, so the determinant can be calculated,
!! and the density element computed
!! \param dima,dimb [in]: number of alpha and beta electrons
subroutine build_abket( dima, dimb )
  use densitywork
  use integrals
  implicit none
  integer dima, dimb,k,l, num_sc_elec

  !      dima=npair+ndocc+nunpd
  !      dimb=npair+ndocc
  num_sc_elec = npair*2
  ! finish forming [ab]ket, so the determinant of it can be computed

  ! first for the sections of [ab]ket that involve spin-coulped
  ! e-. these sections change with each different spin-coupled pair
  ! ket_[ab] map the electron in [ab]bra_docc_un and [ab]bra_npair to
  ! the appropriate OLCAO for this determinant calculation
  if( npair .gt. 0 ) then
     do k = 1, npair
        do l = 1+npair, dima
           aket( l, k ) = abra_docc_un( l-npair, ket_a( k ) )
        end do
        do l = 1+npair, dimb
           bket( l, k ) = bbra_docc_un( l-npair, ket_b( k ) )
        end do
     end do

     do k = 1, dima
        do l = 1, npair
           aket( l, k ) = abra_npair( l, ket_a( k ) )
        end do
     end do

     do k = 1, dimb
        do l = 1, npair
           bket( l, k ) = bbra_npair( l, ket_b( k ) )
        end do
     end do

  endif


  ! recopy over the DOCC and unpaired e- part of [ab]ket, since they were
  ! destroyed in the last determinant calculation
  do k = 1, ndocc

     do l = 1+npair, npair+ndocc
        aket(l, k+npair ) = aket_docc_un( l-npair, k )
        bket( l, k+npair ) = bket_docc_un( l-npair, k )
     end do

     do l = 1+npair+ndocc, npair+nunpd+ndocc
        aket( l, k+npair ) = aket_docc_un( l-npair, k )
     end do

  end do

  do k = 1+ndocc, ndocc+nunpd
     do l = 1+npair, npair+ndocc+nunpd
        aket( l, k+npair ) = aket_docc_un( l-npair, k )
     end do
  end do


end subroutine build_abket

#ifdef USE_LAPACK
!> calculates the determinant of adet
!! \param adet [in]: matrix of (n,n), with leading dimension (max_n)
!! \param max_n [in]: leading dimension of adet
!! \param n [in]: dimension of adet
!! \param d [out]: determinant
!! \param ipvt [in]: workspace for the pivots from dgetrf
subroutine determinant_lapack( adet, max_n, n, d, ipvt )
  use tools, only : dp
  implicit none
  real(dp) adet( max_n, * ), d, det(2)
  real(dp) ten
  integer  max_n, n, i, ierr, ipvt(n)

  call dgetrf(n,n,adet,max_n,ipvt,ierr)
  ! if (ierr .lt. 0)then
  !    print *,'illegal value detected in dgetrf = ',ierr
  !    stop
  ! endif

  if (ierr .gt. 0) then
     d=0.0_dp
     return
  endif
  !
  ! compute the determinant
  !
  det(1) = 1.0_dp
  det(2) = 0.0_dp
  ten = 10.0_dp
  do i = 1, n
     if (ipvt(i) .ne. i) det(1) = -det(1)
     det(1) = adet(i,i)*det(1)
     !        ...exit
     if (det(1) .eq. 0.0_dp) go to 60
     ! do while( dabs(det(1)) .lt. 1.0_dp )
     !    det(1) = ten*det(1)
     !    det(2) = det(2) - 1.0_dp
     ! enddo

10   if (dabs(det(1)) .ge. 1.0_dp) go to 20
     det(1) = ten*det(1)
     det(2) = det(2) - 1.0_dp
     go to 10
20   continue

30   if (dabs(det(1)) .lt. ten) go to 40
     det(1) = det(1)/ten
     det(2) = det(2) + 1.0_dp
     go to 30
40   continue
  enddo
60 continue
  d = det(1) * 10.0**det(2)
end subroutine determinant_lapack
#endif
!> checks the spin of the given two spin orbitals (one in bra and one in ket)
!! and returns logicals denoting the spin as well as locations of the
!! spin orbitals in the bra and ket submatrices ([ab]ket of the full overlap
!! matrix
!!
!! this routine requires bra_[ab] and ket_[ab] to be filled properly.
!! \param spin_orb_in_bra, spin_orb_in_ket [in] : indexes for a spin orbital in the
!!        bra and another in the ket
!! \param bra_alpha, bra_beta, ket_alpha, ket_beta [in] : mapping arrays that
!!        take in an alpha/beta electron and return the location of the electron
!!        in the full bra/ket wavefunction.
!! \param nalpha, nbeta [in] : number of alpha and beta electrons
!! \param orb_in_bra_alpha_index [out] : if spin_orb_in_bra has alpha spin, this holds
!!        the position of spin_orb_in_bra in the alpha submatrix of the overlap
!         matrix
!! \param orb_in_bra_beta_index [out] : if spin_orb_in_bra has beta spin, this holds
!!        the position of spin_orb_in_bra in the beta submatrix of the overlap
!!        matrix
!! \param orb_in_ket_alpha_index [out] : if spin_orb_in_ket has alpha spin, this holds
!!        the position of spin_orb_in_ket in the alpha submatrix of the overlap
!!        matrix
!! \param orb_in_ket_beta_index [out] : if spin_orb_in_ket has beta spin, this holds
!!        the position of spin_orb_in_ket in the beta submatrix of the overlap
!!        matrix
!! \param both_are_alpha_spin, both_are_beta_spin [out] : logicals that flag
!!        whether spin_orb_in_bra and spin_orb_in_ket are both alpha or both beta.
subroutine check_spin_and_locate( spin_orb_in_bra, spin_orb_in_ket, bra_alpha,&
     bra_beta, ket_alpha,ket_beta, nalpha,nbeta, &
     orb_in_bra_alpha_index, orb_in_ket_alpha_index, orb_in_bra_beta_index,&
     orb_in_ket_beta_index, both_are_alpha_spin, both_are_beta_spin )
  implicit none
  logical     orb_in_bra_is_alpha, orb_in_bra_is_beta, &
       orb_in_ket_is_alpha, orb_in_ket_is_beta
  integer     i, spin_orb_in_bra, spin_orb_in_ket, orb_in_bra_alpha_index, orb_in_ket_alpha_index, &
       orb_in_bra_beta_index, orb_in_ket_beta_index, nalpha, nbeta
  integer     bra_alpha(*),bra_beta(*), ket_alpha(*),ket_beta(*)

  logical both_are_alpha_spin, both_are_beta_spin

  orb_in_bra_is_alpha = .false.
  orb_in_ket_is_alpha = .false.

  do i = 1, nalpha
     if ( spin_orb_in_bra .eq. bra_alpha( i ) ) then
        orb_in_bra_is_alpha = .true.
        orb_in_bra_alpha_index = i
     end if
     if ( spin_orb_in_ket .eq. ket_alpha( i ) ) then
        orb_in_ket_is_alpha = .true.
        orb_in_ket_alpha_index = i
     end if
  end do

  orb_in_bra_is_beta = .false.
  orb_in_ket_is_beta = .false.

  do i = 1, nbeta
     if ( spin_orb_in_bra .eq. bra_beta( i ) ) then
        orb_in_bra_is_beta = .true.
        orb_in_bra_beta_index = i
     end if
     if ( spin_orb_in_ket .eq. ket_beta( i ) ) then
        orb_in_ket_is_beta = .true.
        orb_in_ket_beta_index = i
     end if
  end do

  both_are_alpha_spin = orb_in_bra_is_alpha .and. orb_in_ket_is_alpha
  both_are_beta_spin = orb_in_bra_is_beta .and. orb_in_ket_is_beta

end subroutine check_spin_and_locate


!> minimize energy
subroutine minimize_energy( energy,  &
     w, eig, v1, v2, coefflock,int_out, dbl_out )
  use          integrals
  use tools, only: dp
  use xm
  implicit none
  logical      finished,orbconv,orbsetcnv
  integer      num_iter,iset,iorb
  integer ntol
  real(dp)      energy,eprev,eprv_sc,eprv_set,eprv_orb
  real(dp)      etol,  relaxn,cumulx
  real(dp)      tokcal
  parameter  ( tokcal = 627.509469_dp )
  real(dp)     w( hdim, *), eig(*),v1(*),v2(*)
  logical     coefflock(*)
  integer   int_out(2)    !  lengths not
  real(dp)   dbl_out(2)    !   protected

  cumulx  =  0.0_dp
  !     enter cycle of updating the wfn

  num_iter   = 0

  !     tolerance feathering - effective only if there is more than
  !     one orbital optimization group

  do   ntol  = ntol_e_min, ntol_e_max
     etol       = 10.0_dp**( -ntol )

     !     global optimization loop

     finished  = .false.
     do while ( .not. finished )
        eprv_sc  = energy

        !     orbital optimization loop

        orbconv   = .false.
        do while ( .not. orbconv )
           eprv_orb  = energy

           !     orbital subset loop

           do iset = 1, nset

              !     orbital subset optimization loop

              orbsetcnv   = .false.
              do while ( .not. orbsetcnv )
                 eprv_set  = energy
                 num_iter   = num_iter + 1

                 !     list of orbitals to be optimized together in this set

                 do iorb = orbset( 1, iset ), orbset( 2, iset )

                    if  (  dem_gs  )  then

                       !     zeroth-order optimization method:
                       !     'direct energy minimization' (ground states only)

                       call  demgs_opt( iorb,num_iter,cumulx,energy,etol,coefflock )

                    else

                       !     regular first-order derivative method

                       eprev  = energy
                       call  first_order_opt( iorb, w,eig,v1,v2, energy )

                       relaxn  =  ( energy - eprev )*tokcal
                       cumulx  =  cumulx  +  relaxn
                       int_out( 1 ) = num_iter
                       int_out( 2 ) = iorb
                       dbl_out( 1 ) = cumulx
                       dbl_out( 2 ) = relaxn/etol
                       call  xm_print( 'iters', ' ', int_out, dbl_out )
                       call xm_output( 'save', energy,etol )

                    end  if    !  DEM or 1st-order optimization method

                 end  do    !  orbital list for this subset

                 orbsetcnv = abs( energy - eprv_set )*tokcal .lt. etol
                 if ( num_iter .ge. max_iter ) orbsetcnv = .true.
                 if ( orbset( 1, iset ) - orbset( 2, iset ) .eq. 0  .and.   &
                      dem_gs  ) orbsetcnv = .true.
              end do                        !  orbital subset optimization loop

           end do                        !  list of orbital subsets

           orbconv = abs( energy - eprv_orb )*tokcal .lt. etol
           if ( num_iter .ge. max_iter ) orbconv = .true.
           if ( nset .eq. 1 )  then
              orbconv = .true.   ;   eprv_orb  = energy
           end if
        end do                        !  orbital optimization loop

        !     spin coupling optimization

        if ( nspinc .gt. 1 ) then
           call xm_print( 'comment', 'spin optimization;' )
           num_iter   = num_iter + 1

           eprev  = energy
           call  spin_opt( w,eig,v1,v2, energy )

           relaxn  =  ( energy - eprev )*tokcal
           cumulx  =  cumulx  +  relaxn
           int_out( 1 ) = num_iter
           int_out( 2 ) = 0
           dbl_out( 1 ) = cumulx
           dbl_out( 2 ) = relaxn/etol
           call  xm_print( 'iters', ' ', int_out, dbl_out )
           call xm_output( 'save', energy,etol )

           finished = abs( energy - eprv_sc )*tokcal .lt. etol
        else
           finished = .true.
        end if    !  spin opt

        if ( num_iter .ge. max_iter ) finished = .true.
     end do                        !  global optimization loop

  end do                        !  tolerance feathering


  !     print final result

  eprev  =  eprv_orb
  if  (  nspinc  .gt.  1  )  eprev  =  eprv_sc
  if  (  abs( energy - eprev )*tokcal .lt. etol  )  then

     call xm_print( 'title', 'calculation converged;' )
     dbl_out( 1 ) = energy
     call xm_print('parameter', 'total energy;', int_out, dbl_out )
     call xm_output( 'done', energy,etol )

  else  if  (  num_iter .ge. max_iter  )  then
     call xm_print( 'title', 'reached iteration limit;' )
  end if

end subroutine minimize_energy

!>     compute overlap integral <io(1) | jo(1)>, and put
!!     the result in sint
!!     \param io, jo [in]: index of two OLCAO in the integral
#ifdef SIMINT_INT
subroutine  ovint ( io, jo )
  use tools, only: dp
  use         integrals
  use xm
  use valence_simint, only : shell_map, integrals_store
  use SimintFortran
  implicit    none
  integer     io, jo, indf, i,j,k, n
  integer     mini,maxi,minj,maxj, mnshj,mxshj,jmax,ij
  integer     loxi,loxj, iao,jao, mnshi,mxshi
  integer     ic,ia,ii,lit,it,ig,li, jc,ja,jj,ljt,jt,jg,lj
  integer     funmin(7),funmax(7)
  data  funmin/ 1, 2, 5,11,21,36,57/
  data  funmax/ 1, 4,10,20,35,56,84/
  real(dp)    xi,yi,zi,xj,yj,zj
  real(dp)    zero
  parameter ( zero=0.0d+00 )

  integer simint_ret,myrank,master,nproc

  !     initialize weights of a contiguous basis list
  do i  = 1, max_obs
     coeffi( i ) = zero
     coeffj( i ) = zero
  end do

  !     scatter weights into the contiguous basis lists

  do i  = map_orbs( io ), map_orbs( io + 1 ) - 1
     if  ( xpset( i )  .lt.  1  )  then
        indf  =  xpset( i )  +  norbs
        call  ndf2obs( io, indf )
        k  =  0
        do j  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
           k  =  k  +  1
           coeffi( xpnew( k ) ) =   &
                coeffi( xpnew( k ) ) + coeff( j ) * coeff( i )
        end do
     else
        coeffi( xpset( i ) )  =  coeffi( xpset( i ) ) + coeff( i )
     end if
  end do

  do i  = map_orbs( jo ), map_orbs( jo + 1 ) - 1
     if  ( xpset( i )  .lt.  1  )  then
        indf  =  xpset( i )  +  norbs
        call  ndf2obs( jo, indf )
        k  =  0
        do j  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
           k  =  k  +  1
           coeffj( xpnew( k ) ) =   &
                coeffj( xpnew( k ) ) + coeff( j ) * coeff( i )
        end do
     else
        coeffj( xpset( i ) )  =  coeffj( xpset( i ) ) + coeff( i )
     end if
  end do

  sint = zero
!  print *,  "check!", orbas_atnum( io ),orbas_atnum( jo )
  loxi = 0
  do  ic = 1, orbas_atnum( io )
     ia  =   orbas_atset( ic, io )
     it  =  atom_t( ia )
     mnshi  =         map_atom2shell( it )
     mxshi  = mnshi + num_shell_atom( it ) - 1
     do  ii = mnshi, mxshi
        lit   =  ang_mom( ii ) + 1
        mini  =  funmin( lit )
        maxi  =  funmax( lit )
        loxj = 0

        do  jc = 1, orbas_atnum( jo )
           ja  =   orbas_atset( jc, jo )
           jt = atom_t( ja )
           mnshj =         map_atom2shell( jt )
           mxshj = mnshj + num_shell_atom( jt ) - 1

           do  jj = mnshj, mxshj
              ljt    = ang_mom( jj )+1
              minj   = funmin( ljt )
              maxj   = funmax( ljt )
              iao = loxi
              n = 0
              do i = mini, maxi
                 iao = iao + 1
                 jao = loxj
                 do j = minj, maxj
                    n = n + 1
                    jao = jao + 1
                    dij(n) = angn(i) *angn(j)*coeffi( iao )*coeffj( jao )
 !                  print *, iao,coeffi(iao), coeffj(jao)
                 end do
              end do
        ! integrals should contain an array of all the cartesian BF pairs in the
        ! given shells

                 simint_ret = simint_compute_overlap( shell_map(ii-mnshi+1,ia ), shell_map(jj-mnshj+1,ja ), integrals_store )


              ! check each BF here
              ij = 0
              do i = mini, maxi
                 do j = minj, maxj
                    ij = ij + 1
!                    if ( abs( dij( i ) ) .gt. dtol ) then
                       sint = sint + dij(ij)*integrals_store(ij)
!                       print *, "sint", sint, dij(ij),integrals_store(ij)
!                       print *, "sint", i,j
!                    endif
                 enddo
              enddo
!              print *,"cc",ia, ja,ii-mnshi+1,jj-mnshj+1
!              print *, "sint", sint, myrank
              loxj = loxj + shell_size( ang_mom( jj ) )
           end do          !  shells, jj
        end do          !  atoms, ja
           loxi = loxi + shell_size( ang_mom( ii ) )
 !       print *, "simint loxi",loxi
     end do          !  shells, ii
  end do          !  atoms, ia
end
#endif


!>     compute 1-electron core-hamiltonian integrals
!!
!!     this routine assumes the primitive weights come
!!     normalized both to themselves and the CGTO
!!     \param io, jo [in]: index of two OLCAO in the integral
#ifdef SIMINT_INT
subroutine  int1e ( io, jo )
  use tools, only: dp
  use         integrals
  use valence_simint, only : shell_map, integrals_store
  use SimintFortran
  implicit    none
  integer     io, jo, indf, i,j,k, ictr,nn 
  integer     mini,maxi,minj,maxj, mnshj,mxshj,jmax,ij
  integer     loxi,loxj, iao,jao, mnshi,mxshi
  integer     ic,ia,ii,lit,it,ig,li, jc,ja,jj,ljt,jt,jg,lj
  integer     funmin(7),funmax(7) 
  data  funmin/ 1, 2, 5,11,21,36,57/
  data  funmax/ 1, 4,10,20,35,56,84/
  real(dp)    xi,yi,zi,xj,yj,zj, nuchrg
  real(dp)    zero
  parameter ( zero=0.0d+00 )

  integer simint_ret,myrank,master,nproc


  !     initialize weights of a contiguous basis list

  do i  = 1, max_obs
     coeffi( i ) = zero
     coeffj( i ) = zero
  end do

  !     scatter weights into the contiguous basis lists

  ! coeff holds the coefficient of a basis function for a
  ! given basis function index in the orbital list. this should have been
  ! normalized by "normal()" previous to calling this routine.
  !
  do i  = map_orbs( io ), map_orbs( io + 1 ) - 1
     if  ( xpset( i )  .lt.  1  )  then
        indf  =  xpset( i )  +  norbs
        call  ndf2obs( io, indf )
        k  =  0
        do j  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
           k  =  k  +  1
           coeffi( xpnew( k ) ) =   &
                coeffi( xpnew( k ) ) + coeff( j ) * coeff( i )
        end do
     else
        coeffi( xpset( i ) )  =  coeffi( xpset( i ) ) + coeff( i )
     end if
  end do

  do i  = map_orbs( jo ), map_orbs( jo + 1 ) - 1
     if  ( xpset( i )  .lt.  1  )  then
        indf  =  xpset( i )  +  norbs
        call  ndf2obs( jo, indf )
        k  =  0
        do j  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
           k  =  k  +  1
           coeffj( xpnew( k ) ) =   &
                coeffj( xpnew( k ) ) + coeff( j ) * coeff( i )
        end do
     else
        coeffj( xpset( i ) )  =  coeffj( xpset( i ) ) + coeff( i )
     end if
  end do

  hint = zero

  loxi = 0
  do  ic = 1, orbas_atnum( io )
     ia  =   orbas_atset( ic, io )
     xi  =  coords( 1, ia )
     yi  =  coords( 2, ia )
     zi  =  coords( 3, ia )
     it  =  atom_t( ia )
     mnshi  =         map_atom2shell( it )
     mxshi  = mnshi + num_shell_atom( it ) - 1
     do  ii = mnshi, mxshi
        lit   =  ang_mom( ii )
        lit   =  lit + 1
        mini  =  funmin( lit )
        maxi  =  funmax( lit )

        loxj = 0
        do  jc = 1, orbas_atnum( jo )
           ja  =   orbas_atset( jc, jo )
           xj = coords( 1, ja )
           yj = coords( 2, ja )
           zj = coords( 3, ja )
           jt = atom_t( ja )
           mnshj =         map_atom2shell( jt )
           mxshj = mnshj + num_shell_atom( jt ) - 1
           do  jj = mnshj, mxshj
              ljt    = ang_mom( jj )
              ljt    = ljt + 1
              minj   = funmin( ljt )
              maxj   = funmax( ljt )

              jmax = maxj

              iao = loxi
              nn = 0
              do i = mini, maxi
                 iao = iao + 1
                 jao = loxj
                 do j = minj, jmax
                    jao = jao + 1
                    nn = nn + 1
                    dij( nn ) = angn(i)*angn( j )*coeffj( jao )*coeffi( iao )
!                    print *, "sint", coeffj( jao ),coeffi( iao )
                 end do
              end do

              simint_ret = simint_compute_ke( shell_map(ii-mnshi+1,ia ), shell_map(jj-mnshj+1,ja ), integrals_store )

              ! check each BF here
              ij = 0
              do i = mini, maxi
                 do j = minj, maxj
                    ij = ij + 1
!                    if ( abs( dij( i ) ) .gt. dtol ) then
                       hint = hint + dij(ij)*integrals_store(ij)
!                       print *, "hintsim", hint, dij(ij),integrals_store(ij)
!                       print *, "sint", i,j
!                    endif
                 enddo
              enddo

              do ictr = 1, natom
                 nuchrg = nuc_charge( atom_t( ictr ) )
                 if ( abs( nuchrg ) .gt. 1.0d-12 ) then       !  variable tol ?
                    nproc = 1
                   simint_ret = simint_compute_potential( nproc, nuc_charge( atom_t( ictr ) ), &
                        coords(1,ictr), coords(2,ictr), coords(3,ictr), &
                        shell_map(ii-mnshi+1,ia ), shell_map(jj-mnshj+1,ja ), integrals_store )
                    ! check each BF here
                    ij = 0
                    do i = mini, maxi
                       do j = minj, maxj
                          ij = ij + 1
                          !                    if ( abs( dij( i ) ) .gt. dtol ) then
                          hint = hint + dij(ij)*integrals_store(ij)
!                          print *, "hintsim", hint, dij(ij),integrals_store(ij)
                          !                       print *, "sint", i,j
                          !                    endif
                       enddo
                    enddo
                 endif
              enddo


              loxj = loxj + shell_size( ang_mom( jj ) )
           end do          !  shells, jj
        end do          !  atoms, ja
        loxi = loxi + shell_size( ang_mom( ii ) )
     end do          !  shells, ii
  end do          !  atoms, ia
end subroutine
#endif

!>    computes the 2-electron integral <io(1)jo(1)|ko(2)lo(2)>.
!!
!!    the resulting integral is return in gint in the integrals module
!!    \param io,jo,ko,lo [in]:  OLCAOs to integrate over
#ifdef SIMINT_INT
subroutine  int2e ( io,jo,ko,lo )
  use tools, only: dp
  use valence_simint, only : shell_map,work,integrals_store
  use SimintFortran
  use         integrals
  implicit    none
  integer     io,jo,ko,lo, ncomputed
  real(dp)    sum

  integer     ia,ic,it,ii,mnshi,mxshi,ish_beg
  integer     ja,jc,jt,jj,mnshj,mxshj,jsh_beg
  integer     ka,kc,kt,kk,mnshk,mxshk,ksh_beg
  integer     la,lc,lt,ll,mnshl,mxshl,lsh_beg
  integer     i,j,k,indf
  real(dp)    zero,d1,den,faci
  integer    funmin(7), funmax(7),  ij,kl,kao,lao,iao,jao
  data funmin/ 1, 2, 5,11,21,36,57/
  data funmax/ 1, 4,10,20,35,56,84/
  integer    l,m,n, lit,ljt,lkt,llt
  integer    mini,maxi,minj,maxj,mink,maxk,minl,maxl
  type(c_simint_multi_shellpair), target :: s2p_msh, s1p_msh
  parameter ( zero = 0.0d+00 )


  !     initialize weights of a contiguous basis list

  do i  = 1, max_obs
     coeffi( i ) = zero
     coeffj( i ) = zero
     coeffk( i ) = zero
     coeffl( i ) = zero
  end do

  !     scatter weights into the contiguous basis lists
  ! this loops over the set of basis functions that the "io" orbital is
  ! expanded in. map_obs holds the index of the basis functions in the basis set
  ! set that "io" is expanded in.
  do i  = map_orbs( io ), map_orbs( io + 1 ) - 1
     if  ( xpset( i )  .lt.  1  )  then
        indf  =  xpset( i )  +  norbs

        !     map the NDF OBS to the current orbital OBS

        call  ndf2obs( io, indf )
        k  =  0
        do j  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
           k  =  k  +  1
           coeffi( xpnew( k ) ) =   &
                coeffi( xpnew( k ) ) + coeff( j ) * coeff( i )
        end do
     else
        coeffi( xpset( i ) )  =  coeffi( xpset( i ) ) + coeff( i )
     end if
  end do

  do i  = map_orbs( jo ), map_orbs( jo + 1 ) - 1
     if  ( xpset( i )  .lt.  1  )  then
        indf  =  xpset( i )  +  norbs
        call  ndf2obs( jo, indf )
        k  =  0
        do j  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
           k  =  k  +  1
           coeffj( xpnew( k ) ) =   &
                coeffj( xpnew( k ) ) + coeff( j ) * coeff( i )
        end do
     else
        coeffj( xpset( i ) )  =  coeffj( xpset( i ) ) + coeff( i )
     end if
  end do

  do i  = map_orbs( ko ), map_orbs( ko + 1 ) - 1
     if  ( xpset( i )  .lt.  1  )  then
        indf  =  xpset( i )  +  norbs
        call  ndf2obs( ko, indf )
        k  =  0
        do j  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
           k  =  k  +  1
           coeffk( xpnew( k ) ) =   &
                coeffk( xpnew( k ) ) + coeff( j ) * coeff( i )
        end do
     else
        coeffk( xpset( i ) )  =  coeffk( xpset( i ) ) + coeff( i )
     end if
  end do

  do i  = map_orbs( lo ), map_orbs( lo + 1 ) - 1
     if  ( xpset( i )  .lt.  1  )  then
        indf  =  xpset( i )  +  norbs
        call  ndf2obs( lo, indf )
        k  =  0
        do j  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
           k  =  k  +  1
           coeffl( xpnew( k ) ) =   &
                coeffl( xpnew( k ) ) + coeff( j ) * coeff( i )
        end do
     else
        coeffl( xpset( i ) )  =  coeffl( xpset( i ) ) + coeff( i )
     end if
  end do


  gint  = zero
  ish_beg = 1
  do ia = 1, orbas_atnum( io )
     ic    = orbas_atset( ia, io )
     it    = atom_t( ic )
     mnshi =         map_atom2shell( it )
     mxshi = mnshi + num_shell_atom( it ) - 1
     do ii = mnshi, mxshi

        !     screen on LCAO weights per shell

        sum  =  zero
        do i = ish_beg, ish_beg + shell_size( ang_mom( ii ) ) - 1
           sum  =  sum  +  coeffi( i )**2
        end   do
        if  (  sum .gt. dtol  )  then


           jsh_beg = 1
           do ja = 1, orbas_atnum( jo )
              jc    = orbas_atset( ja, jo )
              jt    = atom_t( jc )
              mnshj =         map_atom2shell( jt )
              mxshj = mnshj + num_shell_atom( jt ) - 1
              do jj = mnshj, mxshj

                 sum  =  zero
                 do i = jsh_beg, jsh_beg + shell_size( ang_mom( jj ) ) - 1
                    sum  =  sum  +  coeffj( i )**2
                 end   do
                 if  (  sum .gt. dtol  )  then
                    call simint_initialize_multi_shellpair(s2p_msh)

                    call simint_create_multi_shellpair(1, shell_map(ii-mnshi+1,ic ), &
                         1, shell_map(jj-mnshj+1,jc ), s2p_msh, 2)

                    ksh_beg = 1
                    do ka = 1, orbas_atnum( ko )
                       kc    = orbas_atset( ka, ko )
                       kt    = atom_t( kc )
                       mnshk =         map_atom2shell( kt )
                       mxshk = mnshk + num_shell_atom( kt ) - 1
                       do kk = mnshk, mxshk

                          sum  =  zero
                          do i = ksh_beg, ksh_beg + shell_size( ang_mom( kk ) ) - 1
                             sum  =  sum  +  coeffk( i )**2
                          end   do
                          if  (  sum .gt. dtol  )  then


                             lsh_beg = 1
                             do la = 1, orbas_atnum( lo )
                                lc    = orbas_atset( la, lo )
                                lt    = atom_t( lc )
                                mnshl =         map_atom2shell( lt )
                                mxshl = mnshl + num_shell_atom( lt ) - 1
                                do ll = mnshl, mxshl

                                   sum  =  zero
                                   do i = lsh_beg, lsh_beg + shell_size( ang_mom( ll ) ) - 1
                                      sum  =  sum  +  coeffl( i )**2
                                   end   do
                                   if  (  sum .gt. dtol  )  then

!                                      print *, "simshe", ii,jj,kk,ll
 !                                     print *, kk-mnshk+1,kc,ll-mnshl+1, lc
                                      call simint_initialize_multi_shellpair(s1p_msh)
                                      call simint_create_multi_shellpair(1, shell_map(kk-mnshk+1,kc ), &
                         1, shell_map(ll-mnshl+1,lc ), s1p_msh, 2)

                                      lkt  = ang_mom( kk ) + 1
                                      mink = funmin( lkt )
                                      maxk = funmax( lkt )

                                      llt  = ang_mom( ll ) + 1
                                      minl = funmin( llt )
                                      maxl = funmax( llt )

                                      lit  = ang_mom( ii ) + 1
                                      mini = funmin( lit )
                                      maxi = funmax( lit )

                                      ljt  = ang_mom( jj ) + 1
                                      minj = funmin( ljt )
                                      maxj = funmax( ljt )

                                      kl = 0
                                      kao = 0
                                      do k = mink, maxk
                                         kao = kao + 1
                                         faci = angn(k)
                                         lao = 0
                                         do l = minl, maxl
                                            lao = lao + 1
                                            kl = kl+1
                                            dkl(kl) = faci*angn(l)*coeffk( ksh_beg+kao-1 )*coeffl( lsh_beg+lao-1 )
                                         end do
                                      end do

                                      ij = 0
                                      iao = 0
                                      do i = mini, maxi
                                         iao = iao + 1
                                         faci = angn(i)
                                         jao = 0
                                         do j = minj, maxj
                                            jao = jao + 1
                                            ij = ij + 1
                                            dij(ij) = faci*angn(j)*coeffi( ish_beg+iao-1 )*coeffj( jsh_beg+jao-1 )
                                         end do
                                      end do

                                      ncomputed = simint_compute_eri(s2p_msh, s1p_msh, 0.0d0, work, &
                                           integrals_store)

                                      n = 0
                                      do i = 1,ij
                                         d1 = dij(i)
                                         do k = 1, kl
                                            den = d1*dkl( k )
                                            n = n +1
!                                            if ( abs( den ) .gt. dtol ) then
                                               gint = gint + dkl(k)*dij(i)*integrals_store( n  )
!                                               print *,  "sss", gint, dij(i)*dkl(k), integrals_store( n  )
!                                            end if
                                         end do
                                      end do

 !                                     print *,  "sgint", gint,ii,jj,kk,ll
                                      call simint_free_multi_shellpair(s1p_msh)

                                   end if   !  screen shell LCAO weights
                                   lsh_beg = lsh_beg + shell_size( ang_mom( ll ) )
                                end do   !  shell l
                             end do   !  atom  l

                          end if   !  screen shell LCAO weights
                          ksh_beg = ksh_beg + shell_size( ang_mom( kk ) )
                       end do   !  shell k
                    end do   !  atom  k
                    
                    call simint_free_multi_shellpair(s2p_msh)

                 end if   !  screen shell LCAO weights
                 jsh_beg = jsh_beg + shell_size( ang_mom( jj ) )
              end do   !  shell j
           end do   !  atom  j

        end if   !  screen shell LCAO weights
        ish_beg = ish_beg + shell_size( ang_mom( ii ) )
     end do   !  shell i
  end do   !  atom  i
end subroutine
#endif

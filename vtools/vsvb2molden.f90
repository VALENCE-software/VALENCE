! gdf 30 jun 2022,   VSVB-to-MOLDEN format converter 
!   adapted from vismmp (GAMESS format)  

!! Currently assuming that normalized expansion coefficients are input 
!! This is true if they came from a VSVB calculation
!! Might not be true for a ModelKit/guess wfn

!! Default is cartesian gaussian functions 
!! Note, if this is used with output from the public version of 
!! VALENCE, shell functions must be sorted from CCA order to MOLDEN 
!! order (conversion is from ERIC order, see 'mdord') 
!! MOLDEN order differs from GAMESS only for F shells (?)  
!!  where the ordering seems a bit unsystematic 
!! See,  https://www3.cmbi.umcn.nl/molden/molden_format.html  


      program      vismmp
      implicit     none

      integer      natom,ndocc,totlen,npair,nunpd,natom_t,ndf
      integer      num_sh,num_pr,xpmax,nspinc,nang,nset,nxorb,mxctr

!     molecular geometry in terms of unique atom 'type' (input)

      integer,  allocatable ::  atom_t( : )
      real(8),  allocatable ::  coords( : , : )

!     basis sets for the unique atom types (input)

      integer,  allocatable ::  map_atom2shell( : )
      integer,  allocatable ::  map_shell2prim( : )
      integer,  allocatable ::         ang_mom( : )
      real(8),  allocatable ::      nuc_charge( : )
      real(8),  allocatable ::        exponent( : )
      real(8),  allocatable ::       con_coeff( : )

!     wave function information (input)

      integer,  allocatable ::     obs_atnum( : )
      integer,  allocatable ::     obs_atset( : , : )
      integer,  allocatable ::  map_orbs( : )
      integer,  allocatable ::     xpset( : )
      integer,  allocatable ::     atom_obs( : )
      integer,  allocatable ::     atom_mbs( : )
      integer,  allocatable ::     atom_ndf( : )
      real(8),  allocatable ::     coeff( : )


      integer     norbs,nbasis,nshell
      integer     i,j,k,n,m,ns,np,it,len,loff,con_length,ierr
      integer     io,ib,ia,iat,ish, mpt,indf,jb
      real(8),  allocatable ::  wavfn( : , : )
      real(8),  allocatable ::  pnorm( : )

      character(2)  atom(54)
      data          atom/' H','He','Li','Be',' B',' C',' N',   &
         ' O',' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',' K',  &
       'Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga',   &
       'Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc',   &
       'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe' /
      character(1)  shell_type, ang_type(0:4)
      data ang_type/'S','P','D','F','G'/
      integer     shell_size, shell_count, prim_count
      integer     shell_start(0:4), shell_func, istart,ifinsh
      data        shell_start/1,2,5,11,21/
      real(8)     energy,zero,dum
      parameter ( zero = 0.0d+00 )
      integer     ang_ptr(0:4)
      data ang_ptr/ 0, 1, 4,10,20/
!     integer     gmsord(35) 
!     data gmsord/ 1, 1,2,3, 1,3,6,2,4,5, 1,4,10,2,5,3,7,8,9,6  &
!     ,  1,5,15,2,6,4,9,13,14,3,10,12,7,8,11/
      integer     mdord(35) 
!       from the old VLNC ordering, 
!     data  mdord/ 1, 1,2,3, 1,3,6,2,4,5, 1,4,10,3,2,5,8,9,7,6  &
!     ,  1,5,15,2,6,4,9,13,14,3,10,12,7,8,11/
!       from the LibInt (CCA) ordering,      !! WARNING: F,G shells not done 
      data  mdord/ 1, 1,2,3, 1,4,6,2,3,5, 1,4,10,3,2,5,8,9,7,6  &
      ,  1,5,15,2,6,4,9,13,14,3,10,12,7,8,11/
      real(8)     tmp(15)


!     input the array dims

      read *,  natom,natom_t, npair,nunpd,ndocc, totlen,xpmax,  &
               nspinc, num_sh,num_pr,nang, ndf,nset,nxorb,mxctr

      if (   natom .lt. 1 ) stop '! atom count is 0 !'
      if ( natom_t .lt. 1 ) stop '! atom type count is 0 !'
      if (  num_sh .lt. 1 ) stop '! shell count is 0: please run basdims !'
      if (  num_pr .lt. 1 ) stop '! primitive count is 0: please run basdims !'
      if (  totlen .lt. 1 .or. xpmax .lt. 1 ) stop   &
      '! no storage for orbitals? Please run orbdims!'
      norbs  =  2*npair  +  nunpd  +  ndocc  +  ndf
      if (   norbs .lt. 1 ) stop '! orbital count is 0 !'

!     dummy read of run params 

      read *,  i, i, i, i, i, i,  dum, dum 


!     input the geometry

      allocate( atom_t(    natom ), stat = ierr )
      allocate( coords( 3, natom ), stat = ierr )

      do i = 1, natom
      read *,  atom_t( i ), ( coords( j, i ), j = 1, 3 )
      end  do


!     input the basis set

      allocate( map_atom2shell( natom_t + 1 ),  stat = ierr )
      allocate( map_shell2prim( num_sh  + 1 ),  stat = ierr )
      allocate(        ang_mom( num_sh      ),  stat = ierr )
      allocate(     nuc_charge( natom_t     ),  stat = ierr )
      allocate(       exponent( num_pr      ),  stat = ierr )
      allocate(      con_coeff( num_pr      ),  stat = ierr )

      ns = 1   ;   np = 1
      do    i  =  1,  natom_t
      map_atom2shell( i ) = ns
      read *,  nuc_charge( i ),  nshell
      do    j  =  1,  nshell
      map_shell2prim( ns ) = np
      read *,  ang_mom( ns ), con_length
      if ( con_length .eq. 1 ) then
      read *,  exponent( np )
      con_coeff( np ) = 1.0d+00
      np = np + 1
      else
      do    k  =  1,  con_length
      read *,  exponent( np ), con_coeff( np )
      np = np + 1
      end   do
      end if
      ns = ns + 1
      end   do
      end   do
      map_atom2shell( natom_t + 1 ) = ns
      map_shell2prim( num_sh  + 1 ) = np


!     output the geometry and basis set in MOLDEN format

      write ( *, * ) '[Molden Format] '
      write ( *, * ) '[Atoms] Angs '
      do    i  =  1,  natom 
      write ( *, 1 )  atom( int( nuc_charge( atom_t( i ) ) ) ), i  &
      , int( nuc_charge( atom_t( i ) ) ), ( coords( j, i ), j = 1, 3 )
      end   do

      write ( *, * ) '[GTO] ' 
      do    i  =  1,  natom 
      it = atom_t( i )
      write ( *, 2 )  i 
      do    j  =  map_atom2shell( it ), map_atom2shell( it+1 ) - 1
      write ( *, 3 )  ang_type( ang_mom( j ) ),  map_shell2prim( j+1 ) -map_shell2prim( j ) 
      do    k  =  map_shell2prim( j ), map_shell2prim( j+1 ) - 1
      write ( *, 4 )  exponent( k ),  con_coeff( k )
      end   do
      end   do
      write ( *, * )
      end   do



!     input the orbitals

      allocate( obs_atnum( norbs ),  stat = ierr )
      allocate( obs_atset( mxctr, norbs ),  stat = ierr )
      allocate( map_orbs( norbs + 1 ),  stat = ierr )
      allocate( xpset( totlen ),  stat = ierr )
      allocate( coeff( totlen ),  stat = ierr )


!     input the spin coupled (paired) orbitals

      j  =  1
      do     i = 1, 2*npair
      read *,   obs_atnum( i ),   &
      ( obs_atset( k, i ), k = 1, obs_atnum( i ) ),   n
      map_orbs( i )  =  j
      read *, ( xpset( k ), coeff( k ), k = j, j + n - 1 )
      j  =  j  +  n
      end   do


!     input unpaired orbitals

      do     i = 2*npair + 1, 2*npair + nunpd
      read *,   obs_atnum( i ),   &
      ( obs_atset( k, i ), k = 1, obs_atnum( i ) ),   n
      map_orbs( i )  =  j
      read *, ( xpset( k ), coeff( k ), k = j, j + n - 1 )
      j  =  j  +  n
      end   do


!     input doubly-occupied orbitals

      do     i = 2*npair + nunpd + 1, 2*npair + nunpd + ndocc
      read *,   obs_atnum( i ),   &
      ( obs_atset( k, i ), k = 1, obs_atnum( i ) ),   n
      map_orbs( i )  =  j
      read *, ( xpset( k ), coeff( k ), k = j, j + n - 1 )
      j  =  j  +  n
      end   do


!     input NDF's

      do     i = 2*npair + nunpd + ndocc + 1, 2*npair + nunpd + ndocc + ndf
      read *,   obs_atnum( i ),   &
      ( obs_atset( k, i ), k = 1, obs_atnum( i ) ),   n
      map_orbs( i )  =  j
      read *, ( xpset( k ), coeff( k ), k = j, j + n - 1 )
      j  =  j  +  n
      end   do
      norbs  =  2*npair  +  nunpd  +  ndocc  +  ndf
      map_orbs( norbs  +  1 )  =  j



!     get basis set size 

      nbasis = 0
      do  i = 1,  natom
      it = atom_t( i )
      do  j = map_atom2shell( it ), map_atom2shell( it + 1 ) - 1
      nbasis = nbasis + shell_size( ang_mom( j ) )
      end do
      end do

!     output the orbitals in MOLDEN format...
!     zero the LCAO-weight matrix for the wavefunction, preventing unsightly garbage

!!  for MOLDEN, orbitals could be processed one-at-a-time, avoiding n^2 storage 

      allocate( wavfn( nbasis, nbasis ),  stat = ierr )

      do i = 1, nbasis
      do j = 1, nbasis
      wavfn( j, i ) = zero
      end do
      end do


!     copy the LCAO weights into this matrix

      allocate( atom_mbs( natom + 1 ),  stat = ierr )
      allocate( atom_obs( mxctr + 1 ),  stat = ierr )
      allocate( atom_ndf( mxctr + 1 ),  stat = ierr )

!     first, store atom addresses within the MBS

      ib = 1
      do  iat = 1,  natom
      atom_mbs( iat )  =  ib
      it =  atom_t( iat )
      do  ish = map_atom2shell( it ), map_atom2shell( it+1 ) - 1
      ib = ib + shell_size( ang_mom( ish ) )
      end do  ;  end do
      atom_mbs( natom + 1 )  =  ib


      do io  = 1, norbs - ndf

!     store atom addresses within this OBS

      ib = 1
      do  iat = 1,  obs_atnum( io )
      atom_obs( iat )  =  ib
      it  =  atom_t( obs_atset( iat, io ) )
      do  ish = map_atom2shell( it ), map_atom2shell( it+1 ) - 1
      ib = ib + shell_size( ang_mom( ish ) )
      end do  ;  end do
      atom_obs( obs_atnum( io ) + 1 )  =  ib


!     loop over the weights

      do ib  =  map_orbs( io ), map_orbs( io+1 ) - 1
      if  (  xpset( ib ) .gt.  0   )  then

!     loop OBS to find which atom this label belongs to 
!     (this does not assume strict ordering of labels)

      do ia = 1, obs_atnum( io )
      if  (  xpset( ib ) .ge. atom_obs( ia )  .and.   &
             xpset( ib ) .lt. atom_obs( ia + 1 )  )   &
      mpt  =  atom_mbs( obs_atset( ia, io ) ) - atom_obs( ia )
      end  do
      wavfn( mpt + xpset( ib ), io ) = coeff( ib )

      else      !     process an NDF

      indf  =  xpset( ib )  +  norbs
      jb = 1
      do  iat = 1,  obs_atnum( indf )
      atom_ndf( iat )  =  jb
      it  =  atom_t( obs_atset( iat, indf ) )
      do  ish = map_atom2shell( it ), map_atom2shell( it+1 ) - 1
      jb = jb + shell_size( ang_mom( ish ) )
      end do  ;  end do
      atom_ndf( obs_atnum( indf ) + 1 )  =  jb

      do k  = map_orbs( indf ), map_orbs( indf + 1 ) - 1
      do ia = 1, obs_atnum( indf )
      if  (  xpset( k ) .ge. atom_ndf( ia )  .and.   &
             xpset( k ) .lt. atom_ndf( ia + 1 )  )   &
      mpt  =  atom_mbs( obs_atset( ia, indf ) ) - atom_ndf( ia )
      end  do
      wavfn( mpt + xpset( k ), io )  =   &
      wavfn( mpt + xpset( k ), io )  +  coeff( k ) * coeff( ib )
      end do

      end if      !   regular AO or NDF
      end do      !   AOs/NDFs
      end do      !   orbitals

      deallocate( atom_ndf,  stat = ierr )
      deallocate( atom_obs,  stat = ierr )
      deallocate( atom_mbs,  stat = ierr )


!     re-order shell functions from VALENCE to MOLDEN

      do   i = 1, norbs - ndf
      m = 0
      do   j = 1, natom
      it = atom_t( j )
      do   k  =  map_atom2shell( it ), map_atom2shell( it+1 ) - 1
      len = shell_size( ang_mom( k ) )
      loff  =  ang_ptr( ang_mom( k ) )
      do   n = 1, len
      tmp( n ) = wavfn( m + n, i )
      end do
      do   n = 1, len
      wavfn( m + n, i ) = tmp( mdord( n + loff ) )
      end do
      m = m + len
      end do   ;   end do   ;   end do


!     print in MOLDEN format

      write ( *, * ) ' [MO] ' 
      do   i = 1, norbs - ndf
      write ( *, * ) 'Sym=        a'  
      write ( *, * ) 'Ene=  -1.0000000000'   !< hand insert Eground-Ehole_state for now 
      write ( *, * ) 'Spin=       Alpha'  
      write ( *, * ) 'Occup=       2.0000000000'  
      do   j = 1, nbasis 
      write ( *, 6 ) j, wavfn( j, i ) 
      end   do
      end   do


    1 format(1a3,1x,1i4,1x,1i3,3f18.10)  
    2 format(1i4)               !  website says the trailing '  0' are deprecated 
    3 format(1a1,1x,1i4)    
    4 format(2f20.10)  
    6 format(1i6,1x,1f16.12)


      deallocate( wavfn,  stat = ierr )
      deallocate( coeff,  stat = ierr )
      deallocate( xpset,  stat = ierr )
      deallocate( map_orbs,  stat = ierr )
      deallocate( obs_atset,  stat = ierr )
      deallocate( obs_atnum,  stat = ierr )
      deallocate(      con_coeff,  stat = ierr )
      deallocate(       exponent,  stat = ierr )
      deallocate(     nuc_charge,  stat = ierr )
      deallocate(        ang_mom,  stat = ierr )
      deallocate( map_shell2prim,  stat = ierr )
      deallocate( map_atom2shell,  stat = ierr )
      deallocate( coords, stat = ierr )
      deallocate( atom_t, stat = ierr )

      end  program  vismmp



      integer function shell_size( l_mom )
      implicit    none
      integer     l_mom
      shell_size  =  ( l_mom + 1 )*( l_mom + 2 )/2
      end




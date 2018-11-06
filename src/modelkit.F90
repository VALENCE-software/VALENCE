! gdf 14 aug 2018,  modelkit

!     Assemble a 'model kit' wave function.
!     Current definitions, assumptions and simplifications:
!     To make sense of input LCAO weights, guesses are assumed 
!     to be in the following 'standard' orientation: the bond 
!     axis is the z-axis ('sigma'), the first atom is at the 
!     origin, the second atom lies along the +ve z-axis; the bond 
!     plane is the z,y-plane; the perpendicular direction is x.
!     sigma lone pairs point in the positive z-direction.

!     A core orbital is a 1-atom orbital that does not require 
!     reorientation (core shells are assumed complete and filled 
!     (spherical) so their orientation is irrelevant).
!     A lone pair orbital is a 1-atom orbital that requires 
!     reorientation, involving either p or s,p shells.
!     A bond orbital involves 2 atoms (but allow for 
!     flexibility to increase this).
!     A pi-bond or pi-lone pair involves only p functions.
!     Inevitably, some of the reorientation policies are 
!     'rules of thumb' - but this is all a guess anyway.
!     Add special case 3-center bonds as required.

!     for compatibility with other tools it is assumed that this 
!     tool makes ground state wave functions of DOCC orbitals


      program     modelkit
      implicit    none
      integer     i,j,k, shell_size, bat1,bat2,batm,batx
      integer     it, jb
      integer     nshell,nbasis,ierr,ns
      real(8)     zw,vx,vy,vz, a1,a2,a3,b1,b2,b3,p1,p2,p3,dum
      real(8)     zero,           two
      parameter ( zero = 0.0d+00, two = 2.0d+00 )
      logical     linear, pix,piy,sigma


!     Array dimensions impractical to set by input
!     max_c  =  max. (atomic) centers per bond 
!     max_v  =  max. atomic valency

      integer      max_c,     max_v 
      parameter  ( max_c= 2,  max_v= 20 )


!     wave function and basis set

      integer      natom,ndocc,totlen,npair,nunpd,natom_t
      integer      num_sh,num_pr,xpmax,nspinc,nang,ndf,nset,nxorb,mxctr
      integer      con_length,orbas_atset( max_c ),xplen

      integer,  allocatable ::  map_atom2shell( : )
      integer,  allocatable ::         ang_mom( : )
      integer,  allocatable ::     xpset( : )
      real(8),  allocatable ::     coeff( : )


!     molecular geometry 

      integer     bondat( max_v )
      integer,  allocatable ::  atom_t( : )
      real(8),  allocatable ::  coords( : , : )


!     model-kit 'parts'
!     core orbitals

      integer    ngcore, ncores,  gpmax
      integer,  allocatable ::  gcore_atomt(     : )
      integer,  allocatable ::  gcore_xplen(     : )
      integer,  allocatable ::  gcore_xpset( : , : )
      real(8),  allocatable ::  gcore_coeff( : , : )

!     bond orbitals

      integer    ngbond
      integer,  allocatable ::  gbond_atnum(     : )
      integer,  allocatable ::  gbond_atset( : , : )
      integer,  allocatable ::  gbond_xplen(     : )
      integer,  allocatable ::  gbond_xpset( : , : )
      real(8),  allocatable ::  gbond_coeff( : , : )

!     lone pair orbitals

      integer    nglone
      integer,  allocatable ::  glone_atomt(     : )
      integer,  allocatable ::  glone_xplen(     : )
      integer,  allocatable ::  glone_xpset( : , : )
      real(8),  allocatable ::  glone_coeff( : , : )

!     define bonding orbitals in target molecule

      integer     nbond_orbs, nbond_at
      integer,  allocatable ::   bond_atoms( : , : )
      integer,  allocatable ::    bond_type(     : )

!     define lone pair orbitals in target molecule

      integer     nlone_pair, lpa
      integer,  allocatable ::     lone_atom( : )
      integer,  allocatable ::     lone_type( : )



!     input the array dims

      read *,  natom,natom_t, npair,nunpd,ndocc, totlen,xpmax,  &
               nspinc, num_sh,num_pr,nang, ndf,nset,nxorb,mxctr

      if (   natom .lt. 1 ) stop '! atom count is 0 !'
      if ( natom_t .lt. 1 ) stop '! atom type count is 0 !'
      if (  num_sh .lt. 1 ) stop '! shell count is 0 !'



!     input the geometry

      allocate( atom_t(    natom ), stat = ierr )
      allocate( coords( 3, natom ), stat = ierr )

      do    i  =  1,  natom 
      read *,  atom_t( i ), ( coords( j, i ), j = 1, 3 )
      end   do

      if ( natom .le. 4 ) then
      do    i  =  1,  natom 
      bondat( i ) = i
      end   do
      end  if




!     input the unique atomic basis sets
!     for compatibility, the gaussian primitive data and 
!     nuclear charges are read but not stored

      allocate( map_atom2shell( natom_t + 1 ),  stat = ierr )
      allocate(        ang_mom( num_sh      ),  stat = ierr )

      ns = 1
      do    i  =  1,  natom_t
      map_atom2shell( i ) = ns
      read *,  dum,  nshell
      do    j  =  1,  nshell
      read *,  ang_mom( ns ), con_length
      if ( con_length .eq. 1 ) then
      read *,  dum
      else
      do    k  =  1,  con_length
      read *,  dum, dum
      end   do
      end if
      ns = ns + 1
      end   do
      end   do
      map_atom2shell( natom_t + 1 ) = ns




!     estimate the largest expansion length in order to read the guesses

      gpmax  = 0
      do  it = 1,  natom_t
      nbasis = 0
      do  j = map_atom2shell( it ), map_atom2shell( it + 1 ) - 1
      nbasis = nbasis + shell_size( ang_mom( j ) )
      end do
      gpmax  =  max( gpmax, nbasis )
      end do 
      gpmax  =  gpmax * max_c





!     construct the guess wave function

      xpmax  = 0
      totlen  = 0
      mxctr  =  1

!     input the core orbitals

      ncores   = 0 
      read *,      ngcore
      if  (  ngcore  .gt.  0  )  then
      allocate( gcore_atomt(        ngcore ),  stat = ierr )
      allocate( gcore_xplen(        ngcore ),  stat = ierr )
      allocate( gcore_xpset( gpmax, ngcore ),  stat = ierr )
      allocate( gcore_coeff( gpmax, ngcore ),  stat = ierr )

      do    i  =  1,  ngcore
      read *,   gcore_atomt( i ), gcore_xplen( i )
      read *, ( gcore_xpset( j, i ), gcore_coeff( j, i ),  &
            j = 1, gcore_xplen( i ) )
      end   do



!     match atoms to their guesses (guesses reference atom TYPES)
!  a) populate cores - this is a straight copy
!     get the total core orbital count

      do    i  =  1,  natom 
      do    j  =  1,  ngcore
      if  ( gcore_atomt( j ) .eq. atom_t( i ))  then 
      ncores   =  ncores  +  1
      write (*,2)  1, i, gcore_xplen( j ) 
      write (*,3) ( gcore_xpset( k, j ), gcore_coeff( k, j )  &
      ,      k = 1, gcore_xplen( j ) )

      xpmax  =  max( xpmax, gcore_xplen( j ) )
      totlen  =  totlen  +  gcore_xplen( j )
      end   if 
      end   do 
      end   do

      deallocate( gcore_atomt,  stat = ierr )
      deallocate( gcore_xplen,  stat = ierr )
      deallocate( gcore_xpset,  stat = ierr )
      deallocate( gcore_coeff,  stat = ierr )

      end  if    !  any cores




!     workspace for orbital reorientations

      allocate( xpset( gpmax ),  stat = ierr )
      allocate( coeff( gpmax ),  stat = ierr )



!  b) populate bonds, aligning guesses with the target molecule
!     input the bonds 

      nbond_orbs  =  0
      read *,      ngbond
      if  (  ngbond  .gt.  0  )  then

      allocate( gbond_atnum(        ngbond ),  stat = ierr )
      allocate( gbond_atset( max_c, ngbond ),  stat = ierr )
      allocate( gbond_xplen(        ngbond ),  stat = ierr )
      allocate( gbond_xpset( gpmax, ngbond ),  stat = ierr )
      allocate( gbond_coeff( gpmax, ngbond ),  stat = ierr )

      do    i  =  1,  ngbond
      read *,   gbond_atnum( i ),  &
      ( gbond_atset( j, i ), j = 1, gbond_atnum( i ) ),   &
        gbond_xplen( i )
      read *, ( gbond_xpset( j, i ), gbond_coeff( j, i ),  &
            j = 1, gbond_xplen( i ) )
      mxctr  =  max(  mxctr,  gbond_atnum( i ) )
      end   do


!     input the bonds in the target molecule
!     the 'bond type' must match one of the guess types above

      read *,      nbond_orbs
      allocate( bond_atoms( 2, nbond_orbs ),  stat = ierr )
      allocate( bond_type(     nbond_orbs ),  stat = ierr )

      do    i  =  1,  nbond_orbs 
      read *, bond_atoms( 1, i ), bond_atoms( 2, i ), bond_type( i )
      end   do



      do    i   =  1,  nbond_orbs
      do    jb  =  1,  ngbond
      if (  jb .eq. bond_type( i ) )  then
      bat1 = bond_atoms( 1, i )
      bat2 = bond_atoms( 2, i )
      orbas_atset( 1 )  =  bat1
      orbas_atset( 2 )  =  bat2

      call  orbtype ( gbond_atnum( jb ), gbond_atset( 1, jb )  &
      ,               gbond_xpset( 1, jb ), map_atom2shell  &
      ,               ang_mom, pix,piy,sigma )

      if ( sigma ) then

!     input weight sense: bond vector points FROM bat1 TO bat2

      vx  =  coords( 1, bat2 ) - coords( 1, bat1 )
      vy  =  coords( 2, bat2 ) - coords( 2, bat1 )
      vz  =  coords( 3, bat2 ) - coords( 3, bat1 )

!     Pi-Bond types ...

      else if ( pix .or. piy ) then

      a1  =  coords( 1, bat2 ) -coords( 1, bat1 )
      a2  =  coords( 2, bat2 ) -coords( 2, bat1 )
      a3  =  coords( 3, bat2 ) -coords( 3, bat1 )

!     now find a 3rd non-colinear bonded atom to define a plane
!     try the first bond atom, then the second

      call  get_plane ( bat1, bat2, bond_atoms, nbond_orbs,  &
                        bondat, nbond_at, coords, p1,p2,p3 ) 
      batm = bat1 
      if ( nbond_at .eq. 1 ) then 
      call  get_plane ( bat2, bat1, bond_atoms, nbond_orbs,  &
                        bondat, nbond_at, coords, p1,p2,p3 ) 
      batm = bat2
      end  if


!     for Px, get perpendicular vector

      b1  =  coords( 1, batm ) - p1
      b2  =  coords( 2, batm ) - p2
      b3  =  coords( 3, batm ) - p3
      vx  =  a2*b3 - a3*b2
      vy  =  a3*b1 - a1*b3
      vz  =  a1*b2 - a2*b1

      if ( piy ) then    !  in-plane vector for Py
      b1  =  vx
      b2  =  vy
      b3  =  vz
      vx  =  a2*b3 - a3*b2
      vy  =  a3*b1 - a1*b3
      vz  =  a1*b2 - a2*b1
      end  if
      end  if     ! sigma or pi bond

      call  reorient( gbond_atnum( jb )  &
      ,               gbond_atset( 1, jb ), map_atom2shell, ang_mom   &
      ,               gbond_xpset( 1, jb ), gbond_coeff( 1, jb )   &
      ,               gbond_xplen( jb ), xplen, xpset, coeff, vx,vy,vz )

      write (*,2)  gbond_atnum( jb ),  &
      ( orbas_atset( j ), j = 1, gbond_atnum( jb ) ), xplen
      write (*,3) ( xpset( j ), coeff( j ), j = 1, xplen )

      xpmax  =  max( xpmax, xplen )
      totlen  =  totlen  +  xplen
      end   if    !  matching guess
      end   do    !  bond orbitals
      end   do    !  bonds in target molecule

      deallocate(  bond_type,   stat = ierr )
      deallocate( gbond_atnum,  stat = ierr )
      deallocate( gbond_atset,  stat = ierr )
      deallocate( gbond_xplen,  stat = ierr )
      deallocate( gbond_xpset,  stat = ierr )
      deallocate( gbond_coeff,  stat = ierr )

      end   if    !  any bonds




!  c) populate lone pairs (sometimes the lone pair atom is called the 'host')
!     input lone pairs 

      nlone_pair  =  0
      read *,      nglone
      if  (  nglone  .gt.  0  )  then
      if  (  ngbond  .lt.  1  )  stop '! no lone pairs without bonds !'

      allocate( glone_atomt(        nglone ),  stat = ierr )
      allocate( glone_xplen(        nglone ),  stat = ierr )
      allocate( glone_xpset( gpmax, nglone ),  stat = ierr )
      allocate( glone_coeff( gpmax, nglone ),  stat = ierr )

      do    i  =  1,  nglone
      read *,   glone_atomt( i ), glone_xplen( i )
      read *, ( glone_xpset( j, i ), glone_coeff( j, i ),  &
            j = 1, glone_xplen( i ) )
      end   do


!     input the lone pairs in the target molecule

      read *,      nlone_pair
      allocate( lone_atom( nlone_pair ),  stat = ierr )
      allocate( lone_type( nlone_pair ),  stat = ierr )

      do    i  =  1,  nlone_pair
      read *, lone_atom( i ), lone_type( i )
      end   do


      do    i   =  1,  nlone_pair
      do    jb  =  1,  nglone
      if (  jb .eq. lone_type( i ) )  then
      lpa = lone_atom( i )

      call  orbtype ( 1, glone_atomt( jb )  &
      ,               glone_xpset( 1, jb ), map_atom2shell  &
      ,               ang_mom, pix,piy,sigma )

      call  get_bondats ( lpa, bond_atoms, nbond_orbs,  &
                          bondat, nbond_at ) 


!     orientation cases ...

      if ( sigma ) then

!     a sigma-LP is always opposed and perpendicular to the 
!     point, axis or plane of just the bonded atoms

           if ( nbond_at .eq. 1 ) then

      vx  =  coords( 1, lpa ) -coords( 1, bondat( 1 ) )
      vy  =  coords( 2, lpa ) -coords( 2, bondat( 1 ) )
      vz  =  coords( 3, lpa ) -coords( 3, bondat( 1 ) )

      else if ( nbond_at .eq. 2 ) then

!     if the 3 atoms turn out to be colinear, assume this is 
!     user-error
!     this policy orients toward the host from a point along 
!     the line connecting the 2 bonded atoms, the point is 
!     weighted toward the bonded atom closest to the host

      bat1 = bondat( 1 )
      bat2 = bondat( 2 )
      vx  =  coords( 1, lpa ) -coords( 1, bat1 ) 
      vy  =  coords( 2, lpa ) -coords( 2, bat1 ) 
      vz  =  coords( 3, lpa ) -coords( 3, bat1 ) 
      a1  =  sqrt( vx**2 + vy**2 + vz**2 )
      vx  =  coords( 1, lpa ) -coords( 1, bat2 ) 
      vy  =  coords( 2, lpa ) -coords( 2, bat2 ) 
      vz  =  coords( 3, lpa ) -coords( 3, bat2 ) 
      a2  =  sqrt( vx**2 + vy**2 + vz**2 )
      if ( a1 .lt. a2 ) then 
      batm  =  bat1
      batx  =  bat2
      else
      batm = bat2
      batx = bat1
      end  if
      zw  =  min( a1, a2 )/( a1 + a2 )
      a1  =  coords( 1, batm )
      a2  =  coords( 2, batm )
      a3  =  coords( 3, batm ) 
      a1  =  a1 + zw*( coords( 1, batx ) -coords( 1, batm ) ) 
      a2  =  a2 + zw*( coords( 2, batx ) -coords( 2, batm ) ) 
      a3  =  a3 + zw*( coords( 3, batx ) -coords( 3, batm ) )
      vx  =  coords( 1, lpa ) - a1
      vy  =  coords( 2, lpa ) - a2
      vz  =  coords( 3, lpa ) - a3

      else if ( nbond_at .ge. 3 ) then

!     for >3 bonded atoms, could help the user by averaging 
!     the possible 'normal' vectors ...??

      a1  =  coords( 1, bondat( 1 ) ) -coords( 1, bondat( 2 ) )
      a2  =  coords( 2, bondat( 1 ) ) -coords( 2, bondat( 2 ) )
      a3  =  coords( 3, bondat( 1 ) ) -coords( 3, bondat( 2 ) )
      b1  =  coords( 1, bondat( 1 ) ) -coords( 1, bondat( 3 ) )
      b2  =  coords( 2, bondat( 1 ) ) -coords( 2, bondat( 3 ) )
      b3  =  coords( 3, bondat( 1 ) ) -coords( 3, bondat( 3 ) )
      vx  =  a2*b3 - a3*b2
      vy  =  a3*b1 - a1*b3
      vz  =  a1*b2 - a2*b1

!     for sigma, vector must point FROM bond atoms TO host
!     so invert direction if opposed

      a1  =  coords( 1, lpa ) -coords( 1, bondat( 1 ) )
      a2  =  coords( 2, lpa ) -coords( 2, bondat( 1 ) )
      a3  =  coords( 3, lpa ) -coords( 3, bondat( 1 ) )

      if (  a1*vx + a2*vy + a3*vz  .lt. zero ) then
      vx  =  -vx
      vy  =  -vy
      vz  =  -vz
      end  if

      end  if     ! no. bonded atoms

!     Pi-Lone Pair types ...

      else if ( pix .or. piy ) then

!     a pi-LP is always perpendicular to the plane of the 
!     (bonded+host) atoms

           if ( nbond_at .eq. 1 ) then

      batm = bondat( 1 )
      call  get_plane ( batm, lpa, bond_atoms, nbond_orbs,  &
                        bondat, nbond_at, coords, p1,p2,p3 ) 

!     orientation for Pi-x is perpendicular to the plane

      a1  =  coords( 1, lpa ) -coords( 1, batm )
      a2  =  coords( 2, lpa ) -coords( 2, batm )
      a3  =  coords( 3, lpa ) -coords( 3, batm )
      b1  =  coords( 1, lpa ) -p1
      b2  =  coords( 2, lpa ) -p2
      b3  =  coords( 3, lpa ) -p3
      vx  =  a2*b3 - a3*b2
      vy  =  a3*b1 - a1*b3
      vz  =  a1*b2 - a2*b1

!     y is in the plane (perp. to both axis and x)

      if ( piy ) then
      b1  =  vx
      b2  =  vy
      b3  =  vz
      vx  =  a2*b3 - a3*b2
      vy  =  a3*b1 - a1*b3
      vz  =  a1*b2 - a2*b1
      end  if

      else if ( nbond_at .ge. 2 ) then

      a1  =  coords( 1, lpa ) -coords( 1, bondat( 1 ) )
      a2  =  coords( 2, lpa ) -coords( 2, bondat( 1 ) )
      a3  =  coords( 3, lpa ) -coords( 3, bondat( 1 ) )
      p1  =  coords( 1, bondat( 2 ) )
      p2  =  coords( 2, bondat( 2 ) )
      p3  =  coords( 3, bondat( 2 ) )

!     nbond_at >2 cannot be colinear, so that situation is 
!     either planar or (user) error
!     but, for nbond_at =2, check for linearity

      if ( nbond_at .eq. 2 ) then
      bondat( 3 ) = lpa 
      if ( linear( coords, bondat, 3 ) ) then
      bat1 = bondat( 1 )
      bat2 = bondat( 2 )

!     seek 3rd, non-colinear, bonded atom to define a plane
!     try the first bond atom, then the other

      call  get_plane ( bat1, lpa, bond_atoms, nbond_orbs,  &
                        bondat, nbond_at, coords, p1,p2,p3 ) 
      if ( nbond_at .eq. 1 ) then
      call  get_plane ( bat2, lpa, bond_atoms, nbond_orbs,  &
                        bondat, nbond_at, coords, p1,p2,p3 ) 
      end  if
      end  if     !   linear
      end  if     !   2 bonded atoms

      b1  =  coords( 1, lpa ) -p1
      b2  =  coords( 2, lpa ) -p2
      b3  =  coords( 3, lpa ) -p3
      vx  =  a2*b3 - a3*b2
      vy  =  a3*b1 - a1*b3
      vz  =  a1*b2 - a2*b1 

      if ( piy ) then
      b1  =  vx
      b2  =  vy
      b3  =  vz
      vx  =  a2*b3 - a3*b2
      vy  =  a3*b1 - a1*b3
      vz  =  a1*b2 - a2*b1
      end  if

      end  if     !  no. bonded atoms
      end  if     !  lone pair type

      call  reorient( 1   &
      ,               glone_atomt( jb ), map_atom2shell, ang_mom   &
      ,               glone_xpset( 1, jb ), glone_coeff( 1, jb )   &
      ,               glone_xplen( jb ), xplen, xpset, coeff, vx,vy,vz )

      write (*,2)  1, lpa, xplen
      write (*,3) ( xpset( j ), coeff( j ), j = 1, xplen )

      xpmax  =  max( xpmax, xplen )
      totlen  =  totlen  +  xplen
      end   if    !  matching guess
      end   do    !  lone pair orbitals
      end   do    !  lone pairs in target molecule


      deallocate( lone_atom,    stat = ierr )
      deallocate( lone_type,    stat = ierr )
      deallocate( glone_atomt,  stat = ierr )
      deallocate( glone_xplen,  stat = ierr )
      deallocate( glone_xpset,  stat = ierr )
      deallocate( glone_coeff,  stat = ierr )

      end   if    !  any lone pairs


      deallocate( bond_atoms,   stat = ierr )
      deallocate( coeff,  stat = ierr )
      deallocate( xpset,  stat = ierr )
      deallocate( ang_mom,  stat = ierr )
      deallocate( map_atom2shell,  stat = ierr )
      deallocate( coords, stat = ierr )
      deallocate( atom_t, stat = ierr )




!     print the updated array dims

      ndocc  =  ncores + nbond_orbs + nlone_pair

      open ( unit = 10, file = 'run_params',   &
             form = 'formatted', status = 'unknown' )

      write ( 10, 4 ) natom,natom_t, npair,nunpd,ndocc, &
      totlen,xpmax, nspinc, num_sh,num_pr,nang, ndf,nset,nxorb,mxctr

      close ( unit = 10 )


    2 format(5x,1i5,5x,10i5)
    3 format(4(1i10,1x,1f13.8))
    4 format( 1i6,3i10,2i10,9i10 )

      end program modelkit




!     identify the orientation of the orbital 

      subroutine  orbtype ( atnum, atset, xpset, map_atom2shell  &
      ,                     ang_mom, ox,oy,oz )
      implicit    none
      integer     atnum, atset(*), xpset(*), map_atom2shell(*)
      integer     ang_mom(*)
      integer     i,j,iat,ish,iao,it, shell_size
      logical     ox,oy,oz


      ox = .false.      !  perpendicular to bond axis
      oy = .false.      !  perpendicular to bond axis and x
      oz = .false.      !  oriented along bond axis

!     align input expansion with atom basis sets to 
!     determine AO type

      i = 0
      j = 1
      do  iat = 1,  atnum
      it = atset( iat )
      do  ish = map_atom2shell( it ), map_atom2shell( it+1 ) - 1
      do  iao = 1, shell_size( ang_mom( ish ) )
      i = i + 1
      if ( i .eq. xpset( j ) ) then 

!     rules for determining the orientation based on AO type
!     s,p AOs so far, add other types here

      if ( ang_mom( ish ) .eq. 1 ) then        ! p
      if ( iao .eq. 1 ) ox = .true.
      if ( iao .eq. 2 ) oy = .true.
      if ( iao .eq. 3 ) oz = .true.
      else if ( ang_mom( ish ) .eq. 0 ) then   ! s
      oz = .true.
      end  if

      j = j + 1
      end  if
      end  do
      end  do
      end  do

!     sanity check: only individual ox,oy,oz can be true

      if ( ox .and. oy ) stop '! mixed x,y bond type ?'
      if ( ox .and. oz ) stop '! mixed x,z bond type ?'
      if ( oy .and. oz ) stop '! mixed y,z bond type ?'

!     but one of them must be true

      if ( .not.ox .and. .not.oy .and. .not.oz )  &
      stop '! no bond type ?'
      end 




!     list all the unique atoms bonded to a given atom

      subroutine  get_bondats ( atom, bonds, num_bonds,  &
                  bondat, nbond_at )
      implicit    none
      integer     atom, bonds( 2, *), num_bonds
      integer     bondat(*), nbond_at, i,j,tmp

      nbond_at = 0
      do i = 1, num_bonds
            if (  bonds( 1, i ) .eq. atom ) then
      nbond_at = nbond_at + 1
         bondat( nbond_at ) = bonds( 2, i )
      else  if (  bonds( 2, i ) .eq. atom ) then
      nbond_at = nbond_at + 1
         bondat( nbond_at ) = bonds( 1, i )
      end  if

!     ensure the bonded atoms listed are unique
!     ie. account for multiple bonds

      tmp    =     nbond_at
      do  j  =  1, nbond_at - 1
      if ( bondat( nbond_at ) .eq. bondat( j ) ) tmp = tmp - 1
      end  do
      nbond_at = tmp
      end  do
      end




!     test for 'linearity' amongst a set of atoms

      logical function linear ( coords, atoms, natom )
      implicit    none
      integer     natom, atoms(*), i
      real*8      coords( 3, * )
      real*8      cosavg, vx,vy,vz,a1,a2,a3,b1,b2,b3
      real*8      zero,           one,           tol
      parameter ( zero = 0.0d+00, one = 1.0d+00, tol = 1.0d-06 )


      if ( natom .lt. 3 ) stop '! too few atoms !'

      linear = .false.
      cosavg = zero
      a1  =       coords( 1, atoms( 1 ) )
      a2  =       coords( 2, atoms( 1 ) )
      a3  =       coords( 3, atoms( 1 ) ) 
      vx  =  a1 - coords( 1, atoms( 2 ) )
      vy  =  a2 - coords( 2, atoms( 2 ) )
      vz  =  a3 - coords( 3, atoms( 2 ) ) 
      do  i  = 3,  natom
      b1  =  a1 - coords( 1, atoms( i ) )
      b2  =  a2 - coords( 2, atoms( i ) )
      b3  =  a3 - coords( 3, atoms( i ) )
      cosavg = cosavg +   &
               ( b1*vx + b2*vy + b3*vz ) /  &
         ( sqrt( b1**2 + b2**2 + b3**2 ) *  &
           sqrt( vx**2 + vy**2 + vz**2 ) ) 
      end  do
      cosavg = cosavg / dble( natom - 2 )
      if ( abs( one - cosavg ) .lt. tol ) linear = .true.
      end




!     realign the input orbital to the target orientation
!     currently, only p-shells considered

      subroutine  reorient ( atnum, atset  &
      ,                      map_atom2shell, ang_mom, g_xpset, g_coeff  &
      ,                      g_xplen, xplen, xpset, coeff, vx,vy,vz )
      implicit    none
      integer     atnum, atset(*), map_atom2shell(*), ang_mom(*)
      integer     g_xpset(*), g_xplen, xplen, xpset(*)
      real*8      g_coeff(*), coeff(*)
      real*8      tol,  vx,vy,vz,zw
      parameter ( tol = 1.0d-06 )
      integer     iorbas, nf_orb, ib, iat,ish,iao,it
      integer     sh_ao_1, shell_size


      nf_orb = 0            ! counts functions of the reoriented orbital
      do  ib = 1, g_xplen   ! counts functions of the guess orbital
      if ( g_xpset( ib ) .lt. 1 ) then   ! NDF's not reoriented 
      nf_orb  =  nf_orb  +  1
      coeff( nf_orb )  =  g_coeff( ib )
      xpset( nf_orb )  =  g_xpset( ib )
      else                ! search for the function in the OBS
      iorbas = 0          ! counts functions in the OBS
      do  iat = 1,  atnum
      it =  atset( iat )
      do  ish = map_atom2shell( it ), map_atom2shell( it+1 ) - 1
      sh_ao_1 = iorbas + 1
      do  iao = 1, shell_size( ang_mom( ish ) )
      iorbas = iorbas + 1

      if ( iorbas .eq. g_xpset( ib ) ) then 
      if ( ang_mom( ish ) .eq. 1 )  then

      zw  =  g_coeff( ib )
      zw  =  zw*(vx**2 + vy**2 + vz**2)**( -0.5d+00 )

      if ( abs( zw*vx ) .gt. tol ) then
      nf_orb  =  nf_orb  +  1
      coeff( nf_orb )  =  zw*vx
      xpset( nf_orb )  =  sh_ao_1
      end  if 
      if ( abs( zw*vy ) .gt. tol ) then
      nf_orb  =  nf_orb  +  1
      coeff( nf_orb )  =  zw*vy
      xpset( nf_orb )  =  sh_ao_1 + 1
      end  if 
      if ( abs( zw*vz ) .gt. tol ) then
      nf_orb  =  nf_orb  +  1
      coeff( nf_orb )  =  zw*vz
      xpset( nf_orb )  =  sh_ao_1 + 2
      end  if

      else           !  weight is just copied 
      nf_orb  =  nf_orb  +  1
      coeff( nf_orb )  =  g_coeff( ib )
      xpset( nf_orb )  =  g_xpset( ib )
      end   if    !  realign shell
      end   if    !  guess function located
      end   do    !  shell functions
      end   do    !  atom shells
      end   do    !  bond atoms
      end   if    !  skip NDF's
      end   do    !  guess functions
      xplen  =  nf_orb 
      end




!     find some way to define a plane in order to orient a 'Pi' 
!     bond or lone pair, returning its point (p1,p2,p3)

      subroutine  get_plane ( bat1, bat2, bond_atoms, nbond_orbs,  &
                              bondat, nbond_at, coords, p1,p2,p3 ) 
      implicit    none
      integer     bat1, bat2, bond_atoms( 2, *)
      integer     nbond_orbs, bondat(*), nbond_at, batp
      real*8      coords( 3, *), p1,p2,p3
      real*8      unlikely_linear
      parameter ( unlikely_linear = 0.5432d+00 )
      logical     linear

!     default to an arbitrary reference point
!  is there a smarter choice ? COM ?

      p1 = unlikely_linear
      p2 = unlikely_linear
      p3 = unlikely_linear

      call  get_bondats ( bat1, bond_atoms, nbond_orbs,  &
                          bondat, nbond_at ) 

      if ( nbond_at .gt. 1 ) then
      batp = bondat( 1 )
!     avoid choosing the original bonded atom
      if ( batp .eq. bat2 ) batp = bondat( 2 )
      bondat( 1 ) = bat1
      bondat( 2 ) = bat2
      bondat( 3 ) = batp 
      if ( .not. linear( coords, bondat, 3 ) ) then

!     nbond_at >2 will pass this test anyway

      p1 = coords( 1, batp )
      p2 = coords( 2, batp )
      p3 = coords( 3, batp )
      end  if

!     for nbond_at = 1, default

      end  if
      end





      integer function shell_size( l_mom )
      implicit    none
      integer     l_mom
      shell_size  =  ( l_mom + 1 )*( l_mom + 2 )/2
      end


      module       integrals
        use tools, only: dp
      implicit     none

      integer      nelec,norbs,nalpha,nbeta, hdim,nstore
      integer      max_obs !< largest possible number of basis functions in a OLCAO,
                           !! based on the atoms specified in the OLCAO expansion
      integer      natom !< The number of atoms/point charges in the geometry
      integer      npair !< Number of spin coupled electron/orbital PAIRS
      integer      ndocc !< Number of double-occupied (DOCC) orbitals
      integer      totlen !< Total length of the orbital weight list in wavefunction
      integer      nunpd !< Number of unpaired electrons/orbitals
      integer      ndf !< Number of derived basis functions (LCAO-type)
      integer      nxorb !< Number of orbital excitations
      integer      natom_t !< The number of atom types
      integer      nspinc !< Number of spin couplings
      integer      nang !< Highest angular momentum in the basis set
      integer      xpmax !< Length of the largest orbital expansion
      integer      num_sh !< Number of unique atomic basis set shells
      integer      num_pr !< Number of unique atomic basis set primitives

      integer      nset,mxctr,ntol_e_min, ntol_e_max, max_iter
      integer,  allocatable ::            orbset( :, : )

!     molecular geometry in terms of unique atom 'type' (input)

      integer,  allocatable ::  atom_t( : )
      real(dp),  allocatable ::  coords( : , : )

!     basis sets for the unique atom types

      integer,  allocatable ::  map_atom2shell( : ) !< holds indexing for the basis function shells for a given the atom type.
      integer,  allocatable ::  num_shell_atom( : ) !< holds the number of basis function shells for a given atom type
      integer,  allocatable ::  map_shell2prim( : ) !< holds indexing into the start of the primitive arrays for a given the total shell index.
      integer,  allocatable ::         ang_mom( : ) !< holds the shell angular momentum for a given shell index
      real(dp),  allocatable ::      nuc_charge( : )
      real(dp),  allocatable ::        exponent( : ) !< holds the exponent for a primitive gaussian
      real(dp),  allocatable ::       con_coeff( : ) !< holds the contraction coefficient for a primitive gaussian

!     wave function information

      integer,  allocatable ::     orbas_atnum( : ) !< holds the number of atoms whose basis sets make up a given OLCAO index
      integer,  allocatable ::     orbas_atset( : , : ) !< holds the overall atom index from the input file list for a given
                                                        !!  atom index in a OLCAO expansion and the OLCAO index
      integer,  allocatable ::        map_orbs( : ) !< holds the starting index for the set of basis functions that
                                                    !! a given OLCAO is expanded in in the total OLCAO wavefunction
                                                    !! (based just on how the orbitals are listed in order in the input file)
                                                    !! can be used to index coeff() to loop through the basis function
                                                    !! coefficients for a given orbital
      integer,  allocatable ::           xpset( : ) !< for a given index of a function in the OLCAO wavefunction, returns
                                                    !! the index of the basis function in the list of all basis functions
                                                    !! associated with the atoms that the OLCAO is expanded over
      integer,  allocatable ::            xorb( : )
      integer,  allocatable ::            root( : )
      real(dp),  allocatable ::           coeff( : ) !< holds the coefficient of a basis function
                                                    !! in an OLCAO for a given basis function index in the total OLCAO wavefunction list
                                                    !! (where the list is just based on the order the orbitals are listed in
                                                    !!  the input file)


      integer,  allocatable ::   nxyz( : , : ) !< powers of x,y,z coordinates for s, p, d shells in order
      real(dp),  allocatable ::   angn( : ) !< pnym(ij) = ashl(i) * 1/(ashl(power of x)) * 1/(ashl(power of y)) * 1/(ashl((power of z))
                                           !! for each primitive in each angular momentum i
                                           !! (where power of x + power of y + power of z = i).
                                           !!  index ij walks over all primitives--1 for s, 3 for p, etc. 
      real(dp),  allocatable ::   ashl( : ) !< sqrt*( (2*i-1)!! ) for each given angular momentum i=0,nang
      real(dp),  allocatable ::   ashi( : ) !< 1/ashl(i) for each given angular momentum i=0,nang

!     int1e

      real(dp),  allocatable ::   dij( : )
      real(dp),  allocatable :: coeffi( : )
      real(dp),  allocatable :: coeffj( : )
      integer,  allocatable :: atom_ndf( : , : )
      integer,  allocatable ::  ndf2orb( : )
      integer,  allocatable ::    xpnew( : )

!     int2e workspaces

      real(dp),  allocatable ::   dkl( : )
      real(dp),  allocatable :: coeffk( : )
      real(dp),  allocatable :: coeffl( : )

      real(dp),  allocatable :: schwarz( : )
      real(dp),  allocatable :: eribuf( : )

!     screening
!     ctol  =  charge-cloud cutoff
!     dtol  =  density cutoff
!     itol  =  electron repulsion integral cutoff

      real(dp)   ctol, dtol, itol

      logical     spinopt,store_eri,eri_stored,dem_gs
      real(dp)     sint,hint,gint, enucrep
      real(dp),  allocatable ::   ham( : , : )
      real(dp),  allocatable ::   ovl( : , : )
      real(dp)     ptbnmax,feather

      contains 

      function  shell_size( l_mom )
      integer   shell_size
      integer   l_mom
      shell_size  =  ( l_mom + 1 )*( l_mom + 2 )/2
      end function shell_size


      end module   integrals

      module       densitywork
        use tools, only: dp
      implicit     none

!     density control code storage

      integer   dme_b(2),dme_k(2)
      integer,  allocatable :: pair_sc( : , : , : )
      real(dp),  allocatable ::        coeff_sc( : )

! for each pair of spin couplings, the spin-coupled part of bra_[ab]
! and ket_[ab] are different.
      integer,  allocatable ::  bra_a( : )!< holds the OLCAO
      !                           !! associated with a given alpha e- in the bra
      integer,  allocatable ::  bra_b( : )!< holds the OLCAO
      !                           !! associated with a given beta e- in the bra
      integer,  allocatable ::  ket_a( : ) !< same as bra_a, but with the ket
      integer,  allocatable ::  ket_b( : ) !< same as bra_b, but with the ket
      integer,  allocatable ::  bexch( : )
      integer,  allocatable ::  kexch( : )
      integer,  allocatable ::  bra( : )
      integer,  allocatable ::  ket( : )

!     determinant solver workspaces

      real(dp),  allocatable ::  wdet( : , : ) !< overlap integrals between all
                                              !! spin orbitals 
      real(dp),  allocatable ::  abra_npair( : , : ) !< holds the alpha part of the overlap matrix,
                                                    !!  stored as the overlap of the OLCAO associated
                                                    !!  with each alpha spin-coupled e- and every other spin orbital (npair, nelec)
      real(dp),  allocatable ::  bbra_npair( : , : ) !< same as abra_npair, but for beta e-
      real(dp),  allocatable ::  abra_docc_un( : , : ) !< same as abra_npair, but for the overlap of the OLCAO associated
                                                      !!  with each non-spin coupled alpha e- (ndocc+ unpaired) and all
                                                      !!  spin-coupled e- (ndocc+nunpd, npair*2)
      real(dp),  allocatable ::  bbra_docc_un( : , : ) !< same as bbra_docc_un but for the beta e-
      real(dp),  allocatable ::  aket( : , : ) !< holds the alpha part of the matrix to take the determinant of for the density
      real(dp),  allocatable ::  bket( : , : ) !< same as aket, but the beta part
      real(dp),  allocatable ::  aket_docc_un( : , : ) !< holds the docc and unpaired (non-spin coupled) part of aket
                                                      !! this never changes during the run, and is copied back
                                                      !! into aket for determinant
      real(dp),  allocatable ::  bket_docc_un( : , : ) !< same as aket_docc_un, but for bket

      integer padded_size

      integer,  allocatable ::  ipvt( : ) !< holds the pivots for the call to dgetrf 

      end module   densitywork

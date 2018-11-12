! currently, PRINT_TIMING needs MPI to work (PRINT_COUNTERS will work without MPI)

  !                   ... The XM interface ...
  !              ('XM' stands for eXecution Manager)

  ! motivation
  !     purpose of XM is to hide all the interfacing needed from 
  !     the scientific code, mainly, i/o and parallel interfacing

  !                     ... i/o functionality ...
  !     main output is on stdout
  !      all  input is on stdin

module xm

  implicit none
  integer :: valence_global_communicator
  integer :: nrank
  integer :: irank

contains

  subroutine   xm_getdims ( natom,natom_t, npair,nunpd,ndocc,  &
       totlen,xpmax,nspinc, num_sh,num_pr,nang, ndf,nset,nxorb, mxctr )

    implicit     none
#ifdef VALENCE_MPI
    include "mpif.h"
    integer      ierr, from,comm
#endif
    integer      natom,ndocc,totlen,npair,nunpd, ndf, mxctr
    integer      natom_t,num_sh,num_pr,xpmax,nspinc,nang,nset,nxorb

    integer      nproc, myrank, master

    call xm_inherit( nproc, myrank, master )

    if ( myrank .eq. master ) read *, natom,natom_t, npair,nunpd,ndocc,  &
         totlen,xpmax, nspinc, num_sh,num_pr,nang, ndf,nset,nxorb, mxctr

#ifdef VALENCE_MPI
    from = 0
    comm = valence_global_communicator

    call mpi_bcast( natom, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( ndocc, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( totlen, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( npair, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( nunpd, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( natom_t, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( num_sh, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( num_pr, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( xpmax, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( nspinc, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( nang, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( ndf, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( nset, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( nxorb, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( mxctr, 1, mpi_integer4, from, comm, ierr )
#endif

  end subroutine  xm_getdims




  !     get the main molecule input data

  subroutine  xm_input ( ntol_c, ntol_e_min_in, ntol_e_max_in, ntol_d,  &
       ntol_i, max_iter_in )
    use         integrals
    use         densitywork
    use tools, only: dp
    implicit    none
#ifdef VALENCE_MPI
    include "mpif.h"
    integer     from, comm,ierr
#endif
    !     optimization control

    integer     ntol_c, ntol_e_min_in, ntol_e_max_in, ntol_d, ntol_i
    integer     max_iter_in
    !     parallel

    integer     nproc, myrank, master

    integer     i,j,k,n
    integer     ns, np, nshell, con_length


    call xm_inherit( nproc, myrank, master )

    if ( myrank .eq. master ) then


       !     optimization control 

       read *,  ntol_c, ntol_d, ntol_i, &
            ntol_e_min_in, ntol_e_max_in, max_iter_in, ptbnmax,feather,  &
            ( orbset( 1, i ), orbset( 2, i ), i = 1, nset )


       !     input the cartesian geometry 

       do    i  =  1,  natom
          read *,  atom_t( i ), (coords( j, i ), j=1,3)
       end   do


       !     input the atomic basis set 

       ns = 1
       np = 1
       do    i  =  1,  natom_t
          map_atom2shell( i ) = ns

          read *,  nuc_charge( i ),  nshell
          num_shell_atom( i ) = nshell
          do    j  =  1,  nshell
             map_shell2prim( ns ) = np
             read *,  ang_mom( ns ), con_length

             !     avoid redundant input of unit weight for uncontracted GTO

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

       !     one more entry to enable differencing, to avoid storing
       !     the primitive counts per shell

       map_shell2prim( ns ) = np



       !     information about the N-electron wave function
       !     including spin couplings

       coeff_sc( 1 ) = 1.0d+00
       if ( npair .gt. 0 ) then
          if ( nspinc .eq. 1 ) then
             read *, ( pair_sc( i, 1, 1 ), pair_sc( i, 2, 1 ), i = 1, npair )
          else  if ( nspinc .gt. 1 ) then
             read *, ( coeff_sc( j ), ( pair_sc( i, 1, j ),  &
                  pair_sc( i, 2, j ), i = 1, npair ), j = 1, nspinc )
          end if
       end if

       !     input the excited orbitals and roots

       if( nxorb .gt. 0 ) read *, ( xorb( i ), root( i ), i = 1, nxorb )



       !     input the wave function; 
       !     first, the spin coupled (paired) orbitals

       j  =  1
       do     i = 1, 2*npair
          read *,   orbas_atnum( i ),   &
               ( orbas_atset( k, i ), k = 1, orbas_atnum( i ) ),   n
          map_orbs( i )  =  j
          read *, ( xpset( k ), coeff( k ), k = j, j + n - 1 )
          j  =  j  +  n
       end   do


       !     input unpaired orbitals

       do     i = 2*npair + 1, 2*npair + nunpd
          read *,   orbas_atnum( i ),   &
               ( orbas_atset( k, i ), k = 1, orbas_atnum( i ) ),   n
          map_orbs( i )  =  j
          read *, ( xpset( k ), coeff( k ), k = j, j + n - 1 )
          j  =  j  +  n
       end   do


       !     input doubly-occupied orbitals

       do     i = 2*npair + nunpd + 1, 2*npair + nunpd + ndocc
          read *,   orbas_atnum( i ),   &
               ( orbas_atset( k, i ), k = 1, orbas_atnum( i ) ),   n
          map_orbs( i )  =  j
          read *, ( xpset( k ), coeff( k ), k = j, j + n - 1 )
          j  =  j  +  n
       end   do


       !     input NDF's

       do     i = 2*npair + nunpd + ndocc + 1, 2*npair + nunpd + ndocc + ndf
          read *,   orbas_atnum( i ),   &
               ( orbas_atset( k, i ), k = 1, orbas_atnum( i ) ),   n
          map_orbs( i )  =  j
          read *, ( xpset( k ), coeff( k ), k = j, j + n - 1 )
          j  =  j  +  n
       end   do
       norbs  =  2*npair  +  nunpd  +  ndocc + ndf
       map_orbs( norbs  +  1 )  =  j


    end if    !  myrank = master process


#ifdef VALENCE_MPI
    !     share input data with other processes
    !  bcast the INT's as a single vector?
    from = 0
    comm = valence_global_communicator

    call mpi_bcast( ntol_c, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( ntol_e_min_in, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( ntol_e_max_in, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( ntol_d, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( ntol_i, 1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( max_iter_in, 1, mpi_integer4, from, comm, ierr )

    call mpi_bcast( xorb, nxorb, mpi_integer4, from, comm, ierr )
    call mpi_bcast( root, nxorb, mpi_integer4, from, comm, ierr )
    call mpi_bcast( orbset, 2*nset, mpi_integer4, from, comm, ierr )
    call mpi_bcast( atom_t, natom, mpi_integer4, from, comm, ierr )
    call mpi_bcast( num_shell_atom, natom_t, mpi_integer4, from, comm, ierr )
    call mpi_bcast( map_atom2shell, natom_t, mpi_integer4, from, comm, ierr )
    call mpi_bcast( map_shell2prim, num_sh+1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( ang_mom, num_sh, mpi_integer4, from, comm, ierr )
    call mpi_bcast( pair_sc, 2*npair*nspinc, mpi_integer4, from, comm, ierr )

    norbs = 2*npair + nunpd + ndocc + ndf

    call mpi_bcast( orbas_atnum, norbs, mpi_integer4, from, comm, ierr )
    call mpi_bcast( orbas_atset, norbs*mxctr, mpi_integer4, from, comm, ierr )
    call mpi_bcast( map_orbs, norbs+1, mpi_integer4, from, comm, ierr )
    call mpi_bcast( xpset, totlen, mpi_integer4, from, comm, ierr )

    call mpi_bcast( ptbnmax, 1, mpi_real8, from, comm, ierr )
    call mpi_bcast( feather, 1, mpi_real8, from, comm, ierr )
    call mpi_bcast( nuc_charge, natom_t, mpi_real8, from, comm, ierr )
    call mpi_bcast( coords, 3*natom, mpi_real8, from, comm, ierr )
    call mpi_bcast( exponent, num_pr, mpi_real8, from, comm, ierr )
    call mpi_bcast( con_coeff, num_pr, mpi_real8, from, comm, ierr )
    call mpi_bcast( coeff, totlen, mpi_real8, from, comm, ierr )
    call mpi_bcast( coeff_sc, nspinc, mpi_real8, from, comm, ierr )
#endif

  end subroutine  xm_input


  !     portable word recognition

  logical      function wordmatch ( word1, word2, nchar )
    implicit     none
    character(1) word1(*), word2(*)
    integer      nchar, i

    wordmatch = .true.
    do i = 1, nchar
       if ( word1( i ) .ne. word2( i ) ) wordmatch = .false.
    end do
  end function wordmatch



  !     generic print routine
  ! Don't forget to add ; to denote the end of message
  ! TODO: use trim, adjustl etc as an alternative to left_justify
  subroutine      xm_print ( mode, message, ints, dbls )
    use tools, only: dp
    implicit        none
    character(*)    mode, message
    integer, intent(in), optional :: ints(*)
    real(dp), intent(in), optional :: dbls(*)

    integer         ichanl,  nproc,myrank,master
    parameter     ( ichanl = 6 )     !  stdout, but variable
    integer         maxlen
    parameter     ( maxlen = 120 )
    character(maxlen)   outbuf


    call left_justify( message, outbuf, maxlen )

    !     'mode' selects a basic format 
    !     'error' mode is wanted on all processes
    call xm_inherit( nproc, myrank, master )

    if ( wordmatch( mode, 'error', 5 ) ) then
       if (present(ints)) then
         write (ichanl,1) outbuf, ints( 1 ), myrank
       else
         write (ichanl,11) outbuf, myrank
       end if
1      format(1a40, 1i8, ' from rank ',1i8)
11     format(1a40, ' from rank ',1i8)
    else


       if ( myrank .eq. master ) then

          !     one line of plain text up to maxlen characters

          if ( wordmatch( mode, 'comment', 7 ) ) then
             write (ichanl,2) outbuf
2            format(1x,1a71)

             !     place-holder for fancy announcement

          else if ( wordmatch( mode, 'title', 5 ) ) then
             write (ichanl,3) outbuf
3            format(/,1x,1a71,/)

             !     place-holder for format of a section intro

          else if ( wordmatch( mode, 'header', 6 ) ) then
             write (ichanl,4) outbuf
4            format(1x,1a71,/)

          else if ( wordmatch( mode, 'parameter', 9 ) ) then
             write (ichanl,5) outbuf, dbls(1)
5            format(1x,1a32,2x,1f24.16)

          else if ( wordmatch( mode, 'dimension', 9 ) ) then
             write (ichanl,6) outbuf, ints(1)
6            format(1x,1a32,2x,1i22)

          else if ( wordmatch( mode, 'energy', 6 ) ) then
             write (ichanl,7) ints(1), ints(2), dbls(1), dbls(2)
7            format(1i5,2x,1i4,2(2x,1f22.16))

          else if ( wordmatch( mode, 'iters', 5 ) ) then
             write (ichanl,8) ints(1), ints(2), dbls(1), dbls(2)
8            format(1i5,2x,1i4,6x,1f14.6,6x,1e12.4)

          else if ( wordmatch( mode, 'demgs', 5 ) ) then
             write (ichanl,9) ints(1), ints(2), ints(3), ints(4),  &
                  dbls(1), dbls(2), dbls(3), dbls(4)
9            format(4i5,1f14.6,1e12.4,2f9.4)

          end if         !  non-error modes 
       end if         !  myrank = master
    end if         !  error mode
  end subroutine  xm_print




  !     left-justify a sentence by padding with spaces

  subroutine   left_justify ( sentence, outbuf, maxlen )
    implicit     none
    character(1) sentence(*), outbuf(*)
    integer      maxlen

    integer      count, i
    character(1) space
    data space /' '/

    count = 1
    do while ( sentence( count ) .ne. ';'  .and.  count .le. maxlen  )
       outbuf( count ) = sentence( count )
       count = count + 1
    end do

    !     now pad the output buffer with blank spaces

    do i = count, maxlen
       outbuf( i ) = space
    end do
  end subroutine  left_justify





  subroutine    xm_output ( mode, energy,etol )
    use           integrals
    use           densitywork
    use tools, only: dp
    implicit      none
    character(*)  mode

    integer        i,j, nproc,myrank,master
    real(dp)        energy, etol



    call xm_inherit( nproc, myrank, master )
    if ( myrank .eq. master ) then

       if (  nspinc .gt. 1  ) then

          !     output N-electron wave function 

          open ( unit = 11, file = 'nelecwfn', form = 'formatted'  &
               ,      status = 'unknown' )

          if ( nspinc .eq. 1 ) then
             write ( 11, 3 ) ( pair_sc( i, 1, 1 ), pair_sc( i, 2, 1 ),  &
                  i = 1, npair )
          else  if ( nspinc .gt. 1 ) then
             do  j = 1, nspinc
                write ( 11, 4 ) coeff_sc( j ), ( pair_sc( i, 1, j ),  &
                     pair_sc( i, 2, j ), i = 1, npair )
             end do
          end if
          write ( 11,3) ( xorb( i ), root( i ), i = 1, nxorb )
          write ( 11, * )

          close  ( unit = 11 )
       end  if


       open ( unit = 10, file = 'orbitals', form = 'formatted'  &
            ,      status = 'unknown' )


       !     output spin coupled orbitals

       do     i = 1, 2*npair
          write (10,1)   orbas_atnum( i ),   &
               ( orbas_atset( j, i ), j = 1, orbas_atnum( i ) ),   &
               map_orbs( i + 1 ) - map_orbs( i )
          write (10,2) ( xpset( j ), coeff( j ),  &
               j = map_orbs( i ), map_orbs( i + 1 ) - 1 )
       end   do


       !     output unpaired orbitals

       do     i = 2*npair + 1, 2*npair + nunpd
          write (10,1)   orbas_atnum( i ),   &
               ( orbas_atset( j, i ), j = 1, orbas_atnum( i ) ),   &
               map_orbs( i + 1 ) - map_orbs( i )
          write (10,2) ( xpset( j ), coeff( j ),  &
               j = map_orbs( i ), map_orbs( i + 1 ) - 1 )
       end   do


       !     output doubly-occupied orbitals

       do     i = 2*npair + nunpd + 1, 2*npair + nunpd + ndocc
          write (10,1)   orbas_atnum( i ),   &
               ( orbas_atset( j, i ), j = 1, orbas_atnum( i ) ),   &
               map_orbs( i + 1 ) - map_orbs( i )
          write (10,2) ( xpset( j ), coeff( j ),  &
               j = map_orbs( i ), map_orbs( i + 1 ) - 1 )
       end   do


       !     output the NDF

       write (10,*)
       do     i = 2*npair + nunpd + ndocc + 1, 2*npair + nunpd + ndocc + ndf
          write (10,1)   orbas_atnum( i ),   &
               ( orbas_atset( j, i ), j = 1, orbas_atnum( i ) ),   &
               map_orbs( i + 1 ) - map_orbs( i )
          write (10,2) ( xpset( j ), coeff( j ),  &
               j = map_orbs( i ), map_orbs( i + 1 ) - 1 )
       end   do


       write (10,5) energy
       if ( mode .eq. 'done' ) write (10,7) etol

       !     trailing blank line for enhanced 'readability'

       write (10,*)
       close  ( unit = 10 )
    end  if         !  myrank = master

1   format(5x,1i2,5x,10i4)
2   format(4(1i4,1x,1f13.8))
3   format(1x,       7(2i3,2x))
4   format(1x,1f13.8,7(2i3,2x))
5   format(/,1x,'total energy in atomic units ',1f32.16)
7   format(1x,'converged to ',1e10.2 ,' kCal/mol')
  end subroutine  xm_output





  !          ... parallel functionality and wrappers ...

  !     in terms of MPI, XM currently needs:
  !                 mpi_init
  !                 mpi_finalize
  !                 mpi_comm_rank
  !                 mpi_comm_size
  !                 mpi_bcast
  !                 mpi_allreduce



  subroutine  xm_propagate( comm )
    use timing_flops
    use tools, only: dp
    implicit    none
#ifdef VALENCE_MPI
    include    'mpif.h'
    integer     ierr
    integer     int_out(2)     !  lengths not 
    real(dp)     dbl_out(2)     !   protected
#endif

    integer, intent(in), optional :: comm

#ifdef VALENCE_MPI
    ! if there is an argument, set valence_global_communicator = input
    ! otherwise, set valence_global_communicator to mpi_comm_world
    if( present( comm ) ) then
      valence_global_communicator = comm
    else
       valence_global_communicator = mpi_comm_world
       call mpi_init( ierr )
    endif

    call mpi_comm_rank( valence_global_communicator, irank, ierr )
    call mpi_comm_size( valence_global_communicator, nrank, ierr )
    int_out(1) = nrank
    call xm_print( 'dimension', 'number of processors;',  &
         int_out, dbl_out ) 
#else
    irank = 0
    nrank = 1
#endif

    kernel_time = 0.0_dp
#ifdef PRINT_TIMING
    initial_time = MPI_Wtime()
#endif
#ifdef PRINT_COUNTERS
    count_determinants = 0
#endif
  end subroutine  xm_propagate




  !     ... likewise, for tidy exit

  subroutine  xm_end( comm )
    use timing_flops  
    implicit    none
#ifdef VALENCE_MPI
    include    'mpif.h'
    integer     ierr
#endif

    integer nproc, myrank,master
    integer, intent(in), optional :: comm

#ifdef PRINT_COUNTERS
    integer(8) sum_determinants
#endif

#ifdef PRINT_TIMING
    real(dp) sum_guess_time, sum_kernel_time
    call MPI_Barrier( valence_global_communicator, ierr )
#endif

    call  xm_inherit( nproc, myrank, master )

#ifdef PRINT_TIMING
    final_time = MPI_Wtime() - initial_time
    guess_time = guess_time - initial_time
    if ( myrank .eq. 0) then
       write( *,'(A, F16.6,4X,F16.6,4x,F16.6)' ) &
       & 'Rank 0: Kernel, guess, total timings in seconds:',&
       &  kernel_time, guess_time, final_time
    endif

! get an average, too.
    call MPI_Allreduce( kernel_time, sum_kernel_time, 1, mpi_real8, mpi_sum,valence_global_communicator, ierr )
    call MPI_Allreduce( guess_time, sum_guess_time, 1, mpi_real8, mpi_sum,valence_global_communicator, ierr )
    if ( myrank .eq. 0) then
       write( *,'(A, F16.6,4X,F16.6)' ) & 
       & 'Average over ranks: Kernel, guess timings in seconds:',&
       & sum_kernel_time / nrank, sum_guess_time / nrank
    endif
#endif

#ifdef PRINT_COUNTERS
#ifdef VALENCE_MPI
    call MPI_Allreduce( count_determinants, sum_determinants, 1, mpi_integer8, mpi_sum,valence_global_communicator, ierr )
#else
    sum_determinants = count_determinants
#endif
    if ( myrank .eq. 0) then
       write( *,'(A, I20)' ) 'Average over ranks: Determinants count:', sum_determinants/ nrank
    endif
#endif


#ifdef VALENCE_MPI
    ! if the user input a communicator, don't finalize--they will finalize later
    if( .not. present( comm ) ) call mpi_finalize( ierr )
#endif

  end subroutine  xm_end





  !     share data with nearby process groups to pass on
  !     in regular parallel, this is a broadcast
  !     start global comm, add other communicators later


  subroutine    xm_share ( chtype, buff, len )
    implicit      none
    character(1)  chtype
    integer       buff(*), len

#ifdef VALENCE_MPI
    include    'mpif.h'
    integer     type, from, comm, ierr

    if ( chtype .eq. 'i' ) then
       type = mpi_integer4
    else if ( chtype .eq. 'd' ) then
       type = mpi_real8
    else 
       !     error
    end if
    from = 0
    comm = valence_global_communicator
    call mpi_bcast( buff, len, type, from, comm, ierr )
#endif
  end subroutine  xm_share



  !     equalize data within a group
  !     start with dbl data and global comm

  subroutine  xm_equalize ( buff, datalen )
    use tools, only: dp
    implicit    none
    real(dp)     buff(*)
    integer     datalen

#ifdef VALENCE_MPI
    include    'mpif.h'

    integer     type, oper, comm, ierr
    integer     maxrcv
    parameter  (maxrcv = 16384 )
    real(dp)     rcv( maxrcv )
    integer     i,j,loc,npass,len

    type = mpi_real8
    oper = mpi_sum
    comm = valence_global_communicator

    !     g'sum split over a fixed-length receive buffer

    npass = 1 + ( datalen - 1)/maxrcv
    loc  = 1
    len  = maxrcv
    do i = 1, npass
       if ( i .eq. npass ) len = datalen - maxrcv*(npass-1)
       call mpi_allreduce( buff( loc ), rcv, len, type,  &
            oper, comm, ierr )
       do j = 1, len
          buff( loc + j - 1 ) = rcv( j ) 
       end do
       loc = loc + len
    end do
#endif
  end subroutine  xm_equalize


  subroutine  xm_equalize_scalar ( buff )
    use tools, only: dp
    implicit    none
    real(dp)     buff

#ifdef VALENCE_MPI
    include    'mpif.h'

    integer     type, oper, comm, ierr
    real(dp)     rcv

    type = mpi_real8
    oper = mpi_sum
    comm = valence_global_communicator

    call mpi_allreduce( buff, rcv, 1, type,  &
      oper, comm, ierr )

    buff = rcv
#endif
  end subroutine  xm_equalize_scalar




  !     get rank, process count, etc

  subroutine  xm_inherit ( num_proc, myrank, master )
    implicit    none
#ifdef VALENCE_MPI
    include    'mpif.h'
    integer     ierr
#endif
    integer     num_proc, myrank, master

    master   = 0
#ifdef VALENCE_MPI
    call mpi_comm_rank( valence_global_communicator, myrank,   ierr ) 
    call mpi_comm_size( valence_global_communicator, num_proc, ierr )
#else
    myrank = 0
    num_proc = 1
#endif
  end subroutine  xm_inherit





  !     any parallel clean-up needed

  !     subroutine     xm_abort ( error_message, master_msg )
  subroutine     xm_abort ( error_message )
    use tools, only: dp
    implicit       none
    character(*)  error_message
    
    call xm_print( 'error', error_message)
    call xm_end()
  end subroutine  xm_abort


  subroutine  xm_dtriang ( ij, i, j )
    implicit    none
    integer     ij, i, j, k

    !     decompose a 2D triangular index into its row,col addresses

    i = sqrt( dble( 2*ij ) )
    k = ( i*i + i )/2
    do  while ( k .lt. ij ) 
       i = i + 1
       k = ( i*i + i )/2
    end do
    j  =  ij - ( i*i - i )/2

  end subroutine  xm_dtriang

  subroutine  xm_dtriang8 ( ij, i, j )
    use tools, only: dp
    implicit    none
    integer     i, j
    integer(8) ij,k,i_8, j_8

    !     decompose a 2D triangular index into its row,col addresses

    i_8 = i
    j_8 = j

    i_8 = sqrt( dble( 2*ij ) )
    k = ( i_8*i_8 + i_8 )/2
    do  while ( k .lt. ij ) 
       i_8 = i_8 + 1
       k = ( i_8*i_8 + i_8 )/2
    end do
    j_8  =  ij - ( i_8*i_8 - i_8 )/2

    i = i_8
    j = j_8

  end subroutine  xm_dtriang8



! TODO: Assumes symmetric matrix, make it general with an option.
  subroutine write_matrix(adet,max_n,n,filename)
    use tools, only: dp
    implicit none
    real(dp), intent(in) :: adet( max_n, * )
    integer,  intent(in) :: max_n, n
    character(len=*), intent(in), optional :: filename
    character :: filename2*37,filename3*37
    integer ::  i,j, myiostat
    integer, save :: counter = 0
    character :: charcounter*7
    character, save :: oldfilename*37
    if (present(filename)) then
      filename2 = filename
    else
      filename2 = 'matrix'
    end if 
    counter = counter + 1
    WRITE(charcounter,'(i7)')counter
    filename3=TRIM(ADJUSTL(filename2))//TRIM(ADJUSTL(charcounter))//'.mm'
    if (irank == 0) then
      open(unit=107, file=filename3, status='UNKNOWN', action='readwrite', &
           position='append', iostat=myiostat)
      if (myiostat /= 0) then
         write(6,*) "Warning:  I/O problem."
         return
      end if
      write(107,fmt=5) 'coordinate', ' ', 'real', ' ', 'symmetric'
    5 format('%%MatrixMarket matrix ',11A,1A,8A,1A,20A)
      write(107,*)n,n, n*(n+1)/2
      do j=1,n
         do i=1,j
            write(107,"(i10,2x,i10,8x,ES23.16)",advance="yes")i,j,adet(i,j)
         end do
      end do
      close(107)
      oldfilename = filename2
    end if
  end subroutine write_matrix
  
  subroutine write_determinant(d)
    use tools, only: dp
    implicit none
    real(dp), intent(in) ::  d
    integer :: myiostat
    if (irank == 0) then
      open(unit=109, file='determinant.txt', status='UNKNOWN', action='readwrite', &
           position='append', iostat=myiostat)
      if (myiostat /= 0) then
         write(6,*) "Warning:  IO problem for determinant.txt."
         return
      end if
      write(109,"(ES23.16)",advance="yes") d
      !close(109)
    end if
  end subroutine write_determinant

end module xm





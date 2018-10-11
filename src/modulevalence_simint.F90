#ifdef SIMINT_INT
module valence_simint
  use SimintFortran
  use tools, only: dp
  implicit none
  type(c_simint_shell), allocatable, target :: shell_map(:,:)
  real(dp), allocatable :: work(:)
  real(dp), allocatable :: integrals_store(:)
contains
  subroutine valence_initialize_simint
    use integrals
    use SimintFortran
    use iso_c_binding
    implicit none
    integer i,it,mnshi,mxshi,ii,ierr,max_shells_on_atom
    integer :: worksize, num_cart

    ! initialize the simint library
    call simint_init()

    ! allocate memory for ERIs
    worksize = simint_eri_worksize(0, nang)
    allocate(work(worksize), stat=ierr )

    num_cart = ( (nang+1)*(nang+2) )/2;
    allocate(integrals_store(num_cart*num_cart*num_cart*num_cart), stat=ierr )

    ! allocate memory for shell_map
    max_shells_on_atom = 0
    do i=1,natom
       max_shells_on_atom = max( num_shell_atom( atom_t(i) ), max_shells_on_atom )
    enddo

    allocate( shell_map(max_shells_on_atom,natom + xpmax), stat=ierr )

    ! create the shells
    ! what if you just wanted to print the basis set information?
    do i=1,natom
       !    print *, "coords", coords( 1, i ),coords( 2, i ),coords( 3, i )
       it  =  atom_t( i )
       mnshi  =         map_atom2shell( it )
       mxshi  = mnshi + num_shell_atom( it ) - 1
       !     print *, "number shells", num_shell_atom( it )
       do  ii = mnshi, mxshi
          !        print *, "shell", ii, "ang mom", ang_mom( ii )
          ! print *, "number of cartesians in shell", maxi-mini+1
          ! print *, "init", i,ii-mnshi+1
          ! print *, "init", ang_mom( ii ),  map_shell2prim( ii + 1) - map_shell2prim( ii )
          ! make a new simint_shell here
          call simint_initialize_shell( shell_map(ii-mnshi+1,i) )
          
          call simint_create_shell( map_shell2prim( ii + 1) - map_shell2prim( ii ), ang_mom( ii ), &
               coords( 1, i),coords( 2, i),coords( 3, i), &
               exponent(map_shell2prim( ii )), con_coeff(map_shell2prim( ii )), &
               shell_map(ii-mnshi+1,i) )

          ! do ig  =  map_shell2prim( ii ), map_shell2prim( ii+1 ) -1
          !    ai =  exponent( ig )
          !    ci = con_coeff( ig )
          !    print *, "for", ig, ": primitive coeffs and exponents", ci, ai
          ! enddo

       enddo
    enddo

    !call C_F_POINTER(shell_map(1,1)%alpha, p1, shape=[shell_map(1,1)%nprim])
    !call C_F_POINTER(shell_map(1,1)%coef, p2, shape=[shell_map(1,1)%nprim])
    ! write(*,*) "s_shell1  infoa"
    ! do i = 1, shell_map(1,1)%nprim
    !   write(*,*) p1(i), p2(i)
    ! end do
    
  end subroutine valence_initialize_simint

  subroutine valence_finalize_simint
    use integrals
    use SimintFortran
    implicit none
    integer i,it,mnshi,mxshi,ierr,ii

    do i=1, natom
       it  =  atom_t( i )
       mnshi  =         map_atom2shell( it )
       mxshi  = mnshi + num_shell_atom( it ) - 1
       !     print *, "number shells", num_shell_atom( it )
       do  ii = mnshi, mxshi
          !        print *, "shell", ii, "ang mom", ang_mom( ii )
          call simint_free_shell( shell_map(ii-mnshi+1,i) )
       enddo
    enddo

  deallocate( shell_map, stat = ierr )
  deallocate( work, stat = ierr )
  deallocate( integrals_store, stat=ierr )

  end subroutine valence_finalize_simint

end module valence_simint
#endif

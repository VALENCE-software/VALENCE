module valence_finit
contains
subroutine valence_finalize( comm )
  use xm
  implicit     none
  integer, intent(in), optional :: comm

  ! deallocates basis information and coordinates read from input
  call deallocate_input

#ifndef SIMINT_INT
  call rys_finalize
#endif

  if(present( comm ) ) then
     call xm_end( comm )
  else
     call xm_end
  endif
  
end subroutine valence_finalize

subroutine deallocate_input
  use          densitywork
  use          integrals
  implicit     none
  integer ierr

  deallocate(    coeff,  stat = ierr )
  deallocate(    xpset,  stat = ierr )
  deallocate( orbas_atset,  stat = ierr )
  deallocate( orbas_atnum,  stat = ierr )
  deallocate( map_orbs,  stat = ierr )

  deallocate(   orbset,  stat = ierr )
  deallocate( root,  stat = ierr )
  deallocate( xorb,  stat = ierr )
  deallocate( coeff_sc,  stat = ierr )
  deallocate(  pair_sc,  stat = ierr )

  deallocate(      con_coeff,  stat = ierr )
  deallocate(       exponent,  stat = ierr )
  deallocate(     nuc_charge,  stat = ierr )
  deallocate(        ang_mom,  stat = ierr )
  deallocate( map_shell2prim,  stat = ierr )
  deallocate( num_shell_atom,  stat = ierr )
  deallocate( map_atom2shell,  stat = ierr )
  deallocate( coords, stat = ierr )
  deallocate( atom_t, stat = ierr )

end subroutine deallocate_input
end module valence_finit

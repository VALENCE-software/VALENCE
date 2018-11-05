module tools
    integer, parameter:: dp=kind(0.d0)
    public angs2bohr,get_nuclear_repulsion_energy
contains

    subroutine angs2bohr(natom, coords)
        implicit none
        integer, intent(in) :: natom
        real(dp), intent(inout)  :: coords(:,:)
        integer             :: i
        real(dp), parameter   :: tobohrs = 1.889725987722_dp
        do    i  =  1,  natom
            coords( 1, i )  =  coords( 1, i )*tobohrs
            coords( 2, i )  =  coords( 2, i )*tobohrs
            coords( 3, i )  =  coords( 3, i )*tobohrs
        end   do
    end subroutine angs2bohr

    subroutine get_nuclear_repulsion_energy(natom,atom_t,nuc_charge,coords,nre)
        implicit none
        integer, intent(in) :: natom
        integer, intent(in) :: atom_t(:)
        real(dp), intent(in)  :: nuc_charge(:)
        real(dp), intent(in)  :: coords(:,:)
        real(dp), intent(out) :: nre
        integer :: i,j
        real(dp)  :: rsq, zij
        nre = 0.0_dp
        do    i = 2, natom
            do    j = 1, i - 1
                zij = nuc_charge( atom_t( i ) ) * nuc_charge( atom_t( j ) )
                if (  abs( zij ) .gt. 1.d-12  )  then
                    rsq = ( coords( 1, i ) - coords( 1, j ) )**2  &
                        + ( coords( 2, i ) - coords( 2, j ) )**2  &
                        + ( coords( 3, i ) - coords( 3, j ) )**2
                    if ( rsq .gt. 1.d-5 ) nre = nre + zij * rsq**( -0.5_dp )
                end   if
            end   do
        end   do
    end subroutine get_nuclear_repulsion_energy

end module tools

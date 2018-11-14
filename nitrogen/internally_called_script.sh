#!/bin/bash
# this is called from inside of valence_nitrogen_api.F90 when 
# the run is in parallel

mpirun -f $COBALT_NODEFILE -np 56 -ppn 56 $VALENCE_ROOT/bin/valence < new_input | tee out
echo "new orbitals", $VALENCE_ROOT
cat orbitals
grep "total energy" out | awk '{print $3}' > energy_output
rm out
rm new_input

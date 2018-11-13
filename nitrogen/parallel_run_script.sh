#!/bin/bash

mpirun -f $COBALT_NODEFILE -np 56 -ppn 56 /home/bertoni/VALENCE/bin/valence < new_input | tee out
echo "new orbitals"
cat orbitals
grep "total energy" out | awk '{print $3}' > energy_output
rm out
rm new_input

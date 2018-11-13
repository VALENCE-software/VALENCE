# Using VALENCE with NITROGEN

Download and install NITROGEN from: https://www.colorado.edu/nitrogen/

Nitrogen reads in a shared object which contains the energy functions
in it, so VALENCE needs to be recompiled as a shared object, instead
of as a static library.

This starts with recompiling simint as a shared object: 

For SIMINT, 

you can use 

``` $ SIMINT_EXTRA=-DBUILD_SHARED_LIBS:Bool=True ./install-simint.sh```

For VALENCE, 

``` $ BUILD_SHARED_LIB=true make ```

Since these are now shared objects, make sure to put them in the dynamic
libary path at runtime:

``` $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/user/VALENCE/simint/lib64/:/home/user/VALENCE/lib/ ```


Then write a VALENCE input file with the basis and geometry that you wish to use as an initial starting point:

```
$ cat nh3_atoms.inp 
 4 2 0 0 5    55 11 0 7 18    1 0 1 0 4    

18 18 18 8 8    900 0.0 0.0 1 5

         1           0.00000000           0.00000000           0.00000000
         2           0.00000000           0.00000000           0.99037000
         2           0.89012061           0.00000000          -0.43637193
         2          -0.69982997           0.55004883          -0.43637079


7.0   5
0   6
      4173.5110000              0.0018348        
       627.4579000              0.0139950        
       142.9021000              0.0685870        
        40.2343300              0.2322410        
        12.8202100              0.4690700        
         4.3904370              0.3604550        
0   3
        11.6263580             -0.1149610              
         2.7162800             -0.1691180              
         0.7722180              1.1458520              
1   3
        11.6263580              0.0675800        
         2.7162800              0.3239070        
         0.7722180              0.7408950        
0   1
         0.2120313              
1   1
         0.2120313
1.0    2
0   3
        18.7311370              0.03349460       
         2.8253937              0.23472695       
         0.6401217              0.81375733       
0   1
         0.1612778

      4        1   2   3   4  11
   1   -1.00564968   6    0.03554411   7    0.00031994   8    0.00092481
   9    0.00019507  10    0.00784133  11   -0.00196591  12    0.00783248
  13   -0.00195437  14    0.00783248  15   -0.00195437
      4        1   2   3   4  11
   2   -0.42302885   6   -0.44934419   7    0.00768917   8    0.02222602
   9    0.00492367  10   -0.13526290  11   -0.00523298  12   -0.13516051
  13   -0.00535214  14   -0.13516046  15   -0.00535209
      4        1   2   3   4  11
   3    0.48094532   6   -0.00677858   7    0.28336513   8    0.05309385
   9    0.01146456  10   -0.00434574  11   -0.00616443  12    0.27141267
  13    0.13558036  14   -0.22113769  15   -0.11757427
      4        1   2   3   4  11
   4   -0.56031644   6    0.02282757   7   -0.06185584   8   -0.48752806
   9   -0.03860790  10    0.01463508  11    0.02075973  12    0.01460710
  13    0.02069763  14   -0.18391356  15   -0.08133531
      4        1   2   3   4  11
   5   -0.47591291   6    0.00454225   7   -0.01133034   8   -0.03275102
   9   -0.26898996  10   -0.30124392  11   -0.15212378  12    0.13644987
  13    0.07253393  14    0.13644964  15    0.07253391
```

Then write a NITROGEN input file, as shown in the standard NITROGEN examples.
For a VALENCE run, it is important to include `USE_FORTRAN_PES`, and point to
the libvalence shared object. An example NITROGEN input file is below.

```
$ cat VSVB_STAT.job
#
# optimization and harmonic 
# frequency calculation

JOBTYPE = JOBTYPE_STAT

USE_FORTRAN_PES

# Optimizing MINimum geometry
STAT_TYPE = MIN
GRAD_STEP = 1.0E-03
PES_GRAD_STEP = 1.0E-03
FD_ORDER = 2

# Use Z-matrix coordinate system
COORD = ZMAT
%ZMAT
N
H1  1    R1
H2   1   R2   2 D1
H3   1    R3  3 D2  2 D3
%
# Masses
MASS = N H H H

STAT_MAX_ITER = 5000

# Select PES library
PES = ./lib/libvalence.so

# REF1 reference geometry is used
# as initial geometry
  REF1 = 0.99 0.99 D116 0.99 D116 D-141.8

# Select BFGS optimization
STAT_OPT = BFGS

# Set gradient norm tolerance to 1e-1
STAT_TOL = 1e-1

# Calculate frequencies after optimization
CALC_FREQ
# Print verbose optimization output
STAT_VERBOSE

```

The possibly tricky part is producing internal coordinates to match the cartesian coordinates which
were used in the VALENCE input file.

Then put the appropriate environment variable in the path, and export VALENCE_SCRIPT:

```
 $ export LD_LIBRARY_PATH=/home/user/VALENCE/lib/:/home/user/VALENCE/simint/lib/:$LD_LIBRARY_PATH
 $ export VALENCE_SCRIPT=/home/user/VALENCE/nitrogen/parallel_run_script.sh
```
 run the NIROGEN binary while passing the VALENCE input on stdin: 

`$ ~/path_to_nitrogen/nitrogen/bin/nitrogen VSVB_STAT.job < nh3_atoms.inp `

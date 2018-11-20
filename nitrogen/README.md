# Using VALENCE with NITROGEN

Download and install NITROGEN from: https://www.colorado.edu/nitrogen/

Nitrogen reads in a shared object which contains the energy functions
in it, so VALENCE needs to be recompiled as a shared object, instead
of as a static library.

This starts with recompiling simint as a shared object: 

1. Compile SIMINT as a shared object 

One option is to use 

``` $ SIMINT_EXTRA=-DBUILD_SHARED_LIBS:Bool=True ./install-simint.sh```

But you may need to modify some of the compiler flags, like CC and CXX.

2. Compile VALENCE as a shared object

As with a normal VALENCE build, this can be done in serial or parallel:

Serial:
``` $ BUILD_SHARED_LIB=true make ```

Parallel:
``` $ SEQUENTIAL=false BUILD_SHARED_LIB=true make ```

There are other options that can be changed in the Makefile, like changing
the compiler and print options.

One important option is the option to include -DVALENCE_NITROGEN_PARALLEL in FFLAGS.
If this is set, it builds the NITROGEN/VALENCE API such that 
each energy calculation from VALENCE is run in parallel, by calling
`internally_called_script.sh` from VALENCE.

Since these are now shared objects, make sure to put them in the dynamic
libary path at runtime:

``` $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/user/VALENCE/simint/lib64/:/home/user/VALENCE/lib/ ```

3. Then write a VALENCE input file with the basis and geometry that you wish
to use as an initial starting point. See subdirectories h2o/, nh3/, and ch4/
for example inputs.

``` 
$ ls */*.inp
ch4/methane_rot.inp
h2o/h2o.inp
nh3/nh3_atoms.inp
```

4. Then write a NITROGEN input file, as shown in the standard NITROGEN examples.
For a VALENCE run, it is important to include `USE_FORTRAN_PES`, and point to
the libvalence shared object.  See subdirectories h2o/, nh3/, and ch4/
for example NITROGEN inputs.

```
$ ls */VSVB*
ch4/VSVB_STAT.job
h2o/VSVB_STAT.job
nh3/VSVB_STAT.job
```

The possibly tricky part is producing internal coordinates to match the cartesian coordinates which
were used in the VALENCE input file.

5. Run the calculation. There is a run_script, `run_script.sh`, which can be modified and used to run
VALENCE/NITROGEN jobs.

This script sets the appropriate paths, and exports VALENCE_SCRIPT and VALENCE_ROOT.
The user should modify it to export VALENCE_ROOT and VALENCE_SCRIPT.
  - VALENCE_ROOT is the base directory for VALENCE.
  - VALENCE_SCRIPT is the script that will be invoked by VALENCE at runtime
    if VALENCE is compiled with -DVALENCE_NITROGEN_PARALLEL. 
    This script should call MPI and run a VALENCE energy calculation in parallel.
    `internally_called_script.sh` is an example, and `run_script.sh` currently
    invokes it.

```
$ ./run_script.sh ~/nitrogen_v1.10dev/nitrogen/bin/nitrogen VSVB_STAT.job h2o.inp
```


# How to build VALENCE

VALENCE currently uses the SIMINT library to compute integrals. 
To build SIMINT, please see the instructions given at https://github.com/simint-chem, or feel free to try the instructions below.
(Unfortunately, the code at https://www.bennyp.org/research/simint/ currently lacks one-electron integrals and cannot be used.)
For an example of how to build SIMINT see  install-simint.sh, or see the instructions below.

Requirements:
gcc 8.2.0 or later

a) code generate phase (can be done somewhere different to the library build)

1. go to SIMINT GitHub page
2. click on 'simint generator'
3. click 'clone or download'
4. click copy to clipboard icon
5. git clone [paste URL], enter
6. cd simint_generator
[Follow instructions in GitHub, needs cmake. For example ...]
7. mkdir build; cd build
8. CC=[your c compiler] CXX=[your c++ compiler] cmake ../
9. make
10. cd ..
11. python ./create.py -g build/generator/ostei -l 3 -p 3 outdir

The code is in ./outdir. This can be moved to where you need it, if necessary.

b) library build phase

[where you've placed 'outdir' ...]
12. mkdir outdir/build
13. cd outdir/build
14. CC=[your c compiler] CXX=[your c++ compiler] FC=[your fortran compiler] cmake -DSIMINT_C_FLAGS=-g  -DENABLE_FORTRAN=ON -DSIMINT_VECTOR=scalar -DCMAKE_INSTALL_PREFIX=[Choose a path to install simint, i.e. PATH_TO_VALANCE/simint] ../

Step 14 is given for a scalar processor (-DSIMINT_VECTOR=scalar). Vectorization options are available. If, for example, you are building for Intel Xeon Phi with AVX2 instructions, then use  -DSIMINT_VECTOR=avx2. 

15. make
16. make install
[only the test code uses c]

The SIMINT library path will resemble   [your path choice]/simint/libsimint.a   possibly including directories such as 'lib64' depending on the architecture.


Next, to build VALENCE itself, please edit the Makefile to reflect your choices for the following options:

1.  parallel or sequential build (SEQUENTIAL=true or false). To help with building the parallel code please see  install-mpich.sh.
2.  compiler (COMPILER_VENDOR=INTEL, GNU, ...)
3.  SIMINT base path (SIMINT_BASE)
4.  (optional) MPI library path (MPI_LIBRARY)

... type `make`

You can also try:

- `make test` to run a < 5 minute set of tests to do a quick validation of the build (recommended).
- `make test-large` to run a much longer set of tests to do a more thorough validation of the build. The longer tests may require multiprocessing (it is recommended to build the code with MPI).
- `make doc` to build Doxygen documentation in the doc/ directory, if Doxygen is installed.

Next, cd to 'vtools' and read the README file.



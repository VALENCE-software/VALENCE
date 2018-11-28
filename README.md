[![Build Status](https://travis-ci.com/VALENCE-software/VALENCE.svg?branch=master)](https://travis-ci.com/VALENCE-software/VALENCE)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/VALENCE-software/VALENCE/master?filepath=valence_tutorial.ipynb)
[![codecov](https://codecov.io/gh/VALENCE-software/VALENCE/branch/master/graph/badge.svg)](https://codecov.io/gh/VALENCE-software/VALENCE)
[![DOI](https://zenodo.org/badge/152630099.svg)](https://zenodo.org/badge/latestdoi/152630099)

# VALENCE 
A Massively Parallel Implementation of Variational Subspace Valence Bond 



16 November, 2018



## Contents:

 1. Overview
 2. Introduction to VSVB
 3. About the code
 4. Input 
 5. Output
 6. Examples
 7. Automatic input generation with vtools
 8. How to make spin-coupled orbitals
 9. How to use derived basis functions 
10. Acknowledgements
11. Contact information






## Section 1. Overview


When first downloaded, this repository includes

- .travis.yml
- How-to-build-VALENCE.md
- LICENSE
- Makefile
- README.md
- doc/
- examples/
- install-mpich.sh*
- install-simint.sh*
- nitrogen/
- src/
- testing/
- valence_tutorial.ipynb
- vtools/


Most of the above items are self-explanatory. The directories, /examples/ and /testing/, both contain input files that can be used to validate the binary. /testing/ contains an established set of inputs primarily for internal use. While there is some overlap with /testing/, /examples/ is oriented more toward educating the user about the various functions of VALENCE and features of VSVB (see Section 6 for more details). /vtools/ contains tools for automatic input generation and is described in Section 7. /doc/ contains a write-up of the method and Doxygen-generated documentation. /nitrogen/ is concerned with the interface to NITROGEN (https://www.colorado.edu/nitrogen/) for optimizing molecular geometries and computing vibrational frequencies. The next section contains a brief, practical introduction to VSVB - the method implemented by VALENCE.






## Section 2. Introduction to VSVB


In molecular electronic structure theory, variational subspace valence bond, or VSVB, is a type of generalized valence bond theory using gaussian atomic basis sets. Unlike the majority of methods in the mainstream of quantum chemistry such as Hartree-Fock, MP2, coupled-cluster, DFT, and so on, VSVB is based purely on orbitals that are allowed to overlap with one another. That is, VSVB does not use orthogonal ('molecular') orbitals. The first benefit is that VSVB orbitals tend to be highly local, typically involving just one or two atoms, in contrast to molecular orbitals which are typically delocalized over all the atoms in any given problem. The highly local orbitals have obvious advantages for chemical interpretability and computational scalability. 

At the moment, the only paper on VSVB is:

Graham D. Fletcher, "The variational subspace valence bond method", 
J. Chem. Phys. 142, 134112 (2015).
http://dx.doi.org/10.1063/1.4916743

In addition to a general background in quantum chemistry, it helps to read the above paper and references therein. The 
document doc/notes-vsvb-energy.pdf also contain detailed information about the method.

In general, VSVB orbitals are linear combinations of the atomic orbitals (LCAO) of the basis set. As mentioned above, no consideration is given to how much the VSVB orbitals overlap with one another. This gives the user complete control over how the orbitals are defined. For example, orbitals can be tailored to represent intuitive concepts in chemistry such as bonds and lone-pairs, and calculations can be performed to test these ideas. Typically, the basis set expansion of an orbital together with guesses for the LCAO expansion weights are input to an optimization run where the weight values are refined to minimize the total energy. The user is free to enter any guess because VALENCE always ensures normalization of the orbitals to machine precision, internally. Once obtained, such orbitals can be used repeatedly to build other wave functions where similar orbitals are needed. If orbitals with the desired forms are available, qualitatively correct wave functions can be made 'for free'.
 
For convenience, VSVB orbitals are grouped into three types, in order of input: spin-coupled, unpaired, and double-occupied. This categorization greatly assists efficiency. More importantly, the orbital types serve different chemical functions. The term 'double-occupied' refers to a single spatial orbital used to model a pair of electrons with opposed spins - corresponding to two spin orbitals in the wave function with the same spatial function. This situation can be thought of as analogous to that in closed-shell Hartree-Fock (HF). Indeed, when properly constructed, such a wave function (called a 'VSHF' wave function) has the HF energy, except with overlapping orbitals. In chemistry, double-occupied orbitals are often used to model atomic core electrons and other 'spectator' electrons. The 'unpaired' orbitals are singly occupied with electrons of the same spin, contributing a 'high spin' configuration. VSVB wave functions involving only unpaired and double-occupied orbitals have a single determinant. 

The term 'spin-coupled orbitals', here, refers to pairs of singly occupied orbitals whose electrons are coupled to an overall singlet using a spin function of the Rumer type. The advantage of spin-coupled orbitals is that they greatly extend the applicability of VSVB beyond 'Hartree-Fock' type problems to model bond-breaking/formation and situations where the spatial polarization of charge is critical to reproducing chemical phenomena. Spin-coupled orbitals can be used to recover a significant proportion of the so-called static electron correlation energy. Since the associated spin functions double the number of determinants in the wave function for each pair involved, spin coupled orbitals are typically used judiciously to model key components of the chemistry of interest, such as the separation of an electron pair in a bond or excited state. Occasionally, chemical problems expose different ways to couple the spins of a given group of electrons to the same overall spin, and multiple spin couplings can be incorporated into the VSVB wave function to reflect this, with an additive increase in the number of determinants.

As mentioned above, the tendency for VSVB orbitals to be highly local also brings computational scalability and efficiency. VSVB is characterized by having low memory requirements, negligible I/O and communication overheads (high parallel scalability), high CPU efficiency (in terms of the percentage of peak FLOPs), and moderate complexity (~cubic or even sub-cubic). The major consideration with VSVB is the cost of optimization which typically corresponds to many energy calculations. However, the overall cost in terms of the time-to-solution is greatly mitigated by the ability to re-use previously determined orbitals and by the high concurrency that is possible. 






## Section 3. About the code


VALENCE is the 'number-crunching' code in an overall system for executing VSVB calculations that includes various tools (the 'vtools') for generating input and processing output. Currently, VALENCE can compute energies and optimize single-reference wave functions for the following types of situation-

1. Ground states. 
2. Excited states, including multiple excitations (single-reference).
3. Closed-shell systems. 
4. Open-shell systems. 
5. Bond-breaking/formation, majority or all of the nuclear potential energy surface. 
6. Spin optimization, including resonance. 

The term 'single-reference' is used here to mean a single product of spatial orbitals. Depending on the spin-coupling used, this wave function will have one or more determinants.

In addition, VALENCE has the ability to expand orbitals in terms of other LCAOs, that is, arbitrary combinations of the atomic basis functions may be defined and used as (new) basis functions. A typical use of such 'derived basis functions' (DBF) is to provide degrees-of-freedom adapted to the molecular symmetry ('symmetry-adaptation'). Other uses include making spherical harmonic functions from the cartesian functions, and making a hybridized basis set consisting of sp,sp2, and sp3 'hybrid' functions. 

Currently, VALENCE can optimize two types of linear parameter - orbital weights, and spin-coupling weights - using a first-order method. The term, 'first-order', refers to taking the first derivative of the VSVB energy with respect to the linear weights. The method solves a generalized eigenproblem of the form  HC=SCE  beginning by forming the hamiltonian (H) and overlap matrices (S), where E are the eigenvalues (energies). The eigenvectors (C) contain the updated orbital or spin-coupling weights. The cost of this method is quadratic with the size of the orbital expansion or spin-coupling space. There is currently an option to use a 'direct energy minimization' (DEM) method, though this is still under development. 

VALENCE is written in Fortran-90 and runs in parallel using MPI. To compute integrals, VALENCE currently uses the vectorized integral package, SIMINT (https://github.com/simint-chem). This release includes instructions on how to build SIMINT. Once the SIMINT library is built, building VALENCE itself begins with setting options in a simple Makefile. The Makefile also contains some optimizations for various platforms, and an option to build the sequential form. A binary called 'valence' should result. VALENCE takes an input file on standard input. To run VALENCE sequentially just type,

`./valence < [name of input file]`

See section 5 for example input files. To run VALENCE in parallel, please consult the documentation for your target platform as to how MPI runs are specified.






## Section 4. Input


The 'valence' binary makes no choices regarding the orbitals, it merely processes the wave function specified by the input to compute energies, execute optimization procedures, and so forth. That said, a highly versatile system for defining the orbitals is supported, together with the N-electron form of the wave function, and this is described in this section.

The direct input to VALENCE is a plain-text numerical input file. No attempt is made to improve the 'look' of the input because advanced tools are included in this package to generate input automatically and offer a more user-friendly environment (see 'vtools'). Although the tools continue to be refined and improved in our research group, it is important to understand how the direct input is structured since it offers the highest degree of generality. 

In what follows, a helpful concept is that of the 'orbital basis set' (OBS). In contrast to the more familiar 'molecular' basis set, an OBS spans the basis sets of just those atoms needed to support a given spatial orbital in VSVB. OBS often involve just one or two atoms in order to represent, for example, core electrons and 'lone-pairs', or chemical bonds, respectively. Unlike a molecular basis set, as the molecular system increases in size the component OBS stay approximately the same size. Thus, OBS facilitate a concise and efficient definition of an orbital that is independent of the target molecule, allowing orbitals to be stored and re-combined to make new wave functions.

Broadly, the input to VALENCE consists of specifications for the geometry, basis set, and a guess wave function, structured to facilitate dynamic memory. Thus, all the counts/sizes/dims/lengths, etc, are given on the first line so the program can allocate sufficient memory to read the data that follow. Lines may be separated by blank lines and items on a line by spaces, as permitted by FORTRAN free-format input rules. All data items are of INTeger type, unless specified to be FLOAT/REAL. Note also that VALENCE checks few input errors, but many errors can be avoided by using the 'vtools'.

The input is organized in the order: 
A) sizes, counts, array dims  (1 line)
B) Wave function optimization control  (1 line)
C) Geometry
D) Basis set
E) N-electron Wave function information  (spin-couplings, optional)
F) Orbitals (various kinds, in order: spin-coupled, single-occupied, double-occupied) 
G) Derived basis functions 


Parts (A) through (G) are now described in more detail.

A) sizes, counts, array dims

There are 15 integers on a single line, they are-
Item  1. The number of atoms/point charges in the geometry
Item  2. The number of unique atom types
Item  3. Number of spin-coupled electron/orbital PAIRS
Item  4. Number of unpaired electrons/orbitals
Item  5. Number of double-occupied (DOCC) orbitals
Item  6. Total length of the orbital weight list (array dim)
Item  7. Length of the largest orbital expansion (array dim)
Item  8. Number of spin-couplings
Item  9. Number of unique atomic basis set shells
Item 10. Number of unique atomic basis set primitives
Item 11. Highest angular momentum in the basis set
Item 12. Number of derived basis functions (DBF)
Item 13. Number of orbital optimization groups
Item 14. Number of orbital excitations 
Item 15: Largest atom count of the orbital basis sets



B) Optimization control 

There are 8 or more items on a single line, they are-
Item 1. The charge-cloud screening tolerance (integer), e.g. '5' means 0.00001 or 1d-5 or 10^-5.
Item 2. The density screening tolerance 
Item 3. The integral screening tolerance (Schwarz inequality) 
Item 4: The coarsest wave function convergence tolerance, given in kilocalories per mole (kCal/Mol).
Item 5: The finest wave function convergence tolerance. The optimization will proceed through multiple orbital groups from coarse to fine.
Item 6: Maximum number of iterations.
Item 7: Initial weight perturbation (DEM only)
Item 8: Weight perturbation scalar (DEM only) 
Item 9: The orbital optimization groups as begin/end pairs of orbital labels (in order), e.g. " 1 3  5 6 " shows 2 groups: first group optimizes orbitals 1 through 3, second group is orbitals 5 and 6 (skipping 4)



C) Geometry

A line for each atom/point charge with the layout:

  [Atom type]  X  Y  Z 

eg.   1    0.0   0.0   0.0

The atom type is an integer that addresses the basis set(s) given in the next section. For example, type '1' addresses the first basis set listed, type '2' the second, and so on. 'X,Y,Z' refers to cartesian coordinates in Angstroms (FLOATs).



D) Basis Set

VALENCE recognizes basis sets of cartesian atom-centered contracted gaussians. The basis set for each atom/etc type is given with the layout:

Item 1: The nuclear/point charge (FLOAT)
Item 2: The number of shells (zero or more), NS

There follow NS datasets defining the shells, as follows:

Next line: 
Item 1: The shell angular momentum
Item 2: The number of primitive gaussian functions, NP, in the shell.

There follow NP lines, as follows:

Next line: 
Item 1: The primitive exponent (FLOAT)
If NP>1, 
Item 2: The primitive coefficient (FLOAT)

Note that VALENCE skips input of the redundant unit weight when NP=1. The counts, NS and NP, are iterated until all shells are input for all the atom types. A point charge can be input as a 'nuclear charge' with NS=0. This input is quite general. For example, 'floating' basis set shells can be placed anywhere using a nuclear charge of zero.



E) The 'N-electron' wave function (optional)

If the number of spin-couplings, NC (Item 8 of (A), above), is greater than zero, then the spin-coupling information will be read next. Currently, the code can make 'Rumer' type couplings for the singlet parts of a system. Each spin-coupling is read as follows. If NC=1, a list of the orbital (electron) label pairs defining the (singlet) couplings is given, e.g.  1 2  3 4, means electrons '1 and 2' are singlet coupled, then electrons '3 and 4' are singlet coupled. If NC>1, the pair list is preceeded with the spin-coupling weight (FLOAT). This is repeated for all NC spin-couplings.

Excited states may be entered next, according to Item 14 of (A), as a sequence of {NX,NR} pairs, where NX addresses the orbital to be promoted and NR labels the root of the secular equation (e.g. '0' for lowest/ground, '1' for first excited, etc). It is an error if spin-coupled orbitals are desired but no spin-couplings are input (Item 8 = 0). Without spin-coupled orbitals and/or excited states, this section is null/empty.



F) Spatial wave function

The spatial wave function in terms of the orbitals is now given. In general, the total number of orbitals is given by-

 2*[no. spin-coupled PAIRs] + [no. unpaired] + [no. DOCC]

Whatever orbitals are needed, they must be entered in the order: spin-coupled; then unpaired; then DOCC, as per the intended use. The layout for each orbital is as follows. The first line defines the orbital basis set (OBS):

Item 1.  The number of atoms, NN, whose basis sets make up the OBS.
Item 2.  List of NN atoms.
Item 3.  Total number of AO's, NA.

There follow NA  {AO address, coefficient (FLOAT)}-pairs. The AO addresses lie within the OBS specified by the NN atoms and in the order the atoms are listed. The scheme iterates until all orbitals required by the above formula, based on items 3-5 of (A), are input. 



(G)  Derived basis functions 

DBFs are input using the same OBS format as described in (F) for the main orbital types. The total number of DBFs is preempted by item 12 of (A). See section 7 for more details.




## Section 5. Output


Broadly, the output of VALENCE is the VSVB wave function and its total energy. VALENCE first prints the outcome of a 'guess' energy calculation, preceded by the nuclear repulsion energy in the case of a molecule, optionally followed by an optimization run with a progress update at each iteration. Each optimization step reports the cumulative relaxation in kCal/Mol compared to the 'guess energy'. Also printed is the relaxation obtained at that step divided by the convergence tolerance to indicate how near the optimization is to convergence. Every energy calculation or optimization step prints a file called 'orbitals' which contains the updated orbitals together with the current total energy at the end. If spin-coupled orbitals are used with more than one spin-coupling, an additional file called 'nelecwfn' is produced, containing the updated spin-coupling weights. 






## Section 6. Examples


The current release includes many examples (see directory './examples/') chosen to illustrate the main features and types of calculation that can currently be performed with VALENCE, while also being of a convenient size for verifying correct execution. The examples can be tested using the 'test_all' script and take less than five minutes on a typical processor. In /examples/, simply type-

`./test_all valence`

In this section, five examples are described in detail. Narrative text is contained within brackets where it is helpful to distinguish it from the example input text. It is also helpful to note the following: the magnitude of a linear weight can occasionally be greater than unity; the overall sign of an orbital is arbitrary as the wave function is only determined to within a phase; agreement between total energies from the same run on different hardware is rarely greater than 10 places, and 6 places is more typical.



1) Beryllium atom

This example uses two DOCC orbitals to model the singlet-S ground state of beryllium, with electronic configuration, 1s^2 2s^2, expanded over a basis set of three 's'-type functions. The variational subspace (VS) of the first orbital consists of functions 1 and 3, with function 1 as its unique degree-of-freedom (UDF). The second orbital's VS contains functions 2 and 3, with 2 as the UDF. Thus, function 3 is the 'shared' basis. Other choices exist (in fact there are two, accounting for symmetry), but this is the most efficient choice given the chemical intuition that function 1 resembles a '1s' orbital, function 2 a '2s' orbital, and so on, in accordance with the structure of a typical atomic basis set. Since there are only DOCC orbitals involved, this run optimizes a VSHF wave function.

```
[counts and array dims]
     1   1   0   0     2     4   2   0   3  17   0   0   1   0   1

[run parameters]
  16 16 16    3 3 100  0.0 0.0   1 2


[geometry]
   1   0.0  0.0  0.0


[basis set]
  4.0  3
0   8
      2940.0000000              0.0006800        
       441.2000000              0.0052360        
       100.5000000              0.0266060        
        28.4300000              0.0999930        
         9.1690000              0.2697020        
         3.1960000              0.4514690        
         1.1590000              0.2950740        
         0.1811000              0.0125870        
0   8
      2940.0000000             -0.0001230        
       441.2000000             -0.0009660        
       100.5000000             -0.0048310        
        28.4300000             -0.0193140        
         9.1690000             -0.0532800        
         3.1960000             -0.1207230        
         1.1590000             -0.1334350        
         0.1811000              0.5307670        
0   1
         0.0589000

[orbitals]
      1        1   2
   1    1.0   3    0.0
      1        1   2
   2    1.0   3    0.0
```


[A guess consisting of ones for the UDF and zeros elsewhere may be termed a 'unit guess' by analogy with a unit matrix. The output to the screen is given below.]

```
 guess energy in atomic units          -14.4746666438408660
 orbital optimization                                                   

 (full) first-order method                                              

 cycle  orbital   relaxation(kCal)   (..per orb.)/tol                   
    1     1           -0.038504       -0.3850E+02
    1     2          -36.147608       -0.3611E+05
    2     1          -58.488652       -0.2234E+05
    2     2          -60.233942       -0.1745E+04
    3     1          -61.130090       -0.8961E+03
    3     2          -61.228910       -0.9882E+02
    4     1          -61.280549       -0.5164E+02
    4     2          -61.286071       -0.5522E+01
    5     1          -61.288965       -0.2894E+01
    5     2          -61.289277       -0.3120E+00
    6     1          -61.289441       -0.1637E+00
    6     2          -61.289458       -0.1761E-01

 calculation converged                                                  

 total energy in atomic units          -14.5723376136726923
```



[After six cycles through the orbital list, the final energy (which matches the Hartree-Fock energy) is printed. The 'orbitals' file looks like this:]

```
      1    1   2
   1    1.00065097   3   -0.00375468
      1        1   2
   2    0.48912818   3    0.58002975


 total energy in atomic units             -14.5723376136726923
 converged to   0.10E-02 kCal/mol
```


[There is no 'nelecwfn' file with this run]




2)  H2O/VSHF/cc-VDZ

This example highlights chemically intuitive choices for the UDF, so just the optimized wave function and energy are given.

```
  3 2    0 0 5    30 8 0    7 25 1    0 1 0    3


  16 16 16    6 6 100  0.0 0.0   1 5


  1   0.0000   0.0000   0.1271
  2   0.0000   0.7580  -0.5085
  2   0.0000  -0.7580  -0.5085 


  8.0  5
0   8
    11720.0000000              0.0007100        
     1759.0000000              0.0054700        
      400.8000000              0.0278370        
      113.7000000              0.1048000        
       37.0300000              0.2830620        
       13.2700000              0.4487190        
        5.0250000              0.2709520        
        1.0130000              0.0154580        
0   8
    11720.0000000             -0.0001600        
     1759.0000000             -0.0012630        
      400.8000000             -0.0062670        
      113.7000000             -0.0257160        
       37.0300000             -0.0709240        
       13.2700000             -0.1654110        
        5.0250000             -0.1169550        
        1.0130000              0.5573680        
0   1
        0.3023000                     
1   3
       17.7000000              0.0430180        
        3.8540000              0.2289130        
        1.0460000              0.5087280        
1   1
        0.2753000                     

  1.0   2
0   3
       13.0100000              0.0196850        
        1.9620000              0.1379770        
        0.4446000              0.4781480        
0   1
        0.1220000                     




    3 1 2 3         6
   1    1.00085930   2    0.00438407   6    0.00072731   9    0.00346894
  11    0.00045391  13    0.00045392
    3 1 2 3         8
   2    0.08137627   5    0.30386789   6   -0.39873764   8    0.16245941
   9   -0.27280445  10    0.40165595  11    0.05882459  13   -0.02137198
    3 1 2 3         8
   2    0.08137642   5   -0.30386766   6   -0.39873750   8   -0.16245933
   9   -0.27280444  11   -0.02137131  12    0.40165620  13    0.05882437
    3 1 2 3         6
   2    0.44922797   3    0.56835690   6    0.25068561   9    0.21482539
  11   -0.03046886  13   -0.03046865
    1 1         2
   4    0.63677843   7    0.51530766   

[ The total energy is:     -75.97812747 AU ]
```

Here, the UDF for the orbitals are as follows:
        Orbital     UDF      AO label
     oxygen core:   O 1s        1
      OH(a) bond:   H(a) 1s    10
      OH(b) bond:   H(b) 1s    12
 Sigma lone-pair:   O 3s        3
    Pi lone-pair:  (O 2px       4)

Chemically sensible alternatives for the Sigma lone-pair UDF include the O2s function. Different orbitals would be obtained with this choice but the same energy. The Pi lone-pair has a different symmetry to the other orbitals as it is perpendicular to the plane of the atoms. Since it is the only orbital of its symmetry, the UDF issue is null (hence the parentheses). The in-plane orbitals have Sigma symmetry.

Note that the O 2,3py,z functions are not chemically sensible UDF for the OH bonds. The atoms are placed in the y,z plane to simplify the symmetry definitions, so the 2py,z functions are needed to direct hybridization of the oxygen valence electrons toward the hydrogens for both bonds. So choosing, say, 2py for one bond and 2pz for the other would yield unsymmetric bond orbitals.

Yet another alternative would be to hybridize the O2s,y,z (and/or the O3s,y,z) AOs to give three Osp2 functions, two directed toward their respective hydrogens and one directed in opposition to them, then base the UDF choices on them. The s/p hybridization ratio would need to be optimized for this case.




3)  H2/SC/SZ 

Hydrogen molecule with spin-coupled orbitals and a single-Zeta basis set. This is the simplest example of using spin-coupled orbitals. The (optimized) wave function will dissociate correctly when the interatomic distance is increased. As with H2O, above, just the result is given.

```
 2 1   1 0 0   4 2 1   1 3 0   0 1 0    2


 16 16 16   6 6 100  0.0 0.0   1 2


  1   0.0    0.0   0.0
  1   0.0    0.0   0.770


  1.0   1
0   3
       13.0100000              0.0196850
        1.9620000              0.1379770
        0.4446000              0.4781480


[ spin-coupling information ]
  1  2


    2 1 2         2
   1    0.89622388   2    0.17099192
    2 1 2         2
   1    0.17099207   2    0.89622378


[ The total energy is:     -1.06381067 AU ]
```

The main qualitative difference between this input and the previous ones is the presence of the spin-coupling information in the middle, indicating that electrons 1 and 2 be coupled to a singlet. 

The spatial polarization of the two orbitals toward either atom is evident in the weights in each orbital of the two 1s atomic basis functions. As the atoms are drawn apart, the weight of the local basis function tends to unity, while that of the remote function tends to zero, in each orbital, respectively, leaving two separated hydrogen atoms with a total energy slightly higher than -1.




4)  LiH/SCval/cc-pVDZ

Lithium Hydride singlet-Sigma ground state. This wave function has a double-occupied Li core and a spin-coupled Li-H 'bond'.

```
  2 2   1 0 1   30 10 1   9 27 2   0 1 0  2


  16 16 16    6 6 100  0.0 0.0   1 3


  1   0.0  0.0  0.0
  2   0.0  0.0  1.646



  3.0   6
0   8
     1469.0000000              0.0007660        
      220.5000000              0.0058920        
       50.2600000              0.0296710        
       14.2400000              0.1091800        
        4.5810000              0.2827890        
        1.5800000              0.4531230        
        0.5640000              0.2747740        
        0.0734500              0.0097510        
0   8
     1469.0000000             -0.0001200        
      220.5000000             -0.0009230        
       50.2600000             -0.0046890        
       14.2400000             -0.0176820        
        4.5810000             -0.0489020        
        1.5800000             -0.0960090        
        0.5640000             -0.1363800        
        0.0734500              0.5751020        
0   1
        0.0280500                     
1   3
        1.5340000              0.0227840        
        0.2749000              0.1391070        
        0.0736200              0.5003750        
1   1
        0.0240300                     
2   1
        0.1239000                     

  1.0   3
0   3
       13.0100000              0.0196850        
        1.9620000              0.1379770        
        0.4446000              0.4781480        
0   1
        0.1220000                     
1   1
        0.7270000                     



  1  2   



      2        1   2  10
   2    0.52988354   3    0.24604657   6    0.37767681   9    0.04693388
  10   -0.03994670  12   -0.03994670  15    0.03634386  16    0.03973432
  17    0.17856129  20   -0.00329185
      2        1   2  10
   2    0.06142962   3    0.02236453   6    0.09631637   9   -0.01160377
  10   -0.01504756  12   -0.01504756  15    0.03646087  16    0.71199723
  17    0.24851686  20   -0.01632559
      2        1   2  10
   1    0.99936284   2   -0.00202778   3   -0.00493372   6   -0.01387483
   9    0.00230970  10   -0.00033165  12   -0.00033165  15   -0.00706813
  17    0.00999584  20   -0.00141802

[ The total energy is:     -8.00046053 AU ]
```

The first two orbitals are read as the spin-coupled pair and the third is read as the double-occupied Li 1s core. As in the H2 example above, the presence of the spin-coupled Li-H bond means the wave function will dissociate correctly when the interatomic separation is increased. The bonding here is not strong. The lithium valence electron is polarized toward H.




5)  Be/2SC/cc-VQZ

This is the simplest closed-shell case with two spin-couplings. Convergence is more efficient when the orbitals are obtained with one spin-coupling first, then allowed to relax in the presence of the two couplings. Again, the choice of the first (dominant) spin-coupling is based on the chemical intuition of which electrons are paired to make bonds, lone-pairs, and so forth.


```
  1 1   2 0 0   16 4 2   5 21 0   0 1 0    1


  20 16 16   6 6 100  0.0 0.0  1 4 



  1   0.0    0.0   0.0



  4.0   5
0   9
    14630.0000000              0.0000920        
     2191.0000000              0.0007130        
      498.2000000              0.0037350        
      140.9000000              0.0154680        
       45.8600000              0.0528740        
       16.4700000              0.1456940        
        6.3190000              0.3026810        
        2.5350000              0.4049360        
        1.0350000              0.2223870        
0   9
    14630.0000000             -0.0000170        
     2191.0000000             -0.0001300        
      498.2000000             -0.0006790        
      140.9000000             -0.0028570        
       45.8600000             -0.0098130        
       16.4700000             -0.0286090        
        6.3190000             -0.0637600        
        2.5350000             -0.1172310        
        1.0350000             -0.1212020        
0   1
        0.2528000                     
0   1
        0.1052000                     
0   1
        0.0426100                     


[ spin-couplings and weights ]

    0.28613755  1  2    3  4
    0.02435745  1  3    2  4


      1        1   4
   1    2.04468084   2    1.08510278   3    0.06760375   4   -0.08225021
      1        1   4
   1    0.08795929   2    1.09357982   3    0.05323659   4   -0.06567380
      1        1   4
   2    0.40901809   3    0.66520392   4    0.38644720   5    0.11934164
      1        1   4
   2    0.14161563   3   -0.17596193   4    0.71862322   5    0.47255306


[ The total energy is:     -14.58839772 AU ]
```

This run creates a file called 'nelecwfn' with the spin-coupling information,weights given above as its contents.

The energy with one spin-coupling (1 2)(3 4) is  -14.58808261 AU, so the relaxation with two couplings is modest in this case. However, the impact of multiple spin-couplings can be much greater than this, particularly in systems with multiple near-degenerate orbitals, such as aromatics. In benzene, the Kekule and Dewar structures correspond to different spin-couplings of the six equivalent Pi electrons.






## Section 7. Automatic input generation with vtools


The previous sections describe how the input to VALENCE can be generated 'by hand'. However, it is also possible to generate this input automatically using the 'vtools' software provided in this package. A description of the vtools command syntax is provided in the file,  ./vtools/README.md. 

vtools uses the 'Model Kit' method for building wave functions, employing the analogy that the VSVB wave function can be built from orbitals much in the same way that a model airplane is built from parts. In the case of a wave function, the orbital 'parts' are contained in a repository, and the 'builder' is a binary, called 'modelkit'. modelkit currently processes three types of orbital: core orbitals, bonds, and lone-pairs. While the cores and lone-pairs involve a single atom, the bonds involve two atoms - that is, at present, only two-center bonds are supported but plans are to extend this in the future to allow multi-center orbitals. modelkit follows straightforward rules for orienting the bonds and lone-pairs with respect to the atoms in the molecule, according to the orbital type. The repository uses a 'standard' orientation to encode the orbital type, while obviating the need for storing the atomic coordinates in the orbital information, as follows:
     Z-axis : Sigma orbitals (bonds, lone-pairs)
   X,Y-axes : Pi bonds,lone-pairs (first, second, respectively)
For example, a two-center orbital involving S and Pz functions is recognized as a Sigma-bond, while the one-center counterpart would be a Sigma-lone-pair orbital. A two-center orbital involving Px functions is recognized as a Pi-bond, that with Py functions would be the second Pi-bond. And so on.  

The repository is currently limited to the major orbital types associated with H,C,N,O atoms, and the 6-31G basis set, with work to incorporate more atom types, orbital types, and basis sets on-going. Though such orbitals are strictly 'guesses', an exciting prospect is the use of machine-learning techniques to develop increasingly accurate 'guesses', with the ultimate goal of obviating the need for wave function optimization entirely.

Whenever orbitals are used (especially if they are generated for the first time), we recommend visualizing them in order to check that their form is reasonable from the chemical standpoint. To this end, we are working to incorporate a visualization capability directly into the vtools and expect this to be available in the near future.






## Section 8. How to make spin-coupled orbitals


The following procedure is recommended. It is advisable to begin a pair of spin-coupled (SC) orbitals from the corresponding converged double-occupied (DOCC) orbital. So, first, choose the DOCC orbital of interest from a previously obtained optimized wave function. In a suitable text editor, cut the DOCC orbital from the DOCC list in the input, make a copy, so there are now two identicle orbitals, and paste this pair into the list of SC orbitals. Be sure to adjust the relevant counters among the 'dims' (section 3 (A)). To initialize the subsequent optimization run for the SC pair, use the following method to provide a 'nudge' in the right direction. Locate the first- and second-largest magnitude weights in either of the SC orbitals. In one of the SC pair, increase the largest magnitude weight by 0.1, reduce the second largest by 0.1. In the other SC orbital, do the opposite - decrease the largest magnitude weight by 0.1, increase the second largest by 0.1. The value of 0.1 is just a suggestion, but this has proved to be a reasonable choice. Execute an optimization run for the pair of SC orbitals, keeping the others fixed. If necessary, re-optimize, including any other orbitals that interact significantly with the new SC pair. 






## Section 9. How to use derived basis functions (DBF)


As mentioned in section 2, VALENCE allows the user to define combinations of the atomic basis functions to form new functions over which to expand the electronic orbitals. A typical use of this feature is to symmetry-adapt the basis set to yield more convenient degrees-of-freedom. For example, making spherical harmonic functions from the cartesian functions (useful for transition metals) and hybridized basis sets consisting of sp,sp2, and sp3 'hybrid' functions. 

To use DBF, set the total number required (item 12 in part A of section 3) and append them to the existing orbitals in the OBS format. To reference the DBF in the electronic orbitals, an index less than 1 (that is, zero, or negative) is used to distinguish them from regular AOs. The DBF are indexed beginning with zero for the last one entered and proceeding backwards up the file to -1, -2, ..., 1-N, where N is the number of DBF.

In the following example, four spherical harmonic functions used to model the 3d10 configuration of copper (I) cation with a double-Zeta basis set are defined. The DBF correspond (nominally) to the 3d'z2', 4d'z2', 3dx2-y2, and 4dx2-y2 functions, with indices 0, -1, -2, and -3, respectively. The DBF are used in two of the valence orbitals of Cu+ to provide variational flexbility.

```
     1   1   0   0    14    29   3   0   7  18   2   4   1   0   1

 20 20 20   2 2 50 0.0 0.0   1 5

  1   0.0  0.0  0.0

   29.0   7
0   3
      4134.3020000              0.0631880        
       625.4912000              0.3748450        
       136.9556000              0.6831000        
0   3
       181.4960000             -0.1113200
        39.5743100              0.0944870
        12.1624600              0.9608790
1   3
       181.4960000              0.1430840        
        39.5743100              0.5677560        
        12.1624600              0.4567140        
0   3
        12.3511100             -0.2922230
         4.0496510              0.3429910
         1.2792250              0.8479460
1   3
        12.3511100              0.0277270        
         4.0496510              0.4835240        
         1.2792250              0.5929780        
2   2
        16.7593800              0.2741120        
         4.1789770              0.8446250        
2   1
         0.9943270


    1 1         2
   0    0.57405221  -1    0.62709111
    1 1         2
  -2    0.57407860  -3    0.62706573
    1 1         2
  11    0.57411273  17    0.62703292
    1 1         2
  13    0.57412360  19    0.62702246
    1 1         2
  14    0.57411685  20    0.62702896
    1 1         1
   1    1.00000000
    1 1         1
   2    1.00000000
    1 1         1
   3    1.00000000
    1 1         1
   4    1.00000000
    1 1         1
   5    1.00000000
    1 1         1
   6    1.00000000
    1 1         1
   7    1.00000000
    1 1         1
   8    1.00000000
    1 1         1
   9    1.00000000

    1 1         2
  16    0.86602540  18   -0.86602540
    1 1         2
  10    0.86602540  12   -0.86602540
    1 1         3
  21    1.00000000  16   -0.50000000  18   -0.50000000
    1 1         3
  15    1.00000000  10   -0.50000000  12   -0.50000000
```


The example below is for the ground state of dinitrogen. Bonding and anti-bonding combinations of the two 2Pz functions on each nitrogen provide well-defined UDF for the sigma-bonding orbital without losing quality in the basis set.

```
     2   1   0   0     7    76  12   0   6  15   2   2   0   0   2

 20 20 20  0 0 0 0.0 0.0  

  1   0.0000    0.0000    0.0000  
  1   0.0000    0.0000    1.0784    

  7.0   6
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
0   1
         0.2120313              
1   3
        11.6263580              0.0675800        
         2.7162800              0.3239070        
         0.7722180              0.7408950        
1   1
         0.2120313              
2   1
         0.8000000              


      2        1   2  12
   2    0.46448740   3    0.51058196   9   -0.09072744  10   -0.00618979
  12   -0.00619147  15    0.00259130  17    0.11744429  24    0.00545631
  25   -0.00741184  27   -0.00741435  30    0.02071109   0   -0.13325599
      2        1   2  12
   2   -0.11738937   9    0.00542484  10    0.00741009  12    0.00740997
  15   -0.02070910  17   -0.46450546  18   -0.51057039  24   -0.09076468
  25    0.00619669  27    0.00619689  30   -0.00257915   0   -0.13327727
      2        1   2  12
   2    0.24499689  -1    0.67364713   9    0.12273824  10   -0.01415403
  12   -0.01415444  15    0.04057664  17    0.24502759  24   -0.12270638
  25   -0.01415739  27   -0.01415748  30    0.04057488   0    0.00002984
      2        1   2   6
   4    0.43505172   7    0.24461579  13    0.04819538  19    0.43497306
  22    0.24456183  28   -0.04819615
      2        1   2   6
   5    0.43498162   8    0.24458671  14    0.04819400  20    0.43502387
  23    0.24461060  29   -0.04819252
      2        1   2  12
   1    0.99531880   2    0.02391452   9   -0.00095882  10   -0.00374759
  12   -0.00374774  15   -0.00204212  17   -0.00169077  24    0.00019574
  25   -0.00003327  27   -0.00003336  30   -0.00107139   0    0.00159460
      2        1   2  12
   2    0.00169436   9    0.00019515  10    0.00003176  12    0.00003188
  15    0.00106951  16   -0.99531842  17   -0.02391508  24   -0.00095836
  25    0.00374769  27    0.00374779  30    0.00204220   0    0.00159441

      2        1   2   2
   6    0.60840187  21   -0.60840187
      2        1   2   2
   6    0.87759376  21    0.87759376

```




## Section 10. Acknowledgements


This work was supported by Argonne Leadership Computing Facility, which is a DOE Office of Science User Facility supported under Contract DE-AC02-06CH11357. 
The name, 'VALENCE', was chosen in honor of the famous book by Charles Coulson, an early pioneer of valence bond theory.






## Section 11. Contact information


Please feel free to send questions/comments to any/all members of 'The VALENCE Group, at Argonne':

Graham D. Fletcher 
Computational Science Division 
Argonne National Laboratory 
Lemont, IL, USA
gfletcher@anl.gov

Murat Keceli 
Computational Science Division 
Argonne National Laboratory 
Lemont, IL, USA
keceli@anl.gov

Colleen Bertoni 
Argonne Leadership Computing Facility 
Argonne National Laboratory 
Lemont, IL, USA
bertoni@anl.gov

Michael D'Mello 
Intel Corporation 
425 N. Martingale Road, Suite 1500 
Schaumburg, IL, USA
mdmello@anl.gov







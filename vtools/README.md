# Vtools

Vtools provides tools to prepare input for VALENCE and parse the output. 


## Installation
Vtools is a python script so does not require installation; however,
it depends on Open Babel python bindings so that must be installed on your
system. The simplest way to install Open Babel is through Conda.

If you do not have conda, you can follow the instructions at
https://conda.io/miniconda.html to install miniconda.

With conda you can install Open Babel with:

```
conda install -c openbabel openbabel
```

You also need to copy "ebsel" python module to be able to use
an arbitrary basis set from EMSL Basis Set Exhange database. To do that
you need to clone ebsel module within vtools directory.

```
cd vtools
git clone https://github.com/jaimergp/ebsel.git
```
## Usage
```
vtools.py [-h] [-i INPUT] [-b BASIS] [-d DIRECTORY] [-f FILENAME]
                    [-o ORBITALS] [-v VALENCE] [-m MODELKIT] [--opt OPT]
                    [-A] [-F] [-G] [-M] [-R] [-W]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        INPUT can be a SMILES or InChI string or an xyz file
                        (default: C)
  -b BASIS, --basis BASIS
                        Basis set name or directory path for basis set files
                        (*.basis) and guess orbitals (*.gorb) (default: 631g)
  -d DIRECTORY, --directory DIRECTORY
                        Directory name to write basis set files and guess
                        orbitals. (default: guessdir)
  -f FILENAME, --filename FILENAME
                        Filename for VALENCE input file (default: )
  -o ORBITALS, --orbitals ORBITALS
                        Filename for the orbitals file (default: )
  --opt OPT             Optimization string: First three are thresholds for
                        charge, density and integral, respectively (6 means
                        10**-6) 4th and 5th are coarse and fine tolerances in
                        kcal/mol for orbital optimization 6th is the max
                        number of iterations to achieve self-consistence
                        (default: 10 10 10 6 6 100)
  -v VALENCE            Executable path for VALENCE (default: valence)
  --modelkit MODELKIT   Eexcutable for VALENCE modelkit tool (default: modelkit)
  -A, --allatoms        Uses all atoms and all basis functions for each
                        orbital. (default: False)
  -F, --fullbasis       Use all basis functions for each orbital (default:
                        False)
  -G, --guess           Write *.gorb and *.basis files for a given valence
                        input and orbitals (optional) Use -f to specify the
                        VALENCE input file (required) -o to specify the
                        orbitals file (optional) -d to specify the directory
                        to write files (default = 'guessdir') (default: False)
  -M, --runmodelkit     Run modelkit to adjust weights of guess orbitals
                        (default: False)
  -R, --runvalence      Run VALENCE calculation (default: False)
  -W, --write           Write files (default: False)
````

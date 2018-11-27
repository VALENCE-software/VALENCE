#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals
from . import obtools as ob
from . import iotools as io
import logging
import json
try:
    #from ebsel.ebsel import EMSL_local  #not sure why this line doesn't work with python2.7
    from ebsel.EMSL_local import EMSL_local
    ebsel = True
except:
    ebsel = False

def get_unique_list(arr,keeporder=True):
    """

    """
    uniqueset = set(arr)
    nunique = len(uniqueset)
    if keeporder:
        uniquelist = []
        for x in arr:
            if x not in uniquelist:
                uniquelist.append(x)
    else:
        uniquelist = list(set(arr))
    return uniquelist


def get_number_of_determinants(x):
    """
    Return the number of determinants for the VALENCE calculation.
    Calculated by the formula:
    Number_of_det = 3 * n**4 / 2 - n**3 + 15 * n**2 / 2 - 2 * n
    Parameters
    -----------
    x: int or any object that can return a mol object with obtools get_mol()
    """
    try:
        norbital = x + 2 - 2
    except TypeError:
        norbital = int(ob.get_nelectron(x) / 2)
    return int(3*norbital**4/2-norbital**3+15*norbital**2/2-2*norbital)


def get_unique_atomic_numbers(x):
    """
    Return a list of elements for a given chemical identifier x.
    Parameters:
    ----------
    x, *, can be any chemical identifier recognized by obtools.get_mol()
    i.e. a SMILES or an InChI string or xyz coordinates or xyz filename.
    Returns:
    --------
    list of unique integers corresponding to atomic numbers in the same order
    as given.
    """
    mol = ob.get_mol(x,make3D=True)
    atomicnumbers = []
    elements = []
    for i in range(len(mol.atoms)):
        a = mol.atoms[i]
        if a.atomicnum not in atomicnumbers:
            atomicnumbers.append(a.atomicnum)
    return atomicnumbers


def get_nlonepair(mol,unique=False):
    """

    """
    symbols = ob.get_symbol_list(mol)
    if unique:
        get_unique_list(symbols)
    n = 0
    for s in symbols:
        if s == 'O':
            n += 2
        elif s == 'N':
            n += 1
    return n


def get_lonepair_table(mol):
    """
    Return the lone pair table required by modelkit.
    """
    symbols = ob.get_symbol_list(mol)
    table = ''
    n = 0
    for i,s in enumerate(symbols):
        if s == 'O':
            n += 2
            table += '{} {}\n'.format(i+1,1)
            table += '{} {}\n'.format(i+1,2)
        elif s == 'N':
            n += 1
            table += '{} {}\n'.format(i+1,1)
    table = '{}\n'.format(n) + table
    return table


def get_valence_geo(x):
    """
    Return the coordinates of a molecule in VALENCE format.
    Parameters:
    ----------
    x, *, can be any chemical identifier recognized by obtools.get_mol()
    i.e. a SMILES or an InChI string or xyz coordinates or xyz filename.
    Returns:
    --------
    str, geometry input for VALENCE
    First column is atom type (numbering starts from 1)
    2nd, 3rd, 4th are x, y, z coordinates in angstroms.
    For C2H6:
  1          0.0000      0.0000      0.0000
  1          0.0000      0.0000      1.5045
  2          0.8652     -0.5411     -0.4043
  2         -0.9015     -0.4785     -0.4045
  2          0.0361      1.0198     -0.4048
  2         -0.8665      0.5390      1.9094
  2         -0.0331     -1.0198      1.9090
  2          0.9002      0.4810      1.9088
    """
    mol = ob.get_mol(x,make3D=True)
    atomtype = 0
    atomicnumbers = get_unique_atomic_numbers(mol)
    valencegeo =''
    for i in range(len(mol.atoms)):
        a = mol.atoms[i]
        atomtype = atomicnumbers.index(a.atomicnum) + 1
        valencegeo += '{0:10g} {1:20.8f} {2:20.8f} {3:20.8f}\n'.format(atomtype,a.coords[0],a.coords[1],a.coords[2])
    return valencegeo


def generate_guess_orbital_old(atoms,bondorder,orbtype,nbfdict,fullbasis=False):
    """
    Parameters
    ----------
    atoms: list of integers,  corresponds to atomic numbers required for the guess
    bondorder: int, 0 for atoms, 1 for single bond, 2 for double bond, ...
    orbtype: int, > 0 denotes the kind of orbital
    suffix: sts, suffix for guess orbital filenames
    Returns:
    --------
    str, guess orbital
    int, number of weights in guess orbital
    """
    natom = len(atoms)
    if natom == 1:
        norb = -(-atoms[0] // 2) # Trick to get floor 3/2=2 in python
        nbf = nbfdict[atoms[0]]
        guess = ''
        if fullbasis:
            nweight = nbf
        elif orbtype < norb:
            nweight = 1
        elif orbtype == norb:
            guess = ''
            nweight = nbf - norb
        for i in range(nweight):
            if i == 0:
                weight = 1.0
            else:
                weight = 0.0
            guess += '{} {}\n'.format(orbtype+i,weight)
    elif natom == 2:
        nbf1 = nbfdict[atoms[0]]
        nbf2 = nbfdict[atoms[1]]
        if fullbasis:
            nweight = nbf1 + nbf2
        else:
            nweight = 2
            guess = '{} 0.5\n'.format(nbf1)
            guess += '{} 0.5\n'.format(nbf1+nbf2)
    else:
        print('Not implemented')
        guess = ''
        nweight = 0
    return guess, nweight


def generate_guess_orbital(atoms,orbindex,nbfdict,fullbasis=False):
    """
    Parameters
    ----------
    atoms: list of integers,  corresponds to atomic numbers required for the guess
    orbtype: int, > 0 denotes the kind of orbital
    suffix: sts, suffix for guess orbital filenames
    Returns:
    --------
    str, guess orbital
    int, number of weights in guess orbital
    """
    natom = len(atoms)
    nbf = 0
    guess = ''
    for i in range(natom):
        nbf += nbfdict[atoms[i]]
    if fullbasis:
        nweight = nbf
    else:
        nweight = natom
    for i in range(nweight):
        if i == orbindex-1:
            weight = 1.0
        else:
            weight = 0.0
        guess += '{} {}\n'.format(i+1,weight)
    return guess, nweight


def read_guess_orbital(filename,checkmaxorbital=0):
    """
    Parameters:
    -----------
    filename: str
    checkmaxorbital: int, enter the max orbital number to check validity of guess wfn
    Returns:
    --------
    gwfn: str
    norbital: int
    """
    with open(filename, 'r') as f:
        gwfn = f.read().strip()
    nweight = len(gwfn.split()) // 2
    if checkmaxorbital:
        for i,tmp in enumerate(gwfn.split()):
            if i % 2 == 0:
                orb = int(tmp)
                if orb > checkmaxorbital:
                    print("{} has an orbital ({}) larger than given max orbital {}".format(filename,orb,checkmaxorbital))
    return gwfn, nweight


def get_bonding_table(mol):
    """
    Parameters:
    -----------
    mol: str,type(<class 'pybel.Molecule'>),
    Returns:
    -------
    str: A string formatted as input for VALENCE orient tool.
    >>> print(get_bonding_table('CC'))
    7
    1 2 1
    1 3 2
    1 4 2
    1 5 2
    2 6 2
    2 7 2
    2 8 2
    """
    mol = ob.get_mol(mol,make3D=True)
    nbond = mol.OBMol.NumBonds()
    ntotalbond = 0
    bonds = []
    bondtable = ''
    for i in range(nbond):
        bond = mol.OBMol.GetBondById(i)
        bondorder = bond.GetBondOrder()
        atomidx1 = bond.GetBeginAtomIdx()
        atomidx2 = bond.GetEndAtomIdx()
        atom1 = mol.OBMol.GetAtomById(atomidx1-1)
        atom2 = mol.OBMol.GetAtomById(atomidx2-1)
        z1 = atom1.GetAtomicNum()
        z2 = atom2.GetAtomicNum()
        for orbtype in range(bondorder):
            ntotalbond += 1
            bondtype = (min(z1,z2),max(z1,z2),bondorder,orbtype)
            if bondtype in bonds:
                bondindex = bonds.index(bondtype)
            else:
                bondindex = len(bonds)
                bonds.append(bondtype)
            bondtable += '{} {} {}\n'.format(min(atomidx1,atomidx2),max(atomidx1,atomidx2),bondindex+1)
    bondtable = '{}\n'.format(ntotalbond) + bondtable
    return bondtable

def get_full_guess(mol, basis, nbfdict):
    """
    Return guess orbitals that uses all atoms as centers, and all basis functions.
    """
    mol        = ob.get_mol(mol,make3D=True)
    natom      = len(mol.atoms)
    nelectron  = int(ob.get_nelectron(mol))
    norbital   = int(nelectron / 2)
    fullbasis = True
    if nelectron % 2 == 1:
        print('WARNING: get_full_guess assumes a closed shell system, this molecule has {} electrons, {} orbitals'.format(nelectron, norbital))
    atoms = [0]*natom
    for i in range(natom):
        atom1 = mol.OBMol.GetAtomById(i)
        atoms[i] = atom1.GetAtomicNum()
    atoms = sorted(atoms,reverse=True)
    centers = ''.join(str(a+1)+ ' ' for a in range(natom))
    guess = ''
    nbf = 0
    for i in range(natom):
        nbf += nbfdict[atoms[i]]
    nweight = nbf - norbital + 1
    for i in range(norbital):
        guess += '{}    {}  {}\n'.format(natom, centers, nweight)
        guess += '{} 1.0\n'.format(i+1)
        for j in range(norbital,nbf):
            guess += '{} 0.0\n'.format(j+1)
    return guess


def get_bonding_guess(mol, basis, nbfdict=None, fullbasis=False):
    """
    Return guess orbitals for all bonding orbitals of a given molecule
    """
    mol        = ob.get_mol(mol,make3D=True)
    nbond      = mol.OBMol.NumBonds()
    ncore      = get_number_of_nonbonding_orbitals(mol)
    guess = ''
    for i in range(nbond):
        bond = mol.OBMol.GetBondById(i)
        bondorder = bond.GetBondOrder()
        atomidx1 = bond.GetBeginAtomIdx()
        atomidx2 = bond.GetEndAtomIdx()
        atom1 = mol.OBMol.GetAtomById(atomidx1-1)
        atom2 = mol.OBMol.GetAtomById(atomidx2-1)
        z1 = atom1.GetAtomicNum()
        z2 = atom2.GetAtomicNum()
        bondedatoms = sorted([z1,z2],reverse=True)
        for orbtype in range(bondorder):
            guessfilename = get_guess_filename(bondedatoms,bondorder,orbtype=orbtype+1)
            guessfilename = io.join_path(*[basis,guessfilename])
            if io.check_file(guessfilename):
                orb, nweight = read_guess_orbital(guessfilename)
            elif nbfdict is not None:
                orb, nweight = generate_guess_orbital(bondedatoms,i+orbtype+1+ncore,nbfdict,fullbasis)
            else:
                print("Error in get_bonding_guess, {} not found and nbfdict not provided".format(guessfilename))
            guess += '2    {} {} {}\n{}\n'.format(atomidx1, atomidx2,nweight,orb)
    return guess


def get_core_guess(mol, datadir='', nbfdict=None, fullbasis=False):
    """
    Return guess orbitals for all core orbitals of a given molecule
    """
    mol   = ob.get_mol(mol,make3D=True)
    natom = len(mol.atoms)
    guess = ''
    elists = get_electrons_list(mol)
    for i in range(natom):
        atomi = mol.OBMol.GetAtomById(i)
        zi = atomi.GetAtomicNum()
        natomorbital = int(elists[i][1]/2)
        for orbtype in range(1, natomorbital+1):
            guessfilename = get_guess_filename([zi], bondorder=0, orbtype=orbtype)
            guessfilename = io.join_path(*[datadir,guessfilename])
            if io.check_file(guessfilename):
                orb, nweight = read_guess_orbital(guessfilename)
            elif nbfdict is not None:
                orb, nweight = generate_guess_orbital([zi],orbtype,nbfdict,fullbasis)
            else:
                print("Error in get_core_guess, {} not found and nbfdict not provided".format(guessfilename))
            guess += '1    {} {}\n{}\n'.format(i+1,nweight,orb)
    return guess


def get_modelkit_input(mol,basis=''):
    """
    Parameters:
    -----------
    mol: str,type(<class 'pybel.Molecule'>),
    Returns:
    -------
    str: A string formatted as input for VALENCE orient tool.

    """
    mol        = ob.get_mol(mol,make3D=True)
    nbond      = mol.OBMol.NumBonds()
    uniquelist = get_unique_atomic_numbers(mol)
    natom      = len(mol.atoms)
    natomtype  = len(uniquelist)
    ntotalbond = 0
    bonds      = []
    bondtable  = ''
    guessinput = ''
    coreinput  = get_modelkit_core_orbitals(mol, basis)
    loneinput  = get_modelkit_lonepair_orbitals(mol,basis)
    for i in range(nbond):
        bond = mol.OBMol.GetBondById(i)
        bondorder = bond.GetBondOrder()
        atomidx1 = bond.GetBeginAtomIdx()
        atomidx2 = bond.GetEndAtomIdx()
        atom1 = mol.OBMol.GetAtomById(atomidx1-1)
        atom2 = mol.OBMol.GetAtomById(atomidx2-1)
        z1 = atom1.GetAtomicNum()
        z2 = atom2.GetAtomicNum()
        atomtype1 = uniquelist.index(z1) + 1
        atomtype2 = uniquelist.index(z2) + 1
        for orbtype in range(bondorder):
            ntotalbond += 1
            bondtype = (min(z1,z2),max(z1,z2),bondorder,orbtype)
            guessfilename = get_guess_filename([z1,z2],bondorder,orbtype=orbtype+1)
            guessfilename = io.join_path(*[basis,guessfilename])
            if io.check_file(guessfilename):
                guess, nweight = read_guess_orbital(guessfilename)
            else:
                print('Orient input requires guess orbital file {}'.format(guessfilename))
                raise SystemExit
            if bondtype in bonds:
                bondindex = bonds.index(bondtype)
            else:
                guessinput += '2 {} {} {}\n{}\n'.format(min(atomtype1,atomtype2),max(atomtype1,atomtype2),nweight,guess)
                bondindex = len(bonds)
                bonds.append(bondtype)
            bondtable += '{} {} {}\n'.format(min(atomidx1,atomidx2),max(atomidx1,atomidx2),bondindex+1)
    geoinput    = get_valence_geo(mol)
    basisinput = get_valence_basis(mol, basis=basis)
    nshell, nprimitive, maxangular, nbfdict = analyze_valence_basis(basisinput)
    guessinput = coreinput + '{}\n'.format(len(bonds)) + guessinput + '\n'
    bondtable  = '{}\n'.format(ntotalbond) + bondtable
    lonetable  = get_lonepair_table(mol)
    sizeline   = "{} {} 0 0 0 0 0 0 {} {} {} 0 0 0 0\n".format(natom,natomtype,nshell,nprimitive,maxangular)
    return sizeline + '\n'  + geoinput + '\n' + basisinput + '\n' + guessinput + '\n' + bondtable + '\n' + loneinput + '\n' + lonetable


def set_bond_lengths(mol,r):
    """
    For a given mol, sets bond lengths to values in array r.
    If r is a float, set all bond lengths to r.
    """
    nbond = mol.OBMol.NumBonds()
    for i in range(nbond):
        bondi = mol.OBMol.GetBondById(i)
        bondi.SetLength(r)
    return mol


def get_guess_filename(atoms,bondorder,orbtype=1,suffix='gorb'):
    """
    Parameters
    ----------
    atoms: list of integers,  corresponds to atomic numbers required for the guess
    bondorder: int, -1 for multicenter (3 or more), 0 for atoms, 1 for single bond, 2 for double bond, ...
    orbtype: int, > 0 denotes the kind of orbital
    suffix: str, suffix for guess orbital filenames
    Returns:
    --------
    str, name of the guess orbital file
    """
    atoms = sorted(atoms,reverse=True)
    symbols = []
    ncenter = len(atoms)
    for atom in atoms:
        symbols.append(ob.get_symbol(atom))
    if bondorder == -1:
        if ncenter == 3:
            filename = '{}-{}-{}_{}.{}'.format(symbols[0], symbols[1], symbols[2], orbtype, suffix)
        elif ncenter == 4:
            filename = '{}-{}-{}-{}_{}.{}'.format(symbols[0], symbols[1], symbols[2], symbols[3], orbtype, suffix)
        elif ncenter == 5:
            filename = '{}-{}-{}-{}-{}_{}.{}'.format(symbols[0], symbols[1], symbols[2], symbols[3], symbols[4], orbtype, suffix)
    elif bondorder == 0:
        filename = '{}_{}.{}'.format(symbols[0], orbtype, suffix)
    elif bondorder > 0:
        filename = '{}-{}_{}-{}.{}'.format(symbols[0], symbols[1], str(bondorder), str(orbtype), suffix)
    else:
        print('Bond order can be an integer >= -1')
    return filename


def analyze_valence_basis(basis):
    """
    Returns the shells in the basis set and the number of basis functions
    Parameters:
    -----------
    basis: str
    Returns:
    --------
    nshelltotal, nprimitivetotal, get_max_angular_momentum(allshells), nbfdict
    list of strings: Shells in the basis set
    int: Number of basis functions

    basis = " 1.0  2\n 0   2 \n  1.9620000 0.1379770 \n  0.444 0.4781480 \n 1 1 \n 0.12\n"
    basis += " 6.0  2\n 0   2 \n  1.9620000 0.1379770 \n  0.444 0.4781480 \n 1 1 \n 0.12\n"
    """
    shells = ['S', 'P', 'D', 'F', 'G']
    nfunc  = [1,    3,   6,   10,  15]
    lines = basis.splitlines()
    nbfdict = {}
    nprimitive = 0
    nshell = 0
    maxangular = 0
    z = 0
    for line in lines:
        if line.strip():
            items = line.split()
            if len(items) == 2:
                if '.' in items[1]: # exponent and coeff
                    nprimitive += 1
                else:
                    if '.' in items[0]: # atomic number and number of shells
                        z = int(float(items[0]))
                        nshell  += int(items[1])
                        nbfdict[z] = 0
                    else:
                        m = int(items[0]) # angular momentum for the shell
                        nbfdict[z] += nfunc[m]
                        if m > maxangular:
                            maxangular = m
            elif len(items) == 1:
                nprimitive += 1
            else:
                print('Error in analyze_valence_basis for line {}'.format(line))

    return nshell, nprimitive, maxangular, nbfdict


def analyze_valence_guess(guess):
    """
    Return total and max length for the orbital weights
    """
    if type(guess) is str or 'unicode':
        lines = guess.splitlines()
    elif type(guess) is list:
        if type(guess[0]) is str or 'unicode':
            lines = guess
        else:
            print('Input for analyze_valence_guess should be a str or a list of str not {}'.format(type(guess)))
    else:
        print('Input for analyze_valence_guess should be a str or a list of str not {}'.format(type(guess)))
    norbital = 0
    nweights = []
    ncenters = []
    for line in lines:
        if '.' in line:
            pass #Passing a line for the weights
        elif len(line.split()) > 1:
            try:
                nweight = int(line.split()[-1])
                ncenter = int(line.split()[0])
                norbital += 1
                ncenters.append(ncenter)
                nweights.append(nweight)
            except:
                print('Error in analyze_valence_guess while parsing line: {}'.format(line))
    return norbital, sum(nweights), max(nweights), max(ncenters)


def get_max_angular_momentum(shells):
    """
    Return the maximum angular momentum for given shells
    Parameters:
    -----------
    shells : list of str
    Returns:
    ---------
    int :

    """
    if 'G' in shells:
        mam = 4
    elif 'F' in shells:
        mam = 3
    elif 'D' in shells:
        mam = 2
    elif 'P' in shells:
        mam = 1
    elif 'S' in shells:
        mam = 0
    else:
        print("Shells should contain any of 'S', 'P', 'D', 'F', 'G'")
    return mam


def get_valence_basis(mol, basis='', ebsel=False):
    """
    """
    mol = ob.get_mol(mol,make3D=True)
    uniquelist = get_unique_atomic_numbers(mol)
    basisinput = ''
    for z in uniquelist:
        sym = ob.get_symbol(z)
        basisfile = io.join_path(*[basis,sym+'.basis'])
        if io.check_file(basisfile):
            with open(basisfile, 'r') as f:
                zbasis = f.read().strip()
        elif ebsel:
            try:
                zbasis = fetch_valence_basis(basis, sym)[0]
            except Exception as e:
                print('Cannot fetch basis set {}, exception {}'.format(basisfile, e))
                raise SystemExit
        else:
            print('Cannot find basis file {}'.format(basisfile))
            raise SystemExit
        basisinput = basisinput + '\n' + zbasis
    return basisinput + '\n'


def fetch_valence_basis(basisname, element, maxshell='F'):
    """
    Return the basis set as a string based on VALENCE input format.
    Parameters:
    -----------
    basisname: str,
    element: str, element symbol
    maxshell: max angular momentum
    Returns:
    --------
    str, list of strings, int
    str: Basis set for VALENCE input
    list of strings: Shells in the basis set
    int: Number of basis functions
    """
    shells = ['S', 'P', 'D', 'F', 'G']
    maxshell = maxshell.upper()
    nfunc  = [1,    3,   6,   10,  15]
    maxshellid = shells.index(maxshell)
    basissetdb = EMSL_local(fmt="g94")
    basis = basissetdb.get_basis(basisname,elements=[element])[0]
    lines = []
    nbasis = 0
    nbf    = 0
    num_sh = 0
    num_pr = 0
    basisshells = []
    for line in basis.splitlines():
        if '**' in line:
            pass
        else:
            items = line.split()
            if len(items) == 2:
                if items[0] == items[0].lower(): # Coefficents line
                    lines.append(line)
                else:
                    lines.append('')
            elif len(items) == 3:
                shellid = shells.index(items[0])
                if shellid <= maxshellid:
                    if items[0] == 'S':
                        nbf += 1
                    elif items[0] == 'P':
                        nbf += 3
                    elif items[0] == 'D':
                        nbf += 6
                    elif items[0] == 'F':
                        nbf += 10

                    basisshells.append(items[0])
                    line = '{} {}'.format(str(shellid),items[1])
                    lines.append(line)
                    nbasis += 1
                else:
                    break
    num = float(ob.get_atomno(element))
    lines[0] = '{} {}'.format(str(num), str(nbasis))
    return '\n'.join(lines), basisshells, nbf


def get_valence_input(mol, basis='cc-pvdz', opt='6 6 6 1 0 1', guess=None, fullbasis=False, allatoms=False):
    """
    """
    mol = ob.get_mol(mol,make3D=True)
    geo = get_valence_geo(mol)
    basisinput = get_valence_basis(mol, basis, ebsel)
    nbfdict    = analyze_valence_basis(basisinput)[-1]
    if guess is None:
        if allatoms:
            guess = get_full_guess(mol, basis, nbfdict)
        else:
            coreguess = get_core_guess(mol, basis, nbfdict, fullbasis)
            bondguess = get_bonding_guess(mol, basis, nbfdict, fullbasis)
            guess = coreguess + bondguess
    sizeline = get_size_line(mol,basisinput,guess)
    norbital = analyze_valence_guess(guess)[0]
    optline  = get_opt_line(opt, norbital)
    input = sizeline + '\n' + optline + '\n' + geo + '\n' + basisinput + '\n' + guess
    return input


def get_guess_orbital(mol, basis='cc-pvdz', fullbasis=False, allatoms=False):
    """
    """
    mol = ob.get_mol(mol,make3D=True)
    geo = get_valence_geo(mol)
    basisinput = get_valence_basis(mol, basis, ebsel)
    nbfdict    = analyze_valence_basis(basisinput)[-1]
    if allatoms:
        guess = get_full_guess(mol, basis, nbfdict)
    else:
        coreguess = get_core_guess(mol, basis, nbfdict, fullbasis)
        bondguess = get_bonding_guess(mol, basis, nbfdict, fullbasis)
        guess = coreguess + bondguess
    return guess


def parse_valence_input(filename='valence.inp', orbitalsfile='orbitals'):
    """
    Read <filename> and return all orbitals as a list of dictionaries.
    Each orbital dictionary is composed of:
    indices: a list of N integers corresponding to basis set index (lowest index=1)
    weights: a list of N floats
    atoms: a list of integers corresponding to atom index in the input geometry (lowest index=1)
    TODO: Fix bugs and complete this def.
    """
    if io.check_file(filename):
        lines = io.read_file(filename, aslines=True)
    else:
        print('Can not find {}.'.format(filename))
    readmore = True
    nline = len(lines)
    orbitals = []
    atomlist = []
    coords = [[0., 0., 0.]*natom]
    readinput = True
    readsizeline = True
    readoptline = False
    readgeometry = False
    readbasis = False
    i = 0
    d = {}
    while readinput:
        line = lines[i]
        if line.strip() == '' or line.upper()!= line.lower():
            i += 1
            continue
        else:
            tokens = line.split()
            ntoken = len(tokens)
            if readsizeline and ntoken == 15 and '.' not in line:
                try:
                    natom = int(tokens[0])
                    natomtype = int(tokens[1])
                    nspinpair = int(tokens[2])
                    nunpaired = int(tokens[3])
                    ndocc     = int(tokens[4])
                    nweight   = int(tokens[5])
                    nweightmax= int(tokens[6])
                    nspincoup = int(tokens[7])
                    norbgroup = int(tokens[12])
                    readsizeline = False
                    readoptline  = True
                except:
                    print('error in parsing size line in the given input')
                    return
                nopt = 8 + norbgroup
            elif readoptline and ntoken == nopt:
                chargescreeningtol   = 10.**(int(tokens[0]))
                densityscreeningtol  = 10.**(int(tokens[1]))
                integralscreeningtol = 10.**(int(tokens[2]))
                coarsewftol          = 10.**(int(tokens[3]))
                finewftol            = 10.**(int(tokens[4]))
                maxiter = int(tokens[5])
                initialperturbation = int(tokens[6])
                weightperturbation  = int(tokens[7])
                optgroups = []
                for i in range(norbgroup):
                    optgroups[i] = [int(tokens[8+2*i]), int(tokens[9+2*i])]
                readoptline = False
                readgeo     = True
            elif readgeo and ntoken==4:
                for j in range(natom):
                    line = lines[i+j]
                    tokens = line.split()
                    atoms[j] = (int(tokens[0]))
                    coords[j] = [float(tokens[1]),float(tokens[2]),float(tokens[3])]
                assert natomtype == max(atoms), 'Number of atom types do not match'
                readgeometry = False
                readbasis = True
            elif readbasis and ntoken==2:
                skipline = 0
                for j in range(natomtype):
                    line = lines[i+skipline]
                    tokens = line.split()
                    chargelist[j] = float(tokens[0])
                    nshelllist[j] = int(float[1])
                    for k in range(nshell):
                        line = lines[i+skipline+1]
                        tokens = line.split()



    return d

def parse_valence_basis(basis):
    """
    Parses valance basis to return nuclear charges and the basis for each atom as a list.
    """
    charges = []
    basislist = []
    idx = -1
    for line in basis.splitlines():
        if line.strip() == '' or line.upper()!= line.lower():
            continue
        else:
            tokens = line.split()
            ntoken = len(tokens)
            if ntoken == 2:
                i1 = tokens[0]
                i2 = tokens[1]
                if '.' in i1 and '.' not in i2:
                    charge = float(i1)
                    charges.append(charge)
                    basislist.append(line + '\n')
                    idx += 1
                else:
                    basislist[idx] += line + '\n'
            elif ntoken == 1:
                basislist[idx] += line + '\n'
            else:
                print('Skipping line: ', line)
    return charges, basislist


def charges2symbols(charges):
    """
    Converts nuclear charges to element symbols.
    """
    return [ob.get_symbol(int(charge)) for charge in charges]


def geo2xyz(geo, charges):
    """
    Converts valence geometry to xyz format.
    """
    lines = geo.splitlines()
    natom = len(lines)
    symbols = [ob.get_symbol(int(charge)) for charge in charges]
    xyz = str(natom) + '\n' + 'Generated by geo2xyz in obtools \n'
    for line in lines:
        tokens = line.split()
        xyz  += '{0}  {1}\t{2}\t{3}\n'.format(symbols[int(tokens[0])-1], tokens[1], tokens[2], tokens[3])
    return xyz


def write_basis_files(basis,geo):
    """
    Given the valence basis and geo informations as a text, writes *.basis file
    for each element.
    """
    charges, basislist = parse_valence_basis(basis)
    symbols = charges2symbols(charges)
    for i,symbol in enumerate(symbols):
        filename = '{}.basis'.format(symbol)
        io.write_file(basislist[i],filename)
    return


def write_guess_files(inp, orbitals=None, directory='guessdir'):
    """
    Given input and optionally orbitals (with a filename or text as a string) and a
    directory name, writes *.basis, and *.gorb files to provide initial guesses.
    """
    sizeline, optline, geo, bas, orb = get_input_sections(inp)
    if orbitals:
        orb = io.read(orbitals)
    guessdir = io.get_unique_filename(directory)
    io.mkdir(guessdir)
    io.cd(guessdir)
    charges, basislist = parse_valence_basis(bas)
    for basis in basislist:
        write_basis_files(basis,geo)
    atoms = [int(charge) for charge in charges]
    symbols = charges2symbols(charges)
    xyz = geo2xyz(geo, charges)
    mol = ob.get_mol(xyz)
    orblist = parse_orbitals(orb)
    orbtype = 1
    lastcenteranums = []
    for orbdict in orblist:
        ncenter = orbdict['ncenter']
        if ncenter > 2:
            bondorder = -1
        elif ncenter == 1:
            bondorder = 0
        else:
            bondorder = 1 # better solution reqd
        centeranums   = [mol.atoms[atomtype-1].atomicnum for atomtype in orbdict['atomtypes']]
        weights = orbdict['weights']
        indices = orbdict['indices']
        if centeranums == lastcenteranums:
            orbtype += 1
        else:
            orbtype = 1
        gorbfile = get_guess_filename(centeranums, bondorder, orbtype)
        wtext = ''
        for i,weight in enumerate(weights):
            wtext += '{} {}\n'.format(indices[i],weights[i])
        io.write_file(wtext,gorbfile)
        lastcenteranums = centeranums
    io.cd('../')
    return


def get_input_sections(filename='valence.inp'):
    """
    Read <filename> and return sizeline, optline,
    geolines, basislines, orbitallines.
    """
    if io.check_file(filename):
        lines = io.read_file(filename,aslines=True)
    else:
        print('Can not find {}.'.format(filename))
    readsize = True
    readopt = False
    readgeo = False
    readbas = False
    readorb = False
    sizeline = ''
    optline  = ''
    geolines = ''
    basislines = ''
    orbitallines = ''
    for i,line in enumerate(lines):
        if line.strip() == '' or line.upper()!= line.lower():
            continue # cycle
        else:
            tokens = line.split()
            ntoken = len(tokens)
            if readsize and ntoken == 15 and '.' not in line:
                sizeline = line
                natom = int(tokens[0])
                readsize = False
                readopt  = True
            elif readopt and ntoken >= 8:
                optline = line
                readopt = False
                readgeo = True
            elif readgeo and ntoken==4:
                geolines += line
                if len(geolines.splitlines()) == natom:
                    readgeo = False
                    readbas = True
            elif readbas:
                if ntoken < 3:
                    basislines += line
                else:
                    readbas = False
                    readorb = True
                    orbitallines += line
            elif readorb and ntoken >= 2:
                orbitallines += line
            else:
                print('Unknown line: ', line)
    return sizeline, optline, geolines, basislines, orbitallines


def parse_valence_output(out, d={}):
    """
    Parses output <out> and returns a dictionary
    """
    lines = out.splitlines()
    if len(lines) == 1:
        if io.check_file(out):
            lines = io.read_file(out,aslines=True)

    eguess, enuclear, etotal = float('nan'), float('nan'), float('nan')
    unit = 'au'
    nproc = 1
    niter = 0
    converged = False
    for i,line in enumerate(lines):
        if 'guess energy' in line:
            eguess = float(line.split()[-1])
        elif 'nuclear repulsion' in line:
            enuclear = float(line.split()[-1])
        elif 'total energy' in line:
            etotal = float(line.split()[-1])
        elif 'calculation converged' in line:
            converged = True
            lastiterline = lines[i-2]
            niter = int(lastiterline.split()[0])
        elif 'number of processors' in line:
            nproc = int(line.split()[-1])
    #d['energies'] = {'unit':'au','total':etotal,'nuclear':enuclear,'guess':eguess}
    d['nuclear_energy_au'] = enuclear
    d['guess_energy_au']   = eguess
    d['total_energy_au']   = etotal
    d['converged'] = converged
    d['niter'] = niter
    d['nproc'] = nproc
    return d


def parse_orbitals(inp):
    """
    Read <filename> and return all orbitals as a list of dictionaries.
    Each orbital dictionary is composed of:
    indices: a list of N integers corresponding to the basis set index (lowest index=1)
    weights: a list of N floats
    atoms: a list of integers corresponding to atom index in the input geometry (lowest index=1)
    """
    inp = io.read(inp)
    lines = inp.splitlines()
    readmore = True
    nline = len(lines)
    orbitals = []
    i = 0
    orbline = True
    while i < nline:
        orb = {}
        line = lines[i]
        print(line)
        if not line or line.strip() == '' or line.upper()!= line.lower():
            break
        atoms = []
        tokens = line.split()
        natom = int(tokens[0])
        atoms = [int(float(atomid)) for atomid in tokens[1:1+natom]]
        orb['atomtypes'] = atoms
        orb['ncenter'] = len(atoms)
        nweight = int(tokens[-1])
        weights = []
        indices = []
        nwfound = 0
        while nwfound < nweight:
            i += 1
            line = lines[i]
            tokens = line.split()
            nwfound += len(tokens) / 2
            for j, token in enumerate(tokens):
                if j % 2 == 0 :
                    indices.append(int(token))
                else:
                    weights.append(float(token))
        orb['indices'] = indices
        orb['weights'] = weights
        orbitals.append(orb)
        i += 1
    return orbitals


def write_json(d, filename):
    """
    Given a dictionary d and a string specifying a filename,
    writes the dictionary as a json file
    """
    import os
    import sys
    try:
        with open(filename, 'w') as f:
            json.dump(d, f, ensure_ascii=False,indent=2)
    except Exception as e:
        logging.error('Error in dumping json file')
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logging.error('Exception {}: {}, {}, {} '.format(e, exc_type, fname, exc_tb.tb_lineno))
    return


def get_valence_input_old(x, datadir = '', npair=0,nunpair=0, basisname='cc-pvdz', maxshell='P',nspinc=0,
                   opt={'ntol_c':12,'ntol_d':12,'ntol_i':12,'ntol_e_min':6,'ntol_e_max':6,'max_iter':100,'ptbnmax':0,'feather':0}):
    """
    Item  1. The number of atoms/point charges in the geometry
    Item  2. The number of atom types
    Item  3. Number of spin coupled electron/orbital PAIRS
    Item  4. Number of unpaired electrons/orbitals
    Item  5. Number of double-occupied (DOCC) orbitals
    Item  6. Total length of the orbital weight list (array dim)
    Item  7. Length of the largest orbital expansion (array dim)
    Item  8. Number of spin couplings
    Item  9. Number of unique atomic basis set shells
    Item 10. Number of unique atomic basis set primitives
    Item 11. Highest angular momentum in the basis set
    Item 12. Number of derived basis functions (LCAO-type)
    Item 13. Number of orbital optimization groups
    Item 14. Number of orbital excitations
    Item 15: Largest atom count of the orbital basis sets
    totlen: Sum of the number of orbital weights in guess wave function
    xpmax: Maximum number of orbital weights in a valance bond orbital expansion
    | XYZ       |  ES               |Guess  WFN |          |Basis set         |   0  ?     0  | Guess WFN natom |
      1      2       3    4     5      6      7      8       9      10     11    12  13   14      15
    natom,natom_t, npair,nunpd,ndocc, totlen,xpmax, nspinc, num_sh,num_pr,nang, ndf,nset,nxorb, mxctr
    c2h6
        8   2         0   0     9    70        10      0       7     25      1   0   0    0         2
    h2
        2   1         1   0     0    12        6       1       3      5     1    0   1    0         2
    h2o-vdz
        3   2         0   0     5    30        8       0       7      25    1    0   1    0         3
    """
    if nspinc != 0:
        print("Spin coupling is not allowed in this implementation")
    shells = ['S', 'P', 'D', 'F']
    maxshell = maxshell.upper()
    nang = shells.index(maxshell)
    mol = ob.get_mol(x,make3D=True)
    formula = ob.get_formula(mol)
    sformula = mol.OBMol.GetSpacedFormula().split()
    natom = len(mol.atoms)
    nunpd = 0
    ndf = 0 #??
    nset = 1
    nxorb = 0 # ??
    mxctr = min(natom,2)
    uniquelist   = get_unique_atomic_numbers(mol)
    natom_t = len(uniquelist)
    nelec = ob.get_nelectron(mol)
    mult  = ob.get_multiplicity(mol)
    if nelec % 2 == 1 and nunpair == 0:
        nunpair = 1
    ndocc = int((nelec - npair * 2 - nunpair) / 2)
    num_sh = 0
    num_pr = 0
    nbf = 0
    nbfdict = {}
    nbond = mol.OBMol.NumBonds()
    basis = ''
    unique_basis_shells = []
    print("Molecular formula: {}".format(formula))
    print("Number of atom types: {}".format(natom_t))
    print("Number of electrons: {}".format(nelec))
    print("Spin multiplicity: {}".format(mult))
    print("Number of doubly occupied orbitals: {}".format(ndocc))
    print("Number of unpaired electrons: {}".format(nunpair))
    print("Number of bonds: {}".format(nbond))
    print("Basis set: {}".format(basisname))
    print("Max angular momentum: {}".format(nang))
    print("Unique atomic numbers: {}".format(uniquelist))
    nelecs = [0]*natom
    nolecs = [0]*natom
    for i in range(natom):
        atom = mol.atoms[i]
        nelecs[i] = atom.atomicnum - atom.formalcharge
        nolecs[i] = nelecs[i]
    assert sum(nelecs) == nelec, 'Error in number of electrons {} - {}'.format(sum(nelecs),nelec)
    for z in uniquelist:
        element = ob.get_symbol(z)
        ncount  = int(sformula[sformula.index(element)+1])
        basisfile = io.join_path(*[basisname,element + '.basis'])
        print(basisfile)
        if io.check_file(basisfile,verbose=True):
            print('Reading basis set for {} from {}'.format(element,basisfile))
            z_basis, z_basis_shells, z_nbf = read_valence_basis(basisfile)
        else:
            print('Fetching basis set {} for {} from EMSL'.format(basisname,element))
            z_basis, z_basis_shells, z_nbf = fetch_valence_basis(basisname,element,maxshell=maxshell)
        nbfdict[z] = z_nbf
        unique_basis_shells.append(z_basis_shells)
        # The first line of z_basis should have z and number of shells
        z_num_sh = int(z_basis.splitlines()[0].split()[1])
        z_num_pr = len(z_basis.splitlines()) - 1 - z_num_sh
        print("\t Number of {} atoms: {}".format(element,ncount))
        print("\t Number of shells for {}: {}".format(element,z_num_sh))
        print("\t Number of primitives for {}: {}".format(element,z_num_pr))
        print("\t Number of basis functions for {}: {}".format(element,z_nbf))
        print("\t Basis shells for {}: {}".format(element,z_basis_shells))
        num_sh += z_num_sh
        num_pr += z_num_pr
        basis  += z_basis + '\n'
        nbf += z_nbf * ncount
    print("Total number of shells: {}".format(num_sh))
    print("Total number of primitives: {}".format(num_pr))
    print("Total number of basis functions: {}".format(nbf))
    print("\nVALENCE input\n")
    valence_wfn = ''
    weights = [0.]*(nbf)
    weights[0] = 1.
    numberofweights = []
    nbond = mol.OBMol.NumBonds()
    bonds = []
    for i in range(nbond):
        bond = mol.OBMol.GetBondById(i)
        bondorder = bond.GetBondOrder()
        atomidx1 = bond.GetBeginAtomIdx()
        atomidx2 = bond.GetEndAtomIdx()
        nolecs[atomidx1-1] -= bondorder
        nolecs[atomidx2-1] -= bondorder
        atom1 = mol.OBMol.GetAtomById(atomidx1-1)
        atom2 = mol.OBMol.GetAtomById(atomidx2-1)
        z1 = atom1.GetAtomicNum()
        z2 = atom2.GetAtomicNum()
        bondedatoms = sorted([z1,z2],reverse=True)
        for orbtype in range(bondorder):
            guessfilename = get_guess_filename(bondedatoms,bondorder,orbtype=orbtype+1)
            guessfilename = io.join_path(*[datadir,guessfilename])
            if io.check_file(guessfilename):
                guess, nweight = read_guess_orbital(guessfilename)
            else:
                guess, nweight = generate_guess_orbital(bondedatoms,bondorder,orbtype+1,nbfdict)
            numberofweights.append(nweight)
            valence_wfn += '2    {} {} {}\n{}\n'.format(atomidx1, atomidx2,nweight,guess)

    for i in range(natom):
        atomi = mol.OBMol.GetAtomById(i)
        zi = atomi.GetAtomicNum()
        natomorbital = int(nolecs[i]/2)
        for orbtype in range(1, natomorbital+1):
            guessfilename = get_guess_filename([zi], bondorder=0, orbtype=orbtype)
            guessfilename = io.join_path(*[datadir,guessfilename])
            if io.check_file(guessfilename):
                guess, nweight = read_guess_orbital(guessfilename)
            else:
                guess, nweight = generate_guess_orbital([zi],0,orbtype,nbfdict)
            numberofweights.append(nweight)
            valence_wfn += '1    {} {}\n{}\n'.format(i+1,nweight,guess)

    xpmax = max(numberofweights)
    totlen = sum(numberofweights)
    ntol_c = opt['ntol_c']
    ntol_d = opt['ntol_d']
    ntol_i = opt['ntol_i']
    ntol_e_min = opt['ntol_e_min']
    ntol_e_max = opt['ntol_e_max']
    max_iter = opt['max_iter']
    ptbnmax = opt['ptbnmax']
    feather = opt['feather']
    sizeline = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(natom,natom_t, npair,nunpd,ndocc, totlen,xpmax, nspinc, num_sh,num_pr,nang, ndf,nset,nxorb, mxctr)
    optline = "{} {} {} {} {} {} {} {} 1 {}\n".format(ntol_c, ntol_d, ntol_i, ntol_e_min, ntol_e_max, max_iter, ptbnmax, feather, ndocc)
    geolines = get_valence_geo(mol)
    valenceinput = sizeline + '\n'+ optline + '\n'+ geolines+ '\n' + basis + '\n' + valence_wfn + '\n'
    return valenceinput


def get_electrons_list(mol):
    """
    Return a list of lists for number of electrons and number of non-bonded electons
    for each atom in the mol object.
    """
    mol = ob.get_mol(mol,make3D=True)
    natom = len(mol.atoms)
    elists = [[0,0]]*natom
    for i in range(natom):
        atom = mol.atoms[i]
        ne = atom.atomicnum - atom.formalcharge
        elists[i] = [ne,ne]
    nbond = mol.OBMol.NumBonds()
    for i in range(nbond):
        bond = mol.OBMol.GetBondById(i)
        bondorder = bond.GetBondOrder()
        atomidx1 = bond.GetBeginAtomIdx()
        atomidx2 = bond.GetEndAtomIdx()
        elists[atomidx1-1][1] -= bondorder
        elists[atomidx2-1][1] -= bondorder
    return elists


def get_number_of_bonding_orbitals(mol):
    """
    Return the number of bonding orbitals for the given mol object.
    """
    mol = ob.get_mol(mol,make3D=True)
    natom = len(mol.atoms)
    nbond = mol.OBMol.NumBonds()
    nbonded = 0
    for i in range(nbond):
        bond = mol.OBMol.GetBondById(i)
        bondorder = bond.GetBondOrder()
        nbonded += 2 * bondorder
    return nbonded // 2


def get_number_of_nonbonding_orbitals(mol):
    """
    Return the number of nonbonding orbitals (i.e. one-center orbitals including
    core and lone pair) for the given mol object.
    """
    mol = ob.get_mol(mol,make3D=True)
    nelectron = ob.get_nelectron(mol)
    norbital = nelectron // 2
    nbonded = get_number_of_bonding_orbitals(mol)
    return norbital - nbonded


def get_modelkit_core_orbitals_old(mol,datadir=''):
    """
    """
    mol = ob.get_mol(mol,make3D=True)
    elistall = get_electrons_list(mol)
    uniquelist = []
    for elist in elistall:
        if elist in uniquelist:
            pass
        else:
            uniquelist.append(elist)
    coreinput = ''
    ncore = 0
    for i,t in enumerate(uniquelist):
        ae, be = t
        if be > 0:
            natomorbital = int(be / 2)
            for orbtype in range(1, natomorbital+1):
                ncore += 1
                guessfilename = get_guess_filename([ae], bondorder=0, orbtype=orbtype)
                guessfilename = io.join_path(*[datadir,guessfilename])
                if io.check_file(guessfilename):
                    guess, nweight = read_guess_orbital(guessfilename)
                else:
                    print('Requires guess orbital file {}'.format(guessfilename))
                coreinput += '{} {}\n{}\n'.format(i+1,nweight,guess)
    coreinput = '\n{}\n'.format(ncore) + coreinput
    return coreinput


def get_modelkit_core_orbitals(mol,datadir=''):
    """
    """
    mol = ob.get_mol(mol,make3D=True)
    elistall = get_electrons_list(mol)
    uniquelist = []
    for elist in elistall:
        if elist in uniquelist:
            pass
        else:
            uniquelist.append(elist)
    coreinput = ''
    ncore = 0
    for i,t in enumerate(uniquelist):
        ae, be = t
        if be > 0:
            orbtype = 1
            ncore += 1
            guessfilename = get_guess_filename([ae], bondorder=0, orbtype=orbtype)
            guessfilename = io.join_path(*[datadir,guessfilename])
            if io.check_file(guessfilename):
                guess, nweight = read_guess_orbital(guessfilename)
            else:
                print('Requires guess orbital file {}'.format(guessfilename))
            coreinput += '{} {}\n{}\n'.format(i+1,nweight,guess)
    coreinput = '\n{}\n'.format(ncore) + coreinput
    return coreinput


def get_modelkit_lonepair_orbitals(mol,datadir=''):
    """
    """
    mol = ob.get_mol(mol,make3D=True)
    elistall = get_electrons_list(mol)
    uniquelist = []
    for elist in elistall:
        if elist in uniquelist:
            pass
        else:
            uniquelist.append(elist)
    coreinput = ''
    lonetable = ''
    ncore = 0
    for i,t in enumerate(uniquelist):
        ae, be = t
        if be > 0:
            natomorbital = int(be / 2)
            for orbtype in range(2, natomorbital+1):
                ncore += 1
                guessfilename = get_guess_filename([ae], bondorder=0, orbtype=orbtype)
                guessfilename = io.join_path(*[datadir,guessfilename])
                if io.check_file(guessfilename):
                    guess, nweight = read_guess_orbital(guessfilename)
                else:
                    print('Requires guess orbital file {}'.format(guessfilename))
                coreinput += '{} {}\n{}\n'.format(i+1,nweight,guess)
                lonetable += '{} {}\n'.format(i+1,orbtype-1)
    coreinput = '\n{}\n'.format(ncore) + coreinput
    return coreinput


def get_opt_line(opt, norbital=1):
    """
    Item 1. The charge-cloud screening tolerance (integer)
            e.g. '5' means 0.00001 (10^-5)
    Item 2. The density screening tolerance
    Item 3. The integral screening tolerance (Schwarz inequality)
    Item 4: The coarsest wave function convergence tolerance, given
            in kilocalories per mole (kCal/Mol).
    Item 5: The finest wave function convergence tolerance.
            The optimization will proceed through multiple orbital
            groups from coarse to fine. Should be equal or larger than Item 4.
    Item 6: Maximum number of iterations.
    Item 7: Initial weight perturbation (DEM only)
    Item 8: Weight perturbation scalar (DEM only)
    Item 9: The orbital optimization groups as begin/end pairs of
            orbital labels (in order), e.g. 1 3  5 6
            shows 2 groups: first group optimizes orbitals 1 through
            3, second group is orbitals 5 and 6 (skipping 4)
    """
    if len(opt.split()) == 1 and len(opt.split(',')) == 1:
        if io.check_file(opt):
            opt = io.read_file(opt)
        else:
            print('Error in geo_opt_line, cannot read {}'.format(opt))
    if len(opt.split()) == 6:
        i1, i2, i3, i4, i5, i6 = opt.split()
    elif len(opt.split(',')) == 6:
        i1, i2, i3, i4, i5, i6 = opt.split(',')
    else:
        print('Error in geo_opt_line, cannot split {}'.format(opt))

    line  = '{} {} {} {} {}    {} '.format(i1, i2, i3, i4, i5, i6)
    line += '0.0 0.0 1 {}\n'.format(norbital)
    return line


def get_size_line(mol,basis,guess,npair= 0, nunpair=0,spincoupling=0):
    """
    The first line of the valence input determines the sizes of arrays.
    It consists of 15 integers:
    Item  1. The number of atoms/point charges in the geometry
    Item  2. The number of atom types
    Item  3. Number of spin coupled electron/orbital PAIRS
    Item  4. Number of unpaired electrons/orbitals
    Item  5. Number of double-occupied (DOCC) orbitals
    Item  6. Total length of the orbital weight list (array dim)
    Item  7. Length of the largest orbital expansion (array dim)
    Item  8. Number of spin couplings
    Item  9. Number of unique atomic basis set shells
    Item 10. Number of unique atomic basis set primitives
    Item 11. Highest angular momentum in the basis set
    Item 12. Number of derived basis functions (LCAO-type)
    Item 13. Number of orbital optimization groups
    Item 14. Number of orbital excitations
    Item 15: Largest atom count of the orbital basis sets
    """
    mol = ob.get_mol(mol,make3D=True)
    i1 = len(mol.atoms)
    i2 = len(get_unique_atomic_numbers(mol))
    i3 = npair
    nelec = ob.get_nelectron(mol)
    mult  = ob.get_multiplicity(mol)
    if nelec % 2 == 1 and nunpair == 0:
        nunpair = 1
    ndocc = int((nelec - npair * 2 - nunpair) / 2)
    i4 = nunpair
    i5 = ndocc
    norbital, i6, i7, i15 = analyze_valence_guess(guess)
    i8 = spincoupling
    i9, i10, i11 = analyze_valence_basis(basis)[0:3]
    i12 = 0
    i13 = 1
    i14 = 0
    line  = '{} {} {} {} {}    '.format(i1,i2,i3,i4,i5)
    line += '{} {} {} {} {}    '.format(i6,i7,i8,i9,i10)
    line += '{} {} {} {} {}    '.format(i11,i12,i13,i14,i15)
    line += '\n'
    return line


def get_args():
    """
    Returns args object that contains command line options.
    """
    import argparse
    parser = argparse.ArgumentParser(#formatter_class=argparse.RawDescriptionHelpFormatter,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=
    """
    August 26, 2018

    Generates VALENCE and Modelkit input files.

    """)
    parser.add_argument('-i', '--input', type=str,
                        default='C',
                        help="""INPUT can be a SMILES or InChI string or an xyz file""")
    parser.add_argument('-b', '--basis', type=str,
                        default='631g',
                        help="""Basis set name or directory path for basis set files (*.basis) and guess orbitals (*.gorb)""")
    parser.add_argument('-d', '--directory', type=str,
                        default='guessdir',
                        help="""Directory name to write basis set files and guess orbitals.""")
    parser.add_argument('-f', '--filename', type=str,
                        default='',
                        help="""Filename for the VALENCE input file""")
    parser.add_argument('-o', '--orbitals', type=str,
                        default='',
                        help="""Filename for the orbitals file""")
    parser.add_argument('--opt', type=str,
                        default='10 10 10 6 6 100',
                        help="""Optimization string: First three are thresholds for charge, density and integral, respectively (6 means 10**-6)
                        4th and 5th are coarse and fine tolerances in kcal/mol for orbital optimization
                        6th is the max number of iterations to achieve self-consistence""")
    parser.add_argument('-v', type=str,
                        default='valence',
                        help="""Executable path for VALENCE""")
    parser.add_argument('--modelkit', type=str,
                        default='modelkit',
                        help="""Eexcutable for VALENCE modelkit tool""")
    parser.add_argument('-A', '--allatoms', action='store_true',
                        help='Uses all atoms and all basis functions for each orbital.')
    parser.add_argument('-F', '--fullbasis', action='store_true',
                        help='Use all basis functions for each orbital')
    parser.add_argument('-G', '--guess', action='store_true',
                        help="""Write *.gorb and *.basis files for a given valence input and orbitals (optional)
                        Use -f to specify the VALENCE input file (required)
                        -o to specify the orbitals file (optional)
                        -d to specify the directory to write files (default = 'guessdir')""")
    parser.add_argument('-M', '--runmodelkit', action='store_true',
                        help='Run modelkit to adjust weights of guess orbitals')
    parser.add_argument('-R', '--runvalence', action='store_true',
                        help='Run VALENCE calculation')
    parser.add_argument('-W', '--write', action='store_true',
                        help='Write files')
    return parser.parse_args()


def main():
    from timeit import default_timer as timer
    args = get_args()
    parameters = vars(args)
    mol = parameters['input']
    basis = parameters['basis']
    opt = parameters['opt']
    prefix = parameters['filename']
    fullbasis = parameters['fullbasis']
    allatoms = parameters['allatoms']
    generateguess = parameters['guess']
    valenceoutput = None
    guess      = None
    mkinput    = None

    if generateguess:
        if parameters['orbitals']:
            orbitalsfile = parameters['orbitals']
        else:
            orbitalsfile = None
        write_guess_files(parameters['filename'], orbitals=orbitalsfile, directory=parameters['directory'])
        return 0
    if parameters['runmodelkit']:
        modelkit = parameters['modelkit']
    else:
        modelkit = None
    if io.check_dir(basis):
        print('Basis directory {} found'.format(basis))
    elif ebsel:
        print('Using EMSL Basis Set Exchange library.')
    else:
        print('Basis directory {} not found'.format(basis))
        print('EMSL Basis Set Exchange library not available')
        print('You need to provide basis set function files, i.e. *.basis, or within vtools directory run "git clone https://github.com/jaimergp/ebsel.git"')

    if prefix == '':
        prefix = ob.get_formula(mol)
    if modelkit:
        mkinput = get_modelkit_input(mol, basis)
        if parameters['write']:
            mkinputfile = io.get_unique_filename(prefix + '_mk.inp')
            print('Modelkit input file: {}'.format(mkinputfile))
            io.write_file(mkinput,mkinputfile)
        if io.check_file(modelkit):
            guess = io.run(mkinput, modelkit)
        else:
            print('Cannot find {}. Specify the path to modelkit executable with --modelkit <pathtomodelkit>'.format(modelkit))
    valenceinput = get_valence_input(mol,basis=basis,opt=opt,guess=guess,fullbasis=fullbasis,allatoms=allatoms)
    if parameters['write']:
        valenceinputfile = io.get_unique_filename(prefix + '_valence.inp')
        print('Valence input file: {}'.format(valenceinputfile))
        io.write_file(valenceinput,valenceinputfile)
    else:
        print('VALENCE input: \n {}'.format(valenceinput))
        if mkinput:
            print('Modelkit input: \n {}'.format(mkinput))
    if parameters['runvalence']:
        start = timer()
        try:
            valenceoutput = io.run(valenceinput, parameters['valence'])
        except Exception as e:
            print('Error in running VALENCE: {}'.format(e))
        runtime = timer() - start
        print('VALENCE run time in seconds: {}'.format(runtime))
        results = parse_valence_output(valenceoutput)
        print(results)
        print('VALENCE output: \n {}'.format(valenceoutput))
        if parameters['write']:
            valenceoutputfile = io.get_unique_filename(prefix + '_valence.out')
            jsonfile = io.get_unique_filename(prefix + '.json')
            io.write_file(valenceoutput,valenceoutputfile)
            write_json(results,jsonfile)
        else:
            print('VALENCE output: \n {}'.format(valenceoutput))
    return 0


if __name__ == "__main__":
    main()

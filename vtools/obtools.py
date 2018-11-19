#!/usr/bin/env python
import pybel
import openbabel
import os

"""
Module for simplifying and enhancing the usage of Open Babel.
Open Babel is a tool-box mainly used for cheminformatics.
It allows us to make conversions among different chemical formats,
such as inchi, smiles, xyz, zmat, etc
and extract molecular information such as force field optimized geometry,
spin, bonding information, conformers, etc.
More info on:  openbabel.org
Documentation on: http://openbabel.org/docs/current/index.html
This module is useful for a new user of Open Babel since it
provides information on the functionalities and how to use them
in python.
"""
__updated__ = "2018-01-12"
def get_periodic_table():
    """
    Return the periodic table as a list.
    Includes elements with atomic number less than 55.
    >>> pt = get_periodic_table()
    >>> print(len(pt))
    54
    """
    pt = ['X' ,
          'H' ,'He',
          'Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
          'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar'
          'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
          'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe']
    return pt


def get_symbol(atomno):
    """
    Returns the element symbol for a given atomic number.
    Returns 'X' for atomno=0
    >>> print(get_symbol(1))
    H
    """
    pt = get_periodic_table()
    return pt[atomno]


def get_symbol_list(mol):
    """
    Return the list of element symbols
    >>> print(get_symbol_list('CC'))
    ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    """
    xyz = get_xyz(mol)
    natom = get_natom(mol)
    symbols = ['']*natom
    for i,line in enumerate(xyz.splitlines()[2:]):
        tokens = line.split()
        if len(tokens) == 4:
            symbols[i] = tokens[0]
    return symbols


def get_atomno(symbol):
    """
    Return the atomic number for a given element symbol.
    >>> print(get_atomno('H')
    >>> 1
    """
    pt = get_periodic_table()
    symbol = symbol.capitalize()
    return pt.index(symbol)


def get_format(s):
    """
    Returns the Open Babel format of the given string.
    Note: It is primitive, distinguishes only xyz, smiles, inchi formats.
    >>> print(get_format('C'))
    smi
    >>> print(get_format('InChI=1S/H2O/h1H2'))
    inchi
    >>> print(get_format(get_xyz('C')))
    xyz
    """
    frm = 'unknown'
    lines = s.splitlines()
    n = len(lines)  # number of lines
    if n == 1:
        if s.startswith('InChI'):
            frm = 'inchi'
        else:
            frm = 'smi'
    else:
        try:
            natom = int(lines[0].strip())
            if n >= natom + 2:
                frm = 'xyz'
        except:
            pass

    return frm


def get_mol(s, make3D=False, mult=None):
    """
    Returns open-babel mol object from a given slabel, smiles string
    >>> mol = get_mol('[O][O]')
    >>> print(mol.formula)
    O2
    >>> print(mol.spin)
    3
    >>> mol = get_mol('InChI=1S/H2O/h1H2')
    >>> print(mol.formula)
    H2O
    >>> print(mol.spin)
    1
    """
    import pybel
    if type(s) is pybel.Molecule:
        mol = s
    elif type(s) is str or 'unicode' in str(type(s)):
        if s.endswith('.xyz'):
            mol = next(pybel.readfile('xyz', s))
        elif '_m' in s and len(s.splitlines()) == 1:
            s, mult = s.split('_m')
            mol = set_mult(s,int(mult))
        else:
            frm = get_format(s)
            mol = pybel.readstring(frm, s)
    else:
        print('Incompetible type {} for ob.get_mol'.format(type(s)))
        return None
    if make3D and mol.dim < 3:
        mol.make3D()
    if mult:
        mol.OBMol.SetTotalSpinMultiplicity(int(mult))
    return mol


def get_multiplicity(s):
    """
    Returns the spin multiplicity (2S+1) of the molecule, where S is the
    total spin angular momentum.
    Usage:
    >>> mol = get_mol('C')
    >>> get_multiplicity(mol)
    1
    >>> mols = [get_mol(s) for s in ['[CH3]', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_multiplicity(mol) for mol in mols]
    [2, 1, 3, 1]
    """
    return get_mult(s)


def get_mult(s):
    """
    Returns spin multiplicity as an integer for a given smiles or inchi string
    """
    mult = None
    if type(s) is str:
        if '_m' in s:
            try:
                mult = int(s.split('_m')[-1])
            except:
                pass
    if mult is None:
        mol = get_mol(s,make3D=False)
        mult = mol.spin
    return mult


def get_slabel(s,mult=None):
    """
    slabel is a unique smiles string for labeling species in QTC.
    Composed of two parts 'canonical smiles' and 'multiplicity'
    Canical smiles strings are unique only for a given code.
    QTC uses open babel.
    slabel = s + '_m' + str(mult)
    """
    if '_m' in s:
        s, mult = s.split('_m')
    s = get_smiles(s)
    if not mult:
        mult = get_multiplicity(s)
    return s + '_m' + str(mult)


def get_isomers(s):
    """
    Return a list of stereoisomers for a given smiles, slabel, xyz or mol object.
    Currently only works for cis/trans isomerism with one double bond.
    Generated isomers include chirality but only one of them.
    Depends on xyz generated by open_babel.
    Note:
    For multiple double bonds algorithm is complicated when there is branching.
    """
    import logging
    xyz = get_xyz(s)
    s2  = get_smiles(xyz)
    ndouble  = s2.count('=')
    nslashes = s2.count('/') + s2.count('\\')
    nchiral  = s2.count('@')
    mult = None
    if '_m' in s:
        mult = get_mult(s)
        s2 = get_slabel(s2,mult)
    isomers = [s2]
    if ndouble > 1:
        logging.debug('{0} double bonds in {1}'.format(ndouble,s))
        logging.debug('Can only work with one double bond')
    elif ndouble == 1 and nslashes > 1:
        logging.debug('One double bond and a stereocenter found in {0}'.format(s))
        if nslashes > 3:
            logging.debug('More than 3 slashes {} --> {}'.format(s, s2))
        left,right = s2.split('=')
        newright = ''
        for char in right:
            if char == '/':
                newchar = '\\'
            elif char == '\\':
                newchar = '/'
            else:
                newchar = char
            newright += newchar
        news = left + '=' + newright
        if mult:
            slabel = get_slabel(news,mult)
            news = get_slabel(news,mult)
        isomers.append(news)
    if nchiral > 0:
        print('{0} chiral centers in {1}'.format(nchiral,s))
    return isomers


def write_isomers_list(listfile):
    """
    Writes a file containing the isomers of the species in a given list.
    The filename for the list is required as the input.
    The filename of the new file is returned.
    """
    import iotools as io
    slist = io.read_list(listfile)
    newlist = ''
    for s in slist:
        isomers = get_isomers(s)
        for isomer in isomers:
            newlist += '{}\n'.format(isomer)
    newfilename = listfile.split('.')[0] + '_isomers.txt'
    io.write_file(newlist,filename=newfilename)
    return newfilename


def get_formula(x, hydrogens=True, stoichemetry=True):
    """
    Returns the molecular formula.
    I think it is based on Hill system.
    Print first carbon atom stoichemetry then hydrogens,
    and all the other elements in an alphabetical order.
    Usage:
    >>> get_formula('CCC')
    'C3H8'
    >>> [get_formula(s) for s in ['C', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    ['CH4', 'C2H6O', 'O2', 'H2O']
    >>> get_formula('CCC',hydrogens=False)
    'C3'
    >>> get_formula('CCC',hydrogens=False,stoichemetry=False)
    'C'
    >>> get_formula('CCC',hydrogens=True,stoichemetry=False)
    'CH'
    """
    mol = get_mol(x)
    formula = mol.formula
    if stoichemetry:
        if hydrogens:
            s = formula
        else:
            n = len(formula)
            s = ''
            hdigit = False
            for i in range(n):
                if formula[i].isdigit():
                    if not hdigit:
                        s += formula[i]
                elif formula[i] == 'H':
                    hdigit = True
                else:
                    hdigit = False
                    s += formula[i]
            if s == '':
                s = formula
    else:
        if hydrogens:
            s = ''.join(i for i in formula if not i.isdigit())
        else:
            s = ''.join(i for i in formula if not i.isdigit() and not i == 'H')
            if s == '':
                s = 'H'
    return s


def get_natom(x):
    """
    Return number of atoms in mol.
    >>> mol = get_mol('CC')
    >>> get_natom(mol)
    8
    >>> mols = [get_mol(s) for s in ['[CH3]', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_natom(mol) for mol in mols]
    [4, 9, 2, 3]
    """
    mol = get_mol(x,make3D=True)
    return len(mol.atoms)


def get_natom_heavy(x):
    """
    Return number of heavy atoms (nonHydrogens) in mol.
    >>> mol = get_mol('CC')
    >>> get_natom(mol)
    8
    >>> mols = [get_mol(s) for s in ['[CH3]', 'CCO', '[O][O]', 'InChI=1S/H2O/h1H2']]
    >>> [get_natom_heavy(mol) for mol in mols]
    [1, 3, 2, 1]
    """
    mol = get_mol(x,make3D=True)
    return mol.OBMol.NumHvyAtoms()


def get_nrotor(x):
    """
    Return number of rotors.
    """
    if type(x) == str:
        mol = get_mol(x,make3D=True)
    else:
        mol = get_mol(x)
    return mol.OBMol.NumRotors()

def get_nelectron(x):
    """
    Return number of electrons.
    >>> get_nelectron('C')
    10
    """
    mol = get_mol(x)
    n = 0
    for i in range(get_natom(mol)):
        a = mol.OBMol.GetAtomById(i)
        n += a.GetAtomicNum()
    return n


def get_charge(x):
    """
    Return charge.
    TODO
    """
    mol = get_mol(x, make3D=True)
    return mol.OBMol.GetTotalCharge()


def get_xyz(x):
    """
    Returns coordinates as a string in xyz format.
    Note: xyz coordinates are not deterministic.
    Each run gives a different set of coordinates.
    >>> mol = get_mol('CCCC')
    >>> print(get_xyz(mol).splitlines()[0])
    14
    """
    mol = get_mol(x, make3D=True)
    return mol.write(format='xyz')


def get_molpro_mol(logfile):
    """
    Returns xyz file from molpro logfile.
    """
    import pybel
    return next(pybel.readfile('mpo',logfile))


def get_gaussian_mol(logfile):
    """
    Returns mol file from gaussian logfile.
    """
    import pybel
    return next(pybel.readfile('g09',logfile))


def get_geo(x):
    """
    Returns coordinates-only as a string.
    Note: coordinates are not deterministic.
    Each run gives a different set of coordinates.
    """
    mol = get_mol(x, make3D=True)
    xyz = mol.write(format='xyz').splitlines(True)
    natom = int(xyz[0].strip())
    return ''.join(xyz[2:natom+2])


def get_zmat(x):
    """
    Returns internal coordinates as as string suitable for Gaussian zmat input.
    Note: zmat coordinates are deterministic.
    >>> print(get_zmat('C'))
    C
    H  1  r2
    H  1  r3  2  a3
    H  1  r4  2  a4  3  d4
    H  1  r5  2  a5  3  d5
    Variables:
    r2= 1.0922
    r3= 1.0922
    a3= 109.47
    r4= 1.0922
    a4= 109.47
    d4= 240.00
    r5= 1.0922
    a5= 109.47
    d5= 120.00
    <BLANKLINE>
    """
    mol = get_mol(x, make3D=True)
    return '\n'.join(mol.write('gzmat').splitlines()[5:])


def get_mop(x, keys='pm3 precise nosym threads=1 opt'):
    """
    Returns mopac input as a string.
    Note: For doctest I had to escape newline characters \n as \\n
    Since it gives EOL error.
    >>> xyz = "2\\n \\n H 0. 0. 0.\\n H 0. 0. 0.9\\n  \\n"
    >>> print(get_mop(xyz))
    pm3 precise nosym threads=1 opt
    <BLANKLINE>
    <BLANKLINE>
    H   0.00000 1  0.00000 1  0.00000 1
    H   0.00000 1  0.00000 1  0.90000 1
    <BLANKLINE>
    """
    mol = get_mol(x)
    return mol.write(format='mop', opt={'k': keys})


def get_inchi(x):
    """
    Returns a unique key composed of inchikey and multiplicity
    >>> mol = get_mol('[O][O]')
    >>> get_inchi_key(mol)
    'MYMOFIZGZYHOMD-UHFFFAOYSA-N3'
    """
    mol = get_mol(x)
    return mol.write("inchi").strip()


def get_inchi_key(x, mult=0, extra=''):
    """
    Returns a unique key composed of inchikey and multiplicity
    >>> mol = get_mol('[O][O]')
    >>> get_inchi_key(mol)
    'MYMOFIZGZYHOMD-UHFFFAOYSA-N3'
    """
    mol = get_mol(x)
    if mult == 0:
        mult = mol.spin
    return mol.write("inchikey").strip() + str(mult) + extra


def get_unique_name(x, mult=0, extra=''):
    """
    Returns a unique key composed of inchikey and multiplicity
    >>> mol = get_mol('[O][O]')
    >>> get_unique_name(mol)
    'MYMOFIZGZYHOMD-UHFFFAOYSA-N3'
    """
    mol = get_mol(x, make3D=True)
    if mult == 0:
        mult = mol.spin
    return mol.write("inchikey").strip() + str(mult) + extra


def get_unique_path(x, mult=0, method=''):
    """
    Returns a portable unique path based on inchikey for database directory.
    >>> import os
    >>> if os.path.sep == '/': print(get_unique_path('C',method='pm6'))
    database/C/C/CH4/VNWKTOKETHGBQD-UHFFFAOYSA-N1/pm6
    """
    import iotools as io
    mol = get_mol(x, make3D=True)
    if mult == 0:
        mult = mol.spin
    formula = get_formula(mol)
    formula_noH = get_formula(mol, stoichemetry=True, hydrogens=False)
    elements_noH = get_formula(mol, stoichemetry=False, hydrogens=False)
    uniquekey = get_inchi_key(mol, mult)
    dirs = 'database', elements_noH, formula_noH, formula, uniquekey, method
    return io.join_path(*dirs)


def get_formats():
    """
    Return available write formats in Open Babel.
    >>> for k, v in get_formats().items():
    ...     if 'gaussian z-matrix' in v.lower():
    ...         print(k)
    ...
    ...
    gzmat
    """
    return pybel.outformats


def get_smiles_mult(slabel):
    """
    Splits slabel into smiles and multiplicity and returns them as
    a string and an integer, respectively.
    """
    smi = slabel
    mult = 0
    if '_m' in smi:
        smi, mult = slabel.split('_m')
    if not mult:
        mult = get_multiplicity(smi)
    return smi, int(mult)


def get_smiles_path(x, mult=0, db= 'database'):
    """
    Returns a smiles based path for database directory.
    Note that smiles strings are not unique. Even the
    canonical smiles strings are unique only for the same
    code that generates the smiles string.
    """
    import iotools as io
    if type(x) is pybel.Molecule:
        if mult == 0:
            mult = x.spin
        s = x.write(format='can').strip().split()[0]
        s = s + '_m' + str(mult)
    elif type(x) is str:
        if '_m' in x:
            s = x
        else:
            s = get_smiles(x)
            mult = get_mult(s)
            s = s + '_m' + str(mult)
    formula = get_formula(s)
#    formula_noH = get_formula(s, stoichemetry=True, hydrogens=False)
#    elements_noH = get_formula(s, stoichemetry=False, hydrogens=False)
    s = get_smiles_filename(s)
#    dirs = db, elements_noH, formula_noH, formula, s
    dirs = db, formula, s
    return io.join_path(*dirs)


def get_smiles(x):
    """
    Returns open-babel canonical smiles.
    >>> print(get_smiles('O'))
    O
    >>> print(get_smiles('[H][O][H]'))
    O
    >>> print(get_smiles('O-O'))
    OO
    >>> print(get_smiles('[O]=[O]'))
    O=O
    """
    mol = get_mol(x)
    s = mol.write(format='can').strip().split()[0]
    return s


def get_smiles_filename(x):
    """
    Returns a suitable filename for a given smiles.
    Smiles strings may contain characters not suitable for file names,
    such as \/:*?"<>|(). Not sure if all these characters appear, but here
    they are replaced by an underscore, '_' followed by:
    """
    if type(x) is pybel.Molecule:
        s = x.write(format='can').strip().split()[0]
    elif type(x) is str:
        s = x
    else:
        s = ''
    s = s.replace('[','_b_')
    s = s.replace(']','_d_')
    #s = s.replace('=','_e_')
    s = s.replace(':','_i_')
    s = s.replace('|','_j_')
    s = s.replace('\\','_k_')
    s = s.replace('/','_l_')
    s = s.replace('(','_p_')
    s = s.replace(')','_q_')
    s = s.replace('*','_s_')
    #s = s.replace('#','_x_')
    s = s.replace('<','_v_')
    s = s.replace('>','_y_')
    s = s.replace('?','_z_')

    return s


def smiles2formula(filename):
    import iotools as io
    mols = io.read_list(filename)
    s = ''
    for mol in mols:
        formula = get_formula(mol)
        s += '{0} {1}\n'.format(mol,  formula)
    return s

def get_coordinates_array(xyz):
    """
    Given xyz string, return natom*3 array that contains the
    coordinates of the atoms.
    """
    lines = xyz.splitlines()
    n = int(lines[0])
    coords = [0.]*(n*3)
    i = 0
    for line in lines[2:2+n]:
        coords[i:i+3] = [float(c) for c in line.split()[1:4]]
        i += 3
    return coords


def set_mult(x,mult):
    """
    Sets the total spin multiplicity.
    """
    mol = get_mol(x)
    mult = int(mult)
    assert mult > 0, 'Multiplicity should be a positve integer'
    mol.OBMol.SetTotalSpinMultiplicity(mult)
    return mol


def set_xyz(x,coords):
    """
    Parameters:
    mol : Open babel mol object, or anything that can be converted to mol object with get_mol.
    coords : One-d array of natom floats
    """
    mol = get_mol(x)
    mol.OBMol.SetCoordinates(openbabel.double_array(coords))
    return mol


def fetch_smiles(s):
    """
    Returns the smiles string for a given chemical name.
    Requires cirpy module and internet connection
    >>> fetch_smiles('methane')
    'C'
    """
    try:
        import cirpy
    except:
        r = 'cirpy module not installed, see http://cirpy.readthedocs.io/'
        return
    if cirpy:
        return cirpy.resolve(s,'smiles')
    else:
        return None


def fetch_inchi(s):
    """
    Returns the smiles string for a given chemical name.
    Requires cirpy module and internet connection
    >>> fetch_inchi('methane')
    'InChI=1/CH4/h1H4'
    """
    try:
        import cirpy
    except:
        r = 'cirpy module not installed, see http://cirpy.readthedocs.io/'
        return
    if cirpy:
        r = cirpy.resolve(s,'inchi')
    return r


def fetch_IUPAC_name(s):
    """
    Return IUPAC name for a given smiles or inchi string.
    Requires cirpy module and internet connection
    >>> print(fetch_IUPAC_name('C=O'))
    FORMALDEHYDE
    """
    try:
        import cirpy
    except:
        r = 'cirpy module not installed, see http://cirpy.readthedocs.io/'
        return
    frm = get_format(s)
    if frm == 'smi':
        name = cirpy.resolve(s,'iupac_name',resolvers=['smiles'])
    elif frm == 'inchi':
        name = cirpy.resolve(s,'iupac_name',resolvers=['inchi'])
    elif frm == 'xyz':
        mol = get_mol(s)
        name = cirpy.resolve(mol.write('inchi').strip(),'iupac_name',resolvers=['inchi'])
    else:
        name = None
    return name


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

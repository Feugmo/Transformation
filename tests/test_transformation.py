from transformation import __version__
from transformation.align import transform
from ase.io import read
import numpy as np
import os

from transformation.tools import Molecule

from pathlib import Path
from ase import Atoms

path = Path(__file__).parent.absolute()


# def test_version():
#     assert __version__ == '0.2.0'


def test_reorient():
    assert os.path.exists(os.path.join(path, 'data/test_HEA_83_Transformed.POSCAR'))
    atms = transform(param_file=os.path.join(path, 'params.yaml'), writeoutput=True)
    test_atom = read(os.path.join(path, 'data/test_HEA_83_Transformed.POSCAR'))
    assert test_atom.get_chemical_symbols() == atms.get_chemical_symbols()


def extract_box():
    atms = read(os.path.join(path, 'data/test_HEA_83_Transformed.POSCAR'))
    box = np.zeros((8, 3))
    a = 3.552
    x = 5 / 4
    y = 3 / 4
    z = 3 / 4
    box[0] = np.array([x, -y, -z]) * a  # A
    box[1] = np.array([x, y, -z]) * a  # B
    box[2] = np.array([-x, y, -z]) * a  # C
    box[3] = np.array([-x, -y, -z]) * a  # D
    box[4] = np.array([x, -y, z]) * a  # E
    box[5] = np.array([x, y, z]) * a  # F
    box[6] = np.array([-x, y, z]) * a  # G
    box[7] = np.array([-x, -y, z]) * a  # H

    inside, cell = Molecule.in_box2(atms.positions, box)
    symbols = [atms[i].symbol for i in inside]
    coords = [atms[i].position for i in inside]

    transformed_atoms = Atoms(symbols, positions=coords, cell=cell)
    box = read(os.path.join(path, 'data/sliced_atoms.POSCAR'))
    assert transformed_atoms.get_chemical_symbols() == box.get_chemical_symbols()
#
#

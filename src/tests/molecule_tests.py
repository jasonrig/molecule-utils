import unittest
import os
from io import StringIO

from molutils.util.molecule import Molecule
from molutils.util.periodic_table import lookup_element_by_symbol

PATH_TO_PSI4 = "/opt/psi4/bin/psi4.run"

DIMER_XYZ_FILE = (
    "21\n"
    "molecule_title\n"
    "C	0.00000000	0.00000000	0.00000000\n"
    "N	0.00000000	0.00000000	1.37405500\n"
    "C	1.26815400	0.00000000	1.80843900\n"
    "N	2.07877900	-0.02792200	0.74098800\n"
    "C	1.31304500	-0.01765900	-0.39985000\n"
    "C	-1.17649400	0.19094100	2.22511900\n"
    "C	3.53427200	0.12780900	0.79070400\n"
    "F	1.09782400	2.71748500	1.13271500\n"
    "B	1.51343800	2.90805500	2.48894200\n"
    "F	0.55290600	2.25488800	3.33063400\n"
    "F	1.62575900	4.25335000	2.80352900\n"
    "F	2.76311600	2.22460800	2.65697000\n"
    "H	1.58129100	0.09619200	2.83199800\n"
    "H	1.75047000	-0.01157500	-1.38238800\n"
    "H	-0.91045500	0.02428600	-0.57204000\n"
    "H	-1.96615300	-0.46999100	1.87769000\n"
    "H	-0.90427100	-0.04968500	3.24681200\n"
    "H	-1.47152400	1.23473700	2.17603800\n"
    "H	3.87180400	-0.11397200	1.79244500\n"
    "H	3.97984900	-0.54952400	0.06697900\n"
    "H	3.77751100	1.16445100	0.57797700\n")

WATER_XYZ_FILE = (
    "3\n"
    "\n"
    "O                    0.008165   -0.658276    0.215737\n"
    "H                    1.508729    0.441434    0.259102\n"
    "H                   -1.419465    0.543189    0.243811\n"
)


class MoleculeTest(unittest.TestCase):
    @staticmethod
    def get_sorted_atom_string(molecule):
        return ''.join(sorted([atom[0] for atom in molecule]))

    def test_read(self):
        molecule = Molecule.from_xyz_file(StringIO(DIMER_XYZ_FILE))
        self.assertEqual(len(molecule.atom_list), 21)
        self.assertEqual(len(molecule), 21)

        self.assertEqual(self.get_sorted_atom_string(molecule), 'BCCCCCFFFFHHHHHHHHHNN')

    def test_fragment(self):
        molecule = Molecule.from_xyz_file(StringIO(DIMER_XYZ_FILE))
        # Create two fragments and sort them by number of atoms
        fragments = sorted(molecule.fragment(2), key=len)

        # Check that two fragments are created
        self.assertEqual(len(fragments), 2)

        # The first fragment should have five atoms, BFFF
        self.assertEqual(len(fragments[0]), 5)
        self.assertEqual(self.get_sorted_atom_string(fragments[0]), 'BFFFF')

        # The second fragment should have 16 atoms, CCCCCHHHHHHHHHNN
        self.assertEqual(len(fragments[1]), 16)
        self.assertEqual(self.get_sorted_atom_string(fragments[1]), 'CCCCCHHHHHHHHHNN')

    def test_z_sum(self):
        molecule = Molecule.from_xyz_file(StringIO(WATER_XYZ_FILE))
        self.assertEqual(molecule.get_z_sum(), 10)

    def test_lookup_atom_z(self):
        self.assertEqual(lookup_element_by_symbol('H')[0], 1)
        self.assertEqual(lookup_element_by_symbol('h')[0], 1)
        self.assertIsNone(lookup_element_by_symbol('zz'))

    def test_possible_closed_shell_charges(self):
        ionic_dimer = Molecule.from_xyz_file(StringIO(DIMER_XYZ_FILE))
        fragments = sorted(ionic_dimer.fragment(2), key=len)
        self.assertEqual(ionic_dimer.get_possible_closed_shell_charges(), [0])
        self.assertEqual(fragments[0].get_possible_closed_shell_charges(), [-1, 1])
        self.assertEqual(fragments[1].get_possible_closed_shell_charges(), [-1, 1])

        water = Molecule.from_xyz_file(StringIO(WATER_XYZ_FILE))
        self.assertEqual(water.get_possible_closed_shell_charges(), [0])

    def test_guess_charge(self):
        if not os.path.isfile(PATH_TO_PSI4) or not os.access(PATH_TO_PSI4, os.X_OK):
            self.skipTest("No Psi4 executable found.")
        molecule = Molecule.from_xyz_file(StringIO(DIMER_XYZ_FILE), PATH_TO_PSI4)
        fragments = sorted(molecule.fragment(2), key=len)
        self.assertEqual(molecule.guess_charge(), 0)
        self.assertEqual(fragments[0].guess_charge(), -1)
        self.assertEqual(fragments[1].guess_charge(), 1)

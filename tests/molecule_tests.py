import unittest
import os
from io import StringIO

from molutils.util.molecule import Molecule
from molutils.util.periodic_table import lookup_element_by_symbol
from molutils.util.job_formatters.psi4 import Psi4JobFormatter

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

NITROGEN_ATOM = (
    "1\n"
    "\n"
    "N  0.00 0.00 0.00"
)

PSI4_OUTPUT = (
    "  ==> Geometry <==\n"
    "\n"
    "    Molecular point group: c1\n"
    "    Full point group: C1\n"
    "\n"
    "    Geometry (in Angstrom), charge = 0, multiplicity = 1:\n"
    "\n"
    "       Center              X                  Y                   Z               Mass\n"
    "    ------------   -----------------  -----------------  -----------------  -----------------\n"
    "           F         -2.193629820047    -1.855516068908    -1.383582467737    18.998403224000\n"
    "           C         -1.172015920047    -1.019644636908    -1.218746331737    12.000000000000\n"
    "           C          0.114669063953    -1.517962241908    -1.059929846737    12.000000000000\n"
    "           F          0.325520175953    -2.831149664908    -1.070726500737    18.998403224000\n"
    "           F         -2.628991937047     0.828237203092    -1.313676415737    18.998403224000\n"
    "           C         -1.394261984047     0.350968436092    -1.184005598737    12.000000000000\n"
    "           C         -0.330144528047     1.222797545092    -0.992235587737    12.000000000000\n"
    "           F         -0.545459093047     2.533419542092    -0.947497268737    18.998403224000\n"
    "           C          1.178862928953    -0.646155231908    -0.868311389737    12.000000000000\n"
    "           F          2.407922502953    -1.124703761908    -0.705774696737    18.998403224000\n"
    "           C          0.957567807953     0.724990373092    -0.837317381737    12.000000000000\n"
    "           F          1.975930935953     1.559110990092    -0.661179892737    18.998403224000\n"
    "           H         -1.644566974047    -1.210383778908     2.174557304263     1.007825032070\n"
    "           C         -0.812314959047    -0.527350407908     2.281532461263    12.000000000000\n"
    "           H         -2.042421771047     1.234253206092     2.190249710263     1.007825032070\n"
    "           C         -1.036491811047     0.849083608092     2.292233299263    12.000000000000\n"
    "           C          0.484176730953    -1.023293069908     2.416217104263    12.000000000000\n"
    "           H          0.658329572953    -2.091089033908     2.410131913263     1.007825032070\n"
    "           C          0.035280831953     1.729129498092     2.439058175263    12.000000000000\n"
    "           H         -0.138951482047     2.796453733092     2.446813182263     1.007825032070\n"
    "           C          1.331069796953     1.232787556092     2.574283147263    12.000000000000\n"
    "           H          2.162281353953     1.915280645092     2.685863082263     1.007825032070\n"
    "           C          1.555807218953    -0.143098610908     2.562921850263    12.000000000000\n"
    "           H          2.561079041953    -0.528098246908     2.666538093263     1.007825032070\n"
    "\n"
    "  Running in c1 symmetry.\n"
    "\n"
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

        molecule = Molecule.from_psi4_output(StringIO(PSI4_OUTPUT))
        self.assertEqual(len(molecule.atom_list), 24)
        self.assertEqual(len(molecule), 24)
        self.assertEqual(self.get_sorted_atom_string(molecule), 'CCCCCCCCCCCCFFFFFFHHHHHH')

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

    def test_possible_charges(self):
        ionic_dimer = Molecule.from_xyz_file(StringIO(DIMER_XYZ_FILE))
        fragments = sorted(ionic_dimer.fragment(2), key=len)
        self.assertEqual(ionic_dimer.get_possible_charges(), [0])
        self.assertEqual(fragments[0].get_possible_charges(), [-1, 1])
        self.assertEqual(fragments[1].get_possible_charges(), [-1, 1])

        water = Molecule.from_xyz_file(StringIO(WATER_XYZ_FILE))
        self.assertEqual(water.get_possible_charges(), [0])

        nitrogen_atom = Molecule.from_xyz_file(StringIO(NITROGEN_ATOM))
        self.assertEqual(nitrogen_atom.get_possible_charges(multiplicity=4), [0])
        self.assertEqual(nitrogen_atom.get_possible_charges(multiplicity=1), [-1, 1])

    def test_guess_charge(self):
        if not os.path.isfile(PATH_TO_PSI4) or not os.access(PATH_TO_PSI4, os.X_OK):
            self.skipTest("No Psi4 executable found.")
        molecule = Molecule.from_xyz_file(StringIO(DIMER_XYZ_FILE), PATH_TO_PSI4)
        fragments = sorted(molecule.fragment(2), key=len)
        self.assertEqual(molecule.guess_charge(), 0)
        self.assertEqual(fragments[0].guess_charge(), -1)
        self.assertEqual(fragments[1].guess_charge(), 1)

        nitrogen_atom = Molecule.from_xyz_file(StringIO(NITROGEN_ATOM), PATH_TO_PSI4)
        self.assertEqual(nitrogen_atom.guess_charge(multiplicity=4), 0)
        self.assertEqual(nitrogen_atom.guess_charge(multiplicity=1), -1)

    def test_efp_atom_padding(self):
        molecule = Molecule('', [('Cl', 0.0, 0.0, 0.0)], 1, 1)
        molecule.efp_pad_dummy_atoms()
        self.assertEqual(len(molecule), 3)

        molecule = Molecule('', [('Cl', 0.0, 0.0, 0.0), ('Cl', 0.5, 0.0, 0.0),], 0, 1)
        molecule.efp_pad_dummy_atoms()
        self.assertEqual(len(molecule), 3)

    def test_sapt_job_format(self):
        expected_output = (
            "memory 250 mb\n"
            "molecule molecule_title2 {\n"
            "-1 1\n"
            "  F 1.0978240000 2.7174850000 1.1327150000\n"
            "  B 1.5134380000 2.9080550000 2.4889420000\n"
            "  F 1.6257590000 4.2533500000 2.8035290000\n"
            "  F 2.7631160000 2.2246080000 2.6569700000\n"
            "  F 0.5529060000 2.2548880000 3.3306340000\n"
            "--\n"
            "1 1\n"
            "  C 0.0000000000 0.0000000000 0.0000000000\n"
            "  H -0.9104550000 0.0242860000 -0.5720400000\n"
            "  C 1.3130450000 -0.0176590000 -0.3998500000\n"
            "  H 1.7504700000 -0.0115750000 -1.3823880000\n"
            "  N 0.0000000000 0.0000000000 1.3740550000\n"
            "  C 1.2681540000 0.0000000000 1.8084390000\n"
            "  H 1.5812910000 0.0961920000 2.8319980000\n"
            "  N 2.0787790000 -0.0279220000 0.7409880000\n"
            "  C -1.1764940000 0.1909410000 2.2251190000\n"
            "  H -0.9042710000 -0.0496850000 3.2468120000\n"
            "  H -1.4715240000 1.2347370000 2.1760380000\n"
            "  H -1.9661530000 -0.4699910000 1.8776900000\n"
            "  C 3.5342720000 0.1278090000 0.7907040000\n"
            "  H 3.8718040000 -0.1139720000 1.7924450000\n"
            "  H 3.7775110000 1.1644510000 0.5779770000\n"
            "  H 3.9798490000 -0.5495240000 0.0669790000\n"
            "}\n"
            "\n"
            "set {\n"
            "  guess sad\n"
            "  basis_guess 3-21G\n"
            "  basis cc-pVTZ\n"
            "  scf_type DF\n"
            "  freeze_core True\n"
            "}\n"
            "energy('sapt0')\n"
        )
        molecule = Molecule.from_xyz_file(StringIO(DIMER_XYZ_FILE), PATH_TO_PSI4)
        fragments = sorted(molecule.fragment(2), key=len)
        fragments[0].charge = -1
        fragments[0].multiplicity = 1
        fragments[1].charge = 1
        fragments[1].multiplicity = 1
        job = Psi4JobFormatter(fragments).energy("sapt0")
        self.assertEqual(job, expected_output)

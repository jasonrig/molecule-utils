import math
import os
import re
import subprocess

from .molecule_formatters import MoleculeFormatterMixin
from .periodic_table import lookup_element_by_symbol

DEFAULT_PSI4_EXECUTABLE = "psi4"


class Molecule(MoleculeFormatterMixin):
    def __init__(self, title, atom_list=None, charge=0, multiplicity=1, psi4_path=None):

        if psi4_path is not None:
            self.psi4_path = psi4_path
        elif os.path.isfile(DEFAULT_PSI4_EXECUTABLE) and os.access(DEFAULT_PSI4_EXECUTABLE, os.X_OK):
            self.psi4_path = DEFAULT_PSI4_EXECUTABLE
        else:
            self.psi4_path = None

        self.charge = charge
        self.multiplicity = multiplicity

        self.title = title.strip()
        if not self.title:
            self.title = "input_molecule"

        if atom_list is not None:
            self.atom_list = list(atom_list)
        else:
            self.atom_list = []

    def __iter__(self):
        return self.atom_list.__iter__()

    def __len__(self):
        return self.atom_list.__len__()

    @staticmethod
    def from_xyz_file(xyz_file, psi4_path=None):
        """
        Creates a Molecule object from an XYZ file
        :param xyz_file: the xyz file as a file-like object or path
        """

        def read_xyz_file(fp):
            title = ''
            atom_list = []
            line_number = 0
            for line in fp:
                # Ignore leading blank lines
                if line_number == 0 and len(line.strip()) == 0:
                    continue
                xyz_parts = line.split()
                if line_number == 1:
                    title = line
                elif line_number > 1 and len(xyz_parts) >= 4:
                    atom_list.append((xyz_parts[0], float(xyz_parts[1]), float(xyz_parts[2]), float(xyz_parts[3])))
                line_number += 1
            return Molecule(title, atom_list, psi4_path=psi4_path)

        if isinstance(xyz_file, str):
            with open(xyz_file, 'r') as f:
                return read_xyz_file(f)
        else:
            return read_xyz_file(xyz_file)

    def add_atom(self, label, x, y, z):
        """
        Adds an atom to the molecule
        :param label: the atom label
        :param x: x coordinate
        :param y: y coordinate
        :param z: z coordinate
        """
        self.atom_list.append((label, float(x), float(y), float(z)))

    def merge(self, molecule):
        """
        Merges all atoms from another molecule into this one
        :param molecule: a second molecule
        """
        for atom in molecule:
            self.add_atom(*atom)

    def distance_from(self, molecule):
        """
        Calculates the distance between the two nearest atoms of this and another molecule
        :param molecule: a second molecule
        """
        distance = None
        for a1 in self:
            for a2 in molecule:
                d = math.sqrt(math.pow(a1[1] - a2[1], 2) + math.pow(a1[2] - a2[2], 2) + math.pow(a1[3] - a2[3], 2))
                if distance is None:
                    distance = d
                elif d < distance:
                    distance = d
        return distance

    def fragment(self, n_frags):
        """
        Fragments a molecule based on nearest neighbor classification
        :param n_frags: number of fragments
        """

        # start off with each atom as a separate molecule
        fragments = [Molecule(self.title, [atom], psi4_path=self.psi4_path) for atom in self]

        while len(fragments) > n_frags:
            distance = None
            i = None
            j = None
            for _i in range(len(fragments)):
                for _j in range(_i + 1, len(fragments)):
                    d = fragments[_i].distance_from(fragments[_j])
                    if distance is None:
                        distance = d
                        i = _i
                        j = _j
                    elif d < distance:
                        distance = d
                        i = _i
                        j = _j
            fragments[i].merge(fragments[j])
            del fragments[j]

        for i in range(n_frags):
            fragments[i].title += str(n_frags)

        return fragments

    def get_z_sum(self):
        return sum([lookup_element_by_symbol(atom[0])[0] for atom in self])

    def electron_count(self):
        return self.get_z_sum() + self.charge

    def get_possible_charges(self, lower_range=-1, upper_range=1, multiplicity=1):
        if multiplicity % 2 > 0:
            return [q for q in range(lower_range, upper_range + 1) if (self.get_z_sum() + q) % 2 == 0]
        else:
            return [q for q in range(lower_range, upper_range + 1) if (self.get_z_sum() + q) % 2 > 0]

    def guess_charge(self, lower_range=-1, upper_range=1, multiplicity=1):
        if self.psi4_path is None:
            raise ValueError("Psi4 path must be provided for this method to work")

        possible_charges = self.get_possible_charges(lower_range=lower_range, upper_range=upper_range, multiplicity=multiplicity)
        if len(possible_charges) == 1:
            return possible_charges[0]

        job_template = (
            "{molecule}"
            "set basis STO-3G\n"
            "set reference {reference}\n"
            "set guess sad\n"
            "energy('scf')"
        )

        p = re.compile(r'Total Energy\s+=\s+([-0-9\.]+)')

        self.multiplicity = multiplicity
        lowest_energy_and_charge = None
        for q in possible_charges:
            proc = subprocess.Popen([self.psi4_path, '-i', 'stdin', '-o', 'stdout'], stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE)
            self.charge = q
            proc.stdin.write(
                str.encode(
                    job_template.format(
                        molecule=self.format_psi4(),
                        reference='rhf' if multiplicity == 1 else 'uhf'
                    )
                )
            )
            proc.stdin.close()
            result = proc.stdout.read().decode('utf-8')
            print(result)
            energy_search = p.search(result)
            if energy_search is not None:
                energy = float(energy_search.group(1))
                if lowest_energy_and_charge is None or energy < lowest_energy_and_charge[1]:
                    lowest_energy_and_charge = (q, energy)
            proc.wait()
        if lowest_energy_and_charge:
            self.charge = lowest_energy_and_charge[0]
            return lowest_energy_and_charge[0]
        else:
            return None

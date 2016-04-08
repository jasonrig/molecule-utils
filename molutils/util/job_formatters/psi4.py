from ..molecule import Molecule


class Psi4JobFormatter(object):

    def __init__(self, molecule, basis_set="cc-pVTZ", memory=250, memory_units="mb"):
        self.molecule = molecule
        self.basis_set = basis_set
        self.memory = memory
        self.memory_units = memory_units

    def energy(self, type="scf", molecule=None):
        if molecule is None:
            molecule = self.molecule

        if isinstance(molecule, Molecule):
            molecule_text = molecule.format_psi4()
        else:
            molecule_text = Molecule.format_psi4_group(molecule)

        template = (
            "memory {memory} {memory_units}\n"
            "{molecule}\n"
            "set {{\n"
            "  guess sad\n"
            "  basis_guess 3-21G\n"
            "  basis {basis_set}\n"
            "  scf_type DF\n"
            "  freeze_core True\n"
            "}}\n"
            "energy('{type}')\n"
        )
        return template.format(
            memory=self.memory,
            memory_units=self.memory_units,
            molecule=molecule_text,
            basis_set=self.basis_set,
            type=type
        )
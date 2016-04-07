from ..molecule import Molecule


class Psi4JobFormatter(object):

    def __init__(self, molecule, basis_set="cc-pVTZ", memory="250 mb"):
        self.molecule = molecule
        self.basis_set = basis_set
        self.memory = memory

    def energy(self, type="scf", molecule=None):
        if molecule is None:
            molecule = self.molecule

        if isinstance(molecule, Molecule):
            molecule_text = molecule.format_psi4()
        else:
            molecule_text = Molecule.format_psi4_group(molecule)

        template = (
            "memory {memory}\n"
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
            molecule=molecule_text,
            basis_set=self.basis_set,
            type=type
        )
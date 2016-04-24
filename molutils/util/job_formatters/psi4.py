from ..molecule import Molecule
from .formatter import Formatter


class Psi4JobFormatter(Formatter):

    def __init__(self, molecule, basis_set="cc-pVTZ", memory=250, memory_units="mb"):
        self.molecule = molecule
        self.basis_set = basis_set
        self.memory = memory
        self.memory_units = memory_units

    def format(self, type, method, guess_charge=False):
        if type == "energy":
            if method is None:
                method = "mp2"
            return self.energy(type=method, guess_charge=guess_charge)
        else:
            raise NotImplemented("Calculation type %s not implemented for Psi4 calculations" % type)

    def energy(self, type="scf", guess_charge=False):
        if isinstance(self.molecule, Molecule):
            molecule_text = self.molecule.format_psi4(guess_charge=guess_charge)
        elif len(self.molecule) == 1:
            molecule_text = self.molecule[0].format_psi4(guess_charge=guess_charge)
        else:
            molecule_text = Molecule.format_psi4_group(self.molecule, guess_charge=guess_charge)

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
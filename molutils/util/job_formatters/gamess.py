from ..molecule import Molecule
from .formatter import Formatter


class GamessJobFormatter(Formatter):

    def __init__(self, molecule, basis_set=None, memory_replicated_gb=1, memory_distributed_gb=1):
        self.molecule = molecule
        self.basis_set = basis_set
        if self.basis_set is None:
            self.basis_set = "GBASIS=CCT"
        self.mwords = int(memory_replicated_gb * 125)
        self.memddi = int(memory_distributed_gb * 125)

    def _format_basis_set(self):
        return "$BASIS %s $END" % self.basis_set

    def _format_memory(self):
        return "$SYSTEM MWORDS=%i MEMDDI=%i $END" % (self.mwords, self.memddi)

    def format(self, type, method, guess_charge=False):
        if type == "makefp":
            return self.makefp(guess_charge=guess_charge)
        else:
            raise NotImplementedError("Calculation type %s not implemented for GAMESS calculations" % type)

    def makefp(self, guess_charge=False):
        if isinstance(self.molecule, Molecule):
            molecule_text = self.molecule.efp_pad_dummy_atoms().format_gamess(guess_charge=guess_charge)
        elif len(self.molecule) == 1:
            molecule_text = self.molecule[0].efp_pad_dummy_atoms().format_gamess(guess_charge=guess_charge)
        else:
            raise ValueError("makefp cannot format a molecule group")

        template = (
            " {memory}\n"
            " $CONTRL RUNTYP=MAKEFP $END\n"
            " {basis_set}\n"
            "{molecule}\n"
        )

        return template.format(memory=self._format_memory(),
                               basis_set=self._format_basis_set(),
                               molecule=molecule_text)
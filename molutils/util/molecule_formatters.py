from .periodic_table import lookup_element_by_symbol


class MoleculeFormatterMixin(object):
    def format_psi4(self, guess_charge=False):
        """
        Formats a molecule in a format suitable for psi4
        """

        if guess_charge:
            self.guess_charge()

        template = ("molecule {title} {{\n"
                    "{charge} {multiplicity}\n"
                    "{coordinates}"
                    "}}\n")
        coordinate_template = "  {label} {x:.10f} {y:.10f} {z:.10f}\n"

        coordinate_block = ""
        for atom in self:
            coordinate_block += coordinate_template.format(label=atom[0], x=atom[1], y=atom[2], z=atom[3])

        return template.format(title=self.title, charge=self.charge, multiplicity=self.multiplicity,
                               coordinates=coordinate_block)

    @staticmethod
    def format_psi4_group(molecule_list, guess_charge=False):

        for molecule in molecule_list:
            if guess_charge:
                molecule.guess_charge()

        coordinate_template = "  {label} {x:.10f} {y:.10f} {z:.10f}\n"
        molecule_inner_template = (
            "{charge} {multiplicity}\n"
            "{coordinates}"
        )
        divider = '--\n'
        molecule_outer_template = ("molecule {title} {{\n"
                          "{molecule_inner}"
                          "}}\n")

        molecule_inners = divider.join([molecule_inner_template.format(
            charge=m.charge,
            multiplicity=m.multiplicity,
            coordinates=''.join([coordinate_template.format(label=atom[0], x=atom[1], y=atom[2], z=atom[3]) for atom in m])
        ) for m in molecule_list])
        return molecule_outer_template.format(title=molecule_list[0].title, molecule_inner=molecule_inners)

    def format_gamess(self, guess_charge=False):
        """
        Formats a molecule in a format sutiable for GAMESS
        """

        if guess_charge:
            self.guess_charge()

        scf_type = "RHF"
        if self.multiplicity != 1:
            scf_type = "UHF"

        template = (" $CONTRL ICHARG={charge} MULT={multiplicity} SCFTYP={scf_type} $END\n"
                    " $CONTRL ISPHER=1 $END\n"
                    " $SCF DIRSCF=.TRUE. DIIS=.TRUE. $END\n"
                    " $DATA\n"
                    "C1\n"
                    "{title}\n"
                    "{coordinates}"
                    " $END\n")
        coordinate_template = "{label}   {z_number:.1f}   {x:.10f} {y:.10f} {z:.10f}\n"

        coordinate_block = ""
        for atom in self:
            coordinate_block += coordinate_template.format(label=atom[0], z_number=lookup_element_by_symbol(atom[0])[0],
                                                           x=atom[1], y=atom[2], z=atom[3])

        return template.format(title=self.title, charge=self.charge, multiplicity=self.multiplicity, scf_type=scf_type,
                               coordinates=coordinate_block)

#!/usr/bin/env python3
import argparse
import re
from molutils.util.molecule import Molecule
from molutils.util.job_formatters.psi4 import Psi4JobFormatter
from molutils.util.job_formatters.gamess import GamessJobFormatter


def main(args):
    for file in args.input:
        molecules = Molecule.from_file(file, psi4_path=args.path_to_psi4)
        if args.n_frags > 1:
            molecules = molecules.fragment(args.n_frags)
        else:
            molecules = [molecules]

        # Psi4 calcs
        if args.output_format.lower() == "psi4":
            job_formatter = Psi4JobFormatter(molecules, basis_set=args.basis_set, memory=args.memory, memory_units="Gb")
            _output(job_formatter.format(args.calc_type, args.calc_method, guess_charge=args.guess_charge),
                    file,
                    'inp',
                    args.output_to)

        # GAMESS calcs
        elif args.output_format.lower() == "gamess":
            i = 0
            for m in molecules:
                job_formatter = GamessJobFormatter(m, args.basis_set, args.memory, args.memory_ddi)
                _output(job_formatter.format(args.calc_type, args.calc_method, guess_charge=args.guess_charge),
                        "%i%s" % (i, file),
                        'inp',
                        args.output_to)

        else:
            raise NotImplemented("%s output format not yet implemented" % args.output_format)


def _output(content, input_file_name, output_ext, destination):
    file_name_parts = input_file_name.rsplit('.', 1)
    if len(file_name_parts > 1):
        ext_matcher = re.compile('\\.(%s)$' % file_name_parts[-1])
    else:
        ext_matcher = None
    if destination == "STDOUT":
        print(content)
        return
    elif destination == "AUTO":
        if ext_matcher:
            output_file_name = ext_matcher.sub(output_ext, input_file_name, count=1)
        else:
            output_file_name = "%s.%s" % (input_file_name, output_ext)
        print("Created: %s" % output_file_name)
    else:
        input_file_name = destination

    with open(input_file_name, "w") as f:
        f.write(content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the input XYZ files to process", nargs="+", type=str)
    parser.add_argument("--output_format", help="the software for which an input file should be generated", type=str,
                        choices=['psi4', 'gamess'], default="psi4")
    parser.add_argument("--output_to", help="file name for output to be written,"
                                            "'STDOUT' to print to screen, or 'AUTO' to autogenerate file names",
                        type=str, default="STDOUT")
    parser.add_argument("--calc_type", help="the type of calculation to request", type=str, choices=['energy'],
                        default="energy")
    parser.add_argument("--calc_method", help="the method of calculation (e.g. type=energy, method=mp2)", type=str,
                        default=None)
    parser.add_argument("--basis_set", help="the basis set to use", type=str, default="cc-pVTZ")
    parser.add_argument("--n_frags", help="the number of fragments the XYZ file should be split into", type=int,
                        default=1)
    parser.add_argument("--guess_charge", help="guess the charge", action="store_true", default=False)
    parser.add_argument("--memory", help="memory to use in calculation in GB", type=int, default=1)
    parser.add_argument("--memory_ddi", help="distributed memory to use in GAMESS calculations in GB", type=int,
                        default=1)
    parser.add_argument("--path_to_psi4", help="path to the psi4 executable", default="psi4")
    args = parser.parse_args()
    main(args)

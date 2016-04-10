import argparse
import re
from molutils.util.molecule import Molecule
from molutils.util.job_formatters.psi4 import Psi4JobFormatter


def main(args):
    xyz_ext_matcher = re.compile('\\.(xyz)$')
    if args.output_format.lower() == "psi4":
        input_files = []
        for file in args.input:
            molecules = Molecule.from_xyz_file(file, psi4_path=args.path_to_psi4)

            if args.n_frags > 1:
                molecules = molecules.fragment(args.n_frags)
            else:
                molecules = [molecules]

            if args.guess_charge:
                for m in molecules:
                    m.guess_charge()
            job_formatter = Psi4JobFormatter(molecules, basis_set=args.basis_set, memory=args.memory, memory_units="Gb")
            if args.calc_type == "energy":
                output = job_formatter.energy(type=args.calc_method, molecule=molecules)
            else:
                raise NotImplemented("%s calculation type not implemented" % args.calc_type)
            if args.output_to == "STDOUT":
                print(output)
                continue
            elif args.output_to == "AUTO":
                input_file_name = xyz_ext_matcher.sub('.inp', file, count=1)
                print("Created: %s" % input_file_name)
                if not input_file_name.endswith('.inp'):
                    input_file_name += '.inp'
            else:
                input_file_name = args.output_to

            output_file_name = xyz_ext_matcher.sub('.out', file, count=1)
            if not output_file_name.endswith('.out'):
                output_file_name += '.out'
            input_files.append((input_file_name, output_file_name))

            with open(input_file_name, "w") as f:
                f.write(output)

        run_script = ["#!/bin/bash"]
        for input, output in input_files:
            run_script.append("echo \"Running %s...\"" % input)
            run_script.append(args.path_to_psi4 + " " + input + " " + output)
        run_script = '\n'.join(run_script)
        with open('run.sh', 'w') as f:
            f.write(run_script)
    else:
        raise NotImplemented("%s output format not yet implemented" % args.output_format)


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
                        default="mp2")
    parser.add_argument("--basis_set", help="the basis set to use", type=str, default="cc-pVTZ")
    parser.add_argument("--n_frags", help="the number of fragments the XYZ file should be split into", type=int,
                        default=1)
    parser.add_argument("--guess_charge", help="guess the charge", action="store_true", default=False)
    parser.add_argument("--memory", help="memory to use in calculation in GB", type=int, default=1)
    parser.add_argument("--memory_ddi", help="distributed memory to use in GAMESS calculations in GB", type=int, default=1)
    parser.add_argument("--path_to_psi4", help="path to the psi4 executable", default="psi4")
    args = parser.parse_args()
    main(args)

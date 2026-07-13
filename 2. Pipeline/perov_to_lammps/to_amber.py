import os
import subprocess
from typing import List

class ToAmber:
    def __init__(self, cif_name: str, amount_unique_spacers: int, output_dir: str, charge_per_spacer: List[int], box_lengths: List[float], box_angles: List[float]):
        self.cif_name = cif_name
        self.amount_unique_spacers = amount_unique_spacers
        self.output_dir = output_dir
        self.charge_per_spacer = charge_per_spacer
        self.box_lengths = box_lengths # Nieuw
        self.box_angles = box_angles   # Nieuw

    def parameterise_and_convert_to_mol2(self) -> bool:
        """
        For each unique spacer, first run OpenBabel to convert from XYZ to MOL2 format, then use Antechamber to
        parameterise the molecule using the GAFF2 force field and the BCC charge model. Finally, generate the
        corresponding frcmod file using parmchk2. All files are saved in output_dir.

        Note: This requires that OpenBabel, Antechamber, and Parmchk2 are installed and accessible from the command
        line.
        
        :return: Returns True if all commands executed successfully, otherwise False is returned.
        """
        # Loop over each unique spacer molecule
        for i in range(self.amount_unique_spacers):
            base = f"{self.cif_name}_spacer_{i}"
            xyz = f"{base}.xyz"
            mol2 = f"{base}.mol2"
            gaff_mol2 = f"{base}_gaff.mol2"
            frcmod = f"{base}.frcmod"
            charge = self.charge_per_spacer[i]

            # Use OpenBabel to convert from XYZ to MOL2, then Antechamber to parameterise, then parmchk2 to generate
            # frcmod file, which is needed for tleap (simulation box creation) later on.
            commands = [f"obabel -ixyz {xyz} -omol2 -O {mol2} --gen3d",
                        f"antechamber -i {mol2} -fi mol2 -o {gaff_mol2} -fo mol2 -c bcc -at gaff2 -s 2 -nc {charge}",
                        f"parmchk2 -i {gaff_mol2} -f mol2 -o {frcmod}"]

            # Run commands
            for command in commands:
                print(f"Running in {self.output_dir}: {command}")
                result = subprocess.run(command, shell=True, cwd=self.output_dir, capture_output=True, text=True)
                
                if result.returncode != 0:
                    print(f"Error running '{command}':\n{result.stderr}")
                    return False
                else:
                    print(result.stdout)
                    
        return True

    def build_tleap_input(self, box_lengths: tuple[float, float, float]) -> str:
        """
        Create the tleap input script to build the full system. First, load the necessary parameters and mol2 files,
        then combine them into a single system, define the box dimensions, and finally save the Amber parameter and
        coordinate files, namely system.prmtop and system.crd.

        Note: This requires that tleap (part of AmberTools) is installed and accessible from the command line, and
        expects the inorg.frcmod and all spacer mol2 and frcmod files to be present in output_dir.
        
        :param box_lengths: Tuple containing the box lengths (a, b, c).
        :return: Path to the generated tleap input file.
        """
        tleap_path = os.path.join(self.output_dir, "tleap_build.in")

        lines = ["source leaprc.gaff2\n"]
        lines.append("set default nocenter on\n\n")

        # Add atom types for inorganic elements; since tLeap is built for organice molecules, these atom types need to
        # be manually added to make sure they get included in the final simulation box, since adding them afterwards is
        # not possible, and not defining them throws errors.
        lines.append("""addAtomTypes {
  { "PB" "Pb" "sp3" }
  { "SN" "Sn" "sp3" }
  { "I"  "I"  "sp3" }
  { "BR" "Br" "sp3" }
  { "CL" "Cl" "sp3" }
  { "CU" "Cu" "sp3" }
  { "CS" "Cs" "sp3" }
}\n\n""")

        # Load frcmod for each spacer and inorganics
        for i in range(self.amount_unique_spacers):
            lines.append(f"loadamberparams {self.cif_name}_spacer_{i}.frcmod\n")
        lines.append(f"loadamberparams inorg.frcmod\n\n")

        # Load combined mol2 files for each spacer symmetry
        all_sys_names = ""
        for i in range(self.amount_unique_spacers):
            lines.append(f"sys{i} = loadMol2 {self.cif_name}_combined_spacers_{i}_gaff.mol2\n")
            all_sys_names += f" sys{i} "
        lines.append(f"inorg = loadMol2 {self.cif_name}_inorganics.mol2\n\n")

        # Combine into one block
        lines.append(f"sys = combine {{{all_sys_names} inorg }}\n\n")

        # Add box definition
        a,b,c = box_lengths
        lines.append(f"set sys box {{ {a} {b} {c} }}\n")

        # Save output
        lines.append("saveAmberParm sys system.prmtop system.crd\n")
        lines.append("quit\n")

        # Write tleap input file
        with open(tleap_path, "w") as file:
            file.writelines(lines)

        print(f"Wrote tleap input script: {tleap_path}")
        return tleap_path

    def run_tleap(self, tleap_input: str="tleap_build.in") -> bool:
        """
        Run tleap using the script as generated above to generate the Amber parameter and coordinate files.
        
        :param tleap_input: Path to the tleap input script.
        :return: Returns True if tleap executed successfully, otherwise False is returned.
        """
        command = ("cp ../../inorg.frcmod . && "
                   f"tleap -f {tleap_input}")
        
        print(f"Running tleap in {self.output_dir}: {command}")
        result = subprocess.run(command, shell=True, cwd=self.output_dir, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"TLeap error:\n{result.stderr}")
            return False
        else:
            print(result.stdout)
            print("TLeap completed successfully.\n")
        return True

    def convert_to_lammps(self) -> bool:
        """
        Convert the generated Amber files (system.prmtop and system.crd) to LAMMPS format using InterMol. First, convert
        from ParmEd to GROMACS format, then use InterMol to convert from GROMACS to LAMMPS. The final LAMMPS files
        are named data_<cif_name>.txt and extra_<cif_name>.txt.

        Note: This requires that ParmEd and InterMol are installed and accessible from the command line.

        :return: Returns True if conversion was successful, otherwise False is returned.
        """
        
        a, b, c = self.box_lengths
        alpha, beta, gamma = self.box_angles
        
        command = ("python - << 'PY'\n"
                   "import parmed as pmd\n"
                   "amb = pmd.load_file('system.prmtop', xyz='system.crd')\n"
                   f"amb.box = [{a}, {b}, {c}, {alpha}, {beta}, {gamma}]\n" # Dwing de tricliene box af
                   "amb.save('system.top', format='gromacs', combine='all', parameters='inline')\n"
                   "amb.save('system.gro', combine='all')\n"
                   "PY\n"
                   "python -m intermol.convert --gro_in system.top system.gro --lammps --oname system\n"
                   f"mv system.lmp data_{self.cif_name}.txt\n"
                   f"mv system.input extra_{self.cif_name}.txt")

        result = subprocess.run(command, shell=True, cwd=self.output_dir, capture_output=True, text=True)
        if result.returncode != 0: print(f"Error during LAMMPS conversion:\n{result.stderr}"); return False
        print(result.stdout)
        return True
from openbabel import pybel
from openbabel.pybel import Molecule
from typing import List

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
import os



class FileWriter:
    """
    Class with different kinds of writer functions that are used throughout the pipeline.
    """

    def __init__(self, cif_file: str, output_dir: str) -> None:
        self.cif_file = cif_file
        self.output_dir = output_dir

    def write_to_xyz(self, molecule: Molecule, filename: str) -> bool:
        """
        Write a given molecule to an XYZ file.

        :param molecule: Molecule instance representing the spacer molecule.
        :param filename: Name of the output XYZ file.
        :return: Writes molecule to XYZ and returns True if writing was successful.
        """
        output_path = os.path.join(self.output_dir, filename)
        molecule.to(fmt="xyz", filename=output_path)
        return True

    def write_all_parameterised_symmetries_to_mol2(self, molecule_group: List[Molecule], template_mol2_name: str) -> str | None:
        """
        Use a template mol2 file with correct atom types to create a combined mol2 file for all symmetries of a unique
        spacer molecule.

        :param molecule_group: List of Molecule instances representing different symmetries of a unique spacer.
        :param template_mol2_name: Name of the template mol2 file located in the output directory.
        :return: Combined mol2 file content as a string.
        """
        template_mol2_path = os.path.join(self.output_dir, template_mol2_name)

        try:
            with open(template_mol2_path, 'r') as file:
                template_mol2_text = file.read()
        except Exception as error: print(f"Error reading template file: {error}"); return None

        # Check if required sections exist
        required_sections = ['@<TRIPOS>MOLECULE', '@<TRIPOS>ATOM', '@<TRIPOS>BOND', '@<TRIPOS>SUBSTRUCTURE']
        for section in required_sections:
            if section not in template_mol2_text:
                print(f"Error: {section} section not found in template file"); return None

        # Get sections from template
        sections = template_mol2_text.split('@<TRIPOS>')

        # Extract atom information (types, charges) from template
        atom_section = template_mol2_text.split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0].strip()
        atom_lines = atom_section.split('\n')
        atom_types = []
        atom_charges = []

        for line in atom_lines:
            parts = line.split()
            if len(parts) >= 9:
                atom_types.append(parts[5])  # atom type
                atom_charges.append(float(parts[8]))  # charge

        # Extract bond information from template
        bond_section = template_mol2_text.split('@<TRIPOS>BOND')[1].split('@<TRIPOS>')[0].strip()
        bond_lines = bond_section.split('\n')
        bond_info = []

        for line in bond_lines:
            parts = line.split()
            if len(parts) >= 4:
                bond_info.append((int(parts[1]), int(parts[2]), parts[3]))  # atom1, atom2, bond_type

        # Create combined mol2 file
        num_atoms_per_mol = len(molecule_group[0])
        num_bonds_per_mol = len(bond_info)
        num_molecules = len(molecule_group)
        total_atoms = num_atoms_per_mol * num_molecules
        total_bonds = num_bonds_per_mol * num_molecules

        combined_text = '@<TRIPOS>MOLECULE\n'
        combined_text += 'Combined Molecules\n'
        combined_text += f'{total_atoms} {total_bonds} {num_molecules} 0 0\n'
        combined_text += 'SMALL\nbcc\n\n\n'

        # Add atoms section
        combined_text += '@<TRIPOS>ATOM\n'
        atom_count = 1
        for mol_idx, molecule in enumerate(molecule_group):
            for atom_idx, (site, atom_type, charge) in enumerate(zip(molecule, atom_types, atom_charges)):
                x, y, z = site.coords
                element = str(site.species)
                combined_text += f'{atom_count:7d} {element:<4s} {x:9.4f} {y:9.4f} {z:9.4f} {atom_type:<4s} '
                combined_text += f'{mol_idx + 1:7d} MOL{mol_idx + 1:<3d} {charge:9.6f}\n'
                atom_count += 1

        # Add bonds section
        combined_text += '@<TRIPOS>BOND\n'
        bond_count = 1
        for mol_idx in range(num_molecules):
            atom_offset = mol_idx * num_atoms_per_mol
            for b_idx, (atom1, atom2, bond_type) in enumerate(bond_info):
                new_atom1 = atom1 + atom_offset
                new_atom2 = atom2 + atom_offset
                combined_text += f'{bond_count:7d} {new_atom1:7d} {new_atom2:7d} {bond_type}\n'
                bond_count += 1

        # Add substructure section
        combined_text += '@<TRIPOS>SUBSTRUCTURE\n'
        for mol_idx in range(num_molecules):
            combined_text += f'{mol_idx + 1:7d} MOL{mol_idx + 1:<3d} {mol_idx * num_atoms_per_mol + 1:7d} TEMP 0 **** **** 0 ROOT\n'

        return combined_text

    def write_inorganics_mol2(self, inorganic_structure: Structure, filename: str) -> str:
        """
        Write the inorganic structure to a mol2 file. Note, bonds are not included in this file, since there are no
        bonds in the inorganic framework. The BOND section is still required for the mol2 format, but is empty.

        :param inorganic_structure: Pymatgen Structure instance representing the inorganic framework.
        :param filename: Name of the output mol2 file.
        :return: Path to the written mol2 file.
        """
        # Charges for common inorganic elements
        charge_type = {"Pb": 2.0, "Sn": 2.0, "I": -1.0, "Br": -1.0, "Cl": -1.0, "Cu": 1.0, "Cs": 1.0}
        output_path = os.path.join(self.output_dir, filename)

        atoms = []

        # Loop through each site in the inorganic structure and collect atom data like element, coordinates, and charge,
        # after which they are stored in the atoms list and written to the mol2 file.
        for i, site in enumerate(inorganic_structure, start=1):
            element = site.species_string
            x, y, z = site.coords
            atom_type = element.upper()  # tie this to frcmod types below
            name = f"{element}{i}"
            charge = charge_type.get(element, 0.0)
            atoms.append((i, name, x, y, z, atom_type, 1, "INORG", charge))

        with open(output_path, "w") as file:
            file.write("@<TRIPOS>MOLECULE\n")
            file.write("INORGANIC_BLOCK\n")
            file.write(f"{len(atoms)} 0 1 0 0\n")
            file.write("SMALL\n")
            file.write("USER_CHARGES\n\n")
            file.write("@<TRIPOS>ATOM\n")
            for i, name, x, y, z, atom_type, resid, resname, charge in atoms:
                file.write(
                    f"{i:7d} {name:<4s} {x:10.4f} {y:10.4f} {z:10.4f} {atom_type:<4s} {resid:5d} {resname:<8s} {charge:9.6f}\n")
            file.write("\n@<TRIPOS>BOND\n")
            file.write("\n@<TRIPOS>SUBSTRUCTURE\n")
            file.write(f"{1:7d} INORG    {1:7d} TEMP 0 **** **** 0 ROOT\n")
        return output_path
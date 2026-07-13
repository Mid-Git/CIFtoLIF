import os
import pandas as pd
from pymatgen.core import Molecule, Structure
from typing import List, Tuple, Optional

class SystemAnalyser:
    """
    Class to perform system-wide analysis, such as determining charges based on the unit cell composition
    and calculating dynamic replication factors for LAMMPS simulations.
    """
    
    @staticmethod
    def determine_automatic_charge(inorganic_structure: Structure, all_spacers: List[List[Molecule]]) -> Optional[List[int]]:
        """
        Calculates the net charge of the inorganic framework in the unit cell and deduces
        the charge per organic molecule.
        
        :param inorganic_structure: Pymatgen Structure instance representing the inorganic framework.
        :param all_spacers: A nested list containing all unique spacer molecules and their symmetries.
        :return: A list of integers representing the charge per spacer, or None if automatic detection fails.
        """
        inorg_charges = {
            "Pb": 2, "Sn": 2, "Cu": 2, "Cs": 1,
            "I": -1, "Br": -1, "Cl": -1
        }
        
        total_inorg_charge = 0
        for site in inorganic_structure:
            element = site.species.elements[0].symbol
            if element in inorg_charges:
                total_inorg_charge += inorg_charges[element]
                
        total_organic_charge = -total_inorg_charge
        amount_unique_spacers = len(all_spacers)
        total_nitrogens_in_cell = 0
        
        for spacer_group in all_spacers:
            ref_mol = spacer_group[0]
            num_copies = len(spacer_group)
            n_count = sum(1 for atom in ref_mol.species if atom.symbol == "N")
            total_nitrogens_in_cell += n_count * num_copies

        if total_nitrogens_in_cell > 0 and total_organic_charge == total_nitrogens_in_cell:
            charge_list = []
            print(" -> Nitrogen-heuristic activated")
            print(f"    Total charge ({total_organic_charge:+.0f}) equal to amount of nitrogen atoms")
            
            for spacer_group in all_spacers:
                ref_mol = spacer_group[0]
                n_count = sum(1 for atom in ref_mol.species if atom.symbol == "N")
                charge_list.append(n_count)
                print(f"    Assigned {ref_mol.composition.formula} -> {n_count:+}")
            return charge_list
        
        if amount_unique_spacers == 1:
            total_spacer_molecules = len(all_spacers[0])
            charge_per_spacer = total_organic_charge / total_spacer_molecules
            
            if abs(charge_per_spacer - round(charge_per_spacer)) < 1e-4:
                calculated_charge = int(round(charge_per_spacer))
                spacer_formula = all_spacers[0][0].composition.formula
                
                print(" -> Automatically assigned charge for spacer:")
                print(f"    {spacer_formula} -> {calculated_charge:+}")
                return [calculated_charge]
            else:
                print(f" -> Warning: fractional charge found for ion ({charge_per_spacer}).")
                return None
        else:
            print(" -> Multiple spacers found: no automatic charge detection")
            return None
        
    @staticmethod
    def calculate_replication_factors(box_lengths: List[float]) -> Tuple[List[int], str]:
        """
        Determines the Z-axis (longest side) and calculates dynamic replication factors for LAMMPS
        to scale the simulation box to a realistic size (~25x25x45 Å).
        
        :param box_lengths: List containing the three box lengths [a, b, c].
        :return: A tuple containing the replication multipliers [nx, ny, nz] and the identified Z-axis string ('a', 'b', or 'c').
        """
        lengths = box_lengths
        sorted_lengths = sorted(lengths, reverse=True)
        
        if sorted_lengths[0] < 1.5 * sorted_lengths[1]:
            print(" -> !!! Z-AXIS DETECTION NOT SUFFICIENT !!!")
            print(f"    Longest axis ({sorted_lengths[0]:.2f} Å) and second longest axis ({sorted_lengths[1]:.2f} Å) do not differ enough.")
            
            while True:
                print(f"    a = {lengths[0]:.2f}, b = {lengths[1]:.2f}, c = {lengths[2]:.2f}")
                user_axis = input("    Axis perpendicular to 2D layers (choose 'a', 'b', or 'c'): ").strip().lower()
                if user_axis in ['a', 'b', 'c']:
                    z_axis_idx = ['a', 'b', 'c'].index(user_axis)
                    break
                else:
                    print("    Invalid response. Type a, b or c.")
        else:
            z_axis_idx = lengths.index(max(lengths))
        
        replicates = [1, 1, 1]
        for i in range(3):
            if i == z_axis_idx:
                rep = max(2, round(60.0 / lengths[i]))
                replicates[i] = int(rep)
            else:
                rep = max(4, round(50.0 / lengths[i]))
                replicates[i] = int(rep)
                
        z_axis_str = ['a', 'b', 'c'][z_axis_idx]
        
        print(" -> Z-axis detection summary:")
        print(f"    Longest axis (Z-axis) : Axis {z_axis_str} ({lengths[z_axis_idx]:.2f} Å)")
        print(f"    Calculated replication: {replicates[0]}x{replicates[1]}x{replicates[2]}")
        
        return replicates, z_axis_str

    @classmethod
    def log_z_axis_to_csv(cls, cif_name: str, z_axis_str: str, csv_path: str = "combined_perovskite_data.csv") -> None:
        """
        Logs the identified Z-axis orientation to the specified column in the dataset.
        
        :param cif_name: Name of the CIF file being processed.
        :param z_axis_str: The identified Z-axis ('a', 'b', or 'c').
        :param csv_path: Path to the target CSV file.
        """
        if os.path.exists(csv_path):
            try:
                df = pd.read_csv(csv_path, sep=';')
                
                while len(df.columns) < 17:
                    df[f'Column_{len(df.columns)+1}'] = None
                    
                col_cif = df.columns[0]
                col_z_axis = df.columns[16]
                
                if df.columns[16].startswith('Unnamed') or df.columns[16].startswith('Column_'):
                    df.rename(columns={col_z_axis: 'Z_Axis_Direction'}, inplace=True)
                    col_z_axis = 'Z_Axis_Direction'
                
                row_index = df.index[df[col_cif].astype(str) == str(cif_name)].tolist()
                
                if row_index:
                    df.at[row_index[0], col_z_axis] = z_axis_str
                    df.to_csv(csv_path, sep=';', index=False)
                else:
                    print(f" -> Warning: CIF '{cif_name}' not found in column A of {csv_path}. Z-axis tracking skipped.")
            except Exception as e:
                print(f" -> Error updating {csv_path}: {e}")
        else:
            print(f" -> Warning: File {csv_path} not found. Z-axis tracking skipped.")
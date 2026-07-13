import os
from perov_to_lammps import CIFLoader, MoleculeIdentifier, SystemAnalyser, FileWriter, ToAmber, LammpsFinaliser

def pipeline(cif_path: str, config: dict) -> bool:
    """
    Executes the full conversion pipeline for a single CIF file.

    :param cif_path: Path to the input CIF file.
    :param config: Dictionary containing boolean flags for pipeline steps.
    :return: True if the pipeline completes successfully, False otherwise.
    """
    cif_name = os.path.basename(cif_path).replace(".cif", "")
    output_dir = f"intermediate_files/{cif_name}_files"
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'=' * 60}\n| STARTING PIPELINE FOR CIF: {cif_name:<27} |\n{'=' * 60}")

    loader = CIFLoader(cif_path)
    state = loader.run()
    if not state:
        return False
    organic_structure, inorganic_structure, box_lengths, box_angles = state

    print(" -> Identifying spacer molecules...")
    mol_identifier = MoleculeIdentifier(organic_structure)
    all_spacers = mol_identifier.run()
    all_spacers = mol_identifier.reorder_all_spacer_atom_indices(all_spacers)
    if not all_spacers:
        return False
    amount_unique_spacers = len(all_spacers)

    charge_per_spacer = SystemAnalyser.determine_automatic_charge(inorganic_structure, all_spacers)
    if charge_per_spacer is None:
        charge_per_spacer = []
        for i in range(amount_unique_spacers):
            formula = all_spacers[i][0].composition.formula
            while True:
                user_input = input(f" -> Spacer charge for {formula} is: ").strip()
                if user_input.lower() == "stop":
                    return False
                try:
                    charge_per_spacer.append(int(user_input))
                    break
                except ValueError:
                    print("    Please enter an integer (e.g., -1, 0, +1). Type 'stop' to abort.")

    writer = FileWriter(cif_path, output_dir)
    to_amber = ToAmber(cif_name, amount_unique_spacers, output_dir, charge_per_spacer, box_lengths, box_angles)

    if config.get("write_xyz"):
        print(" -> Writing spacers to XYZ...")
        for i in range(amount_unique_spacers):
            writer.write_to_xyz(all_spacers[i][0], filename=f"{cif_name}_spacer_{i}.xyz")

    if config.get("write_params"):
        print(" -> Parameterising spacers and converting to MOL2...")
        if not to_amber.parameterise_and_convert_to_mol2(): return False

    if config.get("write_parameterised_symmetries_mol2"):
        print(" -> Writing all spacer symmetries...")
        for i in range(amount_unique_spacers):
            template_name = f"{cif_name}_spacer_{i}_gaff.mol2"
            combined_text = writer.write_all_parameterised_symmetries_to_mol2(all_spacers[i], template_name)
            with open(os.path.join(output_dir, f"{cif_name}_combined_spacers_{i}_gaff.mol2"), "w") as f:
                f.write(combined_text)

    if config.get("write_inorganics_mol2"):
        print(" -> Writing inorganics to MOL2...")
        writer.write_inorganics_mol2(inorganic_structure, f"{cif_name}_inorganics.mol2")

    if config.get("write_and_run_tleap"):
        print(" -> Building and running tLeap script...")
        to_amber.build_tleap_input(box_lengths)
        if not to_amber.run_tleap(): return False

    if config.get("write_to_lammps"):
        print(" -> Converting Amber files to LAMMPS...")
        if not to_amber.convert_to_lammps(): return False

    replicate_factors, z_axis_str = SystemAnalyser.calculate_replication_factors(box_lengths)
    SystemAnalyser.log_z_axis_to_csv(cif_name, z_axis_str)

    lammps = LammpsFinaliser(cif_name, output_dir, replicate_factors)

    if config.get("gather_coeffs"):
        organic_pair_coeffs_block = lammps.extract_and_modify_pair_coeff_block()
    if config.get("gather_atom_labels"):
        lammps.get_atom_ids()
    if config.get("write_coeffs_and_input"):
        lammps.coeffs_writer(organic_pair_coeffs_block)
        lammps.input_writer()

    lammps.save_everything_to_folder()
    print(f" -> Pipeline completed successfully for {cif_name}.")
    return True


def main() -> None:
    """
    Main execution loop that iterates over all CIF files in the input directory and logs failures.
    """
    cif_dir = "input_cifs"
    failed_log = "failed_cifs.txt"
    config = {
        "write_xyz": True,
        "write_params": True,
        "write_parameterised_symmetries_mol2": True,
        "write_inorganics_mol2": True,
        "write_and_run_tleap": True,
        "write_to_lammps": True,
        "gather_atom_labels": True,
        "gather_coeffs": True,
        "write_coeffs_and_input": True
    }

    if not os.path.exists(failed_log):
        with open(failed_log, "w") as log:
            log.write("CIF_FILE\tERROR_MESSAGE\n")

    cifs = [f for f in os.listdir(cif_dir) if f.lower().endswith(".cif")]

    for cif in cifs:
        cif_path = os.path.join(cif_dir, cif)
        try:
            if not pipeline(cif_path, config):
                with open(failed_log, "a") as log:
                    log.write(f"{cif}\tPipeline returned False or stopped manually\n")
        except Exception as e:
            print(f" -> CRITICAL ERROR for {cif}: {e}")
            with open(failed_log, "a") as log:
                log.write(f"{cif}\t{str(e).replace(chr(10), ' ')}\n")

if __name__ == "__main__":
    main()
from pymatgen.core import Structure, Molecule
import itertools
import warnings


class CIFLoader:
    """
    Class that is used to load in a CIF, using file_loader(), and then prune all excess H-atoms using prune_H_atoms().
    These two functions can automatically be called by simply calling the class and running the run() function. For
    information about any of the other functions, see the respective documentation.
    """
    # Max. amount of H-atoms per other atom, assuming the other atom has at least one non-H-atom bond.
    max_H_atoms = {"C": 3, "N": 3, "O": 1, "S": 1}

    # Max. allowed distance for H-atom bonding with parents (used in prune_H_atoms())
    max_parent_dist = 1.65

    def __init__(self, cif_file: str) -> None:
        self.cif_file = cif_file  # Location of the current CIF (e.g. "my.cif" or "CIFs/my.cif")
        self.structure = None  # Will containing pymatgen Structure instance based on the CIF.
        self.state = True  # State is used to check whether the loading of a CIF was successful.

    def run(self) -> tuple[Structure, Structure, list, list] | bool:
        """
        Function that automatically runs the file_loader() and prune_H_atoms() functions, without having to explicitly
        Call them first. Saves some code in main.py, and seeing as these functions will not be used seperately, simply
        calling the class and running the run() function will do everything automatically.

        :return: Tuple containing the full structure, inorganic structure, box lengths, and box angles; or False if
        something went wrong.
        """
        state = self.file_loader()

        # Always return False when ANYTHING goes wrong. This is not the most informative, but since a lot of things can
        # go wrong when loading in and parsing CIFs.
        if not state:
            return False

        try:
            # First, extract box lengths and angles, then inorganic atoms, and finally spacers molecules. Note that,
            # excess H-atoms are pruned from the spacer molecules.
            lattice = self.structure.lattice
            self.box_lengths = [lattice.a, lattice.b, lattice.c]
            self.box_angles = [lattice.alpha, lattice.beta, lattice.gamma]
            self.inorg_structure = self.get_inorganics(self.structure)
            self.get_spacers()
            self.prune_H_atoms()
            return self.structure, self.inorg_structure, self.box_lengths, self.box_angles
        except Exception as error: print(f"Error during run() when reading CIF: {error}"); return False

    def file_loader(self) -> bool:
        """
        Function that loads in the CIF, and catches any errors. In case Pymatgen throws an error saying there are too few H-atoms,
        also immediately quit. If it gives an exception error, the formatting of the CIF is broken.

        :return: The read-in CIF structure
        :rtype: pymatgen.core.Structure
        """
        with open(self.cif_file, "r", encoding="utf-8", errors="ignore") as file:
            cif_text = file.read()

        # Add pymatgen error to catch missing H-atoms
        warnings.filterwarnings("error",
                                message=r"Missing elements\s+H\s+from PMG structure composition",
                                category=UserWarning,
                                module=r"pymatgen\.io\.cif")

        # Try initialising structure; stop if there are too few H-atoms or the formatting is broken
        try:
            self.structure = Structure.from_str(cif_text, fmt="cif")
            
        except UserWarning:
            print("Your CIF is probably missing H-atoms; stopping code.")
            return False

        except Exception:
            print("Your CIF is broken; stopping code.")
            return False

        return True

    def get_spacers(self) -> Structure:
        """
        Function that extracts all organic atoms from the CIF structure, leaving only the organic atoms belonging to the
        spacer molecules, which are identified later. The inorganic atoms are filtered out based on the "inorg" set
        defined below.

        :return: The read-in CIF, with only organic atoms.
        """
        inorg = {"Pb", "Sn", "I", "Br", "Cl", "Cu", "Cs"}

        # Go through the entire CIF structure, and keep all organic atoms
        keep = []
        for index, location in enumerate(self.structure):
            if location.species_string not in inorg:
                keep.append(index)

        # Check which indices are NOT in keep, and remove them from them structure.
        remove_index = []
        for index in range(len(self.structure)):
            if index not in keep:
                remove_index.append(index)
        structure_organic_only = self.structure.copy()
        structure_organic_only.remove_sites(remove_index)

        self.structure = structure_organic_only
        return self.structure

    # snippet highlighting the modified function inside init_cif.py

    def prune_H_atoms(self) -> Structure:
        """
        Function that resolves symmetric/duplicate H-atoms in CIF structures.
        If any non-H atom has a fractional occupancy (< 1.0), it raises a RuntimeError
        to allow the pipeline to log it and skip the specific CIF file.

        :return: Pymatgen Structure with fully occupied sites and pruned overlapping hydrogens.
        """
        non_H_atoms = []
        H_atoms_single = []
        H_atoms_double = []
        
        for index, location in enumerate(self.structure):
            element_str = location.species.elements[0].symbol
            occ = location.species.get_el_amt_dict().get(element_str, 0)
            
            if element_str != "H" and occ < 1.0:
                raise RuntimeError(f"Fractional occupation ({occ}) found for atom '{element_str}'. This pipeline only supports disorder in H atoms.")
            
            if element_str != "H":
                non_H_atoms.append(index)
            else:
                species_str = str(self.structure[index].species)
                if species_str == "H0.5":  
                    H_atoms_double.append(index)
                elif species_str == "H1" or species_str == "H": 
                    H_atoms_single.append(index)
                else:
                    raise RuntimeError(f"Unknown fractional H-atom found ({species_str}). Only H1 and H0.5 are supported.")

        keep = set()
        keep.update(non_H_atoms)
        keep.update(H_atoms_single)

        H_double_pairs = list(itertools.combinations(H_atoms_double, 2))
        distances = []
        for atom_1, atom_2 in H_double_pairs:
            distance = self.structure.get_distance(atom_1, atom_2)
            distances.append((atom_1, atom_2, distance))

        discard = set()
        for atom_1, atom_2, distance in distances:
            if distance < 0.5:
                if atom_2 not in discard:
                    keep.add(atom_1)
                    discard.add(atom_2)

        for atom in H_atoms_double:
            if atom not in discard:
                keep.add(atom)

        locations_to_keep = []
        for location in keep:
            locations_to_keep.append(self.structure[location])
            
        self.structure = Structure.from_sites(locations_to_keep)
        return self.structure
    
    def get_inorganics(self, structure: Structure) -> Structure | None:
        """
        Extract inorganic atoms from the CIF structure, based on the "inorg" set defined below, and already translates
        the fractional coordinates to cartesian coordinates.

        :param structure: The full CIF structure.
        :return: The CIF structure with only inorganic atoms, or None if something went wrong.
        """
        inorg = {"Pb", "Sn", "I", "Br", "Cl", "Cu", "Cs"}
        locations_to_keep = []
        for site in structure:
            element = site.species.elements[0].symbol
            if element in inorg:
                frac_coords = site.frac_coords % 1.0
                cart_coords = structure.lattice.get_cartesian_coords(frac_coords)
                new_site = site.__class__(species=site.species,
                                          coords=cart_coords,
                                          lattice=structure.lattice,
                                          coords_are_cartesian=True)
                locations_to_keep.append(new_site)

        if not locations_to_keep:
            print("Warning: no inorganics found in CIF. Check inorg list or CIF.")
            return None

        return Structure.from_sites(locations_to_keep)


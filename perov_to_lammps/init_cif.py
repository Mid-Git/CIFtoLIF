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
from pymatgen.core import Molecule, Structure
import networkx as nx
import numpy as np
from typing import List

class MoleculeIdentifier:
    """
    A class that identifies each unique spacer molecule based on a ``pymatgen.core.Structure`` (which in this case will
    be the structure gathered from the CIF from ``init_cif.py``).
    """
    bond_tolerance = 0.2
    bond_radius = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "F": 0.57, "Cl": 1.02, "Br": 1.20, "I": 1.39, "P": 1.07}
    inorg = {"Pb", "Sn", "I", "Br", "Cl", "Cu", "Cs"}

    def __init__(self, structure: Structure) -> None:
        self.structure = structure
        self.molecules = []

    def run(self) -> List[List[Molecule]]:
        """
        Function that, when called, runs the ``build_graph()`` to build a graph to identify the molecules in the CIF,
        still containing the wrong edges within each graph, after which ``identify_molecules()`` is used to infer
        correct bonds and transform each correctly connected graph (molecule) in a fractional cell.

        :return: A list within a list, where each outer list contains an inner list of all spacer molecule symmetries
        with each uniquely identified Molecule instance.
        """
        self.build_graph()
        self.identify_molecules()
        return self.molecules

    def build_graph(self) -> nx.Graph:
        """
        This functions builds a graph by first adding a node with each atom symbol, and then looping over each pair of
        neighbours, calculating a cutoff, and then adding edges between atoms each time the calculated distance
        in the lattice is lower (or equal) to the cutoff.

        Note that, here, a graph is created based on a cutoff distance (e.g. neighbours can not be further than x Å).
        This means the graph has INCORRECT edges and thus a perception of the molecule, but this doesn't matter because
        it DOES get the correct atoms for a single molecule, and thus correctly identifies the indices of the atoms
        (from the CIF structure) belonging to a single spacer molecule. Therefore, it correctly identifies atoms
        belonging to the same molecule, but NOT the way it's bonded. The correct edges are found in ``build_graph()``.

        :return: Graph of the structure, where each node is an atom and the edges are (incorrect) bonds between atoms
        """
        structure = self.structure
        cif_graph = nx.Graph()  # Initialise a graph structure to start building molecule

        # Add each atom in the structure as a node in the graph
        for index, atom in enumerate(structure):
            cif_graph.add_node(index, element=atom.species_string)

        max_cutoff = 3  # There DEFINITELY will not be any (organic) bonds past 3 Å
        neighbors_all = structure.get_all_neighbors(max_cutoff, include_index=True)

        for atom_i_index, neighbour_list in enumerate(neighbors_all):
            atom_i = structure[atom_i_index].species.elements[0].symbol  # Get element; for bond_radii list

            # Loop over each neighbouring atom based on the neighbouring atoms found using the cutoff
            # Here, neighbour[1] is the distance to the neighbouring atom j (from atom i)
            # Here, neighbour[2] is the index of the neighbouring atom j (from atom i)
            for neighbour in neighbour_list:
                atom_j_index = neighbour[2]
                if atom_j_index <= atom_i_index:  # Make sure pairs aren't repeated
                    continue
                atom_j = structure[atom_j_index].species.elements[0].symbol  # Get element; for bond_radii list

                # Check if both atoms are in the bond_radii list (and thus organic)
                # If so, determine cutoff, and assemble atoms belonging to the same molecule in the graph
                if atom_i in self.bond_radius and atom_j in self.bond_radius:
                    cutoff = self.bond_radius[atom_i] + self.bond_radius[atom_j] + self.bond_tolerance

                    dist_to_neighbour = neighbour[1]
                    if dist_to_neighbour <= cutoff:
                        cif_graph.add_edge(atom_i_index, atom_j_index)
        return cif_graph

    def identify_molecules(self) -> List[List[Molecule]]:
        """
        This function identifies organic molecules based on a pymatgen Structure and the graph from ``build_graph()``;
        the inorganic atoms are not included in the graph. The function uses a BFS approach to go over each neighbour
        in the graph from the previous function (cif_graph), and perform a transformation on each atom if needed. It
        then moves any molecules still outside the (0,1) boundary inside.

        Two transformations are performed in this function. First, per molecule, per atom each image vector is
        determined to make sure all atoms in the molecule are between 0 and 1. After doing that for all molecules,
        it is still possible for one or more of the molecules to lay outside 0 and 1; you start "building" the molecules
        based on a reference atom, but this reference atom will NOT be at (0,0,0), meaning some molecules might fall
        outside 0 and 1 because of the reference atom not being at the origin. This transformation is done as follows:
        Per molecule, the lowest x, y and z of any atom in the molecule is found, this value is then floored and given
        the opposite sign. Adding this value to all respective coordinates (x, y or z) moves the entire molecule inwards
        in the unit cell.

        Note that, it's not actually each neighbour, but rather each atom in said atom's vicinity of 3 Å (see
        ``build_graph()``). This still works, because you don't need to go over each specific neighbour, but rather just
        each atom in a single molecule, thus if the neighbours are correct does not matter.

        :return: A list within a list, where each outer list contains an inner list of all spacer molecule symmetries
        with each uniquely identified Molecule instance.
        """
        structure = self.structure
        cif_graph = self.build_graph()
        molecules = []

        # Loop over all molecules, with the intention to make sure all their coordinates are within the fractional
        # boundaries of 0 and 1.
        for molecule in nx.connected_components(cif_graph):
            molecule = list(molecule)

            shift = {}  # A dictionary of image vectors to know by how much the frac. coordinates should be shifted
            queue = []

            # Get the first atom of the first molecule, add it to shift, and use it as a reference for all other atoms
            first_atom = molecule[0]
            shift[first_atom] = np.array([0, 0, 0])
            queue.append(first_atom)

            # Use a BFS approach: Get an atom, get its neighbours from the cif_graph, calculate the neighbours' image
            # vectors based on the position relative to the initial atom, and then add the neighbours to the queue to do
            # the same for theirs.
            while queue:
                atom_i = queue.pop(0)
                for atom_j in cif_graph.neighbors(atom_i):
                    if atom_j not in shift:
                        dist, image = structure[atom_i].distance_and_image(structure[atom_j])
                        shift[atom_j] = shift[atom_i] + np.array(image)
                        queue.append(atom_j)

            # Now, actually perform the shifts on the coordinates, based on the image vectors gathered by BFS
            shifted_frac_coords = []
            for atom in molecule:
                shifted_coords = structure[atom].frac_coords + shift[atom]
                shifted_frac_coords.append(shifted_coords)
            shifted_coords_np = np.asarray(shifted_frac_coords)

            # The shifts based on the image vectors are all relative to the reference atom. However, this means that
            # the shifts can lie outside the unit cell (e.g. ([0.2,0.5,1.2]), where 1.2 is outside the fractional
            # coordinates, between 0 and 1. To fix this get the outer (np.min) x,y,z coordinates of the enitre molecule,
            # use floor to round the values down, and then use the opposite sign to shift them inwards of the unit cell.
            mins = shifted_coords_np.min(axis=0)
            enitre_mol_shift = -np.floor(mins)

            # Use the shift as calculated above to shift ALL frac. coordinates, and then convert to cartesian
            symbols = []
            coords = []
            for frac_coords, index in zip(shifted_frac_coords, molecule):
                frac_shifted_entirely = frac_coords + enitre_mol_shift
                cart_coords = structure.lattice.get_cartesian_coords(frac_shifted_entirely)
                symbol = structure[index].species.elements[0].symbol
                symbols.append(symbol)
                coords.append(cart_coords)

            mol = Molecule(symbols, coords)
            molecules.append(mol)

        # Based on the structural formulas, create a list of lists, with each list containing all symmetries of a spcaer
        unique_formulas = {}
        for molecule in molecules:
            formula = str(molecule.composition.formula)
            if formula not in unique_formulas:
                unique_formulas[formula] = []
            unique_formulas[formula].append(molecule)

        all_molecules = list(unique_formulas.values())
        self.molecules = all_molecules
        return self.molecules

    def reorder_to_reference(self, ref_mol: Molecule, target_mol: Molecule) -> Molecule:
        """
        This function reorders the atoms of the target molecule (all the symmetry configurations of a unique spacer)
        to match the ordering of the reference molecule (the first symmetry configuration of a unique spacer). Then, the
        atoms of the target molecule are all reordered (using graph isomorphism) to match the ordering of the atoms of
        the reference molecule. EACH symmetry copy needs to have the same ordering of atoms to prevent errors in later
        steps (e.g. when writing to mol2 files or performing parametrisation).

        :param ref_mol: The reference molecule, used to determine the correct ordering of atoms for the target molecule.
        :param target_mol: The target molecule, which atoms will be reordered according to the reference molecule.
        :return: A Molecule instance with atoms that are reordered based on the reference molecule.
        """
        # Build graphs for both, the target and reference molecule(s)
        molecule_graph_ref = self._molecule_to_graph(ref_mol)
        molecule_graph_target = self._molecule_to_graph(target_mol)

        # Define nodes matching based on the element types
        matching_nodes = nx.algorithms.isomorphism.categorical_node_match("element", None)
        matching_graphs = nx.isomorphism.GraphMatcher(molecule_graph_ref, molecule_graph_target, node_match=matching_nodes)
        if not matching_graphs.is_isomorphic():
            print("\n" + "="*40)
            print("--- ISOMORFISM ERROR DEBUG ---")
            print(f"Reference: {len(ref_mol)} atoms, {molecule_graph_ref.number_of_edges()} bonds")
            print(f"Copy:      {len(target_mol)} atoms, {molecule_graph_target.number_of_edges()} bonds")
            
            ref_counts = {}
            for e in [s.symbol for s in ref_mol.species]: ref_counts[e] = ref_counts.get(e, 0) + 1
            tgt_counts = {}
            for e in [s.symbol for s in target_mol.species]: tgt_counts[e] = tgt_counts.get(e, 0) + 1
            
            print(f"Elements Ref:   {ref_counts}")
            print(f"Elements Copy: {tgt_counts}")
            print("="*40 + "\n")
            raise RuntimeError("No isomorphism found between symmetry copies")

        # Reorder target atoms according to reference order
        mapping = matching_graphs.mapping  # Maps the indices of target_mol atoms to ref_mol atoms
        new_symbols = [str(target_mol.species[mapping[i]]) for i in range(len(ref_mol))]
        new_coords = [target_mol.cart_coords[mapping[i]] for i in range(len(ref_mol))]
        return Molecule(new_symbols, new_coords)

    def _molecule_to_graph(self, mol: Molecule) -> nx.Graph:
        """
        Converts a pymatgen Molecule instance into a graph to be later used for reordering atoms based on graph
        isomorphism. 
        """
        single_molecule_graph = nx.Graph()
        
        for i, elem in enumerate(mol.species):
            single_molecule_graph.add_node(i, element=str(elem))
            
        for i in range(len(mol)):
            for j in range(i+1, len(mol)):
                dist = np.linalg.norm(mol.cart_coords[i] - mol.cart_coords[j])
                atom_i = mol.species[i].symbol
                atom_j = mol.species[j].symbol
                
                if atom_i in self.bond_radius and atom_j in self.bond_radius:
                    cutoff = self.bond_radius[atom_i] + self.bond_radius[atom_j] + self.bond_tolerance
                else:
                    cutoff = 2.0  # Unknowns
                    
                if dist <= cutoff:
                    single_molecule_graph.add_edge(i, j)
                    
        return single_molecule_graph

    def reorder_all_spacer_atom_indices(self, all_molecules: List[List[Molecule]]) -> List[List[Molecule]] | bool:
        """
        This function reorders all found spacer molecules (all symmetries of each unique spacer) to have the same atom
        ordering as the first symmetry of each unique spacer. This is needed for later steps in the pipeline.

        :param all_molecules: A list within a list, where each outer list contains an inner list of all spacer molecules
        :return: A list within a list, where each outer list contains an inner list of all spacer molecule symmetries
        where each symmetry is reordered to have the same atom ordering; returns False if reordering fails.
        """
        reordered_groups = []
        try:
            for group in all_molecules:
                ref_molecule = group[0]  # Use first molecule as reference
                new_group = [ref_molecule]  # Start new group with reference molecule (because it's already correct)

                for molecule_symmetry in group[1:]:
                    new_group.append(self.reorder_to_reference(ref_molecule, molecule_symmetry))
                reordered_groups.append(new_group)
            return reordered_groups

        except RuntimeError as e: print(f"Reordering failed: {e}"); return False

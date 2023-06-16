"""
This is the code to calculate the CM, CM_star and Cse metrics.
"""
from collections import Counter
from typing import Dict, Iterator, List, NamedTuple, Tuple

import numpy as np
from rdkit import Chem

# Atom types are defined as (symbol, total degree, non-h degree)
AtomType = Tuple[str, int, int]
# Atoms are defined as (atom index, atom type)
Atom = Tuple[int, AtomType]
AtomDict = Dict[Atom, List[Atom]]


class MolecularComplexity(NamedTuple):
    cm: float
    cm_star: float
    cse: float


def _non_h_items(data: Dict[Atom, any]) -> Iterator[Tuple[Atom, any]]:
    """
    Generator for non-H items from a dictionary where the keys are atom tuples.
    """
    for key, val in data.items():
        if key[1][0] != "H":
            yield key, val


def _collect_atom_paths(neighbors: AtomDict) -> List[List[tuple]]:
    """
    Returns list of atom paths for each atom.

    An atom path is a tuple of atom types.
    """
    atom_paths = []
    for atom, nbs in _non_h_items(neighbors):
        paths = []
        for nb in nbs:
            if nb[1][0] == "H" or neighbors[nb] == [atom]:
                # No second neighbors
                paths.append((atom[1], nb[1]))
            else:
                paths.extend(
                    (atom[1], nb[1], nb2[1]) for nb2 in neighbors[nb] if nb2 != atom
                )

        atom_paths.append(paths)

    return atom_paths


def get_atom_type(atom: Chem.rdchem.Mol) -> AtomType:
    """
    Return a tuple describing the atom type.

    Considers element, total number of connections, and number of non-H connections.
    """
    symbol = atom.GetSymbol()
    degree = atom.GetTotalDegree()
    h_count = atom.GetTotalNumHs(includeNeighbors=True)
    non_h = degree - h_count
    return (symbol, degree, non_h)


def fractional_occurrence(data: list) -> np.ndarray:
    """
    Calculate the fractional occurrence of unique items in the input data.

    Uniqueness determined by collections.Counter.

    Returns:
        np.ndarray: fractional occurrence of unique items
    """
    counter = Counter(data)
    counts = np.array(list(counter.values()))
    return counts / len(data)


def calculate_molecular_complexity(mol: Chem.rdchem.Mol) -> MolecularComplexity:
    """
    This is a function to calculate the molecular complexity metrics described in
    Proudfoot, Bioorganic & Medicinal Chemistry Letters 27 (2017) 2014-2017.
    https://doi.org/10.1016/j.bmcl.2017.03.008

    This function takes an RDKit Mol object, enumerates atom paths,
    and then calculates the complexity environment for each atom CA as

    CA = - Sum (pi*log2(pi)) + log2(N)

    where pi is the fractional occurrence of each path type emanating from
    an atom and N is the total number of paths emanating from that atom.

    Molecular complexity CM can be defined as either the simple sum of the CA,
    or CM* which is the log-sum of the exponentials of the CA.

    CM = Sum (CA)

    CM* = log2(Sum (2**CA))

    Cse = - Sum (qi*log2(qi))

    where qi is the fractional occurrence of an atom environment.
    """
    # get atom types for each atom in the molecule
    atoms = [(atom.GetIdx(), get_atom_type(atom)) for atom in mol.GetAtoms()]

    # create dict with neighbors of each atom
    neighbors = {
        atom: [
            atoms[neighbor.GetIdx()]
            for neighbor in mol.GetAtomWithIdx(atom[0]).GetNeighbors()
        ]
        for atom in atoms
    }

    atom_paths = _collect_atom_paths(neighbors)

    cas = np.zeros(len(atom_paths))
    for i, paths in enumerate(atom_paths):
        total_paths = len(paths)
        pi = fractional_occurrence(paths)
        cas[i] = -np.sum(pi * np.log2(pi)) + np.log2(total_paths)

    cm = np.sum(cas)

    cm_star = np.log2(np.sum(2**cas))

    # sort and concatenate the individual paths to compare the atom environments
    atom_environments = [tuple(sorted(paths)) for paths in atom_paths]

    # Now we can calculate the Cse metric as the fractional occurrence of each atom environment
    qi = fractional_occurrence(atom_environments)
    cse = -np.sum(qi * np.log2(qi))

    return MolecularComplexity(cm, cm_star, cse)


def molecular_complexity(smiles: str) -> MolecularComplexity:
    """This function takes SMILES and returns the CM, CM*, and Cse metrics"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        return calculate_molecular_complexity(mol)
    except Exception:
        return MolecularComplexity(np.nan, np.nan, np.nan)

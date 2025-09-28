#!/usr/bin/env python

"""
Sequence-based structural alignment of two proteins.
"""

import os
from pathlib import Path

import requests
import numpy as np
from Bio.Align import substitution_matrices, PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import PDBParser, FastMMCIFParser, Superimposer, PDBIO, MMCIFIO, Structure, Select, Model, Chain
from Bio.PDB.PDBList import PDBList
from Bio.PDB.Polypeptide import is_aa


class AlignCystalStructures:

    SUBST_MATRICES = list(map(str.lower, substitution_matrices.load()))

    @staticmethod
    def by_sequence(struct1: Structure,
                    struct2: Structure,
                    matrix: str = 'blosum62',
                    gap_open: float = -10.0,
                    gap_extend: float = -0.5,
                    mapping: bool = False) -> tuple[str, str] | tuple[str, str, dict[int, int]]:
        f"""Align the protein sequences of two protein structures using the Needleman-Wunsch algorithm.

        :param struct1: First structure whose sequence to align to that of `struct2`.
        :param struct2: Second structure whose sequence to align to that of `struct1`.
        :param matrix: Substitution penality matrix (one of {AlignCystalStructures.SUBST_MATRICES}).
        :param gap_open: Gap opening value.
        :param gap_extend: Gap extending value.
        :param mapping: should the indices of aligned residues be returned
        :return: the aligned sequences, and if `mapping` is True, the amino acid mapping indices
        """
        if not matrix in AlignCystalStructures.SUBST_MATRICES:
            raise ValueError('Invalid substitution penality matrix')
        matrix = substitution_matrices.load(matrix.upper())
        # Parse sequences with amino acid indices
        resseq1 = extract_pdb_sequence(struct1, with_resseq=True)
        resseq2 = extract_pdb_sequence(struct2, with_resseq=True)
        seq1 = ''.join(list(zip(*resseq1))[1])
        seq2 = ''.join(list(zip(*resseq2))[1])
        # Align sequences
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.substitution_matrix = matrix
        aligner.open_gap_score = gap_open
        aligner.extend_gap_score = gap_extend
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        alns = aligner.align(seq1, seq2)
        aligned1, aligned2 = next(alns)
        # Renumber residues relative to the reference
        if not mapping:
            return aligned1, aligned2
        mapping_ = {}
        new_index1, new_index2 = 0, 0
        for i, (aligned_aa1, aligned_aa2) in enumerate(zip(aligned1, aligned2)):
            if aligned_aa1 == '-' and aligned_aa2 != '-':
                    new_index2 += 1
            elif aligned_aa2 == '-' and aligned_aa1 != '-':
                    new_index1 += 1
            elif aligned_aa1 == '-' and aligned_aa2 == '-':
                    raise ValueError(f'Both aligned sequences contain a gap at index {i}.')
            else:
                assert resseq1[new_index1][1] == aligned_aa1
                assert resseq2[new_index2][1] == aligned_aa2
                mapping_[resseq1[new_index1][0]] = resseq2[new_index2][0]
                new_index1 += 1
                new_index2 += 1
        return aligned1, aligned2, mapping_

    @staticmethod
    def by_coordinates(reference: Structure,
                       mobile: Structure,
                       ref_chain: str = 'A',
                       mobile_chain: str = 'A') -> tuple[str, float]:
        """Align two protein structures on their alpha carbons.

        :param reference: the protein structure to superimpose `mobile` onto
        :param mobile: the protein structure to be superimposed onto `reference`
        :param ref_chain: the chain of `reference` to consider for superimposition
        :param mobile_chain: the chain of `mobile` to consider for superimposition
        :returns: the filename of the aligned protein structure and the RMSD between the aligned structures
        """
        # Parse structures & take only the necessary chain
        try:
            reference = reference[0][ref_chain]
        except KeyError:
            raise Exception('Chain {0} not found in reference structure'.format(ref_chain))
        try:
            mobile_struct = mobile[0][mobile_chain]
        except KeyError:
            raise Exception('Chain {0} not found in mobile structure'.format(mobile_chain))
        # Align sequences to get mapping between residues
        _, _, res_map = AlignCystalStructures.by_sequence(reference, mobile_struct, mapping=True)
        # Identify alpha carbons of residues that align
        ref_ca_list, mobile_ca_list = [], []
        for ref_res in res_map:
            ref_ca_list.append(reference[ref_res]['CA'])
            mobile_ca_list.append(mobile_struct[res_map[ref_res]]['CA'])
        # Superimpose matching residues
        si = Superimposer()
        si.set_atoms(ref_ca_list, mobile_ca_list)
        # Transform & Write Mobile
        si.apply(mobile_struct.get_atoms())
        rmsd = si.rms
        io = PDBIO()
        io.set_structure(mobile_struct)
        m_transformed_name = '{0}_transformed.pdb'.format(mobile.id)
        io.save(m_transformed_name)
        return m_transformed_name, rmsd


def aligned_sequences_identity(seq1: str, seq2: str) -> float:
    """Return the percentage of identical characters between two already aligned sequences."""
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must have same length')
    matches = [x == y for x, y in zip(seq1, seq2)]
    seq_id = (100 * sum(matches)) / len( seq1)
    return seq_id


def extract_pdb_sequence(structure: Structure, with_resseq: bool = False) -> str | list[tuple[int, str]]:
    """Retrieve the 1 letter protein sequence from a PDB structure.

    :param structure: the PDB structure to extract sequence from
    :param with_resseq: should the amino acid indices be returned
    :returns: the protein sequence if `with_resseq` is False, otherwise a list of tuples
    whose first items are amino acid indices and whose second items are the one letter amino acid codes.
    """
    oneletter = lambda res: (res.id[1], protein_letters_3to1.get(res.resname, 'X'))
    seq = [oneletter(res) for res in structure.get_residues() if is_aa(res)]
    if not with_resseq:
        return ''.join(list(zip(*seq))[1])
    return seq


def parse_structure(filepath: str) -> Structure:
    """Parse a crystal structure from a PDB or CIF file."""
    if not os.path.isfile(filepath):
        return IOError('File not found: {0}'.format(filepath))
    if filepath.endswith(('pdb', 'ent')):
        parser = PDBParser()
    elif filepath.endswith('cif'):
        parser = FastMMCIFParser()
    else:
        raise Exception('Format not supported ({0}). Must be .pdb/.ent or .cif'.format(filepath))
    basename = Path(filepath).stem
    return parser.get_structure(basename, filepath)


def write_structure(structure: Structure, filepath: str, pdb_format: str = 'pdb') -> None:
    """Write a structure to a file."""
    if pdb_format not in ('pdb', 'ent', 'cif'):
        raise ValueError('format must be either "pdb" or "cif".')
    if pdb_format in ('pdb', 'ent'):
        io = PDBIO()
    else:
        io = MMCIFIO()
    io.set_structure(structure)
    io.save(filepath)


def download_structure(pdbid: str, mirror: str = 'wwPDB', pdb_format: str = 'pdb', chain_id: str = None, ligand_to_keep: str | list[str] = None) -> str:
    """Download a structure as an mmCIF file from one mirror of the PDB in the current folder.

    :param pdbid: PDB ID of the structure to download.
    :param mirror: Mirror of the PDB ('PDB', 'wwPDB', 'PDBe', or 'PDBj'). Is case-insensitive.
    :param pdb_format:
    :param chain_id: chain to keep (all others are removed). Set to None to keep all.
    :param ligand_to_keep: All ligands but the one(s) specified are dropped from the structure. If None, keep all.
    """
    if pdb_format not in ['pdb', 'cif']:
        raise ValueError('format must be either "pdb" or "cif".')
    mirror = mirror.lower()
    if mirror == 'wwpdb':
        pdb = PDBList(server='https://files.wwpdb.org', verbose=False)
        outname = pdb.retrieve_pdb_file(pdbid, file_format=('mmCif' if pdb_format == 'cif' else 'pdb'), pdir='.')
        _ = os.rename(outname, f'./{pdbid}.{pdb_format}')
        url = None
    elif mirror == 'pdbe':
        url = "https://www.ebi.ac.uk/pdbe/entry-files/download/" + pdbid + f'.{pdb_format}'
    elif mirror == 'pdbj':
        if pdb_format == 'cif':
            url = "https://data.pdbj.org/pub/pdb/data/structures/divided/mmCIF/" + pdbid[1:3] +"/" + pdbid + ".cif"
        else:
            url = "https://data.pdbj.org/pub/pdb/data/structures/divided/pdb/" + pdbid[1:3] +"/pdb" + pdbid + ".ent"
    elif mirror in ['pdb', 'rcsb']:
        url = "http://www.rcsb.org/pdb/files/" + pdbid + f'.{pdb_format}'
    if url is not None:
        content = requests.get(url).text
        with open(f'./{pdbid}.cif', 'w') as oh:
            oh.write(content)
    if chain_id is not None:
        # Read the file
        structure = parse_structure(f'./{pdbid}.{pdb_format}')
        structure = extract_chain_with_ligands(structure, chain_id=chain_id, ligand_to_keep=ligand_to_keep)
        write_structure(structure=structure, filepath=f'./{pdbid}.{pdb_format}', pdb_format=pdb_format)
    return f'./{pdbid}.{pdb_format}'


def get_min_distance(chain, ligand):
    """Calculates the minimum distance between a ligand and a protein chain."""
    min_dist = float('inf')
    for residue in chain:
        for atom in residue:
            for ligand_atom in ligand:
                dist = np.linalg.norm(atom.coord - ligand_atom.coord)
                if dist < min_dist:
                    min_dist = dist
    return min_dist


def extract_chain_with_ligands(structure, chain_id, ligand_to_keep):
    """
    Extracts a specific chain and its associated ligands from a BioPython Structure.

    Args:
        structure (Bio.PDB.Structure.Structure): The input BioPython Structure object.
        chain_id (str): The ID of the chain to extract.
        ligand_to_keep (str, optional): The name of the ligands, lying within 2 angstroms of the protein, that are kept.

    Returns:
        Bio.PDB.Structure.Structure: A new Structure object containing the desired
                                     chain and its ligands.
    """
    new_structure = Structure.Structure("extracted_chain")
    new_model = Model.Model(0)
    new_chain = Chain.Chain(chain_id)
    target_chain_residues = []
    ligands_to_add = set()
    for model in structure:
        # First, get all residues of the target protein chain
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    # We are only interested in standard amino acids for the protein part
                    if residue.id[0] == ' ':
                        target_chain_residues.append(residue.copy())
        # Now, find all ligands in the original structure
        for chain in model:
            for ligand in chain:
                if ligand_to_keep is None:
                    ligands_to_add.add(ligand)
                elif ligand.id[0] == f'H_{ligand_to_keep}':
                    # Calculate the minimum distance between the ligand and the target chain
                    if get_min_distance(target_chain_residues, ligand) <= 4:
                        ligands_to_add.add(ligand)
    # Add the protein residues and the associated ligands to the new chain
    for res in target_chain_residues:
        new_chain.add(res)
    for lig in ligands_to_add:
        try:
            new_chain.add(lig)
        except:
            pass
    new_model.add(new_chain)
    new_structure.add(new_model)
    return new_structure


class ResSelect(Select):

    def __init__(self, ligand):
        self.ligand_id = ligand

    def accept_residue(self, residue):
        if residue.get_resname() == self.ligand_id:
            return 1
        else:
            return 0

class NonHetSelect(Select):

    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

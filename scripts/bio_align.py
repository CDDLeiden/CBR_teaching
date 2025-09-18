#!/usr/bin/env python

"""
Sequence-based structural alignment of two proteins.
"""

import os
from pathlib import Path

import requests
from Bio.Align import substitution_matrices, PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import PDBParser, FastMMCIFParser, Superimposer, PDBIO, Structure, Select
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


def download_structure(pdbid: str, mirror: str = 'wwPDB') -> str:
    """Download a structure as an mmCIF file from one mirror of the PDB in the current folder.

    :param pdbid: PDB ID of the structure to download.
    :param mirror: Mirror of the PDB ('PDB', 'wwPDB', 'PDBe', or 'PDBj'). Is case-insensitive.
    """
    mirror = mirror.lower()
    if mirror == 'wwpdb':
        pdb = PDBList(server='https://files.wwpdb.org', verbose=False)
        pdb.retrieve_pdb_file(pdbid, file_format='mmCif', pdir='.')
        return f'./{pdbid}.cif'
    elif mirror == 'pdbe':
        url = "https://www.ebi.ac.uk/pdbe/entry-files/download/" + pdbid + ".cif"
    elif mirror == 'pdbj':
        url = "https://data.pdbj.org/pub/pdb/data/structures/divided/mmCIF/" + pdbid[1:3] +"/" + pdbid + ".cif"
    elif mirror in ['pdb', 'rcsb']:
        url = "http://www.rcsb.org/pdb/files/" + pdbid + ".cif"
    content = requests.get(url).text
    with open(f'./{pdbid}.cif', 'w') as oh:
        oh.write(content)
    return f'./{pdbid}.cif'


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

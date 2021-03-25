from Bio import SeqIO
from biofile import getaccnum
from helpers import save_file, make_int, make_set
from pathlib import Path
import itertools
import numpy as np

def main():
    ##################################################################
    # From files of each phage caller extract coordinates of viruses #
    # and viral contigs with their accesion numbers.                 #
    ##################################################################

    # Number of nucleotides that detections will be extended by.
    # Default 2kb (2000 nucleotides).
    extend_value = 2000
    # Extract prophages' coordinates from phispy.
    phispy_out = Path(snakemake.input[0], 'prophage_coordinates.tsv')
    phispy = parse_phispy(phispy_out)

    # Extract prophages' coordinates and viral contigs from virsorter.
    virsorter_out = Path(snakemake.input[1], 'Predicted_viral_sequences')
    virsorter = parse_virsorter_prophages(virsorter_out)

    # Contigs feature to be added.
    ###################################################################
    # Extract accesion numbers of whole viral contigs.                #
    # Convert them to coordinates with accesion numbers.              #
    # virsorter_contigs = parse_virsorter_viral_contigs(virsorter_out)#
    # virsorter_contigs = contig_to_coordinates(virsorter_contigs)    #
    # virsorter = virsorter + virsorter_contigs                       #
    ###################################################################

    # Get all unique accesion numbers.
    unique_accnums = get_unique_accnum(phispy, virsorter)
    phispy, virsorter = list(zip(*phispy)), list(zip(*virsorter))
    print(unique_accnums)

    for unique_accnum in unique_accnums:
        odata = [phispy, virsorter, unique_accnum]
        ostarts, oends = union_of_overlapping_detections(*odata)
        ostarts, oends = collapse(ostarts, oends)

        nodata = [phispy, virsorter, unique_accnum, ostarts, oends]
        nostarts, noends = nonoverlapping_detections(*nodata)

        starts, ends = ostarts + nostarts, oends + noends
        starts, ends = extend(starts, ends, extend_value)

        lines = ['\t'.join([str(start), str(end), unique_accnum])
        for start, end in zip(starts, ends)]
        lines = '\n'.join(lines)

        save_file(snakemake.output[0], lines)


def parse_phispy(table):
    """Input: phispy table.
    Output: if contigs: coordinates with contigs ids.
    Output: if complete: coordinates and empty list of contig ids."""
    coordinates = []
    try:
        with open(table, 'r') as f:
            phages = [line.strip().split()[1:4] for line in f.readlines()]
            accesion_numbers, starts, ends = list(zip(*phages))
    except: pass

    starts, ends = make_int(starts), make_int(ends)
    accesion_numbers = [getaccnum(accnum) for accnum in accesion_numbers]
    return starts, ends, accesion_numbers


def parse_virsorter_prophages(fasta_files_path):
    """Input: Path to Predicted_viral_sequences of virsorter.
    Output: coordinates associated conigs"""
    # Prophages on complete genomes.
    phages = list(Path(fasta_files_path).glob('*[54].fasta'))

    phage_files = []
    accesion_numbers = []
    for phage in phages:
        with open(phage, 'r') as file:
            _ = [phage_files.append(line.split('-')[-3:-1]) for line
            in file.readlines() if line[0] == '>']

            file.seek(0)

            _ = [accesion_numbers.append(line) for line
            in file.readlines() if line[0] == '>']

    starts, ends = list(zip(*phage_files))
    starts, ends = make_int(starts), make_int(ends)
    accesion_numbers = [getaccnum(accnum) for accnum in accesion_numbers]
    return starts, ends, accesion_numbers


def parse_virsorter_viral_contigs(fasta_files_path):
    """Input: Path to Predicted_viral_sequences of virsorter.
    Output: list of contig's accesion numbers that are viral."""
    # Whole viral contigs. Categories of virsorter: one and two.
    contigs = list(Path(fasta_files_path).glob('*[12].fasta'))

    # Whole viral contigs.
    contigs_accesion_numbers = []
    for contig in contigs:
        with open(contig, 'r') as file:
            [contigs_accesion_numbers.append(getaccnum(line)) for line
            in file.readlines() if line[0] == '>']

    return contigs_accesion_numbers


def get_unique_accnum(*args):
    accnumbers = [phagecaller[-1] for phagecaller in args if phagecaller]
    accnumbers = list(itertools.chain.from_iterable(accnumbers))
    unique_accnum = set(accnumbers)
    return list(unique_accnum)


def contig_to_coordinates(genome, accesion_numbers):
    ##################################
    # Function needs to be evaluated #
    ##################################
    """Input: contig accesion number
    Output: contig coordaintes (from zero to it's length)"""
    lengths = []
    for accesion_number in accesion_numbers:
        records = SeqIO.parse(genome, 'genbank')
        for record in records:
            record_id = getaccnum(record.id)
            if record_id == accesion_number:
                lengths.append(len(record.seq))

    starts = [0] * len(accesion_numbers)
    ends = lengths
    return starts, ends, accesion_numbers


def get_matrix(phispy, virsorter):
    matrix = [[len(set.intersection(p,v)) for v in virsorter] for p in phispy]
    matrix = np.array(matrix)
    return matrix


def overlapping_indicies(phispy, virsorter, matrix):
    if len(np.nonzero(matrix)) < 2:
        phi_indicies = [0] * len(np.nonzero(matrix)[0])
        vrs_indicies = np.nonzero(matrix)[0].tolist()
    else:
        phi_indicies, vrs_indicies = np.nonzero(matrix)
    return phi_indicies, vrs_indicies


def filter_and_convert(program, unique_accesion_number):
    """Filters detecions for the given unique_accession_number and converts
    coordinates to sets."""
    program = [make_set(start, end) for start, end, accnum in program
    if accnum == unique_accesion_number]
    return program


def collapse(starts, ends):
    """Get unions of all overlapping detections. """
    unions = [make_set(start, end) for start, end in zip(starts, ends)]

    merged = []
    # indicies = []
    for element1 in unions:
        vector = []
        for element2 in unions:
            if set.intersection(element1, element2):
                vector.append(unions.index(element2))
        to_merge = [unions[index] for index in vector]
        merged_overlaping = set.union(*to_merge)
        merged.append(merged_overlaping)

    starts = [min(union) for union in merged]
    ends = [max(union) for union in merged]

    starts = list(dict.fromkeys(starts))
    ends = list(dict.fromkeys(ends))
    return starts, ends


def union_of_overlapping_detections(phispy, virsorter, unique_accesion_number):
    """Get union of overlapping detections from phispy and virsorter."""
    phispy = filter_and_convert(phispy, unique_accesion_number)
    virsorter = filter_and_convert(virsorter, unique_accesion_number)
    matrix = get_matrix(phispy, virsorter)
    phi_indicies, vrs_indicies = overlapping_indicies(phispy, virsorter, matrix)

    unions = [set.union(phispy[i_phi], virsorter[i_vrs]) for i_phi, i_vrs
    in zip(phi_indicies, vrs_indicies)]

    starts = [min(union) for union in unions]
    ends = [max(union) for union in unions]
    return starts, ends


def nonoverlapping_detections \
(phispy, virsorter, unique_accesion_number, ostarts, oends):
    phispy = filter_and_convert(phispy, unique_accesion_number)
    virsorter = filter_and_convert(virsorter, unique_accesion_number)
    union_of_overlapping = [make_set(start, end) for start, end
                            in zip(ostarts, oends)]

    union_of_overlapping = set.union(*union_of_overlapping)

    lonely_detections = []
    for detection in phispy:
        if not set.intersection(detection, union_of_overlapping):
            lonely_detections.append(detection)

    for detection in virsorter:
        if not set.intersection(detection, union_of_overlapping):
            lonely_detections.append(detection)

    starts = [min(union) for union in lonely_detections]
    ends = [max(union) for union in lonely_detections]
    return starts, ends



def extend(starts, ends, extend_value):
    ######################################
    # Take care of start and end of contigs (do not extend beyond).
    ######################################
    starts = [start-extend_value for start in starts]
    ends = [end+extend_value for end in ends]
    return starts, ends



if '__main__' == __name__:
    main()

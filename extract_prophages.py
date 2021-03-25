from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import itertools

def main():
    ####################################################################
    # Handle contigs by trying to open all unique accesion numbers.    #
    # When file loaded iterate over records to find contig.            #
    ####################################################################

    coordinates = load_coordinates(snakemake.input[1])
    unique_accnums = get_unique_accnum(coordinates)

    ########################################
    # Find approriate genome and load it.  #
    ########################################
    for unique_accnum in unique_accnums:
        try: records = load_genome(snakemake.input[0])
        except: continue

    ############################################
    # Load prophage coordinates of the genome. #
    ############################################
    for unique_accnum in unique_accnums:
        filt_coordinates = [(start, end) for start, end, accnum in coordinates
        if accnum == unique_accnum]

        seq_prophages = []
        accnums = []
        for record in records:
            if unique_accnum in record.id:
                for start, end in filt_coordinates:
                    seq_prophages.append(record.seq[start:end+1])
                    accnums.append(unique_accnum)

    records = []
    for i, (seqobject, accnum) in enumerate(zip(seq_prophages, accnums)):
        record = SeqRecord(seqobject, id=accnum, name=f'pp_{i}',
                           description=f'putative_prophage_{i}')
        # records.append(record)

        output = Path(Path.cwd(), snakemake.output[0], f'pp_{i}_{accnum}.fna')
        folder = output.parent.mkdir(exist_ok=True)
        count = SeqIO.write(record, output, 'fasta')


def load_genome(fpath):
    return list(SeqIO.parse(fpath, 'genbank'))


def load_coordinates(fpath):
    with open(fpath) as f:
        lines = [line.strip().split() for line in f.readlines()]
    lines = [(int(start), int(end), accnum) for start, end, accnum in lines]
    return lines

def get_unique_accnum(*args):
    accnumbers = [phagecaller[-1] for phagecaller in args[0] if phagecaller]
    unique_accnum = set(accnumbers)
    return list(unique_accnum)


if __name__ == "__main__":
    main()

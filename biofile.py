import re
import string
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio import Entrez
import random

def getaccnum(anystring):
    """Takes any string and extracts by regex an accesion number."""

    ################################################################
    # Only refseq accesion numbers.                                #
    # ?Add genbank accesion numbers (underscore at third postion)? #
    ################################################################
    upper_letters = string.ascii_uppercase
    digits = string.digits
    pattern = rf'[{upper_letters}]{{1,2}}[{digits}]{{5,8}}'
    accnum = re.search(pattern, anystring)
    if accnum:
        accnum = accnum.group()
        return accnum
    else:
        print('Failed to find accesion number! Random string returned.')
        return ''.join(random.choice(string.ascii_lowercase) for x in range(8))


def get_genome_lenght(genome, ftype='fasta'):
    """Get lenght of the analyzed nucleotide sequence."""
    seq_record = SeqIO.parse(genome, ftype)
    sequence = next(seq_record).seq
    return len(sequence)

def name_as_accnum(fpath, **kwargs):
    """Input: genbank file.
    Output: the same genbank file in the same directory,
    renamed to accession number.
    Output path can be defined."""
    try:
        records = SeqIO.parse(fpath, 'genbank')
    except:
        print('FileNameError\n')
        exit()

    records = list(records)
    accnum = getaccnum(records[0].name)
    if kwargs:
        new_name = Path(kwargs['output'], accnum + '.gbk')
    else:
        new_name = Path(Path(fpath).parent, accnum + '.gbk')

    count = SeqIO.write(records, new_name, 'genbank')
    print(f'Converted {count} record(s).')
    print('Genbank successfully renamed!')


def genbank_to_fasta(gbkpath, **kwargs):
    """Input: genbank file.
    Output: fasta format with appropriate headers."""
    try:
        records = SeqIO.parse(gbkpath, 'genbank')
    except:
        print('FileNameError\n')
        exit()

    records = list(records)
    accnum = getaccnum(records[0].name)

    if kwargs:
        new_name = Path(kwargs['output'], accnum + '.fasta')
    else:
        new_name = Path(Path(gbkpath).parent, accnum + '.fasta')

    SeqIO.convert(gbkpath, "genbank", new_name, "fasta")
    print('Genbank successfully converted!')


def get_cds(genome, coordinates):
    """Get CDS's from genbank at given coordinates. """
    gb = f'/Users/januszkoszucki/{genome}.gb'
    start, end = coordinates

    cds = []
    for rec in SeqIO.parse(gb, "genbank"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    if start < feature.location.start and \
                        end > feature.location.end:
                        if 'product' in list(feature.qualifiers.keys()):
                            cds.append(feature.qualifiers["product"][0])
    return len(cds)


def download_genbank(ids, names):
    Entrez.email = 'j.koszucki@gmail.com'
    cwd = Path.cwd()
    ext = '.gb'

    for id, name in zip(ids,names):
        fname = Path(cwd, name)
        if not fname.exists():
            try:
                record = Entrez.efetch(db='nucleotide',
                                            id=id,
                                            rettype='gbwithparts',
                                            retmode='text')
            except:
                continue


            with open(str(fname) + ext, 'w+') as f:
                _ = [f.write(line) for line in record.readlines()]
        else:
            continue

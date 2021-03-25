from pathlib import Path
from biofile import name_as_accnum, genbank_to_fasta

raw = Path(Path.cwd(), 'data', 'raw').glob('*[gbk|genbank|gb]')

#############################
# ADD ANOTATION WITH PATRIC #
#############################

genbank_out = Path(Path.cwd(),'data', 'processed', 'genbank')
_ = [name_as_accnum(file, output=genbank_out) for file in raw]

gbkfiles = genbank_out.glob('*gbk')
fasta_out = Path(Path.cwd(), 'data', 'processed', 'fasta')
_ = [genbank_to_fasta(gbkfile, output=fasta_out) for gbkfile in gbkfiles]

gbkfiles = genbank_out.glob('*gbk')
fnames = [file.stem for file in gbkfiles]
print('\nUseful input for config file. File names. \n')
print('\n'.join(fnames))

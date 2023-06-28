import functools
import gzip
from pathlib import Path


annotation_path = Path(snakemake.input['annotation'])
output_path = Path(snakemake.output['tx2gene'])

if annotation_path.suffix == '.gz':
    opener = gzip.open
else:
    opener = open

with opener(annotation_path, mode = 'rt') as annotation_file, open(output_path, mode = 'w') as output_file:
    for line in annotation_file:
        if not line.startswith('#'):
            seqid, source, feature_type, start, end, score, strand, phase, attributes = line.rstrip().split("\t")
            if feature_type == 'mRNA':
                attributes = {k: v for k, v in (kv.split('=') for kv in attributes.split(';'))}
                output_file.write(f"{attributes['ID']}\t{attributes['Parent']}\n")

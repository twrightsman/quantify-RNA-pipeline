#!/usr/bin/env python3

from collections import defaultdict
import hashlib
import logging
from pathlib import Path
import sys
from typing import Optional

import Bio.Seq
import Bio.SeqRecord
import Bio.SeqIO
import gffutils
from pysam import FastaFile


def create_unique_feature_id(feature: gffutils.Feature) -> str:
    if (feature.featuretype == 'gene') or (feature.featuretype == 'mRNA'):
        # don't mess up names of potential parents
        return feature.attributes['ID'][0]
    else:
        # features without children can have autoincremented IDs
        return f"autoincrement:{feature.attributes['ID'][0]}"


def main(reference_path: Path, annotation_path: Path, output_path: Path, gffutils_cache_path: Optional[Path]):
    reference = FastaFile(str(reference_path))

    if gffutils_cache_path is not None:
        # cache database to disk
        gffutils_cache_db_path = gffutils_cache_path / str(hashlib.sha256(str(annotation_path).encode('utf-8')).hexdigest())
        logging.debug("using '%s' as gffutils cache path", gffutils_cache_db_path)
        # load cached gffutils database if previously generated
        if gffutils_cache_db_path.exists():
            annotation = gffutils.FeatureDB(str(gffutils_cache_db_path))
        else:
            annotation = gffutils.create_db(
                str(annotation_path),
                str(gffutils_cache_db_path),
                id_spec = create_unique_feature_id,
                merge_strategy = 'error'
            )
    else:
        # do not cache to disk, create fresh database in memory
        annotation = gffutils.create_db(
            str(annotation_path),
            ':memory:',
            id_spec = create_unique_feature_id,
            merge_strategy = 'error'
        )

    transcripts = []
    for gene in annotation.features_of_type('gene'):
        for transcript in annotation.children(gene, featuretype='mRNA'):
            exons = []
            spliced_sequence = []
            # extract exonic sequence
            for exon in annotation.children(transcript, featuretype='exon', order_by='start', reverse=(transcript.strand == '-')):
                sequence = Bio.Seq.Seq(reference.fetch(reference = exon.seqid, start = exon.start - 1, end = exon.end))
                spliced_sequence.append(sequence.reverse_complement() if transcript.strand == '-' else sequence)
                exons.append(exon)
            # join exons
            if spliced_sequence and transcript.id:
                mRNA = Bio.SeqRecord.SeqRecord(
                    seq = Bio.Seq.Seq('').join(spliced_sequence),
                    id = transcript.id,
                      name = '',
                      description = ''
                )
                if all((len(s) > 0 for s in spliced_sequence)):
                    transcripts.append(mRNA)

    if transcripts:
        with open(output_path, 'w') as output_file:
            Bio.SeqIO.write(transcripts, output_file, 'fasta')


if __name__ == "__main__":
    logging.basicConfig(
        format="{asctime} [{module}:{levelname}] {message}",
        style="{",
        level=max(logging.DEBUG, logging.WARNING - (snakemake.params.verbosity * 10)),
    )

    gffutils_cache = Path(snakemake.params['gffutils_cache']) if ('gffutils_cache' in snakemake.params) else None
    if gffutils_cache is not None:
        gffutils_cache.mkdir(exist_ok = True)

    main(
        reference_path = Path(snakemake.input['assembly']),
        annotation_path = Path(snakemake.input['annotation']),
        output_path = Path(snakemake.output['transcriptome']),
        gffutils_cache_path = gffutils_cache
    )


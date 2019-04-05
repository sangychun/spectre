#!/usr/bin/env python

"""SPECtre

Usage: 
    spectre.py [-hc FLOAT] [-t INT] [-m FILE] (-g FILE | -k FILE) -f FILE -b FILE

Options:
  -h, --help          Show this help screen
  -c, --cutoff=FLOAT  Minimum TPM cutoff for active translation [default: 3.0]
  -t, --threads=INT   Number of threads for multiprocessing [default: 1]
  -m, --map=FILE      Path to the UCSC-Ensembl map file [default: None]
  -g, --gtf=FILE      Path to the GTF/GFF annotation file
  -k, --known=FILE    Path to the UCSC knownGene.txt file
  -f, --fasta=FILE    Path to genome FASTA file
  -b, --bam=FILE      Path to BAM file of aligned reads

"""


# Import standard modules:
#import sys
#import re
#import math
#import pickle
#import operator
#import itertools
#import collections

#from functools import partial
#from multiprocessing import Pool
#from multiprocessing.dummy import Pool as ThreadPool

# Import third-party modules:
from docopt import docopt
#import HTSeq as hts
#import pandas as pd
#import numpy as np

# Import SPECtre modules:
#from utils import *
#from anno import *
#from metagene import *
#from scoring import *

def check_inputs():
    

def main():
    if check_inputs():

    """
    # Read in arguments:
    #arguments = docopt(__doc__)
    # Parse parameters:
    #annotation_file = arguments['<annotations>']
    #annotation_type = arguments['<type>']
    #bam_file = arguments['<bam>']
    # Parse options:
    #out_path = arguments['-o']
    #map_file = arguments['-m']
    #tpm_min = arguments['-c']
    #nt = arguments['n']
    # Validate the chromosome formats in the supplied input files:
    validated = check_input_chromosomes(bamfile=bam_file, annofile=annotation_file, annotype=annotation_type)
    if not validated:
        # Validation failure is a fatal error:
        print('Input validation failed, please check that the choromosomes in your BAM file match your reference annotations.')
        sys.exit()
    else:   
        # Count the number of mapped reads:
        n_mapped_reads = count_mapped_reads(infile=bamfile)
        # Parse the read alignment offsets:
        read_offset_positions = parse_custom_offsets(offsets=offsets_file) if offsets_file else None
        # Load BAM file into memory as an HTSeq object:
        alignments = HTSeq.BAM_Reader(bam_file)
        # Adjust the alignments based on the read offset positions:
        offset_alignments = offset_read_alignment_positions(bam=alignments, offsets=read_offset_positions)
        # Build the transcript annotation database:
        anno_db = (load_ensembl_annotations(infile=annotation_file) if annotation_type == 'GTF'
            else load_ucsc_annotations(infile=annotation_file, mapfile=map_file) if (map_file is not None
            and annotation_type == 'UCSC') else load_ucsc_annotations(infile=annotation_file, mapfile=None)
            if (map_file is None and annotation_type == 'UCSC') else None)
        # Extract the coverage over each transcript/region, then normalize coverage to number of mapped reads:
        coverage = calculate_coverage_over_transcript_regions(db=anno_db, bam=bam_file, nreads=n_mapped_reads, threads=nt)
        # Calculate the SPECtre score over each transcript/region:
        coverage = calculate_coherence_over_regions(db=coverage)
        # Calculate the transcripts per million mapped reads (TPM) over each transcript/region:
        coverage = calculate_transcripts_per_million(db=coverage)
        # Build the active/inactive score distributions:
        model = build_translational_probability_model(db=coverage, tpm_minimum=cutoff)
        threshold = binary_search_translational_threshold(df=model)
        # Calculate the posterior probability of translation for each transcript/region:
        coverage = calculate_posterior_probability_by_region(db=coverage, model=model, cutoff=threshold)
        # Pickle results and databases to output directory:
        anno_db.to_pickle(outpath + '/spectre_annotations.pkl')
        coverage.to_picke(outpath + '/spectre_coverage.pkl')
        model.to_pickle(outpath + '/spectre_model.pkl')
        # Pickle basic metrics and temporary variables:
        out_file = open(out_path + '/spectre_metrics.pkl', 'wb')
        pickle.dump({'n_mapped_reads': n_mapped_reads, 'offsets': read_offset_positions, 
            'threshold': threshold}, out_file)
        out_file.close()
    """
    pass

if __name__ == "__main__":
    # Extract input parameters from docstring:
    args = docopt(__doc__)
    # Re-format the annotation file variable:
    afile = args['--gtf'] if args['--gtf'] is not None else args['--known']
    atype = 'Ensembl' if args['--gtf'] is not None else 'UCSC'
    main()

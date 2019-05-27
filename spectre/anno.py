#!/usr/bin/env python

""" Annotation file parser.

This provides spectre a common data structure for UCSC and Ensembl/GENCODE alignments and workflows.

INPUT FILES FORMAT: Transcript annotation formats accepted for use as part of the spectre analytical pipeline are 
limited to the following.

1) Ensembl GFF/GTF (https://ensembl.org/info/website/upload/gff.html):

    Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" 
    columns should be denoted with a '.'

    1. seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
    2. source - name of the program that generated this feature, or the data source (database or project name)
    3. feature - feature type name, e.g. Gene, Variation, Similarity
    4. start - start position of the feature, with sequence numbering starting at 1.
    5. end - end position of the feature, with sequence numbering starting at 1.
    6. score - a floating point value.
    7. strand - defined as + (forward) or - (reverse).
    8. frame - one of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, 
               '1' that the second base is the first base of a codon, and so on.
    9. attribute - a semicolon-separated list of tag-value pairs, providing additional information about each feature.

    NOTE: The 'gene_type' or 'gene_biotype' attribute must be provided for spectre to build its distribution of coding
    and non-coding regions, failure to provide this attribute will result in program termination.

2) a UCSC knownGene table with geneSymbol derived from the kgXref database (refer to the `Addendum` section below for
instructions on how to generate this additional field):

    Fields must be tab-delimited.

    1. name - name of gene
    2. chrom - reference sequence chromosome or scaffold
    3. strand - + or - for strand
    4. txStart - transcript start position (or end position for minus strand item)
    5. txEnd - transcript end position (or start position for minus strand item)
    6. cdsStart - coding region start (or end position if minus strand item)
    7. cdsEnd - coding region end (or start position if minus strand item)
    8. exonCount - number of exons
    9. exonStarts - exon start positions (or end positions for minus strand item)
    10. exonEnds - exon end positions (or start positions for minus strand item)
    11. proteinID - UniProt display ID, UniProt accession, or RefSeq protein ID
    12. alignID - unique identifier (for example, GENCODE transcript ID)
    13. geneSymbol - see ADDENDUM below for instructions on how to include this field

    NOTE: The 'gene_type' will be inferred from a combination of protein identity, and the presence of annotated
    5' and 3'UTRs (to be extracted based on non-matching txStart/txEnd with cdsStart/cdsEnd).

ADDENDUM: The default output of the knownGene database does not include the `geneSymbol` field as required by spectre.
This additional field may be added through the Table Browser interface provided by the UCSC Genome Browser team. Steps
detailed below are specific to the hg38 build of the human genome using the GENCODE v29 database:

    1. First, direct your browser to the UCSC Table Browser (http://genome.ucsc.edu/cgi-bin/hgTables)
    2. Ensure that the following options are selected:
        clade: Mammal
        genome: Human
        assembly: Dec. 2013 (GRCh38/hg38)
        group: Genes and Gene Predictions
        track: GENCODE v29
        table: knownGene
        region: genome
        output format: selected fields from primary and related tables
        output file: <file-name>
    3. Click `get output`
    4. Under `Select Fields from hg38.knownGene`, ensure that all boxes are checked
    5. Under `hg38.kgXref fields`:, select `geneSymbol`
    6. Click `get output` (found under the `Select Fields from hg38.knownGene`)


OUTPUT: Gene and transcript-level information will be parsed into a Pandas DataFrame object, which is then pickled for
downstream use by other portions of the spectre package. The information parsed and tabulated will be structured
as follows:

    0. id - unique identification number, primary key
    1. gene_name - name of gene (eg. TP53, ACTB, etc.)
    2. gene_id - gene identifier (eg. ENSG ID, RefSeq ID)
    3. transcript_id - transcript identifier (eg. ENST ID, UCSC ID), primary key
    4. protein_id - protein name (if available)
    5. gene_type - gene type annotation
    6. transcript_type - transcript type annotation
    7. source - data source or name of program that generated this item
    8. chrom - reference sequence chromosome or scaffold
    9. strand - one of '+' for forward, or '-' for reverse
    10. start - equivalent to UCSC txStart
    11. end - equivalent to UCSC txEnd
    12. cds_start - equivalent to UCSC cdsStart
    13. cds_end - equivalent to UCSC cdsEnd
    14. exon_starts - equivalent to UCSC exonStarts
    15. exon_ends - equivalent to UCSC exonEnds
    16. utr5_starts - 5'UTR start positions (or end positions for minus strand item)
    17. utr5_ends - 5'UTR end positions (or start positions for minus strand item)
    18. cds_starts - coding region start positions (or end positions for minus strand item)
    19. cds_ends - coding region end position (or start positions for minus strand item)
    20. utr3_starts - 3'UTR start positions (or end positions for minus strand item)
    21. utr3_ends - 3'UTR end positions (or start positions for minus strand item)

"""

# Import SPECtre utilities:
from utils import *

# Turn off warnings for chained assignments:
pd.options.mode.chained_assignment = None

##################################################################
# Global variables for transcript annotation database generation #
##################################################################
# Transcript annotation fields:
annotation_fields = ['id','gene_name','gene_id','transcript_id','protein_id','gene_type','transcript_type',
    'source','chrom','strand','start','end','cds_start','cds_end','exon_starts','exon_ends','utr5_starts',
    'utr5_ends','cds_starts','cds_ends','utr3_starts','utr3_ends']

# Transcript annotation data types:
annotation_datatypes = [np.dtype('uint8'),np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_'),
    np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('uint8'),
    np.dtype('U1'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),
    np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32'),
    np.dtype('uint32'),np.dtype('uint32'),np.dtype('uint32')]

# Transcript annotation fields for knownGene-formatted files:
known_gene_fields = ['name','chrom','strand','tx_start','tx_end','cds_start','cds_end','exon_count','exon_starts',
    'exon_ends','protein_id','align_id','gene_name']

# Accepted chromosome names. comment in the chromosome set desired or comment out and create a custom set:
# Homo sapiens:
chroms = ['chr' + str(num) for num in list(range(1,23)) + ['M','X','Y']] + list(range(1,23)) + ['MT','X','Y']

def _initialize_annotation_dataframe():
    """ Instantiates an empty DataFrame for parsed transcript annotation records.

    Creates a Pandas DataFrame object with specified columns and data types. Please refer
    to: https://bit.ly/2nloFkG

    Returns:
        df (:pandas:object:DataFrame): Pre-formatted empty Pandas DataFrame to store annotation records

    """
    try:
        df = pd.DataFrame(index=None)
        for col, datatype in zip(annotation_fields, annotation_datatypes):
            df[col] = pd.Series(dtype=datatype)
        if df is None:
            raise ValueError('Initialization of annotation database failed')
    except ValueError:
        return None
    return df

def _extract_coordinates_from_record(rec=None, region=None):
    try:
        exon_chain = expand_exons_to_chain(list(zip(convert_coordinates(rec['exon_starts']), 
            convert_coordinates(rec['exon_ends']))))
        if exon_chain:
            # Extract the CDS start and end coordinates of the transcript:
            cds_start = int(rec['cds_start'])
            cds_end = int(rec['cds_end'])
            # Parse the coordinates for the targeted region:
            if region == 'utr5_starts':
                coordinates = None if cds_start == cds_end else ','.join([
                    str(start) for start, end in collapse_chain_to_ranges(
                        pos for pos in exon_chain if pos < cds_start)])
            elif region == 'cds_starts':
                coordinates = rec['exon_starts'][:-1] if cds_start == cds_end else ','.join([
                    str(start) for start, end in collapse_chain_to_ranges(pos for pos in exon_chain if (
                        pos >= cds_start and pos <= cds_end))])
            elif region == 'utr3_starts':
                coordinates = None if cds_start == cds_end else ','.join([
                    str(start) for start, end in collapse_chain_to_ranges(
                        pos for pos in exon_chain if pos > cds_end)])
            elif region == 'utr5_ends':
                coordinates = None if cds_start == cds_end else ','.join([
                    str(end) for start, end in collapse_chain_to_ranges(
                        pos for pos in exon_chain if pos < cds_start)])
            elif region == 'cds_ends':
                coordinates = rec['exon_ends'][:-1] if cds_start == cds_end else ','.join([
                    str(end) for start, end in collapse_chain_to_ranges(pos for pos in exon_chain if (
                        pos >= cds_start and pos <= cds_end))])
            elif region == 'utr3_ends':
                coordinates = None if cds_start == cds_end else ','.join([
                    str(end) for start, end in collapse_chain_to_ranges(
                        pos for pos in exon_chain if pos > cds_end)])
            else:
                raise TypeError('Invalid target region')
        else:
            raise ValueError('Expansion of exon chain failed')
    except (TypeError, ValueError):
        return None
    return coordinates

def _add_ensembl_record(record=None, database=None):
    """ Adds an Ensembl-formatted GTF record into the transcript annotation database.

    This function parses an Ensembl-formatted GTF record and extracts relevant information
    into the transcript annotation DataFrame. Records are first scanned for valid attributes
    then parsed according to the type of record encountered. With the exception of coding
    start and end coordinates, only `exon` records are parsed into the transcript anno-
    tation database. Once the entire GTF transcript annotation file has been parsed, then
    CDS start and end coordinate positions are inferred from the annotated coding start and
    end positions.

    """
    def extract_type(rec=None):
        try:
            rtype = (
                'biotype' if ('gene_biotype' in record.attr or 'transcript_biotype' in record.attr)
                else 'type' if ('gene_type' in record.attr or 'transcript_type' in record.attr) else None
                )
            if rtype is None:
                raise KeyError('Type extraction from attributes failed')
        except KeyError:
            return None
        return rtype
    def convert_type(rtype=None):
        try:
            ctype = 'coding' if rtype == 'protein_coding' else 'noncoding'
        except:
            return None
        return ctype
    def parse_record(rec=None, db=None):
        try:
            if rec.type == 'start_codon':
                if rec.iv.strand == '+':
                    db.cds_start[db.transcript_id == rec.attr['transcript_id']] = int(rec.iv.start)
                else:
                    db.cds_end[db.transcript_id == rec.attr['transcript_id']] = int(rec.iv.end)
            elif rec.type == 'stop_codon':
                if rec.iv.strand == '+':
                    db.cds_end[db.transcript_id == rec.attr['transcript_id']] = int(rec.iv.end)
                else:
                    db.cds_start[db.transcript_id == rec.attr['transcript_id']] = int(rec.iv.start)
            elif rec.type == 'exon':
                if not db.exon_starts[db.transcript_id == rec.attr['transcript_id']].values[0]:
                    db.exon_starts[db.transcript_id == rec.attr['transcript_id']] = (str(rec.iv.start))
                    db.exon_ends[db.transcript_id == rec.attr['transcript_id']] = (str(rec.iv.end))
                else:
                    db.exon_starts[db.transcript_id == rec.attr['transcript_id']] += (','+str(rec.iv.start))
                    db.exon_ends[db.transcript_id == rec.attr['transcript_id']] += (','+str(rec.iv.end))
            else:
                # Do not parse other types of records:
                pass
        except KeyError:
            pass
        return db
    try:
        if all([record is not None, database is not None]):
            # Parse records with a valid `transcript_id`:
            if 'transcript_id' in record.attr:
                # Instantiate the transcript annotation record in the database if the encountered
                # record is of type `transcript`, note that the `cds_start` and `cds_end` are set
                # transcript `start` and `end` by default, respectively:
                if not (database.transcript_id.any() == record.attr['transcript_id']) and record.type == 'transcript':
                    database = database.append(dict(zip(annotation_fields, [
                        1, # id
                        record.attr['gene_name'], # gene_name
                        record.attr['gene_id'], # gene_id
                        record.attr['transcript_id'], # transcript_id
                        record.attr['transcript_name'], # protein_id
                        convert_type(record.attr['gene_%s' % (extract_type(record))]), # gene_type
                        record.attr['transcript_%s' % (extract_type(record))], # transcript_type
                        record.attr['transcript_source'], # source
                        record.iv.chrom, # chrom
                        record.iv.strand, # strand
                        record.iv.start, # start
                        record.iv.end, # end
                        record.iv.start, # cds_start
                        record.iv.end, # cds_end
                        str(), # exon_starts
                        str(), # exon_ends
                        str(), # utr5_starts
                        str(), # utr5_ends
                        str(), # cds_starts
                        str(), # cds_ends
                        str(), # utr3_starts
                        str() # utr3_ends
                        ])), ignore_index=True)
                else:
                    database = parse_record(rec=record, db=database)
            else:
                # Skip records without a valid `transcript_id` attribute:
                pass
        else:
            raise ValueError('Missing or invalid database or record input')
    except (KeyError, ValueError):
        return None
    return database

def parse_known_genes(infile=None, chromset=None):
    """ Loads a UCSC-formatted knownGene transcript annotation file into the workspace.

    This function takes as input a UCSC-derived knownGene transcript annotation file and
    parses out the 5'UTR, CDS, and 3'UTR coordinate positions from the provided exon start
    and end coordinates. Non-coding transcripts are defined as those that have identical
    `txStart` and `cdsStart` coordinate positions, or alternatively by the lack of protein
    sequence annotation as denoted by the `proteinID` field.

    Args:
        infile (:py:object:str): File pointer to a knownGene transcript annotation file.

    Returns:
        db (:pandas:object:DataFrame): DataFrame of parsed transcript annotation records.

    """
    def get_transcript_type(rec=None):
        assert rec is not None, 'Missing input annotation record'
        try:
            biotype = 'noncoding' if rec['cds_start'] == rec['cds_end'] else 'coding'
            if not biotype:
                raise ValueError('Parsing of biotype from record failed')
        except ValueError:
            return None
        return biotype
    # Initialize the annotation DataFrame prior to parsing the knownGene annotation file:
    df = _initialize_annotation_dataframe()
    # Parse each knownGene record from the input file:
    with open(infile) as f:
        for i, line in enumerate(f):
            record = dict(zip(known_gene_fields, line.strip().split('\t')))
            df = df.append(dict(zip(annotation_fields, [
                i, # id
                record['gene_name'], # gene_name
                record['name'], # gene_id
                record['align_id'], # transcript_id
                record['protein_id'], # protein_id
                get_transcript_type(rec=record), # gene_type
                get_transcript_type(rec=record), # transcript_type
                'UCSC', # source
                record['chrom'], # chrom
                record['strand'], # strand
                record['tx_start'], # start
                record['tx_end'], # end
                record['cds_start'], # cds_start
                record['cds_end'], # cds_end
                record['exon_starts'][:-1], # exon_starts
                record['exon_ends'][:-1], # exon_ends
                _extract_coordinates_from_record(rec=record, region='utr5_starts'),
                _extract_coordinates_from_record(rec=record, region='utr5_ends'),
                _extract_coordinates_from_record(rec=record, region='cds_starts'),
                _extract_coordinates_from_record(rec=record, region='cds_ends'),
                _extract_coordinates_from_record(rec=record, region='utr3_starts'),
                _extract_coordinates_from_record(rec=record, region='utr3_ends')
                ])), ignore_index=True)
    # Subset the final annotation table to the desired chromosome set:
    if chromset is not None:
        df = df[df.chrom.isin(chromset)]
    # Set the index column:
    df = df.reset_index(drop=True).set_index('id')
    return df

def parse_ensembl_gtf(infile=None, subset=None):
    """ Loads an Ensembl-formatted GTF transcript annotation file into the workspace.

    This function takes as input an Ensembl-derived GTF transcript annotation file and
    parses out the 5'UTR, CDS, and 3'UTR coordinate positions from the provided records.
    Since transcript annotations for each gene are recorded across multiple lines, the
    GTF file must first be parsed in its entirety, and then re-parsed to extract the
    5'UTR, CDS and 3'UTR coordinates.

    Args:
        infile (:py:object:str): File pointer to an Ensembl GTF transcript annotation file.
        subset (:py:object:str): The default set of annotation records to include in the
            final database.

    Returns:
        db (:pandas:object:DataFrame): DataFrame of parsed transcript annotation records.

    """
    # Initialize the annotation DataFrame prior to parsing the transcript annotation file:
    df = _initialize_annotation_dataframe()
    # Load the GFF transcript annotation file into an HTSeq GFF_Reader() object:
    gtf = hts.GFF_Reader(infile)
    for line in gtf:
        df = _add_ensembl_record(record=line, database=df)
    # Check that the first pass through the input GTF transcript annotation file was
    # completed successfully:
    assert len(df) > 0, 'Initial extraction of records from GTF failed'
    # Ensure that the exon start and end coordinates for each transcript are ordered:
    df['exon_starts'] = df['exon_starts'].apply(lambda x: ','.join([str(m) for m in sorted(
        [int(n) for n in x.split(',')])]))
    df['exon_ends'] = df['exon_ends'].apply(lambda x: ','.join([str(m) for m in sorted(
        [int(n) for n in x.split(',')])]))
    return df
    # Subset the transcript annotation table the desired source and chromosome set:
    if subset:
        df = df[df.chrom.isin(subset)]
    # Parse the 5'UTR, CDS and 3'UTR start and end coordinates from the exon coordinates:
    count = 0
    for i, row in df.iterrows():
        df.at[i,'id'] = count
        df.at[i,'utr5_starts'] = _extract_coordinates_from_record(rec=row, region='utr5_starts')
        df.at[i,'utr3_starts'] = _extract_coordinates_from_record(rec=row, region='utr3_starts')
        df.at[i,'cds_starts'] = _extract_coordinates_from_record(rec=row, region='cds_starts')
        df.at[i,'utr5_ends'] = _extract_coordinates_from_record(rec=row, region='utr5_ends')
        df.at[i,'utr3_ends'] = _extract_coordinates_from_record(rec=row, region='utr3_ends')
        df.at[i,'cds_ends'] = _extract_coordinates_from_record(rec=row, region='cds_ends')
        count += 1
    return df.set_index('id')

""" For testing purposes:
ensembl = parse_ensembl_gtf(infile='/home/stonyc/sw/repos/spectre/spectre/data/Homo_sapiens.GRCh38.96.test.gtf')
ucsc = parse_known_genes(infile='/home/stonyc/sw/repos/spectre/spectre/data/GENCODEv29.knownGene.test.txt')
"""

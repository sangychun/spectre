#!/usr/bin/env python

"""Aggregate transcript window parser

This code is part of the spectre distribution and governed by its license. Please see the LICENSE file that should have
been included as part of this package.

This module takes as input a Pandas DataFrame of parsed transcript annotation structures and builds an aggregate tran-
script window database of shared regions across a gene family (ie. transcript isoforms). The purpose of these
aggregated regions is to ensure that the SPECtre scoring distributions are built off of unique sections of
the transcriptome.

For example, where U, E, and I indicate untranslated regions, exons, and introns respectively, the shared aggregate
regions (S) of a family of transcript isoforms may be considered as:

                AUG
POSN    0123456789012345678901234567890123456789012345678901234567890123456789
ISO1    UUUUUUUUEEEEEEEEEEEIIIIIIIEEEEEEEEIIIIIIIIIIEEEEEEEEEEEEEIIIIEEEEEUUUU
ISO2    ....UUUUEEEEEEEEEEEIIIIIIIEEEEEEEEIIIIIIIIIIEEEEEIIIIIIIIIIIIEEEEEUUUU
ISO3    ..UUUUUUEEEEEEEEEEEIIIIIIIEEEEEEEEIIIIEEEEEEEEEEEEEEEEIIIIIIIEEEEEUUUU
META    ....SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS....SSSSS....
SPEC    ........SSSSSSSSSSS.......SSSSSSSS....SSSSSSSSSSSSSSSSSSS....SSSSS....

Where ISO1, ISO2 and ISO3 are annotated isoforms of the same gene, and META defines the shared aggregate regions across
all three isoforms, and SPEC indicate the coding portion of those shared isoform regions. Typically, SPECtre will only
score the shared coding regions of the aggregated metagene to build its scoring distributions.

"""

import multiprocessing as mp
import pandas as pd
import numpy as np
import HTSeq as hts

import itertools
import operator
import datetime

ucsc = pd.read_csv('/home/stonyc/Downloads/knownGene', sep='\t', names=['name','chrom','strand','txStart','txEnd','cdsStart',
    'cdsEnd','exonCount','exonStarts','exonEnds','proteinId','alignId','gene'], header=0)
ucsc = ucsc[ucsc.chrom.isin(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
    'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'])]
ucsc = ucsc.reset_index(drop=True)

test = ucsc[ucsc.chrom == 'chrY']
test = test.reset_index(drop=True)

def check_annotations(gene_types=None):
    if gene_types is None:
        return False
    else:
        return True if (len(list(set(gene_types))) == 1 and 
            list(set(gene_types))[0] == 'coding') else False

def get_minimal_region()


def build_metagene_database(df=None):
    


metagene_fields = ['chrom','start','end','name','status','strand','meta_starts','meta_ends']
metagene_dtypes = [np.dtype('unicode_'),np.dtype('uint32'),np.dtype('uint32'),np.dtype('unicode_'),
    np.dtype('uint8'),np.dtype('unicode_'),np.dtype('unicode_'),np.dtype('unicode_')]

mdf = initialize_annotation_dataframe(columns=metagene_fields, dtypes=metagene_dtypes, index=None)

d = ensembl

starts = list(set(d.cds_start))
for start in starts:
    member_ids = [name for name in list(d.transcript_id[d.cds_start == start])]
    member_types = list(d.gene_type[d.cds_start == start])
    if check_annotations(gene_types=member_types) == True:
        if len(member_ids) > 1:
            print(member_ids)
            chains = list()
            sets = list()
            for member in member_ids:
                print('cds_starts:')
                print(d.cds_starts[d.transcript_id == member].values[0])
                print('cds_ends:')
                print(d.cds_ends[d.transcript_id == member].values[0])
                s = convert_coordinates(d.cds_starts[d.transcript_id == member].values[0])
                e = convert_coordinates(d.cds_ends[d.transcript_id == member].values[0])
                chain = expand_exons_to_chain(exons=list(zip(s,e)))
                sets.append(set(chain))
                chains.append(chain)
            shared = sorted(list(set.union(*sets)))
            blocks = difflib.SequenceMatcher(None, *chains)
            shared_index = blocks.get_matching_blocks()[0].size
            print('union:')
            print(shared)
            print('shared:')
            print(shared[:shared_index])

def _check_for_overlaps(df=None):
    """ Check for overlapping transcripts in the reference genome sequence.

    """

import datetime

def create_genomic_index(df=None):
    """ Fast representation of annotated regions in the genome.

    """
    starts = pd.DataFrame(data={'start': 1}, index=df.cdsStart.tolist())
    ends = pd.DataFrame(data={'end': -1}, index=[i + 1 for i in df.cdsEnd.to_list()])
    exons = pd.merge(starts, ends, how='outer', left_index=True, right_index=True).fillna(0)
    exons['region'] = (exons.pop('end') + exons.pop('start')).cumsum()
    return exons

#ens = create_genomic_index(df=ensembl)

def check_position(pos=None, exons=None):
    query = pd.DataFrame(index=[int(pos)])
    check = pd.merge(exons, query, how='outer', left_index=True, right_index=True).fillna(
        method='ffill').loc[query.index].astype(bool)
    if any(check.region):
        return True
    else:
        return False

def check_overlap(i=None, df=None):
    if (i % 500) == 0:
        print('Multithreaded overlaps checked: [ %s / %s ]' % (i+1, len(df)))
    target = df.at[i,'gene']
    gidx = create_genomic_index(df=df[(df.gene != target) & (df.chrom == list(set(df.chrom[df.gene == target]))[0])])
    mpos = list(df.cdsStart[df.gene == target]) + list(df.cdsEnd[df.gene == target])
    mstart, mend = min(mpos), max(mpos)
    if (check_position(pos=mstart, exons=gidx) & (check_position(pos=mend, exons=gidx))):
        return 1
    else:
        return 0

def check_genes_pool(df=None, nt=None):
    ti = datetime.datetime.now()
    pool = mp.Pool(processes=nt)
    results = [pool.apply_async(check_overlap, args=(i, df,)) for i in list(df.index)]
    df['overlap'] = [res.get() for res in results]
    tf = datetime.datetime.now()
    print('Multiprocess benchmark: %s' % (tf - ti))
    return df

sdf = check_genes_pool(df=test, nt=1)


def check_genes(df=None):
    ti = datetime.datetime.now()
    db = df
    db['overlap'] = 0
    for i, row in df.iterrows():
        target = df.at[i,'gene']
        #print('Creating genomic index for: %s [ %s / %s ]' % (gene, i+1, len(genes)))
        t1 = datetime.datetime.now()
        gidx = create_genomic_index(df=db[(db.gene != target) & (db.chrom == list(set(db.chrom[db.gene == target]))[0])])
        t2 = datetime.datetime.now()
        #print('Check gene (%s) for overlaps to rest of transcriptome..' % (gene))
        mstart = min(list(db.cdsStart[db.gene == target]) + list(db.cdsEnd[db.gene == target]))
        mend = max(list(db.cdsStart[db.gene == target]) + list(db.cdsEnd[db.gene == target]))
        if (check_position(pos=mstart, exons=gidx) == True) or (check_position(
            pos=mend, exons=gidx) == True):
            db.overlap[db.gene == target] = 1
        t3 = datetime.datetime.now()
        if (i % 500) == 0:
            print('Completed check for genomic overlaps with gene: %s [ %s / %s ]' % (target, i+1, len(df)))
        #print('Overlap check completed in: %s' % (t3 - t1))
    tf = datetime.datetime.now()
    print('All %s genes checked for overlap, time for completion: %s' % (len(df), tf - ti))
    print('Average check for each gene = %s' % ((tf-ti)/len(df)))
    return db

odf = check_genes(df=test)
#All 898 genes checked for overlap, time for completion: 0:00:14.208302





def check_denovo_overlaps(df=None):
    """ Check if de novo ORFs overlap an annotated ORF.

    Iterate through the remaining de novo annotated ORFs that have not
    passed previous filters, then check for overlap with accepted ORF
    annotations. If the de novo ORF overlaps with an accepted annotated
    ORF, the de novo ORF should be marked for potential removal from
    the final annotation table.

    Inputs:
    Args:
        df (pandas:object:DataFrame): Merged annotation table of current and
            previous strain

    Returns:
        df (pandas:object:DataFrame): Modified annotation table of current and
            previous strain

    """
    def create_genomic_index(df=None):
        """ Fast representation of annotated regions in the genome.

        This function converts the start and end coordinates of the processed
        GFF records to a representation of genomic region transitions. These
        transitions are then converted into an internal representation for
        fast annotation retrieval query.

        """
        assert df.empty is False, 'Empty input DataFrame provided'
        # Convert the start and end coordinates of the processed GFF records to
        # a representation of genomic region transitions:
        starts = pd.DataFrame(data={'start': 1}, index=df.START.tolist())
        # Increment the to account for 1-based coordinate system used by Prokka:
        ends = pd.DataFrame(data={'end': -1}, index=[i + 1 for i in df.END.tolist()])
        exons = pd.merge(starts, ends, how='outer', left_index=True, 
            right_index=True).fillna(0)
        # Convert the start and end transitions into an interval representation
        # that returns an indicator if a query is contained within the start and
        # end coordinate interval:
        exons['region'] = (exons.pop('end') + exons.pop('start')).cumsum()
        assert exons.empty is False, 'Conversion of exon intervals failed'
        # Return the region intervals:
        return exons
    def check_position(pos=None, exons=None):
        """ Query the interval index if it includes and annotated ORF.

        This function takes an input coordinate position, then performs a quick
        merge with the indexed genome to check if it overalps with an annotated
        ORF.

        """
        assert isinstance(pos, int), 'Invalid `pos` input query'
        assert pos > 0, 'Query `pos` must be a positive integer'
        query = pd.DataFrame(index=[int(pos)])
        check = pd.merge(exons, query, how='outer', left_index=True, right_index=True).fillna(
            method='ffill').loc[query.index].astype(bool)
        if any(check.region):
            return True
        else:
            return False
    t0 = datetime.datetime.now()
    print('[ %s ] Check of overlapping de novo ORFs started...' % (t0), flush=True)
    assert df is not None, 'Resolution of repeat regions failed'
    # Create a genomic coordinate index for annotated ORFs that have passed
    # existing filters, by strand:
    pidx = create_genomic_index(df=df[(df.STATUS == 'PASS') & (df.TYPE.isin(
        ['SNP','INS','DEL']) == False) & (df.STRAND == '+')])
    midx = create_genomic_index(df=df[(df.STATUS == 'PASS') & (df.TYPE.isin(
        ['SNP','INS','DEL']) == False) & (df.STRAND == '-')])
    # Check the start and end coordinates of the de novo ORF for overlap
    # with an acceepted annotated ORF:
    for i, record in df[(df.STATUS == 'FAIL') & (df.CHECK == 'denovo')].iterrows():
        conflict = any([check_position(pos=int(record.START), exons=pidx), check_position(
            pos=int(record.END), exons=pidx)]) if record.STRAND == '+' else any([check_position(
            pos=int(record.START), exons=midx), check_position(pos=int(record.END), exons=midx)])
        # Annotate the de novo ORF:
        if conflict == True:
            df.at[i,'CHECK'] = 'denovo_overlap'
        else:
            df.at[i,'STATUS'] = 'PASS'
    t1 = datetime.datetime.now()
    print('[ %s ] Check of overlapping de novo ORFs finished, module runtime: %s' % (t1, t1-t0), flush=True)
    return df



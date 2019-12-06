"""
######################################################
#                                                    #
#              Pipeline Superenhancer                #
#                                                    #
######################################################

Overview
========

# TODO: add documentation!

Usage
=====


Configuration
-------------


Input files
-----------


Requirements
------------


Pipeline output
===============


Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import cgatcore.experiment as E
from cgatcore import pipeline as P
import PipelineSuperenhancer as SE
import glob
import pandas as pd
from cgat.BamTools import bamtools


# Pipeline configuration
P.get_parameters(
		 ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
		  "../pipeline.yml",
		  "pipeline.yml"],
		 )

PARAMS = P.PARAMS


db = PARAMS['database']['url'].split('./')[1]


def connect():
    '''connect to database.
    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(db)

    if not os.path.exists(PARAMS["annotations_database"]):
        raise ValueError(
                     "can't find database '%s'" %
                     PARAMS["annotations_database"])

    statement = '''ATTACH DATABASE '%s' as annotations''' % \
    (PARAMS["annotations_database"])

    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

# ---------------------------------------------------
# Specific pipeline tasks
    
###############################################################################
##################### get latest gene annotations #############################
###############################################################################
@follows(mkdir("annotations.dir"))
@files(None,"annotations.dir/ensemblGeneset.txt")
def fetchEnsemblGeneset(infile,outfile):
    ''' Get the *latest* gene records using biomart. The aim here is NOT to match
        the great gene set: For that we would only want protein coding genes with
        GO annotations '''

    # statement = '''select count(gene_id) from (select gi.gene_id, gi.gene_name,
    #                       gs.contig, gs.start, gs.end, gs.strand
    #                from gene_info gi
    #                inner join gene_stats gs
    #                on gi.gene_id=gs.gene_id
    #                where gi.gene_biotype="protein_coding")
    #             '''

    # new annotations do not have gene_stats col
    # removed error causing row duplications, now statement retrieves max size transcript per gene_id

    statement = '''select a.gene_id, a.gene_name, b.contig, min(b.start) as start, max(b.end) as end, b.strand
                     from gene_info a 
                     inner join geneset_all_gtf b 
                     on a.gene_id = b.gene_id 
                     where b.gene_biotype = "protein_coding" and b.strand = "+" 
                     group by b.gene_id 
                   union 
                   select a.gene_id, a.gene_name, b.contig, max(b.start) as start, min(b.end) as end, b.strand
                     from gene_info a 
                     inner join geneset_all_gtf b on a.gene_id = b.gene_id 
                     where b.gene_biotype = "protein_coding" and b.strand = "-" group by b.gene_id'''

    anndb = os.path.join(PARAMS["annotations_dir"], "csvdb")

    df = SE.fetch_DataFrame(statement, anndb)
    df.to_csv(outfile, index=False, sep="\t", header=True)

    
@transform(fetchEnsemblGeneset,suffix(".txt"),".load")
def uploadEnsGenes(infile,outfile):
    '''Load the ensembl annotation including placeholder GO ID's'''
    P.load(infile, outfile, options='-i "gene_id" -i "go_id" ')

    
@follows(uploadEnsGenes)
@transform(fetchEnsemblGeneset,
           regex(r"(.*).txt"),
           r"\1.bed")
def getAllGenesBed(infile, outfile):
    '''Make bed file of all Ensembl protein coding genes'''

    tmp_dir = PARAMS["tmpdir"]
    
    statement = f'''tmp=`mktemp -p {tmp_dir}` &&
                   awk 
                     'BEGIN {{OFS="\\t"}} {{print $3,$4,$5,$1,0,$6}}' 
                     {infile} 
                     > $tmp &&
                   tail -n +2 $tmp | 
                     sort -u | 
                     sort -k1,1 -k2,2n 
                     > {outfile} &&
                   rm $tmp'''

    P.run(statement)

    
@follows(getAllGenesBed)
def getGeneLists():
    pass


###############################################################################
##########################  make great promoters  #############################
###############################################################################
@follows(getGeneLists, mkdir("greatBeds.dir"))
@files(uploadEnsGenes, "greatBeds.dir/ens_great.bed")
def greatPromoters(infile, outfile):
    ''' Make great promoters for the genes retrieved from Ensembl'''

    basalup = PARAMS["great_basal_up"]
    basaldown = PARAMS["great_basal_down"]
    maxext = PARAMS["great_max"]
    half = PARAMS["great_half"]
    
    statement = '''select distinct contig, start, end, strand, gene_id from ensemblGeneset'''
    df = SE.fetch_DataFrame(statement, db)
    
    locations = [ [str(r[0]), int(r[1]), int(r[2]),str(r[3]), str(r[4])] for i, r in df.iterrows()]
    
    SE.writeGreat(locations,basalup,basaldown,maxext,outfile,half)


@transform(greatPromoters, suffix(".bed"),".load")
def loadGreatPromoters(infile, outfile):
    '''Load the great promoter regions'''
    P.load(infile, outfile, options='-H "chr,start,end,gene_id" -i "gene_id"')

    
@follows(loadGreatPromoters)
def GreatAnnotation():
    pass

###############################################################################
########################### load interval info  ###############################
###############################################################################
@transform("data.dir/*.bam",
           regex(r"(.*).bam"),
           r"\1.bam.bai")
def indexBAM(infile, outfile):
    '''Index input BAM files'''

    statement = '''samtools index %(infile)s %(outfile)s'''

    P.run(statement)


@follows(indexBAM)   
# @transform("data.dir/*.bed",
#            regex(r"(.*).bed"),
#            r"\1_table.txt")
# def addPeakInfo(infile, outfile):
#     '''Add peak information to table'''

#     df = pd.read_table(infile, sep="\t", header=None)
    
#     o = open(outfile, "w")
#     o.write("\t".join (
#         ["contig","start","end","peak_id","peak_score","width","peak_center"]) + "\n")
    
#     for i, r in df.iterrows():
#         contig, start, end, peak_id, peak_score = r[0:5]

#         width = max(start, end) - min(start, end)
#         peak_center = (start + end)/ 2

#         columns = [str(x) for x in [contig, start, end, peak_id, peak_score, width, peak_center]]

#         o.write("\t".join ( columns ) + "\n")
        
#     o.close()

    
# @transform(addPeakInfo, suffix(".txt"), r".load")
# def loadPeakInfo(infile, outfile):
#     P.load(infile, outfile, options='-i "peak_id"')      



###############################################################################
################ Annotate Peaks with "Regulated" Genes  #######################
###############################################################################
@follows(GreatAnnotation, mkdir("regulated_genes.dir"))
@transform("data.dir/*.bed",
           regex(r"data.dir/(.*).bed"),
           add_inputs("greatBeds.dir/ens_great.bed"),
           r"regulated_genes.dir/\1.GREATassociations.txt")
def regulatedGenes(infiles,outfile):
    '''Annotate input peaks with putative regulated genes'''

    infile, greatPromoters = infiles

    statement = f'''intersectBed 
                      -wa 
                      -wb 
                      -a <(awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$3-$2,$2+$3/2}}' {infile}) 
                      -b {greatPromoters} | 
                    cut 
                      -f1-7,11 
                      > {outfile}'''

    # Filter on nearest peak 2 gene later

    P.run(statement)

    
@transform(regulatedGenes,
           regex(r"(.*).txt"),
           r"\1.load")
def loadRegulatedGenes(infile, outfile):
    '''Load the regulated genes'''

    P.load(infile, outfile, 
           options='-H "contig,pstart,pend,peak_id,peak_score,peak_width,peak_center,gene_id" -i "gene_id"')
 
    
@transform(loadRegulatedGenes,
           suffix(".GREATassociations.load"),
           add_inputs(loadGreatPromoters,uploadEnsGenes),
           ".nearestGene.txt")
def regulatedTables(infiles, outfile):
    '''Make an informative table about peaks and "regulated" genes'''
    
    regulated, great, ensGenes = [ P.toTable(x) for x in infiles ]

    query = f'''select distinct r.contig,
                  r.pstart, r.pend, r.peak_id, r.peak_score,
                  g.gene_id, e.gene_name, e.strand,
                  e.start, e.end
                  from {regulated} as r
                  inner join {great} as g
                     on g.gene_id = r.gene_id 
                  inner join {ensGenes} as e
                     on g.gene_id = e.gene_id'''
    
    dbh = sqlite3.connect(db)
    cc = dbh.cursor()
    sqlresult = cc.execute(query).fetchall()

    tmp = outfile + "_tmp"
    
    o = open(tmp, "w")
    o.write("\t".join ( 
            ["chromosome","peak_start","peak_end","peak_id","peak_score","width","dist2peak",
             "gene_id","gene_name","gene_strand","gene_TSS","gene_start","gene_end"]) + "\n" )

    for r in sqlresult:
        contig, pstart, pend, peak_id, peak_score, gene_id, gene_name, gene_strand = r[0:8]
        gene_start, gene_end = r[8:10]
        
        tss = SE.getTSS(gene_start,gene_end,gene_strand)
        pwidth = max(pstart,pend) - min(pstart,pend)
        ploc = (pstart + pend)/2

        if gene_strand == "+": gstrand = 1
        else: gstrand = 2
        if gstrand==1: tssdist = tss - ploc
        else: tssdist = ploc - tss

        columns = [ str(x) for x in 
                    [  contig, pstart, pend, peak_id, pwidth, tssdist, gene_id, gene_name,
                       gstrand, tss, gene_start, gene_end, peak_score] ]
        o.write("\t".join( columns  ) + "\n")
    o.close()

    # get closest genes 2 peaks, 1 gene per peak
    tmp_dir = PARAMS["tmpdir"]
    
    statement = f'''tmp=`mktemp -p {tmp_dir}` &&
                    tail -n +2 {tmp} | 
                      awk 'BEGIN {{OFS="\\t"}} {{print $4,$5,$6,$7,$8,$9,$10,$11,$12,$1,$2,$3,$13}}' | 
                      sort -k10,10 -k11,11n -k3,3n | 
                      cat | uniq -f9 
                      > $tmp &&
                      mv $tmp {outfile}'''

    P.run(statement)

    
@transform(regulatedTables, suffix(".txt"), ".load")
def loadRegulatedTables(infile,outfile):
    P.load(infile,outfile,
           options='-H"peak_id,width,TSSdist,gene_id,gene_name,gene_strand,TSS,gene_start,gene_end,contig,peak_start,peak_end,peak_score" -i "peakid" ')
           #              1         2        3        4         5        6       7      8          9      10       11        12
           

@follows(loadRegulatedTables)
def peak2gene():
    pass


#################################################################################
###################### Count reads in intervals #################################
#################################################################################
@follows(peak2gene, mkdir("interval_beds"))
@transform("data.dir/*_table.txt",
           regex(r"data.dir/(.*)_table.txt"),
           r"interval_beds/PEAKS_\1.bed")
def movePeaks(infile, outfile):
    '''move peaks to new folder'''

    statement = f'''tail -n +2 {infile} | 
                      cut -f1-6 
                      > {outfile}'''
    
    P.run(statement)
    

@follows(movePeaks)
@transform(loadRegulatedTables,
           regex(r"regulated_genes.dir/(.*).nearestGene.load"),
           r"interval_beds/ENHANCERS_\1.tsv")
def getEnhancers(infile, outfile):
    '''Get peaks outside of TSS +/- 2.5kb'''

    table = os.path.basename(infile)[:-len(".load")].replace(".", "_")
    
    # TSS dist is calculated from the center of the peaks
    # therefore wide peaks in ATAC data may have edges closer to TSS than the cutoff
    # correct TSSdist for peak width:  TSSdist + peak_width/2
    
    query = f'''select contig, peak_start, peak_end, peak_id, peak_score, TSSdist+(width/2) 
                 from {table}''' 

    dbh = sqlite3.connect(db)
    cc = dbh.cursor()
    sqlresult = cc.execute(query).fetchall()
    cc.close()

    # convert to df & write file
    tmp = outfile + "_tmp"
    
    o = open(tmp, "w")
    o.write("\t".join (
        ["contig","start","end","peak_id","peak_score","TSSdist"]) + "\n")

    for r in sqlresult:
        contig, start, end, p_id, peak_score, TSSdist = r[0:6]

        columns = [str(x) for x in [
            contig, start, end, p_id, peak_score, TSSdist] ]

        o.write("\t".join ( columns ) + "\n")
    o.close()

    # use awk to filter out peaks adjacent to TSSs:
    statement = f'''cat <(tail -n +2 {tmp} ) | 
                      awk 'BEGIN {{OFS="\\t"}} {{if ($6 > 2500 || $6 < -2500) print $0}}' - | 
                      sort -k1,1 -k2,2n) |
                      awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$3-$2}}' - 
                      > {outfile} &&
                      rm {tmp} '''
    
    # cols: chr, start, end, p_id, score, width, no. peaks in region
    
    P.run(statement)

    
@follows(getEnhancers)
@transform(loadRegulatedTables,
           regex(r"regulated_genes.dir/(.*).nearestGene.load"),
           r"interval_beds/PROMOTERS_\1.bed")
def getPromoterPeaks(infile, outfile):
    '''Get promoter peaks (TSSs +/- 2.5kb)'''

    table = os.path.basename(infile)[:-len(".load")].replace(".", "_")
    
    # get db data
    query = f'''select contig, peak_start, peak_end, peak_id, peak_score, TSSdist+(width/2) as TSSdist 
                 from  {table} ''' 

    dbh = sqlite3.connect(db)
    cc = dbh.cursor()
    sqlresult = cc.execute(query).fetchall()
    cc.close()

    # convert to df & write file
    tmp = outfile + "_tmp"
    o = open(tmp, "w")
    for r in sqlresult:
        contig, start, end, pid, pscore, TSSdist  = r[0:6]

        columns = [str(x) for x in [contig, start, end, pid, pscore, TSSdist]]

        o.write("\t".join ( columns ) + "\n")
    o.close()

    statement = f'''awk 'BEGIN {{OFS="\\t"}} {{if ($6 < 2500) print $1,$2,$3,$4,$5,$3-$2}}' 
                     <(sed 's/-//g' {tmp} ) 
                     > {outfile}
                     rm {tmp}'''

    P.run(statement)

    
@transform(getPromoterPeaks, suffix(".bed"), r".load")
def loadPromoterPeaks(infile, outfile):
    P.load(infile, outfile, options='-H "contig,start,end,peak_id,peak_score,width" -i "peak_id"' )

    
def filterEnhancersAgainstPromotersGenerator():

    enhancers = glob.glob("interval_beds/ENHANCERS*.tsv")
    promoters = glob.glob("interval_beds/PROMOTERS*.bed")

    if len(enhancers)==0:
        yield []

    if len(promoters)==0:
        yield []

    for pbed in promoters:
        pfile = os.path.basename(pbed).lstrip("PROMOTERS_")[:-len(".bed")]
        for ebed in enhancers:
            outfile = ebed.rstrip(".tsv") + ".bed"
            if pfile in ebed:
                 yield [[ebed, pbed], outfile]

                 
@follows(loadPromoterPeaks)
@files(filterEnhancersAgainstPromotersGenerator)
def filterEnhancersAgainstPromoters(infiles, outfile):
    '''Filter ATAC enhancer peaks against ATAC promoter peaks to ensure there is no overlap.
       Also, filter against genes'''

    enhancers, promoters = infiles
    genes = glob.glob("annotations.dir/ensemblGeneset.bed")
    genes = ' '.join(genes)
    
    statement = '''intersectBed 
                     -v 
                     -a <(sort -k1,1 -k2,2n %(enhancers)s ) 
                     -b <(cat %(promoters)s %(genes)s | sort -k1,1 -k2,2n ) 
                     > %(outfile)s'''

    P.run(statement)

    
@transform(filterEnhancersAgainstPromoters, suffix(".bed"), r".load")
def loadEnhancers(infile, outfile):
    P.load(infile, outfile, options='-H "contig,start,end,peak_id,mean_score,width,peak_no" -i "peak_id"' ) 

    
@follows(loadEnhancers)
@transform(filterEnhancersAgainstPromoters,
           regex(r"interval_beds/ENHANCERS_(.*).bed"),
           r"interval_beds/ENHANCERSmerged_\1.tsv")
def mergeEnhancers(infile, outfile):
    '''Merge enhancer peaks within 12.5kb (default) together using bedtools merge'''

    dist = PARAMS["superenhancer_merge_dist"]
    
    statement = f'''awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,"enhancer_"NR,$5,$3-$2,$4}}'
                      <(mergeBed 
                          -c 4,5 
                          -o count,mean  
                          -d {dist} 
                          -i {infile} ) 
                      > {outfile}'''

    P.run(statement)        

    
@transform(mergeEnhancers,
           regex(r"(.*).tsv"),
           add_inputs("annotations.dir/ensemblGeneset.bed"),
           r"\1.bed")
def filterEnhancersAgainstGenes(infiles, outfile):
    '''subtract any regions spanning genes from enhancer merged peaks'''

    enhancers, genes = infiles
    
    statement = f'''intersectBed 
                      -v 
                      -a <(sort -k1,1 -k2,2n {enhancers}) 
                      -b <(cut -f1-3 {genes} | 
                      sort -k1,1 -k2,2n) 
                      > {outfile}'''

    P.run(statement)

    
@transform(filterEnhancersAgainstGenes, suffix(".bed"), r".load")
def load12kbEnhancers(infile, outfile):
    P.load(infile, outfile, options='-H "contig,start,end,peak_id,mean_score,width,peak_no" -i "peak_id"' ) 


@transform("interval_beds/ENHANCERSmerged_*.bed",
           regex(r"interval_beds/ENHANCERSmerged_(.*).bed"),
           r"interval_beds/REGULATORYFEATURES_\1.bed")
def regulatoryFeatures(infile, outfile):
    '''cat all beds to one file & annotate by feature bedtools 
       multiCov isn't fussy about correctly formatted beds, 
       & just tags on counts as last col'''

    mergedEnh = infile
    Enh = infile.replace("merged", "")
    Prom = infile.replace("ENHANCERSmerged", "PROMOTERS")

    statement = f'''cat
                      <(awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$6,"mergedEnhancerPeaks"}}' {mergedEnh} )
                      <(awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$6,"EnhancerPeak"}}' {Enh} )
                      <(awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$6,"PromoterPeak"}}' {Prom} )
                      > {outfile}'''

    P.run(statement)

    
def generate_scoreIntervalsBAM_jobs():

    # list of bed files & bam files, from which to create jobs
    intervals = glob.glob("interval_beds/REGULATORYFEATURES*.bed")
    bams = glob.glob("data.dir/*.bam")
    
    outDir = "BAM_counts.dir/"

    for interval in intervals:
        ifile = [i.split("/")[-1][:-len(".bed")] for i in [interval] ]
        # iterate over intervals generating infiles & partial filenames

        for bam in bams:
            bfile = [b.split("/")[-1][:-len(".bam")] for b in [bam] ]
            # for each interval, iterate over bams generating infiles & partial filenames
            bedfile = ' '.join(str(x) for x in ifile )
            bamfile = ' '.join(str(x) for x in bfile )
            
            output = outDir + bamfile + ".counts.txt"

            if bamfile in bedfile:
               yield ( [ [interval, bam], output ] )

               
@follows(load12kbEnhancers, mkdir("BAM_counts.dir/"))
@files(generate_scoreIntervalsBAM_jobs)
def scoreIntervalsBAM(infiles, outfile):
    '''Count reads in bed intervals'''

    interval, bam = infiles

    tmp_file = bam.replace(".merge.bam", ".tmp")
    
    if bamtools.isPaired(bam):
        # -p flag specifes only to count paired reads
        options = "-p"
        
    else:
        options = " "

    statement = f'''bedtools multicov 
                       {options} 
                       -q 10 
                       -bams {bam} 
                       -bed <(cut -f1-7 {interval} ) 
                       > {outfile} &&
                     sed -i '1i \contig\\tstart\\tend\\tpeak_id\\tpeak_score\\twidth\\tfeature\\ttotal' 
                       {outfile}'''
        
    P.run(statement)

    
@transform(scoreIntervalsBAM, suffix(".txt"), ".load")
def loadIntervalscoresBAM(infile, outfile):
    P.load(infile, outfile, options='-i "gene_id"')

    
def generator_BAMtotalcounts():
    bams = glob.glob("data.dir/*.bam")
    outDir = "BAM_counts.dir/"

    for bam in bams: 
        bfile = [b.split("/")[-1][:-len(".bam")] for b in [bam] ]
        bfile = ' '.join(str(x) for x in bfile ) # unpack from list
        output = outDir + bfile + ".total_reads.txt"
        
        yield ( [ bam, output ] )

        
@follows(loadIntervalscoresBAM)
@files(generator_BAMtotalcounts)
def BAMtotalcounts(infile, outfile):
    '''Count total reads in BAM for normalisation'''

    if bamtools.isPaired(infile):
        statement = f'''samtools view -f 2 {infile} | wc -l | awk 'BEGIN {{OFS="\\t"}} {{print $0/2}}' > {outfile}''' 
        # count only reads mapped in proper pairs
    else:
        statement = '''samtools view -F 4 {infile} | wc -l  | awk 'BEGIN {{OFS="\\t"}} {{print $0}}' > {outfile}''' 
        # exclude unmapped reads

    P.run(statement)
    

@transform(BAMtotalcounts, suffix(".txt"), r".load")
def loadBAMtotalcounts(infile, outfile):
    P.load(infile, outfile, options='-H "total_reads"')

    
def normaliseBAMcountsGenerator():
    total_reads = glob.glob("BAM_counts.dir/*.total_reads.txt")
    counts = glob.glob("BAM_counts.dir/*.counts.txt")

    if len(total_reads)==0:
        yield []
        
    outdir = "BAM_counts.dir/"

    # generate jobs & match total_reads to counts files
    for interval_count in counts:
        count_table = os.path.basename(interval_count)[:-len(".counts.txt")]

        for read_count in total_reads:
            output = outdir + count_table + ".norm_counts.txt"

            bam_str = os.path.basename(read_count)[:-len(".total_reads.txt")]

            if bam_str in count_table:
                yield ( [ [interval_count, read_count], output ] )

                
@follows(loadBAMtotalcounts)
@files(normaliseBAMcountsGenerator)
def normaliseBAMcounts(infiles, outfile):
    '''normalise BAM counts for file size'''
    
    to_cluster = True
    
    interval_counts, total_reads = infiles
    
    # read counts to dictionary
    name = os.path.basename(total_reads)[:-len(".txt")]
    total_reads = open(total_reads, "r").read().replace("\n", "")
    counts = {}
    counts[name] = total_reads
    counts_str = counts[name]

    # get db table names
    interval_table = os.path.basename(interval_counts)[:-len(".txt")].replace(".", "_")

    # get db data
    query = f'''SELECT contig, start, end, peak_id,
            total, width, feature FROM {interval_table} '''

    dbh = sqlite3.connect(db)
    cc = dbh.cursor()
    sqlresult = cc.execute(query).fetchall()
    cc.close()

    # convert to df & write file
    tmp = outfile + "_tmp"
    
    o = open(tmp, "w")
    o.write("\t".join (
        ["contig","start","end","peak_id","total","width","feature"]) + "\n")

    for r in sqlresult:
        contig, start, end, p_id, total, width, feature = r[0:7]

        
        #pwidth = max(start,end) - min(start,end)

        columns = [str(x) for x in [
            contig, start, end, p_id, total, width, feature ] ]

        o.write("\t".join ( columns ) + "\n")
    o.close()
    
    # normalise counts for sequencing depth using awk & normalise for interval width
    statement = f'''awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5,$5/{counts_str}*1000000,$5/{counts_str}*1000000/$6,$6,$7}}' 
                     <(tail -n +2 {tmp} ) 
                     >{outfile} &&  
                   sed -i '1i \contig\\tstart\\tend\\tpeak_id\\traw_counts\\tRPM\\tRPM_width_norm\\twidth\\tfeature' 
                     {outfile}'''
   
    P.run(statement, job_memory="5G")

    
@transform(normaliseBAMcounts, suffix(".txt"), ".load")
def loadnormaliseBAMcounts(infiles, outfiles):
    P.load(infiles, outfiles)    

    
@follows(loadnormaliseBAMcounts)
def readCounts():
    pass


# @follows(readCounts)
# @transform("BAM_counts.dir/*.norm_counts.load",
#            regex(r"BAM_counts.dir/(.*).norm_counts.load"),
#            r"\1.mergedEnhancerCounts.bed")
# def getSEbeds(infile, outfile):
#     '''Get bed files of 12.5kb enhancer merged peaks with ATAC 
#     read density in 7th column for get-SuperEnhancers.R script'''

#     counts_table = os.path.basename(infile)[:-len(".load")].replace(".", "_")
#     interval_table = "ENHANCERSmerged_" + counts_table.replace("_norm_counts", "")
    
#     # get db data
#     query = '''select a.contig, a.start, a.end, a.peak_id, b.peak_no, a.RPM_width_norm 
#                  from %(counts_table)s a 
#                  inner join %(interval_table)s b 
#                  on a.peak_id= b.peak_id''' % locals()

#     dbh = sqlite3.connect(PARAMS["database"])
#     cc = dbh.cursor()
#     sqlresult = cc.execute(query).fetchall()
#     cc.close()

#     # convert to df & write file
#     tmp = outfile + "_tmp"
    
#     o = open(tmp, "w")

#     for r in sqlresult:
#         contig, start, end, p_id, peak_no, norm_counts = r[0:6]

#         strand = "."

#         columns = [str(x) for x in [
#             contig, start, end, p_id, peak_no, norm_counts ] ]

#         o.write("\t".join ( columns ) + "\n")
#     o.close()

#     statement = '''intersectBed
#                      -wa
#                      -wb
#                      -a <( awk 'BEGIN {OFS="\\t"} {print $1,$2,$3,$4,$5,".",$6}' %(tmp)s )
#                      -b %(greatPromoters)s |
#                    cut -f1-7,11
#                      > %(outfile)s ; checkpoint;
#                      rm %(tmp)s''' # statement incomplete
#     P.run()

    
# @follows(getSEbeds)
# @transform("*.bed",
#            regex(r"(.*).bed"),
#            add_inputs("greatBeds.dir/ens_great.bed"),
#            r"\1_annotated.bed")
# def annotateEnhancers(infiles,outfile):
#     '''Used intersectBed and the great promoters to find regulated Genes'''

#     infile, greatPromoters = infiles

#     statement = '''intersectBed 
#                      -wa 
#                      -wb 
#                      -a <(cut -f1-7 %(infile)s) 
#                      -b %(greatPromoters)s | 
#                      cut -f1-7,11 
#                      > %(outfile)s''' 

#     P.run()


# @transform(annotateEnhancers,
#            regex(r"(.*).bed"),
#            r"\1.load")
# def loadannotateEnhancers(infile, outfile):
#     '''Load the regulated genes'''

#     P.load(infile, outfile, 
#            options='-H "contig,start,end,enhancer_id,no_atac_peaks,strand,ATAC_RPM,gene_id" -i "enhancer_id"')

    
# @transform(loadannotateEnhancers,
#            suffix(r".load"),
#            add_inputs(uploadEnsGenes),
#            r"_table.txt")
# def enhancerTables(infiles, outfile):
#     '''Make an informative table about peaks and "regulated" genes'''
#     regulated, ensGenes = [ P.toTable(x) for x in infiles ]

#     query = '''select distinct r.contig,
#                   r.start, r.end, r.enhancer_id,
#                   r.no_atac_peaks, r.strand, r.ATAC_RPM,
#                   e.gene_name, e.start, e.end, e.strand
#                   from %s as r
#                   inner join %s as e
#                      on r.gene_id = e.gene_id
#                   ''' % (regulated, ensGenes)

#     print(query)
    
#     dbh = sqlite3.connect(PARAMS["database"])
#     cc = dbh.cursor()
#     sqlresult = cc.execute(query).fetchall()

#     tmp = outfile + "_tmp"
    
#     o = open(tmp,"w")
#     o.write("\t".join ( 
#             ["chromosome","start","end","enhancer_id","no_atac_peaks","strand","ATAC_RPM","width","dist2peak",
#              "gene_name","gene_TSS","gene_start","gene_end","gene_strand"]) + "\n" )

#     for r in sqlresult:
#         contig, pstart, pend, enhancer_id, no_atac_peaks, strand, ATAC_RPM, gene_name = r[0:8]
#         gene_start, gene_end, gene_strand = r[8:11]
        
#         if gene_strand == "+": gstrand = 1
#         else: gstrand = 2

#         tss = SE.getTSS(gene_start,gene_end,gene_strand)

#         pwidth = max(pstart,pend) - min(pstart,pend)
#         ploc = (pstart + pend)/2

#         if gstrand==1: tssdist = tss - ploc
#         else: tssdist = ploc - tss

#         columns = [ str(x) for x in 
#                     [  contig, pstart, pend, enhancer_id, no_atac_peaks,
#                        strand, ATAC_RPM, pwidth, tssdist, gene_name] ]
#         o.write("\t".join( columns  ) + "\n")
#     o.close()

#     # get closest genes 2 peaks, 1 gene per peak
#     tmp_file = tmp + "_tmpfile"

    
#     statement = '''tail -n +2 %(tmp)s 
#                 | awk 'BEGIN {OFS="\\t"} {print $4,$5,$6,$7,$8,$9,$10,$1,$2,$3}' 
#                 | sort -k8,8 -k9,9n -k6,6n 
#                 | cat | uniq -f7 > %(tmp_file)s && mv %(tmp_file)s %(outfile)s
#                 && rm %(tmp)s''' % locals()

#     P.run()

# @transform(enhancerTables, suffix(".txt"), ".load")
# def loadenhancerTables(infile,outfile):
#     P.load(infile,outfile,
#            options='-H"enhancer_id,not_atac_peaks,strand,ATAC_RPM,pwidth,tssdist,gene_name,contig,start,end" -i "peak_id" ')




# ---------------------------------------------------
# Generic pipeline tasks
@follows(loadnormaliseBAMcounts)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    

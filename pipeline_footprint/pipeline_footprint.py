"""===========================
Pipeline template
===========================

Overview
========

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
import cgatcore.pipeline as P
import PipelineFootprint as F 
import glob
import pandas as pd

# Pipeline configuration
P.get_parameters(
		 ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
		  "../pipeline.yml",
		  "pipeline.yml"],
		 )

PARAMS = P.PARAMS

db = PARAMS['database']['url'].split('./')[1]
tmp_dir = PARAMS['tmp_dir']

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
@collate("data.dir/*.bam",
       regex(r"(.*)_r[1-9].bam"),
       r"\1_merge.bam")
def mergeBams(infiles, outfile):
    '''merge bams into single file per experimental condition'''

    infiles = ' '.join(infiles)

    statement = f'''samtools cat 
                      {infiles} |
                    samtools sort - 
                      > {outfile}'''
    
    P.run(statement, job_threads=5)

    
@transform(mergeBams,
           suffix(r".bam"),
           r".bam.bai")
def indexBams(infile, outfile):
    '''index merged bams'''

    statement = f'''samtools index -b {infile} {outfile}'''
    
    P.run(statement)


@follows(indexBams, mkdir("cutsites.dir"))
@transform(mergeBams,
           regex(r"data.dir/(.*)_merge.bam"),
           r"cutsites.dir/\1.cutsites.bed.gz")
def getCutSites(infile, outfile):
    '''Seperate forward and reverse strands, 
       offset read start position (F + 4, R - 5),
       and merge back to 1 file'''

    tsv = outfile.replace(".gz", "")
        
    statement = f'''tmp=`mktemp -p {tmp_dir}` &&
                   samtools sort -n {infile} > $tmp &&
                   bedtools bamtobed -i $tmp |
                   awk 'BEGIN {{OFS="\\t"}}
                     {{pshift = $2+4}}
                     {{nshift = $2-5}} 
                     {{if ($6 = "+") print $1,pshift,pshift+1 ;
                     else if ($6 = "-") print $1,nshift,nshift-1 }}' -
                     > {tsv} &&
                   gzip {tsv} &&
                   rm $tmp'''

    P.run(statement, job_memory="6G")

    
@follows(getCutSites)
@transform("data.dir/*.bed",
          regex(r"data.dir/(.*).bed"),
          r"coverage.dir/\1.window.bed")
def offsetPeaks(infile, outfile):
    '''Offset peaks to peak centre +/- n b.p.
       And remove peaks from chr*random, chrun*, chrM'''
    
    summits = PARAMS["peaks_summits"]
    window = PARAMS["peaks_search_range"]
    chrom_sizes = PARAMS["annotations_chrom_sizes"]

    statement = []
    statement.append(f'''tmp=`mktemp -p {tmp_dir}` &&
                        awk 'BEGIN {{OFS="\\t"}} ''' )
    
    if summits == "True":
         statement.append(f'''{{start=$6-1}} {{end=$6+1}}''')
    else:
        statement.append(f'''{{x=($3-$2)/2}} {{centre=$2+x}} {{start=centre-1}} {{end=centre+1}} ''')

    statement.append(f'''{{print $1,int(start),int(end),$4}}'
                          {infile}
                          > $tmp &&
                        slopBed
                          -g {chrom_sizes}
                          -b {window}
                          -i $tmp |
                        grep -v "chr.*random" - | 
                          grep -v "chrUn.*" - |
                          grep -v "chrM" - 
                          > {outfile} &&
                        rm $tmp''')
    
    statement = ' '.join(statement)

    P.run(statement)


@follows(offsetPeaks)
@files(None, "readCounts.tsv")
def getReadCounts(infile, outfile):
    '''get total read counts in merged BAMs for normalisation'''

    statement = f'''tmp=`mktemp -p {tmp_dir}` &&
                    for b in data.dir/*merge.bam ; 
                      do name=`echo $b | 
                        sed 's/data.dir\///' | 
                        sed 's/_merge.bam//'` ; 
                      samtools idxstats $b | 
                        awk -v name=$name 'BEGIN {{OFS="\\t"}}
                        {{sum += $3+$4}} END {{print name,sum}}' - 
                        >> $tmp ;
                      done &&
                    mv $tmp {outfile}'''
    
    P.run(statement, job_memory="5G")

    
@transform(getReadCounts,
           suffix(r".tsv"),
           r".load")
def loadReadCounts(infile, outfile):
    '''load total read counts'''
    
    P.load(infile, outfile, options='-H "sample,total_reads" ')


def coverageBedGenerator():
    '''Create jobs to annotate infiles with motif matches'''

    infiles = glob.glob("coverage.dir/*.window.bed")
    cutsites = glob.glob("cutsites.dir/*.cutsites.bed.gz")

    outdir = "coverage.dir/"

    for c in cutsites:
        for i in infiles:
            c_name = os.path.basename(c).replace(".cutsites.bed.gz", "")
            i_name = os.path.basename(i).replace(".window.bed", "")

            outfile = outdir + '.'.join([i_name, c_name]) + ".bed.gz"

            yield [[i, c], outfile]
            #print([[i, c], outfile])

            
@follows(loadReadCounts, mkdir("coverage.dir"))
@files(coverageBedGenerator)
def coverageBed(infiles, outfile):
    '''Calculate per base coverage of peaks by motif files.
       And correct position to be relative to peak centre, rather than start'''

    bed, cutsite = infiles

#    job_memory = "200G" # jobs keep failing with less memory

    statement = f'''tmp1=`mktemp -p {tmp_dir}` &&
                    tmp2=`mktemp -p {tmp_dir}` &&
                    zcat {cutsite} > $tmp1 &&
                    coverageBed 
                      -d 
                      -a {bed}                      
                      -b $tmp1 
                      > $tmp2 &&
                    awk 'BEGIN {{OFS="\\t"}}
                      {{pwidth=$3-$2}}
                      {{pcentre=pwidth/2}}
                      {{if ($5 >= pcentre) print $1,$2,$3,$4,$5-pcentre,$6;
                      else print $1,$2,$3,$4,"-"pcentre-$5,$6}}' $tmp2 |
                    gzip -c - > {outfile} &&
                    rm $tmp1 $tmp2'''

    P.run(statement, job_memory="200G")

    
@transform("coverage.dir/*.bed.gz",
           regex(r"(.*).bed.gz"),
           r"\1.bed.norm.gz")
def footprintNorm(infile, outfile):
    '''Normalise ATAC cutting frequency'''

    out = outfile.replace(".gz", "")
    script = PARAMS["pipeline_dir"] + "python/footprintNorm.py"
    
    statement = f'''python {script} 
                      --infile {infile} 
                      --outfile {out}
                      --binsize 1 
                      --database {db} &&
                    gzip {out}'''

    P.run(statement, job_memory="200G")
    
    
@transform(footprintNorm,
           regex(r"(.*).bed.norm.gz"),
           r"\1.png")
def plotFootprint(infile, outfile):
    '''Normalise coverage profiles & plot cut site frequency'''

    region_label = os.path.basename(infile).split(".")
    region = region_label[0]
    label = region_label[1]

    window = PARAMS["plot_window"]
    distance = int(window)/2
    smoothed = PARAMS["plot_show_unsmoothed"]
    bandwidth = PARAMS["plot_bandwidth"]
    script = PARAMS["pipeline_dir"] + "R/plotFootprint.R"
    
    statement = f'''Rscript {script}
                      --infiles {infile}
                      --outfile {outfile}
                      --labels {label}
                      --title {region} 
                      -b {bandwidth}
                      --xlims {distance}
                      --unsmoothed {smoothed}'''

    P.run(statement)


def plotFootprintsGenerator():

    footprints = glob.glob("coverage.dir/*.bed.norm.gz") # get all footprint files
    groups = list(set([os.path.basename(x).split(".")[0] for x in footprints])) # get unique footprint regions for plots

    for g in groups:
        infiles = [x for x in footprints if g in os.path.basename(x).split(".")[0]]
        labels = [os.path.basename(x).split(".")[1].replace(".bed.norm.gz", "") for x in infiles]
        outfile = g + ".all_footprints.png"
        
        yield [infiles, outfile, labels]


@follows(coverageBed, plotFootprint)
@files(plotFootprintsGenerator)
def plotFootprints(infiles, outfile, labels):
    '''Plot footprints for all samples over each region'''

    region = ''.join(list(set([os.path.basename(x).split(".")[0] for x in infiles])))
    infiles = ','.join(infiles)
    labels = ','.join(labels)
    
    window = PARAMS["plot_window"]
    distance = int(window)/2
    smoothed = "FALSE" # will be too noisy 
    bandwidth = PARAMS["plot_bandwidth"]
    script = PARAMS["pipeline_dir"] + "R/plotFootprint.R"
    
    statement = '''Rscript {script}
                     --infiles {infiles}
                     --outfile {outfile}
                     --labels {labels}
                     --title {region} 
                     -b {bandwidth}
                     --xlims {distance}
                     --unsmoothed {smoothed}'''

    P.run(statement)

    
    
# ---------------------------------------------------
# Generic pipeline tasks
@follows(plotFootprint, plotFootprints)
def full():
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))    

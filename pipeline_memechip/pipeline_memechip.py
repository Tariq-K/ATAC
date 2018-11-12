##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline template
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_macs2.py config

Input files
-----------
- Bam files & bai indexes

- sample files should be of format: <exp>-<cond>-<treatment>-<sample>.genome.bam
- with matching file for input:     <exp>-<cond>-<treatment>-<WCE>.genome.bam
- if only 1 input is needed for all samples it can be specified in pipeline.ini

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* samtools >= 1.1

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""
import sys, tempfile 
import glob, os 
import subprocess
import sys
import rpy2.robjects as R
import CGAT.Experiment as E
import logging as L
from ruffus import *
import sqlite3
import pybedtools
import PipelineMemechip as S
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineGO as PGO
import CGAT.Biomart as Biomart
import pandas as pd
from CGAT.IndexedFasta import IndexedFasta
import CGAT.Bed as BED
import CGAT.IOTools as IO
import CGAT.FastaIterator as FastaIterator
import CGAT.Database as Database
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns


# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# -----------------------------------------------
# Utility functions

# def connect():
#     '''utility function to connect to database. Use this method to connect to the pipeline database.
#     Additional databases can be attached here as well. Returns an sqlite3 database handle. '''

#     dbh = sqlite3.connect(PARAMS["database"])
#     statement = '''ATTACH DATABASE '%s' as annotations'''% (
#         PARAMS["annotations_database"])
#     cc = dbh.cursor()
#     cc.execute(statement)
#     cc.close()

#     return dbh

def connect():
    '''
    Setup a connection to an sqlite database
    '''

    dbh = sqlite3.connect(PARAMS['database'])

    return dbh


################################################
#####              Meme-ChIP               #####
################################################
#@follows(connect)
@transform("data.dir/*.bed",
           regex(r"(.*).bed"),
           r"\1_meme.bed")
def peakSummit(infile, outfile):
    '''Make bed of peak_summit (+/- 1bp), sorted by peak qvalue'''

    infile = ''.join([x for x in [infile] if "_meme" not in x]) # prevent re-running of this job

    if os.path.isfile(infile):
        df = pd.read_csv(infile, sep="\t", header=None)

        df = df.iloc[:, 0:6] # subset on useful columns
        df.columns = ["contig", "start", "end", "peak_id", "score", "summit"]

        df = df.sort_values("score", ascending=False) # sort by score

        df["summit_start"] = df["summit"] - 1
        df["summit_end"] = df["summit"] + 1

        df = df[["contig", "summit_start", "summit_end", "peak_id", "score"]]

        df.to_csv(outfile, sep="\t", header=None, index=None)

            
def getMemeForegroundBedGenerator():
    
    beds = glob.glob("data.dir/*_meme.bed")

    if len(beds) == 0:
        yield []

    npeaks = [x.strip() for x in str(PARAMS["memechip_npeaks"]).split(",")]
    widths = [y.strip() for y in str(PARAMS["memechip_widths"]).split(",")]

    outdir = "meme.seq.dir/"
    
    for bed in beds:
        int_name = os.path.basename(bed)[:-len("_meme.bed")]
        for n in npeaks:
            for w in widths:
                fgname = ".".join([int_name,n,w,"foreground","bed"])
                outfile = outdir + fgname
                
                yield [ bed, outfile ]

@follows(peakSummit, mkdir("meme.seq.dir/"))
@files(getMemeForegroundBedGenerator)
def getMemeForegroundBed(infile, outfile):
    '''Make meme beds of peak center +/- n bp, for n top scoring peaks'''

    chrom_sizes = PARAMS["annotations_chrom_sizes"]
    
    npeaks = os.path.basename(outfile)[:-len(".foreground.bed")].split(".")[-2]
    width = os.path.basename(outfile)[:-len(".foreground.bed")].split(".")[-1]
    start_width = (int(width) / 2) + 1 # correct for previous offset
    end_width = (int(width) / 2) - 1 # as bed file is chr, peak_centre_minus1, peak_centre_plus1
    
    if npeaks != "all":
        statement = '''slopBed -i <(  head -n %(npeaks)s %(infile)s )
                    -g %(chrom_sizes)s -l %(start_width)s -r %(end_width)s -s > %(outfile)s''' % locals()
        
    else:
        statement = '''slopBed -i  %(infile)s 
                    -g %(chrom_sizes)s -l %(start_width)s -r %(end_width)s -s > %(outfile)s''' % locals()
   
    P.run()

@transform(getMemeForegroundBed, suffix(".bed"), r".load")
def loadMemeForegroundBed(infile, outfile):
    P.load(infile, outfile, options='-H "contig,start,end,peak_id,TFscore" ')
    
@follows(loadMemeForegroundBed)
@transform(getMemeForegroundBed,
           regex(r"(.*).foreground.bed"),
           r"\1.background.bed")
def getMemeBackgroundBed(infile, outfile):
    '''get bed file of peak flanking regions (of equal width to peak) for meme background model'''

    genome_idx = os.path.join(PARAMS["annotations_mm10dir"],"assembly.dir/contigs.tsv")
    
    statement ='''sort -k1,1 -k2,2n %(infile)s 
               | bedtools flank -pct -l 1 -r 1 -g %(genome_idx)s > %(outfile)s'''
               # -pct -l 1 -r 1 will create flanking regions = peak width
    
    P.run()

@follows(getMemeBackgroundBed)
@transform("meme.seq.dir/*.bed",
           regex(r"(.*).bed"),
           r"\1.fasta")
def getMemeSequences(infile, outfile):
    '''Get the peak sequences. The genome sequences are already repeat soft-masked'''
  
    genome_fasta = os.path.join(PARAMS["genome_dir"],PARAMS["genome"]+".fasta")
    genome_idx = os.path.join(PARAMS["annotations_mm10dir"],"assembly.dir/contigs.tsv")

    statement = '''fastaFromBed -fi %(genome_fasta)s -bed %(infile)s -fo %(outfile)s'''
    
    P.run()

@follows(getMemeSequences)
@transform("meme.seq.dir/*background.fasta", suffix(".fasta"), ".bfile")
def getMemeBfiles(infile, outfile):
    '''prepare the meme background model'''

    statement='''fasta-get-markov -m 2 %(infile)s  > %(outfile)s''' % locals()

    P.run()
               
@follows(getMemeSequences, getMemeBfiles, mkdir("meme.chip.dir"))
@transform(getMemeSequences,
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           r"meme.chip.dir/\1.memechip")
def runMemeChIP(infile, outfile):
    '''run MemeChIP'''

    outdir = outfile.replace(".memechip", "")
    bfile= infile.replace(".foreground.fasta", ".background.bfile")

    motifDb =  " -db ".join(P.asList(PARAMS["memechip_motif_db"])) # Meme-Chip needs eac db in list to have "-db" flag

    nmotifs = PARAMS["memechip_nmotif"]
    meme_max_jobs = PARAMS["memechip_meme_maxsize"]

    # Meme options:
    
    # nmeme - The upper bound on the number of sequences that are passed to MEME.
            # This is required because MEME takes too long to run for very large sequence sets.
            # All input sequences are passed to MEME if there are not more than limit.
            # default nmeme = 600
            # specified in pipeline.ini
            
    nmeme_list = PARAMS["memechip_npeaks"].split(",")
    print(nmeme_list)
    
    if "all" in nmeme_list: # get max file length if assessing all peaks
        meme_beds = glob.glob("data.dir/*_meme.bed")

        # if len(meme_beds)==0:
        #     print(len(meme_beds))
        #     yield []

        bed = max([len(pd.read_csv(x, sep="\t")) for x in meme_beds])
        param = max([int(x) for x in nmeme_list if "all" not in x])
        nmeme = max([bed, param])
        
    else: # or max length from npeaks list
        nmeme = max([int(x) for x in nmeme_list])
        # no of sequences actually determined by infile, this just prevents truncation
        
    # ccut - The maximum length of a sequence to use before it is trimmed to a central region of this size.
            # A value of 0 indicates that sequences should not be trimmed.

    # meme-maxsize - Change the largest allowed dataset to be size.
            # default meme-maxsize is 100,000.
            # Fine with the default settings for -nmeme (600) and -ccut (100), largest possible dataset size would be 60000.
            # 600,000, equivalent to max of 600 1000bp seq

    # mode - specified in pipeline.ini
    # oops - One Occurence Per Seqeunce - meme assumes that each sequence only contains 1 occurence of each motif
    #      - fastest and most sensitive, but motifs returned may be "blurry" if any sequences lack them
    # zoops - Zero or One Occurence Per Sequence - meme assumes that each sequence contains at most one occurence of each motif
    #       - takes twice as long as oops, useful when motif may not be present in all sequences, but less sensitive to weak motifs present in all sequences
    # anr - Any Numer of Repetitions - meme assumes that each sequence contains any number of non-overlapping occurences of each motif
    #     - Useful if motifs repeat multiple times within a sequence. If this is the case then will be much more sensitve than oops & zoops.
    #     - Takes 10x more computing power than option 1 and is less sensitive to weak, non-repeated motifs

    mode = PARAMS["memechip_mode"]

    job_memory = "5G"
    
    statement = '''meme-chip
                   -oc %(outdir)s
                   -db %(motifDb)s
                   -bfile %(bfile)s
                   -ccut 0
                   -nmeme %(nmeme)s
                   -meme-mod %(mode)s
                   -meme-minw 5
                   -meme-maxw 30
                   -meme-nmotifs %(nmotifs)s
                   -meme-maxsize %(meme_max_jobs)s
                   %(infile)s
                   > %(outfile)s
                ''' 
    
    P.run()

    
def loadMemeTomTomGenerator():

    meme_tomtom = glob.glob("meme.chip.dir/*memechip")

    if len(meme_tomtom) == 0:
        yield []
        
    for tomtom in meme_tomtom:
        basename = os.path.basename(tomtom)[:-len(".memechip")]
        infile = tomtom.replace(".memechip", "/meme_tomtom_out/tomtom.txt")
        outfile = infile[:-len("tomtom.txt")] + basename + "_Meme_tomtom.load"

        if os.path.isfile(infile): 
            # some empty files are generated w/ just headers, causing errors loading empty tables
            # count lines in file, if file more than 1 line, yield job
            with open(infile, "r") as o:
                n = 0
                for line in o:
                    n = n + 1
                if n > 1:
                    yield [ infile, outfile ]
                    
@follows(runMemeChIP)
@files(loadMemeTomTomGenerator)
def loadMemeTomTom(infile, outfile):
    '''load meme tomtom results'''
    P.load(infile, outfile, options='-H "query_id,target_id,optimal_offset,p_value,e_value,q_value,overlap,query_consensus,targe_consensus,orientation" ')

def loadDremeTomTomGenerator():

    meme_tomtom = glob.glob("meme.chip.dir/*memechip")
    
    if len(meme_tomtom) == 0:
        yield []
            
    for tomtom in meme_tomtom:
        basename = os.path.basename(tomtom)[:-len(".memechip")]
        infile = tomtom.replace(".memechip", "/dreme_tomtom_out/tomtom.txt")
        outfile = infile[:-len("tomtom.txt")] + basename + "_Dreme_tomtom.load"

        if os.path.isfile(infile):
            with open(infile, "r") as o:
                n = 0
                for line in o:
                    n = n + 1
                if n > 1:
                    yield [ infile, outfile ]

@follows(loadMemeTomTom)
@files(loadDremeTomTomGenerator)
def loadDremeTomTom(infile, outfile):
    '''load dreme tomtom results'''
    
    P.load(infile, outfile, options='-H "query_id,target_id,optimal_offset,p_value,e_value,q_value,overlap,query_consensus,targe_consensus,orientation" ')

def summarizeFimoGenerator():

    memes = glob.glob("meme.chip.dir/*/fimo_out*")

    jobs = []
    
    if len(memes) == 0:
        yield []

    n = 0
    for meme in memes:
        n = n + 1
        meme_dir = "/".join(meme.split("/")[0:2]) + "/fimo_out*/fimo.txt"
        outfile = "/".join(meme_dir.split("/")[0:2])+ "/" + meme_dir.split("/")[1] + "_fimo_summary.txt"
        fimo_dirs = glob.glob(str(meme_dir))

        if len(fimo_dirs)==1:
            fimo_dirs = ''.join(fimo_dirs)
        else:
            fimo_dirs = ' '.join(fimo_dirs)
        
        pairs = [ fimo_dirs, outfile ]
        jobs.append(pairs)
        
    unique_jobs = [list(x) for x in set(tuple(x) for x in jobs)] # remove duplicates from list

    for j in unique_jobs:
        infiles = j[0].split(" ")
        outfile = j[1]

        not_empty_infiles = []
        for infile in infiles:
            if os.path.isfile(infile):
                with open(infile, "r") as o:
                    n = 0
                    for line in o:
                        n = n + 1
                if n > 1:
                    not_empty_infiles.append(infile)

        if len(not_empty_infiles)>0:
            yield [not_empty_infiles, outfile]
        
@follows(loadDremeTomTom)    
@files(summarizeFimoGenerator)
def summarizeFimo(infiles, outfile):
    '''summarize fimo results, remove empty files'''

    infiles = ' '.join(infiles)

    tmp_dir = "$SCRATCH_DIR"

    statement = '''tmp=`mktemp -p %(tmp_dir)s`; checkpoint;
                   cat %(infiles)s | grep -v "#" | sort -k1,1n > $tmp; checkpoint;
                   [[ -s $tmp ]] && mv $tmp %(outfile)s || rm $tmp'''

    P.run()
    
@transform(summarizeFimo, suffix(".txt"), ".load")
def loadFimo(infile, outfile):

    P.load(infile, outfile,
           options='-H "pattern_name,sequence_name,start,stop,strand,score,p_value,q_value,matched_sequence" ')
    
@follows(runMemeChIP, loadMemeTomTom, loadDremeTomTom, loadFimo)
def runMemeAnalysis():
    pass

################################################
#####        Find motifs with HOMER        #####
################################################
@follows(runMemeAnalysis, mkdir("homer.chip.dir"))
@transform("meme.seq.dir/*.foreground.fasta",
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           r"homer.chip.dir/\1.homer.log")
def runHomerFindMotifs(infile, outfile):
    '''run Homer findMotifs.pl on fasta sequences'''

    outdir = outfile.replace(".homer.log", "")
    bfile = infile.replace(".foreground.fasta", ".background.fasta")

    statement = '''if [ ! -d %(outdir)s ]; 
                     then mkdir %(outdir)s; 
                     fi;
                   findMotifs.pl 
                     %(infile)s 
                     fasta 
                     %(outdir)s 
                     -fastaBg %(bfile)s
                     &> %(outfile)s''' % locals()

    P.run()

@follows(runHomerFindMotifs, mkdir("homer.genome.dir"))
@transform("meme.seq.dir/*foreground.bed",
           regex(r"meme.seq.dir/(.*).foreground.bed"),
           r"homer.genome.dir/\1.homer.log")
def runHomerFindMotifsGenome(infile, outfile):
    '''Run homer findMotifsGenome.pl with default settings and background'''

    # Homer requires input peaks to have strand (+/- : 0/1) in 6th column
    # Also, sequence length must be specified, default = 200

    tmp_dir = "$SCRATCH_DIR"
    outdir = outfile.replace(".homer.log", "")
    
    statement = '''peaks=`mktemp -p %(tmp_dir)s`; checkpoint;
                  awk 'BEGIN {OFS="\\t"} {print $0,0}' %(infile)s > $peaks; checkpoint; 
                  findMotifsGenome.pl 
                    $peaks
                    mm10
                    %(outdir)s
                    -h
                    -size 200
                    &> %(outfile)s; checkpoint;
                  rm $peaks'''
    
    P.run()

@follows(runHomerFindMotifs, mkdir("motifsCoverage.dir"))
@transform("meme.seq.dir/*.foreground.bed",
           regex(r"meme.seq.dir/(.*).foreground.bed"),
           r"motifsCoverage.dir/\1.motifCoverage.txt")
def annotatePeaks(infile, outfile):
    '''Annotate peaks w/ top discovered motifs for histogram plots'''

    # get original input peak file for motif search
    run = os.path.basename(infile).replace(".foreground.bed", "")
    name = run.split(".")[0]
    peak_file = "data.dir/" + name + ".bed"

    # # select top motifs to search for
    # motifs = "homer.chip.dir/" + run + "/homerResults/motif[1-6].motif"

    # search for all discovered motifs
    motifs = "homer.chip.dir/" + run + "/homerResults/motif*.motif"
    
    statement = '''annotatePeaks.pl 
                     %(peak_file)s
                     mm10 
                     -size 1000 
                     -hist 5 
                     -m %(motifs)s
                     > %(outfile)s'''
    P.run()
                  
@transform(annotatePeaks,
            regex(r"(.*).motifCoverage.txt"),
            r"\1.motifEnrichment.png")
# @transform("test.dir/*txt",
#            regex(r"(.*).motifCoverage.txt"),
#            r"\1.motifEnrichment.png")

def homerMotifEnrichmentPlot(infile, outfile):
    '''Process results from annotatePeaks and plot motif enrichment relative to peaks'''

    statement = '''python motifPlot.py %(infile)s %(outfile)s'''

    print(statement)
    
    P.run()
    
    # #job_memory = "10G"
    
    # homerResults = pd.read_csv(infile, sep="\t", header=0)

    # # limit no. motif to include in plot
    # limit_motifs = 5
    
    # # get plot title
    # title = infile.split("/")[-1]

    # # select total counts for motifs in sequences
    # motifs = homerResults[[x for x in homerResults.columns if "total" in x]]
    # motifs.columns = [x.split(":")[-1].split("/")[0] for x in motifs.columns] # simplify motif names

    # # get position information
    # motifDist = homerResults.iloc[:, 0:1]
    # motifDist.columns = ["dist_from_center"]

    # # rejoin dfs
    # processed_results = pd.concat([motifDist, motifs], axis=1)

    # # melt data
    # processed_results = processed_results.melt(id_vars=["dist_from_center"])
    # processed_results.columns = ["dist_from_center", "motif", "motifs_per_bp_per_peak"]

    # # set plot style
    # matplotlib.style.use(["seaborn-ticks", "seaborn-deep", "seaborn-notebook"])
    # palette = ["r", "g", "b", "c", "m", "y", "k", "w"]

    # n = 0
    # cols = {}
    # for i in processed_results["motif"].unique():
    #     n = n + 1

    #     if limit_motifs:
    #         if n > limit_motifs:
    #             break

    #     color = palette[n]
    #     cols[i] = color

    #     df = processed_results[processed_results["motif"] == i]

    #     plt.plot(df["dist_from_center"], df["motifs_per_bp_per_peak"], label=''.join(df["motif"].unique()))

    # plt.legend()
    # plt.title(title)
    # plt.savefig(outfile, dpi=300, format="png", facecolor='w', edgecolor='w',
    #             orientation='portrait', papertype=None, transparent=False)
    # plt.close() 
    
@follows(runHomerFindMotifs, runHomerFindMotifsGenome, annotatePeaks)
def runHomerAnalysis():
    pass

@follows(runMemeAnalysis, runHomerAnalysis)
def runMotifAnalysis():
    pass

#@follows(runMotifAnalysis)
@files(None, "*.nbconvert.html")
def report(infile, outfile):
    '''Generate html report on pipeline results from ipynb template(s)'''

    templates = PARAMS["report_path"]
    templates = templates.split(",")

    if len(templates)==0:
        print("Specify Jupyter ipynb template path in pipeline.ini for html report generation")
        pass

    for template in templates:
        infile = os.path.basename(template)
        outfile = infile.replace(".ipynb", ".nbconvert.html")
        nbconvert = infile.replace(".ipynb", ".nbconvert.ipynb")
        tmp = os.path.basename(template)
    
        statement = '''cp %(template)s . ; checkpoint;
                   jupyter nbconvert 
                     --to notebook 
                     --allow-errors 
                     --ExecutePreprocessor.timeout=None
                     --ExecutePreprocessor.kernel_name=python3
                     --execute %(infile)s; checkpoint ;
                   jupyter nbconvert 
                     --to html 
                     --ExecutePreprocessor.timeout=None
                     --ExecutePreprocessor.kernel_name=python3
                     --execute %(nbconvert)s; checkpoint;
                   rm %(tmp)s'''

        P. run()




###############################################################################
################ Known motif enrichment analysis (AME) #########################
###############################################################################
@follows(mkdir("meme.ame.dir"), runMemeAnalysis)
@transform(getMemeSequences,
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           r"meme.ame.dir/\1.memeame")
def runAme(infile,  outfile):
    '''Note importance of sorted input sequences...'''
    # FASTA files in meme.seq.dir are sorted by peak score
    
    nfgseq = FastaIterator.count(infile)      
    bfile = infile.replace("foreground.fasta","background.bfile") 
    outdir = outfile.replace(".memeame","")
    
    motifDb =  PARAMS["ame_motif_db"]
    
    job_memory = "10G"
    job_threads = "4"
    
    statement='''ame
                 --verbose 1
                 --oc %(outdir)s
                 --fix-partition %(nfgseq)s
                 --scoring totalhits
                 --bgfile %(bfile)s
                 --bgformat 2
                 --method fisher
                 --length-correct
                 <( cat %(infile)s )
                 <( cat %(motifDb)s )
                 > %(outfile)s
              ''' % locals()
    print(statement)
    P.run()


###############################################################################
############### Get Peaks associated with discovered motifs ###################
###############################################################################
@follows(runAme, mkdir("query_motifs.dir/mast.results.dir"))
@files(None, "query_motifs.dir/mast_motifs.txt")
def filterTFDatabases(infile, outfile):
    '''Filter TF databases for interesting motifs identified by meme-chip
       & s0pecified in pipeline.ini'''

    TFdb = PARAMS["mast_motif_db"]
    TFdb = TFdb.replace(",", " ")
    motifs = PARAMS["mast_motifs"].split(",")

    
    for TF in motifs:
    
        statement = '''awk '/^MOTIF/ {p=0} /%(TF)s/ {p=1} p' <( cat %(TFdb)s ) >> %(outfile)s''' % locals()

        print(statement)
        P.run()

@transform(filterTFDatabases, suffix(".txt"), r".meme")
def addMemeMotifHeader(infile, outfile):
    '''prepend header to query motif .meme file'''

    TFdb = PARAMS["mast_motif_db"].split(",")[0]

    statement = '''head -n9 %(TFdb)s | cat - %(infile)s > tmp && mv tmp %(outfile)s''' % locals()

    P.run()
    
@follows(addMemeMotifHeader)
@transform(getMemeSequences,
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           add_inputs("query_motifs.dir/mast_motifs.meme"),
           r"query_motifs.dir/mast.results.dir/\1.mast.log")
def runMast(infiles, outfile):
    '''Run mast for selected motifs (in pipeline.ini file, output hit coordinated at *.mast.txt'''
    
    job_threads = "4"
    job_memory = "5G"
    
    sequences, motifs = infiles
    
    bfile = sequences.replace(".foreground.fasta", ".background.bfile")
    outdir = outfile[:-len(".mast.txt")]
 
    statement = '''mast 
                   -w              
                   -bfile %(bfile)s
                   -oc %(outdir)s
                   <(cat %(motifs)s )
                   %(sequences)s
                   > %(outfile)s
                '''
    P.run()
    
@follows(runMast)
@transform(getMemeSequences,
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           add_inputs("query_motifs.dir/mast_motifs.meme"),
           r"query_motifs.dir/mast.results.dir/\1/\1.mast.hit_list.txt")
def runMast_HitList(infiles, outfile):
    '''Run mast for selected motifs (in pipeline.ini file, output hit coordinated at *.mast.txt'''
    
    job_threads = "4"
    job_memory = "5G"
    
    sequences, motifs = infiles
    
    bfile = sequences.replace(".foreground.fasta", ".background.bfile")

    head, tail = os.path.split(outfile)
    outdir = head

    statement = '''mast 
                   -hit_list
                   -best
                   -w              
                   -bfile %(bfile)s
                   -oc %(outdir)s
                   <(cat %(motifs)s )
                   %(sequences)s
                   > %(outfile)s
                '''
    P.run()
    
@follows(runMast_HitList)
@transform(filterTFDatabases, suffix(".txt"), r".table.txt")
def motifTable(infile, outfile):
    '''count motifs in .meme db'''

    statement = '''cat %(infile)s | grep MOTIF | tr -s " " | sed 's/ /\\t/g' | awk 'BEGIN {OFS="\\t"} {print $2,$3,NR}' > %(outfile)s'''

    P.run()

@transform(motifTable, suffix(".txt"), r".load")
def loadmotifTable(infile, outfile):
    P.load(infile, outfile, options='-H "motif_id,motif_name,motif_no" ')

@follows(loadmotifTable)
@transform(runMast_HitList, suffix(".txt"), r".table.txt")
def HitList_table(infile, outfile):
    '''make table for csvdb'''

    statement = '''sed 's/\([-+]\)/\\1\\t/' %(infile)s | grep -v ^# | tr -s " " | sed 's/ /\\t/g' | sort -k6,6rn > %(outfile)s'''
    print(statement)
    P.run()

@transform(HitList_table, suffix(".txt"), r".load")
def loadHitList_table(infile, outfile):
    P.load(infile, outfile, options=' -H "sequence_name,strand,motif_no,hit_start,hit_end,score,hit_p_value" ')

@follows(loadHitList_table)
@transform(runMast_HitList,
           regex(r"(.*).mast.hit_list.txt"),
           r"\1.MAST_HitList_results.txt")
def mastHitList_results(infile, outfile):
    '''Make informative table with mast hit list results'''

    hit_list_table = os.path.basename(infile).replace(".", "_")[:-len(".txt")] + "_table"
    motif_id_table = "mast_motifs_table"
    peak_info_table =  "_".join(os.path.basename(infile).split(".")[0:2])
    
    query = '''SELECT a.motif_id, a.motif_name, b.motif_no, b.sequence_name, b.score, 
            b.hit_p_value, c.contig, b.hit_start, b.hit_end, c.start, c.end, b.strand 
            FROM  %(motif_id_table)s a, %(hit_list_table)s b, %(peak_info_table)s c 
            WHERE a.motif_no = b.motif_no AND b.sequence_name = c.peak_id''' % locals()
    
    print(query)
    
    dbh = sqlite3.connect(PARAMS["database"])
    cc = dbh.cursor()
    sqlresult = cc.execute(query).fetchall()

    o = open(outfile,"w")
    o.write("\t".join ( 
            [ "motif_id", "motif_name", "motif_no", "peak_name", "score", "p_value", "contig", "motif_start"
              "motif_end", "peak_start", "peak_end", "strand"]) + "\n" )

    for r in sqlresult:
        motif_id, motif_name, motif_no, peak_name, score, p_value, contig, motif_start = r[0:8]
        motif_end, peak_start, peak_end, strand = r[8:12]

        columns = [ str(x) for x in 
                    [ motif_id, motif_name, motif_no, peak_name, score, p_value, contig, motif_start, motif_end, peak_start, peak_end, strand ] ]
        #               1             2         3       4           5       6      7       8              9         10         11         12
        o.write("\t".join( columns  ) + "\n")
    o.close()

@transform(mastHitList_results, suffix(".txt"), r".bed")
def mastHitList_bed(infile, outfile):
    '''make bedfile of motif locations'''

    # mast reports motif position relative to peak file. Need to return genomic coords of motifs & names of associated peak
    statement = '''awk 'BEGIN {OFS="\\t"} { if ($6 <= 0.005 && $12 = "+") {print $7,$10+$8,$10+$9,$4,$5,$12} 
                   if ($6 <= 0.005 && $12 = "-") {print $7,$11-$8,$11-$9,$4,$5,$12} }' <( tail -n +2 %(infile)s ) > %(outfile)s'''

    P.run()

@follows(mastHitList_bed, mkdir("query_motifs.dir/fimo.results.dir"))
@transform(getMemeSequences,
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           r"query_motifs.dir/fimo.results.dir/\1.fimo.txt")
def runFimo(infiles, outfile):
    '''Run fimo for selected motifs'''
    
    job_threads = "4"
    job_memory = "5G"
    
    sequences = infiles
    motifs = PARAMS["memechip_motif_db"].replace(",", " ")
    
    bfile = sequences.replace(".foreground.fasta", ".background.bfile")
    outdir = outfile[:-len(".fimo.txt")]

    # Parse-genomic-coord flag reportcoords of sequences where available
    
    statement = '''fimo
                   --bgfile %(bfile)s
                   -oc %(outdir)s
                   --parse-genomic-coord 
                   <(cat %(motifs)s )
                   %(sequences)s
                   > %(outfile)s
                '''
    P.run()
    
######################################
# MAST for motifs discovered by MEME #
######################################
@follows(runFimo, mkdir("mast.meme_results.dir"))
@transform(getMemeSequences,
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           add_inputs(r"meme.chip.dir/\1.memechip"),
           r"mast.meme_results.dir/\1.mast.log")
def runMemeMast(infiles, outfile):
    '''Search for meme identified motifs in sequence files'''

    # this searches for de novo meme motifs in "combined.meme" motif file

    job_threads = "2"
    job_memory = "8000M"
    sequences, motifs = infiles

    motifs = motifs.replace(".memechip", "/") + "combined.meme"
    bfile = sequences.replace(".foreground.fasta", ".background.bfile")
    outdir = outfile.replace(".mast.log","")

    # -comp option causes mast to crash, possible solution is here:
    # https://groups.google.com/forum/#!topic/meme-suite/bAiZCzJ36mo
    
    # -hit_list option would output a nice bed file but isn't working...

    # -remcorr removes highly correlated query motifs (lots of these are generated by meme-chip)
    
    statement = '''mast 
                   -bfile %(bfile)s
                   -oc %(outdir)s
                   -w
                   -remcorr
                   %(motifs)s %(sequences)s
                   > %(outfile)s
                '''
    
    P.run()

###############################################################################
########################### Pipeline Targets ##################################
###############################################################################
@follows(runMotifAnalysis, report)
def full():
    pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

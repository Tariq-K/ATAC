"""===========================
Pipeline MEME-ChIP
===========================

Overview
========

This pipeline runs MEME-ChIP and Homer for de novo motif discovery at provided bed files. Motif
discovery settings can be configured in pipeline.yml.

Ame can also be run optionally to perform motif enrichment analysis of known motifs. 

After motif discovery MAST can be run to rank input peaks by specified motifs (either de novo 
motifs discovered by meme, or database entries for known motifs, or both).


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

- Use the target "runMotifAnalysis" for de novo motif discovery using MEME-ChIP and Homer
    * Local background is constructed for MEME-ChIP comprising peak flanking regions
    * Homer is run with both local (results in homer.local.dir) and genomic (results in homer.genome.dir) background
- "runAme" for enrichment of known motifs with Meme
- Motifs of interest can then be speciifed in pipeline.yml and make "runMastAnalysis" target
  to run MAST on input peaks (or "full" target to run both tasks)
- Aternatively, if motif analysis was run elsewhere, MAST can be run on its own 
  to search for defined motif in new peaks etc.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_memechip.py config

Input files
-----------
- bed formatted files
- Pipeline is configured to accept output files from macs2 and pipeline_atac
    * "*.narrowPeak" (MACS2), "*.peaks.bed" or "*.merged.bed" (pipeline_atac)
    * or any bed file, with peak scores in column 5 (for sorting peaks)

Requirements
------------
cgat-flow,
cgat-core,
cgat-apps,

MEME venv:
- use conda to create a python 2 venv with which to run meme BEFORE running pipeline
- e.g. conda create -n meme python=2.7.14 meme=4.11.2
- specify the name of meme venv in pipeline.yml


Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""

from ruffus import *

import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgat.FastaIterator

import PipelineMemechip as S

import sys, tempfile 
import glob, os 
import subprocess
import sys
#import rpy2.robjects as R
import logging as L
import sqlite3
import pybedtools
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns

# -----------------------------------------------
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


################################################
#####              Meme-ChIP               #####
################################################
@follows(connect)
@transform("data.dir/*.peaks.bed",
           regex(r"(.*).peaks.bed"),
           r"\1.meme.bed")
def peakSummit(infile, outfile):
    '''Make bed of peak_summit (+/- 1bp), sorted by peak qvalue'''

    bed_format = PARAMS["infile_format"]
    
    if os.path.isfile(infile):
        df = pd.read_csv(infile, sep="\t", header=None)

        if bed_format == "macs2": 
            df.columns = ["contig", "start", "end", "peak_id", "score",
                          "strand", "FC", "pvalue", "qvalue", "rel_summit_pos"]

            # get peak summit position (start + relative summit position to peak start)
            df["summit"] = df["summit"] = df.apply(lambda x: int(x.start + x.rel_summit_pos), axis=1)
            
        elif bed_format == "deseq2":
            df.columns = ["contig", "start", "end", "peak_id", "score", "summit", "gene_name", "accessibility"]

        elif bed_format == "other":
            # cut to 5 cols & calculate centre of peak
            df = df.iloc[:, [0,1,2,3,4]]
            df.columns = ["contig", "start", "end", "peak_id", "score"]
            df["summit"] = df["summit"] = df.apply(lambda x: int(x.start+((x.end-x.start)/2)), axis=1)            

        else:
            print("Specify input bed file format in pipeline.yml")
            pass
        
        df = df.sort_values("score", ascending=False) # sort by score

        df["summit_start"] = df["summit"] - 1
        df["summit_end"] = df["summit"] + 1

        df = df[["contig", "summit_start", "summit_end", "peak_id", "score"]]

        df.to_csv(outfile, sep="\t", header=None, index=None)

        
def getMemeForegroundBedGenerator():
    
    beds = glob.glob("data.dir/*.meme.bed")

    print(PARAMS["memechip_npeaks"])
    
    if len(beds) == 0:
        pass

    npeaks = [str(x) for x in PARAMS["memechip_npeaks"]]
    widths = [str(x) for x in PARAMS["memechip_widths"]]

    outdir = "meme.seq.dir/"
    
    for bed in beds:
        int_name = os.path.basename(bed)[:-len(".meme.bed")]
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
        statement = f'''slopBed 
                         -i <(head -n {npeaks} {infile} )
                         -g {chrom_sizes}
                         -l {start_width} 
                         -r {end_width} 
                         -s 
                         > {outfile}''' 
        
    else:
        statement = f'''slopBed 
                          -i  {infile} 
                          -g {chrom_sizes} 
                          -l {start_width} 
                          -r {end_width}
                          -s 
                          > {outfile}'''
   
    P.run(statement)

    
@transform(getMemeForegroundBed,
           regex(r"(.*).foreground.bed"),
           r"\1.background.bed")
def getMemeBackgroundBed(infile, outfile):
    '''get bed file of peak flanking regions (of equal width to peak) for meme background model'''

    genome_idx = os.path.join(PARAMS["annotations_mm10dir"],"assembly.dir/contigs.tsv")
    
    statement = f'''sort 
                     -k1,1 
                     -k2,2n 
                     {infile} | 
                   bedtools flank 
                     -pct 
                     -l 1 
                     -r 1 
                     -g {genome_idx} 
                     > {outfile}'''
    
               # -pct -l 1 -r 1 will create flanking regions = peak width

    P.run(statement)
    

@follows(getMemeBackgroundBed)
@transform("meme.seq.dir/*.bed",
           regex(r"(.*).bed"),
           r"\1.fasta")
def getMemeSequences(infile, outfile):
    '''Get the peak sequences. The genome sequences are already repeat soft-masked'''
  
    genome_fasta = os.path.join(PARAMS["genome_dir"],PARAMS["genome"]+".fasta")
    genome_idx = os.path.join(PARAMS["annotations_mm10dir"],"assembly.dir/contigs.tsv")

    statement = f'''fastaFromBed 
                     -fi {genome_fasta} 
                     -bed {infile}
                     -fo {outfile}'''
    
    P.run(statement)
    

@follows(getMemeSequences)
@transform("meme.seq.dir/*background.fasta", suffix(".fasta"), ".bfile")
def getMemeBfiles(infile, outfile):
    '''prepare the meme background model'''

    statement = f'''fasta-get-markov 
                      -m 2 {infile}  
                      > {outfile}''' 

    P.run(statement)
    
               
@follows(getMemeSequences, getMemeBfiles, mkdir("meme.chip.dir"))
@transform("meme.seq.dir/*.foreground.fasta",
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           r"meme.chip.dir/\1.memechip")
def runMemeChIP(infile, outfile):
    '''run MemeChIP'''

    outdir = outfile.replace(".memechip", "")
    bfile = infile.replace(".foreground.fasta", ".background.bfile")
    motifDb =  " -db ".join(PARAMS["memechip_motif_db"]) # Meme-Chip needs each db in list to have "-db" flag
    env = PARAMS["memechip_env"]
    options = PARAMS["memechip_options"]
    
    statement = f'''nmeme=`wc -l {infile} | tr -s "[[:blank:]]" "\\t" | cut -f1` &&
                    meme-chip
                      -oc {outdir}
                      -db {motifDb}
                      -bfile {bfile}
                      -nmeme $nmeme
                      {options}
                      {infile}
                      > {outfile} ''' 

    # python3 version of tomtom does exist but isn't called by meme-chip
    # therefore run in python 2 conda env

    P.run(statement, job_condaenv=env, job_memory="2G", job_threads=5)

    
def loadMemeTomTomGenerator():

    meme_tomtom = glob.glob("meme.chip.dir/*memechip")

    if len(meme_tomtom) == 0:
        pass
        
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

    # cgat default is not to submit P.load() jobs to cluster
    # this is hardcoded in P.load()

    tablename = os.path.basename(outfile).replace(".load", "").replace(".", "_")
    options='-H "query_id,target_id,optimal_offset,p_value,e_value,q_value,overlap,query_consensus,targe_consensus,orientation" '

    statement = []
    statement.append('''cat %(infile)s | ''')
    statement.append(P.build_load_statement(tablename, options=options, retry=True) )
    statement.append(''' > %(outfile)s''')
    statement = ' '.join(statement)
    
    to_cluster = True

    P.run(statement)
    
    
def loadDremeTomTomGenerator():

    meme_tomtom = glob.glob("meme.chip.dir/*memechip")
    
    if len(meme_tomtom) == 0:
        pass
            
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

    to_cluster = True

    tablename = os.path.basename(outfile).replace(".load", "").replace(".", "_")
    options = '-H "query_id,target_id,optimal_offset,p_value,e_value,q_value,overlap,query_consensus,targe_consensus,orientation" '

    statement = []
    statement.append('''cat %(infile)s | ''')
    statement.append(P.build_load_statement(tablename, options=options, retry=True) )
    statement.append(''' > %(outfile)s''')
    statement = ' '.join(statement)
    
    P.run(statement)

    
def summarizeFimoGenerator():

    memes = glob.glob("meme.chip.dir/*/fimo_out*")

    jobs = []
    
    if len(memes) == 0:
        pass

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

    statement = f'''tmp=`mktemp -p {tmp_dir}` &&
                    cat {infiles} | 
                      grep -v "#" | 
                      sort -k1,1n > $tmp &&
                    [[ -s $tmp ]] && mv $tmp {outfile} || rm $tmp'''

    P.run(statement)

    
@transform(summarizeFimo, suffix(".txt"), ".load")
def loadFimo(infile, outfile):
    
    to_cluster = True

    tablename = os.path.basename(outfile).replace(".load", "").replace(".", "_")
    options='-H "pattern_name,sequence_name,start,stop,strand,score,p_value,q_value,matched_sequence" '

    statement = []
    statement.append(f'''cat {infile} | ''')
    statement.append(P.build_load_statement(tablename, options=options, retry=True) )
    statement.append(f''' > {outfile}''')
    statement = ' '.join(statement)
    
    P.run(statement)

    
@follows(runMemeChIP, loadMemeTomTom, loadDremeTomTom, loadFimo)
def runMemeAnalysis():
    pass


################################################
#####        Find motifs with HOMER        #####
################################################
@follows(runMemeAnalysis, mkdir("homer.chip.dir"))
@follows(mkdir("homer.chip.dir"))
@transform("meme.seq.dir/*.foreground.fasta",
           regex(r"meme.seq.dir/(.*).foreground.fasta"),
           r"homer.chip.dir/\1.homer.log")
def runHomerFindMotifs(infile, outfile):
    '''run Homer findMotifs.pl on fasta sequences'''

    outdir = outfile.replace(".homer.log", "")
    bfile = infile.replace(".foreground.fasta", ".background.fasta")

    statement = f'''if [ ! -d {outdir} ]; 
                      then mkdir {outdir}; 
                      fi;
                    findMotifs.pl 
                      {infile} 
                      fasta 
                      {outdir} 
                      -fastaBg {bfile}
                      &> {outfile}''' 

    P.run(statement)

    
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
   
    statement = f'''peaks=`mktemp -p {tmp_dir}` &&
                   awk 'BEGIN {{OFS="\\t"}} {{print $0,0}}' {infile} > $peaks && 
                   findMotifsGenome.pl 
                     $peaks
                     mm10
                     {outdir}
                     -h
                     -size given
                     &> {outfile} &&
                   rm $peaks'''
    
    P.run(statement)
    

@follows(runHomerFindMotifsGenome, mkdir("motifsCoverage.dir"))
@transform("meme.seq.dir/*.foreground.bed",
           regex(r"meme.seq.dir/(.*).foreground.bed"),
           r"motifsCoverage.dir/\1.motifCoverage.txt")
def annotatePeaks(infile, outfile):
    '''Annotate peaks w/ top discovered motifs for histogram plots'''

    # get original input peak file for motif search
    run = os.path.basename(infile).replace(".foreground.bed", "")
    name = run.split(".")[0]
    peak_file = "data.dir/" + name + ".bed"

    # update this so homer run to be used can be specified in pipeline.yml
    motifs = "homer.genome.dir/" + run + "/homerResults/motif*.motif"
    
    statement = f'''annotatePeaks.pl 
                      {peak_file}
                      mm10 
                      -size 1000 
                      -hist 5 
                      -m {motifs}
                      > {outfile}'''
    P.run(statement)
    
                  
@transform("motifsCoverage.dir/*.motifCoverage.txt",
            regex(r"(.*).motifCoverage.txt"),
            r"\1.motifEnrichment.png")
def homerMotifEnrichmentPlot(infile, outfile):
    '''Process results from annotatePeaks and plot motif enrichment relative to peaks'''

    statement = f'''python python/motifPlot.py 
                      {infile} 
                      {outfile}'''

    P.run(statement)
    
        
@follows(runHomerFindMotifs, runHomerFindMotifsGenome, annotatePeaks)
def runHomerAnalysis():
    pass


@follows(runMemeAnalysis)
@files(None, "*.nbconvert.html")
def report(infile, outfile):
    '''Generate html report on pipeline results from ipynb template(s)'''

    templates = PARAMS["report_path"]

    if len(templates)==0:
        print("Specify Jupyter ipynb template path in pipeline.ini for html report generation")
        pass

    for template in templates:
        infile = os.path.basename(template)
        outfile = infile.replace(".ipynb", ".nbconvert.html")
        nbconvert = infile.replace(".ipynb", ".nbconvert.ipynb")
        tmp = os.path.basename(template)
    
        statement = f'''cp {template} .  &&
                        jupyter nbconvert 
                          --to notebook 
                          --allow-errors 
                          --ExecutePreprocessor.timeout=None
                          --ExecutePreprocessor.kernel_name=python3
                          --execute {infile} &&
                        jupyter nbconvert 
                          --to html 
                          --ExecutePreprocessor.timeout=None
                          --ExecutePreprocessor.kernel_name=python3
                          --execute {nbconvert} &&
                        rm {tmp}'''

        P.run(statement)

        
@follows(runMemeAnalysis, runHomerAnalysis)
def runMotifAnalysis():
    pass

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
    
    bfile = infile.replace("foreground.fasta","background.bfile") 
    outdir = outfile.replace(".memeame","")
    
    motifDb =  PARAMS["ame_motif_db"]
    
    statement = f'''ame
                     --verbose 1
                     --oc {outdir}
                     --scoring totalhits
                     --bgfile {bfile}
                     --bgformat 2 
                     --method fisher
                     --length-correction
                     {infile}
                     {motifDb}
                     > {outfile}'''

    P.run(statement, job_memory="10G", job_threads=4)


###############################################################################
############### Get Peaks associated with discovered motifs ###################
###############################################################################
@follows(mkdir("meme.seq.dir"))
@follows(mkdir("query_motifs.dir/mast.results.dir"))
@active_if(bool(PARAMS["mast_meme_motif"]))
@files(None, "query_motifs.dir/denovo_motif.meme")
def getMemeMotif(infile, outfile):
    '''Get meme formatted motif from MEME-CHIP output.
       Specify motif to get in pipeline.ini (path and motif number)'''

    meme_out, motif_no, motif_name = PARAMS["mast_meme_motif"].split(",")
    max_motif = PARAMS["memechip_nmotif"]

    motif_start = "^MOTIF  " + str(motif_no)

    if int(motif_no) < max_motif:
        motif_end = "^MOTIF  " + str(int(motif_no)+1)
    else:
        motif_end = "^SUMMARY" 
        

    statement = f'''sed -n '/^MEME/,/^MOTIF  1/p' {meme_out} | 
                      head -n-1 > {outfile} ; 
                    sed -n '/{motif_start}/,/{motif_end}/p' {meme_out} | 
                      head -n -1 >> {outfile}'''

    P.run(statement)
    
    
@follows(getMemeMotif)
@files(None, "query_motifs.dir/db_motifs.txt")
def filterTFDatabases(infile, outfile):
    '''Filter TF databases for interesting motifs identified by meme-chip
       & specified in pipeline.ini'''

    TFdbs = PARAMS["mast_motif_db"]
    TFdb = ' '.join(TFdbs)
    motifs = PARAMS["mast_motifs"]

    
    for TF in motifs:
    
        statement = f'''awk '/^MOTIF/ {{p=0}} /{TF}/ {{p=1}} p' 
                          <( cat {TFdb} ) 
                          >> {outfile}'''

        P.run(statement)

        
@transform(filterTFDatabases, suffix(".txt"), r".meme")
def addMemeMotifHeader(infile, outfile):
    '''prepend header to query motif .meme file'''

    TFdb = PARAMS["mast_motif_db"][0]

    statement = f'''head -n9 {TFdb} | 
                      cat - {infile} 
                        > tmp && 
                      mv tmp {outfile}'''

    P.run(statement)

    
@follows(addMemeMotifHeader)
@transform("data.dir/*.peaks.bed",
           regex(r"data.dir/(.*).peaks.bed"),
           r"meme.seq.dir/\1.input_sequences.fasta")
def getInputPeakSequences(infile, outfile):
    '''Convert input bed files (all peaks) to fasta sequences for motif searching with MAST & FIMO'''

    genome_fasta = os.path.join(PARAMS["genome_dir"],PARAMS["genome"]+".fasta")

    statement = f'''fastaFromBed 
                      -name 
                      -fi {genome_fasta} 
                      -bed {infile} 
                      -fo {outfile}'''

    P.run(statement)

    
@active_if(PARAMS["mast_background"])
@transform("data.dir/*.peaks.bed",
           regex(r"data.dir/(.*).peaks.bed"),
           r"meme.seq.dir/\1.input_background.fasta")
def getInputPeakBackgroundFASTA(infile, outfile):
    '''get bed file of peak flanking regions (of equal width to peak) for meme background model'''

    genome_idx = os.path.join(PARAMS["annotations_mm10dir"],"assembly.dir/contigs.tsv")
    genome_fasta = os.path.join(PARAMS["genome_dir"],PARAMS["genome"]+".fasta")
    tmp_dir = "$SCRATCH_DIR"
    
    statement = f'''tmp=`mktemp -p {tmp_dir}` &&
                    sort -k1,1 -k2,2n {infile} | 
                      bedtools flank 
                        -pct 
                        -l 1 
                        -r 1 
                        -g {genome_idx} 
                        > $tmp &&
                    fastaFromBed 
                      -name 
                      -fi {genome_fasta} 
                      -bed $tmp 
                      -fo {outfile} &&
                    rm $tmp'''
    
    P.run(statement)
    
    
@active_if(PARAMS["mast_background"])
@transform(getInputPeakBackgroundFASTA, suffix(".fasta"), ".bfile")
def getInputPeakBackgroundModel(infile, outfile):
    '''prepare the meme background model'''

    statement = f'''fasta-get-markov 
                      -m 2 
                      {infile}  
                      > {outfile}'''

    P.run(statement)

    
@follows(getInputPeakSequences, getInputPeakBackgroundModel)
@transform("meme.seq.dir/*.input_sequences.fasta",
           regex(r"meme.seq.dir/(.*).input_sequences.fasta"),
           add_inputs(["query_motifs.dir/*.meme"]),
           r"query_motifs.dir/mast.results.dir/\1.mast.log")
def runMast(infiles, outfile):
    '''Run mast for selected motifs (in pipeline.ini file, output hit coordinated at *.mast.txt'''
    
    sequences, motifs = infiles

    background = PARAMS["mast_background"]

    if background == "custom":
        bfile = sequences.replace(".input_sequences.fasta", ".input_background.bfile")
        b_opt = f'''-bfile {bfile}''' 
    else:
        b_opt = ''' '''
        
    for motif in motifs:

        if len(motif) > 0:
            suffix = os.path.basename(motif).replace(".meme", "") 
            outdir = '.'.join([outfile[:-len(".mast.log")], suffix])

            statement = f'''mast 
                            -w              
                            {b_opt}
                            -oc {outdir}
                            <(cat {motif} )
                            {sequences}
                            &> {outfile} '''

            P.run(statement, job_memory="5G", job_threads=4)

            
@follows(runMast)
@transform("data.dir/*.peaks.bed",
           regex("data.dir/(.*).bed"),
           r"data.dir/\1.load")
def loadPeaks(infile, outfile):
    '''Load input peaks to merge w/ MAST results by peak_id'''

    tablename = os.path.basename(outfile).replace(".load", "").replace(".", "_")
    options='-H "contig,start,end,peak_id,score" '

    statement = []
    statement.append(f'''cat {infile} | cut -f1-5 - |''')
    statement.append(P.build_load_statement(tablename, options=options, retry=True) )
    statement.append(f''' > {outfile}''')
    statement = ' '.join(statement)
    
    to_cluster = True

    P.run(statement)

    
@follows(loadPeaks, mkdir("query_motifs.dir/mast.beds.dir/"))
@transform("query_motifs.dir/mast.results.dir/*/mast.txt",
           regex("query_motifs.dir/mast.results.dir/(.*)/mast.txt"),
           r"query_motifs.dir/mast.results.dir/\1/mast.load")
def loadMast(infile, outfile):
    '''MAST results section1 contains high scoring sequences (for input motifs), 
       ranked by increasing e-value (up to a max of 10)'''
    
    tmp_dir = "$SCRATCH_DIR"
    tablename = "MAST_" + outfile.split("/")[-2].replace(".", "_")
    options='-H "peak_id,e_value,length" '

    statement = []
    statement.append(f'''tmp=`mktemp -p {tmp_dir}` &&
                         sed -n '/^SECTION I:/,/^SECTION II:/p' {infile} | 
                           grep "^[a-zA-Z0-9].*[0-9]$" - | 
                           tr -s "[[:blank:]]" "\\t" 
                           > $tmp &&
                         cat $tmp ''')
    statement.append(P.build_load_statement(tablename, options=options, retry=True))
    statement = ' | '.join(statement) + f''' > {outfile} && rm $tmp'''

    P.run(statement)
    
    
@transform(loadMast,
           regex("query_motifs.dir/mast.results.dir/(.*)/mast.load"),
           r"query_motifs.dir/mast.beds.dir/\1.topMASTpeaks.bed")
def mast_table(infile, outfile):
    '''get mast results as bed file.
       e.g. top scoring peaks for query motifs'''

    mast_table = "MAST_" + infile.split("/")[2].replace(".", "_")
    peak_table = infile.split("/")[2].split(".")[0] + "_peaks"

    query = f'''select a.contig, a.start, a.end, a.peak_id, b.e_value 
               from {peak_table} a, {mast_table} b where a.peak_id = b.peak_id'''

    print(query)
    df = S.fetch_DataFrame(query, db)
    df.to_csv(outfile, header=False, index=None, sep="\t")
    

@follows(mast_table)
@transform("meme.seq.dir/*.input_sequences.fasta",
           regex(r"meme.seq.dir/(.*).input_sequences.fasta"),
           add_inputs(["query_motifs.dir/*.meme"]),
           r"query_motifs.dir/mast.results.dir/\1.mast.hit_list.log")
def runMast_HitList(infiles, outfile):
    '''Run mast for selected motifs (in pipeline.ini file), 
       output hit coordinated at *.mast.txt'''
    
    sequences, motifs = infiles

    background = PARAMS["mast_background"]

    if background == "custom":
        bfile = sequences.replace(".input_sequences.fasta", ".input_background.bfile")
        b_opt = f'''-bfile {bfile}''' 
    else:
        b_opt = ''' '''
         
    for motif in motifs:
        
        suffix = os.path.basename(motif).replace(".meme", "") 
        outdir = '.'.join([outfile[:-len(".mast.hit_list.log")], suffix])
        hlist = outdir + "/mast.hit_list.txt"

        statement = f'''mast 
                          -hit_list
                          -best
                          {b_opt}
                          -w              
                          -oc {outdir}
                          <(cat {motif} )
                          {sequences}
                          > {hlist}
                          2>> {outfile} '''
        
        P.run(statement, job_memory="5G", job_threads=4)
        

@follows(runMast_HitList)
@transform(filterTFDatabases, suffix(".txt"), r".table.txt")
def motifTable(infile, outfile):
    '''count motifs in .meme db'''

    statement = f'''cat {infile} | 
                     grep MOTIF | 
                     tr -s " " | 
                     sed 's/ /\\t/g' | 
                     awk 'BEGIN {{OFS="\\t"}} {{print $2,$3,NR}}' 
                     > {outfile}'''

    P.run(statement)
    

@transform(motifTable, suffix(".txt"), r".load")
def loadmotifTable(infile, outfile):
    P.load(infile, outfile, options='-H "motif_id,motif_name,motif_no" ')

    
@follows(loadmotifTable)
@transform(runMast_HitList,
           regex(r"query_motifs.dir/mast.results.dir/(.*).mast.hit_list.log"),
           r"query_motifs.dir/mast.results.dir/\1.mast.hit_list.touch")
def HitList_table(infile, outfile):
    '''make table for csvdb'''

    head, tail = os.path.split(infile)
    infile1 = head + "/" + tail.split(".")[0] + ".denovo_motif/" + "mast.hit_list.txt"
    infile2 = head + "/" + tail.split(".")[0] +".db_motifs/" + "mast.hit_list.txt"
    infiles = [infile1, infile2]

    for infile in infiles:

        out = infile.replace(".txt", ".table.txt")

        statement = []
        statement.append(f'''grep -v "#" {infile} | 
                               tr -s "[[:blank:]]" "\\t" | 
                               sort -k5,5rn | 
                               sed 's/\([+-]\)\([0-9]\)/\\1\\t\\2/' ''')
        
        if "denovo" in infile:
            ### check what needs to be specified here - it should not be necessary
            ### as motif files can be added into data.dir/ and should be named accordingly
            motif_name = PARAMS["mast_meme_motif"]
            cmd_suffix = f'''sed 's/\(\\t1\\t\)/\\t{motif_name}\\t/' ''' 
            statement.append(cmd_suffix)
            
        if os.path.exists(infile):
            statement = '|'.join(statement) + f''' > {out}'''
            P.run(statement)
            
        else:
            message = f'''{infile} does not exist''' 
            print(message)

    statement = f'''touch {outfile}'''
    
    P.run(statement)

    
@follows(HitList_table)
@transform("query_motifs.dir/mast.results.dir/*/mast.hit_list.table.txt",
           regex(r"query_motifs.dir/mast.results.dir/(.*)/mast.hit_list.table.txt"),
           r"query_motifs.dir/mast.results.dir/\1/mast.hit_list.table.load")
def loadHitList_table(infile, outfile):
    table = "MASThitlist_" + infile.split("/")[2]
    P.load(infile, outfile, tablename=table , options=' -H "peak_id,strand,motif_no,hit_start,hit_end,score,hit_p_value" ')

    
@follows(loadHitList_table)
@transform(runMast_HitList,
           regex(r"(.*).mast.hit_list.log"),
           add_inputs(["query_motifs.dir/*.meme"]),
           r"\1.mast.hit_list.results.touch")
def mastHitList_results(infiles, outfile):
    '''Make informative table with mast hit list results'''

    infile, motifs = infiles

    dfs = []
    n = 0
    for motif in motifs:
        if len(motif)>0:
            n = n + 1
            
            mname = os.path.basename(motif).replace(".meme", "")
            bedfile = os.path.basename(infile).replace(".mast.hit_list.log", "").replace(".", "_")

            hit_list_table = "MASThitlist_" + '_'.join([bedfile, mname])
            peak_table = bedfile + "_peaks"
            motif_id_table = "db_motifs_table"
            
            if mname == "denovo_motif":
                query = f'''select b.contig, b.start + a.hit_start as start, b.start + a.hit_end as end, 
                            a.peak_id, a.score as mast_score, a.strand, a.motif_no as motif_name, a.hit_p_value as p_value 
                            from {hit_list_table} a, {peak_table} b 
                            where a.peak_id = b.peak_id order by mast_score desc'''

            elif mname == "db_motifs":
                query = f'''select b.contig, b.start + a.hit_start as start, b.start + a.hit_end as end, 
                            a.peak_id, a.score as mast_score, a.strand, c.motif_name, a.hit_p_value as p_value 
                            from {hit_list_table} a inner join {peak_table} b 
                            on a.peak_id = b.peak_id inner join {motif_id_table} c 
                            on a.motif_no = c.motif_no order by mast_score desc''' 

            else:
                error_message = f"{motif} file not recognised" 
                print(error_message)

            if n == 1:
                res = S.fetch_DataFrame(query, db)
                res.columns = [ "contig", "motif_start", "motif_end", "peak_id", "score", "strand", "motif_name", "p_value" ]
            else:
                df = S.fetch_DataFrame(query, db)
                df.columns = [ "contig", "motif_start", "motif_end", "peak_id", "score", "strand", "motif_name", "p_value" ]
                res = res.append(df)

    mast_motifs = res["motif_name"].drop_duplicates()

    for m in mast_motifs:
        out_dir = "query_motifs.dir/mast.beds.dir/"
        filename = out_dir + '.'.join([bedfile, m]) + ".MASThitList.bed"
        df = res[res["motif_name"] == m]
        df.to_csv(filename, header=False, index=None, sep="\t")

    statement = f'''touch {outfile}'''

    P.run(statement)
    

@follows(mastHitList_results)
def runMastAnalysis():
    pass


###############################################################################
########################### Pipeline Targets ##################################
###############################################################################
@follows(runMotifAnalysis, runMastAnalysis, runAme)
def full():
    pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

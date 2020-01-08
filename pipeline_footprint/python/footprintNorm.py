#!usr/bin.python
import sys
from argparse import ArgumentParser
import pandas as pd
import os

# add parent dir to python path
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import PipelineFootprint as F

####### Parse commandline arguments
parser = ArgumentParser(prog="normaliseFootprint")
parser.add_argument("--infile", help="Bed file containing per base coverage of ATAC cutsites", required=True)
parser.add_argument("--outfile", help="Name of output file to be written", required=True)
parser.add_argument("--database", help="Name of pipeline database", required=True)
parser.add_argument("--binsize", help="bin size (b.p.)", required=True, type=int)
parser.add_argument("--bam_name", help="name of bam file", default=False)
args = parser.parse_args()

def footprintNorm(infile, outfile, binSize, db, bam_name):
    '''Normalise ATAC cutting frequency'''

    
    # prep df
    motif_cov = pd.read_csv(infile, compression="gzip", sep="\t")
    motif_cov.columns = ["chr", "p_start", "p_end", "peak_id", "position", "motif"]
    motif_cov = motif_cov[["peak_id", "position", "motif"]]
    
    # normalise for sequencing depth

    if bam_name: # get name of db sample name
        sample = bam_name
    else:
        sample = os.path.basename(outfile).split(".")[1]
        
    query = f'''select total_reads from readCounts where sample="{sample}" '''

    total_reads = F.fetch_DataFrame(query, db)["total_reads"].astype(int)
    
    motif_cov["motif"] = motif_cov["motif"].apply(lambda x: x/total_reads)

    # get frequency of motif per bin
    bins = [x for x in range(min(motif_cov["position"])+binSize, max(motif_cov["position"])+binSize, binSize)]
    motif_cov["bins"] = pd.cut(motif_cov["position"], bins)
    motif_cov = motif_cov.groupby(["bins"]).agg({"motif":"mean", "position":"mean"})
    motif_cov.reset_index(inplace=True)
    motif_cov = motif_cov[["position", "motif"]]
    
    motif_cov.to_csv(outfile, index=False)

footprintNorm(args.infile, args.outfile, args.binsize, args.database, args.bam_name)

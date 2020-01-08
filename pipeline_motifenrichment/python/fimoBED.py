#!usr/bin.python
import sys
from argparse import ArgumentParser
import pandas as pd
import os

# import inspect
# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0,parentdir)
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import PipelineMotifenrichment as m

####### Parse commandline arguments
parser = ArgumentParser(prog="fimoBED")
parser.add_argument("--infile", help="loadFimo outfile", required=True)
parser.add_argument("--outfiles", help="Name of results table and bed file", required=True)
parser.add_argument("--db", help="Pipeline database", required=True)
args = parser.parse_args()


def fimoBed(infile, outfiles, db):
    '''Make BED file of FIMO motif locations.
       FIMO filters by p-value. Matches with significant
       q-value exported to bed and results files'''
    
    fimo = ''.join([infile.split("/")[1], "_fimo_results"]).replace(".", "_")
    peaks = ''.join([infile.split("/")[1].split(".")[0], "_foreground"]).replace(".", "_")

    table, bed = outfiles.split(",")

    if "db_motifs" in infile:
        
        motifIDs = ''.join(infile.split("/")[1].split(".")[1])
    
        query = '''select b.chr, b.start + a.start as start, 
                   b.start + a.stop as end, a.sequence_name, 
                   a.score, a.strand, a.pattern_name, c.TF, 
                   a.matched_sequence, a.p_value, a.q_value 
                   from %(fimo)s a, %(peaks)s b, %(motifIDs)s c 
                   where a.sequence_name = b.peak_id 
                   and a.pattern_name = c.pattern_name''' % locals()

    else:
        query = '''select b.chr, b.start + a.start as start, 
                   b.start + a.stop as end, a.sequence_name, 
                   a.score, a.strand, a.pattern_name, 
                   a.matched_sequence, a.p_value, a.q_value 
                   from %(fimo)s a, %(peaks)s b 
                   where a.sequence_name = b.peak_id ''' % locals()
    print(query)
    df = m.fetch_DataFrame(query, db)

#    df = df[df["q_value"] < padj] # filter by FDR
# dont filter by FDR, as this is to stringent for motif enrichment plots
# only motif-sequence match (p-value (FIMO default)

    if "db_motifs" not in infile:
        motifID = infile.split("/")[1].split(".")[1]
        df["TF"] = motifID
        
        df.to_csv(table, sep="\t", index=False)
        df[["chr", "start", "end", "sequence_name", "score", "strand", "TF"]].to_csv(bed, sep="\t", header=None, index=False)

    else:
        df["TF"] = df.apply(lambda x: '_'.join([x.TF, x.pattern_name]), axis=1) # make sure multiple motifs for single TF aren't merged
        motifs = df["TF"].unique()

        if len(motifs) > 1:
            # if multiple TFs were searched for generate seperate BED & summary file for each one
            for motif in motifs:
                res = df[df["TF"] == motif]
                mtable = table.replace("_summary.txt", "_") + motif + "_summary.txt"
                mbed = bed.replace(".bed", "_") + motif + ".bed"
                res.to_csv(mtable, sep="\t", index=False)
                res[["chr", "start", "end", "sequence_name", "score", "strand", "TF"]].to_csv(mbed, sep="\t", header=None, index=False)

        df.to_csv(table, sep="\t", index=False)
        df[["chr", "start", "end", "sequence_name", "score", "strand", "TF"]].to_csv(bed, sep="\t", header=None, index=False)


# run job
fimoBed(args.infile, args.outfiles, args.db)

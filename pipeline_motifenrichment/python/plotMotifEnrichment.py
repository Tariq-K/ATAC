#!usr/bin.python
import sys
from argparse import ArgumentParser
import pandas as pd
import matplotlib
matplotlib.use('agg') # avoids tkinter error
import matplotlib.pyplot as plt
import seaborn as sns
import os

####### Parse commandline arguments
parser = ArgumentParser(prog="plotMotifEnrichment")
parser.add_argument("--infiles", help="Bed file containing per base coverage of motif instances", required=True, nargs="+")
parser.add_argument("--outfile", help="Name of output file to be written", required=True)
parser.add_argument("--motif", default="motif", help="List names of motifs, default is to get this from infile", required=False, nargs="+")
parser.add_argument("--region", default="region", required=False)
parser.add_argument("--norm", default="motif/peak", type=str, help='normalise: "motif/peak" or "motif/bp/peak"', required=False)
parser.add_argument("--window", default="500", type=int, help="no. b.p. either side of peak centre to plot", required=False)
parser.add_argument("--bins", default="10", type=int, help="Size of bins over which to average coverage. default=10", required=False)
parser.add_argument("--dpi", default="300", type=int, required=False)
parser.add_argument("--fast", default="False", type=str, help="Specify --fast True to speed up plots for pipeline", required=False)
parser.add_argument("--gzip", default="False", type=str, help="If files are compressed use --gzip True", required=False)
parser.add_argument("--ylims", help="List of min/max values for y axis", type=float, required=False, nargs="+")
parser.add_argument("--read_no", help="Total no. reads in BAM for ATAC cut site normalisation", type=int, required=False, nargs="+")
args = parser.parse_args()

######## normalise motif coverage & plot
n = 0
colours = ["light red", "windows blue", "dusty purple", "greyish", "amber", "faded green"]
pal = sns.xkcd_palette(colours)

if len(args.infiles) > len(pal):
    extra_colours = sns.xkcd_rgb.keys() # all 954 colours
    pal2 = sns.xkcd_palette(extra_colours)
    pal = pal + pal2
    
sns.set(style="whitegrid") # set seaborn theme

for i in args.infiles:
    n = n +1
    c = n -1

    # read infile
    if args.gzip == "True":
        motif_cov = pd.read_csv(i, compression="gzip", sep="\t")
    else:
        motif_cov = pd.read_csv(i, sep="\t")

    # get motif names
    if args.motif == "motif":
        motif = os.path.basename(i).split(".")[1]
    else:
        motif = args.motif[c]       

    # name cols
    motif_cov.columns = ["chr", "p_start", "p_end", "peak_id", "position", "motif"]

    # get no peaks for normalisation
    no_peaks = len(motif_cov)


    if args.fast == "True":
        # dont correct figure axis or zoom to specified window to speed up plotting
        motif_cov = motif_cov[["peak_id", "position", "motif"]]
        print(motif_cov.head())
        bins = [x for x in range(min(motif_cov["position"])+args.bins, max(motif_cov["position"])+args.bins, args.bins)]
        motif_cov["bins"] = pd.cut(motif_cov["position"], bins)

        # motif frequency normalisation
        # Natoli lab just normalise per per peak - as bin width is constant
        motif_cov = motif_cov.groupby(["bins"]).agg({"motif":"mean", "position":"mean"}) # get mean motifs per bin = motifs/peak

        if args.norm == "motif/bp/peak":
            # or normalise as in HOMER => motif per bp per peak
            motif_cov["motif"] = motif_cov["motif"]/bins # normalise for bin width -> motif/bp/peak
        
        motif_cov.reset_index(inplace=True)

        plt.plot(motif_cov["position"], motif_cov["motif"], color=pal[c], label=motif)
        
    else:
        # get bp position relative to peak centre (rather than peak start)
        motif_cov["pwidth"] = max(motif_cov["position"])
        motif_cov["pcentre"] = motif_cov["pwidth"]/2
        motif_cov["bp_from_centre"] = ""
        motif_cov["bp_from_centre"] = motif_cov.apply(lambda x: x.position - x.pcentre if x.position >= x.pcentre else (x.pcentre - x.position)*-1, axis=1)
        motif_cov["bp_from_centre"].astype(int)

        # subset on peak centre +/- n b.p.
        motif_cov = motif_cov[motif_cov.apply(lambda x: x.bp_from_centre <= args.window and x.bp_from_centre >= -args.window, axis=1)]

        # normalise data
        motif_cov = motif_cov[["peak_id", "bp_from_centre", "motif"]]    
        motif_cov = motif_cov.groupby(["bp_from_centre"]).agg({"motif":"mean"}) # get frequency of motif per bp
        motif_cov["motif"] = motif_cov["motif"]/no_peaks # normalise for no. peaks -> motif/bp/peak
        motif_cov.reset_index(inplace=True)

        plt.plot(motif_cov["bp_from_centre"], motif_cov["motif"], color=palette[c], label=motif)

        # this is very slow... I should use awk to correct axis labels and filter bed by region distance cutoff instead
        # then pipe modified bed to this script (or save as tmp file)

if args.ylims:
    plt.ylim(args.ylims[0], args.ylims[1])
    
plt.xlim(min(motif_cov["position"]), max(motif_cov["position"]))
plt.title(args.region)
plt.legend()
plt.savefig(args.outfile, dpi=args.dpi)
plt.close()



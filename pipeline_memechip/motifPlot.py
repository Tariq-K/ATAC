import sys
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns


def homerMotifEnrichmentPlot(infile, outfile):
    '''Process results from annotatePeaks and plot motif enrichment relative to peaks'''
    
    #job_memory = "10G"

    homerResults = pd.read_csv(infile, sep="\t", header=0)

    # limit no. motif to include in plot
    limit_motifs = 5

    # get plot title
    title = infile.split("/")[-1]

    # select total counts for motifs in sequences
    motifs = homerResults[[x for x in homerResults.columns if "total" in x]]
    motifs.columns = [x.split(":")[-1].split("/")[0] for x in motifs.columns] # simplify motif names

    # get position information
    motifDist = homerResults.iloc[:, 0:1]
    motifDist.columns = ["dist_from_center"]

    # rejoin dfs
    processed_results = pd.concat([motifDist, motifs], axis=1)

    # melt data
    processed_results = processed_results.melt(id_vars=["dist_from_center"])
    processed_results.columns = ["dist_from_center", "motif", "motifs_per_bp_per_peak"]

    # set plot style
    matplotlib.style.use(["seaborn-ticks", "seaborn-deep", "seaborn-notebook"])
    palette = ["r", "g", "b", "c", "m", "y", "k", "w"]

    n = 0
    cols = {}
    for i in processed_results["motif"].unique():
        n = n + 1

        if limit_motifs:
            if n > limit_motifs:
                break

            color = palette[n]
            cols[i] = color

            df = processed_results[processed_results["motif"] == i]

            plt.plot(df["dist_from_center"], df["motifs_per_bp_per_peak"], label=''.join(df["motif"].unique()))

    plt.legend()
    plt.title(title)
    plt.savefig(outfile, dpi=300, format="png", facecolor='w', edgecolor='w',
			 orientation='portrait', papertype=None, transparent=False)
    plt.close()
    

if __name__ == "__main__":    
    homerMotifEnrichmentPlot(sys.argv[1], sys.argv[2])

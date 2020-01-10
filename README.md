# ATAC
Pipelines for ATAC-seq analysis

## pipeline_atac
This is for generic ATAC-seq analysis. Handles mapping, quality control metrics, peak calling, differential accessibility testing with DESeq2, and creating bigWig tracks for visualisation. 

## pipeline_memechip
Runs meme-chip and homer for de novo motif discovery. Also, optionally runs meme suit tools ame and mast to search for instances of known motifs or meme-chip results.

## pipeline_motifenrichment
Fimo scanning for known motifs. Creates per base pair matrix of motif occurences and plots motif frequencies over input intervals. 

## pipeline_superenhancer
Runs ROSE-style super enhancer analysis using ATAC-seq peaks as input

## pipeline_footprint
Requires high sequencing depth. Calculates cut site frequency over intervals (e.g. motif sites)

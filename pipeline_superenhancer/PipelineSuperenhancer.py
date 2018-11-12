import sys
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.Database as DB
import PipelineSuperenhancer as SE
import glob
import pandas as pd
import gzip

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


### Helper Functions

def writeGreat(locations,basalup,basaldown,maxext,outfile,half=False):
    ''' write out a bed file of great promoters from input gene locations
         locations is [contig,gstart,gend,strand,gene_id] '''

    # Gene regulatory domain definition: 
    # Each gene is assigned a basal regulatory domain of a 
    # minimum distance upstream and downstream of the TSS 
    # (regardless of other nearby genes). 
    # The gene regulatory domain is extended in both directions 
    # to the nearest gene's basal domain but no more than the 
    # maximum extension in one direction

    genome = {}
    for location in locations:
        chrom, gstart, gend, strand_int, gid = location
        if strand_int == "-": 
            strand = "minus" 
            tss = gend
        else: 
            strand = "plus"
            tss = gstart
        record = [tss,strand,gid]
        if chrom[3:5]=="NT" or chrom[3:]=="M": continue
        if chrom not in genome: 
            genome[chrom] = [ record ]
        else: genome[chrom].append(record)

    #add the ends of the chromosomes
    contigs = gzip.open(PARAMS["annotations_dir"]+"/assembly.dir/contigs.bed.gz","rt")

    
    nmatched = 0
    for contig_entry in contigs:
        contig, start, end = contig_entry.strip().split("\t")
        
        if contig in genome.keys():
            genome[contig].append([int(end),"end","end"])
            nmatched+=1
    if nmatched < 21:
        raise ValueError("writeGreat: not enough chromosome ends registered")

    #sort the arrays
    for key in genome.keys():
        genome[key].sort( key = lambda entry: entry[0] )
        
    #now we can walk over the regions and make the regulatory regions.

    greatBed = []
   
    for contig in genome.keys():

        locs = genome[contig]
        contig_end = locs[-1][0]
        for i in range(0,len(locs)):

            l,strand,gid = locs[i]

            if strand == "end": continue

            #get the positions of the neighbouring basal domains.

            # - upstream
            if i == 0: frontstop = 0
            else:
                pl, pstrand, pgid = locs[i-1]
                if pstrand == "plus": frontstop = pl + basaldown
                else: frontstop = pl + basalup
            # - downstream
            nl, nstrand, ngid = locs[i+1]
            if nstrand == "plus": backstop = nl - basalup
            else: backstop = nl - basaldown

            # define basal domain
            if strand=="plus":
                basalstart = l - basalup
                basalend = min( l + basaldown, contig_end )
            else:
                basalstart = l - basaldown
                basalend = min( l + basalup, contig_end )

            # upstream extension
            if frontstop > basalstart:
                regstart = basalstart
            else:
                if half == True:
                    frontext = min( maxext, (l - frontstop) / 2 )
                else:
                    frontext = min( maxext, l - frontstop )
                regstart = l - frontext

            # downstream extension
            if backstop < basalend:
                regend = basalend
            else:
                if half == True:
                    backext = min( maxext, ( backstop - l ) / 2 )
                else:
                    backext = min( maxext, backstop - l )
                regend = l + backext

            # greatBed.append(["chr"+contig,str(regstart),str(regend),gid])
            greatBed.append([contig,str(regstart),str(regend),gid])
        
    outfh = open(outfile,"w")
    outfh.write("\n".join(["\t".join(x) for x in greatBed])+"\n")
    outfh.close()


def getTSS(start,end,strand):
    if strand == 1 or strand == "+": tss = start
    elif strand == -1 or strand == "-": tss = end
    else: raise ValueError("getTSS: stand specification not understood")
    return tss

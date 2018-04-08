################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: PipelineGO.py 2877 2010-03-27 17:42:26Z andreas $
#
#   Copyright (C) 2013 Steve Sansom
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
#################################################################################
"""
=======================================================
PipelineChartseq.py - helper functions for this project
=======================================================

:Author: Steve Sansom
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------


Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

"""
# collections, itertools, optparse, os, shutil,sys, tempfile,  shutil, re
# random
#import Database
#import csv
#import IOTools
#import GTF, GFF, Bed
#import numpy
#import pybedtools

import math, glob  
import rpy2.robjects as R
import CGAT.Experiment as E
import logging as L
import sqlite3
import gzip
from bx.intervals.intersection import Interval, IntervalTree
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineGO as PipelineGO
import pandas as pd
import pybedtools

# P.getParameters( 
#      ["pipeline_chartseq/pipeline.ini",
#       "../pipeline.ini",
#       "pipeline.ini" ] )

# PARAMS = P.PARAMS
# PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
#                                        "pipeline_annotations.py" )


try:
    PARAMS = P.getParameters()
except IOError:
    pass

#######################################################################
########### Key helper functions ######################################
#######################################################################

def execute(queries, database=PARAMS["database"], attach=False):
    '''Execute a list of statements sequentially'''

    dbhandle = sqlite3.connect( database )
    cc = dbhandle.cursor()

    if attach:
        for attach_statement in attach:
            cc.execute(attach_statement)

    for statement in queries: cc.execute(statement)
    cc.close()

def fetch(query, database=PARAMS["database"], attach=False):
    '''Fetch all query results and return'''
    dbhandle = sqlite3.connect( database )
    cc = dbhandle.cursor()
    if attach:
        for attach_statement in attach:
            cc.execute(attach_statement)
    sqlresult = cc.execute(query).fetchall()
    cc.close()
    return sqlresult

def fetch_with_names(query, database=PARAMS["database"], attach=False):
    '''Fetch all query results and returns them as a dictionary'''
    dbhandle = sqlite3.connect( database )
    cc = dbhandle.cursor()
    if attach:
        for attach_statement in attach:
            cc.execute(attach_statement)
    sqlresult = cc.execute(query).fetchall()
    data=[]
    # http://stackoverflow.com/questions/4147707/
    # python-mysqldb-sqlite-result-as-dictionary
    field_names = [ d[0] for d in cc.description ]
    data.append( [ name for name in field_names ])
    for record in sqlresult:
        line = [ field for field in record ]
        data.append(line)

    cc.close()
    return data



def write(outfile, lines, header=False):
    ''' expects [[[line1-field1],[line1-field2 ] ],... ]'''
    handle = open(outfile,"w")

    if header:
        handle.write("\t".join([str(title) for title in header])+"\n")

    for line in lines:
        handle.write("\t".join([str(field) for field in line])+"\n")

    handle.close()

##########################################################################
##########################################################################
##########################################################################

##########################################################################
####### Functions for defining promoters #################################
##########################################################################

def makePromoterDictionary(sqlresult,name):
    ''' return a dictionary of promoters from an sql result '''

    annotations = {} #gene and transcript
    #plen = PARAMS["promoter"]
    us = 5000
    ds = 2500

    for r in sqlresult:
        chrom, gstart, gend, strand, gid#, tstart, tend, tid = str(r[0]), int(r[1]), int(r[2]), int(r[3]), str(r[4]), int(r[5]), int(r[6]), str(r[7])
        if chrom[0:2]=="NT" or chrom[0:2]=="MT":
            continue
        if strand==1:
            gps = gstart - us
            gpe = gstart + ds
            #tps = tstart - us
            #tpe = tstart + ds
        else:
            gps = gend - ds
            gpe = gend + us
            #tps = tend - ds
            #tpe = tend + us
        if gid not in annotations:
            annotations[gid]={"chrom":chrom,"start":gps,"end":gpe,"name":name}
    
    return annotations

def makeGreedyPromoterDictionary(sqlresult,name,lim=550000):
    #depreceated, use writeGreat

    annotations = {} #gene and transcript
    #plen = PARAMS["promoter"]
    ds = 550000
    genome = {}
    for r in sqlresult:
        chrom, gstart, gend, strand, gid, tstart, tend, tid = str(r[0]), int(r[1]), int(r[2]), int(r[3]), str(r[4]), int(r[5]), int(r[6]), str(r[7])
        if chrom[0:2]=="NT" or chrom[0:2]=="MT": continue
        if chrom not in genome: genome[chrom] = IntervalTree()
        if strand == 1 : 
            start = gstart 
            end = gend
        else: 
            start = gend
            end = gstart
        genome[chrom].insert(gstart,gend, { "start":start,"end":end, "id":gid,"strand":strand } )
    
    for r in sqlresult:
        chrom, gstart, gend, strand, gid, tstart, tend, tid = str(r[0]), int(r[1]), int(r[2]), int(r[3]), str(r[4]), int(r[5]), int(r[6]), str(r[7])
        #ideally we want lim up and 2.5 down.
        if chrom[0:2]=="NT" or chrom[0:2]=="MT": continue
        if strand==1:
            gps = gstart - lim
            gpe = gstart + ds
            records = genome[chrom].find(gstart-(2*lim),gstart+(2*ds))
            if len(records) == 0:
                annotations[gid]={"chrom":chrom,"start":gps,"end":gpe,"name":name}                
            else:
                #we want to find the closest start.
                upstream, downstream = [], []
                for r in records:
                    if r["id"] == gid: continue
                    s = r["start"]
                    if s < gstart:
                        upstream.append(s)
                    else:
                        downstream.append(s)
                
                if len(upstream) > 0:
                    uss = max(upstream)
                    if uss > gstart:
                        raise ValueError("start greater than start")
                    usep = gstart - uss
                    if usep < 2 * lim:
                        gps = gstart - usep/2
                
                if len(downstream) > 0:
                    dss = min(downstream)
                    dsep = dss - gstart
                    if dsep < 2 * ds:
                        gpe = gstart + dsep/2

                annotations[gid]={"chrom":chrom,"start":gps,"end":gpe,"name":name}                
        else:
            gps = gend - ds
            gpe = gend + lim
            records = genome[chrom].find(gend - (2*ds), gend + (2*lim) )
            if len(records) == 0:
                annotations[gid]={"chrom":chrom,"start":gps,"end":gpe,"name":name}                
            else:
                #we want to find the closest start.
                upstream, downstream = [], []
                for r in records:
                    s = r["start"]
                    if r["id"]== gid:continue
                    if s > gend: downstream.append(s)
                    else: upstream.append(s)
                    
                if len(downstream) > 0:
                    dss = min(downstream)
                    dsep = dss - gend
                    if dsep < 2 * lim:
                        gpe = gend + dsep/2
                
                if len(upstream) > 0:
                    uss = max(upstream)
                    usep = gend - uss
                    if usep < 2 * ds:
                        gps = gend - usep/2
                
                annotations[gid]={"chrom":chrom,"start":gps,"end":gpe,"name":name}                
                
    return annotations

def writeAnnotations(annotations,outfile,ids=False):
    fh = open(outfile,"w")            
    for k,v in annotations.iteritems():
        chrom = "chr" + str(v["chrom"])
        start = str( v["start"] )
        end = str( v["end"] )
        name = str( v["name"] )
        #name = k | because of the way gat works, this needs to identify the experiment
        if not ids:
            entry = [chrom, start, end, name]
        else:
            entry = [chrom, start, end, k]
        fh.write("\t".join(entry)+"\n")
    fh.close()


####
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
        if strand_int == -1: 
            strand = "minus" 
            tss = gend
        else: 
            strand = "plus"
            tss = gstart
        record = [tss,strand,gid]
        if chrom[3:5]=="NT" or chrom[0:2]=="MT": continue
        if chrom not in genome: 
            genome[chrom] = [ record ]
        else: genome[chrom].append(record)

    #add the ends of the chromosomes
    contigs = gzip.open(PARAMS["annotations_dir"]+"/assembly.dir/contigs.bed.gz","r")

    
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

            greatBed.append(["chr"+contig,str(regstart),str(regend),gid])
        
    outfh = open(outfile,"w")
    outfh.write("\n".join(["\t".join(x) for x in greatBed])+"\n")
    outfh.close()

def biomart_iterator( attributes,
                      filters,
                      values,
                      host,
                      biomart, 
                      dataset ):
    ''' stolen from pipeline biomart... '''

    r = R.r
    r.library("biomaRt")
    mart = r.useMart( biomart,
                      dataset=dataset,
                      host=host,
                      path="/biomart/martservice",
                      archive=False )

    values_list = values
    result = r.getBM( attributes=R.StrVector(attributes), 
                      filters=R.StrVector(filters),
                      values=values_list,
                      mart=mart )

    # result is a dataframe.
    # rx returns a dataframe.
    # rx()[0] returns a vector
    for data in zip( *[ result.rx(x)[0] for x in attributes] ):
        yield dict( zip(attributes, data) )


##########################################################################
##########################################################################
##########################################################################

##########################################################################
########### Functions for expression isochores ###########################
##########################################################################


def mergedPromoters(promoters):
    # ensure the list is sorted by contig and start.
    promoters.sort(key = lambda x: ( x[0], x[1]))

    cc = None

    for p in promoters:

        contig, start, end = p[:3]
        tss = p[3]
        iso = p[4]

        if start > end:
            raise ValueError("what is first cannot also be last")

        if cc == None:
            cc, cs, ce, ct, ci = contig, start, end, [tss], [iso]
            continue

        if start > ce or contig != cc:
            region = [cc, cs, ce, ct, ci] 
            cc, cs, ce, ct, ci = contig, start, end, [tss], [iso]            
            yield region
            # yield does not imply continue, hard won knowledge ;)
            continue 

        if end < ce : 
            #add the tss and iso.
            ct.append(tss)
            ci.append(iso)  
            continue
        else: 
            ce = end
            ct.append(tss)
            ci.append(iso)
            continue

    yield region  


def assignIsochores(segment,tss,iso):
    if len(iso) != len(tss): raise ValueError("no tss's does not match number of isochores")

    # the tss's are not necessarily in ascending order
    tss_iso = zip(tss,iso)
    tss_iso.sort(key = lambda k: k[0])

    stss = [x[0] for x in tss_iso]
    siso = [x[1] for x in tss_iso]

    start, end = segment

    if start > end: raise ValueError("bug 2")
    points = [start]

    for n in range(0,len(stss)-1):
        points.append(( stss[n] + stss[n+1] ) / 2)

    points.append(end)

    results = []
    for j in range(0,len(points)-1):
        if j==0:
            s, e = points[j], points[j+1]
        else:
            s, e = points[j] + 1, points[j+1]   
        if s > e + 1: 
            raise ValueError("bug is here 1")
        results.append([s, e, siso[j] ])

    return results


def split_isochores(data):            
    split_isochores = []
    for seg in mergedPromoters(data):
        if seg[2] < seg[1]:
            raise ValueError("end before start")          

        if len(seg[3])==1:
            split_isochores.append( seg[0:3] + [ seg[4][0] ] )

        else:
            
            ais = assignIsochores((seg[1],seg[2]),seg[3],seg[4])
            for ai in ais:
                if ai[1] < ai[0]:
                    raise ValueError("end before start")
                
                split_isochores.append( [seg[0]] + ai )

    return split_isochores
    

def stitch_isochores(split_isochores):
    #now join back together adjacent segments from the same isochore
    ce = None
    final_isochores = []
    for si in split_isochores:
        c, s, e, i = si
        if ce == None:
            cc, cs, ce, ci = c, s, e, i
            continue
        if s - ce <= 1 and i==ci and c==cc:
            ce = e
            continue
        else:
            region = [cc, cs, ce, ci]
            cc, cs, ce, ci = c, s, e, i
            final_isochores.append(region)
    
    return final_isochores



def getTSS(start,end,strand):
    if strand == 1 or strand == "+": tss = start
    elif strand == -1 or strand == "-": tss = end
    else: raise ValueError("getTSS: stand specification not understood")
    return tss

##########################################################################
##########################################################################
##########################################################################

##########################################################################
###### Functions for summarising parameterised GAT runs ##################
##########################################################################

def gutted(gat_result,columns):
    result = {}
    fh = open(gat_result,"r")
    for line in fh:
        fields = line.split("\t")
        if fields[1] == "annotation": 
            field_map = {}
            for i in range(0,len(fields)):
                field_map[fields[i]]=i
                
            continue # actually should construct the mapping
           
        #array = fields[1]
        #obs = fields[2]
        #enr = fields[7]
        #pval = fields[9]
        #qval = fields[10]
        #result[array]=(obs,enr,pval)
        name = fields[field_map["annotation"]]

        result[name] = dict([(x, fields[field_map[x]]) for x in columns])
                                         
    fh.close()
    return result



def summariseGATRunsByAnnotation(infile, outfile):
    '''
       Build a set of per annotation summary tables from a set of paramaterised 
       gat runs.
 
       Expects a parent directory with an informatively named
       subdirectory for each parameter set.
    '''

    gatDir =  infile
    outDir = outfile.split("/")[0]

    counterPatterns = ["segment-overlap","nucleotide-overlap"]

    
    run_sets = glob.glob(gatDir+"/*")

    print run_sets
                  
    gat_runs = {}

    names, segments_names, counter_types, annotations = [], [], [], []

    columns = {"observed":"obs",
               "fold":"fc",
               "pvalue":"pval",
               "qvalue":"qval"}

    for counter_type in counterPatterns:
        print ">>>>"
        print counter_type
        counter_types.append(counter_type)
        gat_runs[counter_type] = {}

        for dir_name in run_sets:
            dir_name_files = glob.glob(dir_name+"/*")

            #print dir_name_files
            for f in dir_name_files:
                basedir, name, filename = f.split("/")
                print filename

                segments = filename.split(".tsv.gz")[0]
                if name not in names: names.append(name)
                if segments not in segments_names: segments_names.append(segments)

                if filename.endswith(counter_type):
                    print "boom"
                    if name not in gat_runs[counter_type]: gat_runs[counter_type][name] = {}
                    gat_runs[counter_type][name][segments] =  gutted(f,columns.keys())

                    for annotation in gat_runs[counter_type][name][segments].keys():
                        if annotation not in annotations: annotations.append(annotation)

    #now we have the lot we need to invert the tree:
    annotations.sort()
    segments_names.sort()
    names.sort(reverse=True)
    #maxtn = max([ len(x) for x in segments_names ])
    #def getpad(n,w,char="."):
    #    out = n
    #    while len(out)<w: out+=char
    #    return out
                  
    def formatData(annotation):

        obs = annotation["observed"]
        fc = "%.2f" % float(annotation["fold"])
        pval = "%.3g" % float(annotation["pvalue"])
        qval = "%.3g" % float(annotation["qvalue"])
        if float(annotation["qvalue"])<0.05: sig="*"
        else: sig="_"
        fcv = float(annotation["fold"])
        if fcv <= 1/1.5:
            enr = "-"
        elif fcv >= 1.5: 
            enr= "+"
        else:
            enr = "~"
                  
        return([obs,fc,pval,qval,sig,enr])


    for counter_type in counter_types:
                      
        for annotation in annotations:
            outfileName = ".".join ([annotation, counter_type, "txt"])
            outfile = outDir + "/" + outfileName

            header = ["segments"]

            for name in names:
                header += [ name + "_" + columns[x] for x in columns.keys() ]

            header += ["enr", "sig"]

            lines = [ ]
            for segments in segments_names:
                line = [segments]
                sig, enr = "", ""
                for name in names:

                    res = formatData(gat_runs[counter_type][name][segments][annotation])
                    line += res[0:4]
                    sig += res[4]
                    enr += res[5]
                line += [enr, sig]

                lines.append(line)

            write(outfile, lines, header)


### were given a track and the paramaterised run names


def summariseGATRunsBySegments(infiles, outfile, params):
    '''Summarise a set of parameterised GAT runs by Segment'''
    

    paramDirs = params
    segments = infiles[0].split("/")[1][:-len(".bed")].replace(".","_")
    goTable = P.toTable(infiles[1])


    names = [ x.split("/")[1] for x in paramDirs ]

    #r, names = {}, []
    r = {}
    counts = {"gc":0,"exprs":0,"noname":0}

    for name in names:
        #if name not in names: names.append( name )
        table = "gatgo_" + name + "_" + segments
        statement = '''select gat.observed, gat.expected, gat.fold, gat.l2fold,
                              gat.pvalue,gat.qvalue,go.go_id,go.description
                       from %s gat
                       inner join %s go
                          on gat.annotation=go.go_id
                    ''' % (table, goTable) 
        results = fetch(statement)

        for result in results:
            
            obs, expect, fold, l2fold, pval, qval, go_id, desc = result[0:8]
            if go_id not in r: r[go_id] = {"desc":desc}
            if "name" not in r[go_id]: r[go_id]["name"] = {}
            r[go_id]["name"][name]={ "obs":obs,"expect":expect, "fold":fold, "l2fold":l2fold, "pval":pval, "qval":qval }
          
    names.sort()

    cols = ["fold","obs","expect","qval","pval"]

    header = ["go_id"]
    for name in names:
        header += [ name + "_" + col for col in cols ]
    header += ["fc","sig","description"]

    def getKey(d):
        p = d["pval"]
        f = d["l2fold"]
        if p<0.05: pchar="*" 
        else: pchar = "~"
        bar = math.log(1.4,2)
        if f > bar: fchar="+"
        elif f < -bar: fchar="-"
        else: fchar="~"
        symbols = {"sig":pchar,"fc":fchar}
        return symbols

    lines, sort = [], []

    min_l2fold = math.log(1.4,2)

    for go in r.keys():

        pvals = [ r[go]["name"][name]["pval"] for name in names ]
        if min(pvals) > 0.05: continue

        l2folds = [ r[go]["name"][name]["l2fold"] for name in names ]

        maxl2fold = max( [ abs(x) for x in l2folds ] )
        if maxl2fold < min_l2fold: continue

        obs = [ r[go]["name"][name]["obs"] for name in names ]
        if min(l2folds) > 0 and max(obs) < 4 : continue

        desc = r[go]["desc"]
        line = [ str(go) ]

        keys={"fc":"","sig":""}
        fold = 0

        for name in names:

            namecols = [ r[go]["name"][name][col] for col in cols ]
            symbols = getKey(r[go]["name"][name])
            keys["fc"] += symbols["fc"]
            keys["sig"] += symbols["sig"]
            fold += r[go]["name"][name]["fold"]
            line += [ "%.3g" % x for x in namecols]

        line += [ keys["fc"], keys["sig"]]
        line += [ str(desc) ]

        lines.append(line)
        sort.append(fold)

    #sort on aggregated fold change
    lines_sort=zip(lines,sort)
    lines_sort.sort(key = lambda s: s[1], reverse=True)
    sorted_lines = [x[0] for x in lines_sort]

    write(outfile, sorted_lines, header)


def multiIntersect(peakset, bedfiles, outfile):

    p = pybedtools.BedTool(peakset)
    all_names = [f[3] for f in p]

    result = pd.DataFrame({"peak":all_names})
    for bedfile in bedfiles:
        b = pybedtools.BedTool(bedfile)
        i = p.intersect(b)
        names = [f[3] for f in i]
        t = [int(x in names) for x in all_names]
        name = bedfile.split("/")[-1][:-len(".bed")]
        result[name] = pd.Series(t)

    result.to_csv(outfile, sep="\t", index = False)



##########################################################################
###### Functions for motif analysis  #####################################
##########################################################################

def count(filename):
    '''FastaIterator.count is generating errors, therefore code for non- .gz files reproduced here'''
    statement = "cat %s | grep -c '>'" % filename

            

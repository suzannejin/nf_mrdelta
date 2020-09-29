#!/usr/bin/env python3

# Python 3.7.2


__description__ = '''

Python version of the Cedric's script for the analysis of homoplasy.

'''



def csv2analyse(df, score, metric, minseq, maxseq, mrdelta, norm, corr, maintain, deltaby, combs):
    '''
    Count homoplasy output
    
    '''

    import collections
    import numpy as np
    import pandas as pd

    # Normalize if required
    if norm:
        if norm == 'lenXnseq':
            df[metric] = df[metric].div(df.nseq, axis = 0).div(df.len, axis = 0)
        elif norm == 'avgLenXnseq':
            df[metric] = df[metric].div(df.nseq, axis = 0).div(df.avgLen, axis = 0)

    # Initialize counts
    namecounts = ["total", "unused", "ndeltaN", "ndeltaP", "ndeltaO", "minScore", "minNumb", "maxScore", "maxNumb"]
    cols = maintain + ["between"] + namecounts
    counts = pd.DataFrame(columns=cols)

    # Group data frame
    index = 0
    g = list(csv.groupby(maintain))  
    for element in g: 
        names = list(element[0])
        #names = dict(zip(keys,element[0]))  # A dictionary with the names of the combinations. Eg { family: seatoxin, bucket: 50, aligner: CLUSTALO }
        comb  = element[1]  # A dataframe with same keys[0:2]. It has as many columns as the original dataframe.

        # Dictionary of counts
        c = dict(zip( maintain + ["between"],  names  + [deltaby] ))
        c["total"]  = comb.shape[0] * (comb.shape[0]-1)  # Total number of deltas, if original data remained untouched
        comb  = filter_unused(comb, minseq, maxseq, score, metric)  
        if comb.shape == 0:
            continue
        delta = compute_delta(comb, score, metric)
        delta = filter_delta(delta, mrdelta)
        c["unused"] = c["total"] - delta.shape[0]  # Number of deltas filtered out
        c["ndeltaN"], c["ndeltaP"], c["ndeltaO"], c["minScore"], c["minNumb"], c["maxScore"], c["maxNumb"] = count_statistics(delta)

        # Merge info to dataframe
        tmp = pd.DataFrame(c, index=[index])
        counts = pd.concat( [counts,tmp] )
        index += 1

    # Reset dtypes
    c = dict(zip(namecounts, ['float64']*len(namecounts)))
    if 'bucket' in maintain:
        c['bucket'] = 'int32'
    counts = counts.astype(c, errors='ignore', copy=False)

    # Set correct correlation: positive or negative
    if corr:
        correlation = corr
    else:
        if counts['ndeltaP'].sum() > counts['ndeltaN'].sum():
            correlation = 'positive'
        else:
            correlation = 'negative'
    if correlation == 'negative':  # Inverse if negative correlation
            tmp = {'minScore':'maxScore', 'maxScore':'minScore', \
                'minNumb' :'maxNumb' , 'maxNumb' :'minNumb'}
            counts = counts.rename(columns=tmp)

    return(counts, correlation)
        

def filter_unused(df, minseq, maxseq, score, metric):
    ''' Remove rows:
            1. with missing values
            2. nseq < minseq or nseq > maxseq
    '''

    import pandas as pd

    # Remove NA
    df = df.dropna(subset=[score, metric])

    # Remove based on the number of sequences
    cond = (df.nseq < minseq) | (df.nseq > maxseq)
    df = df.drop(df[cond].index)
    
    df = df.reset_index(drop=True)
    return(df)


def filter_delta(delta, mrdelta):
    ''' Remove rows from delta if:
            1. both abs(delta) < 0.001
            2. 2*abs(m1-m2)/(m1+m2) < mrdelta threshold given by user
    '''
    
    import pandas as pd

    # Baseline removal
    if mrdelta == -1:
        # Remove if both abs delta == 0
        cond = (abs(delta.delta_score) == 0) & (abs(delta.delta_metric) == 0)
        delta = delta.drop(delta[cond].index)
    else:
        # Remove if both abs delta <= 0.001
        val = 0.001
        cond = (abs(delta.delta_score) <= val) & (abs(delta.delta_metric) <= val)
        delta = delta.drop(delta[cond].index)
        
        # Remove if metrics are 0
        # These are removes so that we won't have the problem of mrdelta = 0/0 = inf
        cond = (delta.metric1 == delta.metric2) & (delta.metric1 == 0)
        delta = delta.drop(delta[cond].index)

    # Remove if mrdelta < threshold
    cond = (delta.mrdelta < mrdelta)
    delta = delta.drop(delta[cond].index)

    delta = delta.reset_index(drop=True)
    return(delta)

                        
def compute_delta(comb, score, metric):
    ''' Compute delta-score and delta-metric.
    
    Input
    -----
    comb   : pandas dataframe. The delta values will be computed between each of its rows. Let's say between row1-row2, row1-row3, [...] 
    score  : string. tc, sp or col
    metric : string. homo, whomo, whomo2, len, ngap, or ngap2

    Output
    ------
    delta  : pandas dataframe. It has six columns: score1, score2, metric1, metric2, delta_score, delta_metric and mrdelta; and as many rows as deltas.
    
    '''

    import numpy as np
    import pandas as pd

    # Initialize
    delta = pd.DataFrame(columns=["score1", "score2", "metric1", "metric2"], dtype='float64')

    vals = { "score" : comb[score].tolist(), "metric": comb[metric].tolist() } # Get values from columns 'score' and 'metric'
    for i in range( len(vals["score"]) ):  # For each row

        # Get score1 and score2
        val_score = vals["score"][i]
        mat_score = vals["score"][:]
        del mat_score[i]

        # Get metric1 and metric2
        val_metric = vals["metric"][i]
        mat_metric = vals["metric"][:]
        del mat_metric[i]

        # Merge dataframe
        tmp = pd.DataFrame(np.array([[val_score] * len(mat_score), \
                                     mat_score, \
                                     [val_metric] * len(mat_metric), \
                                     mat_metric]).transpose(), \
                                     columns=["score1", "score2", "metric1", "metric2"])
        delta = pd.concat( [delta, tmp] )
        
    delta = delta.reset_index(drop=True)

    # Compute deltas
    delta["delta_score"]  = delta.score1  - delta.score2
    delta["delta_metric"] = delta.metric1 - delta.metric2
    
    # Compute mrdelta
    delta['mrdelta'] = 2 * abs(delta.metric1 - delta.metric2) / (delta.metric1 + delta.metric2)
    
    return(delta)


def count_statistics(delta):
    ''' 
    Input
    -----
    delta    : pandas dataframe with the information of the deltas.

    Output
    ------
    ndeltaN  : number of deltas in negative correlation
    ndeltaP  : number of deltas in positive correlation
    ndeltaO  : number of deltas not in negative nor positive correlation
    minScore : sum of score if wrong alignment (according to metric) is always chosen
    minNumb  : counts of wrong alignment considered
    maxScore : sum of score if good alignment (according to metric) is always chosen
    maxNumb  : counts of good alignments considered

    '''

    # Positive correlation between score-metric (the interesting ones)
    # ++, --
    cond = ((delta.delta_score > 0) & (delta.delta_metric > 0)) | ((delta.delta_score < 0) & (delta.delta_metric < 0))
    ndeltaP = delta[cond].shape[0]

    # Negative correlation between score-metric
    # +-, -+
    cond = ((delta.delta_score > 0) & (delta.delta_metric < 0)) | ((delta.delta_score < 0) & (delta.delta_metric > 0))
    ndeltaN = delta[cond].shape[0]

    # Other correlation
    cond = ((delta.delta_score == 0) & (delta.delta_metric != 0)) | ((delta.delta_score != 0) & (delta.delta_metric == 0)) 
    ndeltaO = delta[cond].shape[0]

    # Score if wrong/good alignment (according to metric) is always chosen
    # Supposing positive correlation (otherwise, they will be changed afterwards)
    minScore, maxScore = 0.0, 0.0
    minNumb, maxNumb = 0, 0
    cond1 = delta.metric1 <  delta.metric2
    cond2 = delta.metric1 >= delta.metric2  # So the equal values will be added as wrong and good alignments (50/50)
    df1 = delta[cond1] 
    df2 = delta[cond2]
    # Wrong alignment according to metric
    minScore += df1.score1.sum()
    minScore += df2.score2.sum()
    minNumb += df1.shape[0]
    minNumb += df2.shape[0]
    # Good alignment according to metric
    maxScore += df1.score2.sum()
    maxScore += df2.score1.sum()
    maxNumb += df1.shape[0]
    maxNumb += df2.shape[0]

    return(ndeltaN, ndeltaP, ndeltaO, minScore, minNumb, maxScore, maxNumb)


def calculate_ratio(df, correlation):

    d = {}

    # Calculate
    # Nratio Pratio unusedRatio minacc maxacc deltaacc     
    usedNumber = df["total"] - df["unused"]
    d["unusedRatio"] = df["unused"] / df["total"]
    d["minacc"] = df["minScore"] / df["minNumb"]
    d["maxacc"] = df["maxScore"] / df["maxNumb"]
    d["deltaacc"] = d["maxacc"] - d["minacc"]
    if correlation == 'positive':
        d["Ratio"] = df["ndeltaP"] / usedNumber
    else:
        d["Ratio"] = df["ndeltaN"] / usedNumber

    return(d)


def get_ratio(df, correlation, maintain, deltaby):

    import numpy as np
    import pandas as pd

    original = df.copy()

    # Sum by group
    g = df.groupby(maintain)
    lev = [i for i in range(0, len(maintain) )]
    df = g.sum().reset_index(level=lev)

    # Calculate ratio
    d = calculate_ratio(df, correlation)
    tmp1 = df[maintain]
    tmp2 = pd.DataFrame(d)
    df2 = pd.concat( [tmp1, tmp2], axis=1)

    # Determine nfam
    nfam = []
    for element in g:
        comb = element[1]
        cond = comb.unused < comb.total
        if deltaby != 'family':
            nfam.append( len(comb[cond]["family"].unique()) )
        else:
            nfam.append( np.nan )
    df2["nfam"] = nfam
    return(df2)


def get_global_ratio(df, correlation, deltaby):

    import numpy as np
    import pandas as pd

    # Parse df
    namecounts = ["total", "unused", "ndeltaN", "ndeltaP", "ndeltaO", "minScore", "minNumb", "maxScore", "maxNumb"]
    df2 = df[namecounts].sum()

    # Calculate
    d = calculate_ratio(df2, correlation)
    df2 = pd.DataFrame(d, index=[0])

    # Determina nfam
    cond = df.unused < df.total
    if deltaby != 'family':
        nfam = len(df[cond]["family"].unique())
    else:
        nfam = np.nan
    df2['nfam'] = nfam

    return(df2)


def add_columns(df, add_d, maintain, typeCount):

    n = df.shape[0]

    # Add columns
    d = {}
    for k,v in add_d.items():
        d[k] = [v] * n
    tmp = pd.DataFrame(d)
    df = pd.concat( [tmp, df], axis=1 )

    # Column names
    if typeCount == "counts":
        colcounts = ["total", "unused", "ndeltaN", "ndeltaP", "minScore", "minNumb", "maxScore", "maxNumb"]
        cols = ["mrdelta", "score", "metric"] + maintain + ["between"] + colcounts
    else:
        colcounts = ["Ratio", "unusedRatio", "nfam", "minacc", "maxacc", "deltaacc"]
        cols = ["mrdelta", "score", "metric"] + maintain + ["between"] + ["correlation"] + colcounts

    # Reorder columns
    df = df[cols]

    return(df)


if __name__ == '__main__':

    import argparse
    import collections
    import itertools
    import pandas as pd
    import sys

    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("-csv", type=str, required=True, help=".csv filename")
    app.add_argument("-score", type=str, required=True)#, choices=["tc","sp","col"])
    app.add_argument("-metric", type=str, required=True)#, choices=["homo","whomo","whomo2","len","ngap","ngap2"])
    app.add_argument("-norm", choices=['lenXnseq', 'avgLenXnseq'], default=None, \
                    help="Data normalization. Divide metric by (alignment length or average sequence length) and (number of sequences).")
    app.add_argument("-mrdelta", type=float, default=0.0)
    app.add_argument("-minseq", type=int, default=0)
    app.add_argument("-maxseq", type=int, default=100000)
    app.add_argument("-outdir", type=str, required=True, help="Output directory")
    app.add_argument("-decimal", type=int, default=3)
    app.add_argument("-maintain", type=str, nargs='+', default=["bucket","aligner","family"], \
                    help="Combination to be maintained, given in the order that the counts are calculated.")
    app.add_argument("-deltaby", type=str, default="tree", help="Delta calculated considering this variable")
    app.add_argument("-corr", type=str, choices=['negative','positive'], default=None)
    args = app.parse_args()

    # Read homoplasy.csv 
    csv = pd.read_csv(args.csv, header=0)

    # unused sequences according to family size
    csv = csv[csv.nseq >= args.minseq]
    csv = csv[csv.nseq <= args.maxseq]
    csv = csv.reset_index(drop=True)

    # Combinations: family, bucket, aligner, tree
    combs = {}
    for element in args.maintain:
        combs[element] = csv[element].unique().tolist()
    combs[args.deltaby] = csv[args.deltaby].unique().tolist()

    # Analyze
    # Get counts
    tmp = args.maintain + [args.deltaby]
    counts, correlation = csv2analyse(csv, args.score, args.metric, args.minseq, args.maxseq, args.mrdelta, args.norm, args.corr, \
                         args.maintain, args.deltaby, combs)   # Original combination: same bucket, aligner and/or family, but different tree
                        #order=tmp, bucket=bucket, aligner=aligner, family=family, tree=tree)   # Original combination: same bucket, aligner and/or family, but different tree

    # Summarize counts
    csvs = {}
    csvs["_".join(args.maintain)] = get_ratio(counts, correlation, args.maintain, args.deltaby)
    d = {'base' : {'mrdelta':args.mrdelta, 'score':args.score, 'metric':args.metric} }
    d["_".join(args.maintain)] = {**d['base'], **{'between': args.deltaby, 'correlation':correlation}}
    for i in range(1, len(args.maintain)):
        iter_comb = itertools.combinations(args.maintain, i)  # All possible combinations to stratify analysis
        for j in iter_comb:
            maint = list(j)
            fixed = [x for x in args.maintain if x not in maint ]
            csvs["_".join(maint)] = get_ratio(counts, correlation, maint, args.deltaby)
            tmp = dict(zip(fixed, ['global']*len(fixed) ))
            d["_".join(maint)] = {**d["_".join(args.maintain)], **tmp}
    csvs["global"] = get_global_ratio(counts, correlation, args.deltaby)
    d["global"] = {**d["_".join(args.maintain)], **dict(zip(args.maintain, ['global']*len(args.maintain) ))}

    # Add mrdelta, score, metric and other columns, reorder & print output
    counts = add_columns(counts, d['base'], args.maintain, "counts")
    counts.round(args.decimal).to_csv(args.outdir + "/" + args.score + "_" + args.metric + "_mrdelta" + str(args.mrdelta) + "_counts.tsv", \
                                        sep="\t", index=False, na_rep='NA')
    for element in csvs:
        csvs[element] = add_columns(csvs[element], d[element], args.maintain, None)
        csvs[element].round(args.decimal).to_csv(args.outdir + "/" + \
            args.score + "_" + args.metric + "_mrdelta" + str(args.mrdelta) + "_" + element + ".tsv", 
            sep="\t", index=False, na_rep='NA')

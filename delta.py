##
# CPCantalapiedra 2021

from io import StringIO
import math
import multiprocessing
import sys

import numpy as np
import pandas as pd


from pastml.acr import pastml_pipeline, acr, _validate_input
from pastml.ml import MPPA
from pastml.models.f81_like import F81

##
## ACR functions
##

##
# Compute ACR
def compute_acr_probs(tree, data, columns, threads, data_sep, id_index, verbose = False):
    name_column = None
    root_date = None
    prediction_method = MPPA
    model = F81

    if verbose == True: print("Computing ACR", file=sys.stderr)

    #
    # Validate input

    roots, columns, column2states, name_column, age_label, parameters, rates = \
    _validate_input(tree, columns, name_column, data, data_sep, id_index,                                
                    root_date, copy_only=False, parameters=None, rates=0)

    forced_joint = False
    reoptimise = False
    smoothing = False
    resolve_polytomies = False
    frequency_smoothing = False

    #
    # Compute ACR

    acr_results = acr(forest=roots, columns=columns, column2states=column2states,
                      prediction_method=prediction_method, model=model, column2parameters=parameters,
                      column2rates=rates,
                      force_joint=forced_joint, threads=threads, reoptimise=reoptimise, tau=None if smoothing else 0,
                      resolve_polytomies=resolve_polytomies, frequency_smoothing=frequency_smoothing)

    # marginal probabilities is a <class 'pandas.core.frame.DataFrame'>
    marginal_probs = acr_results[0]['marginal_probabilities']
    # print(marginal_probs, file=sys.stderr)

    return marginal_probs

##
# Computes ACR and remove tips
def get_acr_probs(tree, data_df, columns, threads, data_set, id_index, verbose = False):
    # create file-like object from dataframe
    s_buf = StringIO()
    data_df.to_csv(s_buf)
    s_buf.seek(0)

    # compute ACR
    marginal_probs = compute_acr_probs(tree, s_buf, columns, threads, data_sep, id_index, verbose)

    # remove ACR probabilities of leaves
    marginal_probs = marginal_probs.drop(data_df.index)
    
    return marginal_probs


##
## Entropies functions
##

def compute_entropies(marginal_probs, verbose = False):

    if verbose == True: print("Computing entropies", file=sys.stderr)
    
    k = len(marginal_probs.columns)
    
    for col in marginal_probs:
        marginal_probs.loc[marginal_probs[col]>1/k, col] = marginal_probs.loc[marginal_probs[col]>1/k, col] / (1 - k) - 1 / (1 - k)

    entropies = marginal_probs.sum(axis=1)
    
    np.random.seed(25)
    entropies[entropies == 0] = entropies[entropies == 0] + np.random.rand(1)/10000
    entropies[entropies == 1] = entropies[entropies == 1] - np.random.rand(1)/10000
        
    return entropies

##
## MCMC functions
##

def lpalpha(alpha, beta, x, lambda0):
    return len(x) * (math.lgamma(alpha+beta) - math.lgamma(alpha)) - alpha * (lambda0 - sum(np.log(x)))

def mhalpha(alpha, beta, x, lambda0, se):
    a0 = alpha
    a1 = None
    r = None

    while r is None:
        a1 = math.exp(np.random.normal(math.log(a0), se))
        lpa1 = lpalpha(a1, beta, x, lambda0)
        lpa0 = lpalpha(a0, beta, x, lambda0)
        try:
            r = min(1, math.exp(lpa1 - lpa0))
        except OverflowError:
            if lpa1 >= lpa0:
                r = 1
            else:
                r = float('-inf')

    if np.random.rand(1) < r:
        pass # return a1
    else:
        a1 = a0 # return a0
    
    return a1

def lpbeta(alpha, beta, x, lambda0):
    return len(x) * (math.lgamma(alpha+beta) - math.lgamma(beta)) - beta * (lambda0 - sum(np.log(1-x)))

def mhbeta(alpha, beta, x, lambda0, se):
    b0 = beta
    b1 = None
    r = None

    while r is None:
        b1 = math.exp(np.random.normal(math.log(b0), se))
        lpb1 = lpbeta(alpha, b1, x, lambda0)
        lpb0 = lpbeta(alpha, b0, x, lambda0)
        try:
            r = min(1, math.exp(lpb1 - lpb0))
        except OverflowError:
            if lpb1 >= lpb0:
                r = 1
            else:
                r = float('-inf')

    if np.random.rand(1) < r:
        pass # return b1
    else:
        b1 = b0 # return b0
        
    return b1

def compute_mcmc_parallel(args):
    alpha, beta, x, lambda0, se, sim, thin, burn = args
    return compute_mcmc(alpha, beta, x, lambda0, se, sim, thin, burn)
    
def compute_mcmc(alpha, beta, x, lambda0, se, sim, thin, burn):
    gibbs = []
    # define steps to simulate
    gibbs_range = range(burn, sim+1, thin)
    # print(f"Gibbs sampling: {gibbs_range[0]}, {gibbs_range[1]}, ..., {gibbs_range[-1]}", file=sys.stderr)

    # burn
    for i in range(1, burn):
        alpha = mhalpha(alpha, beta, x, lambda0, se)
        beta = mhbeta(alpha, beta, x, lambda0, se)
        
    # chain
    for i in range(burn, sim+1):
        alpha = mhalpha(alpha, beta, x, lambda0, se)
        beta = mhbeta(alpha, beta, x, lambda0, se)
                
        if i % thin == 0:
            gibbs.append((alpha, beta))
            
    # print(f"Gibbs sampling: {gibbs[0]}, {gibbs[1]}, ..., {gibbs[-1]}", file=sys.stderr)
    
    return gibbs

##
## Compute delta
##

def compute_delta(entropies, threads, verbose = False):
    delta = None
    
    # Based on
    # https://github.com/mrborges23/delta_statistic/blob/master/code.R

    lambda0 = 0.1   #rate parameter of the proposal 
    se      = 0.5   #standard deviation of the proposal
    sim     = 10000 #number of iterations
    thin    = 10    #we kept only each 10th iterate 
    burn    = 100   #100 iterates are burned-in
    mcmc_reps = 2   # number of times MCMC is run

    if verbose == True: print("Computing delta", file=sys.stderr)
    
    threads = threads if threads < 2 else mcmc_reps

    # Single thread
    if threads == 1:
        divs = [b/a for (a, b) in compute_mcmc(np.random.exponential(1), np.random.exponential(1), entropies, lambda0, se, sim, thin, burn)]
        divs.extend([b/a for (a, b) in compute_mcmc(np.random.exponential(1), np.random.exponential(1), entropies, lambda0, se, sim, thin, burn)])

    # Multiprocess
    else:
        
        def mcmc_args(reps):
            for rep in range(1, reps+1):
                yield (np.random.exponential(1), np.random.exponential(1), entropies, lambda0, se, sim, thin, burn)

        divs = None
        with multiprocessing.Pool(processes = threads) as pool:
            vals = pool.imap_unordered(compute_mcmc_parallel, mcmc_args(mcmc_reps))
            divs = [b/a for sublist in vals for (a, b) in sublist]

    # Compute delta
    if divs is None or len(divs) <= 0:
        raise Exception("Could not compute MCMC alpha-beta values.")
    else:
        # mean of all beta/alpha divisions
        delta = sum(divs)/len(divs)

    return delta


##
## p-value functions
##

##
def compute_random_delta_parallel(args):
    tree, df, columns, threads, data_sep, id_index, verbose = args
    return compute_random_delta(tree, df, columns, threads, data_sep, id_index, verbose)

##
def compute_random_delta(tree, df, columns, threads, data_sep, id_index, verbose):
    # random shuffle of input data
    shuffled_df = df.sample(frac = 1, random_state = np.random.RandomState())
    shuffled_df.index = df.index
    
    # compute ACR
    marginal_probs = get_acr_probs(tree, shuffled_df, columns, threads, data_sep, id_index, verbose)
    
    # Entropies
    entropies = compute_entropies(marginal_probs, verbose)
    
    # Delta
    pval_delta = compute_delta(entropies, threads, verbose)

    return pval_delta

##
def compute_pvalue(delta, pval_reps, tree, df, columns, threads, data_sep, id_index, verbose):
    p_value = None
    
    def compute_random_delta_args(reps):
        threads = 1
        for rep in range(1, reps+1):
            yield (tree, df, columns, threads, data_sep, id_index, verbose)
            
    # obtain random deltas
    random_deltas = []
    with multiprocessing.Pool(processes = threads) as pool:
        random_deltas = pool.imap_unordered(compute_random_delta_parallel, compute_random_delta_args(pval_reps))
        random_deltas = [rd for rd in random_deltas] # parse iterator
        # p-value
        p_value = sum([1 for rd in random_deltas if rd > delta]) / len(random_deltas)
        
    return p_value


##
## API
##

def compute_delta_pvalue(tree, data, columns, threads, pval_reps, data_sep, id_index, verbose = False):
    delta = p_value = None
    
    ##
    ## Parse data to File-like object
    ## (needed to shuffle the data and be able to use it
    ## for pastml ACR function which uses pandas read_csv function)
    if verbose == True: print("Parsing input data", file=sys.stderr)
    
    df = pd.read_csv(data, sep=data_sep, index_col=id_index, header=0, dtype=str)

    ##
    ## Ancestral Character Reconstruction
    
    marginal_probs = get_acr_probs(tree, df, columns, threads, data_sep, id_index, verbose)

    ##
    ## Compute entropies
    # Based on
    # https://github.com/mrborges23/delta_statistic/blob/master/code.R
    
    entropies = compute_entropies(marginal_probs, verbose)

    ##
    ## Compute delta
    ## Based on
    # https://github.com/mrborges23/delta_statistic/blob/master/code.R
    delta = compute_delta(entropies, threads, verbose)

    if verbose == True: print(f"delta statistic: {delta}", file=sys.stderr)

    ##
    ## Compute p-value
    
    if pval_reps > 0:
        
        if verbose == True: print(f"Computing p-value: {pval_reps} repetitions", file=sys.stderr)
        
        p_value = compute_pvalue(delta, pval_reps, tree, df, columns, threads, data_sep, id_index, verbose)
        
        if verbose == True: print(f"p-value: {p_value}", file=sys.stderr)
        
    return (delta, p_value)

##
## MAIN
##

if '__main__' == __name__:
    
    ##
    ## Arguments
    
    tree = sys.argv[1]
    data = sys.argv[2]
    columns = sys.argv[3].split(",")
    if len(sys.argv) > 4:
        threads = int(sys.argv[4])
    else:
        threads = 0

    if len(sys.argv) > 5:
        pval_reps = int(sys.argv[5])
    else:
        pval_reps = 0

    if len(sys.argv) > 6:
        data_sep = sys.argv[6]
    else:
        data_sep = ","

    if len(sys.argv) > 7:
        id_index = int(sys.argv[7])
    else:
        id_index = 0
        
    threads = threads if threads > 0 else multiprocessing.cpu_count()
    
    delta, p_value = compute_delta_pvalue(tree, data, columns, threads, pval_reps, data_sep, id_index, verbose = True)
    
    print(f"delta = {delta}")
    if p_value is not None:
        print(f"p-value = {p_value}")
    
    print("Finished.", file=sys.stderr)

## END

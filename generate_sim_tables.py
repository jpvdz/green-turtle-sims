#!/usr/bin/env python

################################ Import modules ################################
import argparse
import msprime
import allel
import numpy as np
import pandas as pd

############################### Define functions ###############################
def read_args():

    '''
    This function parses command line arguments and stores the parameter values
    in respective variables. These are used as input parameters for the
    simulation.
    '''

    # define parser
    parser = argparse.ArgumentParser()

    # add arguments to parser
    parser.add_argument('-i', action = 'store', dest='job_id', type = int,
                        help = 'SLURM array job ID.')
    parser.add_argument('-r', action = 'store', dest='num_sims', type = int,
                        help = 'Number of simulations per parameter set.')
    parser.add_argument('-n', action = 'store', dest='diploid_ne', type = int,
                        help = 'Diploid effective population size.')
    parser.add_argument('-m', action = 'store', dest='migr_rate', type = float,
                        help = 'Migration rate.')
    parser.add_argument('-s', action = 'store', dest='sample_size', type = int,
                        help = 'Number of samples per population.')
    parser.add_argument('-u', action = 'store', dest='mut_rate', type = float,
                        help = 'Mutation rate.')
    parser.add_argument('-b', action = 'store', dest='base_pairs', type = int,
                        help = 'Size of DNA fragment.')
    parser.add_argument('-L', action = 'store', dest='num_loci', type = int,
                        help = 'Number of loci to simulate.')
    parser.add_argument('-o', action='store', dest='output_dir', type = str,
                        help = 'Name of output folder.')

    # parse arguments
    args = parser.parse_args()

    # define new variables containing stored arguments
    job_id = args.job_id
    num_sims = args.num_sims
    diploid_ne = args.diploid_ne
    migr_rate = args.migr_rate
    sample_size = args.sample_size
    mut_rate = args.mut_rate
    base_pairs = args.base_pairs
    num_loci = args.num_loci
    output_dir = args.output_dir

    # return variables
    return job_id, num_sims, diploid_ne, migr_rate, sample_size, mut_rate, \
        base_pairs, num_loci, output_dir

def simulate_rads_2pop(ne, m, n_samples, mu, bp, n_loci):

    '''
    This function simulates data under a two-population model with symmetrical
    migration. The diploid effective population size (ne), migration rate (m),
    diploid sample size (n), mutation rate (m), locus length (bp) and the number
    of loci (n_loci) need to be specified. For each simulated dataset, the
    number of SNPs, Hudson's Fst (count-based) and Weir & Cockerham's theta
    (frequency-based) are calculated and returned.
    '''

    # convert diploid sample size to number of chromosomes
    n_chrom = n_samples * 2

    # configure populations
    pop_configs = [
        msprime.PopulationConfiguration(sample_size=n_chrom, initial_size=ne/2),
        msprime.PopulationConfiguration(sample_size=n_chrom, initial_size=ne/2)
    ]

    # define symmetrical migration matrix
    migr_mat = [
        [0, m],
        [m, 0],
    ]

    # run simulation
    reps = msprime.simulate(population_configurations=pop_configs,
                               migration_matrix=migr_mat,
                               mutation_rate=mu,
                               length=bp,
                               num_replicates=n_loci)

    # create empty matrix to hold genotypes [replicates x chromosomes * 2]
    gen_mat = np.zeros((n_loci, n_chrom * 2), dtype=int)

    n_snps = 0 # count the number of SNPs
    for i, ts in enumerate(reps):
        for variant in ts.variants():
            if (variant.site.id == 0): # keep only the first snp per locus
                gen_mat[n_snps] = variant.genotypes
                n_snps += 1

    # convert matrix -> HaplotypeArray -> GenotypeArray
    h = allel.HaplotypeArray(gen_mat[:n_snps])
    g = h.to_genotypes(ploidy=2) # rows: variants; columns: samples

    # define subpopulations
    subpops = [list(range(0, n_samples)),
              list(range(n_samples, n_samples * 2))]

    # estimate hudson's fst based upon allele counts
    ac1 = g.count_alleles(subpop=subpops[0])
    ac2 = g.count_alleles(subpop=subpops[1])
    num, den = allel.hudson_fst(ac1, ac2)
    hudson_fst = np.sum(num) / np.sum(den)

    # estimate weir and cockerham's theta
    a, b, c = allel.weir_cockerham_fst(g, subpops)
    wc_fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
    wc_fst

    return [n_snps, hudson_fst, wc_fst]

def write_log(output_dir):
    '''
    This function writes the versions of the modules used in the script to
    a .log file.
    '''

    # define log
    log = 'Simulation {0} was run using: {1} ({2}), {3} ({4}), {5} ({6}), {7} ({8})'.format(output_dir,
               msprime.__name__, msprime.__version__,
               allel.__name__, allel.__version__,
               np.__name__, np.__version__,
               pd.__name__, pd.__version__,)

    # define logfile name
    logfile = '{0}/{0}.log'.format(output_dir)

    # write log to file
    with open(logfile, 'w') as f:
        f.write(log)

def main():

    # read command line arguments
    job_id, num_sims, diploid_ne, migr_rate, sample_size, mut_rate, \
        base_pairs, num_loci, output_dir = read_args()

    # define empty lists to hold results
    n_smpl = []
    n_snps = []
    hudson_fst = []
    wc_fst = []

    # run simulations; append results to lists
    for i in range(num_sims):
        snp, huds, wc = simulate_rads_2pop(diploid_ne, migr_rate, sample_size,
            mut_rate, base_pairs, num_loci)
        n_smpl.append(sample_size)
        n_snps.append(snp)
        hudson_fst.append(huds)
        wc_fst.append(wc)

    # convert simulation results to dict
    results = {
        'sample_size' : n_smpl,
        'n_snps' : n_snps,
        'hudson_fst' : hudson_fst,
        'wc_fst' : wc_fst
    }

    # create pandas df from dict and write to csv
    df = pd.DataFrame(results)
    job_num = str(job_id).rjust(2, '0')
    outfile = f"{output_dir}/simulation_table_n{job_num}.csv"
    df.to_csv(outfile, index = False)

    # write log file for first job only
    if job_id == 1:
        write_log(output_dir)

################################# Execute main #################################
if __name__ == '__main__':
    main()

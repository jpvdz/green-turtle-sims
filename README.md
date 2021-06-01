# Assessing minimum sample sizes *in silico*
A collection of scripts used for assessing how many samples are required to reliably estimate genetic differentiation given a certain effective population size and migration rate. The analysis revolves around simulating RAD-like data under a simple two-population island model for a range of sample sizes using `msprime`. For each sample size, genetic differentiation between the two simulated populations is estimated using `scikit-allel`. 

## Dependencies
- Python3 is required to run `generate_sim_tables.py`, along with the `msprime`, `scikit-allel`, `numpy` and `pandas` Python packages.
- `run_sims.sh` is intended to run on a High Performance Computing cluster with SLURM.  
- `summarize_sims.R` requires the `tidyverse`, `Cairo` and `gridExtra` R packages. 

## Instructions
The `run_sims.sh` script is a wrapper for `generate_sim_tables.py` that reads a tab-delimited input file (`sim_r1000.txt`) containing different parameter sets on each line and starts a SLURM array job for each parameter set. It can be executed using:

`sbatch --array=1-15 run_sims.sh -i sim_r1000.txt`

The first column in the input file is the *number of simulations* to perform for a set of parameters, the second the *diploid effective population size*, the third the *migration rate*, the fourth the *sample size*, the fifth the *mutation rate*, the sixth the *size of the DNA fragment* and the seventh column is the *number of loci* to simulate.   

The core of the analysis is `generate_sim_tables.py`, which performs the simulations and estimations. It can be run without the wrapper like this:

`./generate_sim_tables.py -i 1 -r 1000 -n 10000 -m 0.00025 -s 2 -u 1e-8 -b 200 -L 20000 -o sim_r1000`

The `-i` flag sets the job ID and appends this number to the output file in the output directory (here it would be: `sim_r1000/simulation_table_n01.csv`).

To combine different output files, you can simply use:

`cat sim_r1000/simulation_table_n*.csv > simulation_table.csv`

The resulting file can be analyzed using `summarize_sims.sh`, which summarizes the simulation results and generates plots.  

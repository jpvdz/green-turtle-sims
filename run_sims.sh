#!/usr/bin/bash

#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=1gb
#SBATCH --job-name=run_sims
#SBATCH --array=1-1

# load modules
module load Python/3.8.2-GCCcore-9.3.0
module load GSL/2.6-GCC-9.3.0

# read command line arguments
while getopts 'i:' OPTION; do
    case $OPTION in
        i) input_file="$OPTARG"
            ;;
        ?) printf "Usage: %s: [-i] [-o]" ${0##*/} >&2
            exit 2
            ;;
    esac
done
shift $(($OPTIND - 1))

# create output directory using input file name as dirname
output_dir=${input_file%.*}
mkdir -p ${output_dir}

# define empty arrays to hold simulation parameters
num_sims=()
diploid_ne=()
migr_rate=()
num_samples=()
mut_rate=()
base_pairs=()
num_loci=()

# read simulation parameter values
while read -r col1 col2 col3 col4 col5 col6 col7
do
    num_sims+=($col1)
    diploid_ne+=($col2)
    migr_rate+=($col3)
    num_samples+=($col4)
    mut_rate+=($col5)
    base_pairs+=($col6)
    num_loci+=($col7)
done < ${input_file}

# get array job id
job_id=${SLURM_ARRAY_TASK_ID}

# run simulation
python generate_sim_tables.py \
    -i ${job_id} \
    -r ${num_sims[${job_id}-1]} \
    -n ${diploid_ne[${job_id}-1]} \
    -m ${migr_rate[${job_id}-1]} \
    -s ${num_samples[${job_id}-1]} \
    -u ${mut_rate[${job_id}-1]} \
    -b ${base_pairs[${job_id}-1]} \
    -L ${num_loci[${job_id}-1]} \
    -o ${output_dir}

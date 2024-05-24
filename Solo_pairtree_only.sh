#!/usr/bin/bash
#SBATCH --mem=90G
#SBATCH --job-name=Solo_pairtree_RunTime
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=15

# Activate the pairtree environment
source activate pairtree

# Set the concentration parameter from the command line argument
concentration=$2

# Run pairtree cluster-variants
python3 ~/pairtree/bin/clustervars --concentration $concentration --chains 1 --iterations 1000 --parallel 100 $1".complete.ssm" $1".complete.json" $1".filenameout.json"

# Run PairTree
python3 ~/pairtree/bin/pairtree --params $1".filenameout.json" $1".complete.ssm" $1".results.npz"

# Run plotTree
python3 ~/pairtree/bin/plottree --runid $1 $1".comple

#!/bin/bash
#
# Re-run the code locally, to re-create the data and figure.
#
# Usage:
#
#   ./scripts/rerun.sh
#
#SBATCH --partition=gelifes
#SBATCH --time=32:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --job-name=pirex20
#SBATCH --output=example_20.log
#
rm -rf example_20
rm *.png
time Rscript example_20.R
zip -r pirouette_example_20.zip example_20 example_20.R scripts *.png


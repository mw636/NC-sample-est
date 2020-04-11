Copyright:
Marcel Weiß and Sebastian E. Ahnert, 2020

For usage, please cite:
Marcel Weiß and Sebastian E. Ahnert,
"Using small samples to estimate neutral component size and robustness in the genotype-phenotype map of RNA secondary structure",
J. R. Soc. Interface, 2020



########################################

This standalone program returns estimates of neutral component size, extrapolated neutral network size and neutral component robustness, as well as the number of RNA.fold() calls, starting from a given input RNA sequence, sample size, random subsample size and random seed. For the sampling process, accelerated site scanning sampling with random subsampling is used. For the neutral component size estimation, we calculate the correction parameter from the functional relation and parameters derived in the article.

In case of questions, please contact: mw636@cam.ac.uk (Marcel Weiß).



# Computational requirements ########################################

The code is written in Python (assuming version 2.7).

The following Python packages / modules are required:
- random (Python standard library)
- sys (Python standard library)
- os (Python standard library)
- numpy
- tabulate
- viennarna, Python implementation of the ViennaRNA package: https://www.tbi.univie.ac.at/RNA/



# Usage ########################################

The program requires four command line arguments:

1. Input RNA sequence:
This can be a letter string of any length. An error is returned if the sequence includes non-standard nucleotides (i.e. letters not in the RNA alphabet {A,C,G,U}), or if the sequence leads to the undefined phenotype (unbound structure) by calling RNA.fold().

2. Sample size:
This has to be a non-zero positive integer.

3. Random subsample size:
This has to be a non-zero positive integer. An error is returned if the random subsample size is larger than the sample size.

4. Random seed:
This has to be an integer, which is used to initialise the random number generator.

If not existing, the program creates a directory with the name "results". In addition, if not existing, it creates a subdirectory with the name "XX1", where XX1 is the input sequence. The results are written into a .txt file with the name "XX1_estimation_run.XX2.txt", where XX1 is the input sequence and XX2 the number of the estimation run considered for this input sequence. Starting from 1, for every estimation run for this input sequence (no matter if different or the same sample size, random subsample size, or random seed are considered), the number of the estimation run is increased by one.

The .txt file includes the input sequence, the predicted structure, the sample size, the random subsample size, the random seed, the neutral component size estimate, the extrapolated neutral network size estimate, the neutral component robustness estimate and the number of RNA.fold() calls. In the latter case, it is the sum of RNA.fold() calls for the sampling process and the measurement of the one-point mutational neighbourhoods of the random subsample genotypes.



# Example ########################################

This directory includes an example of how to use the program. It can be started by running the script file "example.sh", which calls:

python NC_estimation.py CCAUGGUGGUGGCUGGGGUCAGCCCCACGGUGGUGGCUGG 1000 10 1

i.e. an estimation run for input sequence "CCAUGGUGGUGGCUGGGGUCAGCCCCACGGUGGUGGCUGG" considering a sample size of 1000, random subsample size of 10 and a random seed of 1.

'''
Copyright:
Marcel Weiß and Sebastian E. Ahnert, 2020

For usage, please cite:
Marcel Weiß and Sebastian E. Ahnert,
"Using small samples to estimate neutral component size and robustness in the genotype-phenotype map of RNA secondary structure",
J. R. Soc. Interface, 2020
'''

import random
import sys
import os
import numpy as np
from tabulate import tabulate

import RNA



# parameters ########################################

# command line arguments
input_sequence = sys.argv[1]
sample_size = int(sys.argv[2])
random_subsample_size = int(sys.argv[3])
random_seed = int(sys.argv[4])
random.seed(random_seed)

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet
alphabet_compatible_bp_partners = [["U"],["G"],["C","U"],["A","G"]] # list of (lists of) letters that in principle are able (/ compatible) to form a bp with letter at respective position in RNA alphabet
# alpha function parameters (for NC size estimation correction)
fit_A = 0.68
fit_B = 0.079



# checks ########################################

# sequence length
L = len(input_sequence)

# check if input sequence only includes standard nucleotides ({A,C,G,U})
for letter in input_sequence:
    if letter not in alphabet:
        print "error: input sequence includes non-standard nucleotides"
        quit()

# check if input sequence leads to defined phenotype
(structure, mfe) = RNA.fold(input_sequence)
if structure=="."*L:
    print "error: input sequence leads to undefined phenotype (unbound structure)"
    quit()
ref_structure = structure

# check if random subsample size smaller or equal than sample size
if random_subsample_size>sample_size:
    print "error: random subsample size > sample size"
    quit()



# directory and file structures ########################################

# results directory name
results_directory_name = "results"
results_directory_name += "/%s" % (input_sequence)
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# help functions ########################################

# help function: find structure
def find_structure(sequence):
    (structure, mfe) = RNA.fold(sequence)
    return structure

# help function: find bp partner positions
def find_bp_partner_positions(structure):
    bp_partner_positions = ["x" for site in structure] # list that has at each position the position index of the bp partner, "x" if position is not paired
    for i,symbol_1 in enumerate(structure):
        if symbol_1=='(':
            open_counter = 1
            close_counter = 0
            for j,symbol_2 in enumerate(structure):
                if j>i and symbol_2=='(':
                    open_counter += 1
                if j>i and symbol_2==')':
                    close_counter += 1
                if close_counter==open_counter:
                    bp_partner_positions[i]=j
                    bp_partner_positions[j]=i
                    break
    return bp_partner_positions

# help function: 'accelerated site scanning' neutral mutation for a sequence
def site_scanning_neutral_mutation(sequence,ref_mut_position,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners):
    RNAfold_counter = 0
    success = False
    while success==False:
        # position
        mut_position = ref_mut_position
        # mutation alphabet
        mut_alphabet = [letter for letter in alphabet if letter!=sequence[mut_position]]
        # randomly test mutation alphabet until success
        while len(mut_alphabet)>0:
            # random letter
            mut_letter = random.choice(mut_alphabet)
            # check if mutation affects unpaired site or if it affects paired site if it would lead to compatible bp
            if ref_structure[mut_position]=="." or (ref_structure[mut_position]!="." and mut_letter in alphabet_compatible_bp_partners[alphabet.index(sequence[bp_partner_positions[mut_position]])]):
                mut_sequence = [j for j in sequence] # to preserve "sequence"
                # mutate
                mut_sequence[mut_position] = mut_letter
                # check if neutral
                mut_structure = find_structure("".join(mut_sequence))
                RNAfold_counter += 1
                if mut_structure==ref_structure:
                    success=True
                    break
            # if no success, update mutation alphabet
            mut_alphabet.remove(mut_letter)
        # if no success, go to next position
        if len(mut_alphabet)==0 and success==False:
            ref_mut_position = (ref_mut_position+1) % L
    return mut_sequence, mut_position, RNAfold_counter

# help function: find neutral mutations per site for a sequence
def find_neutral_mut_per_site(sequence,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners):
    seq_neutral_mut_per_site = [0]*L
    RNAfold_counter = 0
    # go through all characters of the sequence
    for i,character in enumerate(sequence):
        # go through all letters of the alphabet
        for letter in alphabet:
            # mutate character to "value" of letter if character != letter
            if character!=letter:
                # check if mutation affects unpaired site or if it affects paired site if it would lead to compatible bp
                if ref_structure[i]=="." or (ref_structure[i]!="." and letter in alphabet_compatible_bp_partners[alphabet.index(sequence[bp_partner_positions[i]])]):
                    mut_sequence = [j for j in sequence] # to preserve "sequence"
                    # mutate
                    mut_sequence[i] = letter
                    # check if neutral
                    mut_structure = find_structure("".join(mut_sequence))
                    RNAfold_counter += 1
                    if mut_structure==ref_structure:
                        # if neutral mutation, add it to the specific site
                        seq_neutral_mut_per_site[i] += 1
    return seq_neutral_mut_per_site, RNAfold_counter

# help function: sample measurement (averaging)
def sample_measurement(sample_seq_neutral_mut_per_site_storage,L):

    sample_avg_neutral_mut_per_site = [0 for l in range(L)]
    sample_SD_neutral_mut_per_site = [0 for l in range(L)]

    for k in range(len(sample_seq_neutral_mut_per_site_storage)):
        for l in range(L):
            sample_avg_neutral_mut_per_site[l] += sample_seq_neutral_mut_per_site_storage[k][l]
    for l in range(L):
        sample_avg_neutral_mut_per_site[l] = float(sample_avg_neutral_mut_per_site[l])/len(sample_seq_neutral_mut_per_site_storage)

    for k in range(len(sample_seq_neutral_mut_per_site_storage)):
        for l in range(L):
            sample_SD_neutral_mut_per_site[l] += (sample_seq_neutral_mut_per_site_storage[k][l]-sample_avg_neutral_mut_per_site[l])**2
    for l in range(L):
        if len(sample_seq_neutral_mut_per_site_storage)>1:
            sample_SD_neutral_mut_per_site[l] = np.sqrt(float(sample_SD_neutral_mut_per_site[l])/(len(sample_seq_neutral_mut_per_site_storage)-1))
        if len(sample_seq_neutral_mut_per_site_storage)==1:
            sample_SD_neutral_mut_per_site[l] = 0

    return sample_avg_neutral_mut_per_site, sample_SD_neutral_mut_per_site

# help function: NC size estimation and extrapolated NN size estimation
def NC_size_and_NN_size_estimation(ref_structure,sample_avg_neutral_mut_per_site,sample_SD_neutral_mut_per_site,L,alpha):
    NC_size_est = 1.0
    for l in range(L):
        if ref_structure[l]==".":
            if (1+sample_avg_neutral_mut_per_site[l]+alpha*sample_SD_neutral_mut_per_site[l])>4:
                NC_size_est = NC_size_est*4.0
            else:
                NC_size_est = NC_size_est*(1+sample_avg_neutral_mut_per_site[l]+alpha*sample_SD_neutral_mut_per_site[l])
        if ref_structure[l]!=".":
            if (1+sample_avg_neutral_mut_per_site[l]+alpha*sample_SD_neutral_mut_per_site[l])>2:
                NC_size_est = NC_size_est*2.0
            else:
                NC_size_est = NC_size_est*(1+sample_avg_neutral_mut_per_site[l]+alpha*sample_SD_neutral_mut_per_site[l])
    NN_size_est = NC_size_est*2**(ref_structure.count("("))
    return NC_size_est, NN_size_est

# help function: NC robustness estimation
def NC_robustness_estimation(ref_structure,sample_avg_neutral_mut_per_site,alphabet,L):
    NC_rob_est = float(sum(sample_avg_neutral_mut_per_site))/L/(len(alphabet)-1)
    return NC_rob_est



# program ########################################

RNAfold_counter_sampling = 1 # initial structure prediction
RNAfold_counter_nb_measurement = 0

bp_partner_positions = find_bp_partner_positions(ref_structure)
    
# sampling (accelerated site scanning) from input sequence
ref_sequence = input_sequence
sample = [ref_sequence]
ref_mut_position = 0
while len(sample) < sample_size:
    mut_sequence, mut_position, RNAfold_counter = site_scanning_neutral_mutation(ref_sequence,ref_mut_position,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners)
    sample.append(mut_sequence)
    ref_sequence = mut_sequence
    ref_mut_position = (mut_position+1) % L
    RNAfold_counter_sampling += RNAfold_counter

# random subsample
random_subsample = random.sample(sample,random_subsample_size)
                
# one-point mutational neighbourhood measurement of random subsample sequences
random_subsample_seq_neutral_mut_per_site_storage = []
for sequence in random_subsample:
    seq_neutral_mut_per_site, RNAfold_counter = find_neutral_mut_per_site(sequence,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners)
    random_subsample_seq_neutral_mut_per_site_storage.append(seq_neutral_mut_per_site)
    RNAfold_counter_nb_measurement += RNAfold_counter

# random subsample measurement (averaging)
random_subsample_avg_neutral_mut_per_site, random_subsample_SD_neutral_mut_per_site = sample_measurement(random_subsample_seq_neutral_mut_per_site_storage,L)

# NC and NN size estimation
alpha = fit_A*(1-np.exp(-fit_B*L))
NC_size_est, NN_size_est = NC_size_and_NN_size_estimation(ref_structure,random_subsample_avg_neutral_mut_per_site,random_subsample_SD_neutral_mut_per_site,L,alpha)

# NC robustness estimation
NC_rob_est = NC_robustness_estimation(ref_structure,random_subsample_avg_neutral_mut_per_site,alphabet,L)



# write results to .txt file ########################################

table = [["input sequence:","%s" % (input_sequence)],
         ["predicted structure:","%s" % (ref_structure)],
         ["sample size:","%s" % (sample_size)],
         ["random subsample size:","%s" % (random_subsample_size)],
         ["random seed:","%s" % (random_seed)],
         ["neutral component size estimate:","%s" % (NC_size_est)],
         ["neutral network size estimate (extrapolated):","%s" % (NN_size_est)],
         ["neutral component robustness estimate:","%s" % (NC_rob_est)],
         ["# RNA.fold() calls:","%s" % (RNAfold_counter_sampling+RNAfold_counter_nb_measurement)]]
i = 1
# find "possible" file name
while os.path.exists(results_directory_name + "%s_estimation_run.%s.txt" % (input_sequence,i))==True:
    i += 1
file_name = results_directory_name + "%s_estimation_run.%s.txt" % (input_sequence,i)
f = open(file_name,"w")
f.write(tabulate(table))
f.close()



print "done!"

quit()

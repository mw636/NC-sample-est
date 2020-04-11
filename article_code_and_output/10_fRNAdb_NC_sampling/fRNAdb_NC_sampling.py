import random
import sys
import os
import numpy as np
from tabulate import tabulate

import RNA



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
input_list_section_total_number = int(sys.argv[2]) # total number of sections the input list is split
input_list_section_number = int(sys.argv[3]) # section number considered in this run (starting from 1)
sample_size_list = [int(sample_size) for sample_size in sys.argv[4][1:-1].split(",")] # list of considered sample sizes
random_subsample_size_list = [int(sample_size) for sample_size in sys.argv[5][1:-1].split(",")] # list of considered random subsample sizes

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet
alphabet_compatible_bp_partners = [["U"],["G"],["C","U"],["A","G"]] # list of (lists of) letters that in principle are able (/ compatible) to form a bp with letter at respective position in RNA alphabet
random_seed = 1 # seed for random number generator
random.seed(random_seed)



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/6_fRNAdb_input_filtering/results/"

# results directory name
results_directory_name = "results"
results_directory_name += "/L.%s" % (L)
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# load and process input ########################################

# load input
input_file_name = input_directory_name + "L.%s_filtered_fRNAdb_input.txt" % (L)
input_IDs = list(np.loadtxt(input_file_name, usecols=(0,), dtype=str, skiprows=2, unpack=True))
input_sequences = list(np.loadtxt(input_file_name, usecols=(1,), dtype=str, skiprows=2, unpack=True))
input_pred_structures = list(np.loadtxt(input_file_name, usecols=(2,), dtype=str, skiprows=2, unpack=True))

# restrict to considered input list section
input_list_section_length = len(input_IDs) / input_list_section_total_number
input_list_section_min_index = (input_list_section_number-1)*input_list_section_length
input_list_section_max_index = (input_list_section_number)*input_list_section_length
if input_list_section_number==input_list_section_total_number:
    input_list_section_max_index=len(input_IDs) # last section can be longer if len(input_IDs) / input_list_section_total_number not an integer

# considered input
input_IDs = input_IDs[input_list_section_min_index:input_list_section_max_index]
input_sequences = input_sequences[input_list_section_min_index:input_list_section_max_index]
input_pred_structures = input_pred_structures[input_list_section_min_index:input_list_section_max_index]



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



# program ########################################

# results lists
results_input_ID_list = []
results_input_sequence_list = []
results_input_pred_structure_list = []
results_sample_size_list = []
results_random_subsample_size_list = []
results_RNAfold_counter_sampling_list = []
results_RNAfold_counter_nb_measurement_list = []
results_random_subsample_avg_neutral_mut_per_site_list = []
results_random_subsample_SD_neutral_mut_per_site_list = []

# go through input
for index,input_ID in enumerate(input_IDs):

    ref_structure = find_structure(input_sequences[index])
    # check if it matches with previously predicted structure
    if ref_structure!=input_pred_structures[index]:
        print "error: predicted structure not matches with previously predicted structure"
        quit()
    
    bp_partner_positions = find_bp_partner_positions(ref_structure)
    
    # go through sample sizes
    for sample_size in sample_size_list:
        
        RNAfold_counter_sampling = 1 # initial structure prediction
        
        # sampling (accelerated site scanning) from fRNAdb sequence
        ref_sequence = input_sequences[index]
        sample = [ref_sequence]
        ref_mut_position = 0
        while len(sample) < sample_size:
            mut_sequence, mut_position, RNAfold_counter = site_scanning_neutral_mutation(ref_sequence,ref_mut_position,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners)
            sample.append(mut_sequence)
            ref_sequence = mut_sequence
            ref_mut_position = (mut_position+1) % L
            RNAfold_counter_sampling += RNAfold_counter
        
        # go through random subsample sizes
        for random_subsample_size in random_subsample_size_list:
        
            if random_subsample_size<=sample_size:
            
                # random subsample
                random_subsample = random.sample(sample,random_subsample_size)
                
                # one-point mutational neighbourhood measurement of random subsample sequences
                random_subsample_seq_neutral_mut_per_site_storage = []
                RNAfold_counter_nb_measurement = 0
        
                for sequence in random_subsample:
                    seq_neutral_mut_per_site, RNAfold_counter = find_neutral_mut_per_site(sequence,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners)
                    random_subsample_seq_neutral_mut_per_site_storage.append(seq_neutral_mut_per_site)
                    RNAfold_counter_nb_measurement += RNAfold_counter

                # random subsample measurement (averaging)
                random_subsample_avg_neutral_mut_per_site, random_subsample_SD_neutral_mut_per_site = sample_measurement(random_subsample_seq_neutral_mut_per_site_storage,L)

                results_input_ID_list.append(input_IDs[index])
                results_input_sequence_list.append(input_sequences[index])
                results_input_pred_structure_list.append(input_pred_structures[index])
                results_sample_size_list.append(sample_size)
                results_random_subsample_size_list.append(random_subsample_size)
                results_RNAfold_counter_sampling_list.append(RNAfold_counter_sampling)
                results_RNAfold_counter_nb_measurement_list.append(RNAfold_counter_nb_measurement)
                results_random_subsample_avg_neutral_mut_per_site_list.append(random_subsample_avg_neutral_mut_per_site)
                results_random_subsample_SD_neutral_mut_per_site_list.append(random_subsample_SD_neutral_mut_per_site)



# write results to .txt file ########################################

table = [[results_input_ID_list[i],
          results_input_sequence_list[i],
          results_input_pred_structure_list[i],
          results_sample_size_list[i],
          results_random_subsample_size_list[i],
          results_RNAfold_counter_sampling_list[i],
          results_RNAfold_counter_nb_measurement_list[i]] for i,x in enumerate(results_input_ID_list)]
for i,x in enumerate(results_input_ID_list):
    table[i].extend(results_random_subsample_avg_neutral_mut_per_site_list[i])
file_name = results_directory_name + "L.%s_fRNAdb_NC_sampling_sample_average_neutral_mutations_per_site_section.%s.txt" % (L,input_list_section_number)
f = open(file_name,"w")
headers=["ID","sequence","pred. structure","sample size","random subsample size","# RNAfold: sampling","# RNAfold: nb meas."] + [""]*L
f.write(tabulate(table,headers=headers))
f.close()

table = [[results_input_ID_list[i],
          results_input_sequence_list[i],
          results_input_pred_structure_list[i],
          results_sample_size_list[i],
          results_random_subsample_size_list[i],
          results_RNAfold_counter_sampling_list[i],
          results_RNAfold_counter_nb_measurement_list[i]] for i,x in enumerate(results_input_ID_list)]
for i,x in enumerate(results_input_ID_list):
    table[i].extend(results_random_subsample_SD_neutral_mut_per_site_list[i])
file_name = results_directory_name + "L.%s_fRNAdb_NC_sampling_sample_SD_neutral_mutations_per_site_section.%s.txt" % (L,input_list_section_number)
f = open(file_name,"w")
headers=["ID","sequence","pred. structure","sample size","random subsample size","# RNAfold: sampling","# RNAfold: nb meas."] + [""]*L
f.write(tabulate(table,headers=headers))
f.close()



print "done!"

quit()

import random
import sys
import os
import numpy as np
import json
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
desired_NC_rank = int(sys.argv[2]) # desired NC rank
sampling_method_list = [str(sampling_method) for sampling_method in sys.argv[3][1:-1].split(",")] # list of considered sampling methods
sample_size_list = [int(sample_size) for sample_size in sys.argv[4][1:-1].split(",")] # list of considered sample sizes

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet
random_seed = 1 # seed for random number generator
random.seed(random_seed)



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/2_NC_analysis/step_3/results/L.%s/" % (L)

# results directory name
results_directory_name = "results"
results_directory_name += "/L.%s" % (L)
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# load input data ########################################

# NC characteristics
input_file_name_1 = input_directory_name + "L.%s_neutral_component_characteristics.txt" % (L)

NC_ranks = list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=int, skiprows=2, unpack=True))
NC_indices = list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=int, skiprows=2, unpack=True))
NC_phen_indices = list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=int, skiprows=2, unpack=True))
NC_phen_seqs = list(np.loadtxt(input_file_name_1, usecols=(3,), dtype=str, skiprows=2, unpack=True))
NC_sizes = list(np.loadtxt(input_file_name_1, usecols=(4,), dtype=int, skiprows=2, unpack=True))

# dictionary: genotype index -> NC index
input_file_name_2 = input_directory_name + "L.%s_dictionary_genotype_index_to_NC_index.json" % (L)

with open(input_file_name_2, 'r') as handle:
    gNC_map = json.load(handle)

# set up NCs = lists of genotype indices (up to and including maximum NC rank)
NCs = [[] for NC_index in NC_indices]
for gen_index,NC_index in gNC_map.items():
    i = NC_indices.index(NC_index)
    NCs[i].append(int(gen_index))



# help functions ########################################

# help function (base transfer): genotype sequence -> genotype index
def gen_seq_to_gen_index(gen_seq,alphabet,L):
    gen_index = 0
    for site_index in range(L):
       gen_index += alphabet.index(gen_seq[site_index])*len(alphabet)**(L-site_index-1)
    return gen_index

# help function (base transfer): genotype index -> genotype sequence
def gen_index_to_gen_seq(gen_index,alphabet,L):
    gen_seq = []
    for site_index in range(L):
        value = gen_index // len(alphabet)**(L-site_index-1)
        letter = alphabet[value]
        gen_seq.append(letter)
        gen_index -= value*len(alphabet)**(L-site_index-1)
    return ''.join(gen_seq)

# help function: 'RW' neutral mutation for a genotype index
def RW_neutral_mutation(gen_index,ref_NC_index,gNC_map,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    success = False
    while success==False:
        # random position
        mut_position = random.randint(0,L-1)
        # mutation alphabet
        mut_alphabet = [letter for letter in alphabet if letter!=gen_seq_array[mut_position]]
        # random letter
        mut_letter = random.choice(mut_alphabet)
        mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
        # mutate
        mut_gen_seq_array[mut_position] = mut_letter
        mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
        # check if neutral
        try:
            mut_NC_index = gNC_map["%s" % mut_gen_index]
            if mut_NC_index==ref_NC_index:
                success=True
        except KeyError:
            pass
    return mut_gen_index

# help function: 'site scanning' neutral mutation for a genotype index
def site_scanning_neutral_mutation(gen_index,ref_mut_position,ref_NC_index,gNC_map,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    success = False
    while success==False:
        # position
        mut_position = ref_mut_position
        # mutation alphabet
        mut_alphabet = [letter for letter in alphabet if letter!=gen_seq_array[mut_position]]
        # randomly test mutation alphabet until success
        while len(mut_alphabet)>0:
            # random letter
            mut_letter = random.choice(mut_alphabet)
            mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
            # mutate
            mut_gen_seq_array[mut_position] = mut_letter
            mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
            # check if neutral
            try:
                mut_NC_index = gNC_map["%s" % mut_gen_index]
                if mut_NC_index==ref_NC_index:
                    success=True
                    break
            except KeyError:
                pass
            # if no success, update mutation alphabet
            mut_alphabet.remove(mut_letter)
        # if no success, go to next position
        if len(mut_alphabet)==0 and success==False:
            ref_mut_position = (ref_mut_position+1) % L
    return mut_gen_index, mut_position



# program ########################################

# go through all sampling methods
for sampling_method in sampling_method_list:


    # go through all sample sizes
    for sample_size in sample_size_list:
    
        # create results subdirectories (if not exist)
        try:
            os.makedirs(results_directory_name + "%s/sample_size.%s/" % (sampling_method,sample_size))
        except OSError:
            pass
    
    
        # go through all NCs
        for index,NC_rank in enumerate(NC_ranks):
        
        
            # only consider desired NC
            if NC_rank==desired_NC_rank:
            
            
                ref_NC_index = NC_indices[index]
                
                # NC = list of NC genotype indices
                NC = NCs[index]


                # sample creation
                sample = []


                if sampling_method=="random":
                
                    sample = random.sample(NC,sample_size)


                if sampling_method=="RW":
                
                    # randomly choose first genotype index
                    ref_gen_index = random.sample(NC,1)[0]
                    sample = [ref_gen_index]
                
                    # fill sample
                    while len(sample) < sample_size:
                        mut_gen_index = RW_neutral_mutation(ref_gen_index,ref_NC_index,gNC_map,alphabet,L)
                        sample.append(mut_gen_index)
                        ref_gen_index = mut_gen_index


                if sampling_method=="site_scanning":
                
                    # randomly choose first genotype index
                    ref_gen_index = random.sample(NC,1)[0]
                    sample = [ref_gen_index]
                    ref_mut_position = 0

                    # fill sample
                    while len(sample) < sample_size:
                        mut_gen_index, mut_position = site_scanning_neutral_mutation(ref_gen_index,ref_mut_position,ref_NC_index,gNC_map,alphabet,L)
                        sample.append(mut_gen_index)
                        ref_gen_index = mut_gen_index
                        ref_mut_position = (mut_position+1) % L
                            

                # write results to .txt file
                table = [[i,sample[i],gen_index_to_gen_seq(sample[i],alphabet,L)] for i in range(sample_size)]
                file_name = results_directory_name + "%s/sample_size.%s/" % (sampling_method,sample_size)
                file_name += "L.%s_%s_sampling_sample_size.%s_sample_genotypes_NC_rank.%s.txt" % (L,sampling_method,sample_size,NC_rank)
                f = open(file_name,"w")
                headers=["#","genotype index","genotype sequence"]
                f.write(tabulate(table,headers=headers))
                f.close()



print "done!"

quit()

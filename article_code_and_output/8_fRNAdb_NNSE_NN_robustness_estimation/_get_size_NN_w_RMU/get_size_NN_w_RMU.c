/*
Copyright (C) 2008 Thomas Jorg, Olivier Martin, Andreas Wagner

If you use the program "get_size_NN", please cite the following source

Jorg,T, Martin, OC, Wagner, A. Neutral network sizes of biological 
RNA molecules can be computed and are not atypically small.
BMC Bioinformatics 2008 (in press) 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See http://www.gnu.org/licenses/ for the license document.


*/

/*
** get_size_NN.c
** 
** Made by (Jorg_Thomas)
** 
** Started on  Wed Jul 25 13:20:49 2007 Jorg_Thomas
** Last update Mon Feb  4 15:16:22 2008 Jorg_Thomas_080103
*/


/*

This program calculates, via a nested Monte Carlo approach, an estimate of 
the size of the neutral set (network) of a given RNA secondary structure. 
It needs the following input given as command line options to run 
(for each item we give the default value used if the corresponding 
option is not set on the command line). Examples of a command line 
will be given below under EXAMPLES:

(i) The RNA secondary structure in dot parenthesis notation is set with 
the -s option (Its length must not exceed the constant NBASE_MAX which 
is currently set to 1000):
(((((((.(((((...))))).......(((((......)))))...))))))) 

(ii) The number of thermalization Monte Carlo sweeps to be 
applied before the start of measurements is set with the -d option:
2000

(iii) The number of Monte Carlo sweeps per measurement is set with
the -n option:
2000

(iv) The total number of measurements is set with the -m option:
10

(v) The seed for the random number generator is set with the -i option:
1

(vi) The temperature at which the RNA folding is performed is set with 
the -t option:
37.0

Thus, in the above example, a total of 2000*10 = 20000 Monte Carlo sweeps 
in which the neutral set size is estimated are performed.

The compilation of the program needs a functional installation of the
Vienna RNA package that can be downloaded from:

http://www.tbi.univie.ac.at/~ivo/RNA/


The compilation is then done using the following command line:

gcc -I${Includes_ViennaRNA} get_size_NN.c -o \
    get_size_NN -L${Libraries_ViennaRNA} -lm -lRNA

where Includes_ViennaRNA and Libraries_ViennaRNA are the directories
where the include and the library files of the Vienna Package are
installed, respectively. Note that it doesn't make much sense to play
around with the compiler switches to speed up the code as the execution
time is essentially dominated by the RNAfold routine of the Vienna package.

EXAMPLES

The example with the default settings can be run using

./get_size_NN > out.hammerhead &

The result of a sample run is given in the file out.hammerhead_sample and 
should correspond to the output you get in your out.hammerhead file.

The command line without setting any options corresponds to 

 ./get_size_NN -s "(((((((.(((((...))))).......(((((......)))))...)))))))" -d 2000 -i 1 -m 10 -n 2000 -t 37.0

(The order in which the options are given is not important). 

In the following example the neutral set of a secondary structure with 100 base pairs is estimated at 37°C 
using 50 measurement bins:
  
./get_size_NN -s "...((..((.(((((((((..(((.(((........))))))....(((((......)))))..........))))))............)))))..))." -m 50


The output consists of some lines echoing the input data and then of L
lines (where L is the number of nucleotides or bases in the RNA sequence)
where the L initial sequences obtained from inverse_fold are given. 
Then a line indicating the number of pairs of parentheses is given and 
finally the N measurement lines with the results for the neutral set 
size are given (in the secon of the above examples there are 50 measurement bins).

NOTE: In the unlikely case that you want to adapt the code for
      longer secondary structures (i.e. L > 1000) you need to
      adjust the line:

      #define NBASE_MAX 1000

      to the maximal length you want to use before compiling the code.

The code contains the following global constants that can be modified:

NBASE_MAX  is the maximal number of the bases in the sequence, and
           corresponds to the maximal length of the secondary structure in 
           dot-parenthesis representation that can be handled by the code. 
MU_NSWEEP  is the number of mutation trials per shell that is performed
           during one exchange Monte Carlo sweep. (We recommend using 1).
COMPATIBLE defines whether or not we restrict the Monte Carlo procedure
           to the compatible set of the given secondary structure. Using 
           COMPATIBLE TRUE produces more accurate results as the 
           configuration space of the compatible set is noticeably
           smaller than the space of all the 4^L sequences and therefore 
           the sampling is better. (We recommend to use TRUE).
WITH_RMU   defines whether or not we calculate R_mu, the neutral set's mutational 
           robustness (average over the neutral set of R_mu(G) where R_mu(G) 
           is the mutational robustness of the genotype G). R_mu(G) is the 
           fraction of the 3*L point mutants of the sequence G 
           that lie on the same neutral set, i.e., that fold into the
           same secondary structure as G itself. The calculation of R_mu
           is computationally rather expensive and if only a determination 
           of the size of the neutral sets is wanted the calculation 
           of R_mu should be switched off, by setting WITH_RMU to FALSE. 
*/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/inverse.h"
#include "ViennaRNA/part_func.h"

#define TRUE 0
#define FALSE -1

// The global constants that can be modified:
// NBASE_MAX  is the maximal number of the bases in the sequence and
//            corresponds to the maximal length of the secondary structure
//            that can be handled by the compiled code. 
// MU_NSWEEP  is the number of mutation trials per shell that is performed 
//            during one exchange Monte Carlo sweep. (We recommend to use 1).
// COMPATIBLE defines whether or not we restrict the Monte Carlo procedure
//            to the compatible set of the given secondary structure. Using 
//            COMPATIBLE TRUE produces more accurate results as the 
//            configuration space of the compatible set is noticeably
//            smaller than the space of all the 4^L sequences and therefore 
//            the sampling is better. (We recommend to use TRUE).
// WITH_RMU   defines whether or not we calculate the mutational robustness
//            R_mu, i.e., the fraction of the 3*L point mutants of a sequence G 
//            that lie on the same neutral set, i.e., that fold into the
//            same secondary structure as G itself. The calculation of R_mu
//            is computationally rather expensive and if only a determination 
//            of the size of the neutral sets is wanted the calculation 
//            of R_mu should be switched off, by setting WITH_RMU to FALSE. 

#define NBASE_MAX 1000
#define MU_NSWEEP 1
#define COMPATIBLE TRUE
#define WITH_RMU TRUE

// Do not modify this!
#define LSTRING_MAX NBASE_MAX + 1

char    target_phenotype_input[LSTRING_MAX]; 
char*   target_phenotype;
char**  list_actual_sequences;
int*    CURRENT_DIST;
int*    MAX_DISTANCE_OF_SHELL;
int*    partner;
double* BALL_RATIOS;
int     LSTRING, NSHELL, DMAX, NBASE;


void random_genotype(char *random_sequence) 
{
  // Produces a random genotype sequence of length NBASE
  int    i, rnd_int;
  char   rnd_base;
  double rnd_num;
 
  // Loop over the whole sequence and fill in random RNA bases
  for (i = 0; i < NBASE; i++)
    {
      // Get an integer random number between 0 and 3
      // and convert it into an RNA base (0->A,1->U,2->G,3->C)
      rnd_num = drand48();
      rnd_int = floor(rnd_num*4);
      
      switch(rnd_int) {
      case 0:
        rnd_base = 'A'; 
        break;
      case 1:
        rnd_base = 'U'; 
        break;
      case 2:
        rnd_base = 'G'; 
        break;
      case 3:
        rnd_base = 'C'; 
        break;
      }

      random_sequence[i] = rnd_base;

    }
  
  // Make sure that we have the correct string termination
  random_sequence[NBASE] = '\0';
  
}


void get_target_genotype(char *target_sequence, char *target_phenotype)
{
  // This function gets a RNA sequence (genotype) that folds into a given
  // RNA secondary structure (phenotype).

  float success, free_energy;     
  
  // We do not give up the search until we have a correct fold
  give_up = 0;
  
  // We need to allocate memory for the folding routing (which is used by inverse_fold)
  initialize_fold(LSTRING);
  
  do {
  // Generate random starting sequences until we have a good folding.
    random_genotype(target_sequence);
    success = inverse_fold(target_sequence, target_phenotype);
    if (success == 0) 
      {
        free_energy = energy_of_struct(target_sequence, target_phenotype);
        printf( "The initial sequence on the NN: %s and its free energy: %f\n", 
                target_sequence, free_energy );
        fflush(stdout);
      }
  } while (success > 0); // If success = 0 we have a perfect match!

  // Realease the memory used by fold
  free_arrays(); 

}


int check_compatible(char *sequence)
{
  // This routine checks whether a given sequence belongs to the compatible
  // set of the given target phenotype
  int  moving_base;
  char this_base;

  // Loop through the given genotype sequence and for paired
  // bases check whether it's an allowed pairing
  for (moving_base = 0; moving_base < NBASE; moving_base++) 
    {
      if (target_phenotype[moving_base] == '('){
        this_base = sequence[moving_base];
        
        switch(this_base) {
        case 'A':
          if (sequence[partner[moving_base]] != 'U') 
            {
              return FALSE;
            }
          break;
        case 'C':
          if (sequence[partner[moving_base]] != 'G') 
            {
              return FALSE;
            }
         break;
        case 'U':
          if ((sequence[partner[moving_base]] == 'C') || (sequence[partner[moving_base]] == 'U')) 
            {
              return FALSE;
            }
          break;
        case 'G':
          if ((sequence[partner[moving_base]] == 'A') || (sequence[partner[moving_base]] == 'G')) 
            {
              return FALSE;
            }
          break;
        }
      }
    }

  return TRUE;

}


void random_point_mutation(char *initial_sequence, char *modified_sequence)
{
  // Here we perform a random mutation at a random position on a given genotype
  // If COMPATIBLE is TRUE we work only with compatible sequences

  int    position, partner_position, rnd_int, base_num, pair_num;
  char   this_base, partner_base;
  double rnd_num;
  
  // Before we start modifying things we create a copy the initial sequence
  // that later on will be mutated.
  strcpy(modified_sequence, initial_sequence);

  // First we choose the position for the mutation
  rnd_num = drand48();
  position = floor(rnd_num*NBASE);

  // Then we modify the base at the given position.
  // For this we first check what we have at that position 
  // translate it into a number... and modify such that
  // we are sure that we don't get back the same genotype.
  this_base = initial_sequence[position];

  // If we restrict ourselves to the compatible set we need 
  // to distinguish between non-paired and paired bases for
  // the update.
  if (COMPATIBLE == TRUE) 
    {
      if ( target_phenotype[position] == '.' )
        {
          // Convert bases to numbers between 0 and 3
          switch(this_base) {
          case 'A':
            base_num = 0;
            break;
          case 'U':
            base_num = 1; 
            break;
          case 'G':
            base_num = 2; 
            break;
          case 'C':
            base_num = 3; 
            break;
          }
          
          // Get a random shift between 1 and 3 (excluding 0).
          rnd_num  = drand48();
          rnd_int  = floor(rnd_num*3);
          base_num = (base_num + rnd_int + 1) % 4;

          // Convert back to bases
          switch(base_num) {
          case 0:
            modified_sequence[position] = 'A'; 
            break;
          case 1:
            modified_sequence[position] = 'U'; 
            break;
          case 2:
            modified_sequence[position] = 'G'; 
            break;
          case 3:
            modified_sequence[position] = 'C'; 
            break;
          }
        } else {
          // Here are the update rules for the paired bases
          partner_position = partner[position];
          partner_base     = initial_sequence[partner_position];

          // Convert base pairs to numbers between 0 and 5
          if ((this_base == 'A') && (partner_base == 'U')) {
            pair_num = 0;
          } else if ((this_base == 'C') && (partner_base == 'G')) {
            pair_num = 1;
          } else if ((this_base == 'G') && (partner_base == 'C')) {
            pair_num = 2;
          } else if ((this_base == 'G') && (partner_base == 'U')) {
            pair_num = 3;
          } else if ((this_base == 'U') && (partner_base == 'A')) {
            pair_num = 4;
          } else if ((this_base == 'U') && (partner_base == 'G')) {
            pair_num = 5;
          } else {
            printf("Oups! Found an impossible pairing....\n");
            printf("%c %c\n", this_base, partner_base);
            exit(1);
          }
          
          // Get a random shift between 1 and 5 (excluding 0).
          rnd_num  = drand48();
          rnd_int  = floor(rnd_num*5);
          pair_num = (pair_num + rnd_int + 1) % 6;
          
          // And convert back into base pairs.
          switch(pair_num) {
          case 0:
            modified_sequence[position] = 'A'; 
            modified_sequence[partner[position]] = 'U'; 
            break;
          case 1:
            modified_sequence[position] = 'C'; 
            modified_sequence[partner[position]] = 'G'; 
            break;
          case 2:
            modified_sequence[position] = 'G';
            modified_sequence[partner[position]] = 'C';  
            break;
          case 3:
            modified_sequence[position] = 'G';
            modified_sequence[partner[position]] = 'U';  
            break;
          case 4:
            modified_sequence[position] = 'U'; 
            modified_sequence[partner[position]] = 'A'; 
            break;
          case 5:
            modified_sequence[position] = 'U'; 
            modified_sequence[partner[position]] = 'G'; 
            break;
          default: 
            printf("Oups! Something went wrong....\n");
            exit(1);
            break;
          }
        }         
    } else {
      // If we don't restrict ourselves to the compatible set
      // we can simply change one base at random... making again
      // sure that we really change it.
      switch(this_base) {
      case 'A':
        base_num = 0;
        break;
      case 'U':
        base_num = 1; 
        break;
      case 'G':
        base_num = 2; 
        break;
      case 'C':
        base_num = 3; 
        break;
      }
      
      // Get a random shift between 1 and 3 (excluding 0).
      rnd_num  = drand48();
      rnd_int  = floor(rnd_num*3);
      base_num = (base_num + rnd_int + 1) % 4;
      
      // And converting back into bases.
      switch(base_num) {
      case 0:
        modified_sequence[position] = 'A'; 
        break;
      case 1:
        modified_sequence[position] = 'U'; 
        break;
      case 2:
        modified_sequence[position] = 'G'; 
        break;
      case 3:
        modified_sequence[position] = 'C'; 
        break;
      }
      
    }
}


int update_sequence(char* actual_sequence, int maxd_in_shell, int distance_old, int shell_number) 
{
  // Here we try to make a random mutation of a given sequence
  char  modified_sequence[LSTRING], modified_phenotype[LSTRING];
  float energy_fold;
  int   distance_new, distance;
  
  // Preparation for folding the sequence 
  initialize_fold(LSTRING);
  
  // First mutate, then fold the mutated sequence and finally get the 
  // base pair distance between the phenotype of the mutated sequence
  // and the target phenotype.  
  random_point_mutation(actual_sequence, modified_sequence);
  energy_fold  = fold(modified_sequence, modified_phenotype);
  distance_new = bp_distance(modified_phenotype, target_phenotype);
  
  // Free the memory of the fold routine.
  free_arrays();

  // If the mutated phenotype lies within a ball of radius r=maxd_in_shell
  // from the target phenotype we accept the mutation.
  if (distance_new <= maxd_in_shell)
    {
      strcpy(actual_sequence, modified_sequence);
      distance = distance_new;
    }
  // If not we reject the mutation and everything remains as before.
  else
    {
      distance = distance_old;
    }

  return distance;  
}


void sweep_one_shell(int shell_number)
{
  // We do MU_NSWEEP mutations on a given sequence
  // per Monte Carlo sweep  
  int imutation;
  
  for (imutation = 0; imutation < MU_NSWEEP; imutation++) {
    CURRENT_DIST[shell_number] = 
      update_sequence(list_actual_sequences[shell_number], 
                      MAX_DISTANCE_OF_SHELL[shell_number],
                      CURRENT_DIST[shell_number],shell_number);
  }
}


double super_sweep()
{
  // This routine does one complete Monte Carlo sweep consisting
  // of simple Metropolis update for each shell and an exchange 
  // Monte Carlo update. It returns the percentage of success of 
  // the exchange trials.
  int ishell;
  int success = 0;
  int incr_success;

  // The "Metropolis" sweep for each shell
  for (ishell = 0; ishell < NSHELL; ishell++)
    {
      sweep_one_shell(ishell);
    }

  //Here we do the exchange MC  
  for (ishell = 1; ishell < NSHELL; ishell++)
    {
      incr_success     = try_swap(ishell-1,ishell);
      success         += incr_success;
    }
  
  // We return the success rate.
  return success/(double) (NSHELL-1);
}


int try_swap(int ishell1, int ishell2) 
{
  // This defines how the swaps between genotypes in 
  // adjacent shells is done.
  int  dist1, dist2;
  char work_sequence[LSTRING];

  dist1 = CURRENT_DIST[ishell1];
  dist2 = CURRENT_DIST[ishell2];

  // Swap the the two genotypes if both of them lie within a radius smaller
  // or equal the smaller of the two radii of the shells.
  if ((dist1 <= MAX_DISTANCE_OF_SHELL[ishell2]) && (dist2 <= MAX_DISTANCE_OF_SHELL[ishell1]))
    {
      strcpy(work_sequence,list_actual_sequences[ishell1]);
      strcpy(list_actual_sequences[ishell1],list_actual_sequences[ishell2]);
      strcpy(list_actual_sequences[ishell2],work_sequence);
      CURRENT_DIST[ishell1] = dist2;
      CURRENT_DIST[ishell2] = dist1;
    }
  return 0;
}


void initialize(int length)
{
  // Initializes observables, allocates memory, defines the initial
  // sequences and the data structure that allows to keep track of the 
  // compatibility condition. Defines the maximal radius within a given
  // shell.
  int  ishell,ibase,jbase,done;
  char* target_sequence;

  // Define the constants related to the length of the secondary structure
  NBASE = length;
  NSHELL = length;
  LSTRING = length + 1;
  DMAX = length;

  // Get memory
  target_sequence = calloc(LSTRING,sizeof(char));
  target_phenotype = calloc(LSTRING,sizeof(char));
  CURRENT_DIST = calloc(NSHELL,sizeof(int));
  MAX_DISTANCE_OF_SHELL = calloc(NSHELL,sizeof(int));
  partner = calloc(NSHELL,sizeof(int));
  BALL_RATIOS = calloc(NSHELL,sizeof(double));

  list_actual_sequences = calloc(NSHELL,sizeof(char*));
  for (ishell = 0; ishell < NSHELL; ishell++)
    {
      list_actual_sequences[ishell] = calloc(LSTRING,sizeof(char));
    }
  
  // Fill in the target_phenotype 
  strncpy(target_phenotype, target_phenotype_input, LSTRING);

  // Get the initial sequences from an inverse fold (cold start condition)
  // and therefore the initial distances are zero.
  for (ishell = 0; ishell < NSHELL; ishell++)
    {
      get_target_genotype(target_sequence, target_phenotype);
      strcpy(list_actual_sequences[ishell],target_sequence);
      CURRENT_DIST[ishell] = 0;
    }
    
  // Define the maximal radius of the shells. Here define the shells
  // such that the innermost shell has radius 0, i.e., its always containing
  // just sequences from the neutral set. With this definition the first volume
  // ratio is trivially equal to 1. We use this definition since it allows for
  // a very simple and straightfoward way calculate the mutational
  // robustness since all sequences from the innermost shell are element of
  // the neutral set.
  for (ishell = 0; ishell < NSHELL; ishell++)
    {
      MAX_DISTANCE_OF_SHELL[ishell] = ishell;
    }
  
  // Initialize the ratios between two successive shells.
  for (ishell = 0; ishell < NSHELL; ishell++)
    {
      BALL_RATIOS[ishell] = 0.0;
    }

  // Get information on the compatibility condition
  for (ibase = 0; ibase < NBASE; ibase++)
    {
      partner[ibase] = -1;
    }
  
  done  = FALSE;
  ibase = -1;
  jbase = NBASE;

  while (done == FALSE)
    {  
      ibase = get_next_left_paren(ibase);
      if (ibase == NBASE) {
        break;
      } 
      jbase          = matching_paren(ibase);
      partner[ibase] = jbase;
      partner[jbase] = ibase;
    }
}


int get_next_left_paren(int ibase)
{
  // Find the next opening parenthesis.
  int moving_base;

  for (moving_base = ibase + 1; moving_base < NBASE; moving_base++)
    {
      if ( target_phenotype[moving_base] == '(' )
        {
          return moving_base;
        }
    }
  return moving_base;
}


int matching_paren(int ibase)
{
  // Find the corresponding closing parenthesis.
  int height = 1;
  int moving_base;

  for (moving_base = ibase + 1; moving_base < NBASE; moving_base++)
    {
      if ( target_phenotype[moving_base] == '(' )
        {
          height++;
        }
      else if (target_phenotype[moving_base] == ')') 
        {
          height--;
          if (height == 0) 
            {
              return moving_base;
            }
        }
    }
  // Check's that the target phenotype is balanced.
  fprintf(stderr,"There's no matching parenthesis\n");
  exit(1);
}

int count_number_of_paren_pairs()
{
  // Counts the total number of pairs of parenthesis,i.e, the number of
  // paired base pairs.
  int ibase; 
  int nleft_paren  = 0;
  int nright_paren = 0;
  
  for (ibase = 0; ibase < NBASE; ibase++)
    {
      if ( target_phenotype[ibase] == '(' )
        {
          nleft_paren++;
        }
      if ( target_phenotype[ibase] == ')' )
        {
          nright_paren++;
        }
    }
  // Here we do just an additional check that the parenthesis are balanced
  if ( nright_paren != nleft_paren ) 
    {
      fprintf(stderr,"Target phenotype is unbalanced!\n");
      exit(1);
    }
  return nleft_paren;
}


void observables()
{
  // Updates the observables BALL_RATIOS[i] that contains the number of
  // times we find the random walk bounded by a shell with radius R(i) within
  // the subshell with radius R(i-1). 
  int ishell;

  if ( CURRENT_DIST[0] == 0) ++BALL_RATIOS[0];
  for (ishell = 1; ishell < NSHELL; ishell++)
    {
      if ( CURRENT_DIST[ishell] <= MAX_DISTANCE_OF_SHELL[ishell-1] )
        {
          ++BALL_RATIOS[ishell];
        }
    }
}

double mutational_robustness(char* sequence)
{
  // Calculates the mutational robustness of a sequence. This
  // corresponds to evaluate the fraction of the 3L neighboring
  // sequences that fold into the same phenotype as the sequence
  // itself.
  int   position, base_num, modified_base_num, distance, r_mu;
  char  this_base;
  char  neighbor_sequence[LSTRING], neighbor_phenotype[LSTRING];
  float energy_fold;

  r_mu = 0;

  // Loop over all the bases and convert into the base into a number from 0
  // to 3.
  for (position = 0; position < NBASE; position++)
    {
      strcpy(neighbor_sequence, sequence);
      this_base = sequence[position];

      switch(this_base) {
      case 'A':
        base_num = 0;
        break;
      case 'U':
        base_num = 1;
        break;
      case 'G':
        base_num = 2;
        break;
      case 'C':
        base_num = 3;
        break;
      }

      // Do all the 3 mutations per base. 
      for (modified_base_num = 0; modified_base_num < 4; modified_base_num++)
        {

          if (modified_base_num == base_num) continue;

          switch(modified_base_num) {
          case 0:
            neighbor_sequence[position] = 'A';
            break;
          case 1:
            neighbor_sequence[position] = 'U';
            break;
          case 2:
            neighbor_sequence[position] = 'G';
            break;
          case 3:
            neighbor_sequence[position] = 'C';
            break;
          }
          
          // Before we use CPU cycles to fold the sequence we check for compatibility
          if (check_compatible(neighbor_sequence) == TRUE)
            {
              // If compatible we fold the sequence and compare the
              // resulting phenotype using the base pair distance.
              initialize_fold(LSTRING);
              energy_fold = fold(neighbor_sequence, neighbor_phenotype);
              distance    = bp_distance(neighbor_phenotype, target_phenotype);
              if (distance == 0) r_mu++;
              free_arrays();
            }
        }
    }
  // Return the fraction of neighbors having the same phenotype
  return (double) r_mu / (double) (3*NBASE);
}


int main( int argc, char **argv ) 
{
  int    npairs, i, itrial, ishell, length;
  int    NSUPER_SWEEPS, NTHERMAL, NTRIALS;
  long   iseed;
  float  mfrac = 0.0;
  double log_NN_size, r_mu;
  char   k;
  int    c;
  int    dvalue = 2000;
  int    ivalue = 1;
  int    mvalue = 10;
  int    nvalue = 2000;
  char   *svalue = "(((((((.(((((...))))).......(((((......)))))...)))))))\0";
  double tvalue = 37.0;

  opterr = 1;
     
  while ((c = getopt (argc, argv, "d:i:m:n:s:t:")) != -1)
      switch (c)
           {
           case 'd':
             printf("d %s\n",optarg);
             dvalue = atoi(optarg);
             break;
           case 'i':
             ivalue = atoi(optarg);
             break;
           case 'm':
             mvalue = atoi(optarg);
             break;
           case 'n':
             nvalue = atoi(optarg);
             break;
           case 's':
             svalue = optarg;
             break;
           case 't':
             tvalue = atof(optarg);
             break;
           case '?':
             
             exit(1);
           default:
             exit(1);
           }

  temperature = tvalue;
  strcpy(target_phenotype_input,svalue); 
  NSUPER_SWEEPS = nvalue;
  NTRIALS = mvalue;
  iseed = ivalue;
  NTHERMAL = dvalue;

  // Check whether the structure consists just of dots and parenthesis
  // and removes the end-of-line character
  for (i=0; i<LSTRING_MAX; i++){
    k = target_phenotype_input[i];
    if ( k != '(' &&  k != ')' && k != '.' &&  k != '\n' && k != '\0' ) {
      fprintf(stderr,"Junk in the structure!\n");
      exit(1);
    }
    if ( k == '\n') {
      target_phenotype_input[i] = '\0';
    }
  }

  // Get length of the secondary structure 
  length = strlen(target_phenotype_input);
 
  // Print out the parameters used in this run
  printf("The structure under investigation is %s having %i bases\n",target_phenotype_input,length);
  printf("The number of thermalization sweeps is %d\n",NTHERMAL);
  printf("One measurement consists of %d sweeps\n",NSUPER_SWEEPS);
  printf("The total number of measurements is %d\n",NTRIALS);
  printf("The seed of the random number generator is %d\n",iseed);
  printf("The temperature at which the folding is done is %lf\n",temperature);

  // Initialize the random number generator.
  srand48(iseed);
  
  // Further initializations.
  initialize(length);

  // Get the number of pairs of parenthesis in the target phenotype
  npairs = count_number_of_paren_pairs();
  printf("Found %d pairs of parenthesis in the target phenotype.\n", npairs);

  // The "thermalization" loop
  for (i = 0; i < NTHERMAL; i++) 
    {
      super_sweep();
    }

  // Loop over the bins over which we do NSUPER_SWEEPS super sweeps. 
  for (itrial = 0; itrial < NTRIALS; itrial++)
    {
      // Reset the observables
      for (ishell = 0; ishell < NSHELL; ishell++)
        {
          BALL_RATIOS[ishell] = 0.0;
        }

      if ( WITH_RMU == TRUE )
        {
          r_mu = 0.0;
        }

      // Do the measurement loop
      for (i = 0; i < NSUPER_SWEEPS; i++)
        {
          mfrac += super_sweep();
          observables();
          if ( WITH_RMU == TRUE )
          {
            r_mu += mutational_robustness(list_actual_sequences[0]);
          }
        }

      // First we set the volume of the total configuration space
      // (4^L without compatibility condition; 4^(L-2*npairs)*6^npairs
      // with compatibility condition)
      if ( COMPATIBLE == TRUE )
        {
          log_NN_size = (NBASE-2*npairs)*log(4.0) + npairs*log(6);
        }
      else
        {
          log_NN_size = NBASE*log(4.0);
        }
        
      // Determine the estimate of the neutral set size using
      // V(0) = [V(0)/V(1)]*[V(1)/V(2)]*.....*[V(N-1)/V[N]]*V[N].
      // Taking the log we are left with the following:
      for (ishell = 0; ishell < NSHELL; ishell++)
        {
          log_NN_size += log(BALL_RATIOS[ishell]/NSUPER_SWEEPS);
        }
      
      // Here we print out the results.
      if (WITH_RMU == TRUE)
        {
          // Normalize r_mu by the number of measurements
          r_mu /= (double) NSUPER_SWEEPS;
          printf("Estimated neutral set size: %f and its log_10: %f and the mutational robustness: %f\n", 
                 exp(log_NN_size), log_NN_size/log(10), r_mu);
        }
      else 
        {
          printf("Estimated neutral set size: %f and its log_10: %f\n", 
                 exp(log_NN_size), log_NN_size/log(10));
        }

      fflush(stdout); 
    }
  return 0;
}

#!/usr/bin/python3
"""
DESCRIPTION:
    Template code for the FIRST Advanced Question of the Hidden Markov Models
    assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via
    corresponding Canvas assignment. Note this script will be graded manually,
    if and only if your "hmm.py" script succesfully implements Baum-Welch
    training! Continuous Feedback will not be available for this script.

AUTHOR:
    <Theodoros Foskolos id: 2768082>
"""
from argparse import ArgumentParser, RawTextHelpFormatter
from hmm_utility import load_tsv
from numpy.random import choice



def parse_args():
    #####################
    # START CODING HERE #
    #####################
    # Implement a simple argument parser (WITH help documentation!) that parses
    # the information needed by main() from commandline. Take a look at the
    # argparse documentation, the parser in hmm_utility.py or align.py
    # (from the Dynamic Programming exercise) for hints on how to do this.
    
    # Positionals
    parser = ArgumentParser()
    parser.add_argument('transition', help='path to a TSV formatted transition matrix')
    parser.add_argument('emission', help='path to a TSV formatted emission matrix')
    
    # Optionals
    parser.add_argument('-n', dest='n_seq', help='integer number of sequences generated (default: 10)', type=int,
                        default=10)  # N: number of sequences that have to be generated
    parser.add_argument('-o', dest='out_dir', help='path to a directory where output fasta file is written and saved\n'
                                                   '  (directory will be made if it does not exist)\n'
                                                   '  (file names and contents depend on algorithm)',
                        default='random_sequences.fasta')
    return parser.parse_args()
    #####################
    #  END CODING HERE  #
    #####################


def generate_sequence(A,E):
    #####################
    # START CODING HERE #
    #####################
    # Implement a function that generates a random sequence using the choice()
    # function, given a Transition and Emission matrix.
    # Look up its documentation online:
    # https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.random.choice.html


    
    # List to store the random sequences
    seq_list = []
    
    # List to store the random states transitions
    state_list = []
    sequence = str()
    
    # Save temp randomised states
    temp_state = str()
    state = str()
    
    # List of states in order use it for the random sequences by using A
    allstates = list(A.keys())
    
    # List for the emission states
    emissionstates = list(E.keys())
    
    # Random generator from allstates and beginning B state, we wont add B state
    temp_state = choice(a=allstates, p=list(A["B"].values()))
    
    # While loop until i reach end state E, the stop it and remove it i dont need it
    while temp_state != "E":
       # Randomised state from the A list using temp state as a starting point, to get the probability values 
       state_list.append(temp_state)
       state = choice(a=allstates, p=list(A[temp_state].values()))
       temp_state = state
    
    # List with possible value numbers
    possible_seq = list(E[list(emissionstates)[0]].keys())  
    
    for x in state_list:
        seq_list.append(choice(a=possible_seq, p=list(E[x].values())))
        
    # Join the sequence list strings, and map to get the list values as string
    sequence = ''.join(map(str,seq_list))        
    #####################
    #  END CODING HERE  #
    #####################
    
    return sequence


def main():
    args = parse_args()
    #####################
    # START CODING HERE #
    #####################

    # Give a number that you want for generating
    N = args.n_seq 
    # Path to save the sequences (fasta)              
    out_file = args.out_dir
    # A-> transition matrix
    # E-> emission matrix
    A = load_tsv(args.transition)    
    E = load_tsv(args.emission)
    # w for write the sequences in fasta
    with open(out_file,'w') as f:  
        # Generates as many sequences as in parse arguments by using the function generate_sequence for the A and the E
        for i in range(N): 
            seq = generate_sequence(A,E)
            # write it as random sequence i as name \n and the sequence
            f.write('>random_sequence_%i\n%s\n' % (i,seq))
        
    #####################
    #  END CODING HERE  #
    #####################
    


if __name__ == "__main__":
    main()

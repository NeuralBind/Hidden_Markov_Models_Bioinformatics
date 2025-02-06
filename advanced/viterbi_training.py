#!/usr/bin/python3

"""
DESCRIPTION:


INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via
    corresponding Canvas assignment. Note this script is graded automatically,
    if and only if your "hmm.py" script succesfully implements Baum-Welch
    training!

AUTHOR:
    <Theodoros Foskolos >
"""

import os.path as op
from os import makedirs
from argparse import ArgumentParser, RawTextHelpFormatter
from hmm_utility import load_fasta, load_tsv, serialize
from hmm import viterbi



def parse_args():
    "Parses inputs from commandline and returns them as a Namespace object."

    parser = ArgumentParser(prog = 'python3 viterbi_training.py',
        formatter_class = RawTextHelpFormatter, description =
        '  Perform Viterbi training, given a set of sequences with A and E priors.\n\n'
        '  Example syntax:\n'
        '    python3 hmm.py seq.fasta A.tsv E.tsv -i 100 -o /viterbi_outputs'
        '    python3 hmm.py baumwelch in.fa priorA priorE -o ./outputs -i 1')

    # Positionals
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('transition', help='path to a TSV formatted transition matrix')
    parser.add_argument('emission', help='path to a TSV formatted emission matrix')

    # Optionals
    parser.add_argument('-o', dest='out_dir',
        help='path to a directory where output files are saved\n'
             '  (directory will be made if it does not exist)')
    parser.add_argument('-i', dest='max_iter', type=int, default=20,
        help='maximum number of iterations (default: 20 )')

    return parser.parse_args()



def train_viterbi(X,A,E):
  
    # Initialize your posterior matrices
    allStates = A.keys()
    emittingStates = E.keys()
    
    # Initialize a new (posterior) Transition and Emission matrix
    new_A = {}
    for k in A:
        new_A[k] = {l:0 for l in A[k]}
    
    new_E = {}
    for k in E:
        new_E[k] = {s:0 for s in E[k]}
        
    # Get the state path of every sequence in X,
    # using the viterbi() function imported from hmm.py
    # Count the transitions and emissions for every state
    # Every sequence
    for seq in X:
        
        # FOR VITERBI ONLY: Trace back the State Path
        pi, P, V = viterbi(seq, A, E)

        # We need  the first state path from B to the first state of the sequence
        new_A["B"][pi[0]] += 1
        # Enumerate the states of each character of the sequence, 
        # and add +1  on the new_A by checking the current char state and the previous
        for i,s in enumerate(pi):
            if i == 0:    
                pass
            else:
                new_A[pi[i-1]][pi[i]] +=1 
        #For the emission states, add 1 on on matrix new_E with s as a state and seq[i] the current emission on the sequence
            new_E[s][seq[i]] += 1

        # For the emission i dont need the B and the E state everything is 0
        # But for the new_A transitions i need to fill the last one going to "E", end state
        new_A[pi[len(pi)-1]]['E'] += 1

   # Normalize your row sums
   # Normalise transition states rows adding equals 1
    for k in allStates:
        sum_An = 0
        for y in new_A[k]:
            sum_An += new_A[k][y]
        for y in new_A[k]:
            if sum_An != 0:  # we cant divide if its zero
                new_A[k][y] = new_A[k][y] / sum_An #our new Ak divided by the sum of An of the proper sequence state
    
    #Normalise emissions rows adding equals to 1
    for i,e in enumerate(emittingStates):
        sum_En = 0 # make a sum total
        for s in new_E[e]:
            sum_En += new_E[e][s]
        if sum_En != 0:  # we cant divide if its zero
            for s in new_E[e]:
                new_E[e][s] = new_E[e][s] / sum_En 

    return new_A, new_E


def main(args = False):
    "Perform Viterbi training, given a set of sequences with A and E priors."
    
    # Process arguments and load specified files
    if not args: args = parse_args()

    set_X, labels = load_fasta(args.fasta) # List of sequences, list of labels
    A = load_tsv(args.transition) # Nested Q -> Q dictionary
    E = load_tsv(args.emission)   # Nested Q -> S dictionary
    i = 0
    i_max = args.max_iter
    

    while i < i_max:
        # Take the new A and E matrices
        A, E = train_viterbi(set_X, A, E)
        i += 1


    
    if args.out_dir:
        makedirs(args.out_dir, exist_ok=True) # Make sure the output directory exists.
        A_path = op.join(args.out_dir,'viterbi_posterior_A')
        with open(A_path,'w') as f: f.write(serialize(A))
        E_path = op.join(args.out_dir,'viterbi_posterior_E')
        with open(E_path,'w') as f: f.write(serialize(E))        



if __name__ == "__main__":
    main()

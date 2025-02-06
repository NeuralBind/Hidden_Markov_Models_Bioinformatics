#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Hidden Markov Models assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <Theodoros Foskolos id: 2768082>
"""

import os.path as op

from os import makedirs
from math import log10
from hmm_utility import parse_args, load_fasta, load_tsv, print_trellis, print_params, serialize



def viterbi(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the most probable state path, the corresponding P(X), and trellis."""

    allStates = A.keys()
    emittingStates = E.keys()
    L = len(X) + 2

    # Initialize
    V = {k:[0] * L for k in allStates} # The Viterbi trellis
    V['B'][0] = 1.

    # Middle columns
    for i,s in enumerate(X):
        for l in emittingStates:
            terms = [V[k][i] * A[k][l] for k in allStates]
            V[l][i+1] = max(terms) * E[l][s]

    # Last column
    for k in allStates:
        term = V[k][i+1] * A[k]['E'] 
        if term > V['E'][-1]:
            V['E'][-1] = term
            pi = k # Last state of the State Path

    # FOR VITERBI ONLY: Trace back the State Path
    l = pi
    i = L-2
    while i:
        i -= 1
        for k in emittingStates:
            if V[k][i] * A[k][l] * E[l][X[i]] == V[l][i+1]:
                pi = k + pi
                l = k
                break

    P = V['E'][-1] # The Viterbi probability: P(X,pi|A,E)
    return(pi,P,V) # Return the state path, Viterbi probability, and Viterbi trellis



def forward(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the Forward probability and corresponding trellis."""

    allStates = A.keys()
    print_params(A, E)
    emittingStates = E.keys()
    L = len(X) + 2

    # Initialize
    F = {k:[0] * L for k in allStates}
    F['B'][0] = 1


    # HINT: The Viterbi and Forward algorithm are very similar! 
    # Adapt the viterbi() function to account for the differences.

    # Middle columns
    for i,s in enumerate(X):
        for l in emittingStates:
            terms = [F[k][i] * A[k][l] for k in allStates]
            F[l][i+1] = sum(terms) * E[l][s]
        
    # Last columns
    for k in allStates:
        term = F[k][i+1] * A[k]['E'] 
        F['E'][-1] += term
  
        


    P = F['E'][-1] # The Forward probability: P(X|A,E)
    return(P,F)



def backward(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the Backward probability and corresponding trellis."""

    allStates = A.keys()
    emittingStates = E.keys()
    L = len(X) + 2

    # Initialize
    B = {k:[0] * L for k in allStates} # The Backward trellis
    for k in allStates:
        B[k][-2] = A[k]['E']
   

    #####################
    # START CODING HERE #
    #####################
    # Remaining columns
    for i in range(L-3,-1,-1):
        #Value of the following column
        s = X[i]
        for k in allStates:
            # Nested, checking the previous states to make the calculations         
            terms = [B[l][i+1] * A[k][l] * E[l][s] for l in emittingStates] # k: L and D state
            B[k][i] = sum(terms)
                 
            

    #####################
    #  END CODING HERE  #
    #####################

    P = B['B'][0] # The Backward probability -- should be identical to Forward!
    return(P,B)



def baumwelch(set_X,A,E):
    """Given a set of sequences X and priors A and E,
    return the Sum Log Likelihood of X given the priors,
    along with the calculated posteriors for A and E."""

    allStates = A.keys()
    emittingStates = E.keys()
    
    # Initialize a new (posterior) Transition and Emission matrix
    new_A = {}
    for k in A:
        new_A[k] = {l:0 for l in A[k]}
    
    new_E = {}
    for k in E:
        new_E[k] = {s:0 for s in E[k]}
        
    # Iterate through all sequences in X
    SLL = 0 # Sum Log-Likelihood
    for X in set_X:
        P,F = forward(X,A,E)  # Save both the forward probability and the forward trellis
        _,B = backward(X,A,E) # Forward P == Backward P, so only save the backward trellis
        SLL += log10(P)

    
        # Inside the for loop: Expectation
        # Calculate the expected transitions and emissions for the sequence.
        # Add the contributions to your posterior matrices.
        # Remember to normalize to the sequence's probability P!
        
        sum_An = []
        for k in allStates:
            SA=0
            for l in allStates:                
                    sum_A=0
                    #Calculating the state transitions
                    #if l is also in emitting states do the calculation
                    if l in emittingStates:  
                        for i in range(0,len(X)):
                            A_prob = F[k][i] * A[k][l] * E[l][X[i]] * B[l][i+1] 
                            sum_A += A_prob 
                    # Emission for silent states = 0 so: the sum will be plus zero 
                    elif k=="B":
                        sum_A += 0 
                    # If it isnt in the emmission states    
                    else: #For column E(end state)
                        #print("IN")
                        A_prob = F[k][len(X)] * A[k][l]  
                        sum_A = A_prob
                    new_A[k][l] += sum_A/P
                    SA += new_A[k][l] #summing the Ak
            #Append denomination for normalising the Ak        
            sum_An.append(SA)
                    
                    
                    
        #Calculate emissions
        for e in emittingStates:               
            for i,s in enumerate(X):    
                # e passes through emission states, while s and i gives us the correct positions given X, then we divide with P following the formula
                new_E[e][s] += (F[e][i+1] * B[e][i+1]) / P
    # Outside the for loop: Maximization
    # Normalize row sums to 1 (except for one row in the Transition matrix!)
    # new_A = ...
    # new_E = ...
    
    #Normalise transition states rows adding equals 1
    for i,k in enumerate(allStates):
        for l in allStates:
            if sum_An[i] != 0:  # we cant divide if its zero
                new_A[k][l] = new_A[k][l] / sum_An[i] #our new Ak divided by the sum of An of the proper sequence state
    
    #Normalise emissions rows adding equals to 1
    for i,e in enumerate(emittingStates):
        sum_En = 0 # make a sum total
        for s in new_E[e]:
            sum_En += new_E[e][s]
        if sum_En != 0:  # we cant divide if its zero
            for s in new_E[e]:
                new_E[e][s] = new_E[e][s] / sum_En 
                



    return(SLL,new_A,new_E)



def main(args = False):
    "Perform the specified algorithm, for a given set of sequences and parameters."
    
    # Process arguments and load specified files
    if not args: args = parse_args()

    cmd = args.command            # viterbi, forward, backward or baumwelch
    verbosity = args.verbosity
    set_X, labels = load_fasta(args.fasta)  # List of sequences, list of labels
    A = load_tsv(args.transition) # Nested Q -> Q dictionary
    E = load_tsv(args.emission)   # Nested Q -> S dictionary
    
    def save(filename, contents):
        if args.out_dir:
            makedirs(args.out_dir, exist_ok=True) # Make sure the output directory exists.
            path = op.join(args.out_dir,filename)
            with open(path,'w') as f: f.write(contents)
        # Note this function does nothing if no out_dir is specified!



    # VITERBI
    if cmd == 'viterbi':
        for j,X in enumerate(set_X): # For every sequence:
            # Calculate the most probable state path, with the corresponding probability and matrix
            Q, P, T = viterbi(X,A,E)

            # Save and/or print relevant output
            label = labels[j]
            save('%s.path' % label, Q)
            save('%s.matrix' % label, serialize(T,X))
            save('%s.p' % label, '%1.2e' % P)
            print('>%s\n Path = %s' % (label,Q))
            if verbosity: print(' Seq  = %s\n P    = %1.2e\n' % (X,P))
            if verbosity >= 2: print_trellis(T, X)
            


    # FORWARD or BACKWARD
    elif cmd in ['forward','backward']:
        if cmd == 'forward':
            algorithm = forward
        elif cmd == 'backward':
            algorithm = backward

        for j,X in enumerate(set_X): # For every sequence:
            # Calculate the Forward/Backward probability and corresponding matrix
            P, T = algorithm(X,A,E)

            # Save and/or print relevant output
            label = labels[j]
            save('%s.matrix' % label, serialize(T,X))
            save('%s.p' % label, '%1.2e' % P)
            if verbosity >= 2:
                print('\n>%s\n P = %1.2e\n' % (label,P))
                print_trellis(T, X)
            elif verbosity: print('>%-10s\tP = %1.2e' % (label,P))



    # BAUM-WELCH TRAINING
    elif cmd == 'baumwelch':
        # Initialize
        i = 1
        i_max = args.max_iter
        threshold = args.conv_thresh

        current_SLL, A, E = baumwelch(set_X,A,E)
        if verbosity: print('Iteration %i, prior SLL = %1.2e' % (i,current_SLL))
        if verbosity >= 2: print_params(A,E)
        
        last_SLL = current_SLL - threshold - 1 # Iterate at least once

        # Iterate until convergence or limit
        while i < i_max and current_SLL - last_SLL > threshold:
            i += 1
            last_SLL = current_SLL

            # Calculate the Sum Log-Likelihood of X given A and E,
            # and update the estimates (posteriors) for A and E.
            current_SLL, A, E = baumwelch(set_X,A,E)

            if verbosity: print('Iteration %i, prior SLL = %1.2e' % (i,current_SLL))
            if verbosity >= 2: print_params(A,E)

        converged = current_SLL - last_SLL <= threshold
        try:
            final_SLL = sum([log10(forward(X,A,E)[0]) for X in set_X])
        except ValueError:
            final_SLL = 0

        # Save and/or print relevant output
        save('SLL','%1.2e\t%i\t%s' % (final_SLL, i, converged))
        save('posterior_A',serialize(A))
        save('posterior_E',serialize(E))
        if verbosity: print('========================================\n')

        if converged:
            print('Converged after %i iterations.' % i)
        else:
            print('Failed to converge after %i iterations.' % i_max)

        if verbosity:
            print('Final SLL: %1.2e' % final_SLL)
            print('Final parameters:')
            print_params(A,E)



if __name__ == '__main__':
	main()

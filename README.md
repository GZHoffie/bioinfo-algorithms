# CS4330 Simple Algorithm Implementations
Simple testing implementation of several algorithms in bioinformatics taught in the NUS Module of CS4330, *Combinatorial Methods in Bioinformatics* by Dr. Ken Sung.

- [Sequence Alignment](./align.py): Needleman-Wunsch algorithm for aligning two DNA sequences. Supports changing of scoring scheme and simple/affine gap penalties.
- [Profile-profile alignment](./MSA.py): do profile-profile alignment in a step of a multiple sequence alignment using the PSSM matrix for speeding up.
- [Partial Order Multiple Sequence Alignment (PO-MSA)](./POG.py): creates partial order graphs (POG) for alignments and align two POGs.
- [Longest Common Subsequence](./LCS.py): calculates the length of longest common subsequence of two strings using dynamic programming.
- [Markov Chain](./markov.py): several functionalities of markov chain, including
   - calculating stationary distribution,
   - calculating probability of observing a set of states
   and Hidden Markov Model (HMM),
   - calculating emit probability given a sequence of observation,
   - forward/backward/Viterbi algorithm,
   - posterior decoding.

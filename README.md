# HMMs
This project is about the implementation of some Hidden Markov Models as part of the A.I. course at KTH.
Four different HMMs are implemented in C++

<In extension, HMM4 (the complete algorithm) could be used to program the player class in games such as Duck Hunt.>

# HMM 1 - NEXT OBSERVATION DISTRIBUTION

Input: three matrices (in this order); transition matrix, emission matrix, and initial state probability distribution. The initial state probability distribution is a row vector encoded as a matrix with only one row. Each matrix is given on a separate line with the number of rows and columns followed by the matrix elements (ordered row by row).
Output: the emission probability distribution on a single line in the same matrix format, including the dimensions.

More details here:
https://kth.kattis.com/problems/kth.ai.hmm1

# HMM2 Probability of Emission Sequence

Input: three matrices; transition matrix, emission matrix, and initial state probability distribution followed by the number of emissions and the sequence of emissions itself. The initial state probability distribution is a row vector encoded as a matrix with only one row. Each matrix is given on a separate line with the number of rows and columns followed by the matrix elements (ordered row by row). It is assumed that there are M different discrete emission types and these are indexed 0 through M-1 in the emission sequence. For example, if there were M=3 possible different emissions (could be the three colours red, green and blue for example), they would be identified by 0, 1 and 2 in the emission sequence.
Output: the prbability of the given sequence as a single scalar.

More details here:
https://kth.kattis.com/problems/kth.ai.hmm2

# HMM3 Estimate Sequence of States

Input: three matrices; transition matrix, emission matrix, and initial state probability distribution followed by the number of emissions and the sequence of emissions itself. The initial state probability distribution is a row vector encoded as a matrix with only one row. Each matrix is given on a separate line with the number of rows and columns followed by the matrix elements (ordered row by row). It is assumed that there are M different discrete emission types and these are indexed 0 through M-1 in the emission sequence. For example, if there were M=3 possible different emissions (could be the three colours red, green and blue for example), they would be identified by 0, 1 and 2 in the emission sequence.
Output: the most probable sequence of states as zero-based indices separated by spaces. Do not output the length of the sequence.

More details here:
https://kth.kattis.com/problems/kth.ai.hmm3

# HMM4 Estimate Model (the complete project algorithm)

Input: a starting guess of the three matrices; transition matrix, emission matrix, and initial state probability distribution followed by the number of emissions and the sequence of emissions itself. The initial state probability distribution is a row vector encoded as a matrix with only one row. Each matrix is given on a separate line with the number of rows and columns followed by the matrix elements (ordered row by row). It is assumed that there are M different discrete emission types and these are indexed 0 through M-1 in the emission sequence. For example, if there were M=3 possible different emissions (could be the three colours red, green and blue for example), they would be identified by 0, 1 and 2 in the emission sequence.
Output: the estimated transition matrix and emission matrix on one line each in the same matrix format as they were given, including the dimensions. Do not output the estimated initial state probability distribution.

More details here:
https://kth.kattis.com/problems/kth.ai.hmm4

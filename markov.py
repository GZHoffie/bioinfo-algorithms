import numpy as np

def DNAToState(dna):
    state = []
    for c in dna:
        if c == 'A':
            state.append(0)
        elif c == 'C':
            state.append(1)
        elif c == 'G':
            state.append(2)
        else:
            state.append(3)
    return state

class MarkovModel:
    def __init__(self, pi, alpha, epsilon=None) -> None:
        self.pi = pi
        self.alpha = alpha
        self.epsilon = epsilon

    def transitionProbability(self, t):
        return np.linalg.matrix_power(self.alpha, t)
    
    def probabilityDistribution(self, t):
        return np.dot(self.pi, self.transitionProbability(t))

    def stationaryDistribution(self):
        return self.probabilityDistribution(10000)
    

    def stateProbability(self, observedStates):
        prob = 1
        previousState = None
        for state in observedStates:
            if previousState is None:
                prob *= self.pi[0, state]
            else:
                prob *= self.alpha[previousState, state]
            
            previousState = state
        
        return prob


    def forward(self, emits):
        numStates = self.alpha.shape[0]
        F = np.zeros((numStates, len(emits)))
        for i in range(len(emits)):
            for q in range(numStates):
                if i == 0:
                    # base case
                    F[q, i] = self.pi[0, q] * self.epsilon[q][emits[i]]
                else:
                    for s in range(numStates):
                        F[q, i] += F[s, i-1] * self.alpha[s, q]

                    F[q, i] *= self.epsilon[q][emits[i]]
        
        return F, np.sum(F[:, -1])
    
    def backward(self, emits):
        numStates = self.alpha.shape[0]
        B = np.zeros((numStates, len(emits)))
        for i in reversed(range(len(emits))):
            for q in range(numStates):
                if i == len(emits) - 1:
                    # base case
                    B[q, i] = 1
                else:
                    for s in range(numStates):
                        B[q, i] += B[s, i+1] * self.alpha[q, s] * self.epsilon[s][emits[i+1]]

        prob = 0
        for q in range(numStates):
            prob += self.pi[0, q] * self.epsilon[q][emits[0]] * B[q, 0]
        return B, prob
    

    def viterbi(self, emits):
        numStates = self.alpha.shape[0]
        V = np.zeros((numStates, len(emits)))
        states = np.zeros((numStates, len(emits)))
        for i in range(len(emits)):
            for q in range(numStates):
                if i == 0:
                    # base case
                    V[q, i] = self.pi[0, q] * self.epsilon[q][emits[i]]
                else:
                    for s in range(numStates):
                        newJointProbability = V[s, i-1] * self.alpha[s, q] * self.epsilon[q][emits[i]]
                        if newJointProbability > V[q, i]:
                            V[q, i] = newJointProbability
                            states[q, i] = s

        
        return V, states
    
    def posteriorDecoding(self, emits):
        F, pX = self.forward(emits)
        B, _ = self.backward(emits)
        return np.multiply(F, B) / pX





if __name__ == "__main__":
    MM = MarkovModel(pi=np.array([[0.5, 0.5]]), alpha=np.array([[0.7, 0.3], [0.3, 0.7]]),
                    epsilon=[np.array([0.4, 0.1, 0.4, 0.1]), np.array([0.25, 0.25, 0.25, 0.25])])
    print(MM.probabilityDistribution(2))
    print(MM.stationaryDistribution())
    print(MM.backward(DNAToState("CTAG")))
    print(MM.posteriorDecoding(DNAToState("CTAG")))

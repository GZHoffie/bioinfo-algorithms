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

class PairHMM:
    def __init__(self, q, p, eta=0.01, gamma=0.1, epsilon=0.2) -> None:
        self.q = q
        self.p = p
        self.eta = eta
        self.gamma = gamma
        self.epsilon = epsilon
    
    def viterbi(self, A, B):
        def delta(x, y):
            return np.log2(self.p[x, y]/self.q[x]/self.q[y]) + np.log2((1-2*self.gamma)/(1-self.eta)/(1-self.eta))

        h = np.log2(self.gamma / self.epsilon)
        s = np.log2(self.epsilon / (1-self.eta))
        c = np.log2((1-2*self.gamma) / (1-self.epsilon))
        print(h, s, c)

        # Base case
        V_M = np.zeros((len(A) + 1, len(B) + 1))
        W_X = np.zeros(V_M.shape)
        W_Y = np.zeros(V_M.shape)

        
        for j in range(V_M.shape[1]):
            W_X[0, j] = float('-inf')
            V_M[0, j] = float('-inf')
            W_Y[0, j] = h - c + j * s
        for i in range(V_M.shape[0]):
            W_Y[i, 0] = float('-inf')
            V_M[i, 0] = float('-inf')
            W_X[i, 0] = h - c + i * s
        
        V_M[0, 0] = -2 * np.log2(self.eta)
        

        for i in range(1, V_M.shape[0]):
            for j in range(1, V_M.shape[1]):
                V_M[i, j] = delta(A[i-1], B[j-1]) + max(V_M[i-1, j-1],
                                                    W_X[i-1, j-1],
                                                    W_Y[i-1, j-1])
                W_X[i, j] = max(V_M[i-1, j] + h - c + s, W_X[i-1, j] + s)
                W_Y[i, j] = max(V_M[i, j-1] + h - c + s, W_X[i, j-1] + s)

        print(V_M)
        print(W_X)
        print(W_Y)



if __name__ == "__main__":
    HMM = PairHMM(q=np.array([0.25, 0.25, 0.25, 0.25]), p=np.array([[0.13, 0.04, 0.04, 0.04],
                                                                    [0.04, 0.13, 0.04, 0.04],
                                                                    [0.04, 0.04, 0.13, 0.04],
                                                                    [0.04, 0.04, 0.04, 0.13]]))

    HMM.viterbi(DNAToState("ACCGA"), DNAToState("AGTTA"))
    

import numpy as np

class LongestCommonSubsequence:
    def __init__(self, A, B) -> None:
        """
        A and B are two lists of indices of distinct sequence.
        """
        self.A = A
        self.B = B
        self.C = np.zeros((len(A) + 1, len(B) + 1))
    
    def delta(self, i):
        """
        return the index of the character in B such that
        B[delta(i)] = A[i]
        """
        return self.B.index(self.A[i-1])+1
    
    def run(self):
        for i in range(len(self.A)):
            for j in range(len(self.B)):
                if j + 1 >= self.delta(i + 1):
                    self.C[i+1, j+1] = max(self.C[i, j+1], 1 + self.C[i, self.delta(i+1)-1])
                else:
                    self.C[i+1, j+1] = self.C[i, j+1]
        
        print(self.C)


def LCS(S, T):
    lcs = np.zeros((len(S) + 1, len(T) + 1))
    for i in range(1, len(S) + 1):
        for j in range(1, len(T) + 1):
            if S[i-1] == T[j-1]:
                lcs[i, j] = 1 + lcs[i-1, j-1]
            else:
                lcs[i, j] = max(lcs[i-1, j], lcs[i, j-1])

    print(lcs) 


if __name__ == "__main__":
    #LCS = LongestCommonSubsequence([1,2,3,4,5,6,7,8], [4,1,3,2,5,7,6,8])
    #LCS.run()
    LCS("ACTGCATGACTGGTA", 
        "TACCAGTCATGCAGT")
    LCS("ACCATATGGTA", 
        "ACCATATGGT")

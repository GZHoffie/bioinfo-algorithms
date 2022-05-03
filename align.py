import numpy as np

class DNAAlign:
    def __init__(self, dna1, dna2, match=2, mismatch=-1, indel=-1) -> None:
        self.dna1 = dna1
        self.dna2 = dna2
        self.score = {"match": match,
                      "mismatch": mismatch,
                      "indel": indel}
        self.m = len(dna1)
        self.n = len(dna2)
    
        self.D = np.zeros((self.m + 1, self.n + 1))
        self.alignment = {"dna1": "", "dna2": ""}
    

    def runDP(self):
        # Initialize
        for j in range(1, self.n + 1):
            self.D[0][j] = 0
        for i in range(1, self.m + 1):
            self.D[i][0] = 0
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                delete = self.D[i-1][j] + self.score["indel"]
                insert = self.D[i][j-1] + self.score["indel"]
                match = self.D[i-1][j-1]
                if self.dna1[i-1] != self.dna2[j-1]:
                    match += self.score["mismatch"]
                else:
                    match += self.score["match"]
                self.D[i][j] = max(insert, delete, match)
        
        return np.max(self.D)

    def backtrack(self):
        i = self.m
        j = self.n
        while i > 0 or j > 0:
            if i > 0 and j > 0 and ((self.dna1[i-1] == self.dna2[j-1] and self.D[i][j] == self.D[i-1][j-1] + self.score["match"])
                or (self.dna1[i-1] != self.dna2[j-1] and self.D[i][j] == self.D[i-1][j-1] + self.score["mismatch"])):
                # Match with same character
                self.alignment["dna1"] = self.dna1[i-1] + self.alignment["dna1"]
                self.alignment["dna2"] = self.dna2[j-1] + self.alignment["dna2"]
                i -= 1
                j -= 1
                
            elif i > 0 and self.D[i][j] == self.D[i-1][j] + self.score["indel"]:
                self.alignment["dna1"] = self.dna1[i-1] + self.alignment["dna1"]
                self.alignment["dna2"] = "-" + self.alignment["dna2"]
                i -= 1
                self.matching = False
            else:
                self.alignment["dna1"] = "-" + self.alignment["dna1"]
                self.alignment["dna2"] = self.dna2[j-1] + self.alignment["dna2"]
                j -= 1
                self.matching = False
            #print(self.alignment)


    def __str__(self):
        self.backtrack()
        return self.alignment["dna1"] + "\n" + self.alignment["dna2"]


if __name__ == "__main__":
    prob = DNAAlign("ACTCGATC", "GTTCGCCT")
    print(prob.runDP())
    print(prob.D)
    print(prob)

    

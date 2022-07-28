import enum
from typing import Dict
import networkx as nx
import numpy as np



def sequence2Graph(index, sequence):
    graph = nx.DiGraph()
    graph.add_node("start")
    previousNode = "start"
    for i, char in enumerate(sequence):
        graph.add_node((index, i, char))
        graph.add_edge(previousNode, (index, i, char))
        previousNode = (index, i, char)
    return graph


def delta(a, b):
    if (a == b):
        return 1
    else:
        return -1


class PartialOrderGraph:
    def __init__(self) -> None:
        self.pog = nx.DiGraph()
        self.currentIndex = 0
    
    def add(self, sequence):
        self.pog = sequence2Graph(self.currentIndex, sequence)
        self.currentIndex += 1
    
    def align(self, sequence):
        indel_cost = -2
        currentList = list(nx.topological_sort(self.pog))
        print(currentList)
        m = len(sequence)
        n = len(currentList) - 1
        V = np.zeros((m + 1, n + 1))
        B = np.zeros(V.shape, dtype='i, i')

        for j in range(1, n + 1):
            V[0, j] = float('-inf')
            for p in self.pog.predecessors(currentList[j]):
                print(j, currentList[j], p)
                V[0, j] = max(V[0, j], V[0, currentList.index(p)] + indel_cost)
        for i in range(1, m + 1):
            V[i, 0] = V[i-1, 0] + indel_cost

        for i in range(1, m + 1):
            for v in range(1, n + 1):
                V_temp = float('-inf')
                V_candidate = V[i-1, v] + indel_cost
                if V_candidate > V_temp:
                    V_temp = V_candidate
                    B[i, v] = (i-1, v)

                for u in self.pog.predecessors(currentList[v]):
                    #print(u, currentList[v])
                    if len(u) == 3:
                        V_candidate = V[i-1, currentList.index(u)] + delta(sequence[i-1], currentList[v][2])
                        if V_candidate > V_temp:
                            V_temp = V_candidate
                            B[i, v] = (i-1, currentList.index(u))
                        
                        V_candidate = V[i, currentList.index(u)] + indel_cost
                        if V_candidate > V_temp:
                            V_temp = V_candidate
                            B[i, v] = (i, currentList.index(u))
                
                V[i, v] = V_temp
        
        print(V)
        print(B)

        return V, B

    def alignGraph(self, graph):
        indel_cost = -2
        currentList = list(nx.topological_sort(self.pog))
        sequenceList = list(nx.topological_sort(graph))
        print(currentList)
        m = len(currentList) - 1
        n = len(sequenceList) - 1
        
        V = np.zeros((m + 1, n + 1))
        B = np.zeros(V.shape, dtype='i, i')

        for v in range(1, n + 1):
            V[0, v] = float('-inf')
            for q in self.pog.predecessors(currentList[v]):
                V[0, v] = max(V[0, v], V[0, q] + indel_cost)
        for u in range(1, m + 1):
            V[u, 0] = float('-inf')
            for p in self.pog.predecessors(currentList[v]):
                V[u, 0] = max(V[u, 0], V[p, 0] + indel_cost)
        for i in range(1, m + 1):
            for v in range(1, n + 1):
                V_temp = 0
                V_candidate = V[i-1, v] + indel_cost
                if V_candidate > V_temp:
                    V_temp = V_candidate
                    B[i, v] = (i-1, v)

                for u in self.pog.predecessors(currentList[v]):
                    print(u, currentList[v])
                    if len(u) == 3:
                        V_candidate = V[i-1, currentList.index(u)] + delta(sequence[i-1], currentList[v][2])
                        if V_candidate > V_temp:
                            V_temp = V_candidate
                            B[i, v] = (i-1, currentList.index(u))
                        
                        V_candidate = V[i, currentList.index(u)] + indel_cost
                        if V_candidate > V_temp:
                            V_temp = V_candidate
                            B[i, v] = (i, currentList.index(u))
                
                V[i, v] = V_temp
        
        print(V)
        print(B)

        return V, B

    


if __name__ == "__main__":
    POG = PartialOrderGraph()
    POG.add("ACGACGTA")
    POG.pog.add_edge((0, 2, 'G'), (1, 3, 'T'))
    POG.pog.add_edge((1, 3, 'T'), (0, 4, 'C'))
    POG.align("CGTAGTA")

    

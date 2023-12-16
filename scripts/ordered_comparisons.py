from abc import ABC, abstractmethod
from itertools import 
class Comparison(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __lt__(self):
        pass

    @abstractmethod
    def __eq__(self):
        pass

    @abstractmethod
    def __hash__(self):
        pass

def ListPairs(L, keep_same = True, keep_different = True):
    pairs = []
    for i in range(len(L)):
        for j in range(i, len(L)):
            if (i == j and keep_same) or (i != j and keep_different):
                pairs.append((L[i], L[j]))
    return pairs

def DictPairs(A, B):
    pairs = []
    for exp1, sams1 in A.items():
        for sam1 in sams1:
            for exp2, sams2 in B.items():
                for sam2 in sams2:
                    pairs.append((exp1, sam1, exp2, sam2))
    return pairs


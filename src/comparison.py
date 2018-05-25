import itertools
import numpy as np
from Bio import pairwise2
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

def generate_kmers(n, ks):
    chosen = np.random.choice(
        np.arange(0, np.sum(4**np.array(ks))-1),
        size=n, replace=False
    )
    iters = itertools.chain(
        *[itertools.product("ATGC", repeat=a) for a in ks]
    )
    ans = []
    i = 0
    for a in iters:
        if i in chosen:
            ans.append("".join(list(a)))
        if len(ans) == n:
            break
        i += 1
    return(ans)

def get_distances(kmers, f):
    return(
        np.array([f(*a) for a in itertools.combinations(kmers, 2)])
    )

def nw_distance(a,b):
    return(
        pairwise2.align.globalxx(a, b)[0][2]
    )

def compute_spearman(a, b):
    s = spearmanr(a,b)
    return(s[0], s[1])

def boxplot_of_distances(dnw,d2,name):
    to_bp = [d2[dnw == a] for a in np.unique(dnw)]
    plt.boxplot(to_bp, notch=True, bootstrap=1000)
    plt.xlabel("Needleman-Wunsch score")
    plt.ylabel(name)
    plt.grid(True)
    plt.show()
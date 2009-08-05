#!/usr/bin/env python
import gc
import sys

from collections import defaultdict

from csc import divisi

from lib import load_map

def main():
    dimmap = load_map.load_maps()

    num_users = len(dimmap['user_to'])
    num_repos = len(dimmap['repos_to'])

    user_counts = dimmap['user_counts']
    del dimmap
    gc.collect()

    A, user_vectors = build_matrix(sys.stdin, user_counts)
    

    svd = divisi.svd.svd_sparse(A, 200)
    print dir(svd)

    print num_users, num_repos


def build_matrix(input, user_counts):
    A = divisi.tensor.DictTensor(2)

    user_vectors = defaultdict(lambda : divisi.tensor.DictTensor(1))
    
    for i, line in enumerate(input):
        user, repos = map(int, line.rstrip().split(':'))
        user_vectors[user][repos] = A[(user, repos)] = 1.0 / user_counts[user]

        if i % 10000 == 0:
            debug("Finished with data line %s" % i)
    return A, user_vectors
    




def debug(s):
    sys.stderr.write("%s\n" % s)

if __name__ == '__main__':
    main()

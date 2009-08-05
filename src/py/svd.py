#!/usr/bin/env python
import gc
import sys
import time

from collections import defaultdict
from heapq import heappush, heapreplace

import numpy

from csc import divisi

from lib import load_map

def main():
    dimmap = load_map.load_maps()

    num_users = len(dimmap['user_to'])
    num_repos = len(dimmap['repos_to'])

    user_counts = dimmap['user_counts']
    del dimmap
    gc.collect()

    datafile = open('../../dat/mapped_data.txt')
    A, user_vectors = build_matrix(datafile, user_counts)
    datafile.close()

    svd = divisi.svd.svd_sparse(A, 250)
    debug("SVD FINISHED!")

    main_vector = numpy.zeros(num_repos)
    num = 0
    svd_matrix = svd.weighted_v.transpose()
    for user in sys.stdin:
        start_time = time.time()
        user = int(user.rstrip())
        try:
            user_vector =  user_vectors[user]
        except KeyError:
            continue
        main_vector[:] = 0

        for key, value in user_vector.iteritems():
            main_vector[key] = value
        print get_suggestions(user, main_vector, svd_matrix)
        num += 1
        t = time.time() - start_time
        if num % 100 == 0:
            debug("%s took %0.3fs %0.2f min remaining." % (
                    user, t, (num_users - num) * t / 60.0)
                  )

def get_suggestions(user, user_vector, svd, top_n=50):
    heap = []
    start = time.time()

    #debug(user_vector.shape)    
    #debug(svd.shape)
    repos = svd * user_vector #.v_dotproducts_with(user_vector)

    #debug("dot products took %0.4fs" % (time.time() - start), 4)

    i = 0
    for repo, value in repos.iteritems():
        i += 1
        if not value:
            continue

        newval = (value, repo)
        if not heap:
            heap.append(newval)
        elif len(heap) == top_n:
            if newval > heap[0]:
                heapreplace(heap, newval)
        else:
            heappush(heap, newval)
        #if i % 10000 == 0:
        #    debug("Finished with repo %s" % i)
    heap.sort(reverse=True)
    return "%s:%s" % (user, ','.join("%s;%0.4f" % (x[1][0], x[0]) for x in heap))

def build_matrix(input, user_counts):
    A = divisi.tensor.DictTensor(2)

    user_vectors = defaultdict(lambda : {})
    
    for i, line in enumerate(input):
        user, repos = map(int, line.rstrip().split(':'))
        user_vectors[user][repos] = A[(user, repos)] = 1.0 / user_counts[user]

        if i % 10000 == 0:
            debug("Finished with data line %s" % i)
    return A, user_vectors


def debug(s, level=0):
    sys.stderr.write("%s\n" % str(s))

if __name__ == '__main__':
    main()

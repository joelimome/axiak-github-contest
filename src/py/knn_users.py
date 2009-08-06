#!/usr/bin/python
import sys
import os
import gc
import time
import random
import math

import psyco
psyco.full()

from collections import defaultdict

from lib import load_map, neighbors

datadir = os.path.join(os.path.dirname(__file__), '..', '..', 'dat')

def main():
    dimmap = load_map.load_maps()

    k = 150
    if len(sys.argv) > 1:
        k = int(sys.argv[0])

    test_users = set(int(line.rstrip()) for line in sys.stdin)

    edgefile = open(os.path.join(datadir, 'edge_counts.dat'))
    neighbor_values = neighbors.get_neighbors(edgefile,
                                              dimmap['user_counts'],
                                              dimmap['user_to'],
                                              test_users,
                                              k=k)
    edgefile.close()
    del dimmap
    gc.collect()

    datafile = open(os.path.join(datadir, 'mapped_data.txt'))
    suggestions = get_suggestions(datafile, neighbor_values, test_users)
    datafile.close()

    print_suggestions(suggestions, top_n=10)

def get_suggestions(datafile, neighbor_values, test_users):
    inverse_map = defaultdict(list)
    neighbor_counts = {}

    suggestions = defaultdict(lambda : defaultdict(float))

    for user, neighbors in neighbor_values.iteritems():
        neighbor_counts[user] = len(neighbors)
        for distance, neighbor in neighbors:
            inverse_map[neighbor].append((user, distance))

    num = 0
    for line in datafile:
        user, repo = map(int, line.rstrip().split(':'))
        try:
            helpers = inverse_map[user]
        except KeyError:
            continue
        for helper, distance in helpers:
            distance = -distance
            if distance == 0:
                distance = 0.3
            else:
                distance = math.sqrt(distance)
            if random.random() < 0.01:
                debug("d:%s" % distance)
            suggestions[helper][repo] += 1.0 / neighbor_counts[helper] / distance
        num += 1
        if num % 5000 == 0:
            debug("get_suggestions: Finished %s" % (num))

    return suggestions

def print_suggestions(suggestions, top_n):
    num = 0
    total_start = time.time()
    for user, repodata in suggestions.iteritems():
        repodata = repodata.items()
        repodata.sort(key=lambda x:x[1], reverse=True)
        print "%s:%s" % (user,
                         ",".join("%s;%0.4f" % x for x in repodata[:top_n]))
        num += 1
        if num % 100 == 0:
            total_time = time.time() - total_start
            remaining = (total_time / float(num)) * (len(suggestions) - num)
            debug ("Finished with %s [%s min remaining]" %
                   (num, remaining / 60.0))

def debug(s):
    sys.stderr.write("%s\n" % str(s))

if __name__ == '__main__':
    main()

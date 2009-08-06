#!/usr/bin/python
import sys
import os

from collections import defaultdict

def order(x, y):
    if x > y:
        return x, y
    else:
        return y, x

def main():
    edge_counts = defaultdict(int)
    repos = defaultdict(list)
    for line in sys.stdin:
        user, repo = map(int, line.rstrip().split(':', 1))
        repos[repo].append(user)

    for users in repos.itervalues():
        for i in xrange(len(users) - 1):
            for j in xrange(i + 1, len(users)):
                edge_counts[order(users[i], users[j])] += 1

    for edge, count in edge_counts.iteritems():
        print "%d:%d:%d" % (edge[0], edge[1], count)


if __name__ == '__main__':
    main()

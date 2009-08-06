#!/usr/bin/python
import sys
import os

import psyco
psyco.full()

from collections import defaultdict

def order(x, y):
    if x > y:
        return x, y
    else:
        return y, x

def main():
    edge_counts = defaultdict(lambda :defaultdict(int))
    repos = defaultdict(list)
    for i, line in enumerate(sys.stdin):
        user, repo = map(int, line.rstrip().split(':', 1))
        repos[repo].append(user)
        if i % 10000 == 0:
            debug("Done with line %s" % i)

    usernum = 0
    for users in repos.itervalues():
        if len(users) < 2:
            continue
        debug("%s - %s" % (len(users), (len(users) * (len(users) + 1)) / 2))
        tightnum = 0
        for i in xrange(len(users) - 1):
            for j in xrange(i + 1, len(users)):
                a, b = users[i], users[j]
                if b < a:
                    a, b = b, a
                edge_counts[a][b] += 1
                tightnum += 1
                if tightnum % 100000 == 0:
                    debug("Done with %s tightnum" % tightnum)

        usernum += 1
        if usernum % 1000 == 0:
            debug("Done with %s users" % usernum)

    for user1, vals in edge_counts.items():
        for user2, count in vals.items():
            print "%d:%d:%d" % (user1, user2, count)


def debug(s):
    sys.stderr.write("%s\n" % str(s))

if __name__ == '__main__':
    main()

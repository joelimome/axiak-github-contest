#!/usr/bin/env python

import sys
import os
import math

from collections import defaultdict
from numpy import *
from scipy import sparse

dir = os.path.dirname(__file__)

imatrix = sparse.lil_matrix
omatrix = sparse.csr_matrix

already_seen = set()

def debug(s):
    sys.stderr.write("%s\n" % str(s))

def create_matrix(datafile=sys.stdin):
    users_watches = defaultdict(lambda : defaultdict(int))
    all_repos = set()

    for line in datafile:
        user, repo = map(int, line.strip().split(':'))
        users_watches[user][repo] = 1
        all_repos.add(repo)
        already_seen.add((user, repo))

    all_repos = list(all_repos)
    all_users = list(users_watches)

    repo_map = dict(((y, x) for x, y in enumerate(all_repos)))
    user_rmap = dict()

    debug((len(users_watches), len(all_repos)))
    matrix = imatrix((len(users_watches), len(all_repos)),
                     dtype=single)
    
    i = 0
    debug("updating..")
    user_map = {}
    for user, watches in users_watches.iteritems():
        user_rmap[i] = user
        user_map[user] = i
        num = 1 / math.sqrt(float(len(watches)))
        for repo in watches:
            matrix[i, repo_map[repo]] = num
        i += 1

    repo_rmap = dict(((y, x) for x, y in repo_map.iteritems()))

    debug("converting...")
    matrix = omatrix(matrix)
    debug("Multiplying...")
    return transpose(matrix) * matrix, matrix, user_map, user_rmap, repo_rmap

def get_recommendations(matrix, vectors, userid, user_rmap, repo_rmap):
    user = user_rmap[userid]

    vector = transpose(vectors[userid])
    debug("multiplying %s.." % userid)
    repos = imatrix(matrix * vector)
    debug("getting data for %s.." % userid)
    repos = [(repo_rmap[i], repos[i,0]) for i in xrange(len(repo_rmap))
             if repos[i,0] > 0.05 and (user, repo_rmap[i]) not in already_seen]
    debug("sorting for %s.." % userid)
    repos.sort(key=lambda x: x[1], reverse=True)
    return "%s:%s" % (user,
                     ','.join(str(repo[0]) for repo in repos[:10]))
    


if __name__ == '__main__':
    datafile = open(os.path.join(dir, '..', 'input', 'data.txt'))
    full, piece, user_map, user_rmap, repo_rmap = create_matrix(datafile)
    datafile.close()

    for line in sys.stdin:
        user = int(line.strip())
        if user not in user_map:
            print "%s:" % user
            continue
        userid = user_map[user]
        print get_recommendations(full, piece, userid, user_rmap, repo_rmap)
        sys.stdout.flush()

    debug('done!')


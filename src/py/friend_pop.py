#!/usr/bin/env python
"""
Using friend's popularity markings...
"""
import sys
import os
import time
import psyco
psyco.full()


from collections import defaultdict

from lib import datainput
from lib import normalize
from lib import loader

datadir = os.path.join(os.path.dirname(__file__), '..', '..', 'input')

def debug(s):
    sys.stderr.write("%s\n" % s)

def get_friend_volume(tests, datafile=None, usermap=lambda x: x):
    test_set = set(tests)
    if not datafile:
        datafile = open(os.path.join(datadir, 'edge_num.dat'))

    volume_data = defaultdict(dict)

    for line in datafile:
        user1, user2, count = map(int, line.strip().split(':'))
        user1 = usermap(user1)
        user2 = usermap(user2)

        if user1 in test_set:
            volume_data[user1][user2] = count

        if user2 in test_set:
            volume_data[user2][user1] = count

    return volume_data

def get_suggestions(user, volume_data, data_info, user_rmap=lambda x: x, answerfilter=lambda: True, repo_rmap=lambda x: x):
    try:
        volume_data = volume_data[user]
    except KeyError:
        debug("No volume info for %s" % user)
        return "%s:" % user_rmap(user)

    repo_scores = defaultdict(int)
    user_repos = data_info['user_repos']

    user_num_repos = len(user_repos.get(user, ()))

    for new_user, weight in volume_data.items():
        try:
            cur_repos = user_repos[new_user]
        except KeyError:
            continue

        if not cur_repos:
            continue

        # Compute the weighting.
        divisor = float((len(cur_repos) + user_num_repos - 2 * weight))

        # There is nothing for this user to contribute.
        if divisor <= 0:
            debug("Exiting due to complete overlap for %s" % user_rmap(user))
            debug("(%s,%s) (%s,%s): %s:: %s" % (
                    user_rmap(user), user_rmap(new_user),
                    len(cur_repos), user_num_repos, weight,
                    divisor))
            continue

        weight /= divisor

        for repo in cur_repos:
            if not answerfilter((user, repo)):
                continue
            repo_scores[repo] += weight

    normalize.normalize_dict(repo_scores)
    repos = repo_scores.items()
    repos.sort(key=lambda x: x[1], reverse=True)
    return "%s:%s" % (
        user_rmap(user),
        ','.join(["%s;%s" % x for x in repos[:100]])
        )

    answer = []
    for repo in repos:
        repo = repo[0]
        if answerfilter((user, repo)):
            answer.append(repo)
        if len(answer) >= 10:
            break

    return "%s:%s" % (
        user_rmap(user),
        ','.join(map(str, answer)),
        )


def main():
    datafile = open(os.path.join(datadir, 'data.txt'))
    data_info = datainput.get_datainfo(datafile)
    datafile.close()

    usermap = data_info['user_map'].__getitem__
    user_rmap = data_info['user_rmap'].__getitem__
    repo_rmap = data_info['repo_rmap'].__getitem__
    answerfilter = lambda x: x not in data_info['user_repo_set']

    tests, test_data = loader.load_tests(usermap=usermap)

    volume_data = get_friend_volume(tests, usermap=usermap)

    for user in tests:
        print get_suggestions(user, volume_data, data_info,
                              user_rmap=user_rmap,
                              answerfilter=answerfilter,
                              repo_rmap=repo_rmap)
        sys.stdout.flush()



if __name__ == '__main__':
    main()

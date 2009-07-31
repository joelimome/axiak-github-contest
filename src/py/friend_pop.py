#!/usr/bin/env python
"""
Using friend's popularity markings...
"""
import sys
import os

from collections import defaultdict

from lib import datainput

datadir = os.path.join(os.path.dirname(__file__), '..', '..', 'input')

def debug(s):
    sys.stderr.write("%s\n" % s)

def load_tests(datafile=sys.stdin, usermap=lambda x: x):
    return [usermap(int(line.strip())) for line in datafile]

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

def get_suggestions(user, volume_data, data_info, user_rmap=lambda x: x, answerfilter=lambda: True):
    try:
        volume_data = volume_data[user]
    except KeyError:
        debug("No volume info for %s" % user)
        return "%s:" % user_rmap(user)

    repo_scores = defaultdict(int)
    user_repos = data_info['user_repos']

    user_num_repos = len(user_repos.get(user, ()))

    for user, weight in volume_data.iteritems():
        try:
            cur_repos = user_repos[user]
        except KeyError:
            continue

        if not cur_repos:
            continue

        # There is nothing for this user to contribute.
        if not (len(cur_repos) + user_num_repos - 2 * weight):
            continue

        # Compute the weighting.
        divisor = float((len(cur_repos) + user_num_repos - 2 * weight)**2)

        weight /= divisor

        for repo in cur_repos:
            if not answerfilter((user, repo)):
                continue
            repo_scores[repo] += weight

    repos = repo_scores.items()
    repos.sort(key=lambda x: x[1], reverse=True)
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

    usermap = lambda x: data_info['user_map'].get(x, -1)
    user_rmap = lambda x: data_info['user_rmap'].get(x, -1)
    answerfilter = lambda x: x not in data_info['user_repo_set']

    tests = load_tests(usermap=usermap)

    volume_data = get_friend_volume(tests, usermap=usermap)

    for user in tests:
        print get_suggestions(user, volume_data, data_info,
                              user_rmap=user_rmap,
                              answerfilter=answerfilter)
        sys.stdout.flush()



if __name__ == '__main__':
    main()

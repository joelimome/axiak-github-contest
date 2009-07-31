#!/usr/bin/env python
"""
Using friend's popularity markings...
"""
import sys
import os
import math
import time

from numpy import *

#import psyco
#psyco.profile()

from collections import defaultdict

from lib import datainput
from lib import normalize
from lib import loader

datadir = os.path.join(os.path.dirname(__file__), '..', '..', 'input')
NUM_LANGUAGES = 40

USE_TEST = True
ALPHA = 0.24

def debug(s, user=None):
    if user == 36 or user is None:
        sys.stderr.write("%s\n" % s)

def get_lang_info(datafile=None, repo_map=lambda x: x):
    repo_lang = {}
    repo_tloc = defaultdict(int)
    lang_repos = defaultdict(list)

    repo_langs = defaultdict(list)

    lang_map = datainput.dictmap()

    if not datafile:
        datafile = open(os.path.join(datadir, 'lang.txt'))

    for line in datafile:
        repo, langinfo = line.strip().split(':', 1)
        try:
            repo = repo_map(int(repo))
        except:
            continue
        lang_info = ndarray(shape=(NUM_LANGUAGES, 1),
                            dtype=float)
        lang_info[...] = 0
        tloc = 0
        for lang in langinfo.split(','):
            lang, num = lang.split(';', 1)
            num = int(num)

            langid, created = lang_map.addrow(lang)
            lang_info[langid, 0] = num
            tloc += num
            repo_langs[repo].append(lang)
            lang_repos[lang].append(repo)

        lang_info = lang_info ** 1.5
        lang_info = matrix(lang_info)

        repo_tloc[repo] = tloc
        norm = linalg.norm(lang_info)
        if norm >= 0.01:
            lang_info /= norm

        #debug(','.join(map(lambda x: "%.2f" % x[0,0], lang_info)))
        repo_lang[repo] = lang_info

    return {
        'repo_lang': repo_lang,
        'repo_tloc': repo_tloc,
        'lang_repos': lang_repos,
        'repo_langs': repo_langs,
        }

def dot_dicts(dict1, dict2):
    sum = 0
    for key, value1 in dict1.items():
        sum += value1 * dict2.get(key, 0)
    return sum

def get_suggestions(user, repo_data, data_info, user_rmap=lambda x: x, answerfilter=lambda: True, input_data=defaultdict(int), repo_rmap=lambda x: x):

    repo_scores = defaultdict(int)
    user_repos = data_info['user_repos']

    if user not in user_repos:
        return "%s:" % user_rmap(user)

    seen_languages = set()
    user_languages = []
    candidate_repos = set()
    user_repo_set = set(user_repos[user])

    for repo in user_repos[user]:
        try:
            lang_info = repo_data['repo_lang'][repo]
        except:
            continue
        user_languages.append(lang_info)
        for lang in repo_data['repo_langs'][repo]:
            if lang not in seen_languages:
                seen_languages.add(lang)
                for repo in repo_data['lang_repos'].get(lang, ()):
                    if not answerfilter((user, filter)):
                        continue
                    if repo not in user_repo_set:
                        candidate_repos.add(repo)

    user_matrix = matrix(ndarray(shape=(len(user_languages),
                                        NUM_LANGUAGES),
                                 dtype=float))

    debug(user_matrix, user_rmap(user))

    for i, vector in enumerate(user_languages):
        user_matrix[i] = transpose(vector)

    repo_scores = {}
    for repo in candidate_repos:
        repo_matrix = matrix(user_matrix)
        curinfo = repo_data['repo_lang'][repo]
        repo_matrix[...] -= transpose(curinfo)

        score = linalg.norm(repo_matrix)
        repo_scores[repo_rmap(repo)] = 1 / score

    del repo_matrix
    del user_matrix

    normalize.normalize_dict(repo_scores)
    if USE_TEST:
        for repo, score in repo_scores.iteritems():
            repo_scores[repo] = \
                score**(ALPHA) * input_data[user][repo] ** (1 - ALPHA)
            if repo_scores[repo]:
                debug("%s <%s,%s>" % (
                    repo_scores[repo],
                    score,
                    input_data[user][repo],),
                      user_rmap(user))


    normalize.normalize_dict(repo_scores)
    repos = repo_scores.items()
    repos.sort(key=lambda x: x[1], reverse=True)
    return "%s:%s" % (
        user_rmap(user),
        ','.join(("%s;%s" % x for x in repos[:100]))
        )



def main():
    datafile = open(os.path.join(datadir, 'data.txt'))
    data_info = datainput.get_datainfo(datafile)
    datafile.close()

    usermap = data_info['user_map'].__getitem__
    user_rmap = data_info['user_rmap'].__getitem__
    repo_rmap = data_info['repo_rmap'].__getitem__
    repo_map = data_info['repo_map'].__getitem__

    answerfilter = lambda x: x not in data_info['user_repo_set']

    tests, test_data = loader.load_tests(usermap=usermap)

    repo_data = get_lang_info(repo_map=repo_map)

    for user in tests:
        print get_suggestions(
            user,
            repo_data,
            data_info,
            user_rmap=user_rmap,
            answerfilter=answerfilter,
            input_data=test_data,
            repo_rmap=repo_rmap,
            )

        sys.stdout.flush()



if __name__ == '__main__':
    main()

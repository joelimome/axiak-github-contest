#!/usr/bin/env python
"""
Using friend's popularity markings...
"""
import sys
import os
import math
import time
#import psyco
#psyco.full()

from collections import defaultdict

from lib import datainput
from lib import normalize
from lib import loader

datadir = os.path.join(os.path.dirname(__file__), '..', '..', 'input')

USE_TEST = True
ALPHA = 0.2

def debug(s):
    sys.stderr.write("%s\n" % s)

def get_lang_info(datafile=None, repo_map=lambda x: x):
    repo_lang = defaultdict(dict)
    repo_tloc = defaultdict(int)
    lang_repos = defaultdict(list)

    if not datafile:
        datafile = open(os.path.join(datadir, 'lang.txt'))

    for line in datafile:
        repo, langinfo = line.strip().split(':', 1)
        try:
            repo = repo_map(int(repo))
        except:
            continue
        lang_info = {}
        tloc = 0
        for lang in langinfo.split(','):
            lang, num = lang.split(';', 1)
            num = int(num)
            lang_info[lang] = num
            tloc += num
            lang_repos[lang].append(repo)

        repo_tloc[repo] = tloc
        normalize.normalize_dict(lang_info, mapfunc=lambda x: x**2, endfunc=lambda x: math.sqrt(x))
        repo_lang[repo] = lang_info

    return {
        'repo_lang': repo_lang,
        'repo_tloc': repo_tloc,
        'lang_repos': lang_repos,
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
        #debug(lang_info)
        for lang in lang_info:
            if lang not in seen_languages:
                seen_languages.add(lang)
                for repo in repo_data['lang_repos'].get(lang, ()):
                    if repo not in user_repo_set:
                        candidate_repos.add(repo)

    repo_scores = {}
    for repo in candidate_repos:
        score = 0
        curinfo = repo_data['repo_lang'][repo]
        for data in user_languages:
            score += dot_dicts(curinfo, data)

        repo_scores[repo] = score

    normalize.normalize_dict(repo_scores)
    if USE_TEST:
        for repo, score in repo_scores.items():
            repo_scores[repo] = \
                ALPHA * score + (1 - ALPHA) * input_data[user][repo];

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

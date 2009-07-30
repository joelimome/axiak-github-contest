#!/usr/bin/env python
import sys
import os
import math
import datetime

from collections import defaultdict

from numpy import *

d = os.path.dirname(__file__)

project_sizes = {}
repo_data = {}
already_seen = set()
user_repos = defaultdict(list)
ALPHA = 0.8

def debug(s):
    sys.stderr.write("%s\n" % s)

def get_data(datafile=sys.stdin):
    
    user_repos = defaultdict()
    for line in datafile:
        user, repo = map(int, line.strip().split(':'))

def project_repos(repos, num_needed=10):
    cur_matrix = matrix(ndarray(shape=(len(repos), len(repo_data['lang_map'])),
                                dtype=float))

    cur_map = {}
    for i, repo in enumerate(repos):
        try:
            repo_lang = repo_data['all_repos'][repo]
        except:
            continue
        for lang, value in repo_lang.items():
            lang_num = repo_data['lang_map'][lang]
            cur_matrix[i, lang_num] = value

    pieces = cur_matrix * repo_data['repo_matrix']

    scores = []
    for i in range(len(repos)):
        score = 0
        for j in range(len(repo_data['lang_map'])):
            score += pieces[i, j]
        repo = repo_data['repo_map'][i]

        scores.append((repo,
                       score * ALPHA + project_sizes.get(repo, 0)))

    scores.sort(key=lambda x:x[1], reverse=True)
    result = []
    for repo, score in scores:
        if (user, repo) not in already_seen:
            result.append(repo)
        if len(result) >= num_needed:
            return result
    return result

def get_data():
    repo_map = repo_data['repo_rmap']
    datafile = open(os.path.join(d, '..', 'input', 'data.txt'))

    for line in datafile:
        user, repo = map(int, line.strip().split(':'))
        already_seen.add((user, repo))
        if repo in repo_map:
            user_repos[user].append(repo)
    datafile.close()

def get_repolangs():
    lang_map = {}
    lang_rmap = {}
    lang_num = 0

    all_repos = {}
    max_value = 0
    langs = open(os.path.join(d, '..', 'input', 'lang.txt'))
    for line in langs:
        lang_info = {}
        repo, rest = line.strip().split(':', 1)
        repo = int(repo)

        total = 0
        for lang in rest.split(','):
            name, lines = lang.split(';')
            if name in lang_info:
                continue
            if name not in lang_map:
                lang_map[name] = lang_num
                lang_rmap[lang_num] = name
                lang_num += 1
            lines = int(lines)
            total += lines ** 2
            lang_info[name] = lines
        if not total:
            continue
        total = math.sqrt(float(total))
        project_sizes[repo] = total
        for lang, value in lang_info.items():
            value /= total
            lang_info[lang] = value
        all_repos[repo] = lang_info
    langs.close()
    repos = list(all_repos)
    repo_map = dict(((i, x) for i, x in enumerate(repos)))
    repo_rmap = dict(((x, i) for i, x in enumerate(repos)))

    repo_matrix = matrix(ndarray(shape=(len(lang_map), len(repos)),
                                 dtype=float))

    for i, repo in enumerate(repos):
        for lang, value in all_repos[repo].items():
            repo_matrix[lang_map[lang], i] = value

    repo_data.update({
            'repo_matrix': repo_matrix,
            'all_repos': all_repos,
            'repo_map': repo_map,
            'repo_rmap': repo_rmap,
            'lang_map': lang_map,
            'lang_rmap': lang_rmap,
            })


# Deprecated
def get_repodata():
    repos_data = defaultdict(lambda: [None, defaultdict(int)])
    lang_repos = defaultdict(lambda: [])

    dates = open(os.path.join(d, '..', 'input', 'repos.txt'))
    for line in dates:
        repo, rest = line.strip().split(':', 1)
        repo = int(repo)
        info = rest.split(',')
        date = datetime.date(*map(int, info[1].split('-')))
        repos_data[repo][0] = int(date.strftime("%s"))
    dates.close()

    langs = open(os.path.join(d, '..', 'input', 'lang.txt'))
    for line in langs:
        lang_info = defaultdict(int)
        repo, rest = line.strip().split(':', 1)
        repo = int(repo)

        total = 0
        for lang in rest.split(','):
            name, lines = lang.split(';')
            if name in lang_info:
                continue
            lines = int(lines)
            total += lines
            lang_info[name] = lines

        for lang, lines in lang_info.items():
            if lines > float(total) / 2.0:
                lang_repos[lang].append(repo)
        repos_data[repo][1] = lang_info
    langs.close()

if __name__ == '__main__':
    get_repolangs()
    get_data()

    big_projects = [x[0] for x in 
                    sorted(project_sizes.items(), 
                           key=lambda x:x[1],
                           reverse=True)][:10]

    for line in sys.stdin:
        user, repos = line.strip().split(':')
        user = int(user)
        repos = map(int, filter(None, repos.split(',')))
        num_needed = 10 - len(repos)
        if not num_needed:
            sys.stdout.write(line)
            continue
        try:
            old_repos = user_repos[user]
        except KeyError:
            old_repos = []

        if not old_repos:
            debug("NO old repos for %s" % user)
            new_repos = []
        else:
            new_repos = project_repos(old_repos, num_needed)
            if not new_repos:
                debug("NO new repos for %s" % user)

        debug("GOT SOME FOR %s!" % user)
        repos.extend(new_repos)
        num_needed = 10 - len(repos)
        repos.extend(big_projects[:num_needed])
        print "%s:%s" % (user, ','.join(map(str, repos)))

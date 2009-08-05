
from collections import defaultdict
#import psyco
#psyco.full()


class dictmap(dict):
    def __init__(self, *args, **kwargs):
        self.__curnum = 0
        return dict.__init__(self, *args, **kwargs)

    def addrow(self, elem):
        if elem in self:
            return self[elem], False
        else:
            id = self.__curnum
            self.__curnum += 1
            self[elem] = id
            return id, True

class dummydictmap(dictmap):
    def addrow(self, elem):
        return elem, False

    def __getitem__(self, item):
        return item

    def __setitem__(self, key, value):
        pass

def get_datainfo(datafile, usermap=False):
    repo_pop = defaultdict(int)

    if usermap:
        user_map = dictmap()
        user_rmap = {}

        repo_map = dictmap()
        repo_rmap = {}
    else:
        user_map = dummydictmap()
        user_rmap = dummydictmap()

        repo_map = dummydictmap()
        repo_rmap = dummydictmap()

    user_repos = defaultdict(list)

    user_repo_set = set()

    user_num = 0
    for line in datafile:
        user, repo = map(int, line.strip().split(':'))

        userid, created = user_map.addrow(user)
        if created:
            user_rmap[userid] = user

        repoid, created = repo_map.addrow(repo)
        if created:
            repo_rmap[repoid] = repo

        repo_pop[repo] += 1

        user_repos[userid].append(repoid)

        user_repo_set.add((userid, repoid))

    return {
        'user_repos': user_repos,
        'user_map': user_map,
        'user_rmap': user_rmap,
        'repo_map': repo_map,
        'repo_rmap': repo_rmap,
        'repo_pop': repo_pop,
        'user_repo_set': user_repo_set,
        }

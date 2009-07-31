
from collections import defaultdict

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

def get_datainfo(datafile):
    repo_pop = defaultdict(int)

    user_map = dictmap()
    user_rmap = {}

    repo_map = dictmap()
    repo_rmap = {}

    user_repos = defaultdict(list)

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

        user_repos[user].append(repo)

    return {
        'user_repos': user_repos,
        'user_map': user_map,
        'user_rmap': user_rmap,
        'repo_map': repo_map,
        'repo_rmap': repo_rmap,
        'repo_pop': repo_pop,
        }

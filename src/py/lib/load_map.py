"""
Load the dimension mappings from the files in the "dat" directory.
"""
import os

file_dir = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'dat')

def load_maps(root_dir=file_dir):
    ucfile = open(os.path.join(file_dir, "usercount.dat"))
    user_counts = dict(map(lambda line: map(int, line.rstrip().split(":", 1)),
                           ucfile))

    ucfile.close()

    user_to = {}
    user_from = {}
    umfile = open(os.path.join(file_dir, "usermap.dat"))
    for line in umfile:
        user, userid = line.rstrip().split(":", 1)
        userid = int(userid)
        user_to[user] = userid
        user_from[userid] = user
    umfile.close()

    repos_to = {}
    repos_from = {}
    rmfile = open(os.path.join(file_dir, "repomap.dat"))
    for line in rmfile:
        repo, repoid = line.rstrip().split(":", 1)
        repoid = int(repoid)
        repos_to[repo] = repoid
        repos_from[repoid] = repo
    rmfile.close()

    del repo, ucfile, umfile, user, line, rmfile, repoid, userid, root_dir

    return locals()

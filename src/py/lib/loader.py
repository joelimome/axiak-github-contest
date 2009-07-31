import sys
from collections import defaultdict

def load_tests(datafile=sys.stdin, usermap=lambda x: x):
    result = []
    input_data = defaultdict(lambda : defaultdict(int))

    for line in datafile:
        if ':' not in line:
            try:
                result.append(usermap(int(line.strip())))
            except:
                pass
            continue
        user, data = line.strip().split(':', 1)
        try:
            user = usermap(int(user))
        except:
            pass
        result.append(user)

        if not data.strip():
            continue

        data = data.split(',')

        for item in data:
            try:
                repo, score = item.split(';', 1)
            except:
                raise Exception("%s %s" % (user, item))
            input_data[user][int(repo)] = float(score)

    return result, input_data

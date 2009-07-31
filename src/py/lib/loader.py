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
        data = data.split(',')
        result.append(user)
        for item in data:
            repo, score = item.split(';', 1)
            input_data[user][int(repo)] = float(score)

    return result, input_data

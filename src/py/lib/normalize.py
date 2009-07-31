def normalize_dict(d):
    max_val = 0
    for key, value in d.iteritems():
        if value > max_val:
            max_val = value

    if not max_val:
        return
    for key, value in d.iteritems():
        d[key] = value / float(max_val)

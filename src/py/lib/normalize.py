def normalize_dict(d, mapfunc=None, endfunc=None):
    """
    Normalize a dictionary.
    By default the normalization is the infty-norm.
    You can use a 2-norm by passing lambda x: x**2
    1-norm by passing lambda x: x
    and so on...
    """
    max_val = 0
    if mapfunc is None:
        for key, value in d.iteritems():
            if value > max_val:
                max_val = value
    else:
        for key, value in d.iteritems():
            max_val += mapfunc(value)

    if not max_val:
        return
    if endfunc:
        max_val = endfunc(max_val)

    for key, value in d.iteritems():
        d[key] = value / float(max_val)

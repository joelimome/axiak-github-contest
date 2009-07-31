#!/usr/bin/env python
"""
Given a list of watches,
generate a list of edges.

E.g.:
Input
  A:1
  B:2
  B:3
  C:1
  C:3
Output:
  A:C
  B:C
"""

import sys
from collections import defaultdict

def generate_edges(infile, outfile):
    doc_info = defaultdict(set)
    all_edges = set()

    for line in infile:
        user, doc = line.strip().split(':', 1)

        for other_user in doc_info[doc]:
            if user == other_user:
                continue
            edge = "%s:%s" % (user, other_user)
            if edge in all_edges:
                continue
            all_edges.add(edge)
            outfile.write("%s\n" % edge)
            outfile.flush()

        doc_info[doc].add(user)



if __name__ == '__main__':
    generate_edges(sys.stdin, sys.stdout)


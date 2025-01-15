#!/usr/bin/python3
# By Mattias Sjö, 2024
#
# Precomputes integer partitions for FORM
# This can be implemented in FORM but makes for much clumsier code
#
# First argument is the integer, second is a single letter specifying the kind of partition:
#   o - ordered: elements of each partition are ordered in descending order (this is the typical mathematical sense)
#   u - unordered: order of elements is free and different orders are treated as distinct partitions
#   f - fixed-length: unordered, zero elements are permitted, number of elements is fixed by third argument
#  (default is o)
# Generates the partitions in lexicographic order using a straightforward recursive algorithm with caching.
# Prints them to the screen and prepares a form file called ipart/(concatenation of arguments).hf

from datetime import datetime
from functools import cache
import sys

# These three implemetations are very similar, but a bit of code repitition is clearer
# than the various conditionals needed to fit all in one function

@cache
def upart(tot):
    if tot <= 0:
        return [[]]
    heads = [n+1 for n in range(tot)]
    return [[head] + tail
                for head in reversed(heads)
                for tail in upart(tot-head)]
@cache
def opart(tot, maxpart):
    if tot <= 0:
        return [[]]
    heads = [n+1 for n in range(min(maxpart,tot))]
    return [[head] + tail
                for head in reversed(heads)
                for tail in opart(tot-head, head)]
@cache
def fpart(tot, length):
    if length <= 0:
        return [[]]
    elif length == 1:
        return [[tot]]
    heads = [n for n in range(tot+1)]
    return [[head] + tail
                for head in reversed(heads)
                for tail in fpart(tot-head, length-1)]



def main(argv):

    if not argv:
        raise ValueError("Integer missing for ipart, should be called as 'ipart INTEGER [ouf]'")

    tot = int(argv[0])
    if tot < 0:
        raise ValueError(f"Integer must be non-negative ({tot} provided)")

    ord = argv[1] if len(argv) >= 2 else 'o'
    if ord not in 'uof':
        raise ValueError("Invalid option to ipart: should be o,u or f")

    if ord == 'f':
        if len(argv) < 3:
            raise ValueError("Missing length specification for fixed-length partition")
        length = int(argv[2])
        if length <= 0:
            raise ValueError(f"Length must be positive ({length} provided)")
    else:
        length = ''


    name = {'u': "Unordered", 'o': "Ordered", 'f': f"Length-{length}"}
    part = {'u': lambda n: upart(n), 'o': lambda n: opart(n,n), 'f': lambda n: fpart(n, length)}

    filename = f"ipart/{tot}{ord}{length}.hf"
    with open(filename, 'w') as out:
        print(f"{name[ord]} partitions of {tot}:")
        print(f"* Generated with options {argv} on {datetime.now().strftime('%c')}", file=out)
        for part in part[ord](tot):
            print(','.join(str(p) for p in part))
            print(f"   + ipart({','.join(str(p) for p in part)})", file=out)
        print(f"Wrote to {filename}")

if __name__ == '__main__':
    main(sys.argv[1:])

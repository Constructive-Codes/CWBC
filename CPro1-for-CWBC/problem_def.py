"""Constant Weight Binary Code definitions."""

PROB_NAME = "CWBC" # Short abbreviated name for problem
FULL_PROB_NAME = "Constant Weight Binary Code" # Full name for problem in Title Caps
PROB_DEF = """A Constant Weight Binary Code CWBC(n,w,d,s) is a set of at least s vectors of length n over the symbols {0,1}, such that (a) each vector has exactly w occurrences of the symbol 1, and (b) each pair of distinct vectors differ in at least d positions.  Given (n,w,d,s), we want to construct a CWBC(n,w,d,s).  For our purposes, n<64, 3<w<n, d is even and 4<=d<=20, and s<1000.  The code should be output as at least s lines with one vector per line, each vector a space-separated list of elements from {0,1}, with the list vectors satisfying constraints (a) and (b)."""
PARAMS = ["n","w","d","s"] # parameter names
DEV_INSTANCES = [[10,4,4,30],[13,5,6,18],[14,6,6,42],[15,7,6,69],[18,6,8,21],[20,7,8,80],[17,7,6,149],[17,8,6,166],[18,6,6,119],[19,6,6,155],[20,8,6,528],[21,8,6,698],[27,5,6,234],[22,9,8,252],[23,6,8,69],[24,6,8,70],[17,7,6,166],[17,8,6,184],[18,6,6,132],[19,6,6,172],[20,8,6,588],[21,8,6,775],[27,5,6,260],[22,9,8,280],[23,6,8,77],[24,6,8,78],[28,6,8,130],[21,9,10,27],[21,10,10,38],[22,8,10,24],[22,9,10,35],[24,8,10,38]]
OPEN_INSTANCES = [[17,7,6,167],[17,8,6,185],[18,6,6,133],[19,6,6,173],[20,8,6,589],[21,8,6,776],[27,5,6,261],[22,9,8,281],[23,6,8,78],[24,6,8,79],[28,6,8,131],[21,9,10,28],[21,10,10,39],[22,8,10,25],[22,9,10,36],[24,8,10,39]]

import utils

def getsol(soltext: str,parmlist: list[int]) -> list:
    """Extract solution array from raw program output soltext."""
    # utils.parsesoltext will do the work; just need to specify how many rows (how many lines of text) the solution array will have
    s = parmlist[3]
    array = utils.parsesoltext(soltext,s)
    return array

def v(soltext: str,parmlist: list[int]) -> bool:
    """Verify that the solution in raw program output soltext is a valid solution with parameters in parmlist."""
    n = parmlist[0]
    w = parmlist[1]
    d = parmlist[2]
    s = parmlist[3]
    array = getsol(soltext,parmlist)
    if len(array)!=s:
        return False
    for r in array:
        if len(r)!=n:
            return False
        for x in r:
            if x<0 or x>1:
                return False

    for i in range(0,s):
        wt = 0
        for k in range(0,n):
            if array[i][k]==1:
                wt += 1
        if wt!=w:
            return False
        for j in range(i+1,s):
            dist = 0
            for k in range(0,n):
                if array[i][k]!=array[j][k]:
                    dist += 1
            if dist<d:
                return False

    return True                 


#!/usr/bin/env python

# Code from Xiaoou

import sys

def main():
    with open(sys.argv[1], 'r') as f, open(sys.argv[2], 'w') as out:
        out.write('gene\tall')
        out.write('\n')
        for line in f:
            exp = line.rstrip().split("\t")
            out.write(exp[0])
            out.write('\t%f' % ts(exp[1:]))
            out.write('\n')


def ts(lst):
    exp = [float(x) for x in lst]
    n = len(exp)
    max_exp = max(exp)
    if(max_exp==0):
        return -0.1
    else:
        total = 0
        for i in exp:
            total += (max_exp - i) * 1.0 / max_exp
        total /= n - 1
        return total


if __name__ == '__main__':
    main()

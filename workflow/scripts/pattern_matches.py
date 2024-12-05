import sys


code = {'A':['A'],
        'C':['C'],
        'G':['G'],
        'T':['T'],
        'R':['A', 'G'],
        'Y':['C', 'T'],
        'S':['G', 'C'],
        'W':['A', 'T'],
        'K':['G', 'T'],
        'M':['A', 'C'],
        'B':['C', 'G', 'T'],
        'D':['A', 'G', 'T'],
        'H':['A', 'C', 'T'],
        'V':['A', 'C', 'G'],
        'N':['A', 'C', 'G', 'T']}


def matches(pattern):
    """generate all k-mers that match a pattern

    Args:
        pattern (str): IUPAC pattern

    Yields:
        str: k-mer
    """
    if len(pattern) == 0:
        yield ''
    else:
        for y in matches(pattern[1:]):
            for x in code[pattern[0]]:
                yield x+y


def reading(file):
    with open(file, "r") as f:
        for line in f:
            types, pattern, rate = line.split("\n")[0].split(" ")
            if types == "deletion" or types == "insertion":
                res = matches(pattern)
                for ogpat in list(res):
                    print(types, ogpat, rate)
            else:
                print(types, pattern, rate)
                

if __name__ == '__main__':
    file = sys.argv[1]
    reading(file)

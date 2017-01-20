#! /usr/bin/env python3
"""subcellular location"""


def get_gene2sublocation(filepath, prefix, genecol=2, sublocationcol=4):
    """
    get gene-subcellular location assos, file content delimiter is '\t'
    :param filepath: subcellular loaction file path, better be from
    the knowledge channel of compartments website
    :param genecol: number of column indicates gene
    :param sublocationcol: number of column indicates sublocation name
    :param prefix: sublocation name prefix
    :return: a dict object, string-set<string>
    """
    gcol = genecol - 1
    scol = sublocationcol - 1
    sublocprefixes = set()
    for pf in prefix:
        sublocprefixes.add(pf)
    gene2subloc = {}
    # printsub = set()
    with open(filepath, mode='r') as rf:
        for line in rf:
            words = line.strip().split("\t")
            gene = words[gcol].strip()
            subloc = words[scol].strip()
            for sublocprefix in sublocprefixes:
                if subloc.lower().startswith(sublocprefix.lower()):
                    if gene not in gene2subloc.keys():
                        gene2subloc[gene] = set()
                    gene2subloc[gene].add(sublocprefix)
                    # printsub.add(subloc)
                    break
    # import pprint
    # pprint.pprint(printsub)
    return gene2subloc


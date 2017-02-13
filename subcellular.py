#! /usr/bin/env python3
"""subcellular location"""
import numpy as np
import files


def weighting_ppi(ppiarray, pnames, p2subloc):
    """
    references:Tang X, Hu X, Yang X, et al. A algorithm for identifying disease genes by incorporating the subcellular
     localization information into the protein-protein interaction networks[C]//Bioinformatics and Biomedicine (BIBM),
     2016 IEEE International Conference on. IEEE, 2016: 308-311.
    :param ppiarray: numpy array
    :param pnames: a list
    :param p2subloc: a dict
    :return: weighted ppi array
    """
    subloc2p = files.invert_dict(p2subloc)
    maxsubloc = 0
    for subloc in subloc2p.keys():
        if len(subloc2p[subloc]) > maxsubloc:
            maxsubloc = len(subloc2p[subloc])
    subloc2score = {}
    for subloc in subloc2p.keys():
        subloc2score[subloc] = len(subloc2p[subloc]) / maxsubloc
    minsublocscore = 1
    for subloc in subloc2score.keys():
        if minsublocscore > subloc2score[subloc]:
            minsublocscore = subloc2score[subloc]

    wppiarray = np.zeros(ppiarray.shape)
    (rlen, clen) = ppiarray.shape
    for i in range(0, rlen):
        for j in range(0, clen):
            if ppiarray[i, j] == 1.0:
                p1, p2 = pnames[i], pnames[j]
                commonsubloc = set()
                if p1 in p2subloc.keys() and p2 in p2subloc.keys():
                    commonsubloc = p2subloc[p1].intersection(p2subloc[p2])
                if len(commonsubloc) == 0:
                    wppiarray[i, j] = minsublocscore
                else:
                    for csubloc in commonsubloc:
                        if wppiarray[i, j] < subloc2score[csubloc]:
                            wppiarray[i, j] = subloc2score[csubloc]
    return wppiarray


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


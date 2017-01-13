#! /usr/bin/env python3
"""similarity normalize"""


def norm_k_nearest_neighbor(sim, k=5):
    """

    :param sim: sim (must not have duplite sim pairs),
    dict, key-value: string-dict<string-value> ({entity1: {entity2: sim, }, }
    :param k: k
    :return: sim,
    dict, key-value: string-dict<string-value> ({entity1: {entity2: sim, }, }
    """
    simlist = {}
    for e1 in sim.keys():
        if e1 not in simlist.keys():
            simlist[e1] = []
        for e2 in sim[e1]:
            if e2 not in simlist.keys():
                simlist[e2] = []
            if e1 != e2:
                simlist[e1].append((e2, sim[e1][e2]))
                simlist[e2].append((e1, sim[e1][e2]))
    ressim = {}
    for e1 in simlist.keys():
        # simlist[e1] = sorted(simlist[e1], key=lambda x: x[0])
        simlist[e1] = sorted(simlist[e1], key=lambda x: x[1], reverse=True)
        count = k
        if len(simlist[e1]) < k:
            count = len(simlist[e1])
        for i in range(0, count):
            e2 = simlist[e1][i][0]
            sval = simlist[e1][i][1]
            if not ((e1 in ressim.keys() and e2 in ressim[e1].keys()) or
                    (e2 in ressim.keys() and e1 in ressim[e2].keys())):
                if e1 not in ressim.keys():
                    ressim[e1] = {}
                ressim[e1][e2] = sval
        for i in range(count, len(simlist[e1])):
            e2 = simlist[e1][i][0]
            sval = simlist[e1][i][1]
            if sval == simlist[e1][count-1][1]:
                if not ((e1 in ressim.keys() and e2 in ressim[e1].keys()) or
                        (e2 in ressim.keys() and e1 in ressim[e2].keys())):
                    if e1 not in ressim.keys():
                        ressim[e1] = {}
                    ressim[e1][e2] = sval
            else:
                break
    return ressim

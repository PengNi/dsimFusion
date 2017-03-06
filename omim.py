#! /usr/bin/env python3
"""OMIM file parser"""
import re


class OmimEntries:
    def __init__(self):
        self._entries = {}

    def readinfo_file(self, filepath):
        """
        'mim2gene.txt'
        :param filepath:
        :return:
        """
        with open(filepath, mode='r') as rf:
            for line in rf:
                if not line.startswith("#"):
                    mimtemp = OmimEntry()
                    words = line.strip('\n').split('\t')
                    mimnumer = words[0].strip()
                    mimtype = words[1].strip()
                    mimentrezid = words[2].strip()
                    mimsymbol = words[3].strip()
                    mimensembles = words[4].strip()

                    mimtemp.setmimnumber(mimnumer)
                    mimtemp.settype(mimtype)
                    if mimentrezid != '':
                        mimtemp.setentrezid(mimentrezid)
                    if mimsymbol != '':
                        mimtemp.setsymbol(mimsymbol)
                    if mimensembles != '':
                        ensembles = set(mimensembles.split(','))
                        mimtemp.setensembleids(ensembles)

                    if mimnumer not in self._entries.keys():
                        self._entries[mimnumer] = mimtemp
        print('mims:', len(self._entries))

    def readmorbidmap(self, filepath, genetype):
        """
        read morbidmap.txt (must after reading 'mim2gene.txt')
        :param filepath: path of 'morbidmap.txt'
        :param genetype: 'omim' or 'entrezid' or 'symbol'
        :return: phenotype2gene assos
        """
        miminfos = self._entries
        rawp2g = {}
        phepattern = re.compile(r'\d\d\d\d\d\d')
        type_pheno = re.compile(r'phenotype')
        type_gene = re.compile(r'gene')
        with open(filepath, mode='r') as rf:
            for line in rf:
                if not line.startswith('#'):
                    words = line.strip('\n').split('\t')
                    phetext = words[0].strip()
                    mimnumber = words[2].strip()

                    ptemp = phepattern.search(phetext)
                    if ptemp:
                        phenotype = ptemp.group()
                    else:
                        phenotype = mimnumber
                    if type_pheno.search(miminfos[phenotype].gettype()) and \
                            not type_gene.search(miminfos[phenotype].gettype()):
                        # if miminfos[phenotype].gettype() != 'phenotype':
                        #     print(phenotype, miminfos[phenotype].gettype(), mimnumber)
                        if phenotype not in rawp2g.keys():
                            rawp2g[phenotype] = set()
                        rawp2g[phenotype].add(mimnumber)
                    # else:
                    #     print(phenotype, miminfos[phenotype].gettype())
                    #     pass
        p2g = {}
        if genetype == 'omim':
            for p in rawp2g.keys():
                for g in rawp2g[p]:
                    if type_gene.search(miminfos[g].gettype()):
                        if p not in p2g.keys():
                            p2g[p] = set()
                        p2g[p].add(g)
            # ---decide to add the following code---
            ps = list(p2g.keys())
            for p in ps:
                p2g[p].discard(p)
                if len(p2g[p]) == 0:
                    del p2g[p]
                    # print(p)
            # ---------------------------------------
        elif genetype == 'entrezid':
            for p in rawp2g.keys():
                for g in rawp2g[p]:
                    if miminfos[g].getentrezid() != '':
                        if p not in p2g.keys():
                            p2g[p] = set()
                        p2g[p].add(miminfos[g].getentrezid())
        elif genetype == 'symbol':
            for p in rawp2g.keys():
                for g in rawp2g[p]:
                    if miminfos[g].getsymbol() != '':
                        if p not in p2g.keys():
                            p2g[p] = set()
                        p2g[p].add(miminfos[g].getsymbol())
        else:
            print('readmorbidmap(): wrong argument!')
        return p2g


class OmimEntry:

    def __init__(self):
        self._mimnumber = ""
        self._type = ""
        self._entrezid = ""
        self._symbol = ""
        self._ensemblids = set()

    def getmimnuber(self):
        return self._mimnumber

    def setmimnumber(self, mimnubmer):
        self._mimnumber = mimnubmer

    def gettype(self):
        return self._type

    def settype(self, mtype):
        self._type = mtype

    def getentrezid(self):
        return self._entrezid

    def setentrezid(self, entrezid):
        self._entrezid = entrezid

    def getsymbol(self):
        return self._symbol

    def setsymbol(self, symbol):
        self._symbol = symbol

    def getensembleids(self):
        return self._ensemblids

    def setensembleids(self, ensembleids):
        self._ensemblids = ensembleids

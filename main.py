#! /usr/bin/env python3
"""main"""
import omim
import files
import similarity_normalize
import subcellular
import mapping
import similarity_module
import experiments


# ------omim pheno 2 geno------
def omim_p2g():
    mims = omim.OmimEntries()
    mims.readinfo_file("data/omim/mim2gene.txt")
    pheno2geno = mims.readmorbidmap('data/omim/morbidmap.txt', 'symbol')
    files.stat_assos(pheno2geno)
    files.write_assos(pheno2geno, 'data/omim/pheno2genesymbol.tsv')


# ------normalize similarity------
def mimminer_analysis():
    mimminersim = files.read_simmatrix('D:/bioinformatics/database/MimMiner/MimMiner_disease_disease_similarity.tsv',
                                       True, False)
    files.stat_sims(mimminersim)
    # knnsim = similarity_normalize.norm_k_nearest_neighbor(mimminersim)
    # files.stat_sims(knnsim)

    d170order = files.read_one_col('D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\'
                                   'PhenotypeID_170.tsv', 1)
    retrivesims = files.retrive_sims(mimminersim, d170order)
    files.stat_sims(retrivesims)
    # knnsim = similarity_normalize.norm_k_nearest_neighbor(retrivesims)
    # files.write_simmatrix(retrivesims, 'data/sims/miminer170_matrix.tsv', True, d170order)


def normalize_similarity():
    sim = files.read_simmatrix("D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\"
                               "dsim_snf_spavgnrwrknn_crstar170.tsv")
    files.stat_sims(sim)
    sim_knn = similarity_normalize.norm_k_nearest_neighbor(sim)
    files.stat_sims(sim_knn)
    d170 = files.read_one_col("D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\"
                              "PhenotypeID_170.tsv", 1)
    # files.write_simmatrix(sim_knn,
    #                       "D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\"
    #                       "similarity_snfknnspavgnrwrknn_digeomim170_originalppimagger_matrix.tsv",
    #                       True, d170, '\t', False, False)
    pass


def normalize_similarity2():
    sims = files.read_sims("D:\\Documents\\workspace\\pyworkspace\\dsimModuleTheory\\"
                           "outputs\\similarity_hpole_birwomim5080_triplet.tsv")
    files.stat_sims(sims)
    knnsims = similarity_normalize.norm_k_nearest_neighbor(sims)
    dnames = files.read_one_col("data/birw_xie/BiRW_phenotype_annotation.txt", 2)
    # dnames = files.read_one_col("D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\"
    #                             "PhenotypeID_170.tsv", 1)
    files.write_simmatrix(sims, "data/sims/similarity_hpole_birwomim5080_matrix.tsv",
                          True, dnames, '\t', False, False)


# ------similarity calculation------
def simcal_module():
    g = similarity_module.read_interactome("data/birw_xie/BiRW_ppi_network_hprdnumber.txt", False, False)
    print(len(g.vs), len(g.es))
    d2g = files.read_assos("data/birw_xie/BiRW_pheno2genonumber_omim2007.txt")
    files.stat_assos(d2g)
    sim = similarity_module.similarity_cal_spavgn(d2g, g)
    files.write_sims(sim, "data/sims/similarity_spavgn_birwdgomim2007_birwhprd.tsv")

    knnsims = similarity_normalize.norm_k_nearest_neighbor(sim)
    dnames = files.read_one_col("data/birw_xie/BiRW_phenotype_annotation.txt", 2)
    # files.write_simmatrix(knnsims, "data/sims/similarity_spavgnknn_birwdgomim2007_birwhprd_5080mat.tsv",
    #                       True, dnames, '\t', False, False)


def simcal_module_multitimes():
    g = similarity_module.read_interactome("data/ppi/tissuespec_magger/original_ppi.txt", False, False)
    print(len(g.vs), len(g.es))

    d2gdirectory = 'D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\data\\'
    d2gfile = d2gdirectory + 'pheno2geno.txt'
    d2g = files.read_assos(d2gfile, True)
    files.stat_assos(d2g)
    for d in d2g.keys():
        print(d, len(d2g[d]))

    d170order = files.read_one_col(d2gdirectory + 'PhenotypeID_170.tsv', 1)

    d2g_tuple = []
    with open(d2gfile, mode='r') as rf:
        next(rf)
        for line in rf:
            words = line.strip().split('\t')
            d2g_tuple.append((words[0], words[1]))
    print(len(d2g_tuple))
    # with open(d2gdirectory + 'pheno2geno_multimark_v2.txt', mode='w') as wf:
    #     for d2gt in d2g_tuple:
    #         wf.write(str(d2gt[0]) + '\t' + str(d2gt[1]) + '\t')
    #         if len(d2g[d2gt[0]]) > 4:
    #             wf.write('1\t1\t1\t1\t1\n')
    #         elif len(d2g[d2gt[0]]) > 3:
    #             wf.write('1\t1\t1\t1\t0\n')
    #         elif len(d2g[d2gt[0]]) > 2:
    #             wf.write('1\t1\t1\t0\t0\n')
    #         elif len(d2g[d2gt[0]]) > 1:
    #             wf.write('1\t1\t0\t0\t0\n')
    #         elif len(d2g[d2gt[0]]) > 0:
    #             wf.write('1\t0\t0\t0\t0\n')
    #         else:
    #             wf.write('0\t0\t0\t0\t0\n')

    # transformdistance = True
    # dgs = set()
    # for d in d2g.keys():
    #     dgs |= set(d2g[d])
    # print("disease genes num:", len(dgs))
    # genegenesim = similarity_module.sim_gene2gene_shortestpath(dgs, g, transformdistance)
    # for i in range(0, len(d2g_tuple)):
    #     dgtuple = d2g_tuple[i]
    #     filenames = d2gdirectory + 'data\\modulesim\\similarity170_matrix' + str(i+1) + '.txt'
    #     d2g_new = {}
    #     for d in d2g.keys():
    #         d2g_new[d] = set()
    #         d2g_new[d].update(d2g[d])
    #     d2g_new[dgtuple[0]].remove(dgtuple[1])
    #     files.stat_assos(d2g_new)
    #     simtemp = similarity_module.similarity_cal_spavgn_multitimescore(d2g_new, g, genegenesim, True)
    #     files.write_simmatrix(simtemp, filenames, True, d170order, '\t', False, False)


def simcal_netsim_multitimes():
    d2gdirectory = 'D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\data\\'
    d2gfile = d2gdirectory + 'pheno2geno.txt'
    d2g = files.read_assos(d2gfile, True)
    files.stat_assos(d2g)
    d170order = files.read_one_col(d2gdirectory + 'PhenotypeID_170.tsv', 1)
    d2g_tuple = []
    with open(d2gfile, mode='r') as rf:
        next(rf)
        for line in rf:
            words = line.strip().split('\t')
            d2g_tuple.append((words[0], words[1]))
    print(len(d2g_tuple))

    sim_g2d_filepre = "D:\\Documents\\workspace\\rworkspace\\rwr\\data\\netsim\\similarity_g2d_matrix"
    for i in range(0, len(d2g_tuple)):
        dgtuple = d2g_tuple[i]

        d2g_new = {}
        for d in d2g.keys():
            d2g_new[d] = set()
            d2g_new[d].update(d2g[d])
        d2g_new[dgtuple[0]].remove(dgtuple[1])
        if len(d2g_new[dgtuple[0]]) == 0:
            del d2g_new[dgtuple[0]]
        files.stat_assos(d2g_new)
        sim_g2d = files.read_simmatrix(sim_g2d_filepre + str(i+1) + '.txt')
        sim = experiments.sim_geneset2geneset_rwr(d2g_new, sim_g2d)
        filenames = d2gdirectory + 'netsim\\similarity170_matrix' + str(i + 1) + '.txt'
        files.write_simmatrix(sim, filenames, True, d170order, '\t', False, False)


# ---------analysis-----------------
def compare_2_version_dgassos():
    d170order = files.read_one_col('D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\data\\'
                                   'PhenotypeID_170.tsv', 1)
    d2g1 = files.read_assos('data/omim/pheno2geneentrez.tsv')
    d2g2 = files.read_assos('D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\data\\'
                            'pheno2geno.txt', True)
    for d in d170order:
        if d in d2g1.keys():
            print(d, d2g1[d], d2g2[d], sep='\t\t\t')
        else:
            print(d, '--', d2g2[d], sep='\t\t\t')
    pass


def ana_dgassos():
    d2g = files.read_assos('D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\'
                           'pheno2geno.txt', True)
    files.stat_assos(d2g)
    for d in d2g.keys():
        print(d, d2g[d])
    print()
    g2d = files.invert_dict(d2g)
    for g in g2d.keys():
        print(g, g2d[g])


# ---------access data from magger's paper---------
# Magger O, Waldman Y Y, Ruppin E, et al. Enhancing the prioritization of disease-causing genes through tissue
# specific protein interaction networks[J]. PLoS Comput Biol, 2012, 8(9): e1002690.
def get_original_ppi_magger():
    directory = 'D:\\bioinformatics\\paper\\disease\\disease-gene-prediction\\Enhancing the ' \
                'Prioritization of Disease-Causing Genes through Tissue Specific Protein Interaction Networks\\'
    fileloc = directory + 'data_s2.txt'
    mappingfile = directory + 'id2index.txt'
    index2id = files.read_mappings(mappingfile, False, '\t', 2, 1)
    network = {}
    with open(fileloc, mode='r') as rf:
        for line in rf:
            if not line.startswith('#'):
                words = line.strip().split(' ')
                n1, n2 = words[0], words[1]
                if not ((n1 in network.keys() and n2 in network[n1]) or
                        (n2 in network.keys() and n1 in network[n2])):
                    if n1 not in network.keys():
                        network[n1] = set()
                    network[n1].add(n2)
    files.stat_network(network)
    orinet = {}
    for index1 in network.keys():
        id1 = index2id[index1]
        orinet[id1] = set()
        for index2 in network[index1]:
            orinet[id1].add(index2id[index2])
    files.stat_network(orinet)
    # write_assos(orinet, directory + 'original_ppi.txt')


# ------subcellular location------
def subcellularloc_analysis():
    prefix = files.read_one_col('data/compartments/localization_class.txt', 1)
    gene2subloc = subcellular.get_gene2sublocation("data/compartments/human_compartment_knowledge_full.tsv",
                                                   prefix)
    files.stat_assos(gene2subloc)
    # files.write_assos(gene2subloc, "data/compartments/genesymbol2sublocation_human_knowledge.tsv")


def get_entrezid2subcellular():
    genesymbol2subloc = files.read_assos("data/compartments/genesymbol2sublocation_human_knowledge.tsv")
    files.stat_assos(genesymbol2subloc)

    genesymbols = list(genesymbol2subloc.keys())
    genesymbol2entrezid = mapping.geneida2geneidb('symbol', 'entrezgene', genesymbols)
    files.stat_assos(genesymbol2entrezid)
    entrezid2subloc = {}
    for symbol in genesymbol2subloc.keys():
        if symbol in genesymbol2entrezid.keys():
            for e in genesymbol2entrezid[symbol]:
                entrezid2subloc[e] = set()
                entrezid2subloc[e].update(genesymbol2subloc[symbol])
    files.stat_assos(entrezid2subloc)
    # files.write_assos(entrezid2subloc, "data/compartments/geneentrezid2sublocation_human_knowledge_mygenemapping.tsv")


# Tang X, Hu X, Yang X, et al. A algorithm for identifying disease genes by incorporating the subcellular
# localization information into the protein-protein interaction networks[C]//Bioinformatics and Biomedicine (BIBM),
# 2016 IEEE International Conference on. IEEE, 2016: 308-311.
def weighting_ppi_by_subcellularloc():
    p2subloc = files.read_assos("data/birw_xie/CRstar_ppi_gene2subloc.tsv")
    files.stat_assos(p2subloc)
    ppiarray = files.read_symmatrix2array("data/birw_xie/CRstar_ppi_network.txt")
    pnames = files.read_one_col("data/birw_xie/CRstar_ppi_annotation.txt", 1, False)
    print(len(pnames))
    wppiarray = subcellular.weighting_ppi(ppiarray, pnames, p2subloc)
    # files.write_symarray2file(wppiarray, "data/birw_xie/CRstar_ppi_network_sublocweighted.txt")


# ------BiRW_Xie------
# Xie M Q, Xu Y J, Zhang Y G, et al. Network-based phenome-genome association prediction by bi-random walk[J].
# PloS one, 2015, 10(5): e0125138.
def get_gene2subcellular():
    genenumber2symbol = files.read_mappings("data/birw_xie/BiRW_ppi_annotation.txt", True)
    files.stat_maps(genenumber2symbol)
    genesymbol2subcellular = files.read_assos("data/compartments/genesymbol2sublocation_human_knowledge.tsv")
    files.stat_assos(genesymbol2subcellular)

    genenumber2subloc = {}
    for gnum in genenumber2symbol.keys():
        if genenumber2symbol[gnum] in genesymbol2subcellular.keys():
            genenumber2subloc[gnum] = set()
            genenumber2subloc[gnum].update(genesymbol2subcellular[genenumber2symbol[gnum]])
    files.stat_assos(genenumber2subloc)
    # files.write_assos(genenumber2subloc, "data/birw_xie/BiRW_ppi_gene2subloc.tsv")


def filter_gene2subcellular():
    entrezids = files.read_one_col("data/birw_xie/CRstar_ppi_annotation.txt", 1)
    entrezid2subloc = files.read_assos("data/compartments/geneentrezid2sublocation_human_knowledge_mygenemapping.tsv")
    files.stat_assos(entrezid2subloc)
    entrezid2subloc_filter = {}
    for e in entrezids:
        if e in entrezid2subloc.keys():
            entrezid2subloc_filter[e] = set()
            entrezid2subloc_filter[e].update(entrezid2subloc[e])
    files.stat_assos(entrezid2subloc_filter)
    files.write_assos(entrezid2subloc_filter, "data/birw_xie/CRstar_ppi_gene2subloc.tsv")


if __name__ == '__main__':
    omim_p2g()
    pass

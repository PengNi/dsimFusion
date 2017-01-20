#! /usr/bin/env python3
"""main"""
import omim
import files
import similarity_normalize


# ------omim pheno 2 geno------
def omim_p2g():
    mims = omim.OmimEntries()
    mims.readinfo_file("data/omim/mim2gene.txt")
    pheno2geno = mims.readmorbidmap('data/omim/morbidmap.txt', 'entrezid')
    files.stat_assos(pheno2geno)
    files.write_assos(pheno2geno, 'data/omim/pheno2geneentrez.tsv')


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
    files.write_simmatrix(retrivesims, 'data/sims/miminer170_matrix.tsv', True, d170order)


def normalize_similarity():
    sim = files.read_simmatrix("D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\"
                               "dsim_snf_mimminerspavgn_crstar170.tsv")
    files.stat_sims(sim)
    sim_knn = similarity_normalize.norm_k_nearest_neighbor(sim)
    files.stat_sims(sim_knn)
    d170 = files.read_one_col("D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\"
                              "PhenotypeID_170.tsv", 1)
    files.write_simmatrix(sim_knn,
                          "D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\"
                          "similarity_snfmimminerspavgn_digeomim170_originalppimagger_matrix.tsv", True, d170)
    pass


# ---------analysis-----------------
def compare_2_version_dgassos():
    d170order = files.read_one_col('D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\'
                                   'PhenotypeID_170.tsv', 1)
    d2g1 = files.read_assos('data/omim/pheno2geneentrez.tsv')
    d2g2 = files.read_assos('D:\\Documents\\workspace\\matlabworkspace\\NoNCRstar-master\\CRstar\\'
                            'pheno2geno.txt', True)
    for d in d170order:
        if d in d2g1.keys():
            print(d, d2g1[d], d2g2[d], sep='\t\t')
        else:
            print(d, '--', d2g2[d], sep='\t\t')
    pass


# ---------access data--------------
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


if __name__ == '__main__':
    pass

# -*- coding: utf-8 -*-
# @Time    : 11/19/20 5:49 PM
# @Author  : Kay
# @Email   : kahy.shen@gmail.com
# @File    : randomNetwork.py
# @Desc    :

import networkx as nx
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
import verify
import Statistics
import numpy as np
# import matplotlib.pyplot as plt
from scipy.stats import exponnorm
import pandas as pd

def func(x, alpha, c):
    return c*pow(x, -alpha)

def funcXSquare(x, alpha, beta, C):
    return C * pow(beta*pow(x, -alpha), np.log(x))


def readDegree(filename):
    with open(filename, "r") as reader:
        data = reader.readlines()
        data = [int(x) for x in data[0].split("\t")]
        # print(data)
    return data

def C(N, m):
    product = 1
    deviation = 1
    for i in range(m):
        product = product * (N-i)
        deviation = deviation * (i+1)
    return product / deviation


def expectedS(N, m, p):
    prob_line = C(N, 2) * p
    nearbyNode = C(N-2, m) * pow(p, 2*m)
    otherNode = pow(1 - p*p, N-2-m)
    return prob_line * nearbyNode * otherNode


def expectedT(N, m, p):
    prob_tri = C(N, 3) * pow(p, 3)
    nearbyNode = C(N-3, m) * pow((3 * pow(p, 2) * (1-p)), m)
    otherNode = pow(1 - 3*p*p + 2*p*p*p, N-3-m)
    return prob_tri * nearbyNode * otherNode

def adjTriDistribution(simplicesDict):
    data = []
    distribution = {}
    for tri in simplicesDict:
        data.append(len(simplicesDict[tri]))
        if len(simplicesDict[tri]) in distribution:
            distribution[len(simplicesDict[tri])] += 1
        else:
            distribution[len(simplicesDict[tri])] = 1
    params = exponnorm.fit(data, 10)
    print(params)
    # mu, std = norm.fit(data)
    return distribution, params


# def visualizeDistribution(distributionList, legendsList, titlename):
#     # maxX = max([max(dis.keys())for dis in distributionList])
#     # x = [for i in range(1, maxX+1)]
#     plt.figure(figsize=(8, 8.5))
#     colors = ['g', 'b', 'r', 'y']
#     marks = ['+', '-', '.', '*']
#     for i in range(len(distributionList)):
#
#         x = sorted(distributionList[i].keys())
#         y = [distributionList[i][k] / sum(distributionList[i].values()) for k in x]
#         # popt, pcov = curve_fit(func, x, y)
#         # alpha = round(popt[0], 2)
#         # c = round(popt[1], 2)
#         plt.plot(x, y, '{}{}-'.format(colors[i % 4], marks[i // 4]), label='{}'.format(legendsList[i]))
#     plt.legend()
#     # plt.savefig("randon_{}.png".format(titlename))
#     plt.savefig("powerLaw_{}.png".format(titlename))
#     plt.close()
def visualizeT(distribution, legends, titlename, filename):
    # maxX = max([max(dis.keys())for dis in distributionList])
    # x = [for i in range(1, maxX+1)]
    plt.figure(figsize=(8, 8.5))
    x = sorted(distribution.keys())
    y = [distribution[k] / sum(distribution.values()) for k in x]
    plt.plot(x, y, 'r+-', label='realTDistribution_No.∆_{}'.format(legends[0]))
    plt.plot(x, exponnorm.pdf(x, *legends[1]), 'g--', label='approximated exponorm with '
                                                        'lambda: {}, mu: {}, sigma:{}'.
             format(round(1 / legends[1][0] / legends[1][2], 2), round(legends[1][1], 2), round(legends[2], 2)))
    plt.legend()
    # plt.savefig("randon_{}.png".format(titlename))
    plt.savefig("{}/powerLaw_{}.png".format(filename, titlename))
    plt.close()

def visualizeS(distribution, legends, titlename, filename):
    # maxX = max([max(dis.keys())for dis in distributionList])
    # x = [for i in range(1, maxX+1)]
    plt.figure(figsize=(8, 8.5))
    colors = ['g', 'b', 'r', 'y']
    marks = ['+', '-', '.', '*']
    x = sorted(distribution.keys())
    y = [distribution[k] / sum(distribution.values()) for k in x]
    popt, pcov = curve_fit(funcXSquare, x, y)
    plt.loglog(x, y, 'r+-', label='{}'.format(legends))
    plt.loglog(x, funcXSquare(x, *popt), 'g--',
             label='fit: a=%5.2f, b=%5.2f, c=%5.2f' % tuple(popt))
    plt.legend()
    # plt.savefig("randon_{}.png".format(titlename))
    plt.savefig("{}/powerLaw_{}.png".format(filename, titlename))
    plt.close()
    return tuple(popt)

def func(x, alpha, c):
    return c * pow(x, -alpha)


def randomPowerLaw(n, m, p):
    G = nx.powerlaw_cluster_graph(n, m, p)
    print("G ready")
    simplex = verify.main(G, "/mnt/simplices/simplicial_complex/data/enron/poweLawRandom/", "{}_{}_{}".format(n, m, p))
    s_simplices = simplex.s_simplices
    t_simplices = simplex.t_simplices
    disT, params = adjTriDistribution(t_simplices)
    legends = [len(simplex.triangles), params]
    print(legends)
    L = [round(1 / legends[1][0] / legends[1][2], 2), round(legends[1][1], 2), round(legends[2], 2)]
    visualizeT(disT, legends, "N_{}_m_{}_p_{}_T".format(n, m, p), "data/enron/poweLawRandom/{}_{}_{}".format(n, m, p))
    disS, invalidParams = adjTriDistribution(s_simplices)
    legend = "No. ∆_{}".format(len(simplex.triangles))
    s_para = visualizeS(disS, legend, "N_{}_m_{}_p_{}_S".format(n, m, p), "data/enron/poweLawRandom/{}_{}_{}".format(n, m, p))
    return L, s_para

def originalenron(i):
    with open("data/enron/enron_edges.txt", "r") as reader:
        data = reader.readlines()
        nodes = []
        lines = [d.split("\n")[0].split("\t") for d in data]
        L = [[int(l[0]), int(l[1])] for l in lines]
        reader.close()
        for l in lines:
            nodes.append(l[0])
            nodes.append(l[1])
    G = nx.Graph()
    G.add_edges_from([(e[0], e[1]) for e in L])
    print("G ready")
    simplex = verify.main(G, "/mnt/simplices/simplicial_complex/data/enron/poweLawRandom/", "{}".format(i))
    s_simplices = simplex.s_simplices
    t_simplices = simplex.t_simplices
    disT, params = adjTriDistribution(t_simplices)
    legends = [len(simplex.triangles), params]
    print(legends)
    visualizeT(disT, legends, "T", "data/enron/poweLawRandom/{}".format(i))
    # disS, invalidParams = adjTriDistribution(s_simplices)
    # legend = "No. ∆_{}".format(len(simplex.triangles))
    # visualizeS(disS, legend, "S", "data/enron/poweLawRandom/{}".format(i))

if __name__ == '__main__':
    # folderName = "data/enron/random/"
    # # degree = readDegree("degree.txt")
    # random = [[4, 0.2], [4, 0.3]]
    # for p in random:
    #     G = nx.powerlaw_cluster_graph(36265, p[0], p[1])
    #     G.remove_edges_from(G.selfloop_edges())
    # # for i in range(1):
    # #     G = nx.configuration_model(degree, create_using=nx.Graph())
    #     verify.main(G, folderName, "{}-{}".format(p[0], p[1]))
    #     Statistics.main("data/enron/random/{}/".format("{}-{}".format(p[0], p[1])))
    #
    #
    # # for i in range(5):
    # #     # m = 4
    # #     # p = 0.15
    # #     # gamma = 1.883
    # #     # G = nx.random_powerlaw_tree(36265, gamma=gamma)
    # #     # G = nx.powerlaw_cluster_graph(36265, m, p)
    # #     G = nx.configuration_model(degree, create_using=nx.Graph)
    # #     clusterCoef = nx.clustering(G)
    # #     Coef = [clusterCoef[c] for c in clusterCoef]
    # #     Coef = sum(Coef) / len(Coef)
    # #     print(Coef)
    # #     degreeDis = nx.degree(G)
    # #     Dis = {}
    # #     for d in degreeDis:
    # #         if d[1] not in Dis:
    # #             Dis[d[1]] = 1
    # #         else:
    # #             Dis[d[1]] += 1
    # #     x = sorted(list(Dis.keys()))
    # #     S = sum([Dis[x] for x in Dis])
    # #     y = [Dis[k] / S for k in x]
    # #
    # #     plt.plot(x, y, 'b-')
    # #     popt, pcov = curve_fit(func, x, y)
    # #
    # #     y2 = [func(i, popt[0], popt[1]) for i in x]
    # #     plt.loglog(x, y, 'g+-')
    # #     plt.loglog(x, y2, 'r--', label="alpha={}, c={}, clusCoef ={}".format(round(popt[0], 4), round(popt[1], 4), round(Coef, 4)))
    # #     plt.legend()
    # #     # plt.savefig('/Users/keshen/PycharmProjects/sk/ISI/simplicial_complex/output/RandomNetWork/enronRand/powerLawCluster/powerLaw-m{}-p{}.png'.format(m, p))
    # #     plt.savefig('/Users/keshen/PycharmProjects/sk/ISI/simplicial_complex/output/RandomNetWork/enronRand/specifyDegree/random{}.png'.format(i))
    # #
    # #     plt.show()
    # # # plt.savefig('degree_distribution.png')
    #####
    """random non powerlaw network"""

    # N = 1000
    # Dis = []
    # Legends = []
    # for p in [0.01, 0.03, 0.05, 0.07]:
    # # p = 0.05
    # # m = 4
    #     randomG = nx.binomial_graph(N, p)
    #     simplex = verify.main(randomG, "/Users/keshen/PycharmProjects/sk/ISI/simplicial_complex/data/Random/", "test")
    #
    #     # tri0adj = [1 for tri in simplex.triadj if len(simplex.triadj[tri]) == 1]
    #     s_simplices = simplex.s_simplices
    #     t_simplices = simplex.t_simplices
    #     dis, mu, std = adjTriDistribution(t_simplices)
    #     Dis.append(dis)
    #     Legends.append("p={},mu={},std={}".format(p, round(mu, 2), round(std, 2)))
    # visualizeDistribution(Dis, Legends, "N_{}".format(N))



    # expectS = expectedS(N, m, p)
    # T = expectedT(N, m-1, p)
    # expectT = N * (N - 1) / 2 * p * (N - 2) * (N - 3) / 2 * 2 * p*p * p*p *(1-p) * pow(1 - 3*p*p + 2*p*p*p, N - 2 - 2)
    # print(T, expectT)
    # print("except T: {}, except S: {}".format(T, expectS))

    # s_2adj = [1 for tri in s_simplices if len(s_simplices[tri]) == m]
    # t_1adj = [1 for tri in t_simplices if len(t_simplices[tri]) == m-1]
    # print(sum(s_2adj), sum(t_1adj))
    # adjTriDistribution(t_simplices)

    """random powerlaw network"""
    L = []
    for n in range(100, 201, 100):
        for m in range(2, 3):
            for p in np.arange(0.0, 0.1, 0.1):
                print(m,n,p)
                l, s_para = randomPowerLaw(n, m, p)
                l.append(s_para[0])
                l.append(s_para[1])
                l.append(s_para[2])
                l.insert(0, "{}_{}_{}".format(n, m, p))
                print(l)
                L.append(l)
    import pandas as pd
    name = ["setting", "lambda", "mu", "sigma", "a", "b", "c"]
    test = pd.DataFrame(columns=name, data=L)
    # test.to_csv('para2para_.csv')




    # N = [3000]
    # M = [5, 10, 15]
    # P = [0, 0.01, 0.03, 0.05]
    # for n in N:
    #     Dis = []
    #     Legends = []
    #     for m in M:
    #         for p in P:
    #             G = nx.powerlaw_cluster_graph(n, m, p)
    #             # simplex = verify.main(G, "/mnt/simplices/simplicial_complex/data/enron/random/", "powerLaw")
    #             simplex = verify.main(G, "/Users/keshen/PycharmProjects/sk/ISI/simplicial_complex/data/test/", "powerLaw")
    #
    #             s_simplices = simplex.s_simplices
    #             t_simplices = simplex.t_simplices
    #             dis, params = adjTriDistribution(t_simplices)
    #             Dis.append(dis)
    #             Legends.append("m_{}_P_{}_No. ∆_{}_lambda_{}_mu_{}_sigma_{}".format(m, p, len(simplex.triangles), round(1 / params[0] / params[2], 2),
    #                                                                              round(params[1], 2), round(params[2], 2)))
    #     visualizeDistribution(Dis, Legends, "N_{}".format(n))


    # for i in range(100000):
    #     G = nx.powerlaw_cluster_graph(4, 2, 0.1)
    #     cliq_list = list(nx.clique.enumerate_all_cliques(G))
    #     traingle_list = [x for x in cliq_list if len(x) == 3]
    #     cnt += (len(traingle_list))
    # print(cnt)
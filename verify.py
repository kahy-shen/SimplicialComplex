# -*- coding: utf-8 -*-
# @Time    : 10/8/20 9:21 PM
# @Author  : Kay
# @Email   : kahy.shen@gmail.com
# @File    : simplices.py
# @Desc    :
import matplotlib

matplotlib.use('Agg')
import numpy as np
import networkx as nx
import copy
import time
import matplotlib.pyplot as plt
import json
import random
from scipy.optimize import curve_fit
import os


def adjency_matrix(uppermatrix, lowermatrix):
    adjcancy_dict = {}
    for key in uppermatrix:
        adjcancy_dict[key] = []
        upper = set(uppermatrix[key])
        for edge in lowermatrix[key]:
            if edge not in upper:
                adjcancy_dict[key].append(edge)
    return adjcancy_dict


def degreeCentrality(G):
    return nx.degree_centrality(G)


def closenessCentrality(G):
    return nx.closeness_centrality(G)


def subgraphCentrality(G):
    return nx.subgraph_centrality(G)


def construct_Graph(adj_dict, G):
    for key in adj_dict.keys():
        edges = [(key, e) for e in adj_dict[key]]
        # print(edges)
        G.add_edges_from(edges)

    return G


def remove_repetitions(candidates):
    can = list(set([frozenset(list(c)) for c in candidates]))
    return [set(list(c)) for c in can]


def func(x, alpha, c):
    return c * pow(x, -alpha)


class simplices:
    def __init__(self, K, G, folderDir, filename):
        self.tri_lower_adj = {}
        self.quads = []  # quads
        self.tri_upper_adj = {}  # line: [other lines composed of quad]
        self.Triangle2Quad = {}
        self.Line2Triangle = {}
        self.graphs = {i: None for i in range(K)}
        self.folderDir = folderDir
        self.filename = filename
        self.G1 = G
        self._initialize()

    def _initialize(self):
        self._read_data()
        self._update()

    def _read_data(self):
        # print("reading")
        self.graphs[0] = self.G1
        print("nodes: ", len(self.G1.nodes()), "lines: ", len(self.G1.edges()))
        cliq_list = list(nx.clique.enumerate_all_cliques(self.G1))
        traingle_list = [x for x in cliq_list if len(x) == 3]
        quad_list = [x for x in cliq_list if len(x) == 4]
        print("nx triangles: ", len(traingle_list), "quads: ", len(quad_list))
        Coef = nx.average_clustering(self.G1)
        degreeDis = nx.degree(self.G1)
        Dis = {}
        for d in degreeDis:
            if d[1] not in Dis:
                Dis[d[1]] = 1
            else:
                Dis[d[1]] += 1
        x = sorted(list(Dis.keys()))
        S = sum([Dis[x] for x in Dis])
        y = [Dis[k] / S for k in x]

        plt.plot(x, y, 'b-')
        popt, pcov = curve_fit(func, x, y)

        y2 = [func(i, popt[0], popt[1]) for i in x]
        plt.loglog(x, y, 'g+-')
        plt.loglog(x, y2, 'r--',
                   label="alpha={}, c={}, clusCoef ={}".format(round(popt[0], 4), round(popt[1], 4), round(Coef, 4)))
        plt.legend()
        # plt.savefig('/Users/keshen/PycharmProjects/sk/ISI/simplicial_complex/output/RandomNetWork/enronRand/powerLawCluster/powerLaw-m{}-p{}.png'.format(m, p))
        plt.savefig("{}{}/degreedistribution.png".format(self.folderDir, self.filename))
        plt.close()


    def _findTriangle(self):
        lines = {frozenset(l) for l in self.G1.edges()}
        self.line_upper_adj = {}  # line: [other lines composed of triangles]
        self.triangles = []  # triangles
        self.line_lower_adj = {}

        for l1 in lines:
            self.line_lower_adj[l1] = []
            self.line_upper_adj[l1] = []

            L = list(l1)
            l_0_nei = {n for n in self.G1.neighbors(L[0])}
            l_1_nei = {n for n in self.G1.neighbors(L[1])}
            for n in l_0_nei:
                self.line_lower_adj[l1].append(frozenset([L[0], n]))
            for n in l_1_nei:
                self.line_lower_adj[l1].append(frozenset([L[1], n]))
            self.line_lower_adj[l1] = list(set(self.line_lower_adj[l1]))
            third_nodes = l_0_nei.intersection(l_1_nei)
            if third_nodes:
                for n in third_nodes:
                    if frozenset([L[0], n]) in lines and frozenset([L[1], n]) in lines:
                        l1_ = set(list(copy.deepcopy(l1)))
                        l1_.add(n)

                        self.line_upper_adj[l1].append(frozenset([L[0], n]))
                        self.line_upper_adj[l1].append(frozenset([L[1], n]))

                        self.triangles.append(l1_)

                        if l1 in self.Line2Triangle:
                            self.Line2Triangle[l1].append(frozenset(list(l1_)))
                        else:
                            self.Line2Triangle[l1] = [frozenset(list(l1_))]
                        if frozenset([L[0], n]) in self.Line2Triangle:
                            self.Line2Triangle[frozenset([L[0], n])].append(frozenset(list(l1_)))
                        else:
                            self.Line2Triangle[frozenset([L[0], n])] = [frozenset(list(l1_))]
                        if frozenset([L[1], n]) in self.Line2Triangle:
                            self.Line2Triangle[frozenset([L[1], n])].append(frozenset(list(l1_)))
                        else:
                            self.Line2Triangle[frozenset([L[1], n])] = [frozenset(list(l1_))]

    def find_quad(self):
        print("No. triangles: ", len(self.triangles))

        tris = {frozenset(list(t)) for t in self.triangles}

        for t in tris:
            self.tri_upper_adj[t] = []
            self.tri_lower_adj[t] = []

            T = list(t)
            t_0_nei = {n for n in self.G1.neighbors(T[0])}
            t_1_nei = {n for n in self.G1.neighbors(T[1])}
            t_2_nei = {n for n in self.G1.neighbors(T[2])}
            for n in t_0_nei.intersection(t_1_nei):
                self.tri_lower_adj[t].append(frozenset([T[0], T[1], n]))
            for n in t_1_nei.intersection(t_2_nei):
                self.tri_lower_adj[t].append(frozenset([T[1], T[2], n]))
            for n in t_0_nei.intersection(t_2_nei):
                self.tri_lower_adj[t].append(frozenset([T[0], T[2], n]))
            self.tri_lower_adj[t] = list(set(self.tri_lower_adj[t]))

            fourth_nodes = t_0_nei.intersection(t_1_nei).intersection(t_2_nei)
            for n in fourth_nodes:
                Q = set(list(copy.deepcopy(t)))
                Q.add(n)
                other_triangles = [frozenset(list(copy.deepcopy(Q) - {n_})) for n_ in Q if n_ != n]
                if other_triangles[0] in tris and other_triangles[1] in tris and other_triangles[2] in tris:
                    for t_ in other_triangles:
                        self.tri_upper_adj[t].append(t_)
                    self.quads.append(Q)

                    if t in self.Triangle2Quad:
                        self.Triangle2Quad[t].append(Q)
                    else:
                        self.Triangle2Quad[t] = [Q]
                    if frozenset([T[0], T[1], n]) in self.Triangle2Quad:
                        self.Triangle2Quad[frozenset([T[0], T[1], n])].append(Q)
                    else:
                        self.Triangle2Quad[frozenset([T[0], T[1], n])] = [Q]
                    if frozenset([T[0], T[2], n]) in self.Triangle2Quad:
                        self.Triangle2Quad[frozenset([T[0], T[2], n])].append(Q)
                    else:
                        self.Triangle2Quad[frozenset([T[0], T[2], n])] = [Q]
                    if frozenset([T[1], T[2], n]) in self.Triangle2Quad:
                        self.Triangle2Quad[frozenset([T[1], T[2], n])].append(Q)
                    else:
                        self.Triangle2Quad[frozenset([T[1], T[2], n])] = [Q]

    def _update(self):
        self._findTriangle()
        self.triangles = remove_repetitions(self.triangles)
        line_adj = adjency_matrix(self.line_upper_adj, self.line_lower_adj)
        # print(line_adj)
        self.G2 = nx.Graph()
        self.G2 = construct_Graph(line_adj, self.G2)
        self.graphs[1] = self.G2
        self.find_quad()
        self.quads = remove_repetitions(self.quads)
        tri_adj = adjency_matrix(self.tri_upper_adj, self.tri_lower_adj)
        self.triadj = tri_adj
        self.G3 = nx.Graph()
        self.G3 = construct_Graph(tri_adj, self.G3)
        self.graphs[2] = self.G3


        print('No. quads:', len(self.quads))

    def special_simplices_S(self):
        # l2t = {(id2user[k] for k in list(key)): [str((id2userfor x_ in list(x))) for x in self.Line2Triangle[key]] for key in self.Line2Triangle}
        l2t = {key: list(set(self.Line2Triangle[key])) for key in self.Line2Triangle}
        self.s_simplices = l2t
        # l2t = {str(key): str(list(set(self.Line2Triangle[key]))) for key in self.Line2Triangle}
        # # print(topL2T)
        # with open('{}{}/simplices_S_name.json'.format(self.folderDir, self.filename), 'w') as outfile:
        #     l2t = json.dumps(l2t, indent=4)
        #     outfile.write(l2t)
        #     outfile.close()

    def special_simplices_T(self):
        simplices2simplices = {}
        adj_matrix = nx.adjacency_matrix(self.graphs[2])
        nodes = [n for n in self.graphs[2].nodes()]
        for i in range(adj_matrix.shape[0]):
            source = nodes[i]
            target = [nodes[t] for t in np.where(adj_matrix[i].toarray() == 1)[1] if t != i]
            simplices2simplices[source] = list(set(target))

        strong_simplices = {}
        for source in simplices2simplices:

            face = set([])
            for target in simplices2simplices[source]:
                # print(target)
                face.add(target.intersection(source))
            if len(face) == 3:
                strong_simplices[source] = simplices2simplices[source]

        simplices2simplices = {key: list(set(simplices2simplices[key])) for key in simplices2simplices.keys()}
        strongtsimplices = {key: list(set(strong_simplices[key])) for key in strong_simplices.keys()}

        self.t_simplices = simplices2simplices
        self.strong_t_simplices = strongtsimplices
        # simplices2simplices = {str(key): str(simplices2simplices[key])
        #                        for key in simplices2simplices.keys()}
        # strong_simplices = {str(key):str(strong_simplices[key])
        #                     for key in strong_simplices.keys()}
        #
        # with open('{}{}/simplices_T_name.json'.format(self.folderDir, self.filename), 'w') as outfile:
        #     s2s = json.dumps(simplices2simplices, indent=4)
        #     outfile.write(s2s)
        #     outfile.close()
        # with open('{}{}/strong_simplices_T_name.json'.format(self.folderDir, self.filename),
        #           'w') as outfile:
        #     s2s = json.dumps(strong_simplices, indent=4)
        #     outfile.write(s2s)
        #     outfile.close()


def main(G, folderName, i):
    # "data/enron/random/"
    path = "{}{}/".format(folderName, i)
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)

    simplex = simplices(3, G, folderName, i)

    # return simplex
    simplex.special_simplices_S()
    simplex.special_simplices_T()
    return simplex




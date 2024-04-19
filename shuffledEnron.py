# -*- coding: utf-8 -*-
# @Time    : 10/8/20 9:21 PM
# @Author  : Kay
# @Email   : kahy.shen@gmail.com
# @File    : simplices.py
# @Desc    :

import numpy as np
import networkx as nx
import copy
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import random

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
    # labels = {n: str(n) for n in G.nodes}
    # nx.draw(G, labels=labels)
    # plt.show()
    # dc = degreeCentrality(G)
    # DC = {'dc': {str(key): dc[key] for key in dc}}
    # cc = closenessCentrality(G)
    # CC = {'cc': {str(key): cc[key] for key in cc}}
    # sc = subgraphCentrality(G)
    # SC = {'sc': {str(key): sc[key] for key in sc}}
    #
    # with open('data/enron/centrality.json', 'a+') as outfile:
    #     dc = json.dumps(DC, indent=4)
    #     cc = json.dumps(CC, indent=4)
    #     sc = json.dumps(SC, indent=4)
    #     outfile.write(dc + "\n")
    #     outfile.write(cc + "\n")
    #     outfile.write(sc + "\n")
    #     outfile.close()
    # # print("dc: ", degreeCentrality(G))
    # # print("cc: ", closenessCentrality(G))
    # # print("sc: ", subgraphCentrality(G))
    return G


def remove_repetitions(candidates):
    can = list(set([frozenset(list(c)) for c in candidates]))
    return [set(list(c)) for c in can]


class simplices:
    def __init__(self, K, filename, id2user, folderDir):
        self.tri_lower_adj = {}
        self.quads = []  # quads
        self.tri_upper_adj = {}  # line: [other lines composed of quad]
        self.Triangle2Quad = {}
        self.Line2Triangle = {}
        self.graphs = {i: None for i in range(K)}
        self.filename = filename
        self.folderDir = folderDir
        self.id2user = id2user
        self._initialize()

    def _initialize(self):
        self._read_data()
        # self._update()

    def _read_data(self):
        # print("reading")
        with open(self.filename, "r") as reader:
            data = reader.readlines()
            nodes = []
            lines = [d.split("\n")[0].split("\t") for d in data]
            self.lines = [[int(l[0]), int(l[1])] for l in lines]
            reader.close()
            for l in lines:
                nodes.append(l[0])
                nodes.append(l[1])
            self.nodes = list(set(nodes))
        # print("nodes")
        G1 = nx.Graph()
        G1.add_edges_from([(e[0], e[1]) for e in self.lines])
        Num = len(G1.edges())
        Edge = []
        for e in list(G1.edges()):
            if e[0] < e[1]:
                Edge.append((e[0], e[1]))
            else:
                Edge.append((e[1], e[0]))
        assert len(Edge) == Num, print("初始就不对")
        Edge = set(Edge)
        for i in range(10 * Num):
            if i % 100000 == 0:
                print(i)
            edgePair = random.sample(Edge, 2)

            idx0V = edgePair[0][1]
            idx1V = edgePair[1][1]
            E0 = sorted([edgePair[0][0], idx1V])
            E1 = sorted([edgePair[1][0], idx0V])
            if idx1V != edgePair[0][0] and idx0V != edgePair[1][0] and (E0[0], E0[1]) not in Edge or ((E1[0], E1[1])) not in Edge:
                for e in edgePair:
                    Edge.remove((e[0], e[1]))

                Edge.add((edgePair[0][0], edgePair[0][1]))
                Edge.add((edgePair[1][0], edgePair[1][1]))
                assert len(Edge) == Num, print("中途出错")
            else:
                i -= 1
        Edge = list(Edge)
        self.G1 = nx.Graph()
        self.G1.add_edges_from([(e[0], e[1]) for e in Edge])
        print("nodes: ", len(self.G1.nodes()), "lines: ", len(self.G1.edges()))
        self.graphs[0] = self.G1
        # self.degreeDis = self.G1.degree([n for n in list(self.G1.nodes())])
        # Dis = {}
        # for d in self.degreeDis:
        #     if d[1] not in Dis:
        #         Dis[d[1]] = 1
        #     else:
        #         Dis[d[1]] += 1
        # x = sorted(list(Dis.keys()))
        # y = [Dis[k] for k in x]
        # plt.plot(x, y)
        # plt.title('degree_distribution'.format(self.folderDir))
        # plt.xlabel('row')
        # plt.ylabel('column')
        # plt.savefig('data/enron/shuffledEnron/shuffled/{}/degree_distribution.png'.format(self.folderDir))
        # # plt.show()

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

        self.G3 = nx.Graph()
        self.G3 = construct_Graph(tri_adj, self.G3)
        self.graphs[2] = self.G3

        dc = degreeCentrality(self.G3)
        DC = {key: dc[key] for key in dc}
        DC = {
            str(frozenset([id2user[k] for k in list(key)])): DC[key] for key in DC.keys()}

        with open('data/enron/shuffledEnron/shuffled/{}/degreeCentrality.json'.format(self.folderDir), 'w') as outfile:
            dc = json.dumps(DC, indent=4)
            outfile.write(dc + "\n")
            outfile.close()

        print('No. quads:', len(self.quads))

    def special_simplices_S(self, id2user):
        # l2t = {(id2user[k] for k in list(key)): [str((id2userfor x_ in list(x))) for x in self.Line2Triangle[key]] for key in self.Line2Triangle}

        l2t = {str(frozenset([id2user[k] for k in list(key)])): [str(frozenset([id2user[x_] for x_ in list(x)])) for x in self.Line2Triangle[key]] for key in self.Line2Triangle}

        l2t_num = {key: len(l2t[key]) for key in l2t.keys()}
        topK_key = sorted(l2t_num.items(), key=lambda x: x[1], reverse=True)[:20]
        topL2T = {k[0]: l2t[k[0]] for k in topK_key}
        # print(topL2T)
        with open('data/enron/shuffledEnron/shuffled/{}/simplices_S_name.json'.format(self.folderDir), 'w') as outfile:
            l2t = json.dumps(l2t, indent=4)
            outfile.write(l2t)
            outfile.close()
        with open('data/enron/shuffledEnron/shuffled/{}/top_20_simplices_S_name.json'.format(self.folderDir), 'w') as outfile:
            topL2T = json.dumps(topL2T, indent=4)
            outfile.write(topL2T)
            outfile.close()
        # return self.Line2Triangle


    def special_simplices_T(self, id2user):
        simplices2simplices = {}
        adj_matrix = nx.adjacency_matrix(self.graphs[2])
        # print(self.graphs[k].nodes(data=True))
        nodes = [n for n in self.graphs[2].nodes()]
        # print(self.graphs[k].nodes(0))
        for i in range(adj_matrix.shape[0]):
            source = nodes[i]
            target = [nodes[t] for t in np.where(adj_matrix[i].toarray() == 1)[1] if t != i]
            simplices2simplices[source] = target

        strong_simplices = {}
        for source in simplices2simplices:

            face = set([])
            for target in simplices2simplices[source]:
                # print(target)
                face.add(target.intersection(source))
            if len(face) == 3:
                strong_simplices[source] = simplices2simplices[source]


        simplices2simplices = {str(frozenset([id2user[k] for k in list(key)])):
                                   [str(frozenset([id2user[x_] for x_ in list(x)])) for x in simplices2simplices[key]]
                               for key in simplices2simplices.keys()}
        strong_simplices = {str(frozenset([id2user[k] for k in list(key)])):
                                   [str(frozenset([id2user[x_] for x_ in list(x)])) for x in strong_simplices[key]]
                               for key in strong_simplices.keys()}
        # print(len(strong_simplices))
        t2t_num = {key: len(simplices2simplices[key]) for key in simplices2simplices.keys()}
        topK_key = sorted(t2t_num.items(), key=lambda x: x[1], reverse=True)[:20]
        topT2T = {k[0]: simplices2simplices[k[0]] for k in topK_key}

        with open('data/enron/shuffledEnron/shuffled/{}/simplices_T_name.json'.format(self.folderDir), 'w') as outfile:
           s2s = json.dumps(simplices2simplices, indent=4)
           outfile.write(s2s)
           outfile.close()
        with open('data/enron/shuffledEnron/shuffled/{}/strong_simplices_T_name.json'.format(self.folderDir), 'w') as outfile:
           s2s = json.dumps(strong_simplices, indent=4)
           outfile.write(s2s)
           outfile.close()
        with open('data/enron/shuffledEnron/shuffled/{}/top_20_simplices_T_name.json'.format(self.folderDir), 'w') as outfile:
           s2s = json.dumps(topT2T, indent=4)
           outfile.write(s2s)
           outfile.close()

        # return simplices2simplices
    def special_simplices_P(self, id2user):
        shortestPathBetweenT = {}
        import random
        nodes = [n for n in self.graphs[2].nodes()]
        randsourcenodeidx = random.sample(range(len(nodes)), 5)
        randtargetnodeidx = random.sample(range(len(nodes)), 5)

        for s in randsourcenodeidx:
            source = nodes[s]
            for t in randtargetnodeidx:
                target = nodes[t]
                if source != target:
                    if nx.has_path(self.graphs[2], source, target):
                        shortestPathBetweenT[frozenset([source, target])] = nx.shortest_path(self.graphs[2], source, target)

        shortestPathBetweenT = {str(frozenset([frozenset([id2user[k] for k in list(list(key)[0])]), frozenset([id2user[k_] for k_ in list(list(key)[1])])])):
                                   [str(frozenset([id2user[x_] for x_ in list(x)])) for x in shortestPathBetweenT[key]]
                               for key in shortestPathBetweenT.keys()}
        # print(shortestPathBetweenT)
        # shortPath_num = {key: len(shortestPathBetweenT[key]) for key in shortestPathBetweenT.keys()}
        # topK_key = sorted(shortPath_num.items(), key=lambda x: x[1], reverse=True)[:20]
        # topT2T = {k[0]: shortestPathBetweenT[k[0]] for k in topK_key}

        with open('data/enron/shuffledEnron/shuffled/{}/simplices_P_name.json'.format(self.folderDir), 'w') as outfile:
            s2s = json.dumps(shortestPathBetweenT, indent=4)
            outfile.write(s2s)
            outfile.close()


def readUserId(filename):
    with open(filename, 'r') as f:
        data = json.loads(f.read())
        f.close()
    return data


if __name__ == '__main__':
    # folder = os.path.exists("data/enron/shuffledEnron/shuffled/{}".format(1))
    # if not folder:
    #     os.makedirs("data/enron/shuffledEnron/shuffled/{}".format(1))
    # user2id = readUserId("data/enron/user.json")
    # id2user = {user2id[key]: key for key in user2id.keys()}
    #
    # simplex = simplices(3, "data/enron/enron_edges.txt", id2user, 1)
    # # simplex = simplices(3, "data/test.txt")
    # #
    # # simplex.special_simplices_S(id2user)
    # # simplex.special_simplices_T(id2user)
    # # simplex.special_simplices_P(id2user)
    for i in range(10):
        folder = os.path.exists("data/enron/shuffledEnron/shuffled/{}".format(i))
        if not folder:
            os.makedirs("data/enron/shuffledEnron/shuffled/{}".format(i))
        user2id = readUserId("data/enron/user.json")
        id2user = {user2id[key]: key for key in user2id.keys()}

        simplex = simplices(3, "data/enron/enron_edges.txt", id2user, i)
        # simplex = simplices(3, "data/test.txt")

        simplex.special_simplices_S(id2user)
        simplex.special_simplices_T(id2user)
        simplex.special_simplices_P(id2user)

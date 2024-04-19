# -*- coding: utf-8 -*-
# @Time    : 12/19/20 11:33 PM
# @Author  : Kay
# @Email   : kahy.shen@gmail.com
# @File    : para2para.py
# @Desc    :
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import exponnorm

def readCsv(filename):
    csvFile = open(filename, "r")
    reader = csv.reader(csvFile)
    csvdata = [item for item in reader]
    return csvdata
def funcXSquare(x, alpha, beta, C):
    return C * pow(beta*pow(x, -alpha), np.log(x))

def visualPCoef(data):
    plt.figure(figsize=(8, 8.5))
    colors = ['g', 'r', 'b', 'y', 'o']
    marks = ['-', '--', '-.', ':', '.', 'o', 'v', '^', '<', '>', 'p']
    X = [line.split('_') for line in data]
    # print(X)
    N_range = list(set([int(x[0]) for x in X]))
    # N_range.remove(0)
    M_range = list(set([int(x[1]) for x in X]))
    p_range = sorted(list(set([round(float(x[2]), 1) for x in X])))
    X_p = {}
    Y_p = {}
    for n in N_range:
        for m in M_range:
            X_p["{}".format(n)] = []
            Y_p["{}".format(n)] = []
            for p in p_range:
                if "{}_{}_{}".format(n, m, p) in data:
                    X_p["{}".format(n)].append(p)
                    Y_p["{}".format(n)].append(float(data["{}_{}_{}".format(n, m, p)]))
    for key in X_p:
        n = int(key.split('_')[0])
        # m = int(key.split('_')[1])
        plt.plot(X_p[key], Y_p[key], '{}'.format(colors[N_range.index(n)]), label=f'N={key}')
    plt.legend()
    plt.xlabel('P')
    plt.ylabel('cluster coefficient')
    # plt.savefig("relationship between triad formation probability and cluster coefficient.png", dpi = 300)


def visualNMPlambda(data, para, flag):
    plt.figure(figsize=(12, 12))
    colors = ['g', 'r', 'b', 'y', 'o']
    marks = ['-', '--', '-.', ':', '.', 'o-', 'v-', '^-', '<-', '>-', 'p-']
    X = [line.split('_') for line in data]
    N_range = list(set([int(x[0]) for x in X]))
    # N_range.remove(0)
    print(N_range)
    M_range = list(set([int(x[1]) for x in X]))
    p_range = sorted(list(set([round(float(x[2]), 1) for x in X])))

    # plt.subplot(1, 2, 1)  # 只变p
    # X_p = {}
    # Y_p = {}
    # for n in N_range:
    #     for m in M_range:
    #         X_p["{}_{}".format(n, m)] = []
    #         Y_p["{}_{}".format(n, m)] = []
    #         for p in p_range:
    #             if "{}_{}_{}".format(n, m, p) in data:
    #                 X_p["{}_{}".format(n, m)].append(p)
    #                 Y_p["{}_{}".format(n, m)].append(float(data["{}_{}_{}".format(n, m, p)]))
    # for key in X_p:
    #     n = int(key.split('_')[0])
    #     m = int(key.split('_')[1])
    #     plt.plot(X_p[key], Y_p[key], '{}{}'.format(colors[N_range.index(n)], marks[M_range.index(m)]), label=key)
    # # plt.legend()
    # plt.xlabel('P')
    # plt.ylabel('Parameter {} in the estimated function of {} simplices distribution'.format(para, flag))
    ######
    plt.subplot(1, 2, 2)  # 只变m
    X_p = {}
    Y_p = {}
    for n in N_range:
        for p in p_range:
            X_p["{}_{}".format(n, p)] = []
            Y_p["{}_{}".format(n, p)] = []
            for m in M_range:
                if "{}_{}_{}".format(n, m, p) in data:
                    X_p["{}_{}".format(n, p)].append(m)
                    Y_p["{}_{}".format(n, p)].append(float(data["{}_{}_{}".format(n, m, p)]))
    for key in X_p:
        n = int(key.split('_')[0])
        p = float(key.split('_')[1])
        plt.plot(X_p[key], Y_p[key], '{}{}'.format(colors[N_range.index(n)], marks[0]), label=f'N={n}')
    plt.legend()
    plt.xlabel('M')
    plt.ylabel('Parameter {} in the estimated function of {} simplices distribution'.format(para, flag))
    plt.show()
    # plt.savefig("relationship between generation parameter and {}.png".format(para), dpi = 300)


if __name__ == '__main__':
    data = readCsv("para2para.csv")
    para2para = []
    for line in data[1:]:
        _ = line[1].split('_')
        for x in line[2:]:
            _.append(float(x))
        para2para.append(_)
    print(data[0])
    para2para = pd.DataFrame(para2para, columns=['N', 'M', 'P', 'clustercoef', 'lambda', 'mu', 'sigma', 'a', 'b', 'c'])

    chosen_N = '40000'
    Ms = ['4', '6', '8']
    Ps = ['0.2', '0.5', '0.8']
    results = pd.DataFrame()
    for M in Ms:
        for P in Ps:
            filtered_df = para2para[(para2para['N'] == chosen_N) &
                                    (para2para['M'] == M) &
                                    (para2para['P'] == P)]            
            if not filtered_df.empty:
                results = pd.concat([results, filtered_df[['N', 'M', 'P', 'clustercoef', 'lambda', 'mu', 'sigma', 'a', 'b', 'c']]], ignore_index=True)
    print(results)

    plt.figure(figsize=(8, 4))
    df_M = results[(results['M'] == '8')]
    color = ['#02B5D8', '#FE4A52', '#73BD52']
    cnt = 0
    for index, row in df_M.iterrows():
        mu = row['mu']
        sigma = row['sigma']
        lambda_ = row['lambda']  
        K = 1 / lambda_
        x = np.linspace(exponnorm.ppf(0.01, K, loc=mu, scale=sigma),
                exponnorm.ppf(0.99, K, loc=mu, scale=sigma), 100)

        pdf = exponnorm.pdf(x, K, loc=mu, scale=sigma)
        
        plt.plot(x, pdf, '-', color=color[cnt], lw=2, label=f"N={row['N']}, M={row['M']}, P={row['P']}")
        cnt += 1
    plt.title('Estimated function of T*-complexes')
    # plt.xlabel('Value')
    # plt.ylabel('Probability Density')
    plt.legend(loc='best')
    plt.grid(True)
    plt.savefig('T-complex-controlled-P.png', dpi=300)
    plt.show()

    x = np.linspace(0.1, 10, 400)

    plt.figure(figsize=(8, 4))
    cnt = 0
    for index, row in df_M.iterrows():
        alpha = row['a']
        beta = row['b']
        C = row['c']  
        y = funcXSquare(x, C, alpha, beta)
        plt.plot(x, y, '-', color=color[cnt], lw=2, label=f"N={row['N']}, M={row['M']}, P={row['P']}")
        cnt += 1
    plt.title('Estimated function of S*-complexes')
    # plt.xlabel('x')
    # plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.savefig('S-complex-controlled-P.png', dpi=300)

    plt.show()

   


    # # for line in data[1:]:
    #     # N, m, p = line[1].split('_')
    #     # print(N, m, round(float(p), 1))
    # # print(data[1][1:-1])
    # dataDict = {line[1]: line[2:]for line in data[1:]}
    # print(dataDict.keys())
    # NMPpara = {key: dataDict[key][0] for key in dataDict}
    # visualPCoef(NMPpara)
    # NMPpara = {key: dataDict[key][1] for key in dataDict}
    # visualNMPlambda(NMPpara, 'lambda', 'T')
    # NMPpara = {key: dataDict[key][2] for key in dataDict}
    # visualNMPlambda(NMPpara, 'mu', 'T')
    # NMPpara = {key: dataDict[key][3] for key in dataDict}
    # visualNMPlambda(NMPpara, 'sigma', 'T')
    # NMPpara = {key: dataDict[key][4] for key in dataDict}
    # visualNMPlambda(NMPpara, 'a', 'S')
    # NMPpara = {key: dataDict[key][5] for key in dataDict}
    # visualNMPlambda(NMPpara, 'b', 'S')
    # NMPpara = {key: dataDict[key][6] for key in dataDict}
    # visualNMPlambda(NMPpara, 'c', 'S')
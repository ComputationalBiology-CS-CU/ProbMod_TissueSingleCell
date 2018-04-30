import csv
import pandas as pd
from os.path import isfile, join
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score


import pandas as pd
from sklearn.decomposition import PCA
from ggplot import *
import time
from sklearn.manifold import TSNE
import csv
from numpy import loadtxt
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
import math

import numpy as np
from scipy.spatial.distance import euclidean

import numpy as np
from matplotlib import pyplot as plt

from sklearn.datasets import make_checkerboard
from sklearn.datasets import samples_generator as sg
from sklearn.cluster.bicluster import SpectralBiclustering
from sklearn.metrics import consensus_score

from fastdtw import fastdtw, dtw

from util import mPlot

import seaborn as sns
sns.set(color_codes=True)

_RELOAD = False

def load_sc_reading_counts():
    DIR_PATH="../data/GSE81547_RAW"
    sc_df = None
    files = [f for f in os.listdir(DIR_PATH) if isfile(join(DIR_PATH, f))]

    for i in range(len(files)):
        #if(i>100):break
        if files[i].endswith(".csv"):
            cur_file = DIR_PATH + '/' + files[i]
            gsm_id=files[i][0:10]
            df_ = pd.read_csv(cur_file, delimiter="\t",header=None, names=[gsm_id])
            if(sc_df is None):
                sc_df = df_
            else:
                sc_df[gsm_id] = df_[gsm_id]
    return sc_df

def gen_srr():
    file_path = "../meta/annocomb_GSE81547.csv"
    df = pd.read_csv(file_path, header=None)
    M = df.shape[0] - 1
    #print(df.shape)

    ctype = df.loc[1:, 27].tolist()  # c
    srr_list = df.loc[1:, 9].tolist()  # m
    person = df.loc[1:, 15].tolist()  # 8
    gsm_list = df.loc[1:, 8].tolist()  #
    # person_label = person.unique()

    c_list, p_list = [], []

    data = {}
    for i in range(M):
        cur_c = ctype[i].split(":")[1].strip()
        cur_p = person[i].split("_")[0].strip()
        c_list.append(cur_c)
        p_list.append(cur_p)
        if (cur_p not in data):
            data[cur_p] = []
        data[cur_p].append((i, cur_p, srr_list[i], c_list[i]))

    # print first 50 srr for each person
    SRR_NUM = 10
    for cur_p in data:
        cur_info = data[cur_p]
        #print(cur_p)
        out_str = ""
        for i in range(len(cur_info)):
            if (i < SRR_NUM):
                out_str += str(cur_info[i][2]) + "\n"
        with open(cur_p + ".txt", "w") as text_file:
            text_file.write(out_str)

    return c_list, srr_list, p_list, gsm_list

def get_data_label(df_sc, c_list, gsm_list):
    df = df_sc.transpose()
    df["y"] = "none"

    for c_type in set(c_list):
        c_i = [i for i, x in enumerate(c_list) if x == c_type]

        target_i = list(set(c_i))
        idx_list = []
        for i in target_i:
            idx_list.append(gsm_list[i])


        df.loc[df.index.isin(idx_list), "y"] = c_type

    df_x = df[[col for col in df.columns if col not in ["y"]]]
    df_y = df["y"]

    return df_x, df_y

def plot_pca(df_x, df_y, title_str):
    #####################################
    # PCA
    #####################################
    true_label = pd.Series(df_y.tolist())
    row_colors, _ = mPlot.map_label2color(true_label)

    pca = PCA(n_components=5)
    pca_result = pca.fit_transform(df_x.values)

    #kmeans = KMeans(n_clusters=4).fit(pca_result)
    #keamns_labels = kmeans.labels_

    print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))

    f2 = plt.figure()
    ax = Axes3D(f2)
    ax.scatter(pca_result[:, 0], pca_result[:, 1], pca_result[:, 2],
               c=row_colors, s=50)
    #plt.scatter(Y_pca[:, 0], Y_pca[:, 1], color=mergeddf['Color'])
    ax.set_xlabel('pca0')
    ax.set_ylabel('pca1')
    ax.set_zlabel('pca2')
    ax.set_title(title_str)
    f2.show()

    return

def plot_tsne(df_x, df_y, title_str):
    #####################################
    # t-SNE
    #####################################
    true_label = pd.Series(df_y.tolist())
    row_color_, row_color, label2color  = mPlot.map_label2color(true_label)

    _label_list, _color_list = [], []
    for _label in label2color:
        _label_list.append(_label)
        _color_list.append(label2color[_label])

    #_label_list = pd.Series(_label_list)
    #_color_list = pd.Series(_color_list)

    tsne = TSNE(n_components=3, verbose=1, perplexity=100, n_iter=20000, learning_rate=100, early_exaggeration=100)
    # tsne = TSNE(n_components=2, verbose=1)

    tsne_results = tsne.fit_transform(df_x.values)

    df_tsne = df_x.loc[:, :].copy()
    df_tsne['x-tsne'] = tsne_results[:, 0]
    df_tsne['y-tsne'] = tsne_results[:, 1]
    df_tsne['true_label'] = row_color

    '''
    fig = plt.figure(figsize=(500, 500))
    ax = fig.add_subplot(111)
    ax.set_title(title_str, fontsize=14)
    ax.set_xlabel("x", fontsize=12)
    ax.set_ylabel("y", fontsize=12)
    ax.scatter(df_tsne['x-tsne'], df_tsne['y-tsne'], c=row_color, s=50)
    plt.show()
    '''

    f2 = plt.figure()
    ax = Axes3D(f2)
    ax.scatter(tsne_results[:, 0], tsne_results[:, 1], tsne_results[:, 2],
               c=row_color_, s=50)
    # plt.scatter(Y_pca[:, 0], Y_pca[:, 1], color=mergeddf['Color'])
    ax.set_xlabel('0')
    ax.set_ylabel('1')
    ax.set_zlabel('2')
    ax.set_title(title_str)
    f2.show()
    plt.show()
    '''
    plt.figure()
    chart2 = ggplot(df_tsne, aes(x='x-tsne', y='y-tsne', color='true_label')) \
             + scale_color_manual(labels = _label_list, values = _color_list) \
             + geom_point(size=60, alpha=0.6) \
             + ggtitle("tSNE true label" + title_str)
    chart2.show()
    '''

    return


def main():
    c_list, srr_list, p_list, gsm_list = gen_srr()
    if (_RELOAD):
        df_sc = load_sc_reading_counts()
        df_sc.to_pickle("save_df_sc")
    else:
        df_sc = pd.read_pickle('save_df_sc')

    df_x, df_y = get_data_label(df_sc, c_list, gsm_list)

    plot_tsne(df_x, df_y, "tSNE")
    #plot_pca(df_x, df_y, "PCA")


if __name__ == "__main__":
    main()
    print("done.")
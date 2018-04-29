import csv
import pandas as pd
from os.path import isfile, join
import os
import allel
from gtfparse import read_gtf

import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

_DEBUG = True
_DEBUG_ROWS = 2000

_RELOAD = False

_GTEX_VCF = ('../data/gtex/gtex_chr22.vcf', 17) #first row start at 17
_SC_VCF = ('../out/out_10.vcf', 30)


# get G2 x M reading counts matrix
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

def normalize_df(in_df, col_name, out_df):
    # sum in_df and add the col into out_df
    df_sum = in_df.sum(axis=1)
    df_sum = df_sum.to_frame(name=col_name)
    if (out_df is None):
        out_df = df_sum
    else:
        out_df[col_name] = df_sum[col_name]
    return out_df

def get_sc_Y(df_sc, c_list, p_list, gsm_list):
    # for each cell type generate size G2x8
    df_list=[]
    for c_type in set(c_list):
        cur_df = None
        for p_name in set(p_list):
            c_i = [i for i, x in enumerate(c_list) if x == c_type]
            p_i = [i for i, x in enumerate(p_list) if x == p_name]

            target_i = list(set().union(c_i, p_i))
            col_list = []
            for i in target_i:
                col_list.append(gsm_list[i])
            #print(col_list)

            df_ = df_sc[df_sc.columns.intersection(col_list)]
            # todo: may need normalization
            # for each cell type, sum all reading counts of that person
            cur_df = normalize_df(df_, p_name, cur_df)

        df_list.append(cur_df)
    G_list = df_sc.index.tolist()

    print("Y2 shape (CxG1x8): %dx%dx%d"%(len(df_list), df_list[0].shape[0], df_list[0].shape[1]))
    return df_list, G_list

_M_DES_FILE = "../data/gtex/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
def get_gtex_tissue_id():
    df_M_des = pd.read_csv(_M_DES_FILE, header=0, sep="\t")
    M1_id_pancreas = df_M_des.loc[df_M_des['SMTS'] == "Pancreas"]['SAMPID'].tolist()
    return M1_id_pancreas

_GTEX_COUNTS_FILE = "../data/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct"
def get_gtex_Y(tissue_id):
    # return GxM, M individule
    if(_DEBUG):
        df_G1_M = pd.read_csv(_GTEX_COUNTS_FILE, nrows=_DEBUG_ROWS, skiprows=[0,1], header=0,
                              delim_whitespace=True)
    else:
        df_G1_M = pd.read_csv(_GTEX_COUNTS_FILE, skiprows=[0, 1], header=0,
                              delim_whitespace=True)
    df_G1_M = df_G1_M.set_index("Description",drop=False)

    M_cols = [col for col in df_G1_M.columns if col in tissue_id]
    df_G1_M_ = df_G1_M[M_cols]

    # group sample by person
    col_p_dict = {}
    for col in df_G1_M_.columns:
        p_name = col[:10]
        if(p_name not in col_p_dict):
            col_p_dict[p_name] = []
        col_p_dict[p_name].append(col)
    df_G1_M1 = None
    for p_name in col_p_dict:
        cols = col_p_dict[p_name]
        df_ = df_G1_M_[cols]
        # todo: need normalization
        df_G1_M1 = normalize_df(df_, p_name, df_G1_M1)

    print("Y1 shape (G1xM1):"+str(df_G1_M1.shape))
    G1_list = df_G1_M["Description"].tolist()

    return df_G1_M1, G1_list

def get_vcf_counts(vcf_file, first_line):
    #callset = allel.read_vcf(vcf_file)
    #gt = allel.GenotypeArray(callset['calldata/GT'])
    #ac = gt.count_alleles()
    vcf2x = {'./.': 0, '0/0':0, '0/1': 1, '0/2': 1, '1/1':2, '1/2':2, '2/2':2, '2/3':2}
    if(_DEBUG):
        df = pd.read_csv(vcf_file, nrows=_DEBUG_ROWS, header=first_line, sep="\t")
    else:
        df = pd.read_csv(vcf_file, header=first_line, sep="\t")

    df = df.set_index("POS", drop=False)
    N_list = df["POS"].tolist()

    col_M = df.columns[9:]
    df_N_M = df[col_M]
    for m in col_M:
        df_N_M[m] = df_N_M[m].apply(lambda x: vcf2x[x.split(":")[0]])

    return df, df_N_M, N_list



def get_gtex_x(vcf_file, first_line, tissue_id):
    df, df_N_M, N_list = get_vcf_counts(vcf_file, first_line)
    col = [col for col in df_N_M.columns if col in tissue_id]
    df = df[col]
    df_N_M = df_N_M[col]
    return df, df_N_M, N_list

def main():
    print("=====preprocessing=====")
    tissue_id = get_gtex_tissue_id()

    print("=====processing X data=====")
    #df_gtex_vcf, x_gtex_N1xM1, N1_list = get_gtex_x(_GTEX_VCF[0], _GTEX_VCF[1], tissue_id)
    df_gtex_vcf, x_gtex_N1xM1, N1_list = get_vcf_counts(_GTEX_VCF[0], _GTEX_VCF[1])
    df_sc_vcf, x_sc_N2x8, N2_list = get_vcf_counts(_SC_VCF[0], _SC_VCF[1])

    N0_list = list(set(N1_list).intersection(set(N2_list)))
    print("N1:%d; N2:%d ==> intersection: N0: %d" % (len(N1_list), len(N2_list), len(N0_list)))

    x_gtex_N0xM1 = x_gtex_N1xM1.loc[N0_list]
    x_sc_N0x8 = x_sc_N2x8.loc[N0_list]

    ###############################################
    print("=====processing Y data=====")
    c_list, srr_list, p_list, gsm_list = gen_srr()
    if(_RELOAD):
        df_sc = load_sc_reading_counts()
        df_sc.to_pickle("save_df_sc")
    else:
        df_sc = pd.read_pickle('save_df_sc')

    y_gtex_G1xM1, G1_list = get_gtex_Y(tissue_id)
    y_sc_CxG2x8, G2_list = get_sc_Y(df_sc, c_list, p_list, gsm_list)

    # only use intersect gene
    G0_list = list(set(G1_list).intersection(set(G2_list)))
    print("G1:%d; G2:%d ==> intersection: G0: %d" % (len(G1_list), len(G2_list), len(G0_list)))

    y_gtex_G0xM1 = y_gtex_G1xM1.loc[G0_list]
    y_sc_CxG0x8 = []
    for sc_y_c in y_sc_CxG2x8:
        y_sc_CxG0x8.append(sc_y_c.loc[G0_list])

    print("=====fitting Y1=beta*X1 =====")
    y1 = y_gtex_G0xM1.iloc[:,:100].transpose()
    x1 = x_gtex_N0xM1.iloc[:,:100].transpose()
    print("y1 shape:" + str(y1.shape))
    print("x1 shape:" + str(x1.shape))

    # Create linear regression object
    regr = linear_model.LinearRegression()

    # Train the model using the training sets
    regr.fit(x1, y1)

    # The coefficients
    print('Coefficients Beta: \n', regr.coef_)
    print('size of beta: '+str(regr.coef_.shape))

    # Make predictions using the testing set
    y1_pred = regr.predict(x1)

    print("Mean squared error: %.2f"
          % mean_squared_error(y1, y1_pred))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % r2_score(y1, y1_pred))



    return

if __name__ == "__main__":
    main()






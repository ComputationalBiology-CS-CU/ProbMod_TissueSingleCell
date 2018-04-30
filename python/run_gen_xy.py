import csv
import pandas as pd
from os.path import isfile, join
import os
import time
import numpy as np

_RELOAD = True
_DEBUG = False
_DEBUG_ROWS = 50000

_OUT_DIR = "../out/"

_GTEX_VCF = ('../data/gtex/gtex_chr22.vcf', 17) #first row start at 17
_SC_VCF = ('../out/out_20.vcf', 30)
_GENE_LEN_FILE = ('../data/gtex/gencode.v19.genes.v6p_model.patched_contigs.gtf', 6)

_GSE_DIR_PATH="../data/GSE81547_RAW"
_GES_ANNO = "../meta/annocomb_GSE81547.csv"

_M_DES_FILE = "../data/gtex/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
_GTEX_COUNTS_FILE = "../data/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct"



df_gene_len_ = pd.read_csv(_GENE_LEN_FILE[0], header=None, skiprows=range(6), sep="\t")
df_gene_len_["gene_id"] = df_gene_len_.iloc[:, 8].apply(lambda x: x.split(";")[0].split('"')[1])  # get name of gene
df_gene_len_["name"] = df_gene_len_.iloc[:, 8].apply(lambda x: x.split(";")[4].split('"')[1])  # get name of gene
df_gene_len_["length"] = df_gene_len_.iloc[:, 4] - df_gene_len_.iloc[:, 3]  # get length of gene
df_gene_len = df_gene_len_.loc[(df_gene_len_[0] == '22') & (df_gene_len_[2] == 'gene')]  # only use chr22
gene2len={}
gene_id2len={}
for index, row in df_gene_len.iterrows():
    gene2len[row["name"]] = row["length"]
    gene_id2len[row["gene_id"]] = row["length"]
#print("gene length size:"+str(df_gene_len.shape))

# get G2 x M reading counts matrix
def load_sc_reading_counts():
    sc_df = None
    files = [f for f in os.listdir(_GSE_DIR_PATH) if isfile(join(_GSE_DIR_PATH, f))]
    for i in range(len(files)):
        #if(i>100):break
        if files[i].endswith(".csv"):
            cur_file = _GSE_DIR_PATH + '/' + files[i]
            gsm_id=files[i][0:10]
            df_ = pd.read_csv(cur_file, delimiter="\t",header=None, names=[gsm_id])
            if(sc_df is None):
                sc_df = df_
            else:
                sc_df[gsm_id] = df_[gsm_id]
    return sc_df

def gen_srr():
    df = pd.read_csv(_GES_ANNO, header=None)
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

def norm_across(df):
    #df = (df - df.mean(axis=1)) / (df.std(axis=1))
    #df = df / df.sum(axis=0)
    #df = df.sub(df.mean(axis=1), axis=0)
    #df = df
    return df


# https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
# Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
# Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
# Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
def norm_RPKM(df, df_len, df_sum):
    #pre_m = df_sum/1000000
    #df = df / pre_m
    df = df.div(df_len, axis = 'index')
    return df

def normalize_df(in_df, data_cols, index_col, col_name, out_df):
    # sum in_df and add the col into out_df
    #df = in_df[data_cols].median(axis=1)
    df = in_df[data_cols]
    #df = (in_df[data_cols] - in_df[data_cols].mean(axis=1)) / (in_df[data_cols].std(axis=1))
    #df = df.to_frame(name=col_name)
    df[index_col] = in_df[index_col]
    df_sum = in_df[data_cols].sum(axis=0)
    if(index_col=="gene_name"):
        df = df[df[index_col].isin(gene2len)]
        df["len"] = df[index_col].apply(lambda x: gene2len[x])
    elif(index_col=="gene_id"):
        df = df[df[index_col].isin(gene_id2len)]
        df["len"] = df[index_col].apply(lambda x: gene_id2len[x])
    else:
        print("ERROR")
    print("[%s]gene has length:%d->%d" % (col_name, len(in_df.index), len(df.index)))
    if(len(df.index)==0):
        return out_df

    df[data_cols] = norm_RPKM(df[data_cols], df["len"], df_sum)
    df[col_name] = df[data_cols].median(axis=1)
    #df[col_name] = df[data_cols].median(axis=1)

    #df[col_name] = df[col_name]/df["len"]

    #df = (df - df.mean()) / (df.std())
    #df = df.to_frame(name=col_name)
    if (out_df is None):
        out_df = pd.DataFrame(index=df.index)
    out_df[col_name] = df[col_name]
    return out_df

def get_sc_Y(df_sc, c_list, p_list, gsm_list):
    # for each cell type generate size G2x8
    df_list = []
    df_list_ctype = []
    df_sc["gene_name"] = df_sc.index
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
            cols = df_sc.columns.intersection(col_list)
            cur_df = normalize_df(df_sc, cols, "gene_name", p_name, cur_df)
        #cur_df = norm_across(cur_df)
        df_list.append(cur_df)
        df_list_ctype.append(c_type)
    G_list = df_sc.index.tolist()

    print("Y2 shape (CxG1x8): %dx%dx%d"%(len(df_list), df_list[0].shape[0], df_list[0].shape[1]))
    return df_list, G_list, df_list_ctype

def get_gtex_tissue_id():
    df_M_des = pd.read_csv(_M_DES_FILE, header=0, sep="\t")
    M1_id_pancreas = df_M_des.loc[df_M_des['SMTS'] == "Pancreas"]['SAMPID'].tolist()
    return M1_id_pancreas

def get_gtex_Y(tissue_id):
    # return GxM, M individule
    if(_DEBUG):
        df_G1_M = pd.read_csv(_GTEX_COUNTS_FILE, nrows=_DEBUG_ROWS, skiprows=[0,1], header=0,
                                  delim_whitespace=True)
    else:
        save_file_name = "save_df_G1_M"
        df_G1_M = pd.read_csv(_GTEX_COUNTS_FILE, skiprows=[0, 1], header=0,
                              delim_whitespace=True)
        #if (_RELOAD):
        #    df_G1_M = pd.read_csv(_GTEX_COUNTS_FILE, skiprows=[0, 1], header=0,
        #                          delim_whitespace=True)
        #    df_G1_M.to_pickle(save_file_name)
        #else:
        #    df_G1_M = pd.read_pickle(save_file_name)

    df_G1_M = df_G1_M.set_index("Description",drop=False)

    #gene_id = df_G1_M["Name"]
    print("S1:"+str(len(df_G1_M.columns)))
    M_cols = [col for col in df_G1_M.columns if col in tissue_id]
    df_G1_M_ = df_G1_M[M_cols]

    print("duplication index:")
    print(df_G1_M_[df_G1_M_.index.duplicated(keep=False)].index)

    # group sample by person
    col_p_dict = {} # map: person id -> [sample from that person]
    for col in df_G1_M_.columns:
        p_name = col.split("-")[1]
        #p_name = col[:10]
        if(p_name not in col_p_dict):
            col_p_dict[p_name] = []
        col_p_dict[p_name].append(col)
    df_G1_M1 = None
    df_G1_M_["gene_id"] = df_G1_M["Name"]
    for p_name in col_p_dict:
        cols = col_p_dict[p_name]
        # todo: need normalization
        df_G1_M1 = normalize_df(df_G1_M_, cols, "gene_id", p_name, df_G1_M1)

    #df_G1_M1 = norm_across(df_G1_M1)
    print("Y1 shape (G1xM1):"+str(df_G1_M1.shape))
    #G1_list = df_G1_M["Description"].tolist()
    G1_list = df_G1_M1.index.tolist()

    return df_G1_M1, G1_list

def get_vcf_counts(vcf_file, first_line):
    #callset = allel.read_vcf(vcf_file)
    #gt = allel.GenotypeArray(callset['calldata/GT'])
    #ac = gt.count_alleles()
    vcf2x = {'./.': np.nan, '0/0':0, '0/1': 1, '0/2': 1, '1/1':2, '1/2':2, '2/2':2, '0/3':2, '1/3':2, '2/3':2}
    if(_DEBUG):
        df = pd.read_csv(vcf_file, nrows=_DEBUG_ROWS, header=first_line, sep="\t")
    else:
        basename = os.path.basename(vcf_file)
        basename = os.path.splitext(basename)[0]
        save_file_name = "save_"+basename
        if (_RELOAD):
            df = pd.read_csv(vcf_file, header=first_line, sep="\t")
            df.to_pickle(save_file_name)
        else:
            df = pd.read_pickle(save_file_name)

    df = df.set_index("POS", drop=False)
    #N_list = df["POS"].tolist()

    col_M = df.columns[9:]
    df_N_M = df[col_M]
    for m in col_M:
        df_N_M[m] = df_N_M[m].apply(lambda x: vcf2x[x.split(":")[0]])
    print("N total:"+str(df_N_M.shape))
    df_N_M = df_N_M.dropna(axis=0, how='any')
    print("N no ./.:" + str(df_N_M.shape))

    N_list = df.loc[df_N_M.index]["POS"].tolist()

    return df, df_N_M, N_list

def get_gtex_x(vcf_file, first_line, list_pname):
    df, df_N_M, N_list = get_vcf_counts(vcf_file, first_line)
    df_N_M2 = pd.DataFrame(index=df_N_M.index)
    print("gtex VCF NXM size:"+str(df_N_M2.shape))
    for i in range(len(list_pname)):
        cur_name = list_pname[i]
        found = False
        for col in df_N_M.columns:
            p_name = col.split("-")[1]
            if(p_name == cur_name):
                df_N_M2[cur_name] = df_N_M[col]
                if(found):
                    print("[WARNNING] Find more than 1 same person id in cvf file")
                found = True
        #if(not found):
        #    print("[WARNNING] Can't find [%s] in cvf file"%(cur_name))
    return df, df_N_M2, N_list

def main():
    print("=====preprocessing=====")
    tissue_id = get_gtex_tissue_id()

    ###############################################
    print("=====processing Y data=====")
    c_list, srr_list, p_list, gsm_list = gen_srr()
    if(_RELOAD):
        df_sc = load_sc_reading_counts()
        df_sc.to_pickle("save_df_sc")
    else:
        df_sc = pd.read_pickle('save_df_sc')

    y_gtex_G1xM1, G1_list = get_gtex_Y(tissue_id)
    y_sc_CxG2x8, G2_list, y_ctype_C = get_sc_Y(df_sc, c_list, p_list, gsm_list)

    # only use intersect gene
    G0_list = list(set(G1_list).intersection(set(G2_list)))
    print("G1:%d; G2:%d ==> intersection: G0: %d" % (len(G1_list), len(G2_list), len(G0_list)))

    # remove duplication
    y_gtex_G0xM1 = y_gtex_G1xM1.loc[G0_list]
    if(y_gtex_G0xM1.shape[0] != len(G0_list)):
        y_gtex_G0xM1 = y_gtex_G0xM1[~y_gtex_G0xM1.index.duplicated(keep='first')]
    y_gtex_G0xM1 = norm_across(y_gtex_G0xM1)
    y_gtex_G0xM1 = y_gtex_G0xM1.sort_index()

    y_sc_CxG0x8 = []
    for sc_y_c in y_sc_CxG2x8:
        cur_df = sc_y_c.loc[G0_list]
        if (cur_df.shape[0] != len(G0_list)):
            cur_df = cur_df[~cur_df.index.duplicated(keep='first')]
        cur_df = norm_across(cur_df)
        cur_df = cur_df.sort_index()
        y_sc_CxG0x8.append(cur_df)

    print("=====processing X data=====")
    list_pname = y_gtex_G1xM1.columns
    df_gtex_vcf, x_gtex_N1xM2, N1_list = get_gtex_x(_GTEX_VCF[0], _GTEX_VCF[1], list_pname)
    df_sc_vcf, x_sc_N2x8, N2_list = get_vcf_counts(_SC_VCF[0], _SC_VCF[1])

    N0_list = list(set(N1_list).intersection(set(N2_list)))
    print("N1:%d; N2:%d ==> intersection: N0: %d" % (len(N1_list), len(N2_list), len(N0_list)))

    x_gtex_N0xM2 = x_gtex_N1xM2.loc[N0_list]
    x_gtex_N0xM2 = x_gtex_N0xM2.sort_index()

    x_sc_N0x8 = x_sc_N2x8.loc[N0_list]
    x_sc_N0x8 = x_sc_N0x8.sort_index()

    # remove duplication in N0
    if (x_gtex_N0xM2.shape[0] != len(N0_list)):
        x_gtex_N0xM2 = x_gtex_N0xM2[~x_gtex_N0xM2.index.duplicated(keep='first')]
    if (x_sc_N0x8.shape[0] != len(N0_list)):
        x_sc_N0x8 = x_sc_N0x8[~x_sc_N0x8.index.duplicated(keep='first')]

    print("=====matching individual label M=====")
    col_M0 = [col for col in y_gtex_G0xM1.columns if col in x_gtex_N0xM2.columns]
    y_gtex_G0xM0 = y_gtex_G0xM1[col_M0]
    x_gtex_N0xM0 = x_gtex_N0xM2[col_M0]
    print("M1:%d; M2:%d ==> intersection: M0: %d" % (len(y_gtex_G0xM1.columns),
                                                     len(x_gtex_N0xM2.columns), len(col_M0)))

    print("====saving x,y data=====")
    np.savez(_OUT_DIR + "/save_npz", y_ctype_C=y_ctype_C)

    y_gtex_G0xM0.to_pickle(_OUT_DIR+"/save_y_gtex_G0xM0")
    x_gtex_N0xM0.to_pickle(_OUT_DIR+"/save_x_gtex_N0xM0")
    x_sc_N0x8.to_pickle(_OUT_DIR+"/save_x_sc_N0x8")
    for i in range(len(y_sc_CxG0x8)):
        y_sc_CxG0x8[i].to_pickle(_OUT_DIR+"/save_y_sc_CxG0x8_"+str(i))
    return

if __name__ == "__main__":
    t0 = time.time()
    main()
    print("done. time:"+str(time.time()-t0)+" s")






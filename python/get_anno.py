import csv
import pandas as pd

file_path="annocomb_GSE81547.csv"

df = pd.read_csv(file_path, header=None)

M = df.shape[0]-1

print(df.shape)

ctype = df.loc[1:,27].tolist() # c
srr = df.loc[1:, 9].tolist() # m
person = df.loc[1:,15].tolist() # 8
#person_label = person.unique()

ctype_label = []
person_label = []

c_list, p_list = [], []
p_id_arr=[]
p_idx_arr = []
data={}
for i in range(M):
    cur_c = ctype[i].split(":")[1].strip()
    cur_p = person[i].split("_")[0].strip()
    c_list.append(cur_c)
    p_list.append(cur_p)
    if(cur_p not in data):
        data[cur_p]=[]
    data[cur_p].append((i,cur_p,srr[i],c_list[i]))

# print first 50 srr for each person
SRR_NUM = 10
for cur_p in data:
    cur_info = data[cur_p]
    print(cur_p)
    out_str=""
    for i in range(len(cur_info)):
        if(i<SRR_NUM):
            out_str+=str(cur_info[i][2])+"\n"
    with open(cur_p+".txt", "w") as text_file:
        text_file.write(out_str)







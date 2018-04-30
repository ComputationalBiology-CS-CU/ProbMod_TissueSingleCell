import pandas as pd
import time

import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score

import seaborn as sns
sns.set(color_codes=True)

_OUT_DIR = "../out/"

def loading_data():
    npzfile = np.load(_OUT_DIR + "/save_npz.npz")
    y_ctype_C = npzfile['y_ctype_C']

    y_gtex_G0xM0 = pd.read_pickle(_OUT_DIR +'save_y_gtex_G0xM0')
    x_gtex_N0xM0 = pd.read_pickle(_OUT_DIR +'save_x_gtex_N0xM0')
    x_sc_N0x8 = pd.read_pickle(_OUT_DIR +'save_x_sc_N0x8')

    y_sc_CxG0x8=[]
    for i in range(7):
        y_sc_CxG0x8.append(pd.read_pickle(_OUT_DIR+"/save_y_sc_CxG0x8_"+str(i)))

    x_sc_N0x8 = x_sc_N0x8.dropna(axis=0, how='any')
    x_gtex_N0xM0 = x_gtex_N0xM0.dropna(axis=0, how='any')

    return y_ctype_C, y_gtex_G0xM0, x_gtex_N0xM0, x_sc_N0x8, y_sc_CxG0x8

def main():
    print("loading data....")
    y_ctype_C, y_gtex_G0xM0, x_gtex_N0xM0, x_sc_N0x8, y_sc_CxG0x8 = loading_data()

    print("=====preprocessing data=====")
    y_sc_G0x8_avg=None
    for y_sc_G0x8 in y_sc_CxG0x8:
        if y_sc_G0x8_avg is None:
            y_sc_G0x8_avg = y_sc_G0x8
        else:
            y_sc_G0x8_avg += y_sc_G0x8
    y_sc_G0x8_avg = y_sc_G0x8_avg/len(y_sc_CxG0x8)
    for i in range(len(y_sc_CxG0x8)):
        y_sc_CxG0x8[i] -= y_sc_G0x8_avg

    print("=====fitting Y_gtex=beta*X_gtex =====")
    y1 = y_gtex_G0xM0.transpose()
    x1 = x_gtex_N0xM0.transpose()
    print("y1 shape:" + str(y1.shape))
    print("x1 shape:" + str(x1.shape))

    regr = linear_model.LinearRegression() # Create linear regression object
    regr.fit(x1, y1) # Train the model using the training sets
    beta = regr.coef_ # The coefficients
    #print('Coefficients Beta: \n', beta)
    print('size of beta (G0xN0): '+str(beta.shape))

    # Make predictions using the testing set
    y1_pred = regr.predict(x1)

    print("Mean squared error: %.2f"
          % mean_squared_error(y1, y1_pred))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % r2_score(y1, y1_pred))

    print("=====Applying beta on sc genotype: Y_hat=beta*X_sc =====")

    print('size of beta (G0xN0): ' + str(beta.shape))
    print('size of X_sc (N0x8): ' + str(x_sc_N0x8.shape))
    Y_hat = beta.dot(x_sc_N0x8)
    print('size of Y_hat (G0x8): ' + str(Y_hat.shape))

    # normalize Y-hat
    #Y_hat = (Y_hat - Y_hat.mean()) / (Y_hat.std())

    print("=====Covariance between cell-type Y_sc and Y_hat =====")
    Y_hat -= y_sc_G0x8_avg
    Y_hat_ = Y_hat.flatten()
    y_hat_8G0x1 = pd.DataFrame(Y_hat_)
    y_Cyx8G0 = pd.DataFrame(columns=y_ctype_C, index=y_hat_8G0x1.index.values)
    for i in range(len(y_ctype_C)):
        c_type = y_ctype_C[i]
        y_G0x8 = y_sc_CxG0x8[i]
        #print("cell type: %s"%(c_type))
        #print("Y_sc shape:"+str(y_G0x8.shape))
        cur_col = y_G0x8.values.flatten()
        y_Cyx8G0[c_type] = cur_col
    y_Cyx8G0['y_hat'] = Y_hat_

    y_CyxCy_cov = y_Cyx8G0.cov()
    y_CyxCy_corr = y_Cyx8G0.corr()
    print(y_CyxCy_cov)
    print(y_CyxCy_corr)
    #sns.clustermap(y_CyxCy_cov, method='complete', metric="euclidean")
    sns.clustermap(y_CyxCy_corr, method='complete', metric="euclidean")

    #plt.pcolor(y_Cx8G0_cov)
    #plt.yticks(np.arange(0.5, len(y_Cx8G0_cov.index), 1), y_Cx8G0_cov.index)
    #plt.xticks(np.arange(0.5, len(y_Cx8G0_cov.columns), 1), y_Cx8G0_cov.columns)

    y_CyxCy_cov[['y_hat']].plot(kind='bar', title="Cov", figsize=(15, 10), legend=True, fontsize=12)
    y_CyxCy_corr[['y_hat']].plot(kind='bar', title="Correlation", figsize=(15, 10), legend=True, fontsize=12)
    plt.show()

    return

if __name__ == "__main__":
    t0 = time.time()
    main()
    print("done. time:"+str(time.time()-t0)+" s")






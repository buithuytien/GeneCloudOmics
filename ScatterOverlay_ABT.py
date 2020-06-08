import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from numpy import errstate, isneginf, array
import time
import argparse
import mplcursors


# HELPER FUNCTION
# import data and store in dataframe
def getDF(file_name, columns=range(4), cutoff_val=0):
    # columns should indicate: rep1_cond1, rep2_cond1, rep1_cond2, rep2_cond2
    # all replicates should be together.
    data = pd.read_csv(file_name)
    df = pd.DataFrame(data)
    df = df.set_index('Unnamed: 0')
    df = df.iloc[:, columns]
    df.index.name = None

    # remove all genes below cutoff value 
    if cutoff_val > 0:
        df = df[df >= cutoff_val].dropna()

    # EC2_0min and EC2_10min
    df_temp = pd.DataFrame(df.iloc[:, 0])
    df_temp = df_temp.join(df.iloc[:, 2])
    df_temp['Names'] = df_temp.index
    #     print(df_temp.head())

    # EC3_0min and EC3_10min
    df_temp4 = pd.DataFrame(df.iloc[:, 1])
    df_temp4 = df_temp4.join(df.iloc[:, 3])
    df_temp4['Names'] = df_temp4.index
    #     print(df_temp4.head())

    # EC2_0min and EC3_0min
    df_temp2 = pd.DataFrame(df.iloc[:, 0])
    df_temp2 = df_temp2.join(df.iloc[:, 1])
    df_temp2['Names'] = df_temp2.index
    #     print(df_temp2.head())

    # EC2_10min and EC3_10min
    df_temp3 = pd.DataFrame(df.iloc[:, 2])
    df_temp3 = df_temp3.join(df.iloc[:, 3])
    df_temp3['Names'] = df_temp3.index
    #     print(df_temp3.head())

    return df, df_temp, df_temp2, df_temp3, df_temp4


# function for data transformation
def logTrans(df_temp, log_base=10):
    with errstate(divide='ignore'):
        result = np.log10(df_temp.iloc[:, 0]) / np.log10(log_base)
        result2 = np.log10(df_temp.iloc[:, 1]) / np.log10(log_base)
    result[isneginf(result)] = 0
    result2[isneginf(result2)] = 0
    df_temp.iloc[:, 0] = result
    df_temp.iloc[:, 1] = result2

    # remove col with neg values
    #     df_temp = df_temp.loc[(df_temp.iloc[:,0]>=0)]
    #     df_temp = df_temp.loc[(df_temp.iloc[:,1]>=0)]
    return (df_temp)


def overlay(df_temp, df_2_3, theta):
    # overlay df_temp onto df_2_3. both are two-column dataframe
    # returns genes from df_temp that does not overlap with any genes from df_2_3
    start_loop = time.time()
    x_val_inc = []
    x_val_exc1 = []
    y_val_inc = []
    y_val_exc1 = []
    exc1_names = []

    for i in range(0, len(df_temp)):
        radius = np.sqrt((df_temp.iloc[i, 0] - df_2_3.iloc[:, 0]) ** 2 + (df_temp.iloc[i, 1] - df_2_3.iloc[:, 1]) ** 2)
        if (radius < (2 * theta)).any():
            x_val_inc.append(df_temp.iloc[i, 0])
            y_val_inc.append(df_temp.iloc[i, 1])
        else:
            x_val_exc1.append(df_temp.iloc[i, 0])
            y_val_exc1.append(df_temp.iloc[i, 1])
            exc1_names.append(df_temp.iloc[i, 2])  # get names

    df_inc = pd.DataFrame(x_val_inc)
    df_inc.columns = ['X']
    df_inc['Y'] = y_val_inc

    df_exc = pd.DataFrame(x_val_exc1)
    df_exc.columns = ['X']
    df_exc['Y'] = y_val_exc1

    df_similar = pd.DataFrame(exc1_names)  # non-overlapping genes
    end_loop = time.time()
    print("1st loop:", end_loop - start_loop)
    return (exc1_names)


def overlayFig(ax4, df_temp, exc1_names, theta):
    # plot overlay figures with orange = non-DE, green = DE genes
    df_exc = df_temp.loc[exc1_names, :]
    df_inc = df_temp.loc[~df_temp["Names"].isin(exc1_names), :]

    # plotting
    ax4.set_xlabel('Log(X)')
    ax4.set_ylabel('Log(Y)')
    pts_inc = plt.scatter(x=df_inc.iloc[:, 0], y=df_inc.iloc[:, 1], c='orange', marker='o', s=2, zorder=1)
    pts_exc = plt.scatter(x=df_exc.iloc[:, 0], y=df_exc.iloc[:, 1], c='seagreen', marker='o', s=2, zorder=2)

    cor_de = (log_base ** df_exc.iloc[:, :2]).corr().iloc[0, 1]  #
    cor_nonde = (log_base ** df_inc.iloc[:, :2]).corr().iloc[0, 1]  #
    n_de = df_inc.shape[0]
    n_nonde = df_exc.shape[0]
    ax4.legend((pts_inc, pts_exc),
               ("Non-DE (" + str(n_nonde) + ")",  # Cor="+str(round(cor_nonde,3)) + ". N="+str(n_nonde),
                "DE (" + str(n_de) + ")"),  # ". Cor="+str(round(cor_de,3)) + ". N="+str(n_de) ),
               fontsize="x-small", loc="lower right", frameon=False
               )
    plt.gca().set_aspect('equal')

    # set axes limits so scale unchanged during dynamic plotting (also set both to be same for comparison)
    #     x_min, x_max, y_min, y_max = ax4.axis()
    x4_min, x4_max, y4_min, y4_max = ax4.axis()
    ax4.set_xlim(-0.1, x4_max)
    ax4.set_ylim(-0.1, y4_max)
    ax4.xaxis.set_ticks(np.arange(0, math.ceil(max(df_temp.iloc[:, 0])), 1))
    ax4.yaxis.set_ticks(np.arange(0, math.ceil(max(df_temp.iloc[:, 1])), 1))

    ### Marker point converter ###
    # set some variables
    r = theta
    N = 1

    # Calculate radius in pixels (wrt each figure):
    rr_pix4 = (ax4.transData.transform(np.vstack([r, r]).T) -
               ax4.transData.transform(np.vstack([np.zeros(N), np.zeros(N)]).T))
    rpix4, _ = rr_pix4.T

    # Calculate and update size in points (wrt each figure):    
    size_pt4 = (2 * rpix4 / fig.dpi * 72) ** 2
    pts_inc.set_sizes(size_pt4)
    pts_exc.set_sizes(size_pt4)

    return (ax4)


def scatterFig(ax, df_temp, theta, color='darkred'):
    # plot actual figures
    ax.set_xlabel(df_temp.columns[0])
    ax.set_ylabel(df_temp.columns[1])
    pts = plt.scatter(x=df_temp.iloc[:, 0], y=df_temp.iloc[:, 1], c=color, marker='o', s=2)
    plt.gca().set_aspect('equal')

    # set axes limits so scale unchanged during dynamic plotting (also set all to be same for comparison)
    x_min, x_max, y_min, y_max = ax.axis()
    ax.set_xlim(-0.1, x_max)
    ax.set_ylim(-0.1, y_max)
    ax.xaxis.set_ticks(np.arange(0, 6, 1))
    ax.yaxis.set_ticks(np.arange(0, 6, 1))

    ### Marker point converter ###
    # set some variables
    r = theta
    N = 1

    # Calculate radius in pixels (wrt each figure):
    rr_pix = (ax.transData.transform(np.vstack([r, r]).T) -
              ax.transData.transform(np.vstack([np.zeros(N), np.zeros(N)]).T))
    rpix, _ = rr_pix.T

    # Calculate and update size in points (wrt each figure):
    size_pt = (2 * rpix / fig.dpi * 72) ** 2
    pts.set_sizes(size_pt)

    return (ax)


# EXECUTION CODE
# command-line interface
parser = argparse.ArgumentParser(
    description='ScatterOverlay. Generates overlapping scatter from 2 scatter plots.'
)
parser.add_argument('-size', type=float, default=0.0125, help='Marker size')
parser.add_argument('-columns', type=str, default='1,2,3,4', help='Columns to be compared')
parser.add_argument('-filename', type=str, default='data.csv', help='CSV file to be read')
args = parser.parse_args()

theta = args.size
colorlist = ['seagreen', 'rebeccapurple', 'darkred', 'orange']
'''
columns = [int(args.columns.split(',')[0]) - 1, int(args.columns.split(',')[1]) - 1,
           int(args.columns.split(',')[2]) - 1, int(args.columns.split(',')[3]) - 1]
'''
columns = [x1_val - 1, y1_val - 1, x2_val - 1, y2_val - 1]
print(x1_val)
print("x1_val")


print('Generating scatter plots...')
start_loop = time.time()

# prepare dataframe
log_base = 10
print(df)
# df.columns = column_names
# print(df)
# df, df_temp, df_temp2, df_temp3, df_temp4 = getDF(file_name, columns, cutoff_val=0)
df_temp = pd.DataFrame(df.iloc[:, columns[0]])
df_temp = df_temp.join(df.iloc[:, columns[2]])
df_temp['Names'] = df_temp.index
#     print(df_temp.head())

# EC3_0min and EC3_10min
df_temp4 = pd.DataFrame(df.iloc[:, columns[1]])
df_temp4 = df_temp4.join(df.iloc[:, columns[3]])
df_temp4['Names'] = df_temp4.index
#     print(df_temp4.head())

# EC2_0min and EC3_0min
df_temp2 = pd.DataFrame(df.iloc[:, columns[0]])
df_temp2 = df_temp2.join(df.iloc[:, columns[1]])
df_temp2['Names'] = df_temp2.index
#     print(df_temp2.head())

# EC2_10min and EC3_10min
df_temp3 = pd.DataFrame(df.iloc[:, columns[2]])
df_temp3 = df_temp3.join(df.iloc[:, columns[3]])
df_temp3['Names'] = df_temp3.index

# log transformation  
 
df_temp = logTrans(df_temp, log_base)  # rep1_cond1, rep1_cond2. Replicate 1, 2 condition scatter
df_temp4 = logTrans(df_temp4, log_base)  # rep2_cond1, rep2_cond2. Replicate 2, 2 condition scatter
df_temp2 = logTrans(df_temp2, log_base)
df_temp3 = logTrans(df_temp3, log_base)

df_2_3 = pd.DataFrame(np.concatenate((df_temp2, df_temp3), axis=0))
df_2_3.columns = ['X', 'Y', 'Names']
df_2_3 = df_2_3.iloc[:, 0:2].astype('float64')
#     print("df_2_3: \n")
#     print(df_2_3.head())

# get DE genes by overlaying
exc1_names = overlay(df_temp, df_2_3, theta)
# other cmbn
exc2_names = overlay(df_temp4, df_2_3, theta)
exc_common = list(set(exc1_names) & set(exc2_names))
# print("exc1_names: " + str(len(exc1_names)))
# print("exc2_names: " + str(len(exc2_names)))
print("DE genes : " + str(len(exc_common)))

df_diff_genes = pd.DataFrame(exc1_names)  # , columns=["combn1"]
df_diff_genes = pd.concat([df_diff_genes, pd.DataFrame(exc2_names), pd.DataFrame(exc_common)], ignore_index=True,
                          axis=1)
df_diff_genes.columns = ["combn1", "combn2", "common"]
#     print(df_diff_genes.tail())


# save to csv file
# new_file_csv = file_name[:-4] + "_log_" + str(log_base) + "_size_" + str(round(theta, 4)) + ".csv"
# df_diff_genes.to_csv(new_file_csv, index=False)
# print("DE genes saved as ", new_file_csv)

# display image
fig = plt.figure(figsize=(6,6) )
fig.subplots_adjust(hspace=0.3, wspace=0.3)
ax = fig.add_subplot(2,2,1); ax = scatterFig(ax, df_temp,  theta, color=colorlist[0] )
ax = fig.add_subplot(2,2,2); ax = scatterFig(ax, df_temp2, theta, color=colorlist[1] )
ax = fig.add_subplot(2,2,3); ax = scatterFig(ax, df_temp3, theta, color=colorlist[2] )
ax = fig.add_subplot(2,2,4); ax = overlayFig(ax, df_temp, exc1_names=exc_common, theta=theta)

mplcursors.cursor().connect(
    "add", lambda sel: sel.annotation.set_text(df_temp["Names"][sel.target.index]))

mplcursors.cursor().connect(
    "add", lambda sel: sel.annotation.set_text(df_temp2["Names"][sel.target.index]))

mplcursors.cursor().connect(
    "add", lambda sel: sel.annotation.set_text(df_temp3["Names"][sel.target.index]))
            

# save image file
# new_file_img = file_name[:-4] + "_log_" + str(log_base) + "_size_" + str(round(theta, 4)) + "_scatter.png"
# fig.savefig(new_file_img , dpi=1200, bbox_inches='tight')

# show image
plt.show()
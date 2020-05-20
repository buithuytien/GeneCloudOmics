import pandas as pd
import numpy as np

# converts R to Python dataframe
def r2py_DataFrame_test():
        ijv = pd.read_csv("testdata_RAW.csv") 
        df_temp = pd.DataFrame(ijv)
        print (df_temp)

        indexnames = list(df_temp.columns.values) 
        indexnames.pop(0)
        print(indexnames)
        
        df_temp = df_temp.set_index('Unnamed: 0')
        df_temp.index.name = None
        df_temp.columns = np.arange(len(df_temp.columns))
        print (df_temp)
        return df_temp, indexnames


# load dataframe straight from ABioTrans
# don't need set 1st column as index as already done in ABT
def r2py_DataFrame(R_df):
        ijv = R_df
        df_temp = pd.DataFrame(ijv)
        print (df_temp)

        indexnames = list(df_temp.columns.values) 
        print(indexnames)
        
        df_temp.columns = np.arange(len(df_temp.columns))
        print (df_temp)
        return df_temp, indexnames


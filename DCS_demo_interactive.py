#!/usr/bin/env python

import numpy as np
import pandas as pd
import scripts.DigitalCellSorter_interactive as DigitalCellSorter
#import scripts.ReadPrepareDataHCApreviewDataset as HCA
import shutil
import os


if __name__ == '__main__':

    AvailableCPUsCount = 11
    data_Folder = 'data'
    geneListToUse = 'geneLists/testing2.xlsx'
    N_samples_for_distribution = 100
    cellTypeNameForSubclustering = None
    clusterIndex = None   

    print('\n======================\nDone loading raw data!\n======================')

    DigitalCellSorter.DigitalCellSorter().Process(df_expr, 
                                                    dataName, 
                                                    saveDir = 'tsne_output/' + dataName + '/', 
                                                    geneListFileName = geneListToUse,
                                                    N_samples_for_distribution = N_samples_for_distribution,
                                                    AvailableCPUsCount = AvailableCPUsCount,
                                                    clusterIndex=clusterIndex,
                                                    clusterName=cellTypeNameForSubclustering,
                                                    n_clusters=n_clusters,
                                                    n_components_pca=n_components_pca,
                                                    perplexity=perplexity,
                                                    temp_indexnames = temp_indexnames)

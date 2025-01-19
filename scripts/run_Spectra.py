from pathlib import Path
import json
from collections import OrderedDict
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os
import anndata

#spectra imports
import dill
import Spectra as spc
from Spectra import Spectra_util as spc_tl
from Spectra import K_est as kst
from Spectra import default_gene_sets

def run_Spectra(features_file, metadata_file, gex_datafile, genesets_json, seuratToAnnDataDir, print_versions = True):

    if print_versions:
        print('scanpy version: ' + sc.__version__)
        print('scipy version: ' + scipy.__version__)
        print('pandas version: ' + pd.__version__)
        print('anndata version: ' + anndata.__version__)

    #set working directory (inherited from the R session)
    os.chdir(working_directory)
    
    genes_df = pd.read_csv(features_file, header=None)
    genelist = list(genes_df[0].values)
    adataObj = sc.read_10x_h5(gex_datafile)
    adataObj = adataObj[:, genelist]
    metadata = pd.read_csv(metadata_file, index_col=0)
    
    #basic preprocessing
    sc.pp.calculate_qc_metrics(adataObj, percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_genes(adataObj, min_cells=20)

    cells = adataObj.obs_names.intersection(metadata.index)
    if len(cells)==0:
        print("No cells in intersection of adata object and metadata supplied in info. Please ensure these objects are correct.")
        exit(1)
    adataObj = adataObj[cells, :]

    with open(genesets_json, 'r') as j:
      genesets = json.loads(j.read())
    
    adataObj.obs["RIRA_Immune_v2.cellclass"] = metadata.loc[:,"RIRA_Immune_v2.cellclass"]
    
    genesets = spc_tl.check_gene_set_dictionary(adataObj,
    genesets,
    obs_key='RIRA_Immune_v2.cellclass',
    global_key='global')
    
    model = spc.est_spectra(adata=adataObj,
    gene_set_dictionary=genesets,
    use_highly_variable=True,
    cell_type_key="RIRA_Immune_v2.cellclass",
    use_weights=True,
    lam=0.1, # varies depending on data and gene sets, try between 0.5 and 0.001
    delta=0.001,
    kappa=None,
    rho=0.001,
    use_cell_types=True,
    n_top_vals=50,
    label_factors=True,
    overlap_threshold=0.2,
    clean_gs = True,
    min_gs_num = 3,
    num_epochs=2 #here running only 2 epochs for time reasons, we recommend 10,000 epochs for most datasets
    )
    
    adataObj.uns['SPECTRA_overlap'].to_csv(os.path.join(seuratToAnnDataDir,"SPECTRA_overlap.csv"))
    
    index_labels = adataObj.uns['SPECTRA_overlap'].index
    gene_weights = pd.DataFrame(adataObj.uns['SPECTRA_factors'],
    index= index_labels,columns=adataObj.var[adataObj.var['spectra_vocab']].index)
    gene_weights.to_csv(os.path.join(seuratToAnnDataDir, "geneweights.csv"))
    
    pd.DataFrame(adataObj.uns['SPECTRA_markers'], index= index_labels).to_csv(os.path.join(seuratToAnnDataDir, "SPECTRA_markers.csv"))
    
    with open(os.path.join(seuratToAnnDataDir,"SPECTRA_L.json"), "w") as fp:
      json.dump(adataObj.uns['SPECTRA_L'] , fp)
    
    pd.DataFrame(adataObj.obsm['SPECTRA_cell_scores'], columns=index_labels, index = adataObj.obs.index).to_csv(os.path.join(seuratToAnnDataDir, "SPECTRA_cell_scores.csv"))
    
    
    

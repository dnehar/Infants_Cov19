import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

from matplotlib import rcParams

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.style.use('ggplot')

sc.settings.set_figure_params(transparent=True)

# function to quantify monocytes co-expressing IL1B and ISGs

def dual_exp(df,gene,thresh = [0,0]):
    gene_1 =  pd.Series((df[:,gene[0]].X>thresh[0]).flatten(),index= df.obs.index)
    gene_2 =  pd.Series((df[:,gene[1]].X>thresh[1]).flatten(),index= df.obs.index)
    name = ‘_’.join(gene)
    fin = []
    for i in gene_1.index:
        if gene_1[i]:
            if gene_2[i]:
                fin += [name]
            else:
                fin += [gene[0]]
        else:
            if gene_2[i]:
                fin += [gene[1]]
            else:
                fin += [‘Neither’]
    fin = pd.Series(fin,index=gene_1.index)
    df.obs[name] = fin.astype(‘category’)
    col = name + ‘_colors’
    trans = {gene[0]:‘#FFFF00’,name:‘#FF0000’,gene[1]:‘#31A354’,‘Neither’:‘#D3D3D3’} ##FFA500
    df.uns[col] = [trans[i] for i in df.obs[name].cat.categories]
    return df

# run function 
dual_exp(CD14_mono, gene=['IL1B','ISG15'],thresh = [0.5,0.5]) #CD14_mono anndata object 

# plot umap #########

#plot bar plot 
matplotlib.style.use('default')
rcParams[‘figure.figsize’] = (3,3)
sc.pl.umap(CD14_mono,
           color=[‘IL1B’,‘IFI27’,‘IL1B_IFI27’],
           size=3,
           ncols=3,
           color_map=‘OrRd’,
           #legend_loc=‘on data’,
           frameon=False) #,  save=‘_ISG_Inflam_coExpression_CD14_mono.pdf’)

#-- proportion table #########


Groups_tab = pd.crosstab(index=CD14_mono.obs[‘CD14_mono_Res0_4’],  # Make a crosstab
                        columns=CD14_mono.obs[‘Patient_groups’], margins=True)               # Name the count column
#-- change index and columns order
col= [ “#99D8C9”,“#9ECAE1"]
MyTab= Groups_tab.div(Groups_tab[“All”], axis=0)
MyTab2 = MyTab.drop(columns=“All”)
MyTab2.plot(kind=“bar”,
            figsize=(7,2.5),
            stacked=True,
            linewidth=1,
            width=0.6,
            fontsize=20,
            color=col_pGroups)
plt.title(“Mono-SCs | Res=0.4 | Patient_groups”, fontsize=18)
plt.ylabel(“Prop. of cells”, fontsize=18)
plt.xlabel(“Clusters”, fontsize=18)
plt.ylim=1.0
plt.legend(loc=‘center left’, bbox_to_anchor=(1, 0.5), fontsize=10)
plt.show()

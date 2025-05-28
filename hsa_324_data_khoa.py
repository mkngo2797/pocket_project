#Analysis of human data pre-mir-324 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

plt.rcdefaults()
folder_path = 'D:/1.miRNA/END_randomized project/dataframes'
#Importing summarized dataframes
human_dfA_nop = pd.read_csv(folder_path+'/human_dfA_nop.txt', sep='\t')
human_dfA_pnk = pd.read_csv(folder_path+'/human_dfA_pnk.txt', sep='\t')
human_dfT_nop = pd.read_csv(folder_path+'/human_dfT_nop.txt', sep='\t')
human_dfT_pnk = pd.read_csv(folder_path+'/human_dfT_pnk.txt', sep='\t')
human_dfG_nop = pd.read_csv(folder_path+'/human_dfG_nop.txt', sep='\t')
human_dfG_pnk = pd.read_csv(folder_path+'/human_dfG_pnk.txt', sep='\t')
human_dfC_nop = pd.read_csv(folder_path+'/human_dfC_nop.txt', sep='\t')
human_dfC_pnk = pd.read_csv(folder_path+'/human_dfC_pnk.txt', sep='\t')


#Checking the impact of nucleotide at both ends
list_dfname=['human_dfA_nop','human_dfA_pnk','human_dfT_nop','human_dfT_pnk','human_dfG_nop','human_dfG_pnk', 'human_dfC_nop', 'human_dfC_pnk']
inputdf=[human_dfA_nop,human_dfA_pnk,human_dfT_nop,human_dfT_pnk,human_dfG_nop,human_dfG_pnk, human_dfC_nop, human_dfC_pnk]
i=0
for dfcheck in inputdf: 
    df = dfcheck.copy()
    df = df.pivot(index='Variant',columns='Cleavage_site',values='Mean_Cleavage_accuracy')
    df = df[[20,21,22,23]]
    df.sort_index(ascending=True, inplace=True)
    df.fillna(0, inplace=True)
    df.sort_values([21],ascending=True,inplace=True)
    
    
    ax = plt.figure(figsize=(6,8))
    mpl.rcParams['axes.linewidth'] = 1 #set the value globally
    mpl.rcParams['axes.spines.right'] = True
    mpl.rcParams['axes.spines.top'] = True
    my_color = sns.light_palette("seagreen", as_cmap=True)
    ax = sns.heatmap(data=df,cmap=my_color,vmax=1,vmin=0,
                      cbar_kws={"shrink":1 },cbar=False)
    ax.tick_params(axis='y', width = 0, length=0)
    ax.tick_params(axis='x', width = 0, length=0)
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(visible=False,size=8)
    plt.yticks(visible=False)
    #plt.title(list_dfname[i])
    #plt.savefig(folder_path+'/figures/sorted_heatmap/'+list_dfname[i], dpi=150, bbox_inches='tight')
    plt.show()
    i+=1


#Drawing boxplot to compare 
df_combine_pnk = pd.concat([human_dfA_pnk, human_dfT_pnk, human_dfC_pnk, human_dfG_pnk])
#extract to compare with fly
#df_combine_pnk.to_csv('D:/1.miRNA/END_randomized project/DCR1-DCL1/rawcount_dcr_pnk/df_hsa_pnk_combine.txt', sep='\t', index=False)

df_DC21_pnk = df_combine_pnk[df_combine_pnk['Cleavage_site']==21]
group_order = ['A','T','G','C']
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_DC21_pnk, x='Group', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=group_order)
ax =sns.stripplot(data=df_DC21_pnk, x='Group', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/5pnt_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()
    
df_combine_nop = pd.concat([human_dfA_nop, human_dfT_nop, human_dfC_nop, human_dfG_nop])
df_DC21_nop = df_combine_nop[df_combine_pnk['Cleavage_site']==21]
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_DC21_nop, x='Group', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=group_order)
ax =sns.stripplot(data=df_DC21_nop, x='Group', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
#plt.savefig(folder_path+'/figures/sorted_heatmap/5pnt_nop_box.png', dpi=150, bbox_inches='tight')
plt.show()  
#Checking 3p last nu
df_3pend_pnk = df_DC21_pnk.copy()
df_3pend_nop = df_DC21_nop.copy()
df_3pend_pnk['3pend']=df_3pend_pnk['Randomized_nts'].str[-1:]
df_3pend_nop['3pend']=df_3pend_nop['Randomized_nts'].str[-1:]

#Boxplot the 1st 3p
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_3pend_pnk, x='3pend', y='Mean_Cleavage_accuracy', color="#bbbbbb", showfliers=False, order=group_order)
ax = sns.stripplot(data=df_3pend_pnk, x='3pend', y='Mean_Cleavage_accuracy', color="#aaaaaa",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/3pnt_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()
#%%
#Boxplot the 3p_2nd
group_order = ['A','T','G','C']
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_DC21_pnk, x='3p_2nd', y='Mean_Cleavage_accuracy', color="#bbbbbb", showfliers=False, order=group_order)
ax = sns.stripplot(data=df_DC21_pnk, x='3p_2nd', y='Mean_Cleavage_accuracy', color="#aaaaaa",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False,size=8)
#plt.savefig(folder_path+'/figures/sorted_heatmap/3p_2nd_DC21_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

#Boxplot for the 3p-3rd nucleotide

group_order = ['A','T','G','C']
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_DC21_pnk, x='3p_op', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=group_order)
ax = sns.stripplot(data=df_DC21_pnk, x='3p_op', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False,size=8)
#plt.savefig(folder_path+'/figures/sorted_heatmap/3p_3rd_DC21_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

ax = plt.figure(figsize=(8,5))
ax = sns.boxplot(data=df_3pend_nop, x='3pend', y='Mean_Cleavage_accuracy', color="#bbbbbb", showfliers=False, order=group_order)
ax = sns.stripplot(data=df_3pend_nop, x='3pend', y='Mean_Cleavage_accuracy', color="#aaaaaa",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],labels='DC21_accuracy')
#plt.savefig(folder_path+'/figures/sorted_heatmap/3pnt_nop_box.png', dpi=150, bbox_inches='tight')
plt.show()
#%%
#Drawing heatmap for 256 variants 
df =df_combine_pnk.copy()
df = df.pivot(index='Variant',columns='Cleavage_site',values='Mean_Cleavage_accuracy')
df = df[[20,21,22,23]]
df.sort_index(ascending=True, inplace=True)
df.fillna(0, inplace=True)
    
ax = plt.figure(figsize=(8,24))  #width and height
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True
my_color = sns.light_palette("seagreen", as_cmap=True)
ax = sns.heatmap(data=df,cmap=my_color,vmax=1,vmin=0,
                      cbar_kws={"shrink":1 },cbar=False)
ax.tick_params(axis='y', width = 0, length=0)
ax.tick_params(axis='x', width = 0, length=0)
plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/unsort_heatmap_256.png', dpi=150, bbox_inches='tight')
plt.show()
#%% Checking number of variants obtained per subgroup
df = df_combine_pnk.copy()
df.drop_duplicates(subset=['Group','Variant'],keep='first',inplace=True)
df['Variant_per_group'] = df.groupby(['Group'])['Variant'].transform('count')
df.drop_duplicates(subset='Group',keep='first',inplace=True)    
df = df[['Group', 'Variant_per_group']]

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

ax = plt.figure(figsize=(4,5))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

ax = sns.barplot(data=df, x="Group",y='Variant_per_group',color='#5AC9A1', order = ['A','T','G','C'])


ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'y', color = 'black', linestyle = '--', linewidth = 0.2)
ax.set_axisbelow(True)
plt.ylim(0,64)
plt.yticks([0,20,40,60,64])
plt.xlim(-1,4)
#plt.xticks([1,2,3,4,5,6])
plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/numberofvariants_pnk.png', dpi=150, bbox_inches='tight')
plt.show()

#%% Plotting reproducibility (Manually change the repeat and subgroup)
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

data = human_dfC_pnk.copy()
sample = 'C-PNKed'
data.drop_duplicates(subset=['Variant'], keep='first', inplace=True)


ax = plt.figure(figsize=(5,5))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

replicate_a = 'rep1'
replicate_b = 'rep2'
replicate_c = 'rep3'

choice1 = replicate_c
choice2 = replicate_a

ax = sns.scatterplot(data=data,x='Global_efficiency_'+choice1,y='Global_efficiency_'+choice2,color='#5AC9A1',s=30)


ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(color = 'black', linestyle = '--', linewidth = 0.2)
plt.ylim(-2.2,2.2)
plt.yticks([-2,-1,0,1,2])
plt.xlim(-2.2,2.2)
plt.xticks([-2,-1,0,1,2])
plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)

plt.savefig(folder_path+f'/figures/sorted_heatmap/Pre-mir-324-end-randomization-reproducibility_{choice1}_{choice2}_{sample}.png', dpi=150, bbox_inches='tight')
plt.show()
#pearson correlation analysis
import numpy as np
data1 = data['Global_efficiency_'+choice1].to_numpy()
data2 = data['Global_efficiency_'+choice2].to_numpy()
r = np.corrcoef(data1, data2)
print (r)

#%% Can skip this cell 
#Checking 2nd structure last nu: Comparing 0-2 and 1-3OVH (Run after next cell)
ax = plt.figure(figsize=(6,8))
ax = sns.boxplot(data=df_DC21_pnk, x='5p_flanking_length', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False)
ax =sns.stripplot(data=df_DC21_pnk, x='5p_flanking_length', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True)
plt.xlabel('')
plt.ylabel('')
ax.set_xticklabels([])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
#plt.savefig(folder_path+'/figures/sorted_heatmap/overhang_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

#Further check according to group and structure 
df_checking = df_DC21_pnk.copy()
df_checking['5pnulength'] = df_checking['Group'] + df_checking['5p_flanking_length'].astype(str)
ax = plt.figure(figsize=(14,8))
x_order = ['A0.0','A1.0', 'T0.0', 'T1.0', 'G0.0', 'G1.0', 'C0.0', 'C1.0']

ax = sns.boxplot(data=df_checking, x='5pnulength', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=x_order)
ax =sns.stripplot(data=df_checking, x='5pnulength', y='Mean_Cleavage_accuracy', color="#2CA87B",size=4, jitter=True, order=x_order)
plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
#ax.set_xticklabels(labels=x_order,rotation=45, horizontalalignment='right')
#plt.savefig(folder_path+'/figures/sorted_heatmap/5pnuandOVH.png', dpi=150, bbox_inches='tight')
plt.show()

#Checking the impact of 3rd nu 3p (opposite) - not a good analysis as affected by 5p nu 

df_checking['3p3rdnulength'] = df_checking['3p_op'] + df_checking['5p_flanking_length'].astype(str)
ax = plt.figure(figsize=(12,8))
x_order = ['A0.0','A1.0', 'T0.0', 'T1.0', 'G0.0', 'G1.0', 'C0.0', 'C1.0']

ax = sns.boxplot(data=df_checking, x='3p3rdnulength', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=x_order)
ax =sns.stripplot(data=df_checking, x='3p3rdnulength', y='Mean_Cleavage_accuracy', color="#2CA87B",size=4, jitter=True, order=x_order)
plt.xlabel('')
plt.ylabel('')
ax.set_xticklabels(labels=x_order,rotation=45, horizontalalignment='right')
#plt.savefig(folder_path+'/figures/sorted_heatmap/5pnuandOVH.png', dpi=150, bbox_inches='tight')
plt.show()

#%%
#List all combines in 5p-3p opposite
df_DC21_pnk = df_DC21_pnk.copy()
df_DC21_pnk['3p_op'] =df_DC21_pnk['Randomized_nts'].str[:1]
df_DC21_pnk['pair']=df_DC21_pnk['Group']+df_DC21_pnk['3p_op']

df_DC21_nop = df_DC21_nop.copy()
df_DC21_nop['3p_op'] =df_DC21_nop['Randomized_nts'].str[:1]
df_DC21_nop['pair']=df_DC21_nop['Group']+df_DC21_nop['3p_op']


df_test=df_DC21_pnk.copy()
df_test=df_test[df_test['5p_flanking_length']==1]
ax = plt.figure(figsize=(13,6))
ax = sns.boxplot(data=df_test, x='pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False)
ax =sns.stripplot(data=df_test, x='pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=4, jitter=True)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/all_1-3_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

df_test=df_DC21_pnk.copy()
df_test=df_test[df_test['5p_flanking_length']==0]
ax = plt.figure(figsize=(8,6))
ax = sns.boxplot(data=df_test, x='pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False)
ax =sns.stripplot(data=df_test, x='pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=4, jitter=True)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9])
#plt.savefig(folder_path+'/figures/sorted_heatmap/all_0-2_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()
#Drawing for nop group
df_test=df_DC21_nop.copy()
df_test=df_test[df_test['5p_flanking_length']==1]
ax = plt.figure(figsize=(13,6))
ax = sns.boxplot(data=df_test, x='pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False)
ax =sns.stripplot(data=df_test, x='pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=4, jitter=True)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
#plt.savefig(folder_path+'/figures/sorted_heatmap/all_1-3_nop_box.png', dpi=150, bbox_inches='tight')
plt.show()

df_test=df_DC21_nop.copy()
df_test=df_test[df_test['5p_flanking_length']==0]
ax = plt.figure(figsize=(8,6))
ax = sns.boxplot(data=df_test, x='pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False)
ax =sns.stripplot(data=df_test, x='pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=4, jitter=True)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9])
#plt.savefig(folder_path+'/figures/sorted_heatmap/all_0-2_nop_box.png', dpi=150, bbox_inches='tight')
plt.show()

#%% Checking impact of combination of 5p with 2nd and 1st nu 
#Must run: calling part
df_DC21_pnk = df_DC21_pnk.copy()
df_DC21_pnk['3p_1st'] =df_DC21_pnk['Randomized_nts'].str[-1]
df_DC21_pnk['3p_2nd'] =df_DC21_pnk['Randomized_nts'].str[1]
df_DC21_pnk['5p3p1st_pair']=df_DC21_pnk['Group']+df_DC21_pnk['3p_1st']
df_DC21_pnk['5p3p2nd_pair']=df_DC21_pnk['Group']+df_DC21_pnk['3p_2nd']

df_testing=df_DC21_pnk.copy()
bar_order = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
ax = plt.figure(figsize=(18,6))
ax = sns.boxplot(data=df_testing, x='5p3p1st_pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=bar_order)
ax = sns.stripplot(data=df_testing, x='5p3p1st_pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order =bar_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/5p3p1st_pair_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

df_testing=df_DC21_pnk.copy()
bar_order = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
ax = plt.figure(figsize=(18,6))
ax = sns.boxplot(data=df_testing, x='5p3p2nd_pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=bar_order)
ax =sns.stripplot(data=df_testing, x='5p3p2nd_pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order =bar_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/5p3p2nd_pair_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

df_testing=df_DC21_pnk.copy()
bar_order = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
ax = plt.figure(figsize=(18,6))
ax = sns.boxplot(data=df_testing, x='pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=bar_order)
ax =sns.stripplot(data=df_testing, x='pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order =bar_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/5p3p3rd_pair_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()


#Save the df for merging with fly Dicer
df_DC21_pnk.to_csv('D:/1.miRNA/END_randomized project/dataframes/human_df_DC21_pnk.txt', sep ='\t', index =False)
#%% Checking the enrichment of 3p strand according to DC21 and 5p nucleotide
#This is a separated cell

df_working = df_DC21_pnk.copy()
df_top =pd.DataFrame()
def extracting_top6(group_5p): 
    df_using =df_working.copy()
    df_using = df_using[df_using['Group']==group_5p]
    threshold = df_using['Mean_Cleavage_accuracy'].quantile(0.80)
    top = df_using[df_using['Mean_Cleavage_accuracy'] >= threshold]
    return(top)

for group in df_working['Group'].unique():
    top_in_group = extracting_top6(group)
    df_top = pd.concat([df_top,top_in_group], axis=0)

ran_top=df_top['Randomized_nts']
#Save ran_top as fasta file                    
with open('D:/1.miRNA/END_randomized project/dataframes/top10_DC21.fasta', 'w') as f:
    for index, sequence in ran_top.items():
        f.write(f'>variant_{index}\n{sequence}\n')                      
top20_list =df_top['Randomized_nts'].tolist()
top20_DC21_list =[sequence.replace('T','U') for sequence in top20_list]

#Draw logo plot

import logomaker as lm
import matplotlib.pyplot as plt

def draw_weblogo (list_sequence, save_fig = 'no'):
    
    
#    list_sequence = ['AACGC', 'ACACT']
    
    # counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat['correct_index'] = counts_mat.index.map(lambda x: x+1)
    counts_mat = counts_mat.set_index('correct_index')
    
    crp_logo  = lm.Logo(counts_mat, 
                        figsize = [2*len(counts_mat), 4],  
                        color_scheme = {'A':'#89CFF0' , 'C':"#5AC9A1" ,  'G':'#9B111E' , 'U':'#3da4dc','T':'#3da4dc'} ,
                        font_name='Arial Rounded MT Bold', zorder = 3)
    
    for _, spine in crp_logo.ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_color('black')
    
    plt.yticks([0,0.5,1], fontsize = 0, color = 'white')
    plt.xticks( fontsize = 0, color = 'white')
    
    if save_fig != 'no':
        # plt.title(save_fig)
        # plt.show()
        plt.savefig(folder_path+'/figures/sorted_heatmap/3pstrand_DC21_top20per.png', bbox_inches="tight", dpi =1500)
    else:
        plt.show()
    return()

draw_weblogo(top20_DC21_list, save_fig ='yes')


#Draw top 20% DC22 variant enrichment (Run after get the df_DC22_pnk)
df_DC22_working = df_DC22_pnk.copy()
df_top =pd.DataFrame()

def extracting_top20(group_5p): 
    df_using =df_DC22_working.copy()
    df_using = df_using[df_using['Group']==group_5p]
    threshold = df_using['Mean_Cleavage_accuracy'].quantile(0.80)
    top = df_using[df_using['Mean_Cleavage_accuracy'] >= threshold]
    return(top)

for group in df_DC22_working['Group'].unique():
    top_in_group = extracting_top20(group)
    df_top = pd.concat([df_top,top_in_group], axis=0)

ran_top=df_top['Randomized_nts']
top20_DC22_list =df_top['Randomized_nts'].tolist()
top20_DC22_list =[sequence.replace('T','U') for sequence in top20_DC22_list]

def draw_weblogo (list_sequence, save_fig = 'no'):
    
    
#    list_sequence = ['AACGC', 'ACACT']
    
    # counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat = lm.alignment_to_matrix(list_sequence, to_type = 'information', pseudocount = 0)
    counts_mat['correct_index'] = counts_mat.index.map(lambda x: x+1)
    counts_mat = counts_mat.set_index('correct_index')
    
    crp_logo  = lm.Logo(counts_mat, 
                        figsize = [2*len(counts_mat), 4],  
                        color_scheme = {'A':'#89CFF0' , 'C':"#5AC9A1" ,  'G':'#9B111E' , 'U':'#3da4dc','T':'#3da4dc'} ,
                        font_name='Arial Rounded MT Bold', zorder = 3)
    
    for _, spine in crp_logo.ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_color('black')
    
    plt.yticks([0,0.5,1], fontsize = 0, color = 'white')
    plt.xticks( fontsize = 0, color = 'white')
    
    if save_fig != 'no':
        # plt.title(save_fig)
        # plt.show()
        plt.savefig(folder_path+'/figures/sorted_heatmap/3pstrand_DC22_top20per.png', bbox_inches="tight", dpi =1500)
    else:
        plt.show()
    return()

draw_weblogo(top20_DC22_list, save_fig ='yes')

#%%Drawing DC22 plot for each subgroup
#This cell for DC22 analysis
import math
def calculate_log2(x):
        return math.log2(x)


df_DC22_pnk = df_combine_pnk[df_combine_pnk['Cleavage_site']==22]
group_order = ['A','T','G','C']
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_DC22_pnk, x='Group', y='Mean_Cleavage_accuracy', color="#f94449", showfliers=False, order=group_order)
ax =sns.stripplot(data=df_DC22_pnk, x='Group', y='Mean_Cleavage_accuracy', color="#d1001f",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/DC22_5pnt_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

#Calculating ratio
df_using_DC21 = df_DC21_pnk.copy()
df_using_DC22 = df_DC22_pnk.copy()
df_using_DC21.rename(columns={'Mean_Cleavage_accuracy':'DC21_accuracy'}, inplace=True)
df_using_DC22.rename(columns={'Mean_Cleavage_accuracy':'DC22_accuracy'}, inplace=True)
df_using_merged = pd.merge(df_using_DC21, df_using_DC22[['Variant', 'DC22_accuracy']], on='Variant', how='right')
df_using_merged['ratio']=df_using_merged['DC22_accuracy']/df_using_merged['DC21_accuracy']
df_using_merged['log2_ratio']=df_using_merged['ratio'].apply(calculate_log2)
df_ratio = df_using_merged[['Variant','Group', 'Randomized_nts', 'ratio', 'log2_ratio']]

#Draw boxplot for log2(ratioDC22/DC21):
group_order = ['A','T','G','C']
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_ratio, x='Group', y='log2_ratio', color="#ca5cdd", showfliers=False, order=group_order)
ax = sns.stripplot(data=df_ratio, x='Group', y='log2_ratio', color="#a000c8",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([-4,-3,-2,-1,0,1,2,3,4])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(folder_path+'/figures/sorted_heatmap/ratio2221_5p_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

#Draw boxplot for log2(ratioDC21/DC22)

df_ratio21on22 = df_using_merged.copy()
df_ratio21on22['ratio'] = df_ratio21on22['DC21_accuracy']/df_ratio21on22['DC22_accuracy']
df_ratio21on22['log2_ratio']=df_ratio21on22['ratio'].apply(calculate_log2)

group_order = ['A','T','G','C']
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_ratio21on22, x='Group', y='log2_ratio', color="#ca5cdd", showfliers=False, order=group_order)
ax = sns.stripplot(data=df_ratio21on22, x='Group', y='log2_ratio', color="#a000c8",size=5, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([-4,-3,-2,-1,0,1,2,3,4])
ax = plt.xticks(visible=False, size=8)
ax = plt.yticks(visible=False, size=8)
plt.savefig(folder_path+'/figures/sorted_heatmap/ratio21on22_5p_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()











#Draw line plot with sns
# Plot the lines for each group using seaborn

#Draw the lineplot for DC21 accuracy 
color_scheme = {'A':'#f08080' , 'C':"#5AC9A1" ,  'G':'lightgrey' , 'U':'#3da4dc','T':'#3da4dc'}
ax = plt.figure(figsize=(24,10))
ax = sns.set(font_scale=2.5)
ax = sns.set_style("whitegrid")
ax = sns.lineplot(data=df_using_merged, x='Randomized_nts', y='DC21_accuracy', hue='Group', linewidth=5, palette=color_scheme)
ax = plt.xlabel('')
ax = plt.ylabel('')
ax = plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
ax = plt.xticks(visible=False,size=8)
ax = plt.yticks(visible=False)
ax = plt.xlim(-1, 64)
ax = plt.title('')
ax = plt.legend().remove()
#plt.savefig(folder_path+'/figures/sorted_heatmap/DC21_accuracy_line.png', dpi=150, bbox_inches='tight')
plt.show()

#Draw the lineplot for ln_ratio
ax = plt.figure(figsize=(24,10))
ax = sns.set(font_scale=2.5)
ax = sns.set_style("whitegrid")
ax = sns.lineplot(data=df_ratio, x='Randomized_nts', y='log2_ratio', hue='Group', linewidth=5)
ax = plt.xlabel('')
ax = plt.ylabel('')
ax = plt.yticks([-4,-3,-2,-1,0,1,2,3,4])
ax = plt.xticks(visible=False,size=8)
ax = plt.yticks(visible=False)
ax = plt.xlim(-1, 64)
ax = plt.title('')
ax = plt.legend().remove()
#plt.savefig(folder_path+'/figures/sorted_heatmap/ln_ratio_line.png', dpi=150, bbox_inches='tight')
plt.show()

#Draw color-based figure 16 combinations last pair for analyse both G/C/OVH

# Define a color palette for each box
color_palette = ['#f94449', '#9B111E', '#89CFF0', '#5AC9A1', '#f94449', '#9B111E', '#89CFF0', '#5AC9A1',
                 '#f94449', '#9B111E', '#89CFF0', '#89CFF0', '#f94449', '#f94449', '#89CFF0', '#5AC9A1']

bar_order = ['AA', 'AG', 'AT', 'AC', 'TT', 'TG', 'TA', 'TC', 'GA', 'GG', 'GT', 'GC', 'CA', 'CT', 'CG', 'CC']
ax = plt.figure(figsize=(18,6))
ax = sns.boxplot(data=df_DC21_pnk, x='pair', y='Mean_Cleavage_accuracy', palette=color_palette, showfliers=False, order=bar_order)
ax = sns.stripplot(data=df_DC21_pnk, x='pair', y='Mean_Cleavage_accuracy', palette=color_palette,size=5, jitter=True, order =bar_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.savefig(folder_path+'/figures/sorted_heatmap/lastpair_colored_G_C_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()


red #f94449
light blue #89CFF0
green  #5AC9A1

#%% Visualize Lee 2023 variantsfor discussion figure
df_lee = df_DC21_pnk[(df_DC21_pnk['3p_op'] =='A') & (df_DC21_pnk['3p_1st']=='T')& (df_DC21_pnk['3p_2nd']=='T')]
group_order = ['A','T','G','C']
ax = plt.figure(figsize=(10,8))
ax = sns.boxplot(data=df_DC21_pnk, x='Group', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=group_order)
ax =sns.stripplot(data=df_DC21_pnk, x='Group', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order=group_order)
ax =sns.stripplot(data=df_lee, x='Group', y='Mean_Cleavage_accuracy', color="#C21807",size=12, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.savefig(folder_path+'/figures/sorted_heatmap/compare_Lee_5pnt.png', dpi=150, bbox_inches='tight')
plt.show()

#%% Redo plotting for figure in dual pocket manuscript

def boxplot_DC21_DC22(): 
    df_using = df_combine_pnk[(df_combine_pnk['Cleavage_site'].isin([21, 22]))]
    group_order = ['A', 'G', 'T', 'C']
    box_color = {21:'#5AC9A1', 22:'#f94449'}
    dot_color = {21:'#358856', 22:'#d1001f'}

    plt.figure(figsize=(12, 8))
    sns.set_style("white")
    ax = sns.boxplot(data=df_using, x='Group', y='Mean_Cleavage_accuracy', hue='Cleavage_site', palette=box_color, showfliers=False, order=group_order, legend=False)
    ax = sns.stripplot(data=df_using, x='Group', y='Mean_Cleavage_accuracy', hue='Cleavage_site',palette=dot_color,size=7, jitter=True, dodge=True, order=group_order, legend=False)
    plt.xlabel('', fontsize=12)
    plt.ylabel('', fontsize=12)
    plt.yticks([0,0.2, 0.4, 0.6, 0.8, 1.0])
    plt.xticks(visible=False,size=8)
    plt.yticks(visible=False)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    
    plt.tick_params(left = True, bottom = True, right = False , labelleft = False,labelbottom = False) 
    plt.savefig(folder_path+'/figures/sorted_heatmap/5pnt_box_DC21_DC22.png', dpi=150, bbox_inches='tight')
    plt.show()
    return
boxplot_DC21_DC22()

#Draw boxplot_DC21_DC22 for the last nucleotide 3p end
def boxplot_DC21_DC22_3p_end(): 
    df_using = df_combine_pnk[(df_combine_pnk['Cleavage_site'].isin([21, 22]))]
    df_using.loc[:, '3p-end'] = df_using['Randomized_nts'].str[-1]
    group_order = ['A', 'G', 'T', 'C']
    box_color = {21:'#5AC9A1', 22:'#f94449'}
    dot_color = {21:'#358856', 22:'#d1001f'}

    plt.figure(figsize=(12, 8))
    sns.set_style("white")
    ax = sns.boxplot(data=df_using, x='3p-end', y='Mean_Cleavage_accuracy', hue='Cleavage_site', palette=box_color, showfliers=False, order=group_order, legend=False)
    ax = sns.stripplot(data=df_using, x='3p-end', y='Mean_Cleavage_accuracy', hue='Cleavage_site',palette=dot_color,size=7, jitter=True, dodge=True, order=group_order, legend=False)
    plt.xlabel('', fontsize=12)
    plt.ylabel('', fontsize=12)
    plt.yticks([0,0.2, 0.4, 0.6, 0.8, 1.0])
    plt.xticks(visible=False,size=8)
    plt.yticks(visible=False)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    
    plt.tick_params(left = True, bottom = True, right = False , labelleft = False,labelbottom = False) 
    plt.savefig(folder_path+'/figures/sorted_heatmap/3pnt_end_box_DC21_DC22.png', dpi=150, bbox_inches='tight')
    
    plt.show()
    return
boxplot_DC21_DC22_3p_end()


#ax =sns.stripplot(data=df_using, x='Group', y='Mean_Cleavage_accuracy', hue='Cleavage_site',color="#2CA87B",size=5, jitter=True, order=group_order)
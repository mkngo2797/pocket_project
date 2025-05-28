import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

path ='D:/1.miRNA/END_randomized project/DCR1-DCL1/rawcount_dcr_pnk/'
#%% Processing from raw_count files to df
#Making control df
df_using = pd.read_csv(path+'control_pnk_rawcount.bed', sep=' ', names =['Sequencing_read', 'Raw_count', 'Group', 'Randomized_nts'])
df_structure = pd.read_csv(path+'Pre-mir-324-end-randomization-structure.bed', sep='\t')
del df_using['Sequencing_read']
df_using['Sum_raw_count'] = df_using.groupby(['Group','Randomized_nts'])['Raw_count'].transform('sum')
del df_using['Raw_count']
df_control = df_structure.merge(df_using, on=['Group','Randomized_nts'], how='left')
df_control['RPM_control']=df_control['Sum_raw_count']/df_control['Sum_raw_count'].sum()*1000000
df_control_A = df_control[df_control['Group']=='A']
df_control_T = df_control[df_control['Group']=='T']
df_control_G = df_control[df_control['Group']=='G']
df_control_C = df_control[df_control['Group']=='C']






#Inputting product reads 
df_clv_A_rep1 = pd.read_csv(path+'DCR_A_RP1_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_A_rep2 = pd.read_csv(path+'DCR_A_RP2_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_A_rep3 = pd.read_csv(path+'DCR_A_RP3_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_T_rep1 = pd.read_csv(path+'DCR_T_RP1_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_T_rep2 = pd.read_csv(path+'DCR_T_RP2_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_T_rep3 = pd.read_csv(path+'DCR_T_RP3_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_G_rep1 = pd.read_csv(path+'DCR_G_RP1_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_G_rep2 = pd.read_csv(path+'DCR_G_RP2_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_G_rep3 = pd.read_csv(path+'DCR_G_RP3_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_C_rep1 = pd.read_csv(path+'DCR_C_RP1_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_C_rep2 = pd.read_csv(path+'DCR_C_RP2_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])
df_clv_C_rep3 = pd.read_csv(path+'DCR_C_RP3_rawcount.bed',sep=' ',names=['Sequencing_read','Raw_count', 'Randomized_nts'])

import math
def merging_clv_ctrl(df_clv, df_con, repeat):
    df1=df_clv.copy()
    df2=df_con.copy()
    df2=df2.merge(df1, on=['Randomized_nts'], how='inner')
    df2.reset_index(inplace=True, drop=True)
    df2['Cleavage_site']=df2['Sequencing_read'].apply(len)
    df2.dropna(inplace=True)
    df2['RPM_CP'+repeat] = df2['Raw_count'] / df2['Raw_count'].sum() * 1000000
    df2['Sum_clv_of_variant'+repeat] = df2.groupby(['Variant'])['RPM_CP'+repeat].transform('sum')
    df2['Positional_efficiency'+repeat] = ((df2["RPM_CP"+repeat]+0.1).apply(math.log2) - (df2["RPM_control"]+0.1).apply(math.log2))
    df2['Global_efficiency'+repeat] = ((df2["Sum_clv_of_variant"+repeat]+0.1).apply(math.log2) - (df2["RPM_control"]+0.1).apply(math.log2))
    df2['Cleavage_accuracy'+repeat] = df2['RPM_CP'+repeat] / df2['Sum_clv_of_variant'+repeat]
    del df2["Sum_clv_of_variant"+repeat]
    del df2["Raw_count"]
    return (df2)
df_A_rep1= merging_clv_ctrl(df_clv_A_rep1, df_control_A, '_rep1')
df_A_rep2= merging_clv_ctrl(df_clv_A_rep2, df_control_A, '_rep2')
df_A_rep3= merging_clv_ctrl(df_clv_A_rep3, df_control_A, '_rep3')

df_T_rep1= merging_clv_ctrl(df_clv_T_rep1, df_control_T, '_rep1')
df_T_rep2= merging_clv_ctrl(df_clv_T_rep2, df_control_T, '_rep2')
df_T_rep3= merging_clv_ctrl(df_clv_T_rep3, df_control_T, '_rep3')

df_G_rep1= merging_clv_ctrl(df_clv_G_rep1, df_control_G, '_rep1')
df_G_rep2= merging_clv_ctrl(df_clv_G_rep2, df_control_G, '_rep2')
df_G_rep3= merging_clv_ctrl(df_clv_G_rep3, df_control_G, '_rep3')

df_C_rep1= merging_clv_ctrl(df_clv_C_rep1, df_control_C, '_rep1')
df_C_rep2= merging_clv_ctrl(df_clv_C_rep2, df_control_C, '_rep2')
df_C_rep3= merging_clv_ctrl(df_clv_C_rep3, df_control_C, '_rep3')

#Merging 3 repeats into one 
from functools import reduce
def merge_dataframe(df_input1, df_input2, df_input3):
    data_frames = [df_input1, df_input2, df_input3]
    df = reduce(lambda  left,right: pd.merge(left,right,on=['Variant','Group','Randomized_nts','Pre_miRNA_sequence','New_define_structure_1',
                                                            'New_define_structure_2', 'concrete_struct', '5p_flanking_length', 'RPM_control',
                                                            'Sequencing_read', 'Cleavage_site'],
                                                how='inner'), data_frames)
    df['Mean_Position_efficiency'] = df[['Positional_efficiency_rep1', 'Positional_efficiency_rep2', 'Positional_efficiency_rep3']].mean(axis=1)
    df['Mean_Global_efficiency'] = df[['Global_efficiency_rep1', 'Global_efficiency_rep2', 'Global_efficiency_rep3']].mean(axis=1)
    df['Mean_Cleavage_accuracy'] = df[['Cleavage_accuracy_rep1', 'Cleavage_accuracy_rep2', 'Cleavage_accuracy_rep3']].mean(axis=1)
    return (df)

df_A=merge_dataframe(df_A_rep1, df_A_rep2, df_A_rep3)
df_T=merge_dataframe(df_T_rep1, df_T_rep2, df_T_rep3)
df_G=merge_dataframe(df_G_rep1, df_G_rep2, df_G_rep3)
df_C=merge_dataframe(df_C_rep1, df_C_rep2, df_C_rep3)

df_dcr_combine = pd.concat([df_A, df_T, df_G, df_C])

#Output the combine dcr df to text file
#df_dcr_combine.to_csv(path+'df_dcr_pnk_combine.txt', sep='\t', index=False)

#%% Processing dcr similar as human enzyme 

#Drawing number of variant obtained
df = df_dcr_combine.copy()
df.drop_duplicates(subset=['Group','Variant'],keep='first',inplace=True)
df['Variant_per_group'] = df.groupby(['Group'])['Variant'].transform('count')
df.drop_duplicates(subset='Group',keep='first',inplace=True)    
df = df[['Group', 'Variant_per_group']]


ax = plt.figure(figsize=(4,5))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

ax = sns.barplot(data=df, x="Group",y='Variant_per_group',color='#FFAE42', order = ['A','T','G','C'])


ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'y', color = 'black', linestyle = '--', linewidth = 0.2)
ax.set_axisbelow(True)
plt.ylim(0,64)
plt.yticks([0,20,40,60,64])
plt.xlim(-1,4)
plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/dcr_numberofvariants_pnk.png', dpi=150, bbox_inches='tight')
plt.show()

#Drawing heatmap for 256 variants 
df_extract=df_dcr_combine[df_dcr_combine['Cleavage_site'].isin([20,21,22,23])]

df =df_extract.copy()
df = df.pivot(index='Variant',columns='Cleavage_site',values='Mean_Cleavage_accuracy')
df = df[[20,21,22,23]]
df.sort_index(ascending=True, inplace=True)
df.fillna(0, inplace=True)
    
ax = plt.figure(figsize=(8,24))  #width and height
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True
my_color = sns.color_palette("YlOrBr", as_cmap=True)
ax = sns.heatmap(data=df,cmap=my_color,vmax=1,vmin=0,
                      cbar_kws={"shrink":1 },cbar=False)
ax.tick_params(axis='y', width = 0, length=0)
ax.tick_params(axis='x', width = 0, length=0)
plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/DCR_unsort_heatmap_256.png', dpi=150, bbox_inches='tight')
plt.show()
#%%
#Mannually check the reproducibility
data = df_C.copy()
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
choice2 = replicate_b

ax = sns.scatterplot(data=data,x='Global_efficiency_'+choice1,y='Global_efficiency_'+choice2,color='#FFAE42',s=30)


ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(color = 'black', linestyle = '--', linewidth = 0.2)
plt.ylim(-2,8)
plt.yticks([-1,1,3,5,7])
plt.xlim(-2,8)
plt.xticks([-1,1,3,5,7])
plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)

plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/'+f'dcr_reproducibility_{choice1}_{choice2}_{sample}.png', dpi=150, bbox_inches='tight')
plt.show()
#pearson correlation analysis
import numpy as np
data1 = data['Global_efficiency_'+choice1].to_numpy()
data2 = data['Global_efficiency_'+choice2].to_numpy()
r = np.corrcoef(data1, data2)
print (r)
#%% Input human processed data 
df_hsa_combine = pd.read_csv(path+'df_hsa_pnk_combine.txt', sep='\t')

#Drawing lineplot for all the variants of fly and human
df_dme_DC22 = df_dcr_combine[df_dcr_combine['Cleavage_site']==22]
data=df_dme_DC22.copy()
data=data.reset_index()
data=data.drop('index', axis=1)

df_hsa_DC22 = df_hsa_combine[df_hsa_combine['Cleavage_site']==22]
data2=df_hsa_DC22.copy()
data2=data2.reset_index()
data2=data2.drop('index', axis=1)


ax = plt.figure(figsize=(40,8))
group_order = ['A','T','G','C']
ax = sns.lineplot(data=data, x='Variant', y='Mean_Cleavage_accuracy', color="#FFAE42")
ax = sns.lineplot(data=data2, x='Variant', y='Mean_Cleavage_accuracy', color="#5AC9A1")

plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=True,size=8)
plt.yticks(visible=True)
#plt.savefig(folder_path+'/figures/sorted_heatmap/fly_5pnt_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()

#%% Boxplot for each nucleotide
df_dme_DC21=df_dcr_combine[df_dcr_combine['Cleavage_site']==21]
df_dme_DC21=df_dme_DC21.reset_index()
df_dme_DC21=df_dme_DC21.drop('index', axis=1)
df_dme_DC21['3p_1st']=df_dme_DC21['Randomized_nts'].str[-1]
df_dme_DC21['3p_2nd']=df_dme_DC21['Randomized_nts'].str[1]
df_dme_DC21['3p_3rd']=df_dme_DC21['Randomized_nts'].str[0]
df_dme_DC21['enzyme']='dcr'
df_dme_DC21['last_pair']=df_dme_DC21['Group']+df_dme_DC21['3p_3rd']


df_hsa_DC21=df_hsa_combine[df_hsa_combine['Cleavage_site']==21]
df_hsa_DC21=df_hsa_DC21.reset_index()
df_hsa_DC21=df_hsa_DC21.drop('index', axis=1)
df_hsa_DC21['3p_1st']=df_hsa_DC21['Randomized_nts'].str[-1]
df_hsa_DC21['3p_2nd']=df_hsa_DC21['Randomized_nts'].str[1]
df_hsa_DC21['3p_3rd']=df_hsa_DC21['Randomized_nts'].str[0]
df_hsa_DC21['enzyme']='hdicer'
df_hsa_DC21['last_pair']=df_hsa_DC21['Group']+df_hsa_DC21['3p_3rd']

df_dme_DC22=df_dme_DC22.reset_index()
df_dme_DC22=df_dme_DC22.drop('index', axis=1)
df_dme_DC22['3p_1st']=df_dme_DC22['Randomized_nts'].str[-1]
df_dme_DC22['3p_2nd']=df_dme_DC22['Randomized_nts'].str[1]
df_dme_DC22['3p_3rd']=df_dme_DC22['Randomized_nts'].str[0]
df_dme_DC22['enzyme']='dcr'
df_dme_DC22['last_pair']=df_dme_DC22['Group']+df_dme_DC22['3p_3rd']

df_hsa_DC22=df_hsa_DC22.reset_index()
df_hsa_DC22=df_hsa_DC22.drop('index', axis=1)
df_hsa_DC22['3p_1st']=df_hsa_DC22['Randomized_nts'].str[-1]
df_hsa_DC22['3p_2nd']=df_hsa_DC22['Randomized_nts'].str[1]
df_hsa_DC22['3p_3rd']=df_hsa_DC22['Randomized_nts'].str[0]
df_hsa_DC22['enzyme']='hdicer'
df_hsa_DC22['last_pair']=df_hsa_DC22['Group']+df_hsa_DC22['3p_3rd']


df_two_species = pd.concat([df_dme_DC21, df_hsa_DC21])
#plotting
mpl.rcParams['axes.linewidth'] = 1
ax = plt.figure(figsize=(10,8))
group_order = ['A','T','G','C']
ax = sns.boxplot(data=df_dme_DC21, x='Group', y='Mean_Cleavage_accuracy', color="#FFAE42", showfliers=False, order=group_order)
ax =sns.stripplot(data=df_dme_DC21, x='Group', y='Mean_Cleavage_accuracy', color="#C36721",size=4, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.tick_params(left = True, bottom = True, right = False , labelleft = False , 
                labelbottom = False) 
plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/dcr_5pnu_box.png', dpi=150, bbox_inches='tight')
plt.show()

#the 1st in 3p in grey color
ax = plt.figure(figsize=(10,8))
group_order = ['A','T','G','C']
ax = sns.boxplot(data=df_dme_DC21, x='3p_1st', y='Mean_Cleavage_accuracy', color="#bbbbbb", showfliers=False, order=group_order)
ax =sns.stripplot(data=df_dme_DC21, x='3p_1st', y='Mean_Cleavage_accuracy', color="#aaaaaa",size=4, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.tick_params(left = True, bottom = True, right = False , labelleft = False , 
                labelbottom = False) 
plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/dcr_3p_1st_box.png', dpi=150, bbox_inches='tight')
plt.show()

#the 2nd in grey
ax = plt.figure(figsize=(10,8))
group_order = ['A','T','G','C']
ax = sns.boxplot(data=df_dme_DC21, x='3p_2nd', y='Mean_Cleavage_accuracy', color="#bbbbbb", showfliers=False, order=group_order)
ax =sns.stripplot(data=df_dme_DC21, x='3p_2nd', y='Mean_Cleavage_accuracy', color="#aaaaaa",size=4, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.tick_params(left = True, bottom = True, right = False , labelleft = False , 
                labelbottom = False) 
plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/dcr_3p_2nd_box.png', dpi=150, bbox_inches='tight')
plt.show()

#the 3rd in 3p
ax = plt.figure(figsize=(10,8))
group_order = ['A','T','G','C']
ax = sns.boxplot(data=df_dme_DC21, x='3p_3rd', y='Mean_Cleavage_accuracy', color="#FFAE42", showfliers=False, order=group_order)
ax =sns.stripplot(data=df_dme_DC21, x='3p_3rd', y='Mean_Cleavage_accuracy', color="#C36721",size=4, jitter=True, order=group_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.tick_params(left = True, bottom = True, right = False , labelleft = False , 
                labelbottom = False) 
plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/dcr_3p_3rd_box.png', dpi=150, bbox_inches='tight')
plt.show()

#Plot overall Cleavage accuracy of human vs fly 
ax = plt.figure(figsize=(10,8))
enzyme_order = ['hdicer', 'dcr']
colors_box = ["#5AC9A1",'#FFAE42' ]
colors_strip =["#2CA87B",'#C36721']
ax = sns.boxplot(data=df_two_species, x='enzyme', y='Mean_Cleavage_accuracy', palette=colors_box, showfliers=False, order=enzyme_order)
ax =sns.stripplot(data=df_two_species, x='enzyme', y='Mean_Cleavage_accuracy', palette=colors_strip,size=4, jitter=True, order=enzyme_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.tick_params(left = True, bottom = True, right = False , labelleft = False , 
                labelbottom = False) 
plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/overall_humanvsfly.png', dpi=150, bbox_inches='tight')
plt.show()


#Calculating p value
import scipy.stats as stats
hdicer_data = df_two_species[df_two_species['enzyme'] == 'hdicer']['Mean_Cleavage_accuracy']
dcr_data = df_two_species[df_two_species['enzyme'] == 'dcr']['Mean_Cleavage_accuracy']

# Perform the Mann-Whitney U test
statistic, p_value = stats.mannwhitneyu(hdicer_data, dcr_data)
print("Mann-Whitney U Test:")
print(f"Statistic: {statistic}")
print(f"p-value: {p_value}")

#%%Check enrichment of 3p strand for DC21 (mannualy set criteria and data)
df_working = df_hsa_DC21.copy()
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

top10_list =df_top['Randomized_nts'].tolist()

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
                        color_scheme = {'A':'#f08080' , 'C':"#5AC9A1" ,  'G':'lightgrey' , 'U':'#3da4dc','T':'#3da4dc'} ,
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
        plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/hsa_top20_3p_logo.png', bbox_inches="tight", dpi =300)
    else:
        plt.show()
    return()

draw_weblogo(top10_list, save_fig ='yes')

#%% coordination 

#Draw box plot for the coordination
bar_order = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
ax = plt.figure(figsize=(18,6))
ax = sns.boxplot(data=df_dme_DC21, x='last_pair', y='Mean_Cleavage_accuracy', color="#FFAE42", showfliers=False, order=bar_order)
ax =sns.stripplot(data=df_dme_DC21, x='last_pair', y='Mean_Cleavage_accuracy', color="#C36721",size=5, jitter=True, order =bar_order)
ax = sns.boxplot(data=df_hsa_DC21, x='last_pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=bar_order)
ax =sns.stripplot(data=df_hsa_DC21, x='last_pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order =bar_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/human_fly_OVL_5p3p3rd_pair_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()



#Draw box plot for the coordination of DC22 accuracy
bar_order = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
ax = plt.figure(figsize=(18,6))
ax = sns.boxplot(data=df_dme_DC22, x='last_pair', y='Mean_Cleavage_accuracy', color="#FFAE42", showfliers=False, order=bar_order)
ax =sns.stripplot(data=df_dme_DC22, x='last_pair', y='Mean_Cleavage_accuracy', color="#C36721",size=5, jitter=True, order =bar_order)
ax = sns.boxplot(data=df_hsa_DC22, x='last_pair', y='Mean_Cleavage_accuracy', color="#5AC9A1", showfliers=False, order=bar_order)
ax =sns.stripplot(data=df_hsa_DC22, x='last_pair', y='Mean_Cleavage_accuracy', color="#2CA87B",size=5, jitter=True, order =bar_order)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/human_fly_OVL_5p3p3rd_pair_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()
#%% Working with df_two species 
df_two_species['5p_enzyme']=df_two_species['Group'] + df_two_species['enzyme']

ax = plt.figure(figsize=(14,8))
#enzyme_order = ['hdicer', 'dcr']
palette_box = {'hdicer':"#5AC9A1",'dcr':'#FFAE42'}
palette_strip ={'hdicer':"#2CA87B",'dcr':'#C36721'}
data=df_two_species.copy()
ax = sns.boxplot(data=data, x='Group', y='Mean_Cleavage_accuracy', showfliers=False, hue ='enzyme', palette=palette_box)
ax = sns.stripplot(data=data, x='Group', y='Mean_Cleavage_accuracy', size=4, jitter=True, hue ='enzyme', palette=palette_strip, dodge=True)
plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.tick_params(left = True, bottom = True, right = False , labelleft = False , 
                labelbottom = False) 
ax = plt.legend().remove()
plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/5p_group_humanvsfly.png', dpi=150, bbox_inches='tight')
plt.show()


#%% Draw the lineplot for DC22 for each replicate of DC22 DCR-1 
df_data=df_dme_DC22.copy()
ax = plt.figure(figsize=(22,12))
ax = sns.lineplot(data=df_data, x='Variant', y='Cleavage_accuracy_rep1', color="#A20021", linewidth=2)
ax = sns.lineplot(data=df_data, x='Variant', y='Cleavage_accuracy_rep2', color="#FED000", linewidth=2)
ax = sns.lineplot(data=df_data, x='Variant', y='Cleavage_accuracy_rep3', color="#00ab41", linewidth=2)

plt.xlabel('')
plt.ylabel('')
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
plt.xticks(visible=True,size=8)
plt.yticks(visible=True)
plt.tick_params(left = True, bottom = True, right = False , labelleft = False , 
                labelbottom = False) 
#plt.savefig(folder_path+'/figures/sorted_heatmap/fly_5pnt_pnk_box.png', dpi=150, bbox_inches='tight')
plt.show()
#%%
#Draw boxplot DC21-DC22 for terminal nucleotides

def boxplot_DC21_DC22(): 
    df_using = df_dcr_combine[(df_dcr_combine['Cleavage_site'].isin([21, 22]))]
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
    plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/dcr_5pnt_box_DC21_DC22.png', dpi=150, bbox_inches='tight')
    return
    plt.show()
    
boxplot_DC21_DC22()

def boxplot_DC21_DC22_3p_end(): 
    df_using = df_dcr_combine[(df_dcr_combine['Cleavage_site'].isin([21, 22]))]
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
    plt.savefig('D:/1.miRNA/END_randomized project/DCR1-DCL1/figures/dcr_3pnt_box_DC21_DC22.png', dpi=150, bbox_inches='tight')
    plt.show()
    return
boxplot_DC21_DC22_3p_end()
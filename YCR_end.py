#Analyzing the interplay between YCR and 5pend nucleotide impact
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
file_path = 'D:/1.miRNA/END_randomized project/YCR_end_analysis/'
#input YCR file as df 
df_ycr= pd.read_excel(file_path+'YCR_fromNAR.xlsx')
df_human = df_ycr[df_ycr['species']=='Hsa']
df_human = df_human.reset_index(drop=True) #data from mirgene db
df_human['premirID']=df_human['premirID'].str[1:]
df_folded = pd.read_csv('D:/1.miRNA/END_randomized project/mirna_check/prefolding_output.txt', sep ='\t')
df_folded['variant']=df_folded['variant'].str[:-4]
df_folded=df_folded[df_folded['variant'].str.contains('Hsa')]

df_human = pd.merge(df_human, df_folded, left_on='premirID', right_on='variant')
df_human['5p_nu']=df_human['preseq'].str[:1]
df_human['clv_assign']=df_human['cleavage_site'].astype(str)
df_human.loc[~df_human['clv_assign'].isin(['21','22']),'clv_assign']='other'
df_human['5p_DC'] = df_human['5p_nu']+df_human['clv_assign']
del df_human['variant']
del df_human['seq']
#Processing data for human pre-miRNA
df_noycr = df_human[(df_human['pos1']=='No')&(df_human['pos2']=='No') &(df_human['pos3']=='No')&(df_human['pos4']=='No')&(df_human['pos5']=='No')]
df_noycr =df_noycr.reset_index(drop=True)
df_haveycr = df_human[~((df_human['pos1']=='No')&(df_human['pos2']=='No') &(df_human['pos3']=='No')&(df_human['pos4']=='No')&(df_human['pos5']=='No'))]
df_haveycr =df_haveycr.reset_index(drop=True)


#Counting for all human pre-miRNAs
df_all_count = df_human.groupby('5p_DC')['premirID'].nunique()
df_all_count = df_all_count.to_frame(name='number').reset_index()
df_all_count['DC'] = df_all_count['5p_DC'].str[1:]
df_all_count['5p'] = df_all_count['5p_DC'].str[:1]



df_data =  df_all_count.copy()

total_pre = df_data.groupby('5p')['number'].sum()
df_data['Proportion']=df_data['number']/df_data['5p'].map(total_pre)

dc_list=['21', '22', 'other']
color = {'other': '#FAF1F1', '21':'#5AC9A1', '22':'#f94449'}
order=['A', 'G', 'U', 'C']
fig, ax = plt.subplots(figsize=(8, 10))
df_data['5p'] = pd.Categorical(df_data['5p'], categories=order, ordered=True)
df_data.groupby(['5p','DC'])['Proportion'].sum().unstack().plot(kind='bar', stacked=True, color=color, ax=ax, legend=False)
plt.xlabel('')
plt.ylabel('')
plt.xticks(rotation=0, fontsize=24)
plt.yticks(fontsize=14)
plt.tick_params(left = True, bottom = True, right = False , labelleft = False, labelbottom = False) 
plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/5p_OVH_hsa_DC_stack_YCR_data.png', dpi=150, bbox_inches='tight')
plt.show()

#counting for no YCR group
df_noycr['5p_DC']=df_noycr['5p_nu']+df_noycr['clv_assign']
df_noycr_count = df_noycr.groupby('5p_DC')['premirID'].nunique()
df_noycr_count = df_noycr_count.to_frame(name='number').reset_index()
df_noycr_count['DC'] = df_noycr_count['5p_DC'].str[1:]
df_noycr_count['5p'] = df_noycr_count['5p_DC'].str[:1]


df_data =  df_noycr_count.copy()

total_pre = df_data.groupby('5p')['number'].sum()
df_data['Proportion']=df_data['number']/df_data['5p'].map(total_pre)

dc_list=['21', '22', 'other']
color = {'other': '#FAF1F1', '21':'#5AC9A1', '22':'#f94449'}
order=['A', 'G', 'U', 'C']
fig, ax = plt.subplots(figsize=(6, 10))
df_data['5p'] = pd.Categorical(df_data['5p'], categories=order, ordered=True)
df_data.groupby(['5p','DC'])['Proportion'].sum().unstack().plot(kind='bar', stacked=True, color=color, ax=ax)
plt.xlabel('')
plt.ylabel('')
plt.xticks(rotation=0, fontsize=24)
plt.yticks(fontsize=14)
#plt.tick_params(left = True, bottom = True, right = False , labelleft = False, labelbottom = False) 
#plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/5p_OVH_hsa_DC_stack.png', dpi=150, bbox_inches='tight')
plt.show()

#Counting for have-YCR group
df_haveycr['5p_DC']=df_haveycr['5p_nu']+df_haveycr['clv_assign']
df_haveycr_count = df_haveycr.groupby('5p_DC')['premirID'].nunique()
df_haveycr_count = df_haveycr_count.to_frame(name='number').reset_index()
df_haveycr_count['DC'] = df_haveycr_count['5p_DC'].str[1:]
df_haveycr_count['5p'] = df_haveycr_count['5p_DC'].str[:1]


df_data =  df_haveycr_count.copy()

total_pre = df_data.groupby('5p')['number'].sum()
df_data['Proportion']=df_data['number']/df_data['5p'].map(total_pre)

dc_list=['21', '22', 'other']
color = {'other': '#FAF1F1', '21':'#5AC9A1', '22':'#f94449'}
order=['A', 'G', 'U', 'C']
fig, ax = plt.subplots(figsize=(6, 10))
df_data['5p'] = pd.Categorical(df_data['5p'], categories=order, ordered=True)
df_data.groupby(['5p','DC'])['Proportion'].sum().unstack().plot(kind='bar', stacked=True, color=color, ax=ax)
plt.xlabel('')
plt.ylabel('')
plt.xticks(rotation=0, fontsize=24)
plt.yticks(fontsize=14)
#plt.tick_params(left = True, bottom = True, right = False , labelleft = False, labelbottom = False) 
#plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/5p_OVH_hsa_DC_stack.png', dpi=150, bbox_inches='tight')
plt.show()


#%% Checking the group of G-DC22 and U-DC21

df_G_DC22 = df_human[(df_human['5p_nu']=='G')&(df_human['clv_assign']=='22')]
print(len(df_G_DC22))
#Now check the distribution based on OVH/YCR/stem-length
df_G_DC22.loc[:, 'overhang'] = df_G_DC22['concrete_struct'].apply(lambda x: x[:2] + x[-2:])

def comparing_ovh_G_DC22(): 
    value_counts = df_G_DC22['overhang'].value_counts().reset_index()
    value_counts.columns = ['overhang', 'count']

    # Plotting using seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x='overhang', y='count', data=value_counts)
    plt.title('Checking_OVH_in_G_DC22')
    plt.xlabel('Type of overhang')
    plt.ylabel('Number of pre-miRNAs')
    plt.xticks(rotation=45)
    plt.yticks([0,5,10,15,20])
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/G_DC22_ovh_check.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
comparing_ovh_G_DC22()

def comparing_YCR_G_DC22():
    df_using = df_G_DC22[['pos1','pos2', 'pos3', 'pos4', 'pos5']]
    yes_counts = df_using.eq('Yes').sum()

    # Plotting using seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x=yes_counts.index, y=yes_counts.values)
    plt.title('Position of YCR appearance')
    plt.xlabel('Positions')
    plt.ylabel('Number of pre-miRNAs')
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/G_DC22_ycr_check.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
comparing_YCR_G_DC22()

def counting_YCR_G_DC22(): 
    df_using = df_G_DC22[['pos1','pos2', 'pos3', 'pos4', 'pos5']]
    all_no_count = 0
    yes_one_pos_count = 0
    yes_multi_pos_count = 0

    for index, row in df_using.iterrows():
        num_yes = sum(value == 'Yes' for value in row)
    
        if num_yes == 0:
            all_no_count += 1
        elif num_yes == 1:
            yes_one_pos_count += 1
        else:
            yes_multi_pos_count += 1
    plot_data = pd.DataFrame({
            'Category': ['NO YCR', 'YCR at one position', 'YCR at Multi Positions'],
            'Count': [all_no_count, yes_one_pos_count, yes_multi_pos_count]})
    plt.figure(figsize=(8, 6))
    sns.barplot(x='Category', y='Count', data=plot_data)
    plt.xlabel('')
    plt.ylabel('Number of pre-miRNAs')
    plt.title('YCR checking in DC22')
    plt.yticks([0,5,10,15,20])
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/G_DC22_ycr_check_counting.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
counting_YCR_G_DC22()

#Checking U_DC21 group

df_U_DC21 = df_human[(df_human['5p_nu']=='U')&(df_human['clv_assign']=='21')]
print(len(df_U_DC21))
#Now check the distribution based on OVH/YCR/stem-length
df_U_DC21.loc[:, 'overhang'] = df_U_DC21['concrete_struct'].apply(lambda x: x[:2] + x[-2:])
def comparing_ovh_U_DC21(): 
    value_counts = df_U_DC21['overhang'].value_counts().reset_index()
    value_counts.columns = ['overhang', 'count']

    # Plotting using seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x='overhang', y='count', data=value_counts)
    plt.title('Checking_OVH_in_U_DC21')
    plt.xlabel('Type of overhang')
    plt.ylabel('Number of pre-miRNAs')
    plt.xticks(rotation=45)
    plt.yticks([0,5,10,15,20])
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/U_DC21_ovh_check.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
comparing_ovh_U_DC21()

def counting_YCR_U_DC21(): 
    df_using = df_U_DC21[['pos1','pos2', 'pos3', 'pos4', 'pos5']]
    all_no_count = 0
    yes_one_pos_count = 0
    yes_multi_pos_count = 0

    for index, row in df_using.iterrows():
        num_yes = sum(value == 'Yes' for value in row)
    
        if num_yes == 0:
            all_no_count += 1
        elif num_yes == 1:
            yes_one_pos_count += 1
        else:
            yes_multi_pos_count += 1
    plot_data = pd.DataFrame({
            'Category': ['NO YCR', 'YCR at one position', 'YCR at Multi Positions'],
            'Count': [all_no_count, yes_one_pos_count, yes_multi_pos_count]})
    plt.figure(figsize=(8, 6))
    sns.barplot(x='Category', y='Count', data=plot_data)
    plt.xlabel('')
    plt.ylabel('Number of pre-miRNAs')
    plt.title('YCR checking in U_DC21')
    plt.yticks([0,5,10,15,20,25,30,35])
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/U_DC21_ycr_check_counting.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
counting_YCR_U_DC21()

def comparing_YCR_U_DC21():
    df_using = df_U_DC21[['pos1','pos2', 'pos3', 'pos4', 'pos5']]
    yes_counts = df_using.eq('Yes').sum()

    # Plotting using seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x=yes_counts.index, y=yes_counts.values, color='green')
    plt.title('Position of YCR appearance')
    plt.xlabel('Positions')
    plt.ylabel('Number of pre-miRNAs')
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/U_DC21_ycr_check.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
comparing_YCR_U_DC21()



#Similar check for G_DC21 and U_DC22 - expected group: 
df_G_DC21 = df_human[(df_human['5p_nu']=='G')&(df_human['clv_assign']=='21')]
print(len(df_G_DC21))
#Now check the distribution based on OVH/YCR/stem-length
df_G_DC21.loc[:, 'overhang'] = df_G_DC21['concrete_struct'].apply(lambda x: x[:2] + x[-2:])

def comparing_ovh_G_DC21(): 
    value_counts = df_G_DC21['overhang'].value_counts().reset_index()
    value_counts.columns = ['overhang', 'count']

    # Plotting using seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x='overhang', y='count', data=value_counts)
    plt.title('Checking_OVH_in_G_DC21')
    plt.xlabel('Type of overhang')
    plt.ylabel('Number of pre-miRNAs')
    plt.xticks(rotation=45)
    plt.yticks([0,5,10,15,20])
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/G_DC21_ovh_check.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
comparing_ovh_G_DC21()

def counting_YCR_G_DC21(): 
    df_using = df_G_DC21[['pos1','pos2', 'pos3', 'pos4', 'pos5']]
    all_no_count = 0
    yes_one_pos_count = 0
    yes_multi_pos_count = 0

    for index, row in df_using.iterrows():
        num_yes = sum(value == 'Yes' for value in row)
    
        if num_yes == 0:
            all_no_count += 1
        elif num_yes == 1:
            yes_one_pos_count += 1
        else:
            yes_multi_pos_count += 1
    plot_data = pd.DataFrame({
            'Category': ['NO YCR', 'YCR at one position', 'YCR at Multi Positions'],
            'Count': [all_no_count, yes_one_pos_count, yes_multi_pos_count]})
    plt.figure(figsize=(8, 6))
    sns.barplot(x='Category', y='Count', data=plot_data)
    plt.xlabel('')
    plt.ylabel('Number of pre-miRNAs')
    plt.title('YCR checking in G_DC21')
    plt.yticks([0,5,10,15,20,25,30,35])
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/G_DC21_ycr_check_counting.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
counting_YCR_G_DC21()

def comparing_YCR_G_DC21():
    df_using = df_G_DC21[['pos1','pos2', 'pos3', 'pos4', 'pos5']]
    yes_counts = df_using.eq('Yes').sum()

    # Plotting using seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x=yes_counts.index, y=yes_counts.values, color='green')
    plt.title('Position of YCR appearance')
    plt.xlabel('Positions')
    plt.ylabel('Number of pre-miRNAs')
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/G_DC21_ycr_check.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
comparing_YCR_G_DC21()


df_U_DC22 = df_human[(df_human['5p_nu']=='U')&(df_human['clv_assign']=='22')]
print(len(df_U_DC22))
#Now check the distribution based on OVH/YCR/stem-length
df_U_DC22.loc[:, 'overhang'] = df_U_DC22['concrete_struct'].apply(lambda x: x[:2] + x[-2:])
def comparing_ovh_U_DC22(): 
    value_counts = df_U_DC22['overhang'].value_counts().reset_index()
    value_counts.columns = ['overhang', 'count']

    # Plotting using seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x='overhang', y='count', data=value_counts)
    plt.title('Checking_OVH_in_U_DC22')
    plt.xlabel('Type of overhang')
    plt.ylabel('Number of pre-miRNAs')
    plt.xticks(rotation=45)
    plt.yticks([0,5,10,15,20])
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/U_DC22_ovh_check.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
comparing_ovh_U_DC22()


def counting_YCR_U_DC22(): 
    df_using = df_U_DC22[['pos1','pos2', 'pos3', 'pos4', 'pos5']]
    all_no_count = 0
    yes_one_pos_count = 0
    yes_multi_pos_count = 0

    for index, row in df_using.iterrows():
        num_yes = sum(value == 'Yes' for value in row)
    
        if num_yes == 0:
            all_no_count += 1
        elif num_yes == 1:
            yes_one_pos_count += 1
        else:
            yes_multi_pos_count += 1
    plot_data = pd.DataFrame({
            'Category': ['NO YCR', 'YCR at one position', 'YCR at Multi Positions'],
            'Count': [all_no_count, yes_one_pos_count, yes_multi_pos_count]})
    plt.figure(figsize=(8, 6))
    sns.barplot(x='Category', y='Count', data=plot_data)
    plt.xlabel('')
    plt.ylabel('Number of pre-miRNAs')
    plt.title('YCR checking in U_DC22')
    plt.yticks([0,5,10,15,20,25,30,35])
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/U_DC22_ycr_check_counting.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
counting_YCR_U_DC22()


def comparing_YCR_U_DC22():
    df_using = df_U_DC22[['pos1','pos2', 'pos3', 'pos4', 'pos5']]
    yes_counts = df_using.eq('Yes').sum()

    # Plotting using seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x=yes_counts.index, y=yes_counts.values, color='green')
    plt.title('Position of YCR appearance')
    plt.xlabel('Positions')
    plt.ylabel('Number of pre-miRNAs')
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/U_DC22_ycr_check.png', dpi=150, bbox_inches='tight')
    plt.show()
    return()
comparing_YCR_U_DC22()
#Plot pie chart to visualize the proportion of YCR pre-miRNA 

# Data for the pie chart. This part drawing with matplotlib manually
import matplotlib.pyplot as plt


def piechart_analysis_pre():
    labels = 'noYCR', 'YCR_else', 'YCR_true'
    sizes = [32, 2, 8]
    explode = (0, 0, 0.1)
    fig, ax = plt.subplots(figsize=(8,8))
    colors = ['#d1d3d4', '#66b3ff', '#2bb673']
    ax.pie(sizes, explode=explode, startangle=90, colors=colors)
    #plt.show()
    plt.savefig('D:/1.miRNA/END_randomized project/mirna_check/figures/U_DC21_ycr_check_piechart.png', dpi=150, bbox_inches='tight')
    return()
piechart_analysis_pre()
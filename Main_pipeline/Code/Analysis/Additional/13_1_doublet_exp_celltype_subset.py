import pandas as pd
import os

os.chdir('/home/sidrah19220/nmd_review/reviewer2/GSVA')
exp=pd.read_csv('normalised_expression_wild.csv')
exp.index=exp['Unnamed: 0']
exp=exp.drop(['Unnamed: 0'],axis=1)
cms=exp.columns
imt=[x for x in cms if x.startswith('IMT')]
df5=exp[imt]
tdf5=df5.T
tdf5.to_csv('ttOSNs.csv')

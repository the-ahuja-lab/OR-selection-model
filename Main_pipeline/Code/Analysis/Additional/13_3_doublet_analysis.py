import pandas as pd
import numpy as np
import scrublet as scr
data = pd.io.parsers.read_csv("1raw_input_for_doublet.csv",sep=",",index_col=0)
indexNamesArr = data.index.values
listOfRowIndexLabels = list(indexNamesArr)
listOfRowIndexLabels
row_names=pd.DataFrame(listOfRowIndexLabels)
scrub = scr.Scrublet(data)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
score=pd.DataFrame(doublet_scores)
value=pd.DataFrame(predicted_doublets)
a=np.concatenate((row_names,score),axis=1)
b=np.concatenate((a,value),axis=1)
final=pd.DataFrame(b)
final.columns = ['Cell','Doublet_Score','Doublet']
print(final)
final.to_csv(r'tdoublet_analysis.csv',index = None, header=True)
exit()
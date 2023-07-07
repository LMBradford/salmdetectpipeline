import pandas as pd

#Read in the previously made summary tables
df_reads = pd.read_csv(snakemake.input[0], sep='\t')
df_hits = pd.read_csv(snakemake.input[1], sep='\t')
df_SSR = pd.read_csv(snakemake.input[2], sep='\t')

#Rename duplicated column names
df_SSR.rename(columns={'unique_reads':'Hits to SSRs','S_enterica_reads':'S.ent hits to SSRs','Non-Salm_hits':'SSRs False Pos'}, inplace=True)


#Merge dataframes together into one master dataframe
df_merged = df_reads.merge(df_hits, on='fileID').merge(df_SSR, on='fileID')

#Master dataframe cleanup and data analysis
#inplace = True means to overwrite the existing df instead of creating a new one
df_merged.drop(['filename_x', 'filename_y', 'Salm reads IDed as Salm'], axis=1, inplace=True) #Remove unnecessary columns

#Calculate true positive rate aka Recall (for PR curves)
#TPR = Tru pos / (True pos + False neg)
# True pos + False neg = total number of Salmonella reads in library
df_merged["Recall_Kr2"] = df_merged["S.ent reads correctly IDed"]/ df_merged["Senterica reads in library"]
df_merged["Recall_SSRs"] = df_merged["S.ent hits to SSRs"]/ df_merged["Senterica reads in library"]



#Calculate precision for PR curves
#Precision = True pos / (True pos + False pos)
df_merged["PrKr2denom"] = df_merged["S.ent reads correctly IDed"] + df_merged["False positives"]
df_merged["Precision Kr2"] = df_merged["S.ent reads correctly IDed"] / df_merged["PrKr2denom"]

df_merged["PrSRRdenom"] = df_merged["S.ent hits to SSRs"] + df_merged["SSRs False Pos"]
df_merged["Precision SSRs"] = df_merged["S.ent hits to SSRs"] / df_merged["PrSRRdenom"]

df_merged.drop(["PrKr2denom", "PrSRRdenom"], axis=1, inplace=True) #Remove denominator columns, no longer needed
df_merged.drop(['filename'], axis=1, inplace=True)

#Some cells have div/0 errors (for SSRs), cells end up blank.
#Replace blank cells with 0 so they plot properly.
df_merged.fillna(0, inplace=True)
#pd.set_option('display.max_columns', 20)
#print(df_merged)

df_merged.to_csv(snakemake.output[0], sep = "\t", index = False)

import pandas as pd #for dataframes

# Set up lists
# Lists are populated with info read from file
# Lists are then pulled into pandas database
FILENAMES, UNIQLINES = [], []
FILEIDS = []
SALMHITCOUNTS = []
sent_hits = []

# Loop through input files defined by snakemake rules
for f in snakemake.input:
  with open(f, 'r') as fileObj:

    #Save just filename, not whole path
    filename = f.split("/")[3]

    #Separate details from filename
    comm = filename.split("_")[0]  
    db = filename.split("_")[1]
    conf = filename.split("_")[2]        
    
    #Combine details into unique fileID  
    fileid = comm +"_"+ db +"_"+ conf    

    #Set up counters
    tlines=0
    sent = 0
    
    #Go through lines to find genome reads came from
    for line in fileObj:
      tlines +=1
      
      if "Salmonella_enterica" in line:
        sent +=1


    #Fill in lists with counts
    FILENAMES.append(filename)
    UNIQLINES.append(tlines)
    SALMHITCOUNTS.append(sent)
    FILEIDS.append(fileid)
    
#Create dataframe with info from lists
df_blhits = pd.DataFrame({
  'filename': FILENAMES,
  'fileID': FILEIDS,
  'unique_reads': UNIQLINES,
  'S_enterica_reads': SALMHITCOUNTS,
  })
 
df_blhits['Non-Salm_hits'] = df_blhits['unique_reads'] - df_blhits['S_enterica_reads']

#Save dataframe as tab-delimited file, with name defined by snakemake rule
df_blhits.to_csv(snakemake.output[0], sep = "\t", index = False)

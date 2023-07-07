import pandas as pd

# Set up lists
# Lists are populated with info read from file
# Lists are then pulled into pandas database
FILENAMES, TOTALLINES, SALM_LINES = [], [], []
FILEIDS = []

for f in snakemake.input:	
  with open(f, 'r') as fileObj:
    namesplit = f.split("/")
    
    db = namesplit[1]
    filename = namesplit[3]

    comm = filename.split("_")[0]
    conf = filename.split("_")[2]

    #Combine details into unique fileID  
    fileid = comm +"_"+ db +"_"+ conf

    #Add file info to external lists
    FILENAMES.append(filename) #Add file name to names list
    FILEIDS.append(fileid) #Add unique fileID to list

#Set up lines to be counted in following loops
    totlines=0
    salmlines_hit=0

    for line in fileObj: 	
      totlines +=1 #Find total number of lines

      col2 = line.split('\t')[1] #Check only column 2 (read IDs) for names

			#Count instances of Salmonella enterica specifically
      if "Salmonella_enterica" in col2: 
        salmlines_hit +=1

		#Add outputs of loops to lists		
    TOTALLINES.append(totlines)
    SALM_LINES.append(salmlines_hit)

    
df_hits = pd.DataFrame({
  'filename': FILENAMES,
  'fileID': FILEIDS, 
  'Tot reads IDed as Salm': TOTALLINES, 
  'Salm reads IDed as Salm': SALM_LINES
  })
 

df_hits["False positives"] = df_hits["Tot reads IDed as Salm"] - df_hits["Salm reads IDed as Salm"]

#print(df_hits)

df_hits.to_csv(snakemake.output[0], sep = "\t", index = False)

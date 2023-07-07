import pandas as pd

# Set up lists
# Lists are populated with info read from file
# Lists are then pulled into pandas database
FILENAMES, TOTALLINES, SALMCOUNTS, TAXCOUNTS, UNCLASSCOUNTS = [], [], [], [], []
ANROUNDS, DATABASES, COMMUNITIES, CONFS, FILEIDS = [], [], [], [], []

#Full taxonomic pathway for the Salmonella genus
#Reads IDed as being in this path are not false IDs, just not very specific
taxpath = ["Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae"]

for f in snakemake.input:	
  with open(f, 'r') as fileObj:
    namesplit = f.split("/")
    
    rnd = namesplit[0]
    db = namesplit[1]
    filename = namesplit[3]

    comm = filename.split("_")[0]
    conf = filename.split("_")[2]

    #Combine details into unique fileID  
    fileid = comm +"_"+ db +"_"+ conf

    #Add file info to external lists
    FILENAMES.append(filename) #Add file name to names list
    ANROUNDS.append(rnd)
    DATABASES.append(db)
    COMMUNITIES.append(comm)
    CONFS.append(conf)
    FILEIDS.append(fileid)
 
    #Set up variables to store counts
    totlines=0
    salmlines_read=0
    taxlines_read=0
    unclass_read=0

    for line in fileObj: 	
      totlines +=1 #Find total number of lines

      col3 = line.split('\t')[2] #Check only column 3 (read assignments) for names

			#Count instances of Salmonella enterica specifically
      if "Salmonella" in col3: 
        salmlines_read +=1
			#Count instances of anything in Salmonella's taxnomic path
      if any(tax in col3 for tax in taxpath): #checks for all strings in taxpath list
        taxlines_read +=1
      #Count instances of read being unclassified
      if "unclassified (taxid 0)" in col3:
        unclass_read +=1
      
		#Add outputs of counting loops to lists		
    TOTALLINES.append(totlines)
    SALMCOUNTS.append(salmlines_read)
    TAXCOUNTS.append(taxlines_read)
    UNCLASSCOUNTS.append(unclass_read)

df_reads = pd.DataFrame({
  'filename': FILENAMES,
  'fileID': FILEIDS,
  'Analysis round': ANROUNDS,
  'Database': DATABASES,
  'Kr2 confidence': CONFS,  
  'Senterica reads in library': TOTALLINES, 
  'S.ent reads correctly IDed': SALMCOUNTS,
  'S.ent reads assigned higher': TAXCOUNTS,
  'Unclassified': UNCLASSCOUNTS})
 
df_reads["S.ent reads incorrectly assigned"] = df_reads["Senterica reads in library"] - df_reads["S.ent reads correctly IDed"] - df_reads["S.ent reads assigned higher"] - df_reads["Unclassified"]

#print(df_reads)

df_reads.to_csv(snakemake.output[0], sep = "\t", index = False)

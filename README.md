# Genebank parser for R   
Genebank files (gb) are challenging to parse in R. 
Using this script you will be able to load them into your environment.   
      
## How do I use it?      
Download the script genebankreader.R into your system and load it using the function source. (e.g., source(path_to_donwloaded_file/genebankreader.R).
Then you call the genebankreader function (e.g. test<-genebankreader("/media/nacho/Data/RDKW20.gb")).   
The results will be a dataframe with one row per CDS and in the columns you will find the name of the CDS, where it is located and the info regardng the gene present in the gb file.


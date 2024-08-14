genebankreader<-function(input){
  gbref<-readLines(input)
  
  out.cds<-list()
  out.gene<-list()
  geneon<-FALSE
  cdson<-FALSE
  pb<-txtProgressBar(max = length(gbref))
  for (i in 1:length(gbref)) {
    setTxtProgressBar(pb,i)
    
    if(length(grep("^     gene", gbref[i]))==1 ){
      geneon<-TRUE
      cdson<-FALSE
      gene.temp<-as.data.frame(gsub(".* ", "", gbref[i]))
      colnames(gene.temp)<-"Position"
      if(exists("cds.temp")) out.cds<-c(out.cds, list(cds.temp))
    }
    
    if(length(grep("^     CDS", gbref[i]))==1 ){
      geneon<-FALSE
      cdson<-TRUE
      cds.temp<-as.data.frame(gsub(".* ", "", gbref[i]))
      colnames(cds.temp)<-"Position"
      if(exists("gene.temp")) out.gene<-c(out.gene, list(gene.temp))
    }
    if(length(grep("translation=", gbref[i]))==1 ) cdson<-FALSE
    
    if(geneon){
      gene.temp$dummy<-gsub("\"" , "", gsub(".*=", "", gbref[i]))
      colnames(gene.temp)[which(colnames(gene.temp)=="dummy")]<- gsub("/","",gsub(" ","",gsub("\"" , "", gsub("=.*", "", gbref[i]))))
    }
    
    if(cdson){
      cds.temp$dummy<-gsub("\"" , "", gsub(".*=", "", gbref[i]))
      colnames(cds.temp)[which(colnames(cds.temp)=="dummy")]<- gsub("/","",gsub(" ","",gsub("\"" , "", gsub("=.*", "", gbref[i]))))
    }
    
  }
  
  cds.factors<-c("codon_start", "locus_tag","note","Position", "product","protein_id","transl_table") 
  
  cds.clean<-list()
  for (i in 1:length(out.cds)) {
    temp<-out.cds[[i]]
    temp<-temp[,-2]
    to.fix<-  c(1:ncol(temp))[-which(colnames(temp)%in% cds.factors)]
    if(length(to.fix)>0){
      to.fix<-to.fix[order(to.fix,decreasing = TRUE)]
      for (i in 1:length(to.fix)) {
        temp[,to.fix[i]-1]<-paste(temp[,to.fix[i]-1], gsub(" ","",temp[,to.fix[i]]) )
      }
      temp<-temp[,-to.fix]
    }
    cds.clean<-c(cds.clean,list(temp))
  }
  
  cds.clean<-do.call(rbind,cds.clean)
  reverse<-grep("complement", cds.clean$Position)
  
  cds.clean$Start<-NA
  cds.clean$End<-NA
  cds.clean$End[c(1:nrow(cds.clean))[-reverse]] <- as.numeric(gsub(".*\\.\\.","",cds.clean$Position[c(1:nrow(cds.clean))[-reverse]] ))
  cds.clean$Start[c(1:nrow(cds.clean))[-reverse]] <- as.numeric(gsub("\\.\\..*","",cds.clean$Position[c(1:nrow(cds.clean))[-reverse]] ))
  
  cds.clean$Position<-gsub(".*\\(", "", gsub("\\)", "", cds.clean$Position))
  cds.clean$Start[c(1:nrow(cds.clean))[reverse]] <- as.numeric(gsub(".*\\.\\.","",cds.clean$Position[c(1:nrow(cds.clean))[reverse]] ))
  cds.clean$End[c(1:nrow(cds.clean))[reverse]] <- as.numeric(gsub("\\.\\..*","",cds.clean$Position[c(1:nrow(cds.clean))[reverse]] ))
  cds.clean$Orientation<-"FW"
  cds.clean$Orientation[reverse]<-"RV"
  return(cds.clean)
}
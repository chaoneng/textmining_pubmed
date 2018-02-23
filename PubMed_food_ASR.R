##pubmed search and download
foodname=read.csv("~/Google 雲端硬碟/Research/FoodinDisease/Code/food_sci_name.csv")
#dim(foodname)
#head(foodname)
#lowfood=tolower(as.matrix(foodname))

lowfood=as.data.frame(foodname)
#dim(lowfood)
#lowfood$name

library(easyPubMed)
library(kulife)
my.query <- paste(lowfood$name,"[MeSH Terms]","OR",lowfood$name,"[All Fields]) AND (2008/01/27[PDat] : 2018/01/23[PDat])")
my.idlist <- get_pubmed_ids(my.query[1])
#my.idlist <- lapply(my.query,get_pubmed_ids)

count<-as.numeric(my.idlist$Count)
#my.seq <-seq(1, as.numeric(my.idlist$Count), by = count-1)

batch.xml <- fetch_pubmed_data(my.idlist, retstart = 1, retmax = 2000)
saveXML(batch.xml,"~/Google 雲端硬碟/Research/FoodinDisease/Code/query3.xml")

##createxml function
#createxml <- function(x) {
#  saveXML(x,paste("~/Google 雲端硬碟/Research/FoodinDisease/Code/",x,".xml",sep = ""))
#}

#batch.xml<-lapply(my.idlist,fetch_pubmed_data)
######

##use xml.r to split
library(tm)
library(XML)

extract_xml <- function(theFile) {
  #library(XML)
  newData <- xmlParse(theFile)
  
  records <- getNodeSet(newData, "//PubmedArticle")
  pmid <- xpathSApply(newData,"//MedlineCitation/PMID", xmlValue)
  authLast <- lapply(records, xpathSApply, ".//Author/LastName", xmlValue)
  
  ## affiliations <- lapply(records, xpathSApply, ".//Author/AffiliationInfo/Affiliation", xmlValue)
  ## affiliations[sapply(affiliations, is.list)] <- NA
  ## affiliations <- sapply(affiliations, paste, collapse = "|")
  
  year <- lapply(records, xpathSApply, ".//PubDate/Year", xmlValue) 
  year[sapply(year, is.list)] <- NA
  year <- unlist(year)
  
  articletitle <- lapply(records, xpathSApply, ".//ArticleTitle", xmlValue) 
  articletitle[sapply(articletitle, is.list)] <- NA
  articletitle <- unlist(articletitle)
  
  journal <- lapply(records, xpathSApply, ".//ISOAbbreviation", xmlValue) 
  journal[sapply(journal, is.list)] <- NA
  journal <- unlist(journal)
  
  abstract <- lapply(records, xpathSApply, ".//Abstract/AbstractText", xmlValue)
  abstract[sapply(abstract, is.list)] <- NA
  abstract <- sapply(abstract, paste, collapse = "|")
  
  theDF <- data.frame(pmid, year, articletitle, journal, abstract,stringsAsFactors = FALSE)
  return(theDF)
}

Angelica = extract_xml("~/Google 雲端硬碟/Research/FoodinDisease/Code/query3.xml")

pmid=Angelica$pmid
ab=Angelica$abstract

#tm , _ % 
docs <- Corpus( VectorSource(ab) ) 

####

##trans word to lower

##clean_corpus function
clean_corpus <- function(corpus){
  # Eliminate extra white spaces
  corpus <- tm_map(corpus, stripWhitespace)
  # Remove punctuations
  corpus <- tm_map(corpus, removePunctuation)
  # Convert the text to lower case
  corpus <- tm_map(corpus, content_transformer(tolower))
  # Remove english common stopwords
  return(corpus)
}

docs<-clean_corpus(docs)

df <- data.frame(text = get("content", docs))

df2 = as.matrix(df)
####

###compare food name
##compa function
compa<-function(x){
  name<-ifelse(grepl(x, df2),x,NA)
  return(name)
}

library(stringr)
t1compa=lapply(str_to_lower(lowfood[,1]),compa)
t2compa=lapply(str_to_lower(lowfood[,2]),compa)

#####

####trans food data frame to matrix
library(data.table)
library(stringi)

###listcom&listcom2 function
listcom<-function(x){
  temp2 <- stri_list2matrix(x, byrow = TRUE)
  temp2[is.na(temp2)] <- ""
  gene=apply(temp2[,1:ncol(temp2)], 1 , paste , collapse = " " )
  return(gene)
}
listcom2<-function(x){
  temp2 <- stri_list2matrix(x, byrow = TRUE)
  temp2[is.na(temp2)] <- ""
  gene=apply(temp2[,1:ncol(temp2)], 1 , paste , collapse = "\t" )
  return(gene)
}

da=data.table(t(sapply(t1compa,c)))
food1=listcom(da)
removspace1<-strsplit(food1, " +")
transfood1=stri_list2matrix(removspace1, byrow = TRUE)

da2=data.table(t(sapply(t2compa,c)))
food2=listcom(da2)
removspace2<-strsplit(food2, " +")
transfood2=stri_list2matrix(removspace2, byrow = TRUE)

transfood1=as.data.frame(transfood1[,-1])
transfood2=as.data.frame(transfood2[,-1])

colnames(transfood1) <- paste('com1food',1:ncol(transfood1))
colnames(transfood2) <- paste('com2food',1:ncol(transfood2))

comfood1=cbind(pmid, as.data.frame(transfood1))
comfood2=cbind(comfood1, as.data.frame(transfood2))

class(comfood2)
dim(comfood2)
head(comfood2)
#run gene

library(pubmed.mineR)
pmids = comfood2$pmid

pubmids1=pmids[1:length(pmids)]
pubtator_result1=lapply(pubmids1,pubtator_function)

#pubmids2=pmids[1001:length(pmids)]
#pubtator_result2=lapply(pubmids2,pubtator_function)

library(stringi)
#library(gdata)
#allpub_result=c(pubtator_result1,pubtator_result2)

allpub_result=pubtator_result1

#library(rlist)
#x = list(pubtator_result1,pubtator_result2)
#allpub_result= list.cbind(x)

genec=sapply(allpub_result , "[",1)
diseasec=sapply(allpub_result , "[",2)
muationc=sapply(allpub_result , "[",3)
chemicalc=sapply(allpub_result , "[",4)
Speciesc=sapply(allpub_result , "[",5)

gene<-listcom2(genec)
disease<-listcom2(diseasec)
chemical<-listcom2(chemicalc)

fieldgene = strsplit(gene, '\t')
strgene <- stri_list2matrix(fieldgene, byrow = TRUE)

fielddisease = strsplit(disease, '\t')
strdisease <- stri_list2matrix(fielddisease,byrow = TRUE)

fieldchemical = strsplit(chemical, '\t')
strchemical <- stri_list2matrix(fieldchemical, byrow = TRUE)

strgene=as.data.frame(strgene)
strdisease=as.data.frame(strdisease)
strchemical=as.data.frame(strchemical)

colnames(strgene) <- paste('gene',1:ncol(strgene))
colnames(strdisease) <- paste('disease',1:ncol(strdisease))
colnames(strchemical) <- paste('chemical',1:ncol(strchemical))


allcombie=cbind(comfood2,strgene,strdisease,strchemical)

#write.csv(allcombie,"D:/foodtextmining.csv")

#foodtest=read.csv("D:/foodtextmining.csv")

allcombie$`gene 1`=sub("No Data", "",allcombie$`gene 1`)

df2=apply(allcombie, 2, function(x) gsub("^$|^ $|^  $", NA, x))

df3=apply(df2, 1, function(x) paste(na.omit(x),collapse=", ") )

df4=strsplit(df3, ",")

fileConn<-file("~/Google 雲端硬碟/Research/FoodinDisease/Code/food_gene_name2.csv")
writeLines(unlist(lapply(df4, paste, collapse=",")),fileConn)
close(fileConn)

####cbind food data query



######

####data mining association rule
library(xlsx)
library(arules)

##
text_demo  <- read.transactions("~/Google 雲端硬碟/Research/FoodinDisease/Code/food_gene_name.csv",format = "basket", sep=",")
head(text_demo)
summary(text_demo)

##################
itemFrequencyPlot(text_demo,topN=20,type="absolute")

rules <- apriori(text_demo, parameter = list(supp = 0.05, conf = 0.8))

options(digits=5)
inspect(rules[1:10])

rules_100=rules[1:300]
###############
library(RColorBrewer)
library(wordcloud)
product_name <- itemLabels(text_demo)
product_cnt <- itemFrequency(text_demo)*1000

#col.pal <- brewer.pal(9, "Blues")
col.pal <- brewer.pal(9, "Dark2")
wordcloud(words = product_name, freq = product_cnt, min.freq = 2, scale = c(3, 0.2), col = col.pal , random.order = FALSE)

##########
library(ggplot2)
library(arulesViz)
plot(rules_100)
plot(rules_100, method="graph", control=list(type="items"))
plot(rules_100, method="grouped")
plot(rules_100, method="paracoord", control=list(reorder=TRUE))

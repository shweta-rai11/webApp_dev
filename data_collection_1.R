library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(annotate)
library(preprocessCore)
library(EnsDb.Hsapiens.v79)

library(DMwR2)

#### Load data ####
###################
# gse38351 <- getGEO("GSE38351") testing
gse93272 <- getGEO("GSE93272")
gse93776 <- getGEO("GSE93776")
#gse68689 <- getGEO("GSE68689")
gse65010 <- getGEO("GSE65010")
gse55457 <- getGEO("GSE55457")
gse77298 <- getGEO("GSE77298")
gse12021 <- getGEO("GSE12021")
#gse1919 <- getGEO("GSE1919")


gse93272 <- gse93272[[1]]
gse93776 <- gse93776[[1]]
#gse68689 <- gse68689[[1]]
gse65010 <- gse65010[[1]]
gse55457 <- gse55457[[1]]
gse77298 <- gse77298[[1]]
gse12021 <- gse12021[[1]]
#gse1919 <- gse1919[[1]]



phen93272 <- pData(gse93272)
phen93776 <- pData(gse93776)
#phen68689 <- pData(gse68689)
phen65010 <- pData(gse65010)
phen55457 <- pData(gse55457)
phen77298 <- pData(gse77298)
phen12021 <- pData(gse12021)
#phen1919 <- pData(gse1919)




feat93272 <- fData(gse93272)
feat93776 <- fData(gse93776)
#feat68689 <- fData(gse68689)
feat65010 <- fData(gse65010)
feat55457 <- fData(gse55457)
feat77298 <- fData(gse77298)
feat12021 <- fData(gse12021)
#feat1919 <- fData(gse1919)


exp93272 <- exprs(gse93272)
exp93776 <- exprs(gse93776)
#exp68689 <- exprs(gse68689)
exp65010 <- exprs(gse65010)
exp55457 <- exprs(gse55457)
exp77298 <- exprs(gse77298)
exp12021 <- exprs(gse12021)
#exp1919 <- exprs(gse1919)


# GSE93272 preprocessing
exp93272 <- as.data.frame(exp93272)
exp93272 <- log2(exp93272)
exp93272 <- cbind(ID=rownames(exp93272), exp93272)
exp93272 <- merge(feat93272[,c(1,2,11)], exp93272, by="ID")
names(exp93272)[names(exp93272) == "Gene Symbol"] <- "Gene_Symbol"

# GSE93776 preprocessing
exp93776 <- as.data.frame(exp93776)
exp93776 <- log2(exp93776)
exp93776 <- cbind(ID=rownames(exp93776), exp93776)
exp93776 <- merge(feat93776[,c(1,2,11)], exp93776, by="ID")
names(exp93776)[names(exp93776) == "Gene Symbol"] <- "Gene_Symbol"
#View(feat1919)
# GSE65010 preprocessing
exp65010 <- as.data.frame(exp65010)
exp65010 <- log2(exp65010)
exp65010 <- cbind(ID=rownames(exp65010), exp65010)
exp65010 <- merge(feat65010[,c(1,2,11)], exp65010, by="ID")
names(exp65010)[names(exp65010) == "Gene Symbol"] <- "Gene_Symbol"

# GSE55457 preprocessing
exp55457 <- as.data.frame(exp55457)
exp55457 <- log2(exp55457)
exp55457 <- cbind(ID=rownames(exp55457), exp55457)
exp55457 <- merge(feat55457[,c(1,2,11)], exp55457, by="ID")
names(exp55457)[names(exp55457) == "Gene Symbol"] <- "Gene_Symbol"

# GSE77298 preprocessing
'''exp77298 <- as.data.frame(exp77298)
exp77298 <- log2(exp77298)
exp77298 <- cbind(ID=rownames(exp77298), exp77298)
exp77298 <- merge(feat77298[,c(1,2,11)], exp77298, by="ID")
names(exp77298)[names(exp77298) == "Gene Symbol"] <- "Gene_Symbol"'''

# GSE12021 preprocessing
exp12021 <- as.data.frame(exp12021)
exp12021 <- log2(exp12021)
exp12021 <- cbind(ID=rownames(exp12021), exp12021)
exp12021 <- merge(feat12021[,c(1,2,11)], exp12021, by="ID")
names(exp12021)[names(exp12021) == "Gene Symbol"] <- "Gene_Symbol"

#GSE1919 PREPROCESSING
'''exp1919 <- as.data.frame(exp1919)
exp1919 <- log2(exp1919)
exp1919 <- cbind(ID=rownames(exp1919), exp1919)
exp1919 <- merge(feat1919[,c(1,2,11)], exp1919, by="ID")
names(exp1919)[names(exp1919) == "Gene Symbol"] <- "Gene_Symbol"'''






dim(exp93272) # 54613   278
dim(exp93776) # 43617   176
dim(exp65010) # 54675    51
dim(exp55457) # 22283    36
#dim(exp77298) # 54675    26
dim(exp12021) # 22283    34
#dim(exp1919) # 12626    18

sum(278+176+51+36 +34 +18)
# total samples are 593 samples 
#View(commonGenes)


exp93272$Gene_Symbol <- sapply(strsplit(exp93272$Gene_Symbol, " /// "), `[`, 1)
exp93776$Gene_Symbol <- sapply(strsplit(exp93776$Gene_Symbol, " /// "), `[`, 1)
exp65010$Gene_Symbol <- sapply(strsplit(exp65010$Gene_Symbol, " /// "), `[`, 1)
exp55457$Gene_Symbol <- sapply(strsplit(exp55457$Gene_Symbol, " /// "), `[`, 1)
#exp77298$Gene_Symbol <- sapply(strsplit(exp77298$Gene_Symbol, " /// "), `[`, 1)
exp12021$Gene_Symbol <- sapply(strsplit(exp12021$Gene_Symbol, " /// "), `[`, 1)
#exp1919$Gene_Symbol <- sapply(strsplit(exp1919$Gene_Symbol, " /// "), `[`, 1)


feat93272 <- fData(gse93272)
feat93776 <- fData(gse93776)
feat65010 <- fData(gse65010)
feat55457 <- fData(gse55457)
#feat77298 <- fData(gse77298)
feat12021 <- fData(gse12021)
#feat1919 <- fData(gse1919)


#exp38351 <- cbind(ID=rownames(exp38351), exp38351) %>% as.data.frame()
#exp93272 <- cbind(ID=rownames(exp93272), exp93272) %>% as.data.frame()
#exp93776 <- cbind(ID=rownames(exp93776), exp93776) %>% as.data.frame()

#exp38351 <- merge(feat38351, exp38351, by = "Gene_Symbol")
#exp93272 <- merge(feat93272, exp93272, by = "Gene_Symbol")
#exp93776 <- merge(feat93776, exp93776, by = "Gene_Symbol")


exp93272 <- exp93272[which(exp93272$Gene_Symbol != ""),]
exp93776 <- exp93776[which(exp93776$Gene_Symbol != ""),]
exp65010 <- exp65010[which(exp65010$Gene_Symbol != ""),]
exp55457 <- exp55457[which(exp55457$Gene_Symbol != ""),]
#exp77298 <- exp77298[which(exp77298$Gene_Symbol != ""),]
exp12021 <- exp12021[which(exp12021$Gene_Symbol != ""),]
#exp1919 <- exp1919[which(exp1919$Gene_Symbol != ""),]



exp93272 <- exp93272[!duplicated(exp93272$Gene_Symbol),]
exp93776 <- exp93776[!duplicated(exp93776$Gene_Symbol),]
exp65010 <- exp65010[!duplicated(exp65010$Gene_Symbol),]
exp55457 <- exp55457[!duplicated(exp55457$Gene_Symbol),]
#exp77298 <- exp77298[!duplicated(exp77298$Gene_Symbol),]
exp12021 <- exp12021[!duplicated(exp12021$Gene_Symbol),]
#exp1919 <- exp1919[!duplicated(exp1919$Gene_Symbol),]



rownames(exp93272) <- exp93272$Gene_Symbol
rownames(exp93776) <- exp93776$Gene_Symbol
rownames(exp65010) <- exp65010$Gene_Symbol
rownames(exp55457) <- exp55457$Gene_Symbol
#rownames(exp77298) <- exp77298$Gene_Symbol
rownames(exp12021) <- exp12021$Gene_Symbol
#rownames(exp1919) <- exp1919$Gene_Symbol



exp93272 <- exp93272[,-c(1,2)]
exp93776 <- exp93776[, -c(1,2)]
exp65010 <- exp65010[,-c(1,2)]
exp55457 <- exp55457[,-c(1,2)]
#exp77298 <- exp77298[, -c(1,2)]
exp12021 <- exp12021[, -c(1,2)]
#exp1919 <- exp1919[, -c(1,2)]



exp93272 <- na.omit(exp93272)
exp93776 <- na.omit(exp93776)
exp65010 <- na.omit(exp65010)
exp55457 <- na.omit(exp55457)
#exp77298 <- na.omit(exp77298)
exp12021 <- na.omit(exp12021)
#exp1919 <- na.omit(exp1919)


write.csv(exp93272, "~/GSE_merging_data/exp93272.csv")
write.csv(exp93776, "~/GSE_merging_data/exp93776.csv")
write.csv(exp65010, "~/GSE_merging_data/exp65010.csv")
write.csv(exp55457, "~/GSE_merging_data/exp55457.csv")
#write.csv(exp77298, "~/blood_paper_work/exp_data/exp77298.csv")
write.csv(exp12021, "~/GSE_merging_data/exp12021.csv")
#write.csv(exp1919, "~/GSE_merging_data/exp1919.csv")

######################

#### Identifying biological classes ####
########################################



for(x in 1:nrow(phen93272)) {
  status <- phen93272$`characteristics_ch1.1`[x]
  if(status == "disease state: healthy control") phen93272$`characteristics_ch1.1`[x] <- "Control" else phen93272$`characteristics_ch1.1`[x] <- "RA"
}
for(x in 1:nrow(phen93272)) {
  status <- phen93272$`characteristics_ch1.6`[x]
  if(status == "gender: F") phen93272$`characteristics_ch1.6`[x] <- "Female" else phen93272$`characteristics_ch1.6`[x] <- "Male"
}

for(x in 1:nrow(phen93272)) {
  status <- phen93272$`characteristics_ch1`[x]
  if(status == "tissue: whole blood") phen93272$`characteristics_ch1.2`[x] <- "Blood" else phen93272$`characteristics_ch1.2`[x] <- "Blood"
}

for(x in 1:nrow(phen93776)) {
  status <- phen93776$`characteristics_ch1`[x]
  if(status == "disease: healthy control") phen93776$`characteristics_ch1`[x] <- "Control" else phen93776$`characteristics_ch1`[x] <- "RA"
}

for(x in 1:nrow(phen93776)) {
  status <- phen93776$`characteristics_ch1.1`[x]
  if(status == "gender: F") phen93776$`characteristics_ch1.1`[x] <- "Female" else phen93776$`characteristics_ch1.1`[x] <- "Male"
}




gender_data1 <- c("Female", "Female", "Female", "Female",
                 "Female", "Female", "Female", "Female",
                 "Female", "Female", "Female", "Female",
                 "Female", "Female", "Female", "Female",
                 "Male","Male","Male","Male",
                 "Male","Male","Male","Male",
                 "Female", "Female", "Female", "Female",
                 "Female", "Female", "Female", "Female",
                 "Female", "Female", "Female", "Female",
                 "Female", "Female", "Female", "Female",
                 "Male","Male","Male","Male",
                 "Male","Male","Male","Male") 
sampleInfo_10 <- phen65010 %>%
  mutate(gender = gender_data1)
sampleInfo_10$gender <- gender_data1
head(sampleInfo_10)

for(x in 1:nrow(sampleInfo_10)) {
  status <- sampleInfo_10$`characteristics_ch1.2`[x]
  if(status == "disease state: healthy control") sampleInfo_10$`characteristics_ch1.2`[x] <- "Control" else sampleInfo_10$`characteristics_ch1.2`[x] <- "RA"
}




for (x in 1:nrow(phen55457)) {
  status <- phen55457$`characteristics_ch1.2`[x]
  
  if (status == "clinical status: normal control") {
    phen55457$`characteristics_ch1.2`[x] <- "Control"
  } else if (status == "clinical status: osteoarthritis") {
    phen55457$`characteristics_ch1.2`[x] <- "OA"
  } else {
    phen55457$`characteristics_ch1.2`[x] <- "RA"
  }
}

for(x in 1:nrow(phen55457)) {
  status <- phen55457$`characteristics_ch1`[x]
  if(status == "gender: female") phen55457$`characteristics_ch1`[x] <- "Female" else phen55457$`characteristics_ch1`[x] <- "Male"
}

for (x in 1:nrow(phen12021)) {
  status <- phen12021$`characteristics_ch1.2`[x]
  
  if (status == "normal control") {
    phen12021$`characteristics_ch1.2`[x] <- "Control"
  } else if (status == "disease: osteoarthritis") {
    phen12021$`characteristics_ch1.2`[x] <- "OA"
  } else {
    phen12021$`characteristics_ch1.2`[x] <- "RA"
  }
}


for(x in 1:nrow(phen12021)) {
  status <- phen12021$`characteristics_ch1`[x]
  if(status == "sex: m") phen12021$`characteristics_ch1`[x] <- "Male" else phen12021$`characteristics_ch1`[x] <- "Female"
}

'''for (x in 1:nrow(phen1919)) {
  status <- phen1919$`source_name_ch1`[x]
  
  if (status == "normal donor") {
    phen1919$`source_name_ch1`[x] <- "Control"
  } else if (status == "osteoarthritis") {
    phen1919$`source_name_ch1`[x] <- "OA"
  } else {
    phen1919$`source_name_ch1`[x] <- "RA"
  }
}

colnames(phen1919)[colnames(phen1919) == "gender:ch1"] <- "Gender"'''


View(phen1919)

dim(phen93272_)
dim(phen93776_)
dim(phen65010_)
dim(phen55457_)
dim(phen12021_)
dim(phen1919_)


phen93272_ <- as.data.frame(cbind(sample = rownames(phen93272), group = phen93272[, c(11, 16)]))
colnames(phen93272_)[colnames(phen93272_) == "group.characteristics_ch1.1"] <- "Disease_state"
colnames(phen93272_)[colnames(phen93272_) == "group.characteristics_ch1.6"] <- "Gender"
phen93272_$platform <- "GPL570"
phen93272_$tissue <- "Blood"
write.csv(phen93272_, "~/GSE_merging_data/phen93272.csv")

phen93776_ <- as.data.frame(cbind(sample=rownames(phen93776), group=phen93776[,c(10,11)]))
colnames(phen93776_)[colnames(phen93776_) == "group.characteristics_ch1"] <- "Disease_state"
colnames(phen93776_)[colnames(phen93776_) == "group.characteristics_ch1.1"] <- "Gender"
phen93776_$platform <- "GPL570"
phen93776_$tissue <- "Blood"
write.csv(phen93776_, "~/GSE_merging_data/phen93776.csv")


phen65010_ <- as.data.frame(cbind(sample=rownames(sampleInfo_10), group=sampleInfo_10[,c(12,35)]))
colnames(phen65010_)[colnames(phen65010_) == "characteristics_ch1.2"] <- "Disease_state"
colnames(phen65010_)[colnames(phen65010_) == "gender"] <- "Gender"
phen65010_$platform <- "GPL570"
phen65010_$tissue <- "Blood"
write.csv(phen65010_, "~/GSE_merging_data/phen65010.csv")

phen55457_ <- as.data.frame(cbind(sample=rownames(phen55457), group=phen55457[,c(10,12)]))
colnames(phen55457_)[colnames(phen55457_) == "characteristics_ch1.2"] <- "Disease_state"
colnames(phen55457_)[colnames(phen55457_) == "Gender"] <- "Gender"
phen55457_$tissue <- "synovial_membrane"
phen55457_$platform <- "GPL96"
write.csv(phen55457_, "~/GSE_merging_data/phen55457.csv")


phen12021_ <- as.data.frame(cbind(sample=rownames(phen12021), group=phen12021[,c(10,12)]))
colnames(phen12021_)[colnames(phen12021_) == "characteristics_ch1.2"] <- "Disease_state"
colnames(phen12021_)[colnames(phen12021_) == "characteristics_ch1"] <- "Gender"
phen12021_$tissue <- "synovial_membrane"
phen12021_$platform <- "GPL96"
#phen12021_$tissue <- "Blood"
write.csv(phen12021_, "~/GSE_merging_data/phen12021.csv")

View(phen1919)
'''phen1919_ <- as.data.frame(cbind(sample=rownames(phen1919), group=phen1919[,c(8,11,12)]))
phen1919_$tissue <- "synovial_membrane"
phen1919_$platform <- "GPL91"
colnames(phen1919_)[colnames(phen1919_) == "source_name_ch1"] <- "Disease_state"
colnames(phen1919_)[colnames(phen1919_) == "characteristics_ch1.1"] <- "Gender1"
colnames(phen1919_)[colnames(phen1919_) == "characteristics_ch1.2"] <- "Gender2"
write.csv(phen1919_, "~/GSE_merging_data/phen1919.xlsx")'''

#read.csv("~/GSE_merging_data/phen1919.csv")
#data <- read.csv("/Users/swetarai/GSE_merging_data/phen1919_new1.csv")
#######merge it

# Load all CSV files
phen93272 <- read.csv("~/GSE_merging_data/phen93272.csv")
phen93776 <- read.csv("~/GSE_merging_data/phen93776.csv")
phen65010 <- read.csv("~/GSE_merging_data/phen65010.csv")
colnames(phen65010)[colnames(phen65010) == "group.characteristics_ch1.2"] <- "Disease_state"
colnames(phen65010)[colnames(phen65010) == "group.gender"] <- "Gender"

phen55457 <- read.csv("~/GSE_merging_data/phen55457.csv")
colnames(phen55457)[colnames(phen55457) == "group.characteristics_ch1.2"] <- "Disease_state"
colnames(phen55457)[colnames(phen55457) == "group.characteristics_ch1"] <- "Gender"
phen12021 <- read.csv("~/GSE_merging_data/phen12021.csv")
colnames(phen12021)[colnames(phen12021) == "group.characteristics_ch1.2"] <- "Disease_state"
colnames(phen12021)[colnames(phen12021) == "group.characteristics_ch1"] <- "Gender"
#phen1919 <- read.csv("~/GSE_merging_data/phen1919_new1.csv")  # Adjust if it's actually .csv



merged_data <- rbind(phen93272,phen93776)
merged_data1 <- rbind(merged_data,phen65010)
merged_data2 <- rbind(merged_data1,phen55457)
merged_data3 <- rbind(merged_data2,phen12021)


dim(merged_data3)



# Write the merged data to a new CSV file
write.csv(merged_data3, "~/GSE_merging_data/merged_phenotypes.csv", row.names = FALSE)



########################################

#### Load expression data ####
##############################

exp93272 <- as.data.frame(fread("~/GSE_merging_data/exp93272.csv"))
rownames(exp93272) <- exp93272$V1
#View(exp93776)
#exp93272 <- exp93272[,-1]

exp93776 <- as.data.frame(fread("~/GSE_merging_data/exp93776.csv"))
rownames(exp93776) <- exp93776$V1
#exp93776 <- exp93776[,-1]

exp65010 <- as.data.frame(fread("~/GSE_merging_data/exp65010.csv"))
rownames(exp65010) <- exp65010$V1
#exp65010 <- exp65010[,-1]

exp55457 <- as.data.frame(fread("~/GSE_merging_data/exp55457.csv"))
rownames(exp55457) <- exp55457$V1
#exp55457 <- exp55457[,-1]

exp12021 <- as.data.frame(fread("~/GSE_merging_data/exp12021.csv"))
rownames(exp12021) <- exp12021$V1
#exp12021 <- exp12021[,-1]

#exp1919 <- as.data.frame(fread("~/GSE_merging_data/exp1919.csv"))
#rownames(exp1919) <- exp1919$V1
#exp1919 <- exp1919[,-1]



intersectSymbol <- Reduce(intersect, 
                          list(
                            rownames(exp93272),
                            rownames(exp93776),
                            rownames(exp65010),
                            rownames(exp55457),
                            rownames(exp12021)
                            
                          )
)

dim(exp93272) # 20959   276
dim(exp93776)
dim(exp65010)
dim(exp55457)
dim(exp12021)
#dim(exp1919)# 20959   174


exp93272 <- exp93272[intersectSymbol,]
exp93776 <- exp93776[intersectSymbol,]
exp65010 <- exp65010[intersectSymbol,]
exp55457 <- exp55457[intersectSymbol,]
exp12021 <- exp12021[intersectSymbol,]
#exp1919 <- exp1919[intersectSymbol,]



exp93272 <- t(exp93272) %>% as.data.frame()
exp93776 <- t(exp93776) %>% as.data.frame()
exp65010 <- t(exp65010) %>% as.data.frame()
exp55457 <- t(exp55457) %>% as.data.frame()
exp12021 <- t(exp12021) %>% as.data.frame()

##############################
#phen38351 <- pData(gse38351)
#phen93272 <- pData(gse93272)
#phen93776 <- pData(gse93776)
#### Load phenotypic data ####
##############################

phen93272 <- as.data.frame(fread("~/GSE_merging_data/phen93272.csv", header = TRUE))
phen93272 <- phen93272[, c("V1", "sample", "Gender", "Disease_state", "platform", "tissue")]
exp93272 <- cbind(sample=rownames(exp93272), exp93272)
#colnames(exp93272)[colnames(exp93272) == "Gene_Symbol"] <- "sample"
exp93272 <- merge(phen93272, exp93272, by="sample")
rownames(exp93272) <- exp93272$sample
dim(exp93272)
table(phen93272$Disease_state)

phen93776 <- as.data.frame(fread("~/GSE_merging_data/phen93776.csv", header = TRUE))
phen93776 <- phen93776[, c("V1", "sample", "Gender", "Disease_state", "platform", "tissue")]
exp93776 <- cbind(sample=rownames(exp93776), exp93776)
#colnames(exp93776)[colnames(exp93776) == "Gene_Symbol"] <- "sample"
exp93776 <- merge(phen93776, exp93776, by="sample")
rownames(exp93776) <- exp93776$sample
dim(exp93776)
table(phen93776$Disease_state)

phen65010 <- as.data.frame(fread("~/GSE_merging_data/phen65010.csv", header = TRUE))

colnames(phen65010)[colnames(phen65010) == "group.gender"] <- "Gender"
colnames(phen65010)[colnames(phen65010) == "group.characteristics_ch1.2"] <- "Disease_state"
phen65010 <- phen65010[, c("V1", "sample", "Gender", "Disease_state", "platform", "tissue")]
exp65010 <- cbind(sample=rownames(exp65010), exp65010)
#colnames(phen65010)[colnames(phen65010) == "group.characteristics_ch1.2"] <- "Disease_state"
#colnames(phen65010)[colnames(phen65010) == "group.gender"] <- "Gender"
exp65010 <- merge(phen65010, exp65010, by="sample")
rownames(exp65010) <- exp65010$sample
View(exp65010)
table(phen65010$Disease_state)


phen55457 <- as.data.frame(fread("~/GSE_merging_data/phen55457.csv", header = TRUE))
colnames(phen55457)[colnames(phen55457) == "group.characteristics_ch1"] <- "Gender"
colnames(phen55457)[colnames(phen55457) == "group.characteristics_ch1.2"] <- "Disease_state"
phen55457 <- phen55457[, c("V1", "sample", "Gender", "Disease_state", "platform", "tissue")]
exp55457 <- cbind(sample=rownames(exp55457), exp55457)

exp55457 <- merge(phen55457, exp55457, by="sample")
rownames(exp55457) <- exp55457$sample
dim(exp55457)
table(phen55457$Disease_state)

phen12021 <- as.data.frame(fread("~/GSE_merging_data/phen12021.csv", header = TRUE))
colnames(phen12021)[colnames(phen12021) == "group.characteristics_ch1"] <- "Gender"
colnames(phen12021)[colnames(phen12021) == "group.characteristics_ch1.2"] <- "Disease_state"
phen12021 <- phen12021[, c("V1", "sample", "Gender", "Disease_state", "platform", "tissue")]
exp12021 <- cbind(sample=rownames(exp12021), exp12021)

exp12021 <- merge(phen12021, exp12021, by="sample")
rownames(exp12021) <- exp12021$sample
dim(exp12021)
table(phen55457$Disease_state)
##############################
dim(exp93272)
dim(exp93776)
dim(exp65010)
dim(exp55457)
dim(exp12021)
#### Merging datasets ####
##########################

exp93272 <- cbind(batch=rep(2, nrow(exp93272)), exp93272)
exp93776 <- cbind(batch=rep(3, nrow(exp93776)), exp93776)
exp65010 <- cbind(batch=rep(4, nrow(exp65010)), exp65010)
exp55457 <- cbind(batch=rep(5, nrow(exp55457)), exp55457)
exp12021 <- cbind(batch=rep(6, nrow(exp12021)), exp12021)

#exp93272$batch <- NULL

colnames(exp93776) <- colnames(exp93272)
colnames(exp65010) <- colnames(exp93272)
colnames(exp55457) <- colnames(exp93272)
colnames(exp12021) <- colnames(exp93272)


# Re-merge after reordering
#mergedPhen <- rbind(phen65010, phen55457, phen93272, phen93776, phen12021)

# Merge the datasets using rbind
mergedExpr <- rbind(exp93272, exp93776, exp65010, exp55457, exp12021)
dim(mergedExpr)

write.csv(mergedExpr, "~/GSE_merging_data/mergeDatasets_raw.csv")

View(mergedExpr)
table(mergedExpr$platform)

unique(mergedExpr$tissue)

######## this is for factor

mergedExpr$Disease_state <- as.factor(mergedExpr$Disease_state)
levels(mergedExpr$Disease_state) <- c("Control", "RA","OA")  # Original levels
levels(mergedExpr$Disease_state) <- c(0, 1, 2)  # 

mergedExpr$Gender <- as.factor(mergedExpr$Gender)
levels(mergedExpr$Gender) <- c("Female", "Male")  # Original levels
levels(mergedExpr$Gender) <- c(0, 1)  # 


mergedExpr$tissue <- as.factor(mergedExpr$tissue)
levels(mergedExpr$tissue) <- c("Blood", "synovial_membrane")  # Original levels
levels(mergedExpr$tissue) <- c(0, 1)  # 

mergedExpr$platform <- as.factor(mergedExpr$platform)
levels(mergedExpr$platform) <- c("GPL570", "GPL96")  # Original levels
levels(mergedExpr$platform) <- c(0, 1)  # 



mergedExpr$Disease_state <- as.numeric(mergedExpr$Disease_state)
mergedExpr$Gender <- as.numeric(mergedExpr$Gender)
mergedExpr$tissue <- as.numeric(mergedExpr$tissue)
mergedExpr$platform <- as.numeric(mergedExpr$platform)

write.csv(mergedExpr, "~/GSE_merging_data/mergeDatasets_with_levels_0_1.csv")
##########################
View(mergedExpr)

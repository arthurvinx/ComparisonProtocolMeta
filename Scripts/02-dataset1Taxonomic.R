# Clean global environment ####

rm(list = ls())

# Load libraries ####

library(data.table)
library(dplyr)
library(vegan)
library(Rmpfr)
library(ggplot2)
library(tidyr)

# Set directory ####

dir <- "../Data/"

# Load files and get general information ####

metadata <- as.data.frame(fread(paste0(dir, "Dataset1/d1Metadata.txt"),
                                stringsAsFactors = F, head = T, sep = "\t"))
colnames(metadata)[3] <- "ranks"
metadata$SK <- sapply(strsplit(metadata$ranks, "; "), "[", 1)
bysk <- metadata %>% group_by(SK) %>% summarise(total = sum(num_reads))
metadata <- metadata[order(metadata$num_reads, decreasing = T),]
rownames(metadata) <- NULL

# Reads after quality control

sum(metadata$num_reads)

# Negative control

metadata$num_reads[1]

# Distinct taxonomy identifiers and descriptions (including the negative control)

length(unique(metadata$taxid))
length(unique(metadata$ranks))

# Grouped by superkingdom

bysk

# Protocol results

taxP <- fread(paste0(dir, "Dataset1/d1ProtocolLCA.txt"), stringsAsFactors = F, head = F, sep = "\t")

# Assigned

nrow(taxP)-sum(grepl("^Unknown$", taxP$V2))

# Not assigned

sum(grepl("^Unknown$", taxP$V2))

# Grouped by superkingdom

table(sapply(strsplit(taxP$V2, ";"), "[", 1))

# Megan results

taxM <- fread(paste0(dir, "Dataset1/d1MeganLCA.txt"), stringsAsFactors = F, head = F)

# Assigned

unname(nrow(taxM)-table(sapply(strsplit(taxM$V2, ";"), "[", 2))[2])

# Not assigned

table(sapply(strsplit(taxM$V2, ";"), "[", 2))[2]

# Grouped by superkingdom

idx <- grepl(";Viruses;", taxM$V2)
sum(idx)

idx <- grepl(";Bacteria;", taxM$V2)
sum(idx)

idx <- grepl(";Archaea;", taxM$V2)
sum(idx)

idx <- grepl(";Eukaryota;", taxM$V2)
sum(idx)

# Remove not assigned

taxP <- taxP[!(grepl("^Unknown$", taxP$V2)),]
taxM <- taxM[!(grepl("^NCBI;Not assigned;$", taxM$V2)),]

# Metrics ####

# Example for perfect metrics

# tp <- 407511 # reads correctly assigned
# fp <- 0 # reads incorrectly assigned
# tn <- 99918 # negative control reads not assigned
# fn <- 0 # reads that should be assigned, but were not assigned

# sen <- tp/(tp+fn)
# spec <- tn/(tn+fp)
# ppv <- tp/(tp+fp)
# npv <- tn/(tn+fn)
# mcc <- ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

metricsP <- function(df, mdf, dRank = "s"){
  tp <- fp <- 0
  tn <- mdf[mdf$refseq=="NegativeCTRL",4]
  if(dRank == "s"){
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 8)
    df$drank <- sapply(strsplit(df$V2, ";"), "[", 7)
  }else if(dRank == "g"){
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 7)
    df$drank <- sapply(strsplit(df$V2, ";"), "[", 6)
  }else{
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 3)
    df$drank <- sapply(strsplit(df$V2, ";"), "[", 2)
  }
  df$drank[is.na(df$drank)] <- "FP"
  df$drank <- gsub("_", " ", df$drank)
  idx <- grepl("^NegativeCTRL", df$V1)
  tn <- tn - sum(idx)
  df <- df[!idx,]
  df$V1 <- gsub("^..._(.*)\\..*$", "\\1", df$V1)
  df$match <- mdf[match(df$V1, mdf$refseq), "drank"]
  idx <- df$drank==df$match
  tp <- sum(idx)
  fp <- sum(!idx)
  fn <- sum(mdf[mdf$refseq!="NegativeCTRL",4]) - (tp + fp)
  sen <- tp/(tp+fn)
  spec <- tn/(tn+fp)
  ppv <- tp/(tp+fp)
  npv <- tn/(tn+fn)
  #mcc <- ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc <- ((mpfr(tp, 128)*mpfr(tn, 128))-(mpfr(fp, 128)*mpfr(fn, 128)))/
    sqrt((mpfr(tp, 128)+mpfr(fp, 128))*(mpfr(tp, 128)+mpfr(fn, 128))*
           (mpfr(tn, 128)+mpfr(fp, 128))*(mpfr(tn, 128)+mpfr(fn, 128)))
  mcc <- as.numeric(mcc)
  res <- round(c(tp, tn, fp, fn, sen, spec, ppv, npv, mcc), digits = 4)
  names(res) <- c("tp", "tn", "fp", "fn", "sen", "spec", "ppv", "npv", "mcc")
  return(res)
}

metricsM <- function(df, mdf, dRank = "s"){
  tp <- fp <- 0
  tn <- mdf[mdf$refseq=="NegativeCTRL",4]
  if(dRank == "s"){
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 8)
  }else if(dRank == "g"){
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 7)
  }else{
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 3)
  }
  idx <- grepl("^NegativeCTRL", df$V1)
  tn <- tn - sum(idx)
  df <- df[!idx,]
  df$V1 <- gsub("^..._(.*)\\..*$", "\\1", df$V1)
  df$match <- mdf[match(df$V1, mdf$refseq), "drank"]
  aux <- strsplit(unlist(df[, 2]), ";")
  df$bool <- vapply(1:nrow(df), function(i){
    return(unlist(df[i, "match"]) %in% aux[[i]])
  }, logical(1))
  tp <- sum(df$bool)
  fp <- sum(!df$bool)
  fn <- sum(mdf[mdf$refseq!="NegativeCTRL",4]) - (tp + fp)
  sen <- tp/(tp+fn)
  spec <- tn/(tn+fp)
  ppv <- tp/(tp+fp)
  npv <- tn/(tn+fn)
  #mcc <- ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc <- ((mpfr(tp, 128)*mpfr(tn, 128))-(mpfr(fp, 128)*mpfr(fn, 128)))/
    sqrt((mpfr(tp, 128)+mpfr(fp, 128))*(mpfr(tp, 128)+mpfr(fn, 128))*
           (mpfr(tn, 128)+mpfr(fp, 128))*(mpfr(tn, 128)+mpfr(fn, 128)))
  mcc <- as.numeric(mcc)
  res <- round(c(tp, tn, fp, fn, sen, spec, ppv, npv, mcc), digits = 4)
  names(res) <- c("tp", "tn", "fp", "fn", "sen", "spec", "ppv", "npv", "mcc")
  return(res)
}

rank <- "s" # Change to "s" for species, "g" for genus, and "p" for phylum ####

resP <- metricsP(taxP, metadata, dRank = rank)
resP

resM <- metricsM(taxM, metadata, dRank = rank)
resM

# Get true positives ####

getPtps <- function(df, mdf, dRank = "s"){
  if(dRank == "s"){
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 8)
    df$drank <- sapply(strsplit(df$V2, ";"), "[", 7)
  }else if(dRank == "g"){
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 7)
    df$drank <- sapply(strsplit(df$V2, ";"), "[", 6)
  }else{
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 3)
    df$drank <- sapply(strsplit(df$V2, ";"), "[", 2)
  }
  df$drank[is.na(df$drank)] <- "FP"
  df$drank <- gsub("_", " ", df$drank)
  idx <- grepl("^NegativeCTRL", df$V1)
  df <- df[!idx,]
  df$V1 <- gsub("^..._(.*)\\..*$", "\\1", df$V1)
  df$match <- mdf[match(df$V1, mdf$refseq), "drank"]
  res <- df[df$drank==df$match,]
  return(res)
}

getMtps <- function(df, mdf, dRank = "s"){
  if(dRank == "s"){
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 8)
  }else if(dRank == "g"){
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 7)
  }else{
    mdf$drank <- sapply(strsplit(mdf$ranks, "; "), "[", 3)
  }
  idx <- grepl("^NegativeCTRL", df$V1)
  df <- df[!idx,]
  df$V1 <- gsub("^..._(.*)\\..*$", "\\1", df$V1)
  df$match <- mdf[match(df$V1, mdf$refseq), "drank"]
  aux <- strsplit(unlist(df[, 2]), ";")
  df$bool <- vapply(1:nrow(df), function(i){
    return(unlist(df[i, "match"]) %in% aux[[i]])
  }, logical(1))
  res <- df[df$bool,]
  return(res)
}

rank <- "s" # change ####

resP <- getPtps(taxP, metadata, dRank = rank)
resM <- getMtps(taxM, metadata, dRank = rank)

aux <- metadata %>% group_by(ranks) %>% summarise(total = sum(num_reads)) %>% arrange(desc(total))
aux <- aux[-1,] # controle negativo
temp <- strsplit(aux$ranks, "; ")
aux$s <- sapply(temp, "[", 8)
aux$g <- sapply(temp, "[", 7)
aux$p <- sapply(temp, "[", 3)
s <- aux %>% group_by(s) %>% summarise(total = sum(total)) %>% arrange(desc(total))
g <- aux %>% group_by(g) %>% summarise(total = sum(total)) %>% arrange(desc(total))
p <- aux %>% group_by(p) %>% summarise(total = sum(total)) %>% arrange(desc(total))
g <- g[-which(g$g=="Unknown genus"),]
p <- p[-which(p$p=="Unknown phylum"),]
rm(aux, temp)

# plots ####

resP <- as.data.frame(table(resP$match))
resM <- as.data.frame(table(resM$match))

if(rank=="s"){
  df <- s
}else if(rank=="g"){
  df <- g
}else{
  df <- p
}
df <- as.data.frame(df)
colnames(df)[2] <- "Metadata"
if(rank=="s"){
  df$Protocol <- resP[match(df$s, resP$Var1), "Freq"]
  df$Megan <- resM[match(df$s, resM$Var1), "Freq"]
}else if(rank=="g"){
  df$Protocol <- resP[match(df$g, resP$Var1), "Freq"]
  df$Megan <- resM[match(df$g, resM$Var1), "Freq"]
}else{
  df$Protocol <- resP[match(df$p, resP$Var1), "Freq"]
  df$Megan <- resM[match(df$p, resM$Var1), "Freq"]
}
df$Protocol[is.na(df$Protocol)] <- 0
df$Megan[is.na(df$Megan)] <- 0
n <- nrow(df)
if(rank=="s"){
  df$s <- NULL
}else if(rank=="g"){
  df$g <- NULL
}else{
  df$p <- NULL
}
df <- df %>% gather("Source", "value")
df$x <- rep(1:n, 3)

# change ####
ylim <- 165000 # 16000 18000 165000
legend <- "Phylum" # Species Genus Phylum

# top15
temp <- df[df$x %in% 1:15, ]
p <- ggplot(temp, aes(x = x, y = value, fill = Source)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(x = legend, y = "Count", title = paste0("Count per ", tolower(legend))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, ylim), breaks = seq.int(0, ylim, 1000))
plot(p)

c <- c("#00ba38", "#619cff", "#f8766d")
temp <- df[df$Source=="Metadata", ]
p <- ggplot(temp, aes(x = x, y = value)) +
  geom_bar(position = "dodge", stat = "identity", colour = c[1], fill = c[1]) +
  labs(x = legend, y = "Count", title = paste0("Metadata ", tolower(legend), " count")) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, ylim), breaks = seq.int(0, ylim, 1000))
plot(p)

temp <- df[df$Source=="Protocol", ]
p <- ggplot(temp, aes(x = x, y = value)) +
  geom_bar(position = "dodge", stat = "identity", colour = c[2], fill = c[2]) +
  labs(x = legend, y = "Count", title = paste0("Protocol ", tolower(legend), " count")) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, ylim), breaks = seq.int(0, ylim, 1000))
plot(p)

temp <- df[df$Source=="Megan", ]
p <- ggplot(temp, aes(x = x, y = value)) +
  geom_bar(position = "dodge", stat = "identity", colour = c[3], fill = c[3]) +
  labs(x = legend, y = "Count", title = paste0("Megan ", tolower(legend), " count")) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, ylim), breaks = seq.int(0, ylim, 1000))
plot(p)

metadata <- metadata[-1,]
metadata <- metadata %>% group_by(ranks) %>% summarise(total = sum(num_reads)) %>% arrange(desc(total))

taxM <- as.data.frame(table(taxM$V2))
taxP <- as.data.frame(table(taxP$V2))

v <- taxP$Freq

H <- diversity(v)
EH <- H/log(specnumber(v))

round(H, digits = 4)
round(EH, digits = 4)

#### Calculate the score 
#!/usr/bin/env Rscript
setwd("/home/zdzhao/mulit_omics_result/result")
library(data.table)
library(reshape2)

### Load heritability and size files for each part
dat.all <- read.table('demo.h2.sum',header=F,sep="\t")
size <- read.table('demo.size.grm',header=F)

#### Merge size file and heritability file, then rename the columns
colnames(size)<-c('V2','V1')
dat1 <- merge(size,dat.all,by='V2')
colnames(dat1) <- c('grm','size','tr','h2','se')

# Create a new table with GRM as rows, trait names as columns, and heritability values as entries
dat2 <- merge(unique(dat1[,1:2]),data.table::dcast(dat1[,c(1,3,4)],grm~tr),by='grm')

### Load the heritability of each part for each trait
dat.all1 <- dcast(dat.all[!grepl('hd',dat.all[,1]),],V1~V2,value.var='V3')

all.dat1 <- dat2[,-(1:2)]/dat2$size
View(all.dat1)
all.dat2 <- dat2[,1:2]
View(all.dat2)
all.dat3 <- cbind(all.dat2,all.dat1)
View(all.dat3)
all.dat3$ave.h2 <- rowMeans(all.dat3[,-(1:2)],na.rm=T)

listpath <- '/home/zdzhao/mulit_omics_result/result/snplist/'
sllist <- list()

### Load all snplist files
for (i in list.files(path=listpath,pattern='*.snplist'))
{sllist[[i]] <- fread(paste0(listpath,i),header=F)
}
# Load the list of all SNPs
faethd <- sllist[['allsnp.snplist']]
colnames(faethd)[1] <- 'SNP'

### Load all snplist files again
for (i in list.files(path=listpath,pattern='*.snplist'))
{sllist[[i]] <- fread(paste0(listpath,i),header=F)
}
# Load the list of all SNPs
faethd <- sllist[['allsnp.snplist']]
colnames(faethd)[1] <- 'SNP'

ihs.list <- c("ihs", "noihs")

ld.list <- c("ld1", "ld2", "ld3", "ld4")

### Process each list in ihs.list
for (fd in c('ihs.list')) {
    grmtypelist <- list()
    for (grmtype in get(fd)) {
        grmtypelist[[grmtype]] <- cbind(sllist[[which(names(sllist) == paste0(grmtype, '.snplist'))]], all.dat3[all.dat3$grm == grmtype,]$ave.h2)
    }
    tmp <- do.call(rbind, grmtypelist)
    colnames(tmp) <- c('SNP', gsub('.list', '', fd))
    setkey(tmp, SNP)
    fdlist[[fd]] <- tmp
}

### Process each list in ld.list
for (fd in c('ld.list')) {
    grmtypelist <- list()
    for (grmtype in get(fd)) {
        grmtypelist[[grmtype]] <- cbind(sllist[[which(names(sllist) == paste0(grmtype, '.snplist'))]], all.dat3[all.dat3$grm == grmtype,]$ave.h2)
    }
    tmp <- do.call(rbind, grmtypelist)
    colnames(tmp) <- c('SNP', gsub('.list', '', fd))
    setkey(tmp, SNP)
    fdlist[[fd]] <- tmp
}

### Merge all data frames
faethd1 <- Reduce(function(...) merge(..., all = T, by = "SNP"), c(list(faethd), fdlist))
faethd1[, faeth.sc := rowMeans(faethd1[,-1])]

### Write the final score to a file
fwrite(faethd1, "demo_score.txt", sep = "\t")

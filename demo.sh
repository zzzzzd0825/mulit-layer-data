###calculate the h2 of target set
#!/bin/bash
for file in *snplist; do
    cd "/home/zdzhao/mulit_omics_result/${file%.snplist}"
    echo "${file%.snplist}" "no_${file%.snplist}" > "${file%.snplist}".mgrm
    gcta64 --bfile /home/zdzhao/1577_cattle/1577 --extract "$file" --make-grm --out "${file%.snplist}"
    gcta64 --bfile /home/zdzhao/1577_cattle/1577 --exclude "$file" --make-grm --out "no_${file%.snplist}"
    gcta64 --mgrm "${file%.snplist}".mgrm --reml --pheno /home/zdzhao/phenotype/ph.txt --mpheno $i --out "tr$i" --reml-alg 2 --threads 40 --reml-maxit 10000
done
### Demo of preparing a document
for file in tr*.hsq; do     sed -n '8,11p' "$file" | awk -v filename="$file" '{print filename, $0 >> "ld.hsq"}'; done
sed -i 's/\s\+/\t/g' ld.hsq 
#####change the names
awk 'BEGIN { FS=OFS="\t" } { 
    if ($2 == "V(G1)/Vp") $2 = "LD1"
    else if ($2 == "V(G2)/Vp") $2 = "LD2"
    else if ($2 == "V(G3)/Vp") $2 = "LD3"
    else if ($2 == "V(G4)/Vp") $2 = "LD4"
    print
}' ld.hsq >ld.h2_sum

###calculate the number of target set
for file in *snplist; do
size=$(wc -l < "$file")
echo -e "${file%.snplist}\t$size"
done > demo.size.grm



####calculate the score 
#!/usr/bin/env Rscript
setwd("/home/zdzhao/mulit_omics_result/result")
library(data.table)
library(reshape2)
###读入遗传力和每个part的size文件
dat.all <- read.table('demo.h2.sum',header=F,sep="\t")
size <- read.table('demo.size.grm',header=F)
####合并size文件和遗传力文件，并对其进行重命名
colnames(size)<-c('V2','V1')
dat1 <- merge(size,dat.all,by='V2')
colnames(dat1) <- c('grm','size','tr','h2','se')
#根据grm为行，性状名字为列，遗传力为数值新建表
dat2 <- merge(unique(dat1[,1:2]),data.table::dcast(dat1[,c(1,3,4)],grm~tr),by='grm')
###载入每个part在每个性状中所占的遗传力
dat.all1 <- dcast(dat.all[!grepl('hd',dat.all[,1]),],V1~V2,value.var='V3')

all.dat1 <- dat2[,-(1:2)]/dat2$size
View(all.dat1)
all.dat2<-dat2[,1:2]
View(all.dat2)
all.dat3<-cbind(all.dat2,all.dat1)
View(all.dat3)
all.dat3$ave.h2 <- rowMeans(all.dat3[,-(1:2)],na.rm=T)

listpath <- '/home/zdzhao/mulit_omics_result/result/snplist/'
sllist <- list()

for (i in list.files(path=listpath,pattern='*.snplist'))
{sllist[[i]] <- fread(paste0(listpath,i),header=F)
}
#载入allsnp的列表
faethd <- sllist[['allsnp.snplist']]
colnames(faethd)[1] <- 'SNP'

for (i in list.files(path=listpath,pattern='*.snplist'))
{sllist[[i]] <- fread(paste0(listpath,i),header=F)
}
#载入allsnp的列表
faethd <- sllist[['allsnp.snplist']]
colnames(faethd)[1] <- 'SNP'


ihs.list<-c("ihs","noihs")

ld.list<-c("ld1","ld2","ld3","ld4")

for (fd in c('ihs.list'))
{grmtypelist <- list()
for (grmtype in get(fd))
{grmtypelist[[grmtype]] <- cbind(sllist[[which(names(sllist)==paste0(grmtype,'.snplist'))]],all.dat3[all.dat3$grm==grmtype,]$ave.h2)
}
tmp <- do.call(rbind,grmtypelist)
colnames(tmp) <- c('SNP',gsub('.list','',fd))
setkey(tmp,SNP)
fdlist[[fd]] <- tmp
}
for (fd in c('ld.list'))
{grmtypelist <- list()
for (grmtype in get(fd))
{grmtypelist[[grmtype]] <- cbind(sllist[[which(names(sllist)==paste0(grmtype,'.snplist'))]],all.dat3[all.dat3$grm==grmtype,]$ave.h2)
}
tmp <- do.call(rbind,grmtypelist)
colnames(tmp) <- c('SNP',gsub('.list','',fd))
setkey(tmp,SNP)
fdlist[[fd]] <- tmp
}

faethd1 <-  Reduce(function(...) merge(...,all=T,by="SNP"),c(list(faethd),fdlist))
faethd1[,faeth.sc:=rowMeans(faethd1[,-1])]
 
fwrite(faethd1, "demo_score.txt", sep = "\t")

### Selection of variants based on thresholds
for i in $(seq 5 10 30); do
    sort -k2 -n demo_score.txt > sorted_scores.txt
    line=$(wc -l < sorted_scores.txt)
    line_bottom_i=$((i * line / 100))
    line_top_i=$(((100 - i) * line / 100))
    awk -v line_i="$line_bottom_i" 'NR < line_i {print $1}' sorted_scores.txt > score_bottom_$i.snplist
    awk -v line_i="$line_top_i" -v line="$line" 'NR > line_i {print $1}' sorted_scores.txt > score_top_$i.snplist
done


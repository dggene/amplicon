vars.tmp <- commandArgs(T)
split.vars  <- vars.tmp[-(1:2)]
input =split.vars[1]
inputindel=split.vars[2]
region=split.vars[3]
out=split.vars[4]
snpbed=split.vars[5]
indelbed=split.vars[6]

library(VariantAnnotation)
#print(sessionInfo())
vcf <- readVcf(input, "hg19")
#print(sessionInfo())
snp <- read.table(snpbed,head=FALSE)[,-2]
colnames(snp)=c('chr','end','RS')

vcfhead=as.data.frame(rowRanges(vcf))[,c(1:3,7)]
geno <- geno(vcf)$GT
ref <- ref(vcf)
alt <- alt(vcf)
dp <- info(vcf)$DP
#print('2')

geno2 <- geno
for (i in 1:nrow(geno)) {
   print(i)
   geno2[i,] <- gsub("0", as.character(ref[i]), geno[i,])
    geno2[i,] <- gsub("/",'', geno2[i,])
   for (j in 1:elementNROWS(alt[i])) {
     geno2[i,] <- gsub(as.character(j),
                       as.character(alt[[i]][j]),
                       geno2[i,])
   }
}

genofinal=cbind(vcfhead,geno2)
genofinal=transform(genofinal,depth=dp,correctgeno=genofinal[,5],check.names=F)

genofinal=merge(genofinal,snp,by.x=c('seqnames','end'),by.y=c('chr','end'),all.x=TRUE)
name=colnames(genofinal)
name[7]=paste(name[5],'cor',sep='')
colnames(genofinal)=name
print('#....$')
print('$....$')

#save.image('test.Rdata')
#indel annotation
vcfindel <- readVcf(inputindel, "hg19")
indel <- read.table(indelbed,head=F)
vcfindelhead=as.data.frame(rowRanges(vcfindel))[,c(1:3,7)]
genoindel <- geno(vcfindel)$GT
refindel <- ref(vcfindel)
altindel <- alt(vcfindel)
dp <- info(vcfindel)$DP
print('$$....$')

#indel
indelfinal <-NULL
region.data <- read.table(region,check.names=F,head=T)
indel.region <- as.character(indel[,10])
indel.dp <- region.data[which(region.data[,1]%in%indel.region),c(1,3)]
indel <- merge(indel,indel.dp,by.x='V10',by.y='Target',all.x=T,sort=F)
indel <-indel[,c(2,3,3,9,10,11,10,4)]
#indel.dp <-rep(0,nrow(indel))
indelfinal= indel[-1,]

indelfinal=as.matrix(indelfinal)
colnames(indelfinal)= name
geno3 <- genoindel
if(nrow(vcfindelhead)!=0){

 for( i in 1:nrow(vcfindelhead)) {
      print(i)
     index=as.numeric(vcfindelhead[i,2])
      for(j in 2:3){

               pos=as.numeric(as.character(indel[j,2]))
               if(index>(pos-10)& index <(pos+10)){
                  if(genoindel[i,]=='0/0')  {indelfinal[j-1,7]=as.character(indel[j,5])}
                  if(genoindel[i,]=='0/1')  {indelfinal[j-1,7]=as.character(indel[j,6])}
                  if(genoindel[i,]=='1/1')  {indelfinal[j-1,7]=as.character(indel[j,7])}
                  indelfinal[j-1,6]=dp[i]
                  geno3[i,] <- gsub("0", as.character(refindel[i]), genoindel[i,])
                  for (k in 1:elementNROWS(altindel[i])) {
                  geno3[i,] <- gsub(as.character(k),as.character(altindel[[i]][k]),geno3[i,])
                    }
                  indelfinal[j-1,5] = geno3[i,]
                  indelfinal[j-1,4] = as.character(refindel[i,])
                   }
                  }
           }
   }

####gst
gst=region.data
gstmup=gst[grep('chr1:110232921',gst[,1]),c(2,9)]
gsttup=gst[grep('chr22:24376736',gst[,1]),c(2,9)]

gstgeno=data.frame(chr=c('chr1','chr22'),start=c('110232921','24376736'),end=c('GSTM1','GSTT1'),ref=c('-','-'),geno=c('present','present'),depth=c(mean(gstmup[1,1],gstmup[1,2]),mean(gsttup[1,1],gsttup[1,2])),geno=c('present','present'),rs=c('DG_001','DG_002'))
gstgeno=as.matrix(gstgeno)

if(gstmup[1,2]<100) {gstgeno[1,5]= 'NULL';gstgeno[1,7]= 'NULL'}
if(gsttup[1,2]<100) {gstgeno[2,5]= 'NULL';gstgeno[2,7]= 'NULL'}
colnames(gstgeno)=name

####ace
ACE=as.matrix(indel[1,])
#ACE=data.frame(indel[1,1:2], indel[1,c(2,8)], indel[1,9],dp=100,indel[1,5],indel[1,3])
#ACE=as.matrix(data.frame(indel[1,1:2], indel[1,c(2,8)], indel[1,9],dp=100,indel[1,5],indel[1,3]))

acegeno=genofinal[which(genofinal$RS%in%'rs4343'),7]
if(acegeno=='GG') {
            ACE[1,5]='DD'
            ACE[1,7]='DD'
            ACE[1,6]=genofinal[which(genofinal$RS%in%'rs4343'),6]
                   }
if(acegeno=='GA') {
            ACE[1,5]='DI'
            ACE[1,7]='DI'
            ACE[1,6]=genofinal[which(genofinal$RS%in%'rs4343'),6]
            }
if(acegeno=='AA') {
            ACE[1,5]='II'
            ACE[1,7]='II'
            ACE[1,6]=genofinal[which(genofinal$RS%in%'rs4343'),6]
            }
colnames(ACE)=name


########
final=rbind(genofinal,ACE,indelfinal,gstgeno)
missgeno=final[is.na(final$depth),]
if(nrow(missgeno)!=0) {missgeno[,7]='NA'}

#depth <10,SNP IS na
allgeno=final[!is.na(final$depth),]
allgeno[,7]=as.character(allgeno[,7])
for(i in 1:nrow(allgeno)){

    if(allgeno[i,7]=='.'){
         allgeno[i,7] ='NA'
         }
    if(as.numeric(allgeno[i,6])<4){
         allgeno[i,7] ='NA'
         }
  }

snpall=rbind(allgeno,missgeno)

#if(snpall[grep('rs676210',snpall[,8]),7]=='GA') snpall[grep('rs676210',snpall[,8]),7]='AG'

write.table(snpall[,c(8,7,6,1:5)],paste(out,'.geno',sep=''),quote=F,row.names=F,sep='\t')

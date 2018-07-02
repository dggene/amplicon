vars.tmp <- commandArgs(T)
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,','))
input =split.vars[1]
inputindel=split.vars[2]
inputgst=split.vars[3]
out=split.vars[4]

library(VariantAnnotation)
#snp annotation
#vcf <- readVcf('d:/Sample10.snp.vcf', "hg19")

vcf <- readVcf(input, "hg19")
snp <- read.table('/DG/programs/beta/scripts/DNA/genetest/adult/dgadultsnp.bed',head=FALSE)[,-2]
colnames(snp)=c('chr','end','RS')
#save.image('test.Rdata')

vcfhead=as.data.frame(rowRanges(vcf))[,c(1:3,7)]
geno <- geno(vcf)$GT
ref <- ref(vcf)
alt <- alt(vcf)
dp <- info(vcf)$DP

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


missgeno=genofinal[is.na(genofinal$depth),]
if(nrow(missgeno)!=0){
   missgeno[,7]='NA'}


#depth <10,SNP IS na
snpgeno=genofinal[!is.na(genofinal$depth),]
snpgeno[,7]=as.character(snpgeno[,7])
for(i in 1:nrow(snpgeno)){
    if(snpgeno[i,7]=='.'){
         snpgeno[i,7] ='NA'}
    if(snpgeno[i,6]<20){
         snpgeno[i,7] ='NA'}
   }
snpall=rbind(snpgeno,missgeno)


#indel annotation
#vcfindel <- readVcf('d:/Sample10.indel.vcf', "hg19")
vcfindel <- readVcf(inputindel, "hg19")

indel <- read.table('/DG/programs/beta/scripts/DNA/genetest/adult/dgadultindel.txt',head=F)
vcfindelhead=as.data.frame(rowRanges(vcfindel))[,c(1:3,7)]
genoindel <- geno(vcfindel)$GT
refindel <- ref(vcfindel)
altindel <- alt(vcfindel)
dp <- info(vcfindel)$DP

indelfinal <-NULL
indelfinal= data.frame(indel[-1,1:2], indel[-1,c(2,8)], indel[-1,9],dp=rep('0',nrow(indel)-1),indel[-1,5],indel[-1,3])
ACE=data.frame(indel[1,1:2], indel[1,c(2,8)], indel[1,9],dp=100,indel[1,5],indel[1,3])
ACE=as.matrix(data.frame(indel[1,1:2], indel[1,c(2,8)], indel[1,9],dp=100,indel[1,5],indel[1,3]))

indelfinal=as.matrix(indelfinal)
colnames(indelfinal)= name
geno3 <- genoindel
if(nrow(vcfindelhead)!=0){ 
                 
 for( i in 1:nrow(vcfindelhead)) {
      print(i)
     index=as.numeric(vcfindelhead[i,2])     
      for(j in 1:3){
              
               pos=as.numeric(as.character(indel[j,2]))
               if(index>(pos-10)& index <(pos+10)){ 
                  if(genoindel[i,]=='0/0')  {indelfinal[j-1,7]=as.character(indel[j,5])}
                  if(genoindel[i,]=='0/1')  {indelfinal[j-1,7]=as.character(indel[j,6])}
                  if(genoindel[i,]=='1/1')  {indelfinal[j-1,7]=as.character(indel[j,7])}
                  indelfinal[j-1,6]=dp[i]
                  geno3[i,] <- gsub("0", as.character(refindel[i]), genoindel[i,])
                  for (k in 1:elementLengths(altindel[i])) {
                  geno3[i,] <- gsub(as.character(k),as.character(altindel[[i]][k]),geno3[i,])
                    }
                  indelfinal[j-1,5] = geno3[i,] 
                  indelfinal[j-1,4] = as.character(refindel[i,]) 
                   }
                  }  
           }            
   }  

#gstm,gstt
gst=read.table(inputgst,head=T)
gstmup=gst[grep('chr1:110232921',gst[,1]),c(2,9)]

gsttup=gst[grep('chr22:24376736',gst[,1]),c(2,9)]

gstgeno=data.frame(chr=c('chr1','chr22'),start=c('110232921','24376736'),end=c('GSTM1','GSTT1'),ref=c('-','-'),geno=c('present','present'),depth=c(mean(gstmup[1,1],gstmup[1,2]),mean(gsttup[1,1],gsttup[1,2])),geno=c('present','present'),rs=c('DG_001','DG_002'))
gstgeno=as.matrix(gstgeno)

if(gstmup[1,2]<100) {gstgeno[1,5]= 'NULL';gstgeno[1,7]= 'NULL'}
if(gsttup[1,2]<100) {gstgeno[2,5]= 'NULL';gstgeno[2,7]= 'NULL'}
colnames(gstgeno)=name

#acegeno=snpall[which(snpall$RS%in%'rs4343'),7]
#if(acegeno=='GG') {ACE[1,5]='DD';ACE[1,7]='DD'}
#if(acegeno=='GA') {ACE[1,5]='DI';ACE[1,7]='DI'}
#if(acegeno=='AA') {ACE[1,5]='II';ACE[1,7]='II'}
#colnames(ACE)=name

#final=rbind(snpall,ACE,indelfinal,gstgeno)

final=rbind(snpall,indelfinal,gstgeno)

write.table(final[,c(8,7,6,1:5)],paste(out,'.geno',sep=''),quote=F,row.names=F,sep='\t') 
   

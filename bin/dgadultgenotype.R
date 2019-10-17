vars.tmp <- commandArgs(T)
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,','))
input =split.vars[1]
inputindel=split.vars[2]
region=split.vars[3]
out=split.vars[4]
snp_bed=split.vars[5]
indel_bed=split.vars[6]

library(VariantAnnotation)
#snp annotation

print(sessionInfo())
save.image('test.Rdata')
vcf <- readVcf(input, "hg19")
#print(sessionInfo())
snp <- read.table(snp_bed,head=FALSE)[,-2]
colnames(snp)=c('chr','end','RS')

vcfhead=as.data.frame(rowRanges(vcf))[,c(1:3,7)]
geno <- geno(vcf)$GT
ref <- ref(vcf)
alt <- alt(vcf)
dp <- info(vcf)$DP
print('2')

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

indel <- read.table(indel_bed,head=F)
vcfindelhead=as.data.frame(rowRanges(vcfindel))[,c(1:3,7)]
genoindel <- geno(vcfindel)$GT
refindel <- ref(vcfindel)
altindel <- alt(vcfindel)
dp <- info(vcfindel)$DP
print('$$....$')

indelfinal <-NULL
region.data <- read.table(region,check.names=F,head=T)
indel.region <- c("chr9:136132901-136132915","chr22:42525081-42525090")
indel.dp <- region.data[which(region.data[,1]%in%indel.region),3]
indelfinal= data.frame(indel[-1,1:2], indel[-1,c(2,8)], indel[-1,9],dp=indel.dp,indel[-1,5],indel[-1,3])


ACE=data.frame(indel[1,1:2], indel[1,c(2,8)], indel[1,9],dp=100,indel[1,5],indel[1,3])
ACE=as.matrix(data.frame(indel[1,1:2], indel[1,c(2,8)], indel[1,9],dp=100,indel[1,5],indel[1,3]))

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

final=rbind(genofinal,ACE,indelfinal)

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

write.table(snpall[,c(8,7,6,1:5)],paste(out,'.geno',sep=''),quote=F,row.names=F,sep='\t') 

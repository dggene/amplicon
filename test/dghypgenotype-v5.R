vars.tmp <- commandArgs(T)
split.vars  <- vars.tmp[-(1:2)]
input =split.vars[1]
inputindel=split.vars[2]
region=split.vars[3]
out=split.vars[4]
snpbed=split.vars[5]
indelbed=split.vars[6]
specialbed=split.vars[7]

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

indelfinal <-NULL

region.data <- read.table(region,check.names=F,head=T)
indel.region <- as.character(indel[,10])
indel.dp <- region.data[which(region.data[,1]%in%indel.region),c(1,3)]
indel <- merge(indel,indel.dp,by.x='V10',by.y='Target',all.x=T )
indel <-indel[,c(2,3,3,9,10,11,10,4)]
#indel.dp <-rep(0,nrow(indel))
indelfinal= indel[-1,]

ACE=as.matrix(indel[1,])
#ACE=data.frame(indel[1,1:2], indel[1,c(2,8)], indel[1,9],dp=100,indel[1,5],indel[1,3])
#ACE=as.matrix(data.frame(indel[1,1:2], indel[1,c(2,8)], indel[1,9],dp=100,indel[1,5],indel[1,3]))



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

snpindel=rbind(genofinal,ACE,indelfinal)

special.geno=read.table(specialbed,head=FALSE)
keysnp=as.character(unique(special.geno$V1))
keysnp.geno=subset(snpindel,RS%in%keysnp)

temp<-apply(keysnp.geno,1,function(x){
    rs=x[8]
    rs.geno=x[7]
    sp=special.geno[special.geno$V1%in%x[8],]
    print(any(c(rs.geno=='.',rs.geno=='',rs.geno==' ',is.na(rs.geno))))
    if(any(c(rs.geno=='.',rs.geno=='',rs.geno==' ',is.na(rs.geno)))){
       x[7]='NA'  
       x[8]=as.character(unique(sp[,4]))
     }else{
        m=grep(rs.geno,sp$V2)
        #print(m)
        if(length(m)>0){
        x[7]=as.character(sp[m,3])
        x[8]=as.character(unique(sp[,4]))
        }else{
        x[7]='*1/*1'
        x[8]=as.character(unique(sp[,4])) 
        }
    }
    return(x)
  })

keysnp_m.geno=data.frame(t(temp))
snpindel=rbind(snpindel,keysnp_m.geno)

rs.merge=c('CYP2C19*1,CYP2C19*2,CYP2C19*3,CYP2C19*17',
  'CYP2C9*1,CYP2C9*2,CYP2C9*3',
  'CYP2D6*1,CYP2D6*2,CYP2D6*4,CYP2D6*6,CYP2D6*10,CYP2D6*41')

getgen=function(x){
  x1=unique(x)
  temp<-lapply(x1,function(i){
    if(sum(x%in%i)>=2){
      x2=rep(i,2)
    }else{
      x2=i
    }
    return(x2)
  })
  x=unlist(temp)
  return(x)
  }

temp<-lapply(rs.merge,function(x){
  x1=unlist(strsplit(x,','))[-1]
  y=as.character(keysnp_m.geno[keysnp_m.geno$RS %in% x1,7])
  g=unlist(strsplit(y,'/'))
  g=sort(as.numeric(gsub('*','',g,fixed=T)))
  print(g)
  g=paste('*',getgen(g),sep='')
  print(g)
  g=paste(tail(g,2),collapse='/') 
  res=data.frame(geno=g,RS=x)
})

special.merge=do.call(rbind,temp)

special.merge=data.frame(seqnames=rep('chrn',nrow(special.merge)),
          end=rep('1',nrow(special.merge)),
          start=rep('1',nrow(special.merge)),
          ref=rep('n',nrow(special.merge)),
          id=rep('n',nrow(special.merge)),
          depth=rep(100,nrow(special.merge)),
          special.merge
          )
colnames(special.merge)=name
snpindel=rbind(snpindel,special.merge)

snp=snpindel[grep('CYP2C9*1,CYP2C9*2,CYP2C9*3',snpindel$RS,fixed=T),]
snp1=snp[,7]
snp2=snpindel[grep('rs9923231',snpindel$RS,fixed=T),7]
snp[,8]='CYP2C9*1,CYP2C9*2,CYP2C9*3|rs9923231'
snp[,7]=paste(snp1,snp2,sep='|')

final=rbind(snpindel,snp)

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

#设置镜像，这部分代码复制粘贴运行就可以
local({
    r <- getOption("repos")  
    r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
    options(repos=r)
}) 
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
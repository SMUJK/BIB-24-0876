rm(list=ls())
library(dplyr)
library(arrow)
library(plink2R)
library(PDSCE)
source("/bib_function.R")
samples<-read.csv("/samples.csv")
for (m in 1:14) {
  tissue<-samples$tis[m]
  every_tissue<-list()
  every_tissue_ivw_i<-list()
  every_tissue_ivw_c<-list()
  for (chromsome in 1:22){
    print(chromsome)
    #loading ADHD gwas data, which are divided into 22 data of 22 chromosomes
    wd2<-paste(c("/gwas_ch",chromsome,".csv"),collapse="")
    gwas.ch<-read.csv(wd2)
    gwas.ch<-subset(gwas.ch,FRQ_A_38691>=0.01 & FRQ_U_186843>=0.01)
    gc()
    #this is a lookup table downloaded from GTEx portal, which is used for referring to the snp id of the GTEx QTL data
    wd3<-paste(c("/lookup_ch",chromsome,".csv"),collapse="")
    lookup.ch<-read.csv(wd3,header = F)
    lookup.ch<-lookup.ch[,c(1,7)]
    colnames(lookup.ch)<-c("variant_id","SNP")
    gc()
    qtl_files<-list.files(paste0("/eqtl/",tissue))
    wd1<-grepl(paste0("chr",chromsome,".par"),qtl_files) %>% qtl_files[.]
    eqtl.ch<-read_parquet(paste0("/eqtl/",tissue,"/",wd1))
    eqtl.ch<-as.data.frame(eqtl.ch)
    eqtl.ch<-subset(eqtl.ch,maf >= 0.01 & pval_nominal<0.05)
    eqtl.ch<-eqtl.ch[is.na(eqtl.ch$slope)==F,]
    gc()
    merged<-merge(eqtl.ch, lookup.ch, by = "variant_id",all=F)
    merged<-subset(merged,SNP!=".")
    smr<-merge(merged,gwas.ch,by="SNP",all=F)
    smr$z1<-smr$slope/smr$slope_se
    smr$z2<-log(smr$OR)/smr$SE
    genes<-unique(smr$phenotype_id)
    
    ld_dir<-"/1000G.EUR."
    geno<-read_plink(paste0(ld_dir,chromsome),impute="avg")
    refer<-scale(geno$bed)
    common_snp<-intersect(unique(smr$SNP),colnames(refer))
    refer<-refer[,common_snp]
    every_chr<-list()
    every_chr_ivw_i<-list()
    every_chr_ivw_c<-list()
    for (ppp in 1:length(genes)){
      #ppp<-1
      if (ppp%%100==0){
        print(paste(chromsome,":",ppp,"/",length(genes)))
      }
      gene<-genes[ppp]
      gdata<-smr[smr$phenotype_id==gene,]
      gdata$tsmr<-(((gdata$z1)^2)*((gdata$z2)^2))/(((gdata$z1)^2)+((gdata$z2)^2))
      gdata$psmr<-1-pchisq(gdata$tsmr,df=1)
      gdata$causal_b<-log(gdata$OR)/gdata$slope
      gdata$causal_se<-sqrt(((gdata$causal_b)^2)*((gdata$slope_se/gdata$slope)^2+(gdata$SE/log(gdata$OR))^2))
      gdata0<-gdata
      snplist<-intersect(gdata$SNP,common_snp)
      if (length(snplist)<=1 | nrow(gdata)==0){
        ivw_result_i<-data.frame(phenotype=gdata0$phenotype_id[1],
                               causal_b=gdata0$causal_b[1],
                               causal_se=gdata0$SE[1]/gdata0$slope[1])
        ivw_result_i$causal_p<-1-pchisq((ivw_result_i$causal_b/ivw_result_i$causal_se)^2,df=1)
      } else {
        gdata<-gdata[gdata$SNP %in% snplist,]
        refer_sub<-refer[,snplist]
        ld.mat <- t(refer_sub) %*% refer_sub / (nrow(refer_sub)-1)
        gdata<-gdata[order(gdata$pval_nominal,decreasing = F),]
        ld.mat<-ld.mat[gdata$SNP,gdata$SNP]
        clump.ivw<-clump(ld.mat,threshold = 0.4)
        if (length(clump.ivw$SNP)==1) {
          ivw_result_i<-data.frame(phenotype=gdata0$phenotype_id[1],
                                   causal_b=gdata0$causal_b[1],
                                   causal_se=gdata0$SE[1]/gdata0$slope[1])
          ivw_result_i$causal_p<-1-pchisq((ivw_result_i$causal_b/ivw_result_i$causal_se)^2,df=1)
          ivw_result_c<-ivw_result_i
        } else {
          gdata_ivw<-gdata[gdata$SNP %in% clump.ivw$SNP,]
          ld_ivw<-clump.ivw$LD[gdata_ivw$SNP,gdata_ivw$SNP]
          #ivw-i
          b_ivw_i<-sum(gdata_ivw$slope*log(gdata_ivw$OR)/((gdata_ivw$SE)^2))/sum((gdata_ivw$slope/gdata_ivw$SE)^2)
          se_ivw_i<-sqrt(1/sum((gdata_ivw$slope/gdata_ivw$SE)^2))
          p_ivq_i<-1-pchisq(((b_ivw_i/se_ivw_i)^2),df=1)
          ivw_result_i<-data.frame(phenotype=gdata_ivw$phenotype_id[1],
                                   causal_b=b_ivw_i,
                                   causal_se=se_ivw_i,
                                   causal_p=p_ivq_i)
        }
        
        
      }
      every_chr_ivw_i[[ppp]]<-ivw_result_i
    }
    every_tissue_ivw_i[[chromsome]]<-do.call(rbind,every_chr_ivw_i)
  }
  out_file<-"/output/"
  single_tissue_result_ivw_i<-do.call(rbind,every_tissue_ivw_i)
  write.csv(single_tissue_result_ivw_i,
            paste0(out_file,"etwas_ivw/",tissue,"_etwas_ivw.csv"),row.names = F)
  rm(every_tissue)
  gc()
}




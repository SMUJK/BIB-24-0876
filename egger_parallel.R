rm(list=ls())
library(doParallel)
source("/bib_function.R")
samples<-read.csv("/samples.csv")
cl <- makeCluster(6) #set number of cluster
registerDoParallel(cl)
for (m in 1:14) {
  tissue<-samples$tis[m]  
  all_chr_results<-foreach(chromsome = 1:22) %dopar% {
    library(dplyr)
    library(arrow)
    library(plink2R)
    library(PDSCE)
    library(doParallel)
    library(purrr)
    library(progress)
    print(chromsome)
    wd2<-paste(c("~/wjk/bib/gwas/gwas_ch",chromsome,".csv"),collapse="")
    gwas.ch<-read.csv(wd2)
    gwas.ch<-subset(gwas.ch,FRQ_A_38691>=0.01 & FRQ_U_186843>=0.01)
    wd3<-paste(c("~/wjk/bib/lookup/lookup_ch",chromsome,".csv"),collapse="")
    lookup.ch<-read.csv(wd3,header = T)
    qtl_files<-list.files(paste0("~/wjk/bib/sqtl/",tissue))
    wd1<-grepl(paste0("chr",chromsome,".par"),qtl_files) %>% qtl_files[.]
    eqtl.ch<-read_parquet(paste0("~/wjk/bib/sqtl/",tissue,"/",wd1))
    eqtl.ch<-as.data.frame(eqtl.ch)
    eqtl.ch<-subset(eqtl.ch,maf >= 0.01 & pval_nominal<0.05)
    eqtl.ch<-eqtl.ch[is.na(eqtl.ch$slope)==F,]
    merged<-left_join(eqtl.ch, lookup.ch, by = "variant_id")
    merged<-merged[!is.na(merged$SNP),]
    smr<-merge(merged,gwas.ch,by="SNP",all=F)
    smr$z1<-smr$slope/smr$slope_se
    smr$z2<-log(smr$OR)/smr$SE
    smr$betay<-log(smr$OR)
    genes<-unique(smr$phenotype_id)
    ld_dir<-"~/wjk/fusion_twas-master/LDREF/1000G.EUR."
    geno<-read_plink(paste0(ld_dir,chromsome),impute="avg")
    refer<-scale(geno$bed)
    common_snp<-intersect(unique(smr$SNP),colnames(refer))
    refer<-refer[,common_snp]
    every_chr_egger<-list()
    every_chr_lda<-list()    
    all_gene_results<-list()
    for (ppp in 1:length(genes)) {
      #ppp<-1
      if (ppp%%100==0){
        print(paste(chromsome,":",ppp,"/",length(genes)))
      }
      gene<-genes[ppp]
      gdata<-smr[smr$phenotype_id==gene,]
      gdata<-gdata[gdata$OR!=1,] #filter out SNP make sure availability of the HEIDI
      snplist<-intersect(gdata$SNP,common_snp)
      
      if (length(snplist)>1 & nrow(gdata)!=0){
        gdata<-gdata[gdata$SNP %in% snplist,]
        refer_sub<-refer[,snplist]
        ld.mat <- t(refer_sub) %*% refer_sub / (nrow(refer_sub)-1)
        gdata<-gdata[order(gdata$pval_nominal,decreasing = F),]
        ld.mat<-ld.mat[gdata$SNP,gdata$SNP]
        aa<-clump(ld.mat,threshold = 0.9)$LD
        if (length(aa)!=1) {
          gdata_lda<-gdata[gdata$SNP%in%rownames(aa),]
          aa<-aa[gdata_lda$SNP,gdata_lda$SNP]
          aa<-0.99*aa+0.01*diag(nrow(aa))
          LDsolvesigma<-solve(aa)
          sigma<-aa
          betaE<-as.vector(LDsolvesigma %*% gdata_lda$slope)
          ##get the conditional estimate for GWAS effect
          betaG<-as.vector(LDsolvesigma %*% log(gdata_lda$OR))
          sigmay<-aa %*% gdata_lda$SE %*% t(gdata_lda$SE)
          solvesigma<-inversesigma(sigmay)
          lda_result<-LDA.MREgger(X=betaE,Y=betaG,W=solvesigma)
          lda_result<-data.frame(phenotype=gdata_lda$phenotype_id[1],
                                 causal_b=lda_result[2,1],
                                 causal_se=lda_result[2,2],
                                 causal_t=lda_result[2,3],
                                 causal_p=lda_result[2,4],
                                 intercept_estimate=lda_result[1,1],
                                 intercept_se=lda_result[1,2],
                                 intercept_t=lda_result[1,3],
                                 intercept_p=lda_result[1,4],
                                 iv_num=nrow(gdata_lda),
                                 chr=chromsome)
        } else {
          lda_result<-data.frame(phenotype=gdata_lda$phenotype_id[1],
                                 causal_b=NA,
                                 causal_se=NA,
                                 causal_t=NA,
                                 causal_p=NA,
                                 intercept_estimate=NA,
                                 intercept_se=NA,
                                 intercept_t=NA,
                                 intercept_p=NA,
                                 iv_num=NA,
                                 chr=NA)
        }
        
        clump.ivw<-clump(ld.mat,threshold = 0.4)
        if (length(clump.ivw$SNP)!=1) {
          #MR-Egger
          gdata_egger<-gdata[gdata$SNP %in% clump.ivw$SNP,]
          mr_egger<-lm(log(gdata_egger$OR)~gdata_egger$slope,weight=(gdata_egger$SE)^-2) %>% summary(.)
          mr_egger_result<-mr_egger$coefficients
          mr_egger_result<-data.frame(phenotype=gdata_egger$phenotype_id[1],
                                      causal_b=mr_egger_result[2,1],
                                      causal_se=mr_egger_result[2,2],
                                      causal_t=mr_egger_result[2,3],
                                      causal_p=mr_egger_result[2,4],
                                      intercept_estimate=mr_egger_result[1,1],
                                      intercept_se=mr_egger_result[1,2],
                                      intercept_t=mr_egger_result[1,3],
                                      intercept_p=mr_egger_result[1,4],
                                      iv_num=nrow(gdata_egger),
                                      rsquare=mr_egger$r.squared,
                                      chr=chromsome)
        } else {
          mr_egger_result<-data.frame(phenotype=gdata_lda$phenotype_id[1],
                                 causal_b=NA,
                                 causal_se=NA,
                                 causal_t=NA,
                                 causal_p=NA,
                                 intercept_estimate=NA,
                                 intercept_se=NA,
                                 intercept_t=NA,
                                 intercept_p=NA,
                                 iv_num=NA,
                                 rsquare=NA,
                                 chr=NA)
        }
        
        result_for_gene<-list(lda_result=lda_result,
                              mr_egger_result=mr_egger_result)
      } else {
        result_for_gene<-NULL
      }
      all_gene_results[[ppp]]<-result_for_gene
    }
    
    #filter out the null genes
    all_gene_results<-discard(all_gene_results, is.null)
    all_gene_result_lda<-lapply(all_gene_results,function(x){
      return(x$lda_result)
    }) %>% do.call(rbind,.)
    all_gene_result_lda<-all_gene_result_lda[!is.na(all_gene_result_lda$causal_b),]
    
    all_gene_result_egger<-lapply(all_gene_results,function(x){
      return(x$mr_egger_result)
    }) %>% do.call(rbind,.)
    all_gene_result_egger<-all_gene_result_egger[!is.na(all_gene_result_egger$causal_b),]
    
    return_list<-list(all_gene_result_lda=all_gene_result_lda,
                      all_gene_result_egger=all_gene_result_egger)
    
    return(return_list)
  }
  
  all_chr_result_gsmr<-lapply(all_chr_results,function(x){
    return(x$all_gene_result_gsmr)
  }) %>% do.call(rbind,.)
  
  all_chr_result_lda<-lapply(all_chr_results,function(x){
    return(x$all_gene_result_lda)
  }) %>% do.call(rbind,.)
  
  all_chr_result_egger<-lapply(all_chr_results,function(x){
    return(x$all_gene_result_egger)
  }) %>% do.call(rbind,.)
  
  
  out_file<-"/output_directory"
  
  write.csv(all_chr_result_egger,
            paste0(out_file,"stwas_mr_egger/",tissue,"_stwas_mr_egger.csv"),row.names = F)
  
  write.csv(all_chr_result_lda,
            paste0(out_file,"stwas_lda/",tissue,"_stwas_lda.csv"),row.names = F)
  
  write.csv(all_chr_result_gsmr,
            paste0(out_file,"stwas_gsmr/",tissue,"_stwas_gsmr.csv"),row.names = F)
}



stopCluster(cl) 

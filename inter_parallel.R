rm(list = ls())
library(arrow)
library(PDSCE)
library(plink2R)
library(doParallel)
library(purrr)
library(dplyr)
source("/bib_function.R")

cl <- makeCluster(6) 
registerDoParallel(cl)
ivw_list<-list()
mr_egger_list<-list()
lda_list<-list()
for (m in 1:22) {
  print(paste("chromosome:",m))
  #loading brainMeta mQTL summary data, these data have been quality controlled based on 0.05 of P-value
  mqtl.ch<-read.csv(paste0("/mqtl_qc-2_ch",m,".csv"))
  mqtl.ch<-mqtl.ch[grepl("cg", mqtl.ch$Probe),]
  location<-data.frame(Probe=mqtl.ch$Probe,Probe_bp=mqtl.ch$Probe_bp)
  probe<-location[!duplicated(location),]
  #loading brainMeta eQTL summary data or sQTL summary data, the data are divided into 22 data based on chromosome
  wd<-paste0("/eqtl_meta_chr",m,".txt")
  eqtl.ch<-read.table(wd,sep = "\t",header = T)
  location<-data.frame(Probe=eqtl.ch$Probe,Probe_bp=eqtl.ch$Probe_bp)
  probe.e<-location[!duplicated(location),]
  probe.e$down<-probe.e$Probe_bp-2000000
  probe.e$up<-probe.e$Probe_bp+2000000
  #loading 1000 Genome data
  ld_dir<-"/1000G.EUR."
  geno<-read_plink(paste0(ld_dir,m),impute="avg")
  refer<-scale(geno$bed)
  result_per_m<-foreach(ppp = 1:nrow(probe.e)) %dopar% {
    source("/bib_function.R")
    library(PDSCE)
    library(purrr)
    library(dplyr)
    print(ppp)
    cpg.vec<-probe[probe.e$down[ppp]<=probe$Probe_bp & probe.e$up[ppp]>=probe$Probe_bp,1]
    part1<-mqtl.ch[mqtl.ch$Probe %in% cpg.vec,]
    part2<-eqtl.ch[eqtl.ch$Probe==probe.e$Probe[ppp],]
    smr<-merge(part1,part2,by="SNP",all=F)
    p_threshold<-0.05/length(unique(smr$Probe.x))
    common_snp<-intersect(unique(smr$SNP),colnames(refer))
    smr<-smr[smr$SNP %in% common_snp,]
    iv_num<-table(smr$Probe.x)
    smr<-smr[smr$Probe.x %in% names(iv_num[iv_num>2]),]
    if (nrow(smr)>=1){
      smr.list<-split(smr,smr$Probe.x)
      res_per_gene<-lapply(smr.list,run_inter)
      ivw_per_gene<-lapply(res_per_gene, function(x){
        return(x$ivw)
      }) %>% discard(., is.null) %>% do.call(rbind,.)
      
      mr_egger_per_gene<-lapply(res_per_gene, function(x){
        return(x$mr_egger)
      }) %>% discard(., is.null) 
      if (length(mr_egger_per_gene)!=0) {
        mr_egger_per_gene<-mr_egger_per_gene %>% do.call(rbind,.) %>% filter(.,!is.nan(causal_p))
      } else {
        mr_egger_per_gene<-NULL
      }
      
      lda_per_gene<-lapply(res_per_gene, function(x){
        return(x$lda)
      }) %>% discard(., is.null)
      if (length(lda_per_gene)!=0) {
        lda_per_gene<-lda_per_gene %>% do.call(rbind,.) %>% filter(.,!is.nan(causal_p))
      } else {
        lda_per_gene<-NULL
      }
      res_per_chr<-list(ivw=ivw_per_gene,
                        mr_egger=mr_egger_per_gene,
                        lda=lda_per_gene)
      res_per_chr<-lapply(res_per_chr,function(tt){
        if (!is.null(tt)) {
          if (nrow(tt)>=1) {
            tt$p_threshold<-p_threshold
            rownames(tt)<-NULL
          } else {
            tt$p_threshold<-numeric(0)
          }
        }
        return(tt)
      })
    } else {
      res_per_chr<-NULL
    }
    write.table(ppp,paste0("/out.txt"))
    return(res_per_chr)
  }
  result_per_m<-result_per_m %>% discard(., is.null)
  all_ivw<-lapply(result_per_m,function(dat){
    return(dat$ivw)
  }) %>% discard(., is.null) %>% do.call(rbind,.)
  
  all_mr_egger<-lapply(result_per_m,function(dat){
    return(dat$mr_egger)
  }) %>% discard(., is.null) %>% do.call(rbind,.)
  
  all_lda<-lapply(result_per_m,function(dat){
    return(dat$lda)
  }) %>% discard(., is.null) %>% do.call(rbind,.)
  
  ivw_list[[m]]<-all_ivw
  mr_egger_list[[m]]<-all_mr_egger
  lda_list[[m]]<-all_lda
}

final_ivw<-do.call(rbind,ivw_list)
final_mr_egger<-do.call(rbind,mr_egger_list)
final_lda<-do.call(rbind,lda_list)

write.csv(final_ivw,"/output/m_e/m_e_ivw.csv",row.names = F)
write.csv(final_mr_egger,"/output/m_e/m_e_mr_egger.csv",row.names = F)
write.csv(final_lda,"/output/m_e/m_e_lda.csv",row.names = F)

stopCluster(cl) 



rm(list=ls())
source("/bib_function.R")
for (m in 1:14) {
  threshold<-0.4
  library(arrow)
  library(dplyr)
  library(PMR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(PDSCE)
  library(ieugwasr)
  samplen<-read.csv("samples.csv")
  n1<-samplen$eup[samplen$tis==tissue]
  wd1<-paste(c("/eqtl/",tissue,"/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_EUR_eQTL_all_associations_Pituitary.v8.EUR.allpairs.chr",m,".parquet"),collapse="")
  eqtl.ch<-read_parquet(wd1)
  eqtl.ch<-as.data.frame(eqtl.ch)
  eqtl.ch<-eqtl.ch[is.na(eqtl.ch$slope)==F,-10]
  eqtl.ch<-subset(eqtl.ch,maf >= 0.01 & pval_nominal<0.05)
  gc()
  wd2<-paste(c("/gwas_ch",m,".csv"),collapse="")
  gwas.ch<-read.csv(wd2)
  gwas.ch<-subset(gwas.ch,FRQ_A_38691>=0.01 & FRQ_U_186843>=0.01)
  gc()
  wd3<-paste(c("lookup/lookup_ch",m,".csv"),collapse="")
  lookup.ch<-read.csv(wd3,header = F)
  lookup.ch<-lookup.ch[,c(1,7)]
  colnames(lookup.ch)<-c("variant_id","SNP")
  gc()
  merged <- left_join(eqtl.ch, lookup.ch, by = "variant_id")
  merged<-subset(merged,SNP!=".")
  smr<-merge(merged,gwas.ch,by="SNP",all=F)
  smr$z1<-smr$slope/smr$slope_se
  smr$z2<-log(smr$OR)/smr$SE
  genes<-unique(smr$phenotype_id)
  #set directory to the folder where PLINK and 1000 Genome data locate
  setwd("/plink")
  pmr.tab<-as.data.frame(matrix(NA,nrow = length(genes),ncol = 7))
  for (ppp in 1:length(genes)) {
    gdata<-subset(smr,phenotype_id==genes[ppp])
    gdata<-gdata[order(gdata$pval_nominal,decreasing = F),]
    write.table(gdata$SNP, file = "output.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    cmd1<-paste(c("plink --bfile 1000G.EUR.QC.",m," --extract output.txt --make-bed  --keep-allele-order --out snpmat"),collapse="")
    a1<-try({suppressWarnings(system(cmd1, intern = T))})
    if (a1[15]=="Error: No variants remaining after --extract.") {
      beta<-gdata$z2[1]/gdata$z1[1]
      beta.se<-abs(gdata$SE[1]/gdata$slope[1])
      t<-(gdata$z1[1]*gdata$z2[1])^2/(gdata$z1[1]^2+gdata$z2[1]^2)
      pv<-1-pchisq(t,df=1)
      pmr.result<-list(causal_effect=beta,
                       causal_pvalue=pv,
                       pleiotropy_effect=NA,
                       pleiotropy_pvalue=NA,
                       sigma_cisSNP=NA,
                       sigma_error_1=NA,
                       sigma_error_2=NA)
    } else {
      invisible(a1)
      a2<-system("plink --bfile snpmat   --r square  --keep-allele-order --write-snplist --out LDMatrixTemp", intern = TRUE)
      invisible(a2)
      refer.ld<-as.matrix(read.table("LDMatrixTemp.ld",header = F))
      snplist<-read.table("LDMatrixTemp.snplist",header=F)[,1]
      colnames(refer.ld)<-snplist
      rownames(refer.ld)<-snplist
      gdata<-merge(data.frame(SNP=snplist),gdata,by="SNP",all.x = T)
      gdata<-gdata[order(gdata$pval_nominal,decreasing = F),]
      r.snp<-match(gdata$SNP,snplist)
      ld.mat<-refer.ld[r.snp,r.snp]
      if (length(snplist)==1) {
        pmr.result<-try({suppressWarnings(PMR_summary_Egger(Zscore_1 = gdata$z1[1],
                                                            Zscore_2 = gdata$z2[1],
                                                            Sigma1sin = matrix(1),
                                                            Sigma2sin = matrix(1),
                                                            samplen1 = n1,
                                                            samplen2 = n2,
                                                            lambda=1, max_iterin = 1000, epsin = 1e-05,
                                                            Heritability_geneexpression_threshold = 0))})
        if (pmr.result[1]=="Error in PMR_summary_Egger_CPP(betaxin, betayin, Sigma1sin, Sigma2sin,  : \n  chol(): decomposition failed\n") {
          
          beta<-gdata$z2[1]/gdata$z1[1]
          beta.se<-abs(gdata$SE[1]/gdata$slope[1])
          t<-(gdata$z1[1]*gdata$z2[1])^2/(gdata$z1[1]^2+gdata$z2[1]^2)
          pv<-1-pchisq(t,df=1)
          pmr.result<-list(causal_effect=beta,
                           causal_pvalue=pv,
                           pleiotropy_effect=NA,
                           pleiotropy_pvalue=NA,
                           sigma_cisSNP=NA,
                           sigma_error_1=NA,
                           sigma_error_2=NA)
        }
        
      } else {
        ld.mat0<-ld.mat 
        gdata0<-gdata
        clump.result<-clump(ld.mat,threshold = threshold)
        ld.mat<-clump.result$LD
        top.snp<-clump.result$SNP
        gdata<-gdata[gdata$SNP%in%top.snp,]
        n2<-round(mean(gdata$Nca+gdata$Nco))
        if (length(top.snp)==1) {
          pmr.result<-try({suppressWarnings(PMR_summary_Egger(Zscore_1 = gdata$z1[1],
                                                              Zscore_2 = gdata$z2[1],
                                                              Sigma1sin = matrix(1),
                                                              Sigma2sin = matrix(1),
                                                              samplen1 = n1,
                                                              samplen2 = n2,
                                                              lambda=1, max_iterin = 1000, epsin = 1e-05,
                                                              Heritability_geneexpression_threshold = 0))},silent = T)
          if (pmr.result[1]=="Error in PMR_summary_Egger_CPP(betaxin, betayin, Sigma1sin, Sigma2sin,  : \n  chol(): decomposition failed\n") {
            clump.result0<-clump(ld.mat0,threshold = 0.8)
            ld.mat0<-clump.result0$LD
            top.snp0<-clump.result0$SNP
            gdata0<-gdata0[gdata0$SNP%in%top.snp0,]
            
            se.x<-as.vector(gdata0$slope_se)
            se.y<-as.vector(gdata0$SE)
            omega<-t(t(se.x*ld.mat0)*se.x)
            beta.x<-as.vector(gdata0$slope)
            beta.y<-as.vector(log(gdata0$OR))
            omega<-pdsoft(omega,lam=0)$sigma
            ivw<-(t(beta.x) %*% solve(omega) %*% beta.y)/(t(beta.x) %*% solve(omega) %*% beta.x)
            ivw.se<-sqrt(1/(t(beta.x) %*% solve(omega) %*% beta.x))
            chisq<-(ivw/ivw.se)^2
            pv<-1-pchisq(chisq,df=1)
            
            pmr.result<-list(causal_effect=ivw,
                             causal_pvalue=pv,
                             pleiotropy_effect=NA,
                             pleiotropy_pvalue=NA,
                             sigma_cisSNP=NA,
                             sigma_error_1=NA,
                             sigma_error_2=NA)
          }
          
        } else {
          w=0
          while (1) {
            if (w>200) {
              print(ppp)
              break
            }
            w=w+1
            pmr.result<-try({suppressWarnings(PMR_summary_Egger(Zscore_1 = gdata$z1,
                                                                Zscore_2 = gdata$z2,
                                                                Sigma1sin = as.matrix(ld.mat),
                                                                Sigma2sin = as.matrix(ld.mat),
                                                                samplen1 = n1,
                                                                samplen2 = n2,
                                                                lambda=1, max_iterin = 1000, epsin = 1e-05,
                                                                Heritability_geneexpression_threshold = 0))},silent = T)
            if (pmr.result[1]!="Error in PMR_summary_Egger_CPP(betaxin, betayin, Sigma1sin, Sigma2sin,  : \n  chol(): decomposition failed\n") {
              break
            } else if (pmr.result[1]=="Error in PMR_summary_Egger_CPP(betaxin, betayin, Sigma1sin, Sigma2sin,  : \n  chol(): decomposition failed\n" & length(ld.mat)==1) {
              #gsmr
              
              clump.result0<-clump(ld.mat0,threshold = 0.8)
              ld.mat0<-clump.result0$LD
              top.snp0<-clump.result0$SNP
              gdata0<-gdata0[gdata0$SNP%in%top.snp0,]
              
              se.x<-as.vector(gdata0$slope_se)
              se.y<-as.vector(gdata0$SE)
              omega<-t(t(se.x*ld.mat0)*se.x)
              beta.x<-as.vector(gdata0$slope)
              beta.y<-as.vector(log(gdata0$OR))
              omega<-pdsoft(omega,lam=0)$sigma
              ivw<-(t(beta.x) %*% solve(omega) %*% beta.y)/(t(beta.x) %*% solve(omega) %*% beta.x)
              ivw.se<-sqrt(1/(t(beta.x) %*% solve(omega) %*% beta.x))
              chisq<-(ivw/ivw.se)^2
              pv<-1-pchisq(chisq,df=1)
              
              pmr.result<-list(causal_effect=ivw,
                               causal_pvalue=pv,
                               pleiotropy_effect=NA,
                               pleiotropy_pvalue=NA,
                               sigma_cisSNP=NA,
                               sigma_error_1=NA,
                               sigma_error_2=NA)
              break
              
              
            } else {
              cut.r<-threshold/(2^w)
              clump.result<-clump(ld.mat,threshold = cut.r)
              ld.mat<-clump.result$LD
              top.snp<-clump.result$SNP
              gdata<-gdata[gdata$SNP%in%top.snp,]
              n2<-round(mean(gdata$Nca+gdata$Nco))
            }
            
          }
          
        }
      }
    }
    pmr.tab[ppp,]<-unlist(pmr.result)
    print(paste("chromosome",m,":",ppp,"/",length(genes)))
  }
  colnames(pmr.tab)<-names(unlist(pmr.result))
  pmr.tab$ENSEMBL<-sub("\\..*", "", genes)
  suppressMessages(suppressWarnings(gene.id <- bitr(sub("\\..*", "", genes), 
                                                    fromType = "ENSEMBL",
                                                    toType = "SYMBOL",
                                                    OrgDb = org.Hs.eg.db)))

  output<-merge(pmr.tab,gene.id,by="ENSEMBL",all.x=T)
  wd4<-paste(c("output/",tissue,"/",tissue,"PMR-Egger_ch",m,".csv"),collapse="")
  write.csv(output,wd4)
  if (m!=22) {
    vars <- ls()
    vars <- vars[vars != "m"] 
    rm(list = vars)
    gc()
  }
}



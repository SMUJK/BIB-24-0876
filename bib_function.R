clump<-function(mat,threshold) {
  input.snp<-colnames(mat)
  mat<-matrix(as.numeric(mat),nrow(mat),nrow(mat))
  colnames(mat)<-input.snp
  rownames(mat)<-input.snp
  ld.mat<-abs(mat)
  top.snp<-c()
  i<-1
  while (1) {
    ld.vec<-ld.mat[1,]
    top.snp[i]<-colnames(ld.mat)[1]
    i=i+1
    left<-ld.vec<threshold
    sum(left)
    if (sum(left)==0) {
      break
    } else if (sum(left)==1) {
      top.snp[i]<-colnames(ld.mat)[left]
      break
    }
    ld.mat<-ld.mat[left,left]
  }
  pos<-colnames(mat)%in%top.snp
  result.mat<-mat[pos,pos]
  ending<-list(LD=result.mat,SNP=top.snp)
  return(ending)
}



LDA.MREgger<-function(X,Y,W){
  bX<-cbind(1,X)
  first<-crossprod(bX,W)%*%bX
  res <- tryCatch({
    bread<-solve(first)
  }, error = function(e) {
    e 
  })
  if (!is.null(res)) {
    bread<-solve(pdsoft(first,lam = 0)$sigma)
  }
  theEsts<-bread%*%crossprod(bX,W%*%Y)
  theresid<-c(Y-theEsts[1]-X*theEsts[2])
  Sig.Est<-c(crossprod(theresid,W%*%theresid))/(length(X)-2)
  finresults<- cbind(theEsts,diag(bread)*Sig.Est)
  TestStat<-theEsts/sqrt(finresults[,2])
  Pvals<-2*pt(abs(TestStat),df = nrow(bX)-2,lower.tail = F)
  return(cbind(finresults,TestStat,Pvals))
}



inversesigma<-function(sigmaG,sigma){
  sigmaG<-as.matrix(sigmaG)
  tmp=sigma*(sigmaG%*%t(sigmaG))
  tmp<-0.99*tmp+0.01*diag(nrow(tmp))
  solvetmp<-solve(tmp)
  re<-sigma%*%solvetmp%*%sigma
  return(re)
}



HEIDI<-function(pos.smr,ld.mat) {
  rin2<-match(rownames(ld.mat),pos.smr$SNP)
  heidi.data<-pos.smr[rin2,]
  heidi.data$causal.se<-sqrt(((heidi.data$slope_se/heidi.data$slope)^2+(heidi.data$SE/log(heidi.data$OR))^2)*((heidi.data$causal_b)^2))
  a<-diag(heidi.data$causal.se)
  mat0<-a %*% ld.mat %*% a
  top<-which.min(heidi.data$pval_nominal)
  if (length(top)>1) {
    top<-top[1]
  }
  mat1<-mat0[-top,-top]
  mat2<-matrix(mat0[-top,top],nrow=(nrow(mat0)-1),ncol=(nrow(mat0)-1))
  mat3<-t(mat2)
  mat4<-matrix(mat0[top,top],nrow=(nrow(mat0)-1),ncol=(nrow(mat0)-1))
  v<-mat1-mat2-mat3+mat4
  d<-heidi.data$causal_b[-top]-heidi.data$causal_b[top]
  chi<-t(d) %*% solve(v) %*% (d)
  heidi.p<-1-pchisq(as.numeric(chi),df=(nrow(mat0)-1))
  return(heidi.p)
}



run_inter<-function(gdata) {
  snp_vec<-gdata$SNP
  rmat<-refer[,snp_vec]
  ld_mat <- t(rmat) %*% rmat / (nrow(rmat)-1)
  clump1<-clump(mat = ld_mat,threshold = 0.4)
  if (length(clump1$SNP)>1) {
    input1<-gdata[gdata$SNP %in% clump1$SNP,]
    ivw_res<-ivw_fun(input1)
    mr_egger_res<-mr_egger_fun(input1)
  } else {
    ivw_res<-NULL
    mr_egger_res<-NULL
  }
  
  clump2<-clump(mat = ld_mat,threshold = 0.9)
  if (length(clump2$SNP)>1) {
    lda_res<-lda_fun(input=gdata,LD=clump2$LD)
  } else {
    lda_res<-NULL
  }
  fun_res<-list(ivw=ivw_res,mr_egger=mr_egger_res,lda=lda_res)
  return(fun_res)
}





ivw_fun<-function(input){
  b<-sum(input$b.x*input$b.y/((input$SE.y)^2))/sum((input$b.x/input$SE.y)^2)
  se<-sqrt(1/sum((input$b.x/input$SE.y)^2))
  p<-1-pchisq(((b/se)^2),df=1)
  ivw_result<-data.frame(exposure=input$Probe.x[1],
                         outcome=input$Probe.y[1],
                         symbol=input$Gene.y[1],
                         chromosome=input$Chr.y[1],
                         causal_b=b,
                         causal_se=se,
                         causal_p=p)
  
  return(ivw_result)
}

mr_egger_fun<-function(input){
  #input<-gdata
  if (var(input$b.x)==0) {
    input$b.x[1]<-input$b.x[1]+0.00001
  }
  mr_egger<-lm(input$b.y ~ input$b.x, weight=(input$SE.y)^-2) %>% summary(.)
  mr_egger_result<-mr_egger$coefficients
  mr_egger_result<-data.frame(exposure=input$Probe.x[1],
                              outcome=input$Probe.y[1],
                              symbol=input$Gene.y[1],
                              chromosome=input$Chr.y[1],
                              causal_b=mr_egger_result[2,1],
                              causal_se=mr_egger_result[2,2],
                              causal_t=mr_egger_result[2,3],
                              causal_p=mr_egger_result[2,4],
                              intercept_b=mr_egger_result[1,1],
                              intercept_se=mr_egger_result[1,2],
                              intercept_t=mr_egger_result[1,3],
                              intercept_p=mr_egger_result[1,4],
                              iv_num=nrow(input),
                              rsquare=mr_egger$r.squared)
  return(mr_egger_result)
}

lda_fun<-function(input,LD) {
  LD<-as.matrix(LD)
  gdata_lda<-input[input$SNP %in% rownames(LD),]
  if (var(gdata_lda$b.x)==0) {
    gdata_lda$b.x[1]<-gdata_lda$b.x[1]+0.00001
  }
  LD<-LD[gdata_lda$SNP,gdata_lda$SNP]
  LD<-0.99*LD+0.01*diag(nrow(LD))
  LDsolvesigma<-solve(LD)
  sigma<-LD
  betaE<-as.vector(LDsolvesigma %*% gdata_lda$b.x)
  ##get the conditional estimate for GWAS effect
  betaG<-as.vector(LDsolvesigma %*% gdata_lda$b.y)
  sigmay<-LD %*% gdata_lda$SE.y %*% t(gdata_lda$SE.y)
  solvesigma<-inversesigma(sigmaG=sigmay,sigma = LD)
  lda_result<-LDA.MREgger(X=betaE,Y=betaG,W=solvesigma)
  lda_result<-data.frame(exposure=input$Probe.x[1],
                         outcome=input$Probe.y[1],
                         symbol=input$Gene.y[1],
                         chromosome=input$Chr.y[1],
                         causal_b=lda_result[2,1],
                         causal_se=lda_result[2,2],
                         causal_t=lda_result[2,3],
                         causal_p=lda_result[2,4],
                         intercept_estimate=lda_result[1,1],
                         intercept_se=lda_result[1,2],
                         intercept_t=lda_result[1,3],
                         intercept_p=lda_result[1,4],
                         iv_num=nrow(gdata_lda))
  return(lda_result)
}



run_ivw<-function(b.x,b.y,SE.y,exposure,chr){
  b<-sum(b.x*b.y/((SE.y)^2))/sum((b.x/SE.y)^2)
  se<-sqrt(1/sum((b.x/SE.y)^2))
  p<-1-pchisq(((b/se)^2),df=1)
  ivw_result<-data.frame(exposure=exposure,
                         chromosome=chr,
                         causal_b=b,
                         causal_se=se,
                         causal_p=p)
  
  return(ivw_result)
}


run_mr_egger<-function(b.x,b.y,SE.y,exposure,chr){
  #input<-gdata
  if (var(b.x)==0) {
    b.x[1]<-b.x[1]+0.00001
  }
  mr_egger<-lm(b.y ~ b.x, weight=(SE.y)^-2) %>% summary(.)
  mr_egger_result<-mr_egger$coefficients
  mr_egger_result<-data.frame(exposure=exposure,
                              chromosome=chr,
                              causal_b=mr_egger_result[2,1],
                              causal_se=mr_egger_result[2,2],
                              causal_t=mr_egger_result[2,3],
                              causal_p=mr_egger_result[2,4],
                              intercept_b=mr_egger_result[1,1],
                              intercept_se=mr_egger_result[1,2],
                              intercept_t=mr_egger_result[1,3],
                              intercept_p=mr_egger_result[1,4],
                              iv_num=length(b.x),
                              rsquare=mr_egger$r.squared)
  return(mr_egger_result)
}


run_lda<-function(b.x,b.y,SE.y,exposure,chr,LD) {
  LD<-as.matrix(LD)
  if (var(b.x)==0) {
    b.x[1]<-b.x[1]+0.00001
  }
  LD<-0.99*LD+0.01*diag(nrow(LD))
  LDsolvesigma<-solve(LD)
  sigma<-LD
  betaE<-as.vector(LDsolvesigma %*% b.x)
  ##get the conditional estimate for GWAS effect
  betaG<-as.vector(LDsolvesigma %*% b.y)
  sigmay<-LD %*% SE.y %*% t(SE.y)
  solvesigma<-inversesigma(sigmaG=sigmay,sigma = LD)
  lda_result<-LDA.MREgger(X=betaE,Y=betaG,W=solvesigma)
  lda_result<-data.frame(exposure=exposure,
                         chromosome=chr,
                         causal_b=lda_result[2,1],
                         causal_se=lda_result[2,2],
                         causal_t=lda_result[2,3],
                         causal_p=lda_result[2,4],
                         intercept_estimate=lda_result[1,1],
                         intercept_se=lda_result[1,2],
                         intercept_t=lda_result[1,3],
                         intercept_p=lda_result[1,4],
                         iv_num=length(b.x))
  return(lda_result)
}



run_3_model<-function(gdata) {
  common_snp<-intersect(gdata$SNP,colnames(refer))
  if (length(common_snp)>1) {
    gdata1<-gdata[gdata$SNP %in% common_snp,]
    snp_vec<-gdata1$SNP
    rmat<-refer[,snp_vec]
    ld_mat <- t(rmat) %*% rmat / (nrow(rmat)-1)
    #mr_egger and ivw
    clump1<-clump(mat = ld_mat,threshold = 0.4)
    if (length(clump1$SNP)>2) {
      input1<-gdata1[gdata1$SNP %in% clump1$SNP,]
      ivw_res<-run_ivw(b.x = input1$b, b.y = input1$b.y, 
                       SE.y = input1$SE.y,exposure = input1$Probe[1],
                       chr = input1$Chr[1])
      mr_egger_res<-run_mr_egger(b.x = input1$b, b.y = input1$b.y, 
                                 SE.y = input1$SE.y,exposure = input1$Probe[1],
                                 chr = input1$Chr[1])
    } else if (length(clump1$SNP)==2) {
      input1<-gdata1[gdata1$SNP %in% clump1$SNP,]
      ivw_res<-run_ivw(b.x = input1$b, b.y = input1$b.y, 
                       SE.y = input1$SE.y,exposure = input1$Probe[1],
                       chr = input1$Chr[1])
      mr_egger_res<-NULL
    } else {
      ivw_res<-NULL
      mr_egger_res<-NULL
    }
    
    #lda
    clump2<-clump(mat = ld_mat,threshold = 0.9)
    if (length(clump2$SNP)>2) {
      input2<-gdata1[gdata1$SNP %in% clump2$SNP,]
      ld<-clump2$LD[input2$SNP,input2$SNP]
      lda_res<-run_lda(b.x=input2$b, b.y=input2$b.y,
                       SE.y=input2$SE.y,exposure=input2$Probe[1],
                       chr = input2$Chr[1],LD=ld)
    } else {
      lda_res<-NULL
    }
    
  } else {
    ivw_res<-NULL
    mr_egger_res<-NULL
    lda_res<-NULL
  }
  return(list(ivw_res=ivw_res,mr_egger_res=mr_egger_res,lda_res=lda_res))
}

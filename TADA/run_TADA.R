#!/histor/sun/wangtao/3_software/anaconda3/envs/linlin_R_4.2/bin/Rscript --vanilla
#PBS -N deeptool
#PBS -l nodes=1:ppn=6
#PBS -j oe
#PBS -q silver

library(openxlsx);library(stringr);library(purrr);library(dplyr)

########## Calculate the BF for SNVs/indels ###################

## General formula

## log Bayes factor for de novo mutations
log.bayes.factor.dn <- function(x.dn, n.dn, mu, gamma.dn, beta.dn){
  marg.lik0 <- dpois(x.dn, 2*n.dn*mu, log=TRUE)
  marg.lik1 <- dnbinom(x.dn, gamma.dn*beta.dn, beta.dn/(beta.dn+2*n.dn*mu), log=TRUE)
  lBF <- marg.lik1 - marg.lik0
  return(lBF = lBF)
}

## log Bayes factor for case–control / inherited variants
bayes.factor.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1, rho0, nu0){
  marglik0.cc <- evidence.null.cc(x.cc, n.cc, rho0, nu0)
  marglik1.cc <- evidence.alt.cc(x.cc, n.cc, gamma.cc, rho1, nu1)
  BF.cn <- marglik1.cc$cn / marglik0.cc$cn
  BF.ca <- marglik1.cc$ca / marglik0.cc$ca
  BF <- BF.cn * BF.ca
  return(BF = BF)
}

## Helper function for bayes.factor.cc — null model
evidence.null.cc <- function(x.cc, n.cc, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x.cc$cn, rho0, nu0/(nu0+n.cc$cn)))
  marglik0.case.log <- log(dnbinom(x.cc$ca, rho0+x.cc$cn, (nu0+n.cc$cn)/(nu0+n.cc$ca+n.cc$cn)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  return(list(cn = exp(marglik0.ctrl.log), ca = exp(marglik0.case.log)))
}

## Helper function for bayes.factor.cc — alternative model
evidence.alt.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1){
  marglik1.ctrl <- dnbinom(x.cc$cn, rho1, nu1/(nu1+n.cc$cn))
  marglik1.case <- dnbinom(x.cc$ca, rho1+x.cc$cn, (nu1+n.cc$cn)/(nu1+n.cc$ca*gamma.cc+n.cc$cn))
  marglik1 <- marglik1.ctrl * marglik1.case
  return(list(cn = marglik1.ctrl, ca = marglik1.case))
}

#### Load input data (1 metadata + 3 SNV/indel datasets)
TADA_sup_info <- read.csv('./TADA_sup_info_mut_hgnc.csv',header=TRUE)
TADA_sup_info1 <- TADA_sup_info[,!names(TADA_sup_info) %in% c("gene.y")]
TADA_sup_info_for_annotation <- TADA_sup_info1[c("hgnc_symbol","ensembl_gene_id")]

#### Filter SNVs/indels located on SRGs 
gene_sensory <- read.csv('./genecards_sensory_hgnc_supple.csv',header = TRUE)
sensory_gene <- list(gene_sensory$Gene_hgnc)
sensory_gene_vector <- as.vector(unlist(sensory_gene))

## Load annotated SNVs/indels
dn_SNV <- read.csv("./dn_SNV_info_hgnc.csv",header = TRUE)
in_SNV <- read.csv("./in_SNV_info_hgnc.csv",header = TRUE)
cc_SNV <- read.csv("./cc_SNV_info_hgnc.csv",header = TRUE)

dn_SNV_SRG <- dn_SNV[dn_SNV$hgnc_symbol %in% sensory_gene_vector,]
in_SNV_SRG1 <- in_SNV[in_SNV$hgnc_symbol %in% sensory_gene_vector,]
cc_SNV_SRG1 <- cc_SNV[cc_SNV$hgnc_symbol %in% sensory_gene_vector,]
head(dn_SNV)
print(nrow(cc_SNV_SRG1))

## Remove NA rows and convert columns to numeric
in_SNV_SRG <- na.omit(in_SNV_SRG1)
in_SNV_SRG$in.t.ptv <- as.numeric(in_SNV_SRG$in.t.ptv)
in_SNV_SRG$in.u.ptv <- as.numeric(in_SNV_SRG$in.u.ptv)
in_SNV_SRG$in.t.misb <- as.numeric(in_SNV_SRG$in.t.misb)
in_SNV_SRG$in.u.misb <- as.numeric(in_SNV_SRG$in.u.misb)
in_SNV_SRG$in.t.misa <- as.numeric(in_SNV_SRG$in.t.misa)
in_SNV_SRG$in.u.misa <- as.numeric(in_SNV_SRG$in.u.misa)

cc_SNV_SRG <- na.omit(cc_SNV_SRG1)
cc_SNV_SRG$case.ptv <- as.numeric(cc_SNV_SRG$case.ptv)
cc_SNV_SRG$control.ptv <- as.numeric(cc_SNV_SRG$control.ptv)
cc_SNV_SRG$case.misb <- as.numeric(cc_SNV_SRG$case.misb)
cc_SNV_SRG$control.misb <- as.numeric(cc_SNV_SRG$control.misb)
cc_SNV_SRG$case.misa <- as.numeric(cc_SNV_SRG$case.misa)
cc_SNV_SRG$control.misa <- as.numeric(cc_SNV_SRG$control.misa)

## Define BF calculation function for de novo SNVs
BF_DN_SNV <- function(count_case, count_con, n_case, n_con, mut, gamma.dn, beta.dn=.2){
  if(n_con > 0){
    BF_con <- exp(sapply(1:length(mut), function(i)
      log.bayes.factor.dn(x.dn = count_con[i], n.dn = n_con, mu = mut[i],
                          beta.dn = beta.dn, gamma.dn = gamma.dn[i])))
  } else {
    BF_con <- rep(1, length(mut))
  }
  BF_con[BF_con<1] <- 1
  
  BF_case <- exp(sapply(1:length(mut), function(i)
    log.bayes.factor.dn(x.dn = count_case[i], n.dn = n_case, mu = mut[i],
                        beta.dn = beta.dn, gamma.dn = gamma.dn[i])))
  BF_case[BF_case<1] <- 1
  BF <- pmax(BF_case / BF_con, 1)
  BF[which(count_case == 0 & count_con == 0)] <- 1
  BF[is.na(BF)] <- 1
  return(BF)
}

beta.dn <- 0.2
dn_n_case <- 18440
dn_n_con <- 5714

## Compute BF for de novo SNVs (ptv, misb, misa)
BF_dn_ptv <- BF_DN_SNV(
  count_case = dn_SNV_SRG$dn.ptv,
  count_con = dn_SNV_SRG$dn.ptv.sib,
  dn_n_case, dn_n_con,
  mut = dn_SNV_SRG$mu_lof,
  gamma.dn = dn_SNV_SRG$prior.dn.ptv,
  beta.dn = beta.dn)

BF_dn_misb <- BF_DN_SNV(
  count_case = dn_SNV_SRG$dn.misb,
  count_con = dn_SNV_SRG$dn.misb.sib,
  dn_n_case, dn_n_con,
  mut = dn_SNV_SRG$mu_misb,
  gamma.dn = dn_SNV_SRG$prior.dn.misb,
  beta.dn = beta.dn)

BF_dn_misa <- BF_DN_SNV(
  count_case = dn_SNV_SRG$dn.misa,
  count_con = dn_SNV_SRG$dn.misa.sib,
  dn_n_case, dn_n_con,
  mut = dn_SNV_SRG$mu_misa,
  gamma.dn = dn_SNV_SRG$prior.dn.misa,
  beta.dn = beta.dn)

dn_SNV_SRG$BF_dn_ptv <- BF_dn_ptv
dn_SNV_SRG$BF_dn_misb <- BF_dn_misb 
dn_SNV_SRG$BF_dn_misa <- BF_dn_misa
dn_SNV_SRG_selected <- dn_SNV_SRG[c("hgnc_symbol","BF_dn_ptv","BF_dn_misb","BF_dn_misa")]

## Define BF calculation function for inherited / case-control SNVs
BF_CC_SNV <- function(count_case, count_con, n_case, n_con, mut, gamma.cc, nu=5000){
  rho.in = nu * sum(count_con, na.rm=T) / (2*length(mut)*n_con)
  rho.in = as.numeric(rho.in) * mut / mean(as.numeric(gsub(",","",mut)))
  BF = sapply((1:length(mut)), function(i){
    bayes.factor.cc(
      x.cc = data.frame(ca=count_case[i], cn=count_con[i]),
      n.cc = data.frame(ca=n_case, cn=n_con),
      gamma.cc = gamma.cc[i],
      rho1 = rho.in[i], nu1 = nu, rho0 = rho.in[i], nu0 = nu)
  })
  BF[mut == 0] <- 1
  BF[which((count_case + count_con)==0)] <- 1
  return(BF)
}

## Inherited SNVs BF
in_n_case <- 15036
in_n_control <- 5492

BF_in_ptv <- BF_CC_SNV(in_SNV_SRG$in.t.ptv, in_SNV_SRG$in.u.ptv,
                       in_n_case, in_n_control,
                       mut = in_SNV_SRG$mu_lof,
                       gamma.cc = in_SNV_SRG$prior.in.ptv)

BF_in_misb <- BF_CC_SNV(in_SNV_SRG$in.t.misb, in_SNV_SRG$in.u.misb,
                        in_n_case, in_n_control,
                        mut = in_SNV_SRG$mu_misb,
                        gamma.cc = in_SNV_SRG$prior.in.misb)

BF_in_misa <- BF_CC_SNV(in_SNV_SRG$in.t.misa, in_SNV_SRG$in.u.misa,
                        in_n_case, in_n_control,
                        mut = in_SNV_SRG$mu_misa,
                        gamma.cc = in_SNV_SRG$prior.in.misa)

in_SNV_SRG$BF_in_ptv <- BF_in_ptv
in_SNV_SRG$BF_in_misb <- BF_in_misb
in_SNV_SRG$BF_in_misa <- BF_in_misa
in_SNV_SRG_selected <- in_SNV_SRG[c("hgnc_symbol","BF_in_ptv","BF_in_misb","BF_in_misa")]

## Case–control SNV BF
cc_n_case <- 5591
cc_n_control <- 8597

BF_cc_ptv <- BF_CC_SNV(cc_SNV_SRG$case.ptv, cc_SNV_SRG$control.ptv,
                       cc_n_case, cc_n_control,
                       mut = cc_SNV_SRG$mu_lof,
                       gamma.cc = cc_SNV_SRG$prior.cc.ptv)

BF_cc_misb <- BF_CC_SNV(cc_SNV_SRG$case.misb, cc_SNV_SRG$control.misb,
                        cc_n_case, cc_n_control,
                        mut = cc_SNV_SRG$mu_misb,
                        gamma.cc = cc_SNV_SRG$prior.in.misb)

BF_cc_misa <- BF_CC_SNV(cc_SNV_SRG$case.misa, cc_SNV_SRG$control.misa,
                        cc_n_case, cc_n_control,
                        mut = cc_SNV_SRG$mu_misa,
                        gamma.cc = cc_SNV_SRG$prior.in.misa)

cc_SNV_SRG$BF_cc_ptv <- BF_cc_ptv
cc_SNV_SRG$BF_cc_misb <- BF_cc_misb
cc_SNV_SRG$BF_cc_misa <- BF_cc_misa
cc_SNV_SRG_selected <- cc_SNV_SRG[c("hgnc_symbol","BF_cc_ptv","BF_cc_misb","BF_cc_misa")]

## Merge in & cc results
SNV_SRG_in_cc <- merge(x=in_SNV_SRG_selected,y=cc_SNV_SRG_selected,by="hgnc_symbol",all.x=TRUE,all.y=TRUE)

## Merge in + cc + dn
SNV_SRG_in_cc_dn <- merge(x=SNV_SRG_in_cc,y=dn_SNV_SRG_selected,by="hgnc_symbol",all.x=TRUE,all.y=TRUE)
SNV_SRG_in_cc_dn[is.na(SNV_SRG_in_cc_dn)] <- 1.0

head(SNV_SRG_in_cc_dn)
write.csv(SNV_SRG_in_cc_dn,file="./SNV_SRG_BF.csv",row.names = FALSE)

################################# CNV Section ###########################
######################## de novo CNV ###################################

library(GenomicRanges);library(openxlsx);library(stringr)

#### Set parameters
loeuf_threshold <- 0.6
cnv_size_range <- 1
beta.dn <- 0.2

## Number of cases for de novo CNV
dn_CNV_case <- 13786
dn_CNV_con <- 5098

## Load CNV table
CNV <- read.csv("./CNV_merge_note_20241023.csv",header = TRUE)
CNV_SRG <- CNV[CNV$hgnc_genes %in% sensory_gene_vector,]
CNV_SRG_info <- merge(CNV_SRG,TADA_sup_info1,by.x="hgnc_genes",by.y="hgnc_symbol",all.x=TRUE)

## Remove rows with missing gene info
CNV_SRG_info_non_na <- CNV_SRG_info[!is.na(CNV_SRG_info$gene.x), ]

## Filter CNVs with LOEUF ≥ 0.6
CNV_SRG_info_non_na1 <- CNV_SRG_info_non_na[which(CNV_SRG_info_non_na$LOEUF >= 0.6),]
CNV_SRG_info_non_na1$num_genes <- 1

## Extract de novo CNVs
CNV_SRG_dn <- CNV_SRG_info_non_na1[which(CNV_SRG_info_non_na1$Modified_Inheritence == "de novo"),]

### Compute del_dup adjustment factor
del_dup_adj <- table(CNV_SRG_dn[,c("case.control","call")])
del_dup_adj <- (del_dup_adj[2,1]/del_dup_adj[1,1])/(del_dup_adj[2,2]/del_dup_adj[1,2])

## CNV BF for de novo CNVs
BF_DN_CNV <- function(cnv_size_range, cnv_use, info, prior, beta.dn=0.2, n, del_dup_adj=1){
  raw_lBF <- lapply(cnv_size_range, function(i, cnv_use, n, del_dup_adj, prior){
    cnv_sub <- cnv_use[cnv_use$num_genes==i, ]
    if(dim(cnv_sub)[1]>0){
      ## Compute relative risk based on prior
      cnv_sub$relative_risk <- (prior[match(cnv_sub$gene_symbol,info$gene_gencodeV33)-1])/(del_dup_adj+1)
      cnv_sub$relative_risk[is.na(cnv_sub$relative_risk)] = 1
      genes_combo <- table(cnv_sub$gene_symbol)
      marg.lik0 <- dpois(genes_combo, n*cnv_sub$mut[match(names(genes_combo), cnv_sub$gene_symbol,nomatch=1)], log=TRUE)
      marg.lik1 <- dnbinom(genes_combo,
                           cnv_sub$relative_risk[match(names(genes_combo), cnv_sub$gene_symbol,nomatch=1)]*beta.dn,
                           beta.dn/(beta.dn+n*cnv_sub$mut[match(names(genes_combo), cnv_sub$gene_symbol,nomatch=1)]),
                           log=TRUE)
      lBF <- marg.lik1 - marg.lik0
    }
  }, cnv_use, n=n, del_dup_adj=del_dup_adj, prior=prior)
  raw_lBF <- as.data.frame(raw_lBF)
  sum_lBF <- by(raw_lBF$Freq, raw_lBF$Var1, sum)
  BF <- rep(1, length(info$gene_gencodeV33))
  BF[match(names(sum_lBF), info$gene_gencodeV33,nomatch=1)] <- exp(as.numeric(sum_lBF))
  return(BF)
}

############# de novo CNV — deletion
## Extract deletions
del_use <- CNV_SRG_dn[CNV_SRG_dn$call=="DEL",]

## Predict mutation rate for deletion CNV
mut.pred.del <- sapply(cnv_size_range, function(x){
  max(1, sum(del_use$num_genes==x & del_use$case.control=="control")) /
    dn_CNV_con / (length(TADA_sup_info1$prior.dn.ptv) - x)
})
del_use$mut <- mut.pred.del[del_use$num_genes]

## Compute BF for deletion CNVs
BF_dn_del_prob <- BF_DN_CNV(cnv_size_range=cnv_size_range,
                            cnv_use=del_use[del_use$case.control=="case",],
                            info=TADA_sup_info1, prior=TADA_sup_info1$prior.dn.ptv,
                            beta.dn=beta.dn, n=dn_CNV_case)
BF_dn_del_sib <- BF_DN_CNV(cnv_size_range=cnv_size_range,
                           cnv_use=del_use[del_use$case.control=="control",],
                           info=TADA_sup_info1, prior=TADA_sup_info1$prior.dn.ptv,
                           beta.dn=beta.dn, n=dn_CNV_con)
BF_dn_del <- pmax(1, BF_dn_del_prob/BF_dn_del_sib)
keep(BF_dn_del, ~ .x >100000)
BF_dn_del[BF_dn_del>100]
TADA_sup_info1$BF_dn_del <- BF_dn_del

############ de novo CNV — duplication
dup_use <- CNV_SRG_dn[CNV_SRG_dn$call=="DUP",]
mut.pred.dup <- sapply(cnv_size_range, function(x){
  max(1, sum(dup_use$num_genes==x & dup_use$case.control=="control")) /
    dn_CNV_con / (length(TADA_sup_info1$prior.dn.ptv) - x)
})
dup_use$mut <- mut.pred.dup[dup_use$num_genes]

## Compute BF for duplication CNVs
BF_dn_dup_prob <- BF_DN_CNV(cnv_size_range=cnv_size_range,
                            cnv_use=dup_use[dup_use$case.control=="case",],
                            info=TADA_sup_info1, prior=TADA_sup_info1$prior.dn.ptv,
                            beta.dn=beta.dn, n=dn_CNV_case,del_dup_adj=del_dup_adj)
BF_dn_dup_sib <- BF_DN_CNV(cnv_size_range=cnv_size_range,
                           cnv_use=dup_use[dup_use$case.control=="control",],
                           info=TADA_sup_info1, prior=TADA_sup_info1$prior.dn.ptv,
                           beta.dn=beta.dn, n=dn_CNV_con,del_dup_adj=del_dup_adj)
BF_dn_dup <- pmax(1, BF_dn_dup_prob/BF_dn_dup_sib)
TADA_sup_info1$BF_dn_dup <- BF_dn_dup

#################################### Case-control + inherited CNV ######################################

## General formula for cc/in CNV BF
bayes.factor.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1, rho0, nu0){
  marglik0.cc <- evidence.null.cc(x.cc, n.cc, rho0, nu0)
  marglik1.cc <- evidence.alt.cc(x.cc, n.cc, gamma.cc, rho1, nu1)
  BF.cn <- marglik1.cc$cn / marglik0.cc$cn
  if(is.na(BF.cn)){BF.cn <- 1}
  BF.ca <- marglik1.cc$ca / marglik0.cc$ca
  BF <- BF.cn * BF.ca
  return(BF = BF)
}

## Helper function for cc/in CNV BF — null model
evidence.null.cc <- function(x.cc, n.cc, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x.cc$cn, rho0, nu0/(nu0+n.cc$cn)))
  marglik0.case.log <- log(dnbinom(x.cc$ca, rho0+x.cc$cn, (nu0+n.cc$cn)/(nu0+n.cc$ca+n.cc$cn)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  return(list(cn = exp(marglik0.ctrl.log), ca = exp(marglik0.case.log)))
}

## BF for case–control / inherited CNVs
BF_CC_CNV <- function(cnv_size_range, cnv_use, n_case, n_con, info, prior, del_dup_adj=1, nu=5000){ 
  n.cc.cnv <- data.frame(ca=n_case, cn=n_con)
  raw_lBF <- do.call(rbind,lapply(cnv_size_range, function(i, cnv_use){
    cnv_sub <- cnv_use[cnv_use$num_genes==i,]
    ## Compute relative risk
    cnv_sub$relative_risk <- (prior[match(cnv_use$gene_symbol,info$gene_gencodeV33)]-1)/(del_dup_adj+1) 
    rho.ptv = nu * sum(cnv_sub$case.control=="control") / ((2*length(info$gene_gencodeV33)-i)*n.cc.cnv[2])
    gene_combos <- unique(cnv_sub$gene_symbol)
    
    ## Compute per-gene BF
    lBF_cc <- as.numeric(do.call(rbind,lapply(gene_combos, function(x){
      tot_lbf <- log(bayes.factor.cc(
        x.cc = data.frame(
          ca=sum(cnv_sub$gene_symbol==x & cnv_sub$case.control=="case"), 
          cn=sum(cnv_sub$gene_symbol==x & cnv_sub$case.control=="control")),
        n.cc = n.cc.cnv,
        gamma.cc = cnv_sub$relative_risk[cnv_sub$gene_symbol==x][1], 
        rho1 = as.numeric(rho.ptv), nu1 = nu, 
        rho0 = as.numeric(rho.ptv), nu0 = nu))
      return(tot_lbf)
    })))
    names(lBF_cc) <- gene_combos
    
    lBF_cc <- do.call(rbind,lapply(1:length(lBF_cc),function(x){
      evidence <- data.frame(gene = names(lBF_cc)[x], lBF = lBF_cc[x])
      return(evidence)
    }))
  },cnv_use=cnv_use))
  
  sum_lBF <- by(raw_lBF$lBF, raw_lBF$gene, sum)
  BF <- rep(1,length(info$gene_gencodeV33))
  BF[match(names(sum_lBF),info$gene_gencodeV33,nomatch=1)] <- exp(as.numeric(sum_lBF))
  return(BF)
}

## Extract cc CNVs
CNV_SRG_cc <- CNV_SRG_info_non_na1[which(CNV_SRG_info_non_na1$Modified_Inheritence=="Unknown"),]

cc_in_CNV_case <- 15581
cc_in_CNV_con <- 6107

## Compute cc CNV BF (deletion + duplication)
BF_cc_del <- BF_CC_CNV(cnv_size_range, cnv_use=CNV_SRG_cc[CNV_SRG_cc$call=="DEL",],
                       n_case=cc_in_CNV_case, n_con=cc_in_CNV_con,
                       info=TADA_sup_info1, prior=TADA_sup_info1$prior.cc.ptv,
                       del_dup_adj=1, nu=5000)
TADA_sup_info1$BF_cc_del <- BF_cc_del

BF_cc_dup <- BF_CC_CNV(cnv_size_range, cnv_use=CNV_SRG_cc[CNV_SRG_cc$call=="DUP",],
                       n_case=cc_in_CNV_case, n_con=cc_in_CNV_con,
                       info=TADA_sup_info1, prior=TADA_sup_info1$prior.cc.ptv,
                       del_dup_adj=del_dup_adj, nu=5000)
TADA_sup_info1$BF_cc_dup <- BF_cc_dup

## Extract inherited CNVs
CNV_SRG_in <- CNV_SRG_info_non_na1[which(CNV_SRG_info_non_na1$Modified_Inheritence == "Inherited"),]

## Compute inherited CNV BFs
BF_in_del <- BF_CC_CNV(cnv_size_range, cnv_use=CNV_SRG_in[CNV_SRG_in$call=="DEL",],
                       n_case=cc_in_CNV_case, n_con=cc_in_CNV_con,
                       info=TADA_sup_info1, prior=TADA_sup_info1$prior.cc.ptv,
                       del_dup_adj=1, nu=5000)
TADA_sup_info1$BF_in_del <- BF_in_del

BF_in_dup <- BF_CC_CNV(cnv_size_range, cnv_use=CNV_SRG_in[CNV_SRG_in$call=="DUP",],
                       n_case=cc_in_CNV_case, n_con=cc_in_CNV_con,
                       info=TADA_sup_info1, prior=TADA_sup_info1$prior.cc.ptv,
                       del_dup_adj=del_dup_adj, nu=5000)
TADA_sup_info1$BF_in_dup <- BF_in_dup

#################### Filter SRGs from info table ####################

### Extract SRGs present in CNV_SRG
CNV_gene_id_list <- list(CNV_SRG$hgnc_genes)
CNV_gene_id_vector <- as.vector(unlist(CNV_gene_id_list))

info_cnv_SRGs <- TADA_sup_info1[which(TADA_sup_info1$hgnc_symbol %in% CNV_gene_id_vector),]
info_cnv_SRGs_remian <- c("hgnc_symbol","BF_dn_del","BF_dn_dup","BF_cc_del","BF_cc_dup","BF_in_del","BF_in_dup")
info_cnv_SRGs1 <- info_cnv_SRGs[,colnames(info_cnv_SRGs) %in% info_cnv_SRGs_remian]
head(info_cnv_SRGs1)
write.csv(info_cnv_SRGs1,file="./CNV_SRG_BF.csv",row.names = FALSE)

#################### Compute integrated BF (SNV + CNV) ####################

########### General functions
Bayesian.FDR <- function(BF, pi0) {
  i.order = order(BF, decreasing = TRUE)
  BF = BF[i.order]
  pi <- 1 - pi0
  q <- pi*BF/(1-pi+pi*BF)
  q0 <- 1 - q
  FDR = cumsum(q0)/(1:length(BF))
  FDR[i.order] = FDR 
  return(FDR = FDR)
}

q2p <- function(x){
  ec <- x
  ec <- sort(ec)
  ec <- ec*(1:length(ec))/(length(ec))/max(ec)
  return(p = ec[names(x)])
}

### Merge SNV and CNV BF tables
CNV_SRG_fi <- read.csv('./CNV_SRG_BF.csv',header=TRUE) 
SNV_SRG_fi <- read.csv('./SNV_SRG_BF.csv',header=TRUE)

CNV_SNV_BF_merge <- merge(x=CNV_SRG_fi,y=SNV_SRG_fi,by="hgnc_symbol",all.x=TRUE,all.y=TRUE)
CNV_SNV_BF_merge[is.na(CNV_SNV_BF_merge)] <- 1.0

CNV_SNV_BF_merge_an <- merge(CNV_SNV_BF_merge,TADA_sup_info_for_annotation,by="hgnc_symbol",all.x=TRUE)
write.csv(CNV_SNV_BF_merge_an,file="./CNV_SNV_SRG_BF.csv",row.names = FALSE)

CNV_SNV_BF_merge <- read.csv('/SNV_CNV_SRG_BF.csv',header=TRUE)

###### Compute combined BF 
BF_asd <- cbind(
  pmax(CNV_SNV_BF_merge$BF_dn_ptv * CNV_SNV_BF_merge$BF_in_ptv * CNV_SNV_BF_merge$BF_cc_ptv, 1),
  pmax(CNV_SNV_BF_merge$BF_dn_misb * CNV_SNV_BF_merge$BF_in_misb * CNV_SNV_BF_merge$BF_cc_misb, 1),
  pmax(CNV_SNV_BF_merge$BF_dn_misa * CNV_SNV_BF_merge$BF_in_misa * CNV_SNV_BF_merge$BF_cc_misa, 1),
  pmax(CNV_SNV_BF_merge$BF_dn_del * CNV_SNV_BF_merge$BF_in_del * CNV_SNV_BF_merge$BF_cc_del, 1),
  pmax(CNV_SNV_BF_merge$BF_dn_dup * CNV_SNV_BF_merge$BF_in_dup * CNV_SNV_BF_merge$BF_cc_dup, 1)
)

# q value
qval_asd <- Bayesian.FDR(apply(BF_asd, 1, prod), pi0 = 1 - 0.05)
names(qval_asd) <- CNV_SNV_BF_merge$hgnc_symbol
sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd < x))

qval_asd1 <- do.call(rbind,lapply(1:length(qval_asd),function(x){
  bb <- data.frame(hgnc_symbol = names(qval_asd)[x], qval_asd = qval_asd[x])
  return(bb)
}))
rownames(qval_asd1) <- NULL

# p value
ec_threshold <- .05/length(qval_asd)
pval_asd <- q2p(qval_asd)
sum(pval_asd <= ec_threshold)

pval_asd1 <- do.call(rbind,lapply(1:length(pval_asd),function(x){
  aa <- data.frame(hgnc_symbol = names(pval_asd)[x], pval_asd = pval_asd[x])
  return(aa)
}))
rownames(pval_asd1) <- NULL

SRG_SNV_CNV_BF_FDR <- merge(x=CNV_SNV_BF_merge,y=qval_asd1,by="hgnc_symbol")

SRG_SNV_CNV_BF_FDR_pval <- merge(x=SRG_SNV_CNV_BF_FDR,y=pval_asd1,by="hgnc_symbol")
SRG_SNV_CNV_BF_FDR_pval_an <- merge(SRG_SNV_CNV_BF_FDR_pval,TADA_sup_info_for_annotation,by="hgnc_symbol",all.x=TRUE)
write.csv(SRG_SNV_CNV_BF_FDR_pval_an,file="./SNV_CNV_SRG_BF_FDR_pval_20250801.csv",row.names = FALSE)

## Annotate ASD risk group from SFARI
SNV_CNV_SRG_BF_FDR_pval_20240822 <- read.csv("./SNV_CNV_SRG_BF_FDR_pval.csv",header=TRUE)
SFARI_ASD_gene <- read.csv("./SFARI_ASD_gene_note.csv",header=TRUE)
SFARI_ASD_gene_need <- SFARI_ASD_gene[,c("hgnc_symbol","ensembl.ID","ASD.risk.group")]

SNV_CNV_SRG_BF_FDR_pval_20240822_ASD_note <- merge(
  SNV_CNV_SRG_BF_FDR_pval_20240822, SFARI_ASD_gene_need,
  by="hgnc_symbol",all.x=TRUE)
write.csv(SNV_CNV_SRG_BF_FDR_pval_20240822_ASD_note,file="./SNV_CNV_SRG_BF_FDR_pval_ASD_note.csv",row.names=FALSE)

## Add sensory type annotation
BF_info <- read.csv("./SNV_CNV_SRG_BF_FDR_pval_ASD_note.csv",header=TRUE)
sensory_info <- read.csv("genecards_hgnc_sensory_info.csv",header=TRUE)

BF_sensory <- merge(sensory_info, BF_info,
                    by.x="Gene_hgnc", by.y="hgnc_symbol", all.y=TRUE)
write.csv(BF_sensory,file="./BF_sensory.csv",row.names=FALSE)

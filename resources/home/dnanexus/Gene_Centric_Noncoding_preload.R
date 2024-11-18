Gene_Centric_Noncoding_preload <- function(chr,gene_name,category=c("all_categories","downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),
                                           genofile,obj_nullmodel,
                                           dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                           dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                           dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                           dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                           rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                           QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                           Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                           Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                           SPA_p_filter=TRUE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  if(category=="all_categories")
  {
    results <- noncoding_preload(chr,gene_name,genofile,obj_nullmodel,
                                 dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                 dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                 dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                 dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                 rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                 rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                 SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="downstream")
  {
    results <- downstream(chr,gene_name,genofile,obj_nullmodel,
                          rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                          rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                          QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                          Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                          SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="upstream")
  {
    results <- upstream(chr,gene_name,genofile,obj_nullmodel,
                        rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                        rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                        SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="UTR")
  {
    results <- UTR(chr,gene_name,genofile,obj_nullmodel,
                   rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                   rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                   QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                   Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                   SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="promoter_CAGE")
  {
    results <- promoter_CAGE_preload(chr,gene_name,genofile,obj_nullmodel,
                                     dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                     rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                     rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                     SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="promoter_DHS")
  {
    results <- promoter_DHS_preload(chr,gene_name,genofile,obj_nullmodel,
                                    dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                    rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                    rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                    SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="enhancer_CAGE")
  {
    results <- enhancer_CAGE_preload(chr,gene_name,genofile,obj_nullmodel,
                                     dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                     rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                     rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                     SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="enhancer_DHS")
  {
    results <- enhancer_DHS_preload(chr,gene_name,genofile,obj_nullmodel,
                                    dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                    rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                    rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                    SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  return(results)
}


noncoding_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                              dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                              dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                              dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                              dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                              rare_maf_cutoff=0.01,rv_num_cutoff=2,
                              rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                              QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                              Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                              Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                              SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno

  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA))
  {
    use_SPA <- obj_nullmodel$use_SPA
  }else
  {
    use_SPA <- FALSE
  }

  #####################################
  #   Gene Info
  ## get SNV id
  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant")
  {
    SNVlist <- filter == "PASS"
  }

  if(variant_type=="SNV")
  {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }

  if(variant_type=="Indel")
  {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }

  variant.id <- seqGetData(genofile, "variant.id")

  rm(filter)
  gc()

  ########################################
  #   Downstream
  ## downstream SNVs
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  is.in <- (GENCODE.Category=="downstream")&(SNVlist)
  variant.id.downstream <- variant.id[is.in]

  seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)

  rm(variant.id.downstream)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

  rm(GENCODE.Info)
  gc()

  rm(variant_gene_num)
  gc()

  Gene <- as.character(unlist(GENCODE.Info.split))

  rm(GENCODE.Info.split)
  gc()

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results_downstream <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- rep(NA,4)
    results_temp[3] <- "downstream"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results_downstream <- rbind(results_downstream,results_temp)
  }

  if(!is.null(results_downstream))
  {
    if(!use_SPA)
    {
      colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
      colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_downstream)[(dim(results_downstream)[2]-1):dim(results_downstream)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
      colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_downstream)[dim(results_downstream)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  ########################################
  #   Upstream

  is.in <- (GENCODE.Category=="upstream")&(SNVlist)
  variant.id.upstream <- variant.id[is.in]

  seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)

  rm(variant.id.upstream)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

  rm(GENCODE.Info)
  gc()

  rm(variant_gene_num)
  gc()

  Gene <- as.character(unlist(GENCODE.Info.split))

  rm(GENCODE.Info.split)
  gc()

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results_upstream <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- rep(NA,4)
    results_temp[3] <- "upstream"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results_upstream <- rbind(results_upstream,results_temp)
  }

  if(!is.null(results_upstream))
  {
    if(!use_SPA)
    {
      colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
      colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_upstream)[(dim(results_upstream)[2]-1):dim(results_upstream)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
      colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_upstream)[dim(results_upstream)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  ########################################################
  #                UTR

  is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
  variant.id.UTR <- variant.id[is.in]

  rm(GENCODE.Category)
  gc()

  seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)

  rm(variant.id.UTR)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")

  rm(GENCODE.Info)
  gc()

  # Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[seq(1,length(z),2)]))
  Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))

  rm(GENCODE.Info.split)
  gc()

  variant.id.SNV <- seqGetData(genofile, "variant.id")

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results_UTR <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- rep(NA,4)
    results_temp[3] <- "UTR"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results_UTR <- rbind(results_UTR,results_temp)
  }

  if(!is.null(results_UTR))
  {
    if(!use_SPA)
    {
      colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
      colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_UTR)[(dim(results_UTR)[2]-1):dim(results_UTR)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
      colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_UTR)[dim(results_UTR)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  #############################################
  #   Promoter-CAGE

  ### Gene
  is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromCAGE[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results_promoter_CAGE <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- dfPromCAGEVarGene.SNV[1,1:4]
    results_temp[3] <- "promoter_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results_promoter_CAGE <- rbind(results_promoter_CAGE,results_temp)
  }

  if(!is.null(results_promoter_CAGE))
  {
    if(!use_SPA)
    {
      colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
      colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_promoter_CAGE)[(dim(results_promoter_CAGE)[2]-1):dim(results_promoter_CAGE)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
      colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_promoter_CAGE)[dim(results_promoter_CAGE)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  ##################################################
  #       Promoter-DHS

  ### Gene
  is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromrOCRs[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results_promoter_DHS <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- dfPromrOCRsVarGene.SNV[1,1:4]
    results_temp[3] <- "promoter_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results_promoter_DHS <- rbind(results_promoter_DHS ,results_temp)
  }

  if(!is.null(results_promoter_DHS))
  {
    if(!use_SPA)
    {
      colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
      colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_promoter_DHS)[(dim(results_promoter_DHS)[2]-1):dim(results_promoter_DHS)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
      colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_promoter_DHS)[dim(results_promoter_DHS)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  ###########################################
  #        Enhancer-CAGE

  ### Gene
  is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerCAGE[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results_enhancer_CAGE <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- dfHancerCAGEVarGene.SNV[1,1:4]
    results_temp[3] <- "enhancer_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results_enhancer_CAGE <- rbind(results_enhancer_CAGE,results_temp)
  }

  if(!is.null(results_enhancer_CAGE))
  {
    if(!use_SPA)
    {
      colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
      colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_enhancer_CAGE)[(dim(results_enhancer_CAGE)[2]-1):dim(results_enhancer_CAGE)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
      colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_enhancer_CAGE)[dim(results_enhancer_CAGE)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  ##################################################
  #       Enhancer-DHS

  ### Gene
  is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerrOCRs[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results_enhancer_DHS <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- dfHancerrOCRsVarGene.SNV[1,1:4]
    results_temp[3] <- "enhancer_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results_enhancer_DHS <- rbind(results_enhancer_DHS,results_temp)
  }

  if(!is.null(results_enhancer_DHS))
  {
    if(!use_SPA)
    {
      colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
      colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_enhancer_DHS)[(dim(results_enhancer_DHS)[2]-1):dim(results_enhancer_DHS)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
      colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_enhancer_DHS)[dim(results_enhancer_DHS)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  ############################################
  #           results

  results_noncoding <- list(upstream=results_upstream,downstream=results_downstream,UTR=results_UTR,
                            promoter_CAGE=results_promoter_CAGE,promoter_DHS=results_promoter_DHS,
                            enhancer_CAGE=results_enhancer_CAGE,enhancer_DHS=results_enhancer_DHS)

  return(results_noncoding)

}


promoter_CAGE_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                                  dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                  Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                  SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno

  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA))
  {
    use_SPA <- obj_nullmodel$use_SPA
  }else
  {
    use_SPA <- FALSE
  }

  ### Gene
  is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromCAGE[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- dfPromCAGEVarGene.SNV[1,1:4]
    results_temp[3] <- "promoter_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }


    results <- rbind(results,results_temp)
  }

  if(!is.null(results))
  {
    if(!use_SPA)
    {
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  return(results)
}


promoter_DHS_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                                 dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                 QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                 Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                 Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                 SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno

  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA))
  {
    use_SPA <- obj_nullmodel$use_SPA
  }else
  {
    use_SPA <- FALSE
  }

  ### Gene
  is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromrOCRs[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- dfPromrOCRsVarGene.SNV[1,1:4]
    results_temp[3] <- "promoter_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results <- rbind(results,results_temp)
  }


  if(!is.null(results))
  {
    if(!use_SPA)
    {
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  return(results)
}


enhancer_CAGE_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                                  dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                  Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                  SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno

  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA))
  {
    use_SPA <- obj_nullmodel$use_SPA
  }else
  {
    use_SPA <- FALSE
  }

  ### Gene
  is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerCAGE[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- dfHancerCAGEVarGene.SNV[1,1:4]
    results_temp[3] <- "enhancer_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results <- rbind(results,results_temp)
  }

  if(!is.null(results))
  {
    if(!use_SPA)
    {
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  return(results)
}


enhancer_DHS_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                                 dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                 QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                 Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                 Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                 SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno

  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA))
  {
    use_SPA <- obj_nullmodel$use_SPA
  }else
  {
    use_SPA <- FALSE
  }

  ### Gene
  is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerrOCRs[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- NULL
  if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
  {
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
  }

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  pvalues <- 0
  if(n_pheno == 1)
  {
    if(!use_SPA)
    {
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                           rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
    }else
    {
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                      rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,
                                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else
  {
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
  }

  results <- c()
  if(inherits(pvalues, "list"))
  {
    results_temp <- dfHancerrOCRsVarGene.SNV[1,1:4]
    results_temp[3] <- "enhancer_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant

    if(!use_SPA)
    {
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else
    {
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }

    results <- rbind(results,results_temp)
  }

  if(!is.null(results))
  {
    if(!use_SPA)
    {
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
    }else
    {
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }

  seqResetFilter(genofile)

  return(results)
}


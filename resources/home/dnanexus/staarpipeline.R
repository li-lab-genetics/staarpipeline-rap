args <- commandArgs(TRUE)
### mandatory
outfile <- args[1]
test.type <- args[2] # Null, Single, Gene_Centric_Coding, Gene_Centric_Coding_incl_ptv, Gene_Centric_Noncoding, ncRNA, Sliding_Window, SCANG
### optional
pheno.file <- args[3]
grm.file <- args[4]
nullobj.file <- args[5]
agds.file <- args[6]
annotation.name.catalog.file <- args[7]
phenotype <- args[8]
pheno_id <- args[9]
covariates <- args[10]
het_vars <- args[11]
random_time_slope <- args[12]
user_cores <- as.numeric(args[13])
min.mac <- as.numeric(args[14])
max.maf <- as.numeric(args[15])
min.rv.num <- as.numeric(args[16])
max.rv.num <- as.numeric(args[17])
max.rv.num.prefilter <- as.numeric(args[18])
sliding_window_length <- as.numeric(args[19])
QC_label <- args[20]
variant_type <- args[21]
geno_missing_imputation <- args[22]
Annotation_dir <- args[23]
Use_annotation_weights <- args[24]
Annotation_name <- args[25]

test.type.vals <- c("Null", "Single", "Gene_Centric_Coding", "Gene_Centric_Coding_incl_ptv", "Gene_Centric_Noncoding", "ncRNA", "Sliding_Window", "SCANG")
if(!test.type %in% test.type.vals) stop("Error: test.type must be Null, Single, Gene_Centric_Coding, Gene_Centric_Coding_incl_ptv, Gene_Centric_Noncoding, ncRNA, Sliding_Window, or SCANG")
if(Use_annotation_weights == "YES") {
  Use_annotation_weights = TRUE
} else {
  Use_annotation_weights = FALSE
}
Annotation_name <- unlist(strsplit(Annotation_name,","))

if(test.type == "Null") {
  cat("Only fit the null model, no association tests will be performed...\n")
  cat("Output null model object:", paste0(outfile, ".Rdata"), "\n")
  cat("Input null model object:", nullobj.file, "\n")
  cat("Phenotype file:", pheno.file, "\n")
  cat("Relatedness matrix file:", grm.file, "\n")
  cat("Phenotype:", phenotype, "\n")
  cat("ID:", pheno_id, "\n")
  cat("Covariates:", covariates, "\n")
  cat("Grouping variable in heteroscedastic linear mixed models:", het_vars, "\n")
  cat("Time variable in random slope longitudinal models:", random_time_slope, "\n")
  cat("The following arguments will be ignored:\n")
  cat("\tGenotype and functional annotation all-in-one GDS (AGDS) file:", agds.file, "\n")
  cat("\tName and the corresponding channel name in the AGDS file:", annotation.name.catalog.file, "\n")
  cat("\tUser requested running the analysis on", user_cores, "cores\n")
  cat("\tMinimum minor allele count to be included for single variant test:", min.mac, "\n")
  cat("\tMaximum minor allele frequency to be included for variant-set test:", max.maf, "\n")
  cat("\tMinimum number of variants of analyzing a given variant-set:", min.rv.num, "\n")
  cat("\tMaximum number of variants of analyzing a given variant-set:", max.rv.num, "\n")
  cat("\tMaximum number of variants before extracting the genotype matrix:", max.rv.num.prefilter, "\n")
  cat("\tSliding window size (bp) to be used in sliding window test:", sliding_window_length, "\n")
  cat("\tChannel name of the QC label in the AGDS file:", QC_label, "\n")
  cat("\tVariants included in the analysis:", variant_type, "\n")
  cat("\tMethod of handling missing genotypes:", geno_missing_imputation, "\n")
  cat("\tChannel name of the annotations in the AGDS file:", Annotation_dir, "\n")
  cat("\tUse annotations as weights or not:", Use_annotation_weights, "\n")
  cat("\tAnnotations used in STAAR:", Annotation_name, "\n")
  user_cores <- 1
} else {
  cat("Output file prefix:", outfile, "\n")
  cat("Type of test:", test.type, "\n")
  if(agds.file == "NO_AGDS_FILE") stop("Error: AGDS genotype file can only be missing when test.type is Null: fitting the null model only")
  cat("Null model object:", nullobj.file, "\n")
  cat("Genotype and functional annotation all-in-one GDS (AGDS) file:", agds.file, "\n")
  cat("Name and the corresponding channel name in the AGDS file:", annotation.name.catalog.file, "\n")
  cat("User requested running the analysis on", user_cores, "cores\n")
  cat("Minimum minor allele count to be included for single variant test:", min.mac, "\n")
  cat("Maximum minor allele frequency to be included for variant-set test:", max.maf, "\n")
  cat("\tMinimum number of variants of analyzing a given variant-set:", min.rv.num, "\n")
  cat("\tMaximum number of variants of analyzing a given variant-set:", max.rv.num, "\n")
  cat("\tMaximum number of variants before extracting the genotype matrix:", max.rv.num.prefilter, "\n")
  cat("Sliding window size (bp) to be used in sliding window test:", sliding_window_length, "\n")
  cat("Channel name of the QC label in the AGDS file:", QC_label, "\n")
  cat("Variants included in the analysis:", variant_type, "\n")
  cat("Method of handling missing genotypes:", geno_missing_imputation, "\n")
  cat("Channel name of the annotations in the AGDS file:", Annotation_dir, "\n")
  cat("Use annotations as weights or not:", Use_annotation_weights, "\n")
  cat("Annotations used in STAAR:", Annotation_name, "\n")
  cat("The following arguments will be ignored:\n")
  cat("\tPhenotype file:", pheno.file, "\n")
  cat("\tRelatedness matrix file:", grm.file, "\n")
  cat("\tPhenotype:", phenotype, "\n")
  cat("\tID:", pheno_id, "\n")
  cat("\tCovariates:", covariates, "\n")
  cat("\tGrouping variable in heteroscedastic linear mixed models:", het_vars, "\n")
  cat("\tTime variable in random slope longitudinal models:", random_time_slope, "\n")
}

if(user_cores > 1) Sys.setenv(MKL_NUM_THREADS = 1)
suppressMessages(library(gdsfmt))
suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(STAAR))
suppressMessages(library(MultiSTAAR))
suppressMessages(library(SCANG))
suppressMessages(library(STAARpipeline))
suppressMessages(library(parallel))

if(nullobj.file == "NO_NULL_OBJ") {
  if(pheno.file == "NO_PHENO_FILE") stop("Error: pheno.file must be specified if no null model object is provided")
  if(phenotype == "NO_PHENO") stop("Error: phenotype must be specified if no null model object is provided")
  pheno <- read.csv(pheno.file, as.is=T)
  if(grm.file == "NO_GRM_FILE") {
    warning("Note: grm.file for relatedness not specified... Assuming individuals are uncorrelated...")
    kins <- NULL
  } else {
    kins <- get(load(grm.file))
    if(!is.null(kins) && !class(kins) %in% c("matrix", "list")) {
      if(is.null(attr(class(kins), "package"))) stop("Error: relatedness matrix (matrices) must be a matrix or a list.")
      else if(attr(class(kins), "package") != "Matrix") stop("Error: if \"kins\" is a sparse matrix, it must be created using the Matrix package.")
    }
  }
  pheno.value <- unique(pheno[,phenotype])
  is.binary <- FALSE
  if(length(pheno.value)==2) is.binary <- all(sort(pheno.value)==c(0,1))
  family <- if(is.binary) binomial(link="logit") else gaussian(link="identity")
  if(is.binary) {
    cat("Running analysis for a binary trait using the logistic mixed model, the following arguments will be ignored:\n")
    cat("\tGrouping variable in heteroscedastic linear mixed models:", het_vars, "\n")
  }
  formula <- if(covariates=="NO_COVAR") paste(phenotype, "~ 1") else paste(phenotype, "~", paste(unlist(strsplit(covariates, ",")), collapse = " + "))
  nullobj <- fit_nullmodel(as.formula(formula), data=pheno, kins=kins, id=pheno_id, random.slope=if(random_time_slope=="NO_RANDOM_TIME_SLOPE") NULL else random_time_slope, groups=if(het_vars=="NO_HET_VARS") NULL else het_vars, family=family, verbose=TRUE)
  rm(pheno, kins); gc()
} else {
  cat("Null model object already exists, the following arguments will be ignored:\n")
  cat("\tPhenotype file:", pheno.file, "\n")
  cat("\tRelatedness matrix file:", grm.file, "\n")
  cat("\tPhenotype:", phenotype, "\n")
  cat("\tID:", pheno_id, "\n")
  cat("\tCovariates:", covariates, "\n")
  cat("\tGrouping variable in heteroscedastic linear mixed models:", het_vars, "\n")
  cat("\tTime variable in random slope longitudinal models:", random_time_slope, "\n")
  ls.tmp <- ls()
  nullobj <- get(load(nullobj.file))
  if(class(nullobj) == "GENESIS.nullMixedModel") {
    nullmod$sample.id <- row.names(nullmod$model.matrix)
    nullobj <- genesis2staar_nullmodel(nullmod)
    rm(nullmod, pheno); gc()
  }
}
if(test.type == "Null") {
  save(nullobj, file = paste0(outfile, ".Rdata"))
} else if(test.type == "Gene_Centric_Coding") {
  if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when test.type is Gene_Centric_Coding")
  Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
  cat("Performing gene-centric test for coding functional categories, the following arguments will be ignored:\n")
  cat("\tMinimum minor allele count to be included for single variant test:", min.mac, "\n")
  cat("\tSliding window size (bp) to be used in sliding window test:", sliding_window_length, "\n")
  rm(list=setdiff(ls(), c("outfile", "nullobj", "agds.file", "max.maf", "min.rv.num", "max.rv.num", "max.rv.num.prefilter", "QC_label", "variant_type", "geno_missing_imputation", "Annotation_dir", "Annotation_name_catalog", "Use_annotation_weights", "Annotation_name", "user_cores"))); gc()

  genofile <- seqOpen(agds.file)

  ## Chr
  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  chr <- CHR[1]
  rm(CHR); gc()
  ## genes info
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  sub_seq_num <- dim(genes_info_chr)[1]

  ## array_id
  sub_seq_id <- 1:sub_seq_num

  gene_centric_coding_dnanexus <- function(genes_info_chr,kk,chr,genofile,obj_nullmodel,rare_maf_cutoff,
                                           rv_num_cutoff,rv_num_cutoff_max,rv_num_cutoff_max_prefilter,
                                           QC_label,variant_type,geno_missing_imputation,
                                           Annotation_dir,Annotation_name_catalog,
                                           Use_annotation_weights,Annotation_name,silent)
  {
    gene_name <- genes_info_chr[kk,1]
    results <- try(Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,
                                       rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent),silent=TRUE)
    return(results)
  }

  out <- mclapply(sub_seq_id,gene_centric_coding_dnanexus,genes_info_chr=genes_info_chr,chr=chr,genofile=genofile,obj_nullmodel=nullobj,rare_maf_cutoff=max.maf,
                  rv_num_cutoff=min.rv.num,rv_num_cutoff_max=max.rv.num,rv_num_cutoff_max_prefilter=max.rv.num.prefilter,
                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=TRUE,mc.cores=user_cores)

  rm(list=setdiff(ls(), c("out", "outfile"))); gc()
  save(out, file = paste0(outfile, ".Rdata"))
  } else if(test.type == "Gene_Centric_Coding_incl_ptv") {
    if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when test.type is Gene_Centric_Coding")
    Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
    cat("Performing gene-centric test for coding functional categories including protein-truncating variants, the following arguments will be ignored:\n")
    cat("\tMinimum minor allele count to be included for single variant test:", min.mac, "\n")
    cat("\tSliding window size (bp) to be used in sliding window test:", sliding_window_length, "\n")
    rm(list=setdiff(ls(), c("outfile", "nullobj", "agds.file", "max.maf", "min.rv.num", "max.rv.num", "max.rv.num.prefilter", "QC_label", "variant_type", "geno_missing_imputation", "Annotation_dir", "Annotation_name_catalog", "Use_annotation_weights", "Annotation_name", "user_cores"))); gc()

    genofile <- seqOpen(agds.file)

    ## Chr
    CHR <- as.numeric(seqGetData(genofile, "chromosome"))
    chr <- CHR[1]
    rm(CHR); gc()
    ## genes info
    genes_info_chr <- genes_info[genes_info[,2]==chr,]
    sub_seq_num <- dim(genes_info_chr)[1]

    ## array_id
    sub_seq_id <- 1:sub_seq_num

    gene_centric_coding_incl_ptv_dnanexus <- function(genes_info_chr,kk,chr,genofile,obj_nullmodel,rare_maf_cutoff,
                                                      rv_num_cutoff,rv_num_cutoff_max,rv_num_cutoff_max_prefilter,
                                                      QC_label,variant_type,geno_missing_imputation,
                                                      Annotation_dir,Annotation_name_catalog,
                                                      Use_annotation_weights,Annotation_name,silent)
    {
      gene_name <- genes_info_chr[kk,1]
      results <- try(Gene_Centric_Coding(chr=chr,gene_name=gene_name,category="all_categories_incl_ptv",genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,
                                         rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                         QC_label=QC_label,variant_type="variant",geno_missing_imputation=geno_missing_imputation,
                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                         Use_annotation_weights=FALSE,Annotation_name=Annotation_name,silent=silent),silent=TRUE)
      return(results)
    }

    out <- mclapply(sub_seq_id,gene_centric_coding_incl_ptv_dnanexus,genes_info_chr=genes_info_chr,chr=chr,genofile=genofile,obj_nullmodel=nullobj,rare_maf_cutoff=max.maf,
                    rv_num_cutoff=min.rv.num,rv_num_cutoff_max=max.rv.num,rv_num_cutoff_max_prefilter=max.rv.num.prefilter,
                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=TRUE,mc.cores=user_cores)

    rm(list=setdiff(ls(), c("out", "outfile"))); gc()
    save(out, file = paste0(outfile, ".Rdata"))
} else if(test.type == "Gene_Centric_Noncoding") {
  if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when test.type is Gene_Centric_Noncoding")
  Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
  cat("Performing gene-centric test for noncoding functional categories, the following arguments will be ignored:\n")
  cat("\tMinimum minor allele count to be included for single variant test:", min.mac, "\n")
  cat("\tSliding window size (bp) to be used in sliding window test:", sliding_window_length, "\n")
  rm(list=setdiff(ls(), c("outfile", "nullobj", "agds.file", "max.maf", "min.rv.num", "max.rv.num", "max.rv.num.prefilter", "QC_label", "variant_type", "geno_missing_imputation", "Annotation_dir", "Annotation_name_catalog", "Use_annotation_weights", "Annotation_name", "user_cores"))); gc()

  suppressMessages(library(GenomicFeatures))
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  source("Gene_Centric_Noncoding_preload.R")

  genofile <- seqOpen(agds.file)

  ## Chr
  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  chr <- CHR[1]
  rm(CHR)
  ## genes info
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  sub_seq_num <- dim(genes_info_chr)[1]

  ## array_id
  sub_seq_id <- 1:sub_seq_num

  gene_centric_noncoding_dnanexus <- function(genes_info_chr,kk,chr,genofile,obj_nullmodel,
                                              dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                              dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                              dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                              dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                              rare_maf_cutoff,rv_num_cutoff,
                                              rv_num_cutoff_max,rv_num_cutoff_max_prefilter,
                                              QC_label,variant_type,geno_missing_imputation,
                                              Annotation_dir,Annotation_name_catalog,
                                              Use_annotation_weights,Annotation_name,silent){
    gene_name <- genes_info_chr[kk,1]
    results <- try(Gene_Centric_Noncoding_preload(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                  dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                                  dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                                  dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                                  dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                                  rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                  rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent),silent=TRUE)
    return(results)
  }

  #########################################################
  #             Promoter_CAGE
  #########################################################

  varid <- seqGetData(genofile, "variant.id")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

  #Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
  CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
  CAGEBvt <- CAGEAnno!=""
  CAGEidx <- which(CAGEBvt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[CAGEidx])
  seqSetFilter(genofile,promGobj,intersect=TRUE)
  CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
  ##obtain variants info
  CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
  CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
  CAGEvref <- as.character(seqGetData(genofile,"$ref"))
  CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
  dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)

  rm(varid)
  gc()

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
  variant.id.SNV.PromCAGE <- variant.id[SNVlist]

  dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
  dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
  dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
  dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)

  seqResetFilter(genofile)

  rm(dfPromCAGEVarGene)
  gc()

  #########################################################
  #             Promoter_DHS
  #########################################################

  varid <- seqGetData(genofile, "variant.id")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

  # Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
  rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
  rOCRsBvt <- rOCRsAnno!=""
  rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[rOCRsidx])

  seqSetFilter(genofile,promGobj,intersect=TRUE)
  rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
  ## obtain variants info
  rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
  rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
  rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
  rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
  dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

  rm(varid)
  gc()

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
  variant.id.SNV.PromrOCRs <- variant.id[SNVlist]

  dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
  dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
  dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
  dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

  seqResetFilter(genofile)

  rm(dfPromrOCRsVarGene)
  gc()

  #########################################################
  #             Enhancer_CAGE
  #########################################################

  varid <- seqGetData(genofile, "variant.id")

  #Now extract the GeneHancer with CAGE Signal Overlay
  genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  genehancer <- genehancerAnno!=""

  CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
  CAGE <- CAGEAnno!=""
  CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
  CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])

  # variants that covered by whole GeneHancer without CAGE overlap.
  genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
  enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
  enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
  enhancervpos <- as.numeric(seqGetData(genofile,"position"))
  enhancervref <- as.character(seqGetData(genofile,"$ref"))
  enhancervalt <- as.character(seqGetData(genofile,"$alt"))
  dfHancerCAGEVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

  rm(varid)
  gc()

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
  variant.id.SNV.HancerCAGE <- variant.id[SNVlist]

  dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
  dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
  dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
  dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)

  seqResetFilter(genofile)

  rm(dfHancerCAGEVarGene)
  gc()

  #########################################################
  #             Enhancer_DHS
  #########################################################

  varid <- seqGetData(genofile, "variant.id")

  #Now extract the GeneHancer with rOCRs Signal Overlay
  genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  genehancer <- genehancerAnno!=""

  rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
  rOCRs <- rOCRsAnno!=""
  rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
  rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
  # variants that covered by whole GeneHancer without rOCRs overlap.

  genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
  enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
  enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
  enhancervpos <- as.numeric(seqGetData(genofile,"position"))
  enhancervref <- as.character(seqGetData(genofile,"$ref"))
  enhancervalt <- as.character(seqGetData(genofile,"$alt"))
  dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

  rm(varid)
  gc()

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
  variant.id.SNV.HancerrOCRs <- variant.id[SNVlist]

  dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
  dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
  dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
  dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)

  seqResetFilter(genofile)

  rm(dfHancerrOCRsVarGene)
  gc()

  out <- mclapply(sub_seq_id,gene_centric_noncoding_dnanexus,genes_info_chr=genes_info_chr,chr=chr,genofile=genofile,obj_nullmodel=nullobj,
                  dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                  dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                  dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                  dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                  rare_maf_cutoff=max.maf,rv_num_cutoff=min.rv.num,
                  rv_num_cutoff_max=max.rv.num,rv_num_cutoff_max_prefilter=max.rv.num.prefilter,
                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=TRUE,mc.cores=user_cores)

  rm(list=setdiff(ls(), c("out", "outfile"))); gc()
  save(out, file = paste0(outfile, ".Rdata"))
} else if(test.type == "ncRNA") {
  if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when test.type is ncRNA")
  Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
  cat("Performing ncRNA test, the following arguments will be ignored:\n")
  cat("\tMinimum minor allele count to be included for single variant test:", min.mac, "\n")
  cat("\tSliding window size (bp) to be used in sliding window test:", sliding_window_length, "\n")
  rm(list=setdiff(ls(), c("outfile", "nullobj", "agds.file", "max.maf", "min.rv.num", "max.rv.num", "max.rv.num.prefilter", "QC_label", "variant_type", "geno_missing_imputation", "Annotation_dir", "Annotation_name_catalog", "Use_annotation_weights", "Annotation_name", "user_cores"))); gc()

  genofile <- seqOpen(agds.file)

  ## Chr
  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  chr <- CHR[1]
  rm(CHR)
  ## genes info
  ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]
  sub_seq_num <- dim(ncRNA_gene_chr)[1]

  ## array_id
  sub_seq_id <- 1:sub_seq_num

  gene_centric_ncRNA_dnanexus <- function(ncRNA_gene_chr,kk,chr,genofile,obj_nullmodel,rare_maf_cutoff,
                                          rv_num_cutoff,rv_num_cutoff_max,rv_num_cutoff_max_prefilter,
                                          QC_label,variant_type,geno_missing_imputation,
                                          Annotation_dir,Annotation_name_catalog,
                                          Use_annotation_weights,Annotation_name,silent)
  {
    gene_name <- ncRNA_gene_chr[kk,2]
    results <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,
                         rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent),silent=TRUE)
    return(results)
  }

  out <- mclapply(sub_seq_id,gene_centric_ncRNA_dnanexus,ncRNA_gene_chr=ncRNA_gene_chr,chr=chr,genofile=genofile,obj_nullmodel=nullobj,rare_maf_cutoff=max.maf,
                  rv_num_cutoff=min.rv.num,rv_num_cutoff_max=max.rv.num,rv_num_cutoff_max_prefilter=max.rv.num.prefilter,
                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=TRUE,mc.cores=user_cores)

  rm(list=setdiff(ls(), c("out", "outfile"))); gc()
  save(out, file = paste0(outfile, ".Rdata"))
}else if(test.type == "Sliding_Window") {
  if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE" && Use_annotation_weights) stop("Error: Annotation name catalog file cannot be missing when test.type is Sliding_Window and Use_annotation_weights is YES")
  Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
  cat("Performing sliding window test, the following arguments will be ignored:\n")
  cat("\tMinimum minor allele count to be included for single variant test:", min.mac, "\n")
  rm(list=setdiff(ls(), c("outfile", "nullobj", "agds.file", "max.maf", "min.rv.num", "max.rv.num", "max.rv.num.prefilter", "sliding_window_length", "QC_label", "variant_type", "geno_missing_imputation", "Annotation_dir", "Annotation_name_catalog", "Use_annotation_weights", "Annotation_name", "user_cores"))); gc()

  genofile <- seqOpen(agds.file)

  ## Chr
  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  chr <- CHR[1]
  rm(CHR)

  ## jobs_num
  jobs_num <- matrix(rep(0,66),nrow=22)

  filter <- seqGetData(genofile, QC_label)
  SNVlist <- filter == "PASS"

  position <- as.numeric(seqGetData(genofile, "position"))
  position_SNV <- position[SNVlist]

  jobs_num[chr,1] <- chr
  jobs_num[chr,2] <- min(position[SNVlist])
  jobs_num[chr,3] <- max(position[SNVlist])

  colnames(jobs_num) <- c("chr","start_loc","end_loc")
  jobs_num <- as.data.frame(jobs_num)

  ## start_loc and end_loc
  start_loc <- jobs_num$start_loc[chr]
  end_loc <- start_loc + (sliding_window_length/2)*20 - 1

  ## sub-sequence num
  sub_seq_num <- ceiling((jobs_num$end_loc[chr] - jobs_num$start_loc[chr] + 1)/20/(sliding_window_length/2))
  sub_seq_id <- 1:sub_seq_num

  sliding_window_dnanexus <- function(kk,chr,start_loc,end_loc,sliding_window_length,genofile,obj_nullmodel,rare_maf_cutoff,
                                      rv_num_cutoff,rv_num_cutoff_max,rv_num_cutoff_max_prefilter,
                                      QC_label,variant_type,geno_missing_imputation,
                                      Annotation_dir,Annotation_name_catalog,
                                      Use_annotation_weights,Annotation_name,silent)
  {
    start_loc_sub <- start_loc + (sliding_window_length/2)*20*(kk-1)
    end_loc_sub <- end_loc + (sliding_window_length/2)*20*(kk-1) + (sliding_window_length/2)

    end_loc_sub <- min(end_loc_sub,jobs_num$end_loc[chr])

    results <- try(Sliding_Window(chr=chr,start_loc=start_loc_sub,end_loc=end_loc_sub,sliding_window_length=sliding_window_length,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,
                                  rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,type="multiple",
                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent),silent=TRUE)

    return(results)
  }

  out <- mclapply(sub_seq_id,sliding_window_dnanexus,chr=chr,start_loc=start_loc,end_loc=end_loc,sliding_window_length=sliding_window_length,genofile=genofile,obj_nullmodel=nullobj,rare_maf_cutoff=max.maf,
                  rv_num_cutoff=min.rv.num,rv_num_cutoff_max=max.rv.num,rv_num_cutoff_max_prefilter=max.rv.num.prefilter,
                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=TRUE,mc.cores=user_cores)

  rm(list=setdiff(ls(), c("out", "outfile"))); gc()
  save(out, file = paste0(outfile, ".Rdata"))
}else if(test.type == "SCANG") {
  if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE" && Use_annotation_weights) stop("Error: Annotation name catalog file cannot be missing when test.type is SCANG and Use_annotation_weights is YES")
  Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
  cat("Performing dynamic window test (SCANG), the following arguments will be ignored:\n")
  cat("\tMinimum minor allele count to be included for single variant test:", min.mac, "\n")
  cat("\tMinimum number of variants of analyzing a given variant-set:", min.rv.num, "\n")
  cat("\tMaximum number of variants of analyzing a given variant-set:", max.rv.num, "\n")
  cat("\tMaximum number of variants before extracting the genotype matrix:", max.rv.num.prefilter, "\n")
  cat("\tSliding window size (bp) to be used in sliding window test:", sliding_window_length, "\n")
  rm(list=setdiff(ls(), c("outfile", "nullobj", "agds.file", "max.maf", "QC_label", "variant_type", "geno_missing_imputation", "Annotation_dir", "Annotation_name_catalog", "Use_annotation_weights", "Annotation_name", "user_cores"))); gc()

  genofile <- seqOpen(agds.file)

  ## Chr
  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  chr <- CHR[1]
  rm(CHR)

  ## jobs_num
  jobs_num <- matrix(rep(0,66),nrow=22)

  filter <- seqGetData(genofile, QC_label)
  SNVlist <- filter == "PASS"

  position <- as.numeric(seqGetData(genofile, "position"))
  position_SNV <- position[SNVlist]

  jobs_num[chr,1] <- chr
  jobs_num[chr,2] <- min(position[SNVlist])
  jobs_num[chr,3] <- max(position[SNVlist])

  colnames(jobs_num) <- c("chr","start_loc","end_loc")
  jobs_num <- as.data.frame(jobs_num)

  nullobj <- staar2scang_nullmodel(nullobj)

  ## sub-sequence num
  sub_seq_num <- ceiling((jobs_num$end_loc[chr] - jobs_num$start_loc[chr] + 1)/0.5e6)
  sub_seq_id <- 1:sub_seq_num

  scang_dnanexus <- function(kk,chr,genofile,obj_nullmodel,rare_maf_cutoff,
                             QC_label,variant_type,geno_missing_imputation,
                             Annotation_dir,Annotation_name_catalog,
                             Use_annotation_weights,Annotation_name,silent)
  {
    start_loc <- (kk-1)*0.5e6 + jobs_num$start_loc[chr]
    end_loc <- min(start_loc + 0.5e6 - 1, jobs_num$end_loc[chr])

    results <- try(Dynamic_Window_SCANG(chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,
                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent),silent=TRUE)

    return(results)
  }

  out <- mclapply(sub_seq_id,scang_dnanexus,chr=chr,genofile=genofile,obj_nullmodel=nullobj,rare_maf_cutoff=max.maf,
                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=TRUE,mc.cores=user_cores)

  rm(list=setdiff(ls(), c("out", "outfile"))); gc()
  save(out, file = paste0(outfile, ".Rdata"))
}else {
  cat("Performing single variant test, the following arguments will be ignored:\n")
  cat("\tMaximum minor allele frequency to be included for variant-set test:", max.maf, "\n")
  cat("\tMinimum number of variants of analyzing a given variant-set:", min.rv.num, "\n")
  cat("\tMaximum number of variants of analyzing a given variant-set:", max.rv.num, "\n")
  cat("\tMaximum number of variants before extracting the genotype matrix:", max.rv.num.prefilter, "\n")
  cat("\tSliding window size (bp) to be used in sliding window test:", sliding_window_length, "\n")
  cat("\tChannel name of the annotations in the AGDS file:", Annotation_dir, "\n")
  cat("\tUse annotations as weights or not:", Use_annotation_weights, "\n")
  cat("\tAnnotations used in STAAR:", Annotation_name, "\n")
  rm(list=setdiff(ls(), c("outfile", "nullobj", "agds.file", "min.mac", "QC_label", "variant_type", "geno_missing_imputation", "user_cores"))); gc()

  genofile <- seqOpen(agds.file)

  ## Chr
  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  chr <- CHR[1]
  rm(CHR)

  ## jobs_num
  jobs_num <- matrix(rep(0,66),nrow=22)

  filter <- seqGetData(genofile, QC_label)
  SNVlist <- filter == "PASS"

  position <- as.numeric(seqGetData(genofile, "position"))
  position_SNV <- position[SNVlist]

  jobs_num[chr,1] <- chr
  jobs_num[chr,2] <- min(position[SNVlist])
  jobs_num[chr,3] <- max(position[SNVlist])

  colnames(jobs_num) <- c("chr","start_loc","end_loc")
  jobs_num <- as.data.frame(jobs_num)

  ## start_loc and end_loc
  start_loc <- jobs_num$start_loc[chr]
  end_loc <- start_loc + 0.5e6 - 1

  ## sub-sequence num
  sub_seq_num <- ceiling((jobs_num$end_loc[chr] - jobs_num$start_loc[chr] + 1)/0.5e6)
  sub_seq_id <- 1:sub_seq_num

  individual_analysis_dnanexus <- function(kk,chr,start_loc,end_loc,genofile,obj_nullmodel,mac_cutoff,subset_variants_num,
                                           QC_label,variant_type,geno_missing_imputation)
  {
    start_loc_sub <- start_loc + 0.5e6*(kk-1)
    end_loc_sub <- end_loc + 0.5e6*(kk-1)

    end_loc_sub <- min(end_loc_sub,jobs_num$end_loc[chr])

    results <- try(Individual_Analysis(chr=chr,start_loc=start_loc_sub,end_loc=end_loc_sub,genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=mac_cutoff,subset_variants_num=subset_variants_num,
                                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation))
    return(results)
  }

  subset_variants_num <- 5e3
  out <- mclapply(sub_seq_id,individual_analysis_dnanexus,chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,obj_nullmodel=nullobj,mac_cutoff=min.mac,subset_variants_num=subset_variants_num,
                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,mc.cores=user_cores)

  rm(list=setdiff(ls(), c("out", "outfile"))); gc()
  save(out, file = paste0(outfile, ".Rdata"))
}


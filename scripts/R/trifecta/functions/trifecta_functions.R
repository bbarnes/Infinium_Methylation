
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               TripleCrown::
#                         CG# Database Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Command Line Options Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )

# Tidyverse Core Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("stringr", quietly = TRUE) ))
suppressWarnings(suppressPackageStartupMessages( 
  base::require("readr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("glue",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("furrr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("purrr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("magrittr", quietly = TRUE) ) )

# Parallel Processing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

# Plotting Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("ggplot2", quietly = TRUE) ) )

# Matrix Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("matrixStats", quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("scales",      quietly = TRUE) ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Common Human Abbreviations for String Variables::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
TAB2 <- "\t\t"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Local Run Time Defaults::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

trifecta_options = function( pars, args,
                             
                             vb=0, vt=4, tc=1, tt=NULL,
                             fun_tag='trifecta_options') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  eflag <- FALSE
  wflag <- FALSE
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb>=vt+2) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   args={args}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Parse Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts <- base::list()
  
  opts$run_name     <- NULL
  opts$out_path     <- NULL
  
  opts$design_key <- 'Seq_ID'
  opts$design_seq <- 'Forward_Sequence'
  opts$design_prb <- 'Probe_Key'
  opts$score_key  <- 'UnMethyl_Final_Score,Methyl_Final_Score'
  opts$strand_fr_key <- 'Methyl_Allele_FR_Strand'
  opts$strand_co_key <- 'Methyl_Allele_CO_Strand'
  opts$strand_tb_key <- 'Methyl_Allele_TB_Strand'
  opts$cpg_cnt_key   <- 'Methyl_Underlying_CpG_Count'
  
  # Probe Filter Parameters::
  opts$min_prb_score <- "0.3,0.2"
  opts$min_cpg_rank  <- "3,2,1,0"
  opts$min_scr_rank  <- "0.6,0.5,0.4,0.3"
  opts$strandFR      <- 'F,R'
  opts$strandCO      <- 'C,O'
  opts$pick_best     <- TRUE
  opts$unique_cpgs   <- TRUE
  
  opts$sample_max   <- 0
  opts$platform     <- "EPIC"
  
  # Genomic Reference Parameters
  opts$ref_path     <- NULL
  opts$ref_file     <- NULL
  opts$ref_build    <- NULL
  opts$ref_source   <- NULL
  opts$ref_species  <- NULL
  
  opts$temp_off_csv  <- NULL
  opts$canonical_csv <- file.path( pars$aux_path, "canonical_cgn_top_grp.csv.gz")
  opts$temp_off_csv  <- file.path( pars$aux_path, "template_offset.csv.gz")
  
  if ( is.null(opts$temp_off_csv) ||
       !file.exists(opts$temp_off_csv) ) opts$temp_off_csv <- 
    file.path( pars$aux_path, "triplecrown/template_offset.csv.gz")
  
  opts$Rscript      <- "Rscript"
  opts$single       <- FALSE
  opts$cluster      <- FALSE
  opts$parallel     <- FALSE
  opts$track_time   <- TRUE
  opts$clean        <- FALSE
  opts$reload       <- 0
  opts$verbose      <- 3
  
  if (opts$verbose > 0)
    cat(glue::glue("[{pars$prgm_tag}]: Starting; {pars$prgm_tag}.{RET2}"))
  
  if (pars$run_mode == 'RStudio') {
    
    pars$local_run_type <- NULL
    
    # Read Options and Test Cases from pre-defined function
    # loc_dat <- lighthoof_local_defaults( opts = opt,
    #                                      pars = par,
    #                                      vb = opts$verbose )
    # opt <- loc_dat$opt
    # par <- loc_dat$par
    
    pars$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R")
    opts$top_path <- pars$top_path
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Params Defaults::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( is.null(pars$version) ) pars$version <- '1'
    if ( is.null(pars$version_key) ) pars$version_key <- "v"
    opts$version  <- paste0(pars$version_key,pars$version)
    
    opts$out_path  <- file.path( pars$top_path, 'scratch' )
    opts$imp_path  <- file.path( pars$top_path, 'data/improbe' )
    opts$ann_path  <- file.path( pars$top_path, 'data/annotation' )
    opts$man_path  <- file.path( pars$top_path, 'data/manifests' )
    opts$idat_path <- file.path( pars$top_path, 'data/idats' )
    
    opts$bsmap_opt <- "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
    opts$bsmap_dir <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz"
    opts$bsmap_exe <- "bsmapz"
    
    # opts$improbe_path <- file.path( opts$top_path, "tmp/improbe-output-10.v2.tsv.gz" )
    # opts$improbe_path <- file.path( opts$top_path, "Projects/Envonik/GRCr00.improbeDesignOutputs.100k.bed.gz" )
    opts$improbe_path <- file.path( opts$top_path, "Projects/Envonik/region-designs/GRCr00.improbeDesignOutputs.bed.gz" )
    if (  is.null(pars$run_name) ) opts$run_name = "Chicken"
    
    opts$improbe_path <- file.path( opts$top_path, "Projects/Envonik/region-designs/Cho.improbeDesignOutputs.bed.gz" )
    if (  is.null(pars$run_name) ) opts$run_name = "Cho"
    
    opts$improbe_path <- file.path( opts$top_path, "Projects/Envonik/region-designs/Crawfish.improbeDesignOutputs.bed.gz" )
    if (  is.null(pars$run_name) ) opts$run_name = "Crawfish"
    
    if ( is.null(opts$run_name) ) opts$run_name = "Evonik"
    
    if ( !is.null(pars$run_name) ) opts$run_name = pars$run_name
    if ( par$prgm_tag == "stable_embark_add_on_cgn" ) opts$run_name = "Embark"
    
    opts$pick_best   <- TRUE
    opts$unique_cpgs <- TRUE
    
    # opts$ref_source <- paste( "NCBI", sep = ',' )
    # opts$ref_source <- paste( "UCSC", sep = ',' )
    opts$ref_source <- paste( "Both", sep = ',' )
    
    opts$run_name <- paste(opts$run_name,opts$ref_source,opts$version, sep='-')
    if ( opts$ref_source == "UCSC" ) {
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$chr_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "hg38.fa.gz",
                              "hg19.fa.gz",
                              "hg18.fa.gz",
                              "mm10.fa.gz",
                              "galGal5.fa.gz",
                              "criGriChoV1.fa.gz",
                              "pvirPacbioDovetail.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 "Gallus_gallus",
                                 "ChineseHamster",
                                 "MarbledCrayfish",
                                 sep = ',' )
      
      opts$ref_build <- paste( "hg38",
                               "hg19",
                               "hg18",
                               "mm10",
                               "galgal5",
                               "criGrichov1",
                               "GRCc00",
                               sep = ',' )
      
      
    } else if ( opts$ref_source == "NCBI" ) {
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$chr_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "GRCh38.genome.fa.gz",
                              "GRCh37.genome.fa.gz",
                              "GRCh36.genome.fa.gz",
                              "GRCm38.genome.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 sep = ',' )
      
      opts$ref_build <- paste( "GRCh38",
                               "GRCh37",
                               "GRCh36",
                               "GRCm10",
                               sep = ',' )
      
    } else if ( opts$ref_source == "Both" ) {
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "GRCh38.genome.fa.gz",
                              "GRCh37.genome.fa.gz",
                              "GRCh36.genome.fa.gz",
                              "GRCm38.genome.fa.gz",
                              
                              "galGal5.fa.gz",
                              "criGriChoV1.fa.gz",
                              "pvirPacbioDovetail.fa.gz",
                              
                              "canFam3.genome.fa.gz",
                              "canFam6.genome.fa.gz",
                              
                              # "hg38.fa.gz",
                              # "hg19.fa.gz",
                              # "hg18.fa.gz",
                              # "mm10.fa.gz",
                              
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 
                                 "Gallus_gallus",
                                 "ChineseHamster",
                                 "MarbledCrayfish",
                                 
                                 "Canis_lupus_familiaris",
                                 "Canis_lupus_familiaris",

                                 # "Homo_sapiens",
                                 # "Homo_sapiens",
                                 # "Homo_sapiens",
                                 # "Mus_musculus",
                                 
                                 sep = ',' )
      
      opts$ref_build <- paste( "GRCh38",
                               "GRCh37",
                               "GRCh36",
                               "GRCm10",
                               
                               "galgal5",
                               "criGrichov1",
                               "GRCc00",
                               
                               "canFam3",
                               "canFam6",
                               
                               # "hg38",
                               # "hg19",
                               # "hg18",
                               # "mm10",
                               
                               sep = ',' )
      
    } else {
      eflag <- TRUE
      errs_mssg <- glue::glue("Unsupported default ref_source = {opts$ref_source}")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    # Define tru-cgn-map csv::
    opts$cgn_tru_csv <- file.path( 
      pars$top_path, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-top/bins-100/map.csv.gz" )
    
    # opts$single       = TRUE
    opts$single       = FALSE
    
    opts$cluster      = FALSE
    opts$parallel     = TRUE
    opts$track_time   = TRUE
    opts$clean        = FALSE
    opts$reload       = 1
    opts$verbose      = 3
    
    opts$workflow = "parse_genome,parse_chromosome,update_database"
    
  } else if (pars$run_mode == 'Command_Line') {
    
    options_list <- list(
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Run Time Version Options:: 
      #                       Platform, Genome Build, etc
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--run_name"), type="character", default=opts$run_name, 
        help=glue::glue(
          "Build run name.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--workflow"), type="character", default=opts$workflow, 
        help=glue::glue(
          "Workflow(s) to be executed.",
          "{RET}{TAB2} e.g. parse_genome,parse_chromosome,update_database, etc.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--top_path"), type="character", default=opts$top_path,
        help=glue::glue(
          "Top directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--out_path"), type="character", default=opts$out_path,
        help=glue::glue(
          "Build output directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                 Options:: Loci Variation/EWAS-Swifthoof
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--improbe_path"), type="character", default=opts$improbe_path,
        help=glue::glue(
          "Path to improbe design output file.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--design_key"), type="character", default=opts$design_key,
        help=glue::glue(
          "Improbe Design Key Name (Seq_ID).",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--design_seq"), type="character", default=opts$design_seq,
        help=glue::glue(
          "Improbe Design Sequence Name (Forward_Sequence, Top_Sequence).",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--design_prb"), type="character", default=opts$design_prb,
        help=glue::glue(
          "Improbe Probe Design Type Name (Probe_Key).",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--design_srd"), type="character", default=opts$design_srd,
        help=glue::glue(
          "Improbe Design Sequence Strand Letters. If null will guess base ",
          "on Design Sequence Name (FR,TB).",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--score_key"), type="character", default=opts$score_key,
        help=glue::glue(
          "Score Column Names. If multiple columns the min of all of them will be used.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--strand_fr_key"), type="character", default=opts$strand_fr_key,
        help=glue::glue(
          "Strand Forward/Opposite column name.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--strand_co_key"), type="character", default=opts$strand_co_key,
        help=glue::glue(
          "Strand Converted/Opposite column name.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--strand_tb_key"), type="character", default=opts$strand_tb_key,
        help=glue::glue(
          "Strand Top/Bot column name.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--cpg_cnt_key"), type="character", default=opts$cpg_cnt_key,
        help=glue::glue(
          "Underlying CpG count column name.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--sample_max"), type="integer", default=opts$sample_max, 
        help=glue::glue(
          "Maximum Number of Samples to Process.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Program Specific Options:
      #                           improbe Parameters
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--min_prb_score"), type="character", default=opts$min_prb_score, 
        help=glue::glue(
          "Minimum Probe Score for each Converted/Opposite strand.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--min_cpg_rank"), type="character", default=opts$min_cpg_rank, 
        help=glue::glue(
          "Infinium I/II Probe Underlying CpG Count Rank String.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--min_scr_rank"), type="character", default=opts$min_scr_rank, 
        help=glue::glue(
          "Infinium I/II Probe Score Rank String.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--strandFR"), type="character", default=opts$strandFR, 
        help=glue::glue(
          "Target Strand Forward/Reverse to design (F,R).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--strandCO"), type="character", default=opts$strandCO, 
        help=glue::glue(
          "Target Strand Converted/Opposite to design (C,O).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--pick_best"), action="store_true", default=opts$pick_best, 
        help=glue::glue(
          "Boolean variable to only return the best scoring probe for each strand.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--unique_cpgs"), action="store_true", default=opts$unique_cpgs, 
        help=glue::glue(
          "Boolean variable to only select unique cpg in the genome.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Program Specific Options:
      #                            Refernce Genome
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--ref_file"), type="character", default=opts$ref_file, 
        help=glue::glue(
          "Reference Genome fasta file name(s).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_build"), type="character", default=opts$ref_build, 
        help=glue::glue(
          "Reference Genome build names(s).",
          "{RET}{TAB2} e.g. GRch38,GRCh37,GRCh36,GRCm37.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_source"), type="character", default=opts$ref_source, 
        help=glue::glue(
          "Reference Source(s).",
          "{RET}{TAB2} e.g. UCSC, NCBI.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_species"), type="character", default=opts$ref_species, 
        help=glue::glue(
          "Reference Specie(s).",
          "{RET}{TAB2} e.g. Homo_sapiens, Mus_musculus, Rattus_norvegicus, ",
          "SARS-CoV-2.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--canonical_csv"), type="character", default=opts$canonical_csv, 
        help=glue::glue(
          "CSV file(s) containg canonical (already defined) CG#/Top Sequence ",
          "Template definitions.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--temp_off_csv"), type="character", default=opts$temp_off_csv, 
        help=glue::glue(
          "CSV file containg Top Sequence Template start and end offsets from ",
          "CG# genomic postion.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Run Time Mode Options::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # optparse::make_option(
      #   c("--Rscript"), type="character", default=opts$Rscript,
      #   help=glue::glue(
      #     "Rscript executable path.",
      #     "{RET}{TAB2}[ default = %default ]" ),
      #   metavar="character"),
      
      # Process Parallel/Cluster Parameters::
      optparse::make_option(
        c("--single"), action="store_true", default=opts$single, 
        help=glue::glue(
          "Boolean variable to run a single sample on a single-core.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--parallel"), action="store_true", default=opts$parallel, 
        help=glue::glue(
          "Boolean variable to run parallel on multi-core.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--cluster"), action="store_true", default=opts$cluster,
        help=glue::glue(
          "Boolean variable to run jobs on cluster.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      # Run=time Options::
      optparse::make_option(
        c("--track_time"), action="store_true", default=opts$track_time,
        help=glue::glue(
          "Boolean variable tack run times.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--clean"), action="store_true", default=opts$clean, 
        help=glue::glue(
          "Boolean variable to run a clean build.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--reload"), type="integer", default=opts$reload, 
        help=glue::glue(
          "Integer value to reload intermediate files (for testing).{RET}",
          "{TAB2} Zero indicates no-reloading, higher numbers more reloads.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer"),
      
      # Verbosity level::
      optparse::make_option(
        c("-v", "--verbose"), type="integer", default=opts$verbose, 
        help=glue::glue(
          "Verbosity level: 0-5 (5 is very verbose).",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer")
    )
    
    opt_parser = optparse::OptionParser( option_list = options_list )
    opts = optparse::parse_args( opt_parser )
    
    if ( is.null(opts$top_path) && !is.null(pars$top_path) )
      opts$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R")
    
  } else {
    stop( glue::glue("{RET}[{pars$prgm_tag}]: ERROR: Unrecognized run_mode = ",
                     "'{pars$run_mode}'!{RET2}") )
  }
  
  if ( opts$clean ) opts$reload <- -1
  
  
  errs_mssg <- glue::glue("improbe_path='{opts$improbe_path}' does NOT exist")
  if ( !file.exists(opts$improbe_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  ret_cnt <- length(opts)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

# End of file

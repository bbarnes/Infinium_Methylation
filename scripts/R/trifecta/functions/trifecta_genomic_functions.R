
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Trifecta:: Genomic
#             improbe (Infinium Methylation Probe) Design Functions::
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
#                         Trifecta Static Tables::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Move these to /dat/auxilary in git repository...
#

#
# Seq_Tag Index Map::
#
seq_tag_srd_tib <- dplyr::tribble(
  ~Strand_TB, ~Strand_CO, ~Next_Base, ~Tag_Idx,
  "T",        "C",        "A",         0,
  "T",        "C",        "T",         1,
  "T",        "C",        "C",         2,
  "T",        "C",        "G",         3,
  # ----- ----- ----- ----- ----- ----- -----
  "T",        "O",        "A",         4,
  "T",        "O",        "T",         5,
  "T",        "O",        "C",         6,
  "T",        "O",        "G",         7,
  # ----- ----- ----- ----- ----- ----- -----
  "B",        "C",        "A",         8,
  "B",        "C",        "T",         9,
  "B",        "C",        "C",        10,
  "B",        "C",        "G",        11,
  # ----- ----- ----- ----- ----- ----- -----
  "B",        "O",        "A",        12,
  "B",        "O",        "T",        13,
  "B",        "O",        "C",        14,
  "B",        "O",        "G",        15,
  # ----- ----- ----- ----- ----- ----- -----
) %>% dplyr::mutate(
  dplyr::across(Strand_TB, as.character ),
  dplyr::across(Strand_CO, as.character ),
  dplyr::across(Next_Base, as.character ),
  dplyr::across(Tag_Idx,   as.integer ) ) # %>% dplyr::mutate( Idx = dplyr::row_number() )

# TBD:: Tables like this and file format col_specs should be defined upfront
#  and passed to to functions!!!
#
prb_desIdx_tib <- dplyr::tribble(
  ~Prb_Type, ~Infinium, ~AlleleAB,
  "cg",      1,         "A",
  "ch",      1,         "A",
  "rs",      1,         "A",
  "cg",      2,         "A",
  "ch",      2,         "A",
  "rs",      2,         "A",
  # ----- ----- ----- ----- -----
  "cg",      1,         "B",
  "ch",      1,         "B",
  "rs",      1,         "B",
  "cg",      2,         "B",
  "ch",      2,         "B",
  "rs",      2,         "B",
  # ----- ----- ----- ----- -----
) %>% dplyr::mutate(
  dplyr::across(Prb_Type, as.character),
  dplyr::across(Infinium,  as.integer),
  dplyr::across(AlleleAB, as.character) ) %>%
  dplyr::mutate( Idx = dplyr::row_number() )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#            Load Pre-Defined Probe Look-up Metrics Tables
#
#  TBD:: This should really live somewhere else...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

umd_target_tib <- dplyr::tribble(
  ~Prb_Type, ~Infinium, ~Strand_BS, ~Strand_FR, ~Strand_CO, ~PrbEnd_3p, ~UMD_Type,
  "cg",      1,         "--",       "F",        "C",        "G",        "M",
  "cg",      1,         "-+",       "F",        "O",        "C",        "M",
  "cg",      1,         "+-",       "R",        "C",        "G",        "M",
  "cg",      1,         "++",       "R",        "O",        "C",        "M",
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "cg",      1,         "--",       "F",        "C",        "A",        "U",
  "cg",      1,         "-+",       "F",        "O",        "T",        "U",
  "cg",      1,         "+-",       "R",        "C",        "A",        "U",
  "cg",      1,         "++",       "R",        "O",        "T",        "U",
  
) %>% dplyr::mutate(
  dplyr::across(Prb_Type, as.character),
  dplyr::across(Infinium,  as.integer),
  dplyr::across(Strand_BS, as.character),
  dplyr::across(Strand_FR, as.character),
  dplyr::across(Strand_CO, as.character),
  dplyr::across(PrbEnd_3p, as.character),
  dplyr::across(UMD_Type,  as.character)
)

# probe_offset_tib <- dplyr::tribble(
#   ~Srd_Key, ~FR1_Key, ~FR2_Key, ~FR1,   ~FR2,   ~TB,    ~CO,    ~Des,   ~Offset, ~Length, ~NxbPos,
#   TC1,     FC1     RC1     F     R     T     C     1          0     50     -1
#   BC1,     RC1     FC1     R     F     B     C     1        -49     50      1
#   TO1,     FO1     RO1     F     R     T     O     1        -50     50    -50
#   BO1,     RO1     FO1     R     F     B     O     1          1     50      0
#   TC2,     FC2     RC2     F     R     T     C     2          1     50      0
#   BC2,     RC2     FC2     R     F     B     C     2          0     50      1
#   TO2,     FO2     RO2     F     R     T     O     2        -51     50      0
#   BO2,     RO2     FO2     R     F     B     O     2          2     50      1
# )

top_offset_tib <- dplyr::tribble(
  ~Prb_Type, ~Infinium, ~Strand_BS, ~Strand_TB, ~Strand_CO, ~Top_Offset, ~Nxb_Offset,
  "rs",       1,        "--",       "T",        "C",        50,          0,
  "rs",       1,        "-+",       "T",        "O",        -1,          0,
  "rs",       1,        "+-",       "B",        "C",         0,          0,
  "rs",       1,        "++",       "B",        "O",        49,          0,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "rs",       2,        "--",       "T",        "C",        50,          0,
  "rs",       2,        "-+",       "T",        "O",         0,          0,
  "rs",       2,        "+-",       "B",        "C",        -1,          0,
  "rs",       2,        "++",       "B",        "O",        49,          0,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "ch",       1,        "--",       "T",        "C",        49,          0,
  "ch",       1,        "-+",       "T",        "O",         0,          0,
  "ch",       1,        "+-",       "B",        "C",         0,          0,
  "ch",       1,        "++",       "B",        "O",        49,          0,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "ch",       2,        "--",       "T",        "C",        50,          0,
  "ch",       2,        "-+",       "T",        "O",        -1,          0,
  "ch",       2,        "+-",       "B",        "C",        -1,          0,
  "ch",       2,        "++",       "B",        "O",        50,          0,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "cg",       1,        "--",       "T",        "C",         0,         -1,
  "cg",       1,        "-+",       "T",        "O",       -50,        -50,
  "cg",       1,        "+-",       "B",        "C",       -49,          1,
  "cg",       1,        "++",       "B",        "O",         1,          0,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "cg",       2,        "--",       "T",        "C",         1,          0,
  "cg",       2,        "-+",       "T",        "O",       -51,          0,
  "cg",       2,        "+-",       "B",        "C",         0,          1,
  "cg",       2,        "++",       "B",        "O",         2,          1,
  
) %>% dplyr::mutate(
  dplyr::across(Prb_Type,   as.character),
  dplyr::across(Infinium,   as.integer),
  dplyr::across(Strand_BS,  as.character),
  dplyr::across(Strand_TB,  as.character),
  dplyr::across(Strand_CO,  as.character),
  dplyr::across(Top_Offset, as.integer),
  dplyr::across(Nxb_Offset, as.integer)
)

# template_offset_tib
# dns_pad dns_len din_len ups_len ups_pad tmp_len ext_len
# <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#   1      -2     -60       2      60       2     122     126

tmp_offset_tib <- dplyr::tribble(
  ~dns_pad, ~dns_len, ~din_len, ~ups_len, ~ups_pad, ~tmp_len, ~ext_len,
  -2,       -60,      2,        60,       2,        122,      126,
) %>% dplyr::mutate(
  dplyr::across( dns_pad, as.integer ),
  dplyr::across( dns_len, as.integer ),
  dplyr::across( din_len, as.integer ),
  dplyr::across( ups_len, as.integer ),
  dplyr::across( ups_pad, as.integer ),
  dplyr::across( tmp_len, as.integer ),
  dplyr::across( ext_len, as.integer )
)

cpg_offset_13092021_tib <- dplyr::tribble(
  ~Prb_Type, ~Infinium, ~Strand_BS, ~Strand_FR, ~Strand_CO, ~CG_Offset,
  "rs",       1,         "--",       "F",        "C",        50,
  "rs",       1,         "-+",       "F",        "O",        -1,
  "rs",       1,         "+-",       "R",        "C",         0,
  "rs",       1,         "++",       "R",        "O",        49,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "rs",       2,         "--",       "F",        "C",        50,
  "rs",       2,         "-+",       "F",        "O",         0,
  "rs",       2,         "+-",       "R",        "C",        -1,
  "rs",       2,         "++",       "R",        "O",        49,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "ch",       1,         "--",       "F",        "C",        49,
  "ch",       1,         "-+",       "F",        "O",         0,
  "ch",       1,         "+-",       "R",        "C",         0,
  "ch",       1,         "++",       "R",        "O",        49,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "ch",       2,         "--",       "F",        "C",        50,
  "ch",       2,         "-+",       "F",        "O",        -1,
  "ch",       2,         "+-",       "R",        "C",        -1,
  "ch",       2,         "++",       "R",        "O",        50,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "cg",       1,         "--",       "F",        "C",        48,
  "cg",       1,         "-+",       "F",        "O",         0,
  "cg",       1,         "+-",       "R",        "C",         0,
  "cg",       1,         "++",       "R",        "O",        48,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "cg",       2,        "--",        "F",        "C",        49,
  "cg",       2,        "-+",        "F",        "O",        -1,
  "cg",       2,        "+-",        "R",        "C",        -1,
  "cg",       2,        "++",        "R",        "O",        49,
  
) %>% dplyr::mutate(
  dplyr::across(Prb_Type, as.character),
  dplyr::across(Infinium,  as.integer),
  dplyr::across(Strand_BS, as.character),
  dplyr::across(Strand_FR, as.character),
  dplyr::across(Strand_CO, as.character),
  dplyr::across(CG_Offset, as.integer)
)

cpg_offset_tib <- dplyr::tribble(
  ~Probe_Type, ~Infinium_Design, ~Strand_BS, ~Strand_FR, ~Strand_CO, ~CG_Offset,
  "rs",        1,                "--",       "F",        "C",        50,
  "rs",        1,                "-+",       "F",        "O",        -1,
  "rs",        1,                "+-",       "R",        "C",         0,
  "rs",        1,                "++",       "R",        "O",        49,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "rs",        2,                "--",       "F",        "C",        50,
  "rs",        2,                "-+",       "F",        "O",         0,
  "rs",        2,                "+-",       "R",        "C",        -1,
  "rs",        2,                "++",       "R",        "O",        49,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "ch",        1,                "--",       "F",        "C",        49,
  "ch",        1,                "-+",       "F",        "O",         0,
  "ch",        1,                "+-",       "R",        "C",         0,
  "ch",        1,                "++",       "R",        "O",        49,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "ch",        2,                "--",       "F",        "C",        50,
  "ch",        2,                "-+",       "F",        "O",        -1,
  "ch",        2,                "+-",       "R",        "C",        -1,
  "ch",        2,                "++",       "R",        "O",        50,
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "cg",        1,                "--",       "F",        "C",        48,
  "cg",        1,                "-+",       "F",        "O",        -1, #  0 ~ -1 is 100% correct!
  "cg",        1,                "+-",       "R",        "C",         0,
  "cg",        1,                "++",       "R",        "O",        49, # 48
  # ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  "cg",        2,               "--",        "F",        "C",        49,
  "cg",        2,               "-+",        "F",        "O",        -2, # -1  ~ B02:: -2 is 100% correct!
  "cg",        2,               "+-",        "R",        "C",        -1,
  "cg",        2,               "++",        "R",        "O",        50, # 49  ~ T02 (50 was looking best so far...)
  
) %>% dplyr::mutate(
  dplyr::across(Probe_Type, as.character),
  dplyr::across(Infinium_Design,  as.integer),
  dplyr::across(Strand_BS, as.character),
  dplyr::across(Strand_FR, as.character),
  dplyr::across(Strand_CO, as.character),
  dplyr::across(CG_Offset, as.integer)
)




# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Local Run Time Defaults::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

trifecta_genomic_options = function( pars, args,
                                
                                vb=0, vt=4, tc=1, tt=NULL,
                                fun_tag='trifecta_genomic_options') {
  
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
  opts$pos_path     <- NULL
  opts$imp_path     <- NULL
  
  opts$ref_path     <- NULL
  opts$ref_file     <- NULL
  opts$ref_build    <- NULL
  opts$ref_source   <- NULL
  opts$ref_species  <- NULL
  
  opts$canonical_csv <- file.path( pars$aux_path, "canonical_cgn_top_grp.csv.gz")
  opts$temp_off_csv  <- file.path( pars$prgm_aux_path, "template_offset.csv.gz" )
  
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
    
    opts$out_path  <- file.path( pars$top_path, 'scratch' )
    opts$imp_path  <- file.path( pars$top_path, 'data/improbe' )
    opts$ann_path  <- file.path( pars$top_path, 'data/annotation' )
    opts$man_path  <- file.path( pars$top_path, 'data/manifests' )
    opts$idat_path <- file.path( pars$top_path, 'data/idats' )
    
    opts$run_name = "chicago-v.1.1"
    
    opts$ref_source <- paste( "NCBI", sep = ',' )
    # opts$ref_source <- paste( "UCSC", sep = ',' )
    
    opts$run_name <- paste(opts$run_name,opts$ref_source, sep='-')
    if ( opts$ref_source == "UCSC" ) {
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$chr_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "hg19.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 sep = ',' )
      
      opts$ref_build <- paste( "hg19",
                               sep = ',' )

    } else if ( opts$ref_source == "NCBI" ) {
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$chr_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "GRCh37.genome.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 sep = ',' )
      
      opts$ref_build <- paste( "GRCh37",
                               sep = ',' )
      
    } else {
      errs_mssg <- glue::glue("Unsupported default ref_source = {ref_source}")
      if ( !is.null(files) ) eflag <- FALSE %in% lapply(files,file.exists) %>% unlist()
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    
    opts$single       = TRUE
    # opts$single       = FALSE
    
    opts$cluster      = FALSE
    opts$parallel     = TRUE
    opts$track_time   = TRUE
    opts$clean        = FALSE
    opts$reload       = 1
    opts$verbose      = 3
    
    opts$workflow = "parse_genome,parse_chromosome,update_database"
    
    opts$imap_csv <- NULL
    if ( dir.exists( opts$imp_path ) ) opts$imap_csv <- 
      file.path( opts$imp_path, 'scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-top/bins-100/map.csv.gz' )

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
        c("--out_path"), type="character", default=opts$out_path,
        help=glue::glue(
          "Build output directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #               Trifecta Genomics Specific input parameters::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--pos_path"), type="character", default=opts$pos_path,
        help=glue::glue(
          "Target posistion path(s). (e.g. bed file)",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),

      optparse::make_option(
        c("--imap_csv"), type="character", default=opts$imap_csv,
        help=glue::glue(
          "Improbe cgn database split mapping file path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                     Other common input parameters::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--ref_path"), type="character", default=opts$ref_path, 
        help=glue::glue(
          "Reference Genome directory path(s).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--chr_path"), type="character", default=opts$chr_path, 
        help=glue::glue(
          "Reference Genome directory path(s) for individual chromosomes.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
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
      
      optparse::make_option(
        c("--Rscript"), type="character", default=opts$Rscript,
        help=glue::glue(
          "Rscript executable path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
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
    
  } else {
    stop( glue::glue("{RET}[{pars$prgm_tag}]: ERROR: Unrecognized run_mode = ",
                     "'{pars$run_mode}'!{RET2}") )
  }
  
  ret_cnt <- length(opts)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

# End of file

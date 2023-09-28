
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                imAssign::
#                    ( formally known as triplecrown ) 
#                         cg# Database Generation
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Options, Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("doParallel", quietly = TRUE) ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Set Run Environment:: RStudio/Command-Line
#                            Source All Functions
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par  <- list()
args <- commandArgs(trailingOnly = FALSE)

par$src_path <- NULL
par$run_mode <- args[1]
par$date_str <- Sys.Date() %>% as.character()
par$prgm_dir <- 'imAssign'
par$prgm_tag <- 'imAssign'
par$verbose  <- 3
local_paths  <- c( 
  "/Users/bbarnes/Documents/tools/imSuite/scripts/R",
  "/home/jsommer/methylation/imAssign/scripts/R"
)

if ( par$run_mode == "RStudio" ) {
  for (path in local_paths) if ( dir.exists(path) ) par$src_path <- path
} else {
  par$exe_path <- 
    base::normalizePath( base::substring( args[grep("--file=", args)], 8 ) )
  par$scr_path <- base::dirname( par$exe_path )
  par$src_path <- base::dirname( par$scr_path )
  par$run_mode <- 'Command_Line'
}
stopifnot( length(par$src_path) > 0 )
stopifnot( dir.exists(par$src_path) )

par$fun_path <- file.path( par$src_path, "functions" )
stopifnot( dir.exists(par$fun_path) )

src_ret <- base::lapply( base::list.files( 
  par$fun_path, pattern = "\\.R$", full.names = TRUE ), source )
par <- source_functions( pars = par, rcpp = 0, vb = par$verbose )
par <- params_check( pars = par, args = args, 
                     prgm_aux_check = TRUE, vb = par$verbose )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$version <- 0
par$version <- 1
# par$version <- 2
# par$version <- 3

# par$run_name <- "sesame"

opt <- NULL
# opt <- triplecrown_options( pars = par, args = args, vb = par$verbose )
opt <- imAssign_options( pars = par, args = args, vb = par$verbose )
vb  <- opt$verbose

vt  <- 0
tc  <- 0

p0  <- vb >= vt + 0
p1  <- vb >= vt + 1
p2  <- vb >= vt + 2
p4  <- vb >= vt + 4
p8  <- vb >= vt + 8

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c( 'run_mode', 
               'src_path', 'scr_path', 'exe_path', 'prgm_dir', 'prgm_tag' )
opt_reqs <- c( 'out_path', 'ref_path', 'ref_file', 'ref_build', 'ref_species',
               'Rscript', 'verbose' )

#
# TBD:: Update docker defaults after building new docker branch::
# TBD:: Add auxilary files to to program_init check::
#
opt$rcpp <- 2
prgm_dat <- NULL
prgm_dat <- program_init( name = par$prgm_tag,
                          opts = opt, opt_reqs = opt_reqs,
                          pars = par, par_reqs = par_reqs, 
                          rcpp = opt$rcpp,
                          vb = opt$verbose, vt=3, tc=0 )

opt <- NULL
opt <- prgm_dat$opt
par <- NULL
par <- prgm_dat$par

opt_tib <- NULL
opt_tib <- prgm_dat$opt_tib
par_tib <- NULL
par_tib <- prgm_dat$par_tib

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Initialize Run Objects
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tt <- timeTracker$new()

# For accumulating failed probes and why::
#  TBD:: Add this to timeTracker()
#
error_ledgar <- NULL

pmssg <- glue::glue("[{par$prgm_tag}]:")
pwarn <- glue::glue("{RET}[{par$prgm_tag}]: Warning:")
perrs <- glue::glue("{RET}[{par$prgm_tag}]: ERROR:")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Pre-defined Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Load Auxiliary CG# Database Files::
temp_off_tib <- safe_read( file = opt$temp_off_csv, 
                           use_spec = TRUE, 
                           spec_prefernce = "rds", 
                           vb = opt$verbose, vt = 9, tt = tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         0.1.1 Load Pre-Defined::
#                                imGenomes
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

work_list <- NULL
work_list <- 
  file_list( path = opt$workflow,
             dir_only = TRUE,
             paths_exists = FALSE,
             vb = opt$verbose )

#
# TBD:: Seems like we need to remove the last two file_list() calls below if
#. opt$ref_source != "Both"
#
gen_path <- NULL
if ( opt$ref_source == "Both" ) {
  gen_path <- 
    file_list( path = opt$ref_path, 
               subs = opt$ref_species,
               unique = FALSE,
               dir_only = TRUE, 
               ret_type = COM,
               subs_exists = FALSE,
               vb = opt$verbose )
} else {
  gen_path <- 
    file_list( path = opt$ref_path, 
               subs = opt$ref_species,
               unique = FALSE,
               dir_only = TRUE, 
               ret_type = COM,
               subs_exists = FALSE,
               vb = opt$verbose ) %>%
    file_list( subs = opt$ref_source,
               unique = FALSE,
               dir_only = TRUE, 
               ret_type = COM,
               subs_exists = FALSE,
               vb = opt$verbose ) %>%
    file_list( subs = opt$ref_build,
               unique = FALSE,
               dir_only = TRUE, 
               ret_type = COM,
               subs_exists = FALSE,
               paths_exists = FALSE,
               vb = opt$verbose )
}

dir_add_str <- NULL
if ( opt$ref_source == "NCBI") dir_add_str <- "Fasta"

ref_seqs <- NULL
ref_seqs <-
  file_list( path = gen_path, 
             subs = paste0("Sequence/WholeGenome",dir_add_str),
             file = opt$ref_file,
             unique = FALSE,
             suffix = c("\\.gz", "\\.fa$", "\\.genome$"),
             ret_type = "list",
             subs_exists = FALSE,
             paths_exists = FALSE,
             files_exists = FALSE,
             vb = opt$verbose )

chr_dirs <- NULL
chr_dirs <- 
  file_list( path = gen_path, 
             subs = "Sequence/Chromosomes",
             names = names(ref_seqs),
             unique = FALSE,
             dir_only = TRUE, 
             ret_type = "list",
             subs_exists = FALSE,
             paths_exists = FALSE,
             vb = opt$verbose )

opt$ref_builds <- safe_mkdir( file.path( opt$out_path, "Genome_Builds" ) )

#
# Off Target "Skip" Genomes::
#
off_keys <- c()
off_keys <- c( "GRCh38","GRCh37","GRCh36","GRCm38","galGal5","criGriChoV1","pvirPacbioDovetail", "canFam6" )
off_keys <- c( "GRCh38","GRCh37","GRCh36","GRCm38","galGal5","criGriChoV1","pvirPacbioDovetail", "canFam3" )
off_keys <- c( "GRCh38","GRCh37","GRCh36","GRCm38","galGal5","criGriChoV1","pvirPacbioDovetail" )

off_keys <- c( "GRCh38","GRCh37","GRCh36","GRCm38","galGal5","criGriChoV1","pvirPacbioDovetail", "canFam6" )
off_keys <- c( "GRCh37","GRCh36","GRCm38","galGal5","criGriChoV1","pvirPacbioDovetail", "canFam3", "canFam6" )

#
# TBD:: Allow all sources for reference genomes...
#  - Fix the Genome_Builds Fasta issue
#

# tmp_list <- readr::read_csv( file.path( opt$top_path, "tmp/ChineseHamster.chr-list.txt" ), col_names = FALSE ) %>% purrr::set_names( c("ID") )

# opt$parallel <- FALSE
# opt$reload <- 0

opt$single <- FALSE
opt$single <- TRUE

# TBD:: Determine if we should parallize this (maybe only when cluster is
#   detected???)
#
chr_csvs <- list()
if ( !is.null( work_list[["parse_genome"]] ) ) {
  if ( opt$verbose > 0 ) cat(glue::glue(
    "{pmssg} Work List = {work_list['parse_genome']}...{RET}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Parse and Prep all Genomes::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  for ( ref_tag in names(ref_seqs) ) {
    # ref_tag <- "GRCh38"
    # ref_tag <- names(ref_seqs)[1]
    # ref_seqs[[ref_tag]] <- "/Users/bbarnes/Documents/data/imGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/GRCh38.genome.fa.gz"
    # ref_seqs[[ref_tag]] <- "/Users/bbarnes/Documents/data/imGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/GRCh38.genome.fa"
    # chr_dirs[[ref_tag]] <- "/Users/bbarnes/Documents/data/imGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes"
    
    if ( ref_tag %in% off_keys ) {
      if ( opt$verbose > 0 ) cat(glue::glue(
        "{pmssg}{TAB} Skipping Reference Genome: '{ref_tag}'...{RET}"))
      next
    }
    if ( opt$verbose > 0 ) cat(glue::glue(
      "{pmssg}{TAB} Processing Reference Genome: '{ref_tag}'...{RET}"))
    
    ret_dir <- safe_mkdir( file.path( opt$ref_builds, ref_tag ) )
    
    # First we want to split genome into individual chromosomes, we'll check if
    #  they already exist (passing file time stamps) and write out a list file
    #  that contains the order in which they should be assigned cg#'s 
    #
    # [Done]:: check if chr_tib can be loaded and validated from existing values
    #   Load file if exist and check each file has .gzi/.fai and is bgzipped...
    #
    
    #
    # TBD:: Implement much faster file checking...
    #
    chr_csv <- NULL
    chr_csv <- genome_to_chrom( file = ref_seqs[[ref_tag]], 
                                # skip_vec = tmp_list$ID,
                                
                                chr_key  = "Chr_Key",
                                fas_key  = "Chr_Fas",
                                head_key = "Fasta_Header",
                                
                                build   = ref_tag,
                                source  = opt$ref_source, 
                                # validate_only = TRUE,
                                validate_only = FALSE,
                                
                                out_dir = chr_dirs[[ref_tag]],
                                pre_tag = NULL,
                                
                                # max = 1,
                                reload = opt$reload,
                                reload_min = 10,
                                parallel = opt$parallel,
                                ret_type = "file",
                                
                                vb=vb,vt=vt,tc=tc, tt=tt )
    
    chr_csvs[[ref_tag]] <- chr_csv
    
    # TBD:: Same download/chromosome partitioning needs to happen to dbSNP
    #   vcf files::
    #
    #   - all_commnon_YYYYMMDD.vcf.gz
    #        [ https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/ ]
    #        [ rsync -a -P rsync://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153Common.bb ]
    #   - clinvar_YYYYMMDD.vcf.gz 
    #        [ https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/ ]
    #        [ rsync -a -P rsync://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153ClinVar.bb ]
    #
    
    if ( p1 )
      cat(glue::glue("{pmssg}{TAB} Done. Processing Reference: '{ref_tag}'.{RET2}"))
    
    if ( opt$single ) break
    
    # if ( ref_tag == "GRCm38" ) break
    # if ( ref_tag == "galGal5" ) break
    # if ( ref_tag == "criGriChoV1" ) break
  }
  if ( p0 ) cat(glue::glue(
    "{pmssg} Done. Work List = {work_list['parse_genome']}.{RET2}") )
}

# opt$parallel <- FALSE

chr_max <- 2
chr_max <- 0
top_csvs <- list()
if ( !is.null( work_list[["parse_chromosome"]] ) ) {
  if ( opt$verbose > 0 ) cat(glue::glue(
    "{pmssg} Work List = {work_list['parse_chromosome']}...{RET}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Extract Forward and Top Sequences::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  for ( ref_tag in names(chr_csvs) ) {
    
    chr_csv <- NULL
    chr_csv <- chr_csvs[[ref_tag]]
    
    if ( opt$verbose > 0 ) cat(glue::glue(
      "{pmssg}{TAB} Processing Reference Chromosomes: '{ref_tag}'...{RET}") )
    
    ret_dir <- safe_mkdir( file.path( opt$ref_builds, ref_tag ) )
    
    top_csv <- NULL
    top_csv <- chroms_to_tops( x = chr_csv,
                               
                               di_nuc = c("CG","YG","CR"), 
                               # di_nuc = c("CG"),
                               chr_key = "Chr_Key",
                               fas_key = "Chr_Bgz",
                               top_key = "Top_Fas",
                               
                               off_tib = temp_off_tib,
                               
                               build   = ref_tag,
                               out_dir = ret_dir,
                               run_tag = opt$run_name,
                               pre_tag = chr_csv,
                               
                               max = chr_max,
                               reload = opt$reload,
                               reload_min = 1,
                               parallel = opt$parallel,
                               # parallel = FALSE,
                               ret_type = "file",
                               
                               vb = opt$verbose + 2, tt = tt )
    
    top_csvs[[ref_tag]] <- top_csv
    
    #
    # Scratch for search chr2 with [CYSMBHV][GRSKBDV] instead of CG
    #
    #  chr2_fas <- "/Users/bretbarnes/Documents/data/imGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes/chr2.fa.gz"
    #  chr2_bio <- Biostrings::readDNAStringSet(filepath = chr2_fas, format = "fasta")
    #  chr2_cpg <- Biostrings::matchPattern( pattern = "CG", subject = chr2_bio[[1]] )
    #  [dope] chr2_all <- Biostrings::matchPattern( pattern = "[CYSMBHV][GRSKBDV]", subject = chr2_bio[[1]] )
    #  chr2_yG <- Biostrings::matchPattern( pattern = "YG", subject = chr2_bio[[1]] )
    #  chr2_Cr <- Biostrings::matchPattern( pattern = "CR", subject = chr2_bio[[1]] )
    #
    # tru_top_tib %>% dplyr::mutate(sub_seq = stringr::str_sub(fwd, 61, 62)) %>% dplyr::filter(sub_seq == "YG")
    
    if ( opt$verbose > 0 )
      cat(glue::glue("{pmssg}{TAB} Done. Processing Reference: '{ref_tag}'.{RET2}"))
    
    if ( opt$single ) break
  }
  
  if ( opt$verbose > 0 ) cat(glue::glue(
    "{pmssg} Done. Work List = {work_list['parse_chromosome']}.{RET2}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
# LEFT OFF HERE::
#
#. Question:
#.  1. Is it better to use map (i.e. hash) on seq or binary?
#.     - Data: data/improbe/dbSNP_Core8.chrMap.v1/design-input/GRCh38_GRCh37_GRCh36_GRCm10_GRCr00_GRCk01_GRGa02_CFAM31.cgn-set.csv.gz
#.     - Load line by line and create map/array
#        1. Based on character seq
#.       2. Based on unit8 sequence
#
#
#
#
#                               Above Works
#
# TBD: Investigation:
#.  1. What does a 2-bit 122mer look like? Unsigned Integer?
#
# TBD:: Workflow
#.  1. Split parse_chromosome into binary TOP files and Bed Files (wtih top int, strand)
#.  2. After investigation can we simply on the fly assign uint rather than cgn?
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# 122-mer::
#  8-52[CG]52-8 = (8+52+2+52+8) = 122 (120) [ bytes = 120/8 = 15 ]
# 106-mer (104-mer):
#  8-52[CG]52-8 = (52+2+52) = 106 (104) [ bytes = 104/8 = 13 ]
#


#
# Need a wrapper function for dna_to_bit()
#
#  bit1 <- dna_to_bit( can_tib$Top_Seq, 122 )

can_tib <- NULL
can_tib <- safe_read( file = opt$canonical_csv,
                      use_spec = TRUE, 
                      spec_prefernce = "rds", 
                      vb = opt$verbose, tt = tt ) %>%
  dplyr::mutate(Can_Cgn=as.integer(Can_Cgn)) %>% 
  dplyr::arrange(Can_Cgn)

Rcpp::sizeof( can_tib$Top_Seq[1] )


#
# Skip Step:: For Dingos Design Round 1::
#

# top_csvs[["canFam3"]] <- "/Users/bretbarnes/Documents/scratch/triplecrown/SNP-Both-v3/Genome_Builds/canFam3/chroms_to_tops/chromosome-order-list.csv.gz"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Internal Loading Test of Full CGN Database::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

old_school_dat <- FALSE

if ( !old_school_dat ) {
  #
  # New Master File below, but need to determine the correct next step...
  #
  max_cgn_tsv <- file.path( opt$top_path, "data/improbe/dbSNP_Core8.chrMap.v1/design-input/GRCh38_GRCh37_GRCh36_GRCm10_GRCr00_GRCk01_GRGa02_CFAM31.cgn-set.csv.gz" )
  
} else {
  max_cgn_cols <- cols(
    Can_Cgn = col_integer(),
    Top_Seq = col_character()
  )
  
  max_cgn_tsv <- file.path( opt$top_path, "data/improbe/scratch/cgnDB/dbSNP_Core7.chrMap.v2/design-input/GRCh38_GRCh37_GRCh36_GRCm10_GRCr00_GRCk01_GRGa02.cgn-set.cgn-top.tsv.gz" )
  max_cgn_tib <- readr::read_tsv( file = max_cgn_tsv, 
                                  col_names = names(max_cgn_cols$cols), 
                                  col_types = max_cgn_cols )
}

#
# Extended version with product categories::
#
if ( FALSE ) {
  # 1,1,1,0, 0,0,0,0
  max_cgn_cols <- cols(
    Can_Cgn = col_integer(),
    Top_Seq = col_character(),
    Src_Bol1 = col_skip(),
    Src_Bol2 = col_skip(),
    Src_Bol3 = col_skip(),
    Src_Bol4 = col_skip(),
    Src_Bol5 = col_skip(),
    Src_Bol6 = col_skip(),
    Src_Bol7 = col_skip(),
    Src_Bol8 = col_skip()
  )
  
  max_cgn_csv <- file.path( opt$top_path, "data/improbe/scratch/cgnDB/dbSNP_Core7.chrMap.v2/design-input/GRCh38_GRCh37_GRCh36_GRCm10_GRCr00_GRCk01_GRGa02.cgn-set.csv.gz" )
  max_cgn_tib <- readr::read_csv( file = max_cgn_csv, 
                                  col_names = names(max_cgn_cols$cols), 
                                  col_types = max_cgn_cols )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                   END::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


chr_max <- 3
chr_max <- 2
opt$single <- TRUE
compare_true <- FALSE

if ( !is.null( work_list[["update_database"]] ) ) {
  
  if ( opt$verbose > 0 ) cat(glue::glue(
    "{pmssg} Work List = {work_list['update_database']}...{RET}") )
  
  cgn_path <- file.path( opt$out_path, "update_database")
  can_path <- file.path( cgn_path, "canonical" )
  can_path <- safe_mkdir( can_path, vb = opt$verbose )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Initialize Canoncial CG# <=> Top Sequence::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  can_tib <- NULL
  can_tib <- safe_read( file = opt$canonical_csv,
                        use_spec = TRUE, 
                        spec_prefernce = "rds", 
                        vb = opt$verbose, tt = tt ) %>%
    dplyr::mutate(Can_Cgn=as.integer(Can_Cgn)) %>% 
    dplyr::arrange(Can_Cgn)
  
  # cgn_max <- 5 # TEST VALUE
  cgn_max <- 0 # RUN VALUE
  cgn_dbs <- NULL
  cgn_dbs <- init_cgDb( top_vec = can_tib$Top_Seq,
                        cgn_vec = can_tib$Can_Cgn,
                        path = can_path,
                        name = opt$run_name,
                        all = TRUE, # all = FALSE
                        sep = "b",  # sep = "t",
                        cgn_max = cgn_max, 
                        trim = TRUE,
                        
                        vb=vb,vt=vt+10,tc=tc )
  
  #
  # Quick Tibble Conversion for easy calculations::
  #
  cgn_dbs_tib <- cgn_dbs %>% 
    tibble::as_tibble( rownames = "CGN" ) %>% 
    dplyr::mutate( CGN = CGN %>% as.integer() ) %>%
    dplyr::rename( Top_Seq = value ) 
  
  # Historic Numbers::
  #  Mapped:  1242920
  #  UnMapd: 47101282
  map_top_cnt <- cgn_dbs_tib %>% dplyr::filter( Top_Seq != "" ) %>% base::nrow()
  if ( opt$verbose > 0 ) cat(glue::glue(
    "{pmssg} Mapped Top to CGN Count={map_top_cnt}.{RET2}") )

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Truth Data Set::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # TBD:: Truth Set 2:: Use spike v2 Chicago order to validate top/cgn assignment...
  #

  if ( compare_true ) {
    tru_min_idx <- 0
    cgn_tru_map <- safe_read( file = opt$cgn_tru_csv, vb = opt$verbose )
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Load all Top Sequences and Assign CG#s::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  for ( ref_tag in names(top_csvs) ) {
    
    if ( opt$verbose > 0 )
      cat(glue::glue("{pmssg} Assigning Reference: '{ref_tag}'...{RET}"))
    
    # cgn_ref_path <- file.path( cgn_dbs_path, ref_tag )
    # safe_mkdir( cgn_ref_path, vb = opt$verbose )
    
    ret_dir <- safe_mkdir( file.path( opt$ref_builds, ref_tag ) )
    safe_mkdir( ret_dir, vb = opt$verbose )
    
    # NOTE:: The code below is for the case where there were different docker
    #   builds for runs on the cluster. We'll ignore this for now...
    #
    # if ( is.null( opt$top_path ) || !dir.exists( opt$top_path ) )
    #   opt$top_path <- opt$out_path
    # top_dir <- file.path( opt$top_path, "Genome_Builds", ref_tag )
    
    #
    #
    # TBD:: LEFT OFF HERE!!!
    #  - Use db.print() capabilities
    #  - Write compare tru-partitions to print() results for each batch...
    #
    #  - Either use new validated CGN for Chicago or look up from tru-list
    #    - improbe on top sequences
    #    - bsmap on selecteed designs against hg19
    #
    #
    top_tib <- NULL
    top_csv <- top_csvs[[ref_tag]]
    
    if ( !file.exists(top_csv) ) {
      wflag <- TRUE
      warn_mssg <- glue::glue("top_csv = '{top_csv}' does not exist...",
                              "Skipping {ref_tag}. This should not happen")
      if ( wflag ) cat(glue::glue("{pwarn} {warn_mssg}!{RET2}"))
      wflag <- FALSE
      
      # eflag <- TRUE
      # errs_mssg <- glue::glue("top_csv = '{top_csv}' does not exist...",
      #                         "Skipping {ref_tag}. This should not happen")
      # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      # if ( eflag ) return(NULL)
      
    } else {
      top_tib <- safe_read( file = top_csv, 
                            type = "csv", 
                            vb = opt$verbose, tt = tt )
      
      tib_len = top_tib %>% base::nrow()
      
      for ( tib_idx in c(1:tib_len) ) {
        
        top_fas <- top_tib$Top_Fas[tib_idx]
        top_chr <- top_tib$Chr_Key[tib_idx]
        
        cgn_chr_name <- paste(opt$run_name, ref_tag, top_chr, sep='.' )
        # cgn_chr_path <- file.path( cgn_ref_path, top_chr )
        cgn_chr_path <- file.path( ret_dir, top_chr )
        safe_mkdir( cgn_chr_path, vb = opt$verbose )
        
        top_bld_tib <- NULL
        top_bld_tib <- safe_read( file = top_fas, vb = vb )
        # top_bld_tib <- read_fast_files( paths_vec = top_fas, sep = COM, vb = 10 )
        
        #
        #
        # TBD:: Test printing after cleaning()
        #   - [Done]: cgn_init
        #   - And push_tops()
        #
        # TBD:: Test < vs. <= loop_size in push_tops()
        # TBD:: Load True Map:: opts$cgn_tru_csv
        # TBD:: Load/Define current sum.csv/idx.csv
        #    - use min_cgn from sum to calculate files to load
        #      - keep data in memory???
        #    - bgzip idx.csv???
        #
        # TBD:: Functionalize all of this to allow for pre=loading (template-2)
        # TBD:: Functionalize all of this to allow for pre=loading (template-2)
        # TBD:: Functionalize all of this to allow for pre=loading (template-2)
        # TBD:: Functionalize all of this to allow for pre=loading (template-2)
        #
        
        # chr1 length(cgn_dbs) = 96688394
        
        # OLD CALL VERSION BELOW::
        #
        # cgn_dbs = push_tops( top_vec = top_bld_tib$Top_Sequence,
        #                      cdb_vec = cgn_dbs,
        #                      
        #                      path = cgn_chr_path,
        #                      name = cgn_chr_name,
        #                      write_idx = TRUE,
        #                      write_seq = FALSE,
        #                      sep = COM,
        #                      
        #                      max = cgn_max,
        #                      clean = TRUE,
        #                      vb = opt$verbose, vt = 2 )
        
        # NEW CALL VERSION::
        #  NOTE:: Need to RcppSource: scripts/R/Rcpp/infinium_bisulfite_functions.cpp
        #
        cgn_dbs = push_tops2( top_vec = top_bld_tib$Top_Sequence,
                              cgn_vec = cgn_dbs,
                              
                              path = cgn_chr_path,
                              name = cgn_chr_name,
                              write_idx = TRUE,
                              write_seq = FALSE,
                              sep = COM,
                              
                              cgn_max = cgn_max,
                              vb=vb, vt=3 )
        
        #
        # Quick Tibble Conversion for easy calculations::
        #
        cgn_dbs_tib <- cgn_dbs %>% 
          tibble::as_tibble( rownames = "CGN" ) %>% 
          dplyr::mutate( CGN = CGN %>% as.integer() ) %>%
          dplyr::rename( Top_Seq = value ) 
        
        # Historic Numbers::
        #  Mapped:  1242920
        #  UnMapd: 47101282
        map_top_cnt <- cgn_dbs_tib %>% dplyr::filter( Top_Seq != "" ) %>% base::nrow()
        if ( opt$verbose > 0 ) cat(glue::glue(
          "{pmssg} Mapped Top to CGN Count={map_top_cnt}.{RET2}") )
        
        
        
        
        
        
        if ( compare_true ) {
          #
          # Expected output::
          #
          idx_cols = cols(
            cgn = col_integer(),
            top = col_character()
          )
          cgn_chr_sum = file.path( cgn_chr_path, paste0(cgn_chr_name,".sum.csv") )
          cgn_chr_idx = file.path( cgn_chr_path, paste0(cgn_chr_name,".idx.csv") )
          
          if ( is.null(cgn_chr_sum) || !file.exists(cgn_chr_sum) ) {
            eflag <- TRUE
            errs_mssg <- glue::glue("cgn_chr_sum = '{cgn_chr_sum}' does not exist...",
                                    "Skipping {ref_tag}/{top_chr}. This should not happen")
            if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
            if ( eflag ) next
            # if ( eflag ) return(NULL)
          }
          
          if ( is.null(cgn_chr_idx) || !file.exists(cgn_chr_idx) ) {
            eflag <- TRUE
            errs_mssg <- glue::glue("cgn_chr_idx = '{cgn_chr_idx}' does not exist...",
                                    "Skipping {ref_tag}/{top_chr}. This should not happen")
            if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
            if ( eflag ) next
            # if ( eflag ) return(NULL)
          }
          
          # Read Candidate data::
          cgn_chr_sum_tib <- safe_read( file = cgn_chr_sum, vb = opt$verbose, tt = tt )
          cgn_chr_idx_tib <- safe_read( file = cgn_chr_idx, spec_col = idx_cols, 
                                        has_head = FALSE, spec_prefernce = "col", 
                                        use_spec = TRUE, over_write = TRUE, type = "csv",
                                        vb = opt$verbose, tt = tt )
          
          tru_max_idx <- cgn_tru_map %>% dplyr::filter( cgn_chr_sum_tib$min_cgn >= Range_Beg & cgn_chr_sum_tib$min_cgn <= Range_End ) %>% head(n=1) %>% dplyr::pull(Index)
          
          cgn_tru_vec <- cgn_tru_map %>% dplyr::filter( Index > tru_min_idx & Index <= tru_max_idx ) %>% dplyr::pull(Path)
          
          cgn_tru_tib <- lapply( cgn_tru_vec, safe_read, spec_col = idx_cols, has_head = FALSE, spec_prefernce = "col", use_spec = TRUE,  vb = opt$verbose+3 ) %>% dplyr::bind_rows()
          
          
          #
          # CGN Join::
          #
          cgn_join_tib <- cgn_tru_tib %>% 
            dplyr::inner_join(cgn_chr_idx_tib, by="cgn", suffix=c("_ref","_can"))
          
          cgn_join_tib %>% filter(top_ref != top_can) %>% print()
          # cgn_join_tib %>% filter(top_ref != top_can) %>% head() %>% as.data.frame()
          
          #
          # TOP Join::
          #
          top_join_tib <- cgn_tru_tib %>% 
            dplyr::inner_join(cgn_chr_idx_tib, by="top", suffix=c("_ref","_can"))
          
          top_join_tib %>% filter(cgn_ref != cgn_can) %>% print()
          
          
          #
          # Scratch for comparison::
          #
          # top_join_tib %>% dplyr::filter( cgn_ref >= 2290514 & cgn_ref <= 2290521)
          # cgn_join_tib %>% dplyr::filter( cgn >= 2290514 & cgn <= 2290521 ) %>% as.data.frame()
          
          
          # TBD:: Should double check this now that it is fixed!!!!
          #
          #  -- Chr 1 example below passes now!!!
          # cgn_tru_tib %>% dplyr::filter( cgn >= 2290514 & cgn <= 2290521 ) %>% as.data.frame()
          # cgn_chr_idx_tib %>% dplyr::filter( cgn >= 2290514 & cgn <= 2290521 ) %>% as.data.frame()
          
          #  -- Chr 2 FAILURE example below
          # TBD:: Looks like the sorting isn't producing the correct order. Probably need to sort
          #   on TOP not Forward...
          #
          # cgn_tru_tib %>% dplyr::filter( cgn >= 2518175 & cgn <= 2518183 ) %>% as.data.frame()
          # cgn_chr_idx_tib %>% dplyr::filter( cgn >= 2518175 & cgn <= 2518183 ) %>% as.data.frame()
          
          
          #
          # Set new min index::
          #
          # tru_min_idx <- tru_max_idx
          
        } # END:: if ( compare_true )
        
        if (tib_idx >= chr_max) break
      } # END:: for ( tib_idx in c(1:tib_len) )
    }
    
    if ( opt$single ) break
  } # END:: for ( ref_tag in names(top_csvs) )
  
  if ( opt$verbose > 0 ) cat(glue::glue(
    "{pmssg} Done. Work List = {work_list['update_database']}.{RET2}") )
}











#
# Basic Validation Testing::
#
if (FALSE) {
  
  cgn_dbs_tib <- cgn_dbs %>% 
    tibble::as_tibble( rownames = "CGN" ) %>% 
    dplyr::mutate( CGN = CGN %>% as.integer() ) %>%
    dplyr::rename( Top_Seq = value ) 
  
  man_csv <- file.path( opt$top_path, "data/CustomContent/McMaster/manifest/latest-10082021/GMELMethylation1_20051285_A1.csv.gz" )
  man_tib <- safe_read( man_csv, vb=vb ) %>%
    dplyr::filter( !is.na(MAPINFO) ) %>%
    dplyr::filter( !is.na( Name ) ) %>%
    dplyr::mutate( Fwd_Sequence = shear_brac(Forward_Sequence) %>% stringr::str_to_upper(),
                   CGN = Name %>% stringr::str_remove("_.*$") %>%
                     stringr::str_remove("^[a-zA-Z]+") %>%
                     stringr::str_remove("^0+") %>%
                     as.integer(),
                   CGN = CGN + 1 ) %>%
    dplyr::select( CGN, Fwd_Sequence, dplyr::everything() ) %>% 
    dplyr::filter( !is.na(CGN) )
  
  #
  # TBD::
  #  - Understand the off one CGN issues
  #  - Re-run with everything...
  #
  
  man_top_tib <- fwd2tops_cpp( seq_vec_r = man_tib$Fwd_Sequence,
                               pos_vec_r = man_tib$MAPINFO, uc = TRUE, return_source = TRUE ) %>%
    tibble::as_tibble()
  
  cgn_map_tib <- cgn_dbs_tib %>% 
    dplyr::inner_join( man_top_tib, by=c("Top_Seq"="Top_Sequence") )
  
  cgn_map_tib %>% dplyr::inner_join(
    man_tib, by=c("Fwd_Sequence"), suffix=c("_cgn", "_man")
  ) %>% head() %>% as.data.frame()
  
  cgn_map_tib %>% dplyr::inner_join(
    man_tib, by=c("Fwd_Sequence"), suffix=c("_cgn", "_man")
  ) %>% dplyr::filter( CGN_cgn == CGN_man )
  
  
  # test_top_seq <- "CGGAGGGAGCCCCACAATGATTGAGATATTCTGAGCCAGCAGGCCCTCCC,GCGTTGGAGACAGCTGAGAGGCGGTTGATAAAATCTAATTGCCCCATCGATCCAGCAGAG[CG]GAGGGAGCCCCACAATGATTGAGATATTCTGAGCCAGCAGGCCCTCCCCTGTGCCTTCAC"
  test_top_seq <- "CGGAGGGAGCCCCACAATGATTGAGATATTCTGAGCCAGCAGGCCCTCCC,GCGTTGGAGACAGCTGAGAGGCGGTTGATAAAATCTAATTGCCCCATCGATCCAGCAGAGCGGAGGGAGCCCCACAATGATTGAGATATTCTGAGCCAGCAGGCCCTCCCCTGTGCCTTCAC"
  test_cgn_num <- 582671
  
  cgn_dbs_tib %>% dplyr::filter( Top_Seq == test_top_seq )
  cgn_dbs_tib %>% dplyr::filter( CGN == test_cgn_num )
  
}

if (FALSE) {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Scratch Space for Expand_Seqs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (FALSE) {
    
    #
    # Scratch space for splitting improbe-validation database::
    #
    imp_tsv <- "/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-probe-designs.sort-unique.wHeader.tsv.gz"
    # imp_tib <- safe_read(imp_csv, type="tsv", vb = opt$verbose, tt=tt)
    
    # Even this doesn't work....
    # imp_tib <- read_fast_files( paths_vec = c(imp_tsv), 
    #                             sep = TAB, vb = 10 )
    
    # To Do::
    #  - toFile(dir, name); # Write db and summary
    #
    
    test_seqs <- c(
      "AAAYAAA",
      "AAARAAA",
      
      "ttyttRt",
      "ccrccyc",
      
      "AAsttyttRt",
      "AAmccrccyc",
      
      "HAAsttyttRt",
      "DAAmccrccyc",
      
      'zzzzzzz'
    )
    
    expand_seqs_cpp(test_seqs)
    
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

# End of file

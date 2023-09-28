
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   Scratch Space for Docker Workhorse::
#                      Finger Printing Ivestigation
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
suppressWarnings(suppressPackageStartupMessages( 
  base::require("sesame",  quietly = TRUE) ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Set Run Environment:: RStudio/Command-Line
#                            Source All Functions
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par  <- list()
args <- commandArgs(trailingOnly = FALSE)

par$src_path <- NULL
par$run_mode <- args[1]
par$date_str <- Sys.Date() %>% as.character()
par$prgm_dir <- 'stable'
par$prgm_tag <- 'stable_analysis_EPICv2_docker_snps'
par$verbose  <- 3
local_paths  <- c(
  "/Users/bbarnes/Documents/tools/imSuite/scripts/R",
  # "/Users/bbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
  # "/Users/bretbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
  # "/illumina/scratch/darkmatter/tools/Workhorse-Unstained/scripts/R",
  "/repo/Workhorse-Unstained/scripts/R" )

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
par <- source_functions( pars = par, vb = par$verbose )
par <- params_check( pars = par, args = args, 
                     prgm_aux_check = FALSE, vb = par$verbose )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$run_name <- "EPICv2"

# Multiple copies of v0 (bk and bk24)
par$version <- 0
# par$version <- 1
# par$version <- 2
# par$version <- 3
# par$version <- 4
# par$version <- 5

opt <- NULL
# opt <- swifthoof_sesame_options( pars = par, args = args, vb = par$verbose )
opt <- stable_options( pars = par, args = args, vb = par$verbose )
vb  <- opt$verbose
vt  <- 0
tc  <- 0

p0  <- vb > vt + 0
p1  <- vb > vt + 1
p2  <- vb > vt + 2
p4  <- vb > vt + 4
p8  <- vb > vt + 8

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c( 'run_mode',
               'src_path', 'scr_path', 'exe_path', 'prgm_dir', 'prgm_tag' )
opt_reqs <- c( 'out_path', # 'ref_path', 'ref_file', 'ref_build', 'ref_species',
               'Rscript', 'verbose' )

opt$rcpp <- 2
opt$rcpp <- 0
#
# TBD:: Return Rcpp value instead of boolean::
#
prgm_dat <- program_init( name = par$prgm_tag,
                          opts = opt, opt_reqs = opt_reqs,
                          pars = par, par_reqs = par_reqs, 
                          rcpp = opt$rcpp,
                          vb=vb, vt=3, tc=0 )

opt <- prgm_dat$opt
par <- prgm_dat$par
opt_tib <- prgm_dat$opt_tib
par_tib <- prgm_dat$par_tib

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Initialize Run Objects
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tt <- timeTracker$new()

pmssg <- glue::glue("[{par$prgm_tag}]:")
pwarn <- glue::glue("{RET}[{par$prgm_tag}]: Warning:")
perrs <- glue::glue("{RET}[{par$prgm_tag}]: ERROR:")

success = TRUE

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Pre-processing::
#
#                        Load Truth Sample Sheets::
#                         Load Auto Sample Sheets::
#                           Load Calls:: raw/ind
#                             Merging VCFs::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

v1_man_tib <- NULL
v1_man_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/manifest/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz" )
v1_man_tib <- readr::read_csv( file = v1_man_csv, show_col_types = FALSE ) %>%
  dplyr::mutate( Full_ID = Probe_ID )

v2_man_tib <- NULL
v2_man_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v4/docker/man/EPICv2/EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
v2_man_tib <- readr::read_csv( file = v2_man_csv, show_col_types = FALSE )

v3_man_tib <- NULL
v3_man_tib <- dplyr::bind_rows(
  v1_man_tib %>% dplyr::filter( !Probe_ID %in% v2_man_tib$Probe_ID ) %>% dplyr::mutate( BP_Group = "EPICv1", Probe_Group_BP = "v1" ),
  v2_man_tib %>% dplyr::filter(  Probe_ID %in% v1_man_tib$Probe_ID ) %>% dplyr::mutate( BP_Group = "EPICv3", Probe_Group_BP = "v3" ),
  # v1_man_tib %>% dplyr::filter(  Probe_ID %in% v2_man_tib$Probe_ID ) %>% dplyr::mutate( BP_Group = "EPICv3", Probe_Group_BP = "v3" ),
  v2_man_tib %>% dplyr::filter( !Probe_ID %in% v1_man_tib$Probe_ID ) %>% dplyr::mutate( BP_Group = "EPICv2", Probe_Group_BP = "v2" )
) %>% 
  dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>%
  dplyr::mutate( Infinium = DESIGN,
                 Probe_Group_Mask = FALSE ) %>%
  dplyr::select( Probe_ID, Probe_Type, Infinium, Probe_Group_Mask,
                 BP_Group, Probe_Group_BP, dplyr::everything() )

# [TBD]:
#  - Truth Comparisons Sample Calling
#  - Comparisons Beta (r2/dB) & SNVs (gt) x EPICv1/v2
#

opt$single <- TRUE
opt$single <- FALSE

if ( opt$single ) {
  v1_dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/EPICv1/chip-206203800149" )
  v2_dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/EPICv2/chip-206891110001" )
} else {
  v1_dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/EPICv1" )
  v2_dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/EPICv2" )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Load Truth Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

true_ssh_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/formatted/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz" )
true_ssh_tib <- NULL
true_ssh_tib <- readr::read_csv( file = true_ssh_csv, show_col_types = FALSE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Load Auto Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# [TBD]: Rename Titration Samples with values and same name...
#
v2_ssh_list <- NULL
v2_ssh_list <- file_list( path    = v2_dat_path, 
                          prefix  = v2_dat_path,
                          suffix  = "_EPIC_A1_AutoSampleSheet.csv.gz", 
                          pattern = "_EPIC_A1_AutoSampleSheet.csv.gz$",
                          recursive = TRUE )

v2_ssh_tib <- NULL
v2_ssh_tib <- v2_ssh_list %>% 
  lapply( readr::read_csv, show_col_types = FALSE ) %>% 
  dplyr::bind_rows() %>%
  dplyr::mutate( cg_calls_pass_perc_1 )

v2_sam_tib <- NULL
v2_sam_tib <- dplyr::right_join(
  true_ssh_tib %>% dplyr::select( Sentrix_Name, Sample_Base ) %>%
    dplyr::mutate( Sample_Base = Sample_Base %>% stringr::str_to_upper() ),
  v2_ssh_tib %>% 
    dplyr::select( Sentrix_Name, 
                   AutoSample_R2_Key_1,AutoSample_dB_Key_1,
                   AutoSample_R2_Val_1,AutoSample_dB_Val_1, 
                   AutoSample_dB_Cnt_1 ),
  by=c("Sentrix_Name")
)

# v2_sam_tib %>% dplyr::filter( Sample_Base == AutoSample_R2_Key_1 )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Calls:: raw
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# [TBD]: Functionalize all of this so we can pass a list of files to generate 
#. a set of plots...
#

v2_raw_list <- NULL
v2_raw_list <- file_list( path    = v2_dat_path, 
                          prefix  = v2_dat_path,
                          suffix  = "_EPIC_A1_raw.call.dat.csv.gz", 
                          pattern = "_EPIC_A1_raw.call.dat.csv.gz$",
                          recursive = TRUE )

v2_raw_tab <- NULL
v2_raw_tab <- v2_raw_list %>%
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows( .id = "Sentrix_Name" )

v2_raw_key_list <- NULL
v2_raw_key_list <- v2_raw_tab %>%
  tidyr::pivot_longer( cols = c( pvals_pOOBAH,pvals_PnegEcdf,betas ), 
                       names_to  = c("Key"), 
                       values_to = c("Val") ) %>% split(.$Key)

# This can now be passed to use consolidate_matricies_wrapper()
v2_raw_tab <- NULL
v2_raw_tab <- v2_raw_key_list %>% 
  lapply( function(x) { 
    x %>% dplyr::select(-Key) %>% 
      tidyr::pivot_wider( id_cols = c(Probe_ID), 
                          names_from = c(Sentrix_Name), 
                          values_from = c(Val) )
  })

#
# Convert to Matricies::
#
v2_raw_beta_mat <- NULL
v2_raw_beta_mat <- v2_raw_tab$betas %>% 
  tibble::column_to_rownames( var = "Probe_ID") %>%
  as.data.frame() %>%
  as.matrix()

v2_raw_poob_mat <- NULL
v2_raw_poob_mat <- v2_raw_tab$pvals_pOOBAH %>% 
  tibble::column_to_rownames( var = "Probe_ID") %>%
  as.data.frame() %>%
  as.matrix()

v2_raw_negs_mat <- NULL
v2_raw_negs_mat <- v2_raw_tab$pvals_PnegEcdf %>% 
  tibble::column_to_rownames( var = "Probe_ID") %>%
  as.data.frame() %>%
  as.matrix()

cond_poob_tib <- NULL
cond_poob_tib <- consolidate_matricies(
  betas = v2_raw_beta_mat,
  detps = v2_raw_poob_mat,
  clean = NULL,
  
  sum_stats = TRUE,
  all_stats = TRUE,
  
  vb=vb,vt=vt,tc=tc, tt=tt )

#
# Plot by Sample Name::
#
v2_sam_sum <- NULL
v2_sam_sum <- v2_sam_tib %>% 
  dplyr::group_by( Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) v2_sam_sum %>% print( n=base::nrow(v2_sam_sum) )

plot_sample_vec <- NULL
plot_sample_vec <- v2_sam_tib$Sample_Base %>% unique() %>% as.vector()
plot_sample_vec <- c("NA12873","HELA","JURKAT","MCF7","RAJI")

analysis_val <- "CellLine"
samp_name <- "swap"
detp_keys <- c("poob","negs")
beta_keys <- c("ind","raw")
prod_keys <- c("EPICv2")
pval_mins <- c(1.0, 0.05 )
# pval_mins <- c(1.0, 0.05, 0.01 )
dB_mins   <- c( 0.2 )
# dB_mins   <- c( 0.2, 0.1 )

# [TBD]: Investigate idn instead of ind...

# [TBD]: To use this correctly we need to have EPICv1/EPICv2 together...
# [TBD]: Load all samples in a single pass instead of v1/v2

v3_man_sum <- NULL
v3_man_sum <- v3_man_tib %>% 
  dplyr::group_by( Probe_Group_Mask,BP_Group,Probe_Group_BP ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) v3_man_sum %>% print( n=base::nrow(v3_man_sum) )

for ( detp_key in detp_keys ) {
  for ( beta_key in beta_keys ) {
    for ( prod_key in prod_keys ) {
      for ( pval_min in pval_mins ) {
        for ( dB_min in dB_mins ) {
          cat(glue::glue("{pmssg} Params: {detp_key}, {beta_key}, {prod_key}, {pval_min}, {dB_min}...{RET}"))
          
          
          cat(glue::glue("{pmssg} Done.{RET2}"))
        }
      }
    }
  }
}
for ( sample_base in plot_sample_vec ) {
  
  # [TBD]: Loop over these variables::
  detp_key <- detp_keys[1]
  beta_key <- beta_keys[1]
  prod_key <- prod_keys[1]
  pval_min <- pval_mins[2]
  dB_min   <- dB_mins[1]
  
  # Need to extract the correct Sentrix_Names
  sentrix_vec <- NULL
  sentrix_vec <- v2_sam_tib %>% dplyr::filter( Sample_Base == sample_base ) %>% dplyr::pull( Sentrix_Name )
  
  sample_red <- 0
  sample_red <- 3
  betas_ncol <- base::ncol( v2_raw_tab$betas %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) )
  detps_ncol <- base::ncol( v2_raw_tab$pvals_pOOBAH %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) )
  
  betas_ncol <- base::ncol( v2_raw_tab$betas %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) ) - sample_red
  detps_ncol <- base::ncol( v2_raw_tab$pvals_pOOBAH %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) ) - sample_red
  
  
  plot_ret <- NULL
  plot_ret <- plot_beta_gg( betas = v2_raw_tab$betas %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) %>% dplyr::select( dplyr::all_of( 1:betas_ncol) ),
                            detps = v2_raw_tab$pvals_pOOBAH %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) %>% dplyr::select( dplyr::all_of( 1:detps_ncol) ),
                            manifest = v3_man_tib, 
                            
                            ids_key = "Probe_ID", prb_key = "Probe_Type", inf_key = "Infinium",
                            
                            dB_min   = dB_min,
                            pval_min = pval_min,
                            sub_per  = 1,
                            
                            top_tag = paste0( "Sample=",sample_base ),
                            sub_tag = paste0( glue::glue("Workflow=raw, pval={pval_min}, dB={dB_min}") ),
                            par_tag = paste( analysis_val,samp_name,
                                             detp_key,beta_key,prod_key,sample_base, sep="-" ),
                            
                            dpi_val = 320,
                            
                            out_dir = file.path( opt$out_path ),
                            run_tag = paste0( opt$run_name ),
                            
                            reload     = opt$reload,
                            reload_min = 10,
                            reload_pre = NULL,
                            
                            # ret_data   = FALSE,
                            ret_data   = TRUE,
                            parallel   = FALSE,
                            
                            vb=vb,vt=vt,tc=tc, tt=tt )
  
  # Old Consolidation Method::
  #
  if ( FALSE ) {
    plot_poob_ret <- NULL
    plot_poob_ret <- consolidate_matricies_wrapper( 
      betas = v2_raw_beta_mat[ ,sentrix_vec ],
      detps = v2_raw_poob_mat[ ,sentrix_vec ],
      manifest = v3_man_tib,
      
      pval_min = pval_min,
      dB_min   = dB_min,
      plot_con = TRUE,
      
      sam_rgq = NULL,
      top_tag = paste0( "Sample=",sample_base ),
      sub_tag = paste0( glue::glue("Workflow=raw, pval={pval_min}, dB={dB_min}") ),
      par_tag = paste( analysis_val,samp_name,
                       detp_key,beta_key,prod_key,sample_base, sep="-" ),
      
      dpi_val = 320,
      
      out_dir = file.path( opt$out_path ),
      run_tag = paste0( opt$run_name ),
      
      reload     = opt$reload,
      reload_min = 10,
      reload_pre = NULL,
      
      ret_data   = FALSE,
      parallel   = FALSE,
      
      vb=vb,vt=vt,tc=tc, tt=tt )
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Calls:: ind
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

v2_ind_list <- NULL
v2_ind_list <- file_list( path    = v2_dat_path, 
                          prefix  = v2_dat_path,
                          suffix  = "_ind.call.dat.csv.gz", 
                          pattern = "_ind.call.dat.csv.gz$",
                          recursive = TRUE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Merging VCFs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

##fileformat=VCFv4.0
##fileDate=20230710
##reference=hg19
##INFO=<ID=PVF,Number=1,Type=Float,Description="Pseudo Variant Frequency">
##INFO=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=GS,Number=1,Type=Integer,Description="Genotyping score from 7 to 85">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

# Bret's Suggestions::
#
##INFO=<ID=PVF,Number=1,Type=Integer,Description="Pseudo Variant Frequency * 100">
##INFO=<ID=RGT,Number=1,Type=Integer,Description="Reference Genotype">
##INFO=<ID=AGT,Number=1,Type=Integer,Description="Allele Genotype">
##INFO=<ID=GTS,Number=1,Type=Integer,Description="Genotyping score from 7 to 85">


vcf_cols <- NULL
vcf_cols <- readr::cols(
  Chromosome = readr::col_character(),
  Coordinate = readr::col_integer(),
  SNP_ID     = readr::col_character(),
  REF        = readr::col_character(),
  ALT        = readr::col_character(),
  Qual       = readr::col_character(),
  Filter     = readr::col_character(),
  Info_Str   = readr::col_character()
)

v2_vcf_list <- NULL
v2_vcf_list <- file_list( path    = v2_dat_path, 
                          prefix  = v2_dat_path,
                          suffix  = "_EPIC_A1_raw.snps.vcf", 
                          pattern = "_EPIC_A1_raw.snps.vcf$",
                          recursive = TRUE )

v2_snp_tab <- NULL
v2_snp_tab <- v2_vcf_list %>% 
  lapply( readr::read_tsv, 
          col_names=names(vcf_cols$cols), 
          col_types=vcf_cols, skip=6 ) %>% 
  dplyr::bind_rows( .id = "Sentrix_Name" ) %>%
  tidyr::separate( Info_Str, into = c("PVF","GT","GS"), sep=";", remove = TRUE ) %>%
  tidyr::unite( Target_ID, Chromosome,Coordinate,SNP_ID,REF,ALT, 
                sep=";", remove = TRUE) %>%
  dplyr::mutate(
    Filter = dplyr::case_when(
      Filter == "PASS" ~ 1.0,
      Filter == "FAIL" ~ 0.0,
      TRUE ~ NA_real_ ) %>% as.integer(),
    PVF = PVF %>% stringr::str_remove("^PVF=") %>% as.double(),
    PVF = as.integer( PVF * 1000 ),
    GT  = GT  %>% stringr::str_remove( "^GT="),
    GTS = GS  %>% stringr::str_remove( "^GS=") %>% as.integer(),
    RGT = GT  %>% stringr::str_remove("/.*$") %>% as.integer(),
    AGT = GT  %>% stringr::str_remove("^.*/") %>% as.integer()
  ) %>% dplyr::select( -GT ) %>% clean_tib()

# Pivot Wider::
# v2_snp_wid_tib <- NULL
# v2_snp_wid_tib <- v2_snp_tab %>% 
#   tidyr::pivot_wider( id_cols = c(Target_ID), 
#                       names_from = c(Sentrix_Name), 
#                       values_from = c( Qual,Filter,PVF,GTS,RGT,AGT ) )
# 
# v2_snp_wid_tib %>% head() %>% as.data.frame()

# Piviot Longer + Split::
v2_snp_log_tib <- NULL
v2_snp_log_tib <- v2_snp_tab %>%
  tidyr::pivot_longer( cols = c( Qual,Filter,PVF,GTS,RGT,AGT ), 
                       names_to  = c("Key"), 
                       values_to = c("Val") )


# cols=c(Target_ID,Sentrix_Name)
v2_key_list <- NULL
v2_key_list <- v2_snp_tab %>%
  tidyr::pivot_longer( cols = c( Qual,Filter,PVF,GTS,RGT,AGT ), 
                       names_to  = c("Key"), 
                       values_to = c("Val") ) %>% split(.$Key)

vet_mat_list <- NULL
vet_mat_list <- v2_key_list %>% 
  lapply( function(x) { 
    x %>% dplyr::select(-Key) %>% 
      tidyr::pivot_wider( id_cols = c(Target_ID), 
                          names_from = c(Sentrix_Name), 
                          values_from = c(Val) )
  })

#
# [Done]: Implement coding below:
##INFO=<ID=PVF,Number=1,Type=Integer,Description="Pseudo Variant Frequency * 100">
##INFO=<ID=RGT,Number=1,Type=Integer,Description="Reference Genotype">
##INFO=<ID=AGT,Number=1,Type=Integer,Description="Allele Genotype">
##INFO=<ID=GTS,Number=1,Type=Integer,Description="Genotyping score from 7 to 85">
#
# [TBD]: Make weighted score (for plotting): [-1,1] * QUAL/PVF/GS | FILTER
#
# vet_mat_list$GTA %>% tibble::column_to_rownames( var = "Target_ID" ) %>% as.matrix() %>% dim()

# Need to look up plotting script...
#  Something like: plotBetasByType()
#  consolidate_matricies_wrapper()
#  beta_panel_gg()



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                             INVESTIGATION::
#
#                        dbSNP 151 Sesame Files::
#                                 Anno S/I
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  
  annoS_rds <- NULL
  annoS_rds <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoS.rds" )
  annoI_rds <- NULL
  annoI_rds <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoI.rds" )
  
  annoS_dat <- NULL
  annoS_dat <- readr::read_rds( file = annoS_rds )
  annoI_dat <- NULL
  annoI_dat <- readr::read_rds( file = annoI_rds )
  
  annoS_tib <- NULL
  annoS_tib <- annoS_dat %>% as.data.frame() %>% 
    tibble::rownames_to_column( var = "Loci_ID" ) %>% 
    tibble::as_tibble() %>%
    dplyr::rename(
      Chromosome_CpG_hg19 = seqnames,
      Beg_CpG_hg19 = start,
      End_CpG_hg19 = end,
      SNP_ID = rs
    ) %>% 
    dplyr::select( Chromosome_CpG_hg19,strand,
                   Loci_ID,SNP_ID, 
                   designType,U,REF,ALT ) %>%
    clean_tib()
  
  annoI_tib <- NULL
  annoI_tib <- annoI_dat %>% as.data.frame() %>% 
    tibble::rownames_to_column( var = "Loci_ID" ) %>% 
    tibble::as_tibble() %>%
    dplyr::rename(
      Chromosome_CpG_hg19 = seqnames,
      Beg_CpG_hg19 = start,
      End_CpG_hg19 = end,
      SNP_ID = rs
    ) %>% 
    dplyr::select( Chromosome_CpG_hg19,strand,
                   Loci_ID,SNP_ID, 
                   designType,In.band,REF,ALT ) %>%
    clean_tib()

}

if ( FALSE ) {
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Check Sesame New Function::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  # sesameData::sesameData_annoProbes()
  
  anno_grab_dat <- NULL
  anno_grab_dat <- sesameData::sesameData_annoProbes( 
    Probe_IDs = annoI_tib$Loci_ID, 
    collapse = FALSE, 
    chooseOne = FALSE, 
    platform = "EPIC" )
  
  # [TBD]:: Ask Wanding about this method...
    
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file

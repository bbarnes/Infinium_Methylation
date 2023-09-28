
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

# install.packages("VennDiagram")
suppressWarnings(suppressPackageStartupMessages( 
  base::require("VennDiagram",  quietly = TRUE) ) )
# library("VennDiagram")

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
par$prgm_tag <- 'stable_analysis_EPICv2_docker_Cancer151'
par$verbose  <- 3
local_paths  <- c(
  "/Users/bbarnes/Documents/tools/imSuite/scripts/R",
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
par$version <- 1
par$version <- 2
par$version <- 3
par$version <- 4
par$version <- 5
par$version <- 6
par$version <- 7
par$version <- 8
par$version <- 9
par$version <- 10
par$version <- 11
par$version <- 12 # Before 151 spike in trifect SNPs
par$version <- 13 # After 151 spike in trifect SNPs

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
#                  Pre-processing:: Load Truth Sample Sheets
#
# 1. [EPICv1]: Load Truth Sample Sheet
# 2. [EPICv2_MV]: Load Truth Sample Sheet
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# 151 Probes:
#
add_151_tib <- NULL
add_151_csv <- file.path( "/Users/bbarnes/mylocal/Documents/Projects.new/EPIC_v2/docker/dat/VA_SNP_Selection/Add_151_Trifecta_SNV_Probes.csv" )
add_151_tib <- readr::read_csv( file = add_151_csv, show_col_types = FALSE )

true_v1_top_dir <- file.path( opt$top_path, "Projects/EPIC_v2/GSIBIOINFO-638" )
true_v2_top_dir <- file.path( opt$top_path, "Projects/EPIC_v2/docker/dat/EPICv2_VA_MVP-Data-2023_09_05" )

# EPICv2_MSA Truth::
#
true_v2_ssh_tib <- NULL
true_v2_ssh_csv <- file.path( true_v2_top_dir, "SampleSheet-20230829-EPICv2_VA_MVP-SampleRun.csv.gz" )
true_v2_ssh_tib <- readr::read_csv( file = true_v2_ssh_csv, show_col_types = FALSE, skip = 7 ) %>%
  clean_tib() %>%
  dplyr::mutate(
    Sentrix_Name = paste( Sentrix_ID,Sentrix_Position, sep="_"),
    Sample_Base = Sample_Name %>% stringr::str_remove("_.*$") %>% stringr::str_to_upper(),
    Sample_ID = Sample_Base %>% stringr::str_sub(1,1),
    Sample_Group = Sample_Group %>% stringr::str_remove(" "),
    Sample_Class = dplyr::case_when(
      Sample_Group == "CellLine" ~ "r2",
      Sample_Group == "Coriell"  ~ "r2",
      Sample_Group == "MeTritration" ~ "me",
      Sample_Group == "NegativeMouse" ~ "nn",
      TRUE ~ NA_character_
    ),
    Concentration = 250 %>% as.integer(),
    Sample_Source = "EPICv2_MVP"
  ) %>% dplyr::select( Sentrix_Name,Sample_Base,Sample_ID,Sample_Class,
                       Concentration,Sample_Group,Sample_Source )

# EPICv1 Truth::
#
true_v1_ssh_tib <- NULL
true_v1_ssh_csv <- file.path( true_v1_top_dir, "SampleSheets/formatted/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz" )
true_v1_ssh_tib <- readr::read_csv( file = true_v1_ssh_csv, show_col_types = FALSE ) %>%
  clean_tib() %>%
  dplyr::mutate(
    Sample_Base = Sample_Base %>% stringr::str_to_upper(),
    Sample_ID = Sample_Base %>% stringr::str_sub(1,1),
    Sample_Source = "Alpha",
    Sample_Class = dplyr::case_when(
      Sample_Group == "CellLine" ~ "r2",
      Sample_Group == "Coriell"  ~ "r2",
      Sample_Group == "MeTritration" ~ "me",
      Sample_Group == "NegativeMouse" ~ "nn",
      TRUE ~ NA_character_
    )
  ) %>% dplyr::select( Sentrix_Name,Sample_Base,Sample_ID,Sample_Class,
                       Concentration,Sample_Group,Sample_Source ) %>% 
  dplyr::filter( Sample_Base %in% true_v2_ssh_tib$Sample_Base )

#
# Reduce to only HELA/RAJI::
#
true_v2_ssh_tib <- true_v2_ssh_tib %>% 
  dplyr::filter( Sample_Base %in% true_v1_ssh_tib$Sample_Base )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  Pre-processing:: Load Auto Sample Sheets
#
# 1. [EPICv1]: Load Auto Sample Sheet
# 2. [EPICv2_MV]: Load Auto Sample Sheet
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$min_call_rate <- 75

data_v1_top_dir <- file.path( opt$top_path, "Projects/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/EPICv1" )
data_v2_top_dir <- file.path( opt$top_path, "Projects/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2.1/EPICv2_MVP" )

auto_v1_ssh_tib <- NULL
auto_v1_ssh_tib <- file_list( path    = data_v1_top_dir, 
                              prefix  = data_v1_top_dir, 
                              suffix  = "_EPIC_B4_AutoSampleSheet.csv.gz", 
                              pattern = "_EPIC_B4_AutoSampleSheet.csv.gz$",
                              recursive = TRUE ) %>%
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows() %>%
  dplyr::select( Sentrix_Name, 
                 detect_version,
                 cg_calls_pass_perc_1,
                 AutoSample_R2_Key_1,AutoSample_dB_Key_1,
                 AutoSample_R2_Val_1,AutoSample_dB_Val_1, 
                 AutoSample_dB_Cnt_1 ) %>%
  dplyr::inner_join( true_v1_ssh_tib, by=c("Sentrix_Name") ) %>%
  dplyr::mutate(
    Sample_Requeue = dplyr::case_when( 
      cg_calls_pass_perc_1 >= opt$min_call_rate ~ FALSE, 
      cg_calls_pass_perc_1 <  opt$min_call_rate ~ TRUE,
      TRUE ~ NA )
  ) %>%
  dplyr::group_by( Sample_Base,detect_version ) %>% 
  dplyr::mutate(
    Rep_Num = dplyr::row_number(),
    PPP_Int = cg_calls_pass_perc_1 %>% as.integer(),
    Plot_ID = paste( detect_version,Sample_ID,Concentration,Rep_Num,PPP_Int, sep="_")
  ) %>% dplyr::ungroup() %>%
  dplyr::select( Sentrix_Name,Plot_ID, dplyr::everything() )

auto_v2_ssh_tib <- NULL
auto_v2_ssh_tib <- file_list( path    = data_v2_top_dir, 
                              prefix  = data_v2_top_dir, 
                              suffix  = "_EPIC_A1_AutoSampleSheet.csv.gz", 
                              pattern = "_EPIC_A1_AutoSampleSheet.csv.gz$",
                              recursive = TRUE ) %>%
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows() %>%
  dplyr::select( Sentrix_Name, 
                 detect_version,
                 cg_calls_pass_perc_1,
                 AutoSample_R2_Key_1,AutoSample_dB_Key_1,
                 AutoSample_R2_Val_1,AutoSample_dB_Val_1, 
                 AutoSample_dB_Cnt_1 ) %>%
  dplyr::inner_join( true_v2_ssh_tib, by=c("Sentrix_Name") ) %>%
  dplyr::mutate(
    Sample_Requeue = dplyr::case_when( 
      cg_calls_pass_perc_1 >= opt$min_call_rate ~ FALSE, 
      cg_calls_pass_perc_1 <  opt$min_call_rate ~ TRUE,
      TRUE ~ NA )
  ) %>%
  dplyr::group_by( Sample_Base,detect_version ) %>% 
  dplyr::mutate(
    Rep_Num = dplyr::row_number(),
    PPP_Int = cg_calls_pass_perc_1 %>% as.integer(),
    Plot_ID = paste( detect_version,Sample_ID,Concentration,Rep_Num,PPP_Int, sep="_")
  ) %>% dplyr::ungroup() %>%
  dplyr::select( Sentrix_Name,Plot_ID, dplyr::everything() )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Pre-processing:: Load Sample VCFs
#
# 1. [EPICv1]: Load Sample VCFs
# 2. [EPICv2_MV]: Load Sample VCFs
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$gts_min <- 20

# gtc_vec <- NULL
# gtc_vec <- c( 0, 20 )

data_v1_vcf_list <- NULL
data_v1_vcf_list <- file_list( path    = data_v1_top_dir, 
                               prefix  = data_v1_top_dir, 
                               suffix  = "_EPIC_.*_raw.snps.vcf", 
                               pattern = "_raw.snps.vcf$",
                               recursive = TRUE )

data_v2_vcf_list <- NULL
data_v2_vcf_list <- file_list( path    = data_v2_top_dir, 
                               prefix  = data_v2_top_dir,
                               suffix  = "_EPIC_.*_raw.snps.vcf", 
                               pattern = "_raw.snps.vcf$",
                               recursive = TRUE )

auto_v3_ssh_tib <- NULL
auto_v3_ssh_tib <- dplyr::bind_rows( auto_v1_ssh_tib,auto_v2_ssh_tib )

snvs_v3_gtc <- NULL
snvs_v3_gtc <- c( data_v1_vcf_list[auto_v1_ssh_tib$Sentrix_Name],
                  data_v2_vcf_list[auto_v2_ssh_tib$Sentrix_Name] ) %>%
  analyze_snvs( run_tag="All",
                ssh_tib = auto_v3_ssh_tib,
                sub_vec = NULL, # auto_v1_ssh_tib$Sentrix_Name,
                gts_min = opt$gts_min, 
                fval = 5, 
                jcnt = 0, 
                uval = TRUE,
                outDir = file.path( opt$out_path ),
                write_out = FALSE, 
                plot_heat = FALSE )

snvs_v1_gtc <- NULL
snvs_v1_gtc <- data_v1_vcf_list[auto_v1_ssh_tib$Sentrix_Name] %>% # head() %>%
  analyze_snvs( run_tag="EPICv1",
                ssh_tib = auto_v1_ssh_tib,
                sub_vec = NULL, # auto_v1_ssh_tib$Sentrix_Name,
                gts_min = opt$gts_min, 
                fval = 5, 
                jcnt = 0, 
                uval = TRUE,
                outDir = file.path( opt$out_path ),
                write_out = FALSE, 
                plot_heat = FALSE )
snvs_v1_gtc_tib <- snvs_v1_gtc %>% tibble::as_tibble()

snvs_v2_gtc <- NULL
snvs_v2_gtc <- data_v2_vcf_list[auto_v2_ssh_tib$Sentrix_Name] %>% # head() %>%
  analyze_snvs( run_tag="EPICv2_MVP",
                ssh_tib = auto_v2_ssh_tib,
                sub_vec = NULL, # auto_v1_ssh_tib$Sentrix_Name,
                gts_min = opt$gts_min, 
                fval = 5, 
                jcnt = 0, 
                uval = TRUE,
                outDir = file.path( opt$out_path ),
                write_out = FALSE, 
                plot_heat = FALSE )
snvs_v2_gtc_tib <- snvs_v2_gtc %>% tibble::as_tibble()

#
# Calculate performance
#
per_snvs_tib <- NULL

samp_v1_ssh_tibs<- auto_v1_ssh_tib %>% split(.$Sample_Base )
for ( sample in names(samp_v1_ssh_tibs) ) {
  cur_ssh <- samp_v1_ssh_tibs[[sample]]
  cut_mat <- snvs_v1_gtc[ , cur_ssh$Plot_ID ]
  cur_per <- gtc_snv_to_performance( cut_mat ) %>%
    dplyr::mutate( Match_Per = round( 100 * Match_Cnt/Total_Cnt, 3 ),
                   Added_Probe = dplyr::case_when(
                     Probe_ID %in% add_151_tib$IlmnID ~ TRUE,
                     TRUE ~ FALSE ),
                   Sample = sample,
                   Platform = "v1" )
  per_snvs_tib <- per_snvs_tib %>% dplyr::bind_rows( cur_per )
}

samp_v2_ssh_tibs<- auto_v2_ssh_tib %>% split(.$Sample_Base )
for ( sample in names(samp_v2_ssh_tibs) ) {
  cur_ssh <- samp_v2_ssh_tibs[[sample]]
  cut_mat <- snvs_v2_gtc[ , cur_ssh$Plot_ID ]
  cur_per <- gtc_snv_to_performance( cut_mat ) %>%
    dplyr::mutate( Match_Per = round( 100 * Match_Cnt/Total_Cnt, 3 ),
                   Added_Probe = dplyr::case_when(
                     Probe_ID %in% add_151_tib$IlmnID ~ TRUE,
                     TRUE ~ FALSE ),
                   Sample = sample,
                   Platform = "v2" )
  per_snvs_tib <- per_snvs_tib %>% dplyr::bind_rows( cur_per )
}

samp_v3_ssh_tibs<- auto_v3_ssh_tib %>% split(.$Sample_Base )
for ( sample in names(samp_v3_ssh_tibs) ) {
  cur_ssh <- samp_v3_ssh_tibs[[sample]]
  cut_mat <- snvs_v3_gtc[ , cur_ssh$Plot_ID ]
  cur_per <- gtc_snv_to_performance( cut_mat ) %>%
    dplyr::mutate( Match_Per = round( 100 * Match_Cnt/Total_Cnt, 3 ),
                   Added_Probe = dplyr::case_when(
                     Probe_ID %in% add_151_tib$IlmnID ~ TRUE,
                     TRUE ~ FALSE ),
                   Sample = sample,
                   Platform = "v3" )
  per_snvs_tib <- per_snvs_tib %>% dplyr::bind_rows( cur_per )
}

per_snvs_sum <- NULL
per_snvs_sum <- per_snvs_tib %>% 
  dplyr::group_by( Added_Probe,Sample,Platform ) %>%
  dplyr::summarise( Tot = n(),
                    Avg = mean(Match_Per, na.rm=TRUE ),
                    Med = median( Match_Per, na.rm=TRUE ),
                    Sds = sd( Match_Per, na.rm=TRUE ), 
                    Mad = mad( Match_Per, na.rm=TRUE ), 
                    .groups = "drop" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Script for Screen Probe Analytical Performance
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

# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("microbenchmark", quietly = TRUE) ) )

# Load Plotting Packages::
# suppressWarnings(suppressPackageStartupMessages(
#   base::require("ggplot2", quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages(
#   base::require("GGally", quietly = TRUE) ) )

# BiocManager::install("sesame")

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
par$prgm_tag <- 'stable_imSesameCpp_MSAv1'
par$verbose  <- 3
local_paths  <- c( 
  "/Users/bbarnes/Documents/tools/imSuite/scripts/R"
)

# local_paths  <- c( 
#   "/Users/bbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
#   "/Users/bretbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
#   "/illumina/scratch/darkmatter/tools/Workhorse-Unstained/scripts/R",
#   "/repo/Workhorse-Unstained/scripts/R" )

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
par <- source_functions( pars = par, rcpp = FALSE, vb = par$verbose )
par <- params_check( pars = par, args = args, 
                     prgm_aux_check = FALSE, vb = par$verbose )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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

# par$run_name <- "EPICv1"
# par$run_name <- "Embarkv1"
# par$run_name <- "FAILv1"
# par$run_name <- "COREv1"

par$run_name <- "MSAv03"
par$run_name <- "MSAv10"

opt <- NULL
opt <- imProbeQC_options( pars = par, args = args, vb = par$verbose )
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
opt_reqs <- c( 'out_path', 
               'Rscript', 'verbose' )

opt$rcpp <- 0
opt$rcpp <- 2
opt$rcpp <- 3
prgm_dat <- program_init( name = par$prgm_tag,
                          opts = opt, opt_reqs = opt_reqs,
                          pars = par, par_reqs = par_reqs, 
                          rcpp = opt$rcpp,
                          vb = opt$verbose, vt=3, tc=0 )

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

success = TRUE;


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Pre-processing:: Search for Idats (MSA v.1.0)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa10_idat_list <- NULL
msa10_idat_path <- file.path( opt$top_path, "data/pre-idats/MSA" )
msa10_idat_list <- sesame::searchIDATprefixes( dir.name = msa10_idat_path )

epic2_idat_list <- NULL
epic2_idat_path <- file.path( opt$top_path, "data/idats/idats_EPIC_v2-20220912-Alpha_subset" )
epic2_idat_list <- sesame::searchIDATprefixes( dir.name = epic2_idat_path )

# A tibble: 655 × 2
idat_path_tib <- NULL
idat_path_tib <- c( msa10_idat_list,epic2_idat_list ) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column( var = "Sentrix_Name") %>% 
  magrittr::set_names( c("Sentrix_Name", "Sentrix_Path") ) %>%
  tibble::as_tibble() %>%
  dplyr::distinct( Sentrix_Name, .keep_all = TRUE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Pre-processing:: Sample Sheet (MSA.v.1.0)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa_ssh_tib <- NULL
msa_ssh_tib <- format_msa_sshs( 
  paths = c( file.path( opt$top_path, "data/pre-idats/MSA/sample_sheets/rep0"),
             file.path( opt$top_path, "data/pre-idats/MSA/sample_sheets/uhm0"),
             file.path( opt$top_path, "data/pre-idats/MSA/sample_sheets/uhm1") ),
  
  out_dir    = file.path( opt$out_path, "sample_sheets" ),
  run_tag    = opt$run_name,
  reload     = opt$reload,
  reload_min = 10,
  ret_data   = FALSE,
  parallel   = opt$parallel,
  write_out  = FALSE,
  
  vb=vb,vt=vt+1,tc=tc+1, tt=tt )

msa_ssh_sum <- NULL
msa_ssh_sum <- msa_ssh_tib %>%
  print_sum( vec = c("Sample_Group","Source_Key","Source_Name","Sample_Base",
                     "Sample_Input","Sample_Titration","User_Format"),
             vb=vb,vt=vt+1,tc=tc+1, tt=tt )

# [TBD]: Update Source_Key
# [TBD]: Add Platform (or can we just preeict it), look at Rccp output...

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Pre-processing:: Sample Sheets (EPICv1/v2)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epi_ssh_tib <- NULL
epi_ssh_tib <- format_epic_ssh_alpha( 
  file = file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/formatted/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz" ),
  write_out  = FALSE,
  vb=vb,vt=vt+3,tc=tc+1, tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Pre-processing:: Sample Sheets (Combine All)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# [TBD]: Where did these go???
#
# A tibble: 25 × 8
ssh_mis_tib <- NULL
ssh_mis_tib <- dplyr::bind_rows( epi_ssh_tib,msa_ssh_tib ) %>% 
  dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
  dplyr::anti_join( idat_path_tib, 
                    by=c("Sentrix_Name") ) # %>% as.data.frame()

# Join and Match::
ssh_tib <- NULL
ssh_tib <- dplyr::bind_rows( epi_ssh_tib,msa_ssh_tib ) %>% 
  dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
  dplyr::inner_join( idat_path_tib, by=c("Sentrix_Name") ) %>%
  dplyr::mutate(
    Sample_Base = dplyr::case_when(
      Sample_Group == "MeTitration" ~ paste0( Sample_Base,Sample_Titration ),
      TRUE ~ Sample_Base ),
    Sample_Group = dplyr::case_when(
      Sample_Base == "NA12873" ~ "CellLine",
      TRUE ~ Sample_Group )
  ) %>%
  dplyr::arrange( Sample_Titration,Sample_Input,Sample_Base )

# [DONE]: Fix this 
#   22 NA12873      Coriell          10
ssh_tib %>% dplyr::group_by( Sample_Base,Sample_Group ) %>% dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)

ssh_tib %>% dplyr::group_by( Sample_Base,Sample_Group ) %>% dplyr::summarise( Count=n(), .groups = "drop" ) %>% 
  dplyr::mutate( Char2 = paste0( stringr::str_sub(Sample_Base,1,2),stringr::str_sub(Sample_Group,1,2) ) ) %>% print(n=1000)

ssh_sum0 <- NULL
ssh_sum0 <- ssh_tib %>% 
  print_sum( vec = c("Sample_Group","Source_Key","Sample_Base"),
             vb=vb,vt=vt+3,tc=tc+1, tt=tt )

ssh_sum1 <- NULL
ssh_sum1 <- ssh_tib %>%
  print_sum( vec = c("Sample_Base","Source_Key",
                     "Sample_Input","Sample_Titration") ,
             vb=vb,vt=vt+3,tc=tc+1, tt=tt )

ssh_sum2 <- NULL
ssh_sum2 <- ssh_tib %>%
  print_sum( vec = c("Sample_Group","Sample_Base",
                     "Sample_Input","Sample_Titration") ,
             vb=vb,vt=vt+3,tc=tc+1, tt=tt )

ssh_sum3 <- NULL
ssh_sum3 <- ssh_tib %>%
  print_sum( vec = c("Sample_Group","Sample_Base",
                     "Sample_Input","Sample_Titration",
                     "Source_Key","Source_Name" ),
             vb=vb,vt=vt+3,tc=tc+1, tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#             [TBD] Pre-processing:: Build Sample Experiments
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssh_tib$Source_Key %>% unique()
ssh_tib$Sample_Input %>% unique()
ssh_tib$Sample_Titration %>% unique()

# Only T0,T1,R0,E2
# Only >= 50 ng
Sample_Input_min <- 50
Source_Key_vec <- c( "T0","T1","R0","E2" )
base_ssh_tib <- NULL
base_ssh_tib <- ssh_tib %>% 
  dplyr::filter( Source_Key %in% Source_Key_vec ) %>% 
  dplyr::filter( Sample_Input >= Sample_Input_min ) %>%
  dplyr::arrange( Sample_Titration,Sample_Input,Sample_Base )

#
# - Sample_Group = [  MeTitration ] Titration
uhm_base_ssh_tib <- NULL
uhm_base_ssh_tib <- base_ssh_tib %>% 
  dplyr::filter( Sample_Group == "MeTitration" ) %>%
  dplyr::filter( Source_Key != "E2" ) %>%
  dplyr::filter( Source_Key != "T0" ) %>%
  dplyr::mutate( Sample_Base = paste0( Sample_Base,Sample_Titration ) )

#
# - Sample_Group = [ !MeTitration ] Replicate
rep_base_ssh_tib <- NULL
rep_base_ssh_tib <- base_ssh_tib %>% dplyr::filter( Sample_Group != "MeTitration" )

#
# - Sample_Group = [ !MeTitration ] Replicate:: CellLine
cel_base_ssh_tib <- NULL
cel_base_ssh_tib <- base_ssh_tib %>% 
  dplyr::filter( Sample_Group != "MeTitration" ) %>%
  dplyr::filter( Sample_Group == "CellLine" )

# rep_base_ssh_tib %>% dplyr::group_by( Sample_Group,Sample_Base,Sample_Input ) %>% summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)

all_ssh_listUHM1 <- NULL
all_ssh_listUHM1 <- dplyr::bind_rows(
  uhm_base_ssh_tib %>% 
    dplyr::filter( Sample_Titration ==  0 | Sample_Titration ==  50 ) %>% 
    dplyr::mutate( Exp_Str = "T1_UH" ) %>%
    dplyr::filter( Source_Key == "T1" & Sample_Input >= 50 ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  uhm_base_ssh_tib %>% 
    dplyr::filter( Sample_Titration ==  0 | Sample_Titration == 100 ) %>% 
    dplyr::mutate( Exp_Str = "T1_UM" ) %>%
    dplyr::filter( Source_Key == "T1" & Sample_Input >= 50 ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  uhm_base_ssh_tib %>% 
    dplyr::filter( Sample_Titration == 50 | Sample_Titration == 100 ) %>% 
    dplyr::mutate( Exp_Str = "T1_HM" ) %>%
    dplyr::filter( Source_Key == "T1" & Sample_Input >= 50 ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() )
  
) %>% split( .$Exp_Str )

all_ssh_listREP0 <- NULL
all_ssh_listREP0 <- dplyr::bind_rows(
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "HELA" ) %>% 
    dplyr::filter( Source_Key == "R0") %>% 
    dplyr::filter( Sample_Input == 50 ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"50", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "HELA" ) %>% 
    dplyr::filter( Source_Key == "R0") %>% 
    dplyr::filter( Sample_Input == 100 ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"100", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "HELA" ) %>% 
    dplyr::filter( Source_Key == "R0") %>% 
    dplyr::filter( Sample_Input == 200 ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"200", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "JURKAT" ) %>% 
    dplyr::filter( Source_Key == "R0") %>% 
    dplyr::filter( Sample_Input == 50 ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"50", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "JURKAT" ) %>% 
    dplyr::filter( Source_Key == "R0") %>% 
    dplyr::filter( Sample_Input == 100 ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"100", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "JURKAT" ) %>% 
    dplyr::filter( Source_Key == "R0") %>% 
    dplyr::filter( Sample_Input == 200 ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"200", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() )
  
) %>% split( .$Exp_Str )

all_ssh_listREP1 <- NULL
all_ssh_listREP1 <- dplyr::bind_rows(
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "HELA" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50 & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste0(Source_Key,"_",stringr::str_sub(Sample_Base,1,4),Sample_Input) ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "HELA" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste0(Source_Key,"_",stringr::str_sub(Sample_Base,1,4),Sample_Input) ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),

  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "JURKAT" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste0(Source_Key,"_",stringr::str_sub(Sample_Base,1,4),Sample_Input) ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "JURKAT" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste0(Source_Key,"_",stringr::str_sub(Sample_Base,1,4),Sample_Input) ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "K562" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste0(Source_Key,"_",stringr::str_sub(Sample_Base,1,4),Sample_Input) ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "K562" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste0(Source_Key,"_",stringr::str_sub(Sample_Base,1,4),Sample_Input) ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "MCF7" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste0(Source_Key,"_",stringr::str_sub(Sample_Base,1,4),Sample_Input) ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  # LEFT OFF REPLACING HERE::
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "MCF7" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"250", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA11922" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"50", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA11922" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"250", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA12752" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"50", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA12752" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"250", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA12877" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"50", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA12877" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"250", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA12878" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"50", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA12878" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"250", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA12882" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"50", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA12882" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"250", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA1879" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 50  & User_Format == "48x1") %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"50", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() ),
  cel_base_ssh_tib %>% 
    dplyr::filter( Sample_Base == "NA1879" ) %>% 
    dplyr::filter( Source_Key == "T1" & Sample_Input == 250 & User_Format == "48x1" ) %>%
    dplyr::mutate( Exp_Str = paste(Source_Key,Sample_Group,Sample_Base,"250", sep="_") ) %>%
    dplyr::select( Exp_Str, Sentrix_Name,Sample_Base,Sample_Input,Sample_Titration,Sample_Group,User_Format, dplyr::everything() )
  
) %>% split( .$Exp_Str )

all_ssh_list <- NULL
all_ssh_list <- c( all_ssh_listUHM1,all_ssh_listREP0,all_ssh_listREP1 )

# all_ssh_list %>% lapply( function(x) { print( x, n=1000 ) })
# all_ssh_list %>% lapply( function(x) { base::nrow(x) }) %>% dplyr::bind_rows( .id = "Experiment_Str" )

# rep_base_ssh_tib %>% dplyr::mutate( Exp_Str = paste("replication",Sample_Base,Sample_Input, sep="_") ) %>% split(.$Exp_Str )
# rep_ssh_list <- dplyr::bind_rows(
#   rep_base_ssh_tib %>% dplyr::mutate( Exp_Str = paste("replication",Sample_Base, sep="_") ) %>% split(.$Exp_Str )
# )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Pre-processing:: Previous Products
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_man_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_manifests.sub8.rds")
rev_man_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/rev_manifests.sub8.rds")
neg_ctl_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_negative_ctls.rds" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Negative Controls
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

neg_ctls_tib  <- NULL
neg_ctls_tib  <- readr::read_rds( neg_ctl_rds ) %>%
  dplyr::mutate( Probe_ID = paste0(Probe_ID,"_0") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Manifest (ALL)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_man_tib  <- NULL
all_man_tib  <- readr::read_rds( all_man_rds )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#              Pre-processing:: Manifest (MSA10) OLD METHOD::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# A tibble: 301,805 × 10
# msa_man_csv <- file.path( opt$out_path, "manifests", "msa_man.csv.gz")
msa_man_tib <- NULL
msa_man_csv <- "/Users/bbarnes/Documents/scratch/stable_imSesameCpp_MSAv1/MSAv10-NCBI-v4/manifests/msa_man.csv.gz"
msa_man_tib <- readr::read_csv( file = msa_man_csv, show_col_types = FALSE )

if ( FALSE ) {
  org_neg_tib  <- NULL
  org_neg_tib  <- readr::read_rds( neg_ctl_rds ) %>%
    dplyr::mutate( Probe_ID = paste0(Probe_ID,"_0") )
  
  ssh_all_tib <- NULL
  ssh_all_tib <- ssh_rep0_tib %>%
    # dplyr::bind_rows( ssh_uhm0_tib ) %>%
    dplyr::bind_rows( ssh_uhm1_tib )
  
  ssh_all_sum0 <- ssh_all_tib %>%
    print_sum( vec = c("Source_Key","Source_Name","Sample_Base",
                       "Sample_Input","Sample_Titration","User_Format") )
  
  ssh_all_sum1 <- ssh_all_tib %>%
    print_sum( vec = c("Sample_Base",
                       "Sample_Input","Sample_Titration","User_Format") )
  
  ssh_all_csv <- file.path( opt$out_path, "MAS.v.1.0_sample_sheet.csv.gz" )
  readr::write_csv( x = ssh_all_tib, file = ssh_all_csv )
  
  v03_ctl_tib <- NULL
  v03_ctl_tib <- v03_man_tib %>%
    dplyr::filter( !Probe_Type == "cg" ) %>%
    dplyr::filter( !Probe_Type == "ch" ) %>%
    dplyr::filter( !Probe_Type == "rs" )
  
  v03_ctl_sum <- NULL
  v03_ctl_sum <- print_sum( tib = v03_ctl_tib, vec = c("Probe_Type") )
  
  v03_man_tib <- NULL
  v03_man_csv <- file.path( opt$top_path, "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_SS_BP123_A1.csv.gz" )
  v03_man_tib <- readr::read_csv( file = v03_man_csv, show_col_types = FALSE ) %>%
    clean_tib() %>%
    dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("^ZZ", "ctl_"),
                   Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_neg", "ctl_Negative_"),
                   Name = Probe_ID %>% stringr::str_remove("_[^_]+$"),
                   M = dplyr::case_when(
                     is.na(M) ~ 0.0,
                     TRUE ~ M
                   ),
                   U = U %>% as.integer(),
                   M = M %>% as.integer(),
                   col = dplyr::case_when(
                     is.na(col) ~ "2",
                     TRUE ~ col
                   ),
                   Manifest = "MSA03",
                   Manifest_Version = "V3",
                   Annotation = BP,
                   Locus_Name = Probe_ID %>% stringr::str_remove("_.*$"),
                   # Probe_Type = Probe_ID %>% stringr::str_sub(1,2),
                   SNP_Probe_ID = NA ) %>%
    dplyr::rename( Chromosome = CHR,
                   Coordinate = MAPINFO ) %>%
    dplyr::select( Probe_ID, U, M, col, Name,
                   Species, Manifest, Manifest_Version, Annotation,
                   Chromosome, Coordinate, Probe_Type, Locus_Name,
                   SNP_Probe_ID, AlleleA_ProbeSeq )
  
  
  # A tibble: 4,077 × 3
  v03_neg_tib <- NULL
  v03_neg_tib <- v03_ctl_tib %>%
    dplyr::filter( Probe_Type == "NEGATIVE" ) %>%
    dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") ) %>%
    dplyr::mutate(
      # Probe_ID = Probe_ID %>% stringr::str_remove("^ctl_"),
      Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_Negative_", "neg"),
      Address = U,
      Sequence = AlleleA_ProbeSeq
    ) %>%
    dplyr::select( Address, Probe_ID, Sequence ) %>%
    dplyr::arrange( Probe_ID ) %>%
    dplyr::distinct( Address, .keep_all = TRUE )
  
  all_neg_tib <- NULL
  all_neg_tib <- dplyr::bind_rows( org_neg_tib,v03_neg_tib ) %>%
    dplyr::arrange( Sequence) %>%
    dplyr::distinct( Address, .keep_all = TRUE )
  
  v03_ctl_tib <- NULL
  v03_ctl_tib <- v03_man_tib %>%
    dplyr::filter( !Probe_Type == "cg" ) %>%
    dplyr::filter( !Probe_Type == "ch" ) %>%
    dplyr::filter( !Probe_Type == "rs" )
  
  v03_ctl_sum <- NULL
  v03_ctl_sum <- print_sum( tib = v03_ctl_tib, vec = c("Probe_Type") )
  
  v10_man_tib <- NULL
  v10_man_csv <- file.path( opt$top_path, "data/pre-idats/MSA/MSA-48v1-0-Post-PQC_2B.body.csv.gz" )
  v10_man_tib <- readr::read_csv( file = v10_man_csv, show_col_types = FALSE ) %>% clean_tib()
  
  msa_man_csv <- file.path( opt$out_path, "manifests", "msa_man.csv.gz")
  msa_man_tib <- NULL
  msa_man_tib <- dplyr::bind_rows(
    v10_man_tib %>%
      dplyr::mutate(
        Probe_ID = IlmnID,
        U = dplyr::case_when(
          is.na(AddressA_ID) ~ 0,
          TRUE ~ AddressA_ID ),
        M = dplyr::case_when(
          is.na(AddressB_ID) ~ 0,
          TRUE ~ AddressB_ID ),
        col = dplyr::case_when(
          is.na(col) & M==0 ~ '2',
          U != 0 & M != 0 ~ col,
          TRUE ~ NA_character_
        ),
        Manifest = "MSA",
        Annotation = dplyr::case_when(
          CHR == "chrY" ~ "ChrY",
          TRUE ~ ""
        ),
        Chromosome = CHR
      ) %>%
      dplyr::select( Probe_ID, U,M, col, Manifest, Annotation, Chromosome ),
    v03_ctl_tib %>%
      dplyr::mutate(
        Annotation = Probe_Type,
        Chromosome = "chr0"
      ) %>%
      dplyr::select( Probe_ID, U,M, col,
                     Species,
                     Manifest, Manifest_Version,
                     Annotation, Chromosome )
  ) %>%
    dplyr::mutate(
      Name = Probe_ID %>% stringr::str_remove("_.*$"),
      Species = "homo_sapiens",
      Manifest = "MSA",
      Manifest_Version = "V1"
    ) %>%
    dplyr::select( Probe_ID, U,M, col, Name,
                   Species, Manifest, Manifest_Version,
                   Annotation, Chromosome )
  readr::write_csv( x = msa_man_tib, file = msa_man_csv )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Manifest (MSA10)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa_all_ctls_tib <- NULL
msa_all_ctls_tib <- format_msa_ctls( 
  neg_man_csv = file.path( opt$top_path, "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_SS_BP123_A1.csv.gz" ),
  
  out_dir    = file.path( opt$out_path, "manifests/inputs" ),
  run_tag    = opt$run_name,
  reload     = opt$reload,
  reload_min = 10,
  ret_data   = FALSE,
  parallel   = opt$parallel,
  write_out  = FALSE,
  
  vb=vb,vt=vt+1,tc=tc+1, tt=tt )

#
# All Controls (EPIC + MSA):: Rccp Format
#
msa_neg_ctls_tib <- NULL
msa_neg_ctls_tib <- msa_ctls_to_negs( msa_all_ctls_tib ) %>%
  dplyr::bind_rows( neg_ctls_tib ) %>%
  dplyr::arrange( Address ) %>%
  dplyr::distinct( Address, .keep_all = TRUE )

msa_sel_ctls_tib <- NULL
msa_sel_ctls_csv <- file.path( opt$top_path, "data/pre-idats/MSA/MSA-48v1-0-Post-PQC_2B.ctls-wHeader.csv.gz" )
msa_sel_ctls_tib <- readr::read_csv( file = msa_sel_ctls_csv, show_col_types = FALSE ) %>%
  clean_tib()

#
# Controls Only (Sesame/Rccp Format)::
#
msa_ses_ctls_tib <- NULL
msa_ses_ctls_tib <- msa_all_ctls_tib %>% 
  dplyr::filter( U %in% msa_sel_ctls_tib$Address ) %>% 
  dplyr::select( Probe_ID,U,M,col,Name, Species,Manifest,Manifest_Version,
                 Annotation,Chromosome,Probe_Type ) %>%
  dplyr::mutate( Probe_Type = "ct" )

#
# Analytical + Controls (Sesame/Rccp Format)::
#
msa_man_full_tib <- NULL
msa_man_full_tib <- format_msa_mans( 
  v10_man_csv = file.path( opt$top_path, "data/pre-idats/MSA/MSA-48v1-0-Post-PQC_2B.body.csv.gz" ), 
  
  out_dir    = file.path( opt$out_path, "manifests/inputs" ),
  run_tag    = opt$run_name,
  reload     = opt$reload,
  reload_min = 10,
  ret_data   = FALSE,
  parallel   = opt$parallel,
  write_out  = FALSE,
  
  vb=vb,vt=vt+1,tc=tc+1, tt=tt ) %>%
  # dplyr::bind_rows( msa_ses_ctls_tib )
  dplyr::bind_rows( 
    # msa_ses_ctls_tib %>% dplyr::select( Probe_ID,U,M,col,Name,Species,Manifest,Manifest_Version,Annotation,Chromosome,Probe_Type ) %>%
    msa_ses_ctls_tib %>% dplyr::select( Probe_ID,U,M,col,Name,Species,Manifest,Manifest_Version,Annotation,Chromosome ) %>%
      dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("ctl_Negative_","neg_") ) )
  # dplyr::bind_rows( 
  #   msa_ses_ctls_tib %>% dplyr::mutate( 
  #     Probe_ID = Probe_ID %>% stringr::str_replace("ctl_Negative_","neg_") ) )

msa_man_body_sum <- NULL
msa_man_body_sum <- msa_man_full_tib %>% 
  # print_sum( vec = c("Manifest","Manifest_Version","Probe_Type"),
  print_sum( vec = c("Manifest","Manifest_Version"),
             vb=vb,vt=vt+1,tc=tc+1, tt=tt )

#
# Interesting Summaries to be Reviewed... Probably Not working now...
#
if ( FALSE ) {
  
  msa_man_tib %>% 
    dplyr::mutate( Probe_Type = stringr::str_sub( Probe_ID, 1,2) ) %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  all_ctl_tib <- NULL
  all_ctl_tib <- all_man_tib %>% 
    dplyr::filter( Probe_Type != "cg" ) %>% 
    dplyr::filter( Probe_Type != "ch" ) %>% 
    dplyr::filter( Probe_Type != "rs" ) %>% 
    dplyr::filter( Probe_Type != "nv" )
  
  all_ctl_sum <- NULL
  all_ctl_sum <- all_ctl_tib %>%
    dplyr::group_by( Manifest,Manifest_Version,Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
}

#
# Investigate difference: msa_man_tib vs. msa_man_full_tib
#

if ( FALSE ) {
  msa_man_cnt <- msa_man_tib %>% base::nrow()
  msa_man_unq <- msa_man_tib %>% dplyr::distinct( U ) %>% base::nrow()
  msa_man_full_unq <- msa_man_full_tib %>% base::nrow()
  msa_man_full_cnt <- msa_man_full_tib %>% dplyr::distinct( U ) %>% base::nrow()
  
  if ( p0 ) cat(glue::glue("{pmssg}      msa_man_cnt='{msa_man_cnt}'.{RET2}"))
  if ( p0 ) cat(glue::glue("{pmssg}      msa_man_unq='{msa_man_unq}'.{RET2}"))
  if ( p0 ) cat(glue::glue("{pmssg} msa_man_full_unq='{msa_man_full_unq}'.{RET2}"))
  if ( p0 ) cat(glue::glue("{pmssg} msa_man_full_cnt='{msa_man_full_cnt}'.{RET2}"))
  
  msa_man_tib %>% dplyr::anti_join( msa_man_full_tib, by=c("U") )
  msa_man_full_tib %>% dplyr::anti_join( msa_man_tib, by=c("U") )
  
  msa_man_tib %>% dplyr::group_by( Annotation ) %>% dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  msa_man_full_tib %>% dplyr::group_by( Annotation ) %>% dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  
  # Qs: How is Annotation used?
  # Qs: Why are there missing controls
  # Qs: Can we write/read bgz or just write sdf.rds..
  #
  # Conclusion: Use msa_man_tib for now and fix later...
}

# v1_body_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv1/rank_historic_probes/EPICv2-UCSC-v5.rank_historic_probes.csv.gz" )
# v2_body_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv2/rank_historic_probes/EPICv2-UCSC-v5.rank_historic_probes.csv.gz" )
# v1_body_tib <- NULL
# v1_body_tib <- readr::read_csv( file = v1_body_csv, show_col_types = FALSE )
# v2_body_tib <- NULL
# v2_body_tib <- readr::read_csv( file = v2_body_csv, show_col_types = FALSE )

v2_dock_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/manifest/EPICv2/EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
v2_dock_tib <- NULL
v2_dock_tib <- readr::read_csv( file = v2_dock_csv, show_col_types = FALSE ) %>%
  clean_tib()
v2_dock_sum <- NULL
v2_dock_sum <- v2_dock_tib %>% print_sum( vec = c("Probe_Type"),
                                          vb=vb,vt=vt+1,tc=tc+1, tt=tt )

#
# [TBD]: Load actual EPICv2 manifest and match by sequence...
#
# MSA EPIC Overlap::
if ( FALSE ) {
  
  msa_man_full_tib %>% dplyr::filter( Name %in% v2_dock_tib$Probe_ID )
  msa_man_full_tib %>% dplyr::filter( Probe_ID %in% v2_dock_tib$Full_ID )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Run-Time Variables
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tc <- 0
vt <- 0

# vb <- 5 # Two level deep (p3)
# vb <- 4 # One level deep (p2)
# vb <- 3 # Standard
# vb <- 2 # Light
# vb <- 1 # Min
vb <- 0 # None

success <- FALSE

opt$min_pval <- 0.05
opt$min_beta <- 0.3
opt$max_beta <- 0.7
opt$min_perO <- 0.75
opt$min_perI <- 0.05

# Mask Methods:: csJjKk
# workflow_vec = c( "dc",
#                    "ds",
#                    "dJ",
#                    "dj",
#                    "dK",
#                    "dk" )
workflow_vec = c( "cd",
                  "sd",
                  "Jd",
                  "jd",
                  "Kd",
                  "kd" )
workflow_vec = c( "id" )

workflow_vec = c( "i" )

opt$return_df <- 0
opt$return_df <- 1
opt$return_df <- 2
opt$return_df <- 3
# opt$return_df <- 4

opt$max_sam <- 4
opt$max_rep <- 8

#
# TBD:: Scatter plots with random sub-selection::
#
opt$max_sam <- 10000
opt$max_rep <- 6

# opt$max_sam <- 3
# opt$max_rep <- 3
# opt$single <- TRUE
opt$single <- TRUE
opt$single <- FALSE
opt$min_dB <- 0.2

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Data-processing:: Standard Sesame::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

exp_idat_lst <- NULL
# exp_idat_lst <- c( msa10_idat_list,epic2_idat_list )
exp_idat_lst <- msa10_idat_list
exp_man_tib  <- msa_man_tib
# exp_ctl_tib  <- msa_man_tib %>% dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_") )
exp_ctl_tib  <- msa_ses_ctls_tib
exp_out_dir  <- safe_mkdir( file.path(opt$out_path, "sdf/core") )

# opt$parallel <- FALSE
if ( FALSE ) {
  cur_sdfs <- NULL
  if ( opt$parallel ) {
    cur_sdfs <- mclapply( exp_idat_lst,
                          prefix_to_sdf, 
                          platform = "EPIC", 
                          manifest = exp_man_tib, 
                          controls = exp_ctl_tib, 
                          out_dir  = exp_out_dir,
                          run_tag  = opt$run_name,
                          reload   = opt$reload, 
                          reload_min = 10,
                          parallel = opt$parallel, 
                          vb=vb, vt=vt, tc=0, tt = tt )
  } else {
    cur_sdfs <- lapply( exp_idat_lst,
                        prefix_to_sdf, 
                        platform = "EPIC", 
                        manifest = exp_man_tib, 
                        controls = exp_ctl_tib, 
                        out_dir  = exp_out_dir,
                        run_tag  = opt$run_name,
                        reload   = opt$reload, 
                        reload_min = 10,
                        parallel = opt$parallel, 
                        vb=vb, vt=vt, tc=0, tt = tt )
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Data-processing:: Titraiton Sample List
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

exp_auto_ssh_tib <- NULL

uhm_dB_tib <- NULL
uhm_data_tib <- NULL
negs_sdf_tab <- NULL
rm_outliers_vec <- c(TRUE, FALSE)
for ( rm_outliers in rm_outliers_vec ) {
  for ( exp_name in names(all_ssh_list) ) {
    all_beta_tibs <- NULL
    all_poob_tibs <- NULL
    all_negs_tibs <- NULL

    # all_beta_tibs[[2]] <- NULL
    # all_poob_tibs[[2]] <- NULL
    # all_negs_tibs[[2]] <- NULL
    
    # exp_name <- "titration_0_100"
    # exp_name <- "uhm_0_50"
    
    rm_outliers_str <- "rm-outliers"
    if ( !rm_outliers ) rm_outliers_str <- "ok-outliers"
    
    cur_out_path <- NULL
    cur_out_path <- safe_mkdir( dir = file.path(opt$out_path, "analysis",rm_outliers_str, paste0(exp_name,"_",rm_outliers_str) ) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Data-processing:: Sentrix Sample List
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # [TBD]: Group by Sample_Base
    ssh_sel_tibs <- NULL
    ssh_sel_tibs <-     all_ssh_list[[exp_name]] %>% 
      dplyr::group_by( Sample_Base ) %>%
      dplyr::mutate( Rep_Idx = dplyr::row_number(),
                     Rep_Key = paste( Sample_Base,Rep_Idx, sep="_") ) %>%
      dplyr::ungroup()

    ssh_sel_list <- NULL
    ssh_sel_list <- ssh_sel_tibs %>% 
      split( .$Sentrix_Name ) # %>% head( n=opt$max_rep )
    
    # sentrix_name <- names(ssh_sel_list)[1]
    # idat_path_tib
    
    # [TBD]: Screen this earlier...
    if ( ssh_sel_list %>% length() < 3 ) next

    for ( sentrix_name in names(ssh_sel_list) ) {
      
      prefix      <- ssh_sel_list[[sentrix_name]]$Sentrix_Path
      src_key     <- ssh_sel_list[[sentrix_name]]$Source_Key
      ng_input    <- ssh_sel_list[[sentrix_name]]$Sample_Input
      sample_base <- ssh_sel_list[[sentrix_name]]$Sample_Base
      if ( p1 ) cat(glue::glue("{pmssg} Current prefix({exp_name}) = '{prefix}'{RET}"))
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                       Data-processing:: Loading
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      # prefix <- "/Users/bbarnes/Documents/data/idats/idats_MSAv03_48x1_Alpha/207545400008/207545400008_R04C01"
      # prefix <- "/Users/bbarnes/Documents/data/idats/EPIC/203010460029/203010460029_R12C01"
      # prefix <- "/Users/bbarnes/Documents/data/idats/EPIC/201502830033/201502830033_R02C01"
      # all_man_tib  <- NULL
      # all_man_tib  <- readr::read_rds( all_sub_rds ) # %>% dplyr::bind_rows( msa_man_tib )
      # all_man_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% dplyr::group_by( Manifest, Manifest_Version ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
      
      cur_out_name <- NULL
      cur_out_name <- paste0(sentrix_name,".",exp_name,"_",rm_outliers_str)
      cur_negs_csv <- NULL
      cur_negs_csv <- file.path( cur_out_path, paste0(cur_out_name,".negs.sdf.csv.gz") )
      
      # msa_neg_ctls_tib %>% dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("ctl_Negative_","neg_") ) %>% tail()
      # msa_neg_ctls_tib %>% dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("ctl_Negative_","neg_") ) %>% dplyr::pull( Probe_ID ) %>% as.vector()
      
      # [TBD]: Update Output Format: Include [ Pval,Beta ]
      # [TBD]: Should it the min or max pval?
      # [TBD]: Is swapping happening?
      # [TBD]: Second round of p-value filtering after outlier removal
      # [TBD]: Compare against zero filtering
      # [TBD]: Compare against direct Sesame loading...
      #
      # [TBD]: Calculate Genotype in Rccp
      #   - Use Actual Bead Count in the model. 
      #   - Basically Comare sum(R)/sum(GR) vs. sum(G)/sum(GR) compare against col (next base)
      #   - Test on known EPICv2 calls...
      #
      
      # Seame Code: /Users/bbarnes/Documents/tools/programs/sesame
      #
      # [TBD]: set_genotype()
      # [TBD]: Swap and Genotyping (Fingerprinting) needs to happen at the same time in the same function...
      #        GT = Genotype Call 0=0/0, 1=0/1, 2=1/1, 3=NA | col
      #        GS = Genotype Score
      #             [TBD]: Need boost binomalial distribution...
      #             [TBD]: Weighted getAFTypeIbySumAlleles() by 1-pval...
      #
      # A tibble: 301,803 × 16
      # Probe_ID           UG    UR    MG    MR col   mask  mask_0 Beta_0 Pval_0 UG_Det_0 UR_Det_0 MG_Det_0 MR_Det_0 Beta_1 Probe_Type [GT] [GS]
      # <chr>           <int> <int> <int> <int> <chr> <lgl> <lgl>   <dbl>  <dbl>    <dbl>    <dbl>    <dbl>    <dbl>  <dbl> <chr> <int> <dbl>
      # 1 cg06185909_TC11  3214   590  2398   527 G     FALSE FALSE   0.427      0   0         0       0         0      0.427 cg    X    S
      # 2 cg18975462_BC11  2382   931   274   397 G     FALSE FALSE   0.103      0   0         0       0.252     0.273  0.103 cg    X    S
      # 3 cg20516119_TC11   231  1250   189   370 R     FALSE FALSE   0.228      0   0.654     0       0.902     0.423  0.228 cg    X    S
      # 4 cg10149399_BC11   297   773   294  3122 R     FALSE FALSE   0.802      0   0.0677    0       0.0850    0      0.802 cg    X    S
      # 
      # [Wrong Manifest] A tibble: 293,649 × 8 [ rm_outliers == TRUE ]
      idat_tib <- NULL
      idat_tib <- read_idat_pair_rcpp( 
        prefix_path      = prefix,
        # prefix_path      = "/Users/bbarnes/Documents/data/pre-idats/MSA/MSA.replicats/07102023_MSAv1_Cancercell_ReRun/207675480005/207675480005_R02C01",
        output_path      = cur_out_path, # file.path( opt$out_path, cur_out_name ),
        workflow_vec     = workflow_vec %>% as.vector(),
        
        pval_add_vec     = msa_neg_ctls_tib$Address %>% as.vector(), 
        # msa_man_tib
        addU_man_vec     = msa_man_tib$U %>% as.vector(),
        addM_man_vec     = msa_man_tib$M %>% as.vector(),
        
        cgns_man_vec     = msa_man_tib$Probe_ID %>% as.vector(), 
        cols_man_vec     = msa_man_tib$col %>% as.vector(),
        keys_man_vec     = msa_man_tib$Manifest %>% as.vector(),
        anns_man_vec     = msa_man_tib$Annotation %>% as.vector(), 
        chrs_man_vec     = msa_man_tib$Chromosome %>% as.vector(),
        
        # pval_add_vec     = msa_neg_ctls_tib %>% dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("ctl_Negative_","neg_") ) %>% dplyr::pull( Address )
        # addU_man_vec     = msa_man_full_tib$U %>% as.vector(),
        # addM_man_vec     = msa_man_full_tib$M %>% as.vector(),
        # 
        # cgns_man_vec     = msa_man_full_tib$Probe_ID %>% as.vector(), 
        # cols_man_vec     = msa_man_full_tib$col %>% as.vector(),
        # keys_man_vec     = msa_man_full_tib$Manifest %>% as.vector(),
        # anns_man_vec     = msa_man_full_tib$Annotation %>% as.vector(), 
        # chrs_man_vec     = msa_man_full_tib$Chromosome %>% as.vector(),
        
        min_pval         = opt$min_pval,
        min_beta         = opt$min_beta,
        max_beta         = opt$max_beta,
        min_perO         = opt$min_perO, 
        min_perI         = opt$min_perI, 
        read_bgz         = FALSE,
        write_bgz        = FALSE, 
        
        rm_pval_outliers = rm_outliers,
        # rm_pval_outliers = FALSE,
        
        # return_df        = opt$return_df,
        return_df        = 10,
        
        vb=vb,vt=vt,tc=tc ) %>%
        # dplyr::mutate(across(where(is.factor), as.character) ) %>%
        dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
      
      # [TBD]: Build everything upfront and then reload as needed
      #        - This will require write_bgz()
      # [TBD]: Run without masking...
      # [TBD]: Run raw Sesame side by side

      # Quick Summaries::
      # [TBD]: Functionalize this stuff...
      if ( FALSE ) {
        
        # Negative Controls that fail:
        idat_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") ) %>% dplyr::filter( mask )
        
        sn <- 0.0001
        neg_beta_den_ggg <- NULL
        neg_beta_den_ggg <- idat_tib %>% 
          dplyr::filter( Probe_Type == "ct" ) %>% 
          dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") ) %>%
          ggplot2::ggplot( aes( x=Beta_1, fill=mask) ) +
          ggplot2::geom_density( alpha=0.2 )
        
        plot_tab <- NULL
        plot_tab <- idat_tib %>% 
          # dplyr::filter( Probe_Type == "ct" ) %>%
          dplyr::filter( Probe_Type == "nn" |
                           Probe_Type == "cg" | 
                           # Probe_Type == "ch" | 
                           # Probe_Type == "rs" |
                           Probe_ID %>% stringr::str_starts("ctl_Negative_") |
                           # Probe_ID %>% stringr::str_starts("ctl_BS_Conversion_I_") |
                           # Probe_ID %>% stringr::str_starts("ctl_BS_Conversion_II_") |
                           Probe_Type == "nn" ) %>%
          dplyr::mutate( 
            # Pval_Plot = as.integer( Pval_0 * 10 ),
            Pval_Plot = dplyr::case_when(
              Pval_0 <= 0.01 ~ "pval <= 01",
              Pval_0 <= 0.05 ~ "pval <= 05",
              Pval_0 <= 0.10 ~ "pval <= 10",
              Pval_0 <= 0.20 ~ "pval <= 20",
              Pval_0 <= 0.30 ~ "pval <= 30",
              Pval_0 >  0.30 ~ "pval >  30",
              
              TRUE ~ NA_character_ ),
            Plot_Name = dplyr::case_when(
              Probe_ID %>% stringr::str_starts("ctl_BS_Conversion_I_") ~ "B1",
              Probe_ID %>% stringr::str_starts("ctl_BS_Conversion_II_") ~ "B2",
              Probe_Type == "ct" ~ Probe_ID %>% stringr::str_remove("ctl_") %>% stringr::str_sub(1,2),
              TRUE ~ Probe_Type ),
            int_sum = UG+UR+MG+MR
          ) %>%
          # dplyr::filter( int_sum < 10000 ) %>%
          # dplyr::filter( int_sum < 2000 ) %>%
          # dplyr::filter( int_sum < 1000 ) %>%
          tidyr::pivot_longer( cols = c(UG,UR,MG,MR), 
                               names_to = "Color", 
                               values_to = "Intensity", 
                               values_drop_na = TRUE ) %>%
          dplyr::filter( Intensity != 0 ) 
        
        prb_ints_den_ggg <- plot_tab %>%
          # ggplot2::ggplot( aes( x=log(Intensity+sn), fill=Color) ) +
          ggplot2::ggplot( aes( x=Intensity, fill=Color ) ) +
          ggplot2::geom_density( alpha=0.2 ) +
          ggplot2::facet_grid( rows = vars(Pval_Plot), # cols = vars(mask),
                               cols = vars(Plot_Name) )
        prb_ints_den_ggg
        
        plot_tab %>% dplyr::group_by( Plot_Name,Pval_Plot ) %>%
          dplyr::summarise( Count=n(), .groups = "drop" )
        
        
        # Count Stats::
        idat_sum1 <- idat_tib %>% print_sum( vec = c("Probe_Type"),
                                             vb=vb+1,vt=vt,tc=tc+1, tt=tt )
        idat_sum2 <- idat_tib %>% print_sum( vec = c("Probe_Type","mask"),
                                             vb=vb+1,vt=vt,tc=tc+1, tt=tt )
        idat_sum3 <- idat_tib %>% print_sum( vec = c("Probe_Type","mask_0"),
                                             vb=vb+1,vt=vt,tc=tc+1, tt=tt )
      }

      # How to measure outlier removal::
      # 
      negs_sdf <- NULL
      negs_sdf <- idat_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% 
        dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") )
      readr::write_csv( x = negs_sdf, file = cur_negs_csv )
      
      negs_sdf_tab <- negs_sdf_tab %>% 
        dplyr::bind_rows(
          dplyr::mutate( negs_sdf,
                         Rm_Outliers = rm_outliers,
                         Sample = exp_name ) %>%
            dplyr::select( Rm_Outliers,Sample, dplyr::everything() )
        )
      
      # negs_sdf %>%
      #   dplyr::group_by( mask ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
      
      #
      # TBD:: Return select columns with appropriate data types in Rcpp...
      #   - This will allow the code above to be removed!
      #
      
      if ( is.null(idat_tib) ) {
        stop( glue::glue("{perrs} Failed to load idats: prefix='{prefix}'{RET}") )
        break
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                       Data-processing:: Summary
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      idat_sum <- NULL
      idat_sum <- print_sum( tib = idat_tib, vec = c("mask","Probe_Type","col"),
                             vb=vb,vt=vt,tc=tc, tt=tt )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Data-processing:: Parse Control Data
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # data_tib <- NULL
      # data_tib <- idat_tib
      # data_tib <- idat_tib %>% 
      #   dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") )
      
      ctls_tib <- NULL
      ctls_tib <- idat_tib %>% 
        dplyr::filter(  stringr::str_starts( Probe_ID, pattern = "ct") )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Data-processing:: Workflow Pipeline
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      work_sdf <- NULL
      # work_sdf <- data_tib %>%
      work_sdf <- idat_tib %>%
        # dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") ) %>%
        as.data.frame() %>%
        sesame::SigDF( platform = "EPIC", 
                       ctl = msa_all_ctls_tib
                       # ctl = ctls_tib
                       # ctl = NULL
        ) %>%
        mutate_sdf_simple( 
          steps = "DB", 
          negs_min = 1.0, 
          poob_min = 1.0, 
          vb=vb+3,vt=vt,tc=tc )
      
      rep_key <- NULL
      rep_key <- ssh_sel_list[[sentrix_name]]$Rep_Key
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Data-processing:: Extract Beta Values
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      beta_tib <- NULL
      beta_tib <- mutate_sdf_simple( 
        sdf = work_sdf,
        # steps = "v", # No Masking...
        steps = "V", # With masking...
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc ) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column( var = "Probe_ID") %>% 
        tibble::as_tibble() %>% 
        magrittr::set_names( value = c("Probe_ID",rep_key) )
      # magrittr::set_names( value = c("Probe_ID","Beta") )

      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Data-processing:: Extract pooBAH P-Values
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      poob_tib <- NULL
      poob_tib <- mutate_sdf_simple( 
        sdf = work_sdf,
        steps = "O", 
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc ) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column( var = "Probe_ID") %>% 
        tibble::as_tibble() %>% 
        magrittr::set_names( value = c("Probe_ID",rep_key) )
      # magrittr::set_names( value = c("Probe_ID","Poob") )
      
      poob_mat <- NULL
      poob_mat <- as.matrix( poob_tib %>% tibble::column_to_rownames( var = "Probe_ID" ) )
      
      # [TBD]: Read Sample Sheet, but this should be done in an analysis function
      #   so everything is separated: build -> analyze...
      cur_auto_ssh_tib <- NULL
      cur_auto_ssh_csv <- file.path( cur_out_path, paste0(sentrix_name,".sample_sheet.csv") )
      cur_auto_ssh_tib <- readr::read_csv( file = cur_auto_ssh_csv, show_col_types = FALSE ) %>%
        clean_tib() %>%
        dplyr::mutate( Remove_Outliers = rm_outliers, 
                       Experiment = exp_name,
                       Pass_Poob_Cnt = which(  poob_mat <= 0.05 ) %>% length(),
                       Fail_Poob_Cnt = which( !poob_mat <= 0.05 ) %>% length(),
                       Pass_Poob_Per = round( 100 * Pass_Poob_Cnt/ (Fail_Poob_Cnt+Pass_Poob_Cnt), 3 )
        ) %>% dplyr::select( Experiment, dplyr::everything() )
      
      exp_auto_ssh_tib <- exp_auto_ssh_tib %>% 
        dplyr::bind_rows( cur_auto_ssh_tib )
      
      if ( cur_auto_ssh_tib %>% dplyr::pull( Pass_Poob_Per ) <= 70 ) next
      
      #
      # Aggergate Data Sets if they pass...
      # [TBD]: Write all of this out and then select base on Call Rates...
      #

      if ( is.null( all_beta_tibs[[sample_base]] ) ) {
        all_beta_tibs[[sample_base]] <- beta_tib
      } else {
        all_beta_tibs[[sample_base]] <- all_beta_tibs[[sample_base]] %>% 
          inner_join( beta_tib, by = c("Probe_ID") )
      }
      
      if ( is.null( all_poob_tibs[[sample_base]] ) ) {
        all_poob_tibs[[sample_base]] <- poob_tib
      } else {
        all_poob_tibs[[sample_base]] <- all_poob_tibs[[sample_base]] %>% 
          inner_join( poob_tib, by = c("Probe_ID") )
      }
    }
    if ( p0 ) cat(glue::glue("{pmssg} Finished Processing Idat Pairs sample='{exp_name}'.{RET2}"))
    
    #
    # LEFT OFF HERE::
    #  - [TBD]: Incorporate: {ng_input,rep_key}, could just pull from the sample sheet...
    #
    #  - [TBD]: Need to implement dual input matrices
    #     - This will allow Titration efforts
    #
    # [TBD]: Incorporate weight score
    # [TBD]: Add sample names
    # [TBD]: Reduce output
    # [TBD]: Add EPICv2
    # [TBD]: Create Comparison Groups Sheet
    #
    # [Done]: Add More Stats...
    #
    # [TBD]: Split all_beta_tibs/all_poob_tibs into titration if titration experiment...
    #
    
    # [TBD]: Add Poob Call Rate to sample sheets above...
    # pval_pas_cnt <- which( as.matrix( all_poob_tibs[[1]] %>% tibble::column_to_rownames( var = "Probe_ID" ) ) <= 0.05 ) %>% length()
    # pval_mis_cnt <- which( as.matrix( all_poob_tibs[[1]] %>% tibble::column_to_rownames( var = "Probe_ID" ) ) >  0.05 ) %>% length()
    # pval_pas_per <- base::round( 100 * pval_pas_cnt / ( pval_pas_cnt+pval_mis_cnt ), 3 )
    
    cur_min_dB  <- opt$min_dB
    cur_cmp_str <- "lte"
    cur_is_abs  <- TRUE
    
    all_beta_tibs2 <- NULL
    all_poob_tibs2 <- NULL
    all_negs_tibs2 <- NULL
    exp_name2 <- NULL
    
    # [TBD]: Find a better way of spliting this up...
    # Pretty Sure this isn't needed anymore...
    # [TBD]: Remove code below...
    if ( exp_name %>% stringr::str_starts("uhm_") ) {
      ssh_uhm_list <- NULL
      ssh_uhm_list <- ssh_sel_tibs %>% split( .$Sample_Titration )
      exp_name2 <- ssh_uhm_list[[2]] %>% dplyr::distinct( Sample_Base ) %>% dplyr::pull( Sample_Base )
      
      all_beta_tibs2 <- all_beta_tibs[[2]]
      all_poob_tibs2 <- all_poob_tibs[[2]]
      # all_negs_tibs2 <- all_negs_tibs[[2]]
      
      if ( exp_name == "uhm_0_100" ) {
        cur_min_dB  <- 0.5
        cur_cmp_str <- "gt"
        cur_is_abs  <- FALSE
      }
      if ( exp_name == "uhm_0_50" ) {
        cur_min_dB  <- 0.2
        cur_cmp_str <- "gt"
        cur_is_abs  <- FALSE
      }
      if ( exp_name == "uhm_50_100" ) {
        cur_min_dB  <- 0.1
        cur_cmp_str <- "gt"
        cur_is_abs  <- FALSE
      }
    }
    
    # [TBD]: Update the usage of the name_strA/B. Currently doesn't really do antyhing...
    cur_dB_tib <- NULL
    cur_dB_tib <- calc_dBs( beta_tibA = all_beta_tibs[[1]],
                            pval_tibA = all_poob_tibs[[1]],
                            name_strA = exp_name, # names(all_beta_tibs[[1]])[1]
                            
                            beta_tibB = all_beta_tibs2,
                            pval_tibB = all_poob_tibs2,
                            name_strB = exp_name2, # names(all_beta_tibs[[2]])[1]

                            min_dB  = cur_min_dB,
                            cmp_str = cur_cmp_str,
                            is_abs  = cur_is_abs,
                            
                            # [TBD]: Change the output directory to not have replicate,
                            #        use experiment output directory instead...
                            # out_dir    = file.path( opt$out_path, "replicate" ),
                            out_dir    = cur_out_path,
                            run_tag    = exp_name, 
                            reload     = opt$reload,
                            reload_min = 10, 
                            # ret_data   = TRUE,
                            ret_data   = FALSE,
                            parallel   = opt$parallel,
                            
                            vb=vb+2,vt=vt+1,tc=tc+1, tt=tt )
    
    if ( FALSE ) {
      
      cur_dB_tib2 <- cur_dB_tib
      
      cur_dB_tib %>% head() %>% as.data.frame()
      cur_dB_tib %>% dplyr::filter( dB_pas_call ) %>% base::nrow()
      
      cur_dB_tib %>% dplyr::filter( dB_med > 0.1 ) %>% base::nrow() /cur_dB_tib %>% base::nrow()
      # cur_dB_tib2 %>% dplyr::filter( dB_med > 0.1 ) %>% base::nrow() /cur_dB_tib %>% base::nrow()
      
      cur_dB_tib  %>% dplyr::filter( db_pas_per > 70 ) %>% base::nrow()
      cur_dB_tib  %>% dplyr::filter( dB_pas_call ) %>% base::nrow()
      # cur_dB_tib2 %>% dplyr::filter( dB_pas_call ) %>% base::nrow()
      
      den_wS_ggg <- NULL
      den_wS_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=dB_med ) ) +
        ggplot2::geom_density()
      den_wS_ggg

    }
    
    #
    # General Test Case Plotting::
    #
    if ( FALSE ) {
      box_pas_wS_ggg <- NULL
      box_pas_wS_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes(x=db_pas_per,y=wS_per, group = db_pas_per) ) + 
        ggplot2::geom_boxplot( varwidth = TRUE )
      box_pas_wS_ggg
      
      pnt_pas_wS_ggg <- NULL
      pnt_pas_wS_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes(x=db_pas_per,y=wS_per) ) + 
        ggplot2::geom_point()
      pnt_pas_wS_ggg
      
      den_wS_ggg <- NULL
      den_wS_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per ) ) +
        ggplot2::geom_density()
      den_wS_ggg
      
      pnt_wS_avg_ggg <- NULL
      pnt_wS_avg_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_avg ) ) +
        ggplot2::geom_point()
      pnt_wS_avg_ggg
      
      d2d_wS_avg_ggg <- NULL
      d2d_wS_avg_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_avg ) ) +
        ggplot2::geom_density2d()
      d2d_wS_avg_ggg
      
      pnt_wS_med_ggg <- NULL
      pnt_wS_med_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
        ggplot2::geom_point()
      pnt_wS_avg_ggg
      
      d2d_wS_med_ggg <- NULL
      d2d_wS_med_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
        ggplot2::geom_density2d()
      d2d_wS_med_ggg
      
      d2d_wS_med_ggg2 <- NULL
      d2d_wS_med_ggg2 <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
        ggplot2::geom_density2d() +
        ggplot2::facet_grid(
          rows = vars(db_pas_per)
        )
      d2d_wS_med_ggg2
    }
    
    # [TBD]: Need to plot EPIC correlation against db_pas_per & wS_per, which is better
    # cur_dat$ret_tib %>% dplyr::filter( is.na(wS_per) )
    
    # Bind Data
    # [TBD]: Add Annotation data...
    uhm_data_tib <- uhm_data_tib %>% dplyr::bind_rows( 
      cur_dB_tib %>% dplyr::mutate( 
        Remove_Outliers = rm_outliers, 
        Experiment = exp_name ) %>% 
        dplyr::select( Remove_Outliers,Experiment, dplyr::everything() ) )

    if ( opt$single ) break
  }
  if ( opt$single ) break
}

if ( FALSE ) {
  
  # all_data_csv <- file.path( opt$out_path, paste0( opt$run_name,"_all_data.csv.gz" ) )
  # readr::write_csv( x = uhm_data_tib, file = all_data_csv )
  
  all_data_rds <- file.path( opt$out_path, paste0( opt$run_name,"_all_data.rds" ) )
  readr::write_rds( x = uhm_data_tib, file = all_data_rds, compress = "gz" )
  
  exp_auto_ssh_csv <- file.path( opt$out_path, paste0( opt$run_name, "_exp_auto_ssh.csv.gz") )
  readr::write_csv( x = exp_auto_ssh_tib, file = exp_auto_ssh_csv )
  
}

# Basic Plotting
if ( FALSE ) {

  # Pass_Poob_Per
  exp_auto_ssh_sum <- NULL
  exp_auto_ssh_sum <- exp_auto_ssh_tib %>%
    dplyr::group_by( Experiment,Remove_Outliers ) %>%
    dplyr::summarise(
      
      cr_cnt = n(),
      cr_min = min(Pass_Poob_Per),
      cr_avg = mean(Pass_Poob_Per),
      cr_med = median(Pass_Poob_Per),
      cr_sds = sd(Pass_Poob_Per),
      cr_mad = mad(Pass_Poob_Per),
      cr_max = max(Pass_Poob_Per),
      
      .groups = "drop" )
  
  exp_auto_ssh_sum2 <- NULL
  exp_auto_ssh_sum2 <- exp_auto_ssh_tib %>%
    dplyr::filter( Pass_Poob_Per > 85 ) %>%
    dplyr::group_by( Experiment,Remove_Outliers ) %>%
    dplyr::summarise(
      
      cr_cnt = n(),
      cr_min = min(Pass_Poob_Per),
      cr_avg = mean(Pass_Poob_Per),
      cr_med = median(Pass_Poob_Per),
      cr_sds = sd(Pass_Poob_Per),
      cr_mad = mad(Pass_Poob_Per),
      cr_max = max(Pass_Poob_Per),
      
      .groups = "drop" )
  
  uhm_data_sum <- NULL
  uhm_data_sum <- uhm_data_tib %>% 
    print_sum( vec = c("Experiment","Remove_Outliers","dB_pas_call"),
               vb=vb+2,vt=vt+1,tc=tc+1, tt=tt )
  
  den_pas_dB_ggg <- NULL
  den_pas_dB_ggg <- uhm_data_tib %>% 
    # dplyr::mutate( Group_Key = Experiment %>% stringr::str_sub(1,2) ) %>%
    ggplot2::ggplot( aes(x=db_pas_per,
                         # y=wS_per, 
                         group = db_pas_per) ) + 
    ggplot2::geom_density( alpha=0.2 ) +
    ggplot2::facet_grid( rows = vars(Remove_Outliers), 
                         cols = vars(Experiment) )
  den_pas_dB_ggg
  
  
  box_pas_wS_ggg <- NULL
  box_pas_wS_ggg <- uhm_data_tib %>% 
    ggplot2::ggplot( aes(x=db_pas_per,y=wS_per, group = db_pas_per) ) + 
    ggplot2::geom_boxplot( varwidth = TRUE ) +
    ggplot2::facet_grid( rows = vars(Remove_Outliers), 
                         cols = vars(Experiment) )
  box_pas_wS_ggg
  
  # exp_auto_ssh_tib 
  #
  # [TBD]: Figure out how to add the size of each group to the plot
  # [TBD]: Try using the weight aes...
  if ( FALSE ) {
    box_pas_wS_ggg2 <- NULL
    box_pas_wS_ggg2 <- uhm_data_tib %>% 
      ggplot2::ggplot( aes(x=db_pas_per,y=wS_per, group = db_pas_per) ) + 
      ggplot2::geom_boxplot( outlier.alpha = 0.2, varwidth = TRUE ) +
      # notch doesn't really do anything...
                             # notch = TRUE ) +
      ggplot2::facet_grid( rows = vars(Remove_Outliers), 
                           cols = vars(Experiment) )
    
    # This geom_text stuff doesn't work...
      # ggplot2::geom_text(
      #   stat = "db_pas_per" # aes(label = after_stat("db_pas_per") )
      # )
    box_pas_wS_ggg2
  }
}
























if ( FALSE ) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #.                    Summarize and Plot Results::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  rep_data_sum <- NULL
  rep_data_sum <- rep_data_tib %>% 
    print_sum( vec = c("Sample","Rm_Outliers"), vb=vb+3,vt=vt,tc=tc )
  
  if ( par$run_name == "EPICv1" ||
       par$run_name == "COREv1" ||
       par$run_name == "FAILv1" ) {
    
    # rep_data_tib %>% dplyr::filter( Probe_ID %in% epic_ses_dat$mask )
    
    epic_mask_tib <- NULL
    epic_mask_tib <- tibble::tibble( Probe_ID = epic_ses_dat$mask, Masked = TRUE )
    
    rep_mask_tib <- NULL
    rep_mask_tib <- rep_data_tib %>% 
      dplyr::left_join( epic_mask_tib, by=c("Probe_ID") ) %>% 
      dplyr::mutate(
        Masked = dplyr::case_when(
          is.na(Masked) ~ FALSE,
          TRUE ~ TRUE )
      )
    
    rep_mask_sum <- NULL
    rep_mask_sum <- rep_mask_tib %>% 
      dplyr::group_by( Masked ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    
    # Quick Validation::
    #
    # rep_mask_tib %>% 
    #   dplyr::group_by( Masked ) %>% 
    #   dplyr::summarise( Count=n(), .groups = "drop" )
    
    sn <- 0.000001
    
    # Apply Cutoffs/Log
    rep_den_mask_gg <- NULL
    rep_den_mask_gg <- rep_mask_tib %>% 
      dplyr::filter( fin_scr >= 0.9 ) %>%
      # ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
      ggplot2::ggplot( aes( x=log(fin_scr) + sn, color = Sample, fill = Sample ) ) +
      ggplot2::geom_density( alpha = 0.2 )  +
      ggplot2::facet_grid( rows = vars(Masked) )
    
    # Conclusion:: Score does show a difference for Replicate score vs. third
    #  party analysis...
    #
    rep_den_mask_pdf <- file.path( opt$out_path, paste(opt$run_name,"rep_den_mask.pdf", sep=".") )
    rep_den_mask_gg <- NULL
    rep_den_mask_gg <- rep_mask_tib %>% 
      dplyr::filter( fin_scr >= 0.9 ) %>%
      ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
      # ggplot2::ggplot( aes( x=log(fin_scr) + sn, color = Sample, fill = Sample ) ) +
      ggplot2::geom_density( alpha = 0.2 )  +
      ggplot2::facet_grid( rows = vars(Masked),
                           cols = vars(Sample) )
    ggplot2::ggsave( filename = rep_den_mask_pdf, 
                     device = "pdf", width = 7, height = 7, dpi = 320 )
    
  } else {
    
    # Looks good::
    rep_den_pdf <- file.path( opt$out_path, paste(opt$run_name,"rep_density.pdf", sep=".") )
    rep_den_gg <- NULL
    rep_den_gg <- rep_data_tib %>% 
      ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
      ggplot2::geom_density( alpha = 0.2 ) +
      ggplot2::facet_grid( cols = vars(Sample),
                           rows = vars(Rm_Outliers) )
    
    ggplot2::ggsave( filename = rep_den_pdf, 
                     device = "pdf", 
                     width = 7, 
                     height = 7, 
                     dpi = 320 )
    
    #
    # Plot Negative w/w0 Outlier Removal
    #
    negs_sdf_ggg <- negs_sdf_tab %>% 
      ggplot2::ggplot( aes(x=UG, fill=Rm_Outliers) ) + 
      ggplot2::geom_density( alpha=0.2 ) +
      ggplot2::facet_grid( cols = vars(Sample),
                           rows = vars(Rm_Outliers) )
    
    # More or less worthless::
    rep_2den_gg <- NULL
    rep_2den_gg <- rep_data_tib %>% 
      ggplot2::ggplot( aes( x=fin_scr, y=med_dbs ) ) +
      ggplot2::geom_density2d() +
      ggplot2::facet_grid( rows = vars(Sample) )
    
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file

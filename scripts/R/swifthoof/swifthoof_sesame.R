
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                Swifthoof:: 
#                Analyze idats <= Sesame => to pval/beta
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
  base::require("doParallel", quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("sesame", quietly = TRUE) ) )

# test_idats_dir <- "/Users/bretbarnes/Documents/data/idats/idats_AKE-MVP-Failed-v1"
# open_ses_dat <- sesame::openSesame( sesame::searchIDATprefixes( test_idats_dir, recursive = TRUE ), platform = "EPIC" )
# open_ses_ssets <- sesame::openSesame( sesame::searchIDATprefixes( test_idats_dir, recursive = TRUE ), platform = "EPIC" )
# open_ses_dat %>% as.data.frame() %>% tibble::as_tibble( rownames = "Probe_ID" ) %>% dplyr::summarise( dplyr::across( -Probe_ID , ~ sum(is.na(.x)) /865918 ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Set Run Environment:: RStudio/Command-Line
#                            Source All Functions
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par  <- list()
args <- commandArgs(trailingOnly = FALSE)

par$src_path <- NULL
par$run_mode <- args[1]
par$date_str <- Sys.Date() %>% as.character()
par$prgm_dir <- 'swifthoof'
par$prgm_tag <- 'swifthoof_sesame'
par$verbose  <- 3
local_paths  <- c( 
  "/Users/bbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
  "/Users/bretbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
  "/illumina/scratch/darkmatter/tools/Workhorse-Unstained/scripts/R",
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

par$run_name <- "EPIC_v2"
par$run_name <- "GSIBIOINFO-597"
par$run_name <- "GSIBIOINFO-638"
par$run_name <- "Embark_v3"

par$version <- 0
par$version <- 1
# par$version <- 2

opt <- NULL
opt <- swifthoof_sesame_options( pars = par, args = args, vb = par$verbose )
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Pre-processing:: Initialize Parameter Vectors
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NOT DEFINED YET::
#
# pval_vec     <- split_str_to_vec(opt$pval)
# min_pval_vec <- split_str_to_vec(opt$minPval)
# min_perc_vec <- split_str_to_vec(opt$minPerc)
# workflow_vec <- split_str_to_vec(opt$workflow)

gst_manifest_vec <- NULL
gst_manifest_vec <- stringr::str_split( 
  string = opt$genomestudio, 
  pattern = ",", simplify = FALSE ) %>% 
  base::unlist() %>% as.vector()

ses_manifest_vec <- NULL
ses_manifest_vec <- stringr::str_split( 
  string = opt$sesame, 
  pattern = ",", simplify = FALSE ) %>% 
  base::unlist() %>% as.vector()

manifest_key_vec <- NULL
manifest_key_vec <- stringr::str_split( 
  string = opt$manifest_name, 
  pattern = ",", simplify = FALSE ) %>% 
  base::unlist() %>% as.vector()

manifest_key_cnt <- manifest_key_vec %>% length()
manfiest_max_cnt <- max( length(gst_manifest_vec),length(ses_manifest_vec) )
if ( manfiest_max_cnt != manifest_key_cnt )
  cat(glue::glue("{pmssg} {BRK}{RET}",
                 "{pmssg} {S15}   FAILED Manifest Counts NOT equal!{RET}",
                 "{pmssg} {S25} max({manfiest_max_cnt}) != key({manifest_key_cnt}) {RET}",
                 "{pmssg} {BRK}{RET}"))
stopifnot( manfiest_max_cnt == manifest_key_cnt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Pre-processing:: Load Manifests/Controls
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Load Controls to be added::
man_ctls_tib <- NULL
man_ctls_tib <- readr::read_rds( file = opt$controls )
# man_ctls_tib <- readr::read_rds( file.path( opt$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" ) )

# Load both EPIC v1/v2 Analytical Manifests::
#   TBD:: Replace with global manifest and manifest identification in c++
manifest_tibs <- base::list()
if ( opt$run_name %>% stringr::str_starts("GSIBIOINFO-638") ) {
  for ( ii in c(1:manfiest_max_cnt) ) {
    manifest_tibs[[manifest_key_vec[ii]]] <- load_epicv2_manifests( 
      sesame_csv = ses_manifest_vec[ii], 
      genome_csv = gst_manifest_vec[ii],
      name = manifest_key_vec[ii],
      ctls = man_ctls_tib,
      out_dir = opt$out_path, 
      run_tag = opt$run_name, 
      reload = opt$reload, 
      reload_min = 0, 
      parallel = opt$parallel,
      vb=vb,vt=vt+2,tc=tc )
  }
}

# if ( opt$run_name %>% stringr::str_starts("Embark_v3") ) {
#   manifest_tibs <- readr::read_csv( file = opt$sesame[1], show_col_types = FALSE ) %>%
#     dplyr::bind_rows( man_ctls_tib )
# }

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Load Sample Sheets
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# A tibble: 327 × 19 (all ACA sample sheets)
# A tibble: 168 × 21 (unique)
sample_sheet_tib <- NULL

if ( opt$run_name %>% stringr::str_starts("GSIBIOINFO-638") ) {
  sample_sheet_tib <- build_epicv2_sample_sheet( 
    sheet_str = opt$sample_sheets, 
    idat_path = opt$idat_path, 
    sentrix_unq = TRUE,
    out_dir = opt$out_path, 
    run_tag = opt$run_name, 
    reload = opt$reload, 
    reload_min = 10, 
    parallel = opt$parallel, 
    vb=vb,vt=vt+4,tc=tc )
} else if ( opt$run_name %>% stringr::str_starts("Embark_v3") ) {
  idat_tib <- sesame::searchIDATprefixes( 
    dir.name = file.path( opt$idat_path ), 
    recursive = TRUE ) %>% 
    as.data.frame() %>% 
    tibble::as_tibble( rownames = "Sentrix_Name" ) %>% 
    purrr::set_names( c("Sentrix_Name","Prefix") )
  
  sheet_rm_vec <- c("206712840012_R01C01","206712840015_R01C01","206712840015_R04C02","206712840015_R06C02")
  
  sample_sheet_tib <- readr::read_csv( file = opt$sample_sheets[1], show_col_types = FALSE ) %>%
    dplyr::inner_join( idat_tib, by=c("Sentrix_Name") ) %>% 
    dplyr::filter( !Sentrix_Name %in% sheet_rm_vec )
  
  unq_man_key <- sample_sheet_tib$Manifest_Key %>% unique()
  manifest_tibs[[unq_man_key]] <- NULL
  manifest_tibs[[unq_man_key]] <- readr::read_csv( 
    file = opt$sesame[1], show_col_types = FALSE ) %>%
    dplyr::bind_rows( man_ctls_tib )
}
# sample_sheet_tib %>% dplyr::group_by(Sample_Group,Sample_Base) %>% dplyr::summarise( Count=n() ) %>% print(n=1000)

sample_sheet_list <- NULL
sample_sheet_list <- sample_sheet_tib %>% split(.$Sentrix_Name)

manifest_sheets <- NULL
manifest_sheets <- sample_sheet_tib %>% split(.$Manifest_Key )

# Print Important Data::
if ( p8 ) sample_sheet_tib %>% 
  dplyr::select( Sheet_Prep,Sheet_Proc,
                 Sentrix_Name, Sample_Group,
                 Chip_Name,Chip_Version,
                 Sample_Base,Concentration,Sample_Name ) %>% 
  print(n=base::nrow(sample_sheet_tib))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#               Pre-processing:: Initialize All Sesame Pre-Game...
#                            Build Sesame SDF::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

validate_mclapply1()
# library(bettermc)

if ( TRUE ) {
  
  # raw_time <- base::system.time({
  #   all_sdfs2 <- NULL
  #   all_sdfs2 <- file_list( 
  #     path = file.path( opt$out_path,"sdf/raw/prefix_to_sdf"), 
  #     suffix = ".prefix_to_sdf.rds$", 
  #     pattern = ".prefix_to_sdf.rds$" ) %>% 
  #     bettermc::mclapply( readr::read_rds )
  # })
  
  sdfs <- NULL
  load_sdfs_time <- base::system.time({
    rds_list <- file_list( 
      path = file.path( opt$out_path,"sdf/raw/prefix_to_sdf"), 
      suffix = ".prefix_to_sdf.rds$", 
      pattern = ".prefix_to_sdf.rds$" )
    
    sdfs <- NULL
    sdfs <- foreach::foreach( sentrix_name = names(rds_list), 
                              .inorder=TRUE, 
                              .final = function(x) setNames(x, names(rds_list)) ) %dopar% {
                                readr::read_rds( file = rds_list[[sentrix_name]] )
                              }
  })
  tt$addTime( load_sdfs_time, paste( "load_sdfs_time", sep="." ) )
  
} else {
  #  TBD:: Replace prefix_to_sdf() with c++
  if ( opt$parallel ) {
    make_sdfs_time <- base::system.time({
      sdfs <- NULL
      sdfs <- foreach::foreach( sn = names(sample_sheet_list), .inorder=TRUE, .final = function(x) setNames(x, names(sample_sheet_list)) ) %dopar% {
        prefix_to_sdf( prefix   = sample_sheet_list[[sn]][["Prefix"]],
                       platform = sample_sheet_list[[sn]][["Chip_Name"]], 
                       manifest = manifest_tibs[[ sample_sheet_list[[sn]][["Manifest_Key"]] ]],
                       out_dir  = file.path( opt$out_path,"sdf/raw" ),
                       run_tag  = sample_sheet_list[[sn]][["Sentrix_Name"]],
                       reload   = opt$reload, 
                       reload_min = 0, 
                       parallel = opt$parallel, 
                       vb=vb,vt=vt+2,tc=tc,tt=tt )
      }
    })
    tt$addTime( make_sdfs_time, paste( "make_sdfs_time_par", sep="." ) )
    
  } else {
    
    make_sdfs_time <- base::system.time({
      sdfs <- NULL
      sdfs <- sample_sheet_list %>% # head(n=3) %>%
        lapply( function(x) {
          prefix_to_sdf( prefix   = x[["Prefix"]],
                         platform = x[["Chip_Name"]],
                         manifest = manifest_tibs[[ x[["Manifest_Key"]] ]],
                         out_dir  = file.path( opt$out_path,"sdf/raw" ),
                         run_tag  = x[["Sentrix_Name"]],
                         reload   = opt$reload,
                         reload_min = 0,
                         parallel = opt$parallel,
                         vb=vb,vt=vt+2,tc=tc,tt=tt )
        } )
    })
    tt$addTime( make_sdfs_time, paste( "make_sdfs_time_lin", sep="." ) )
    
  }
}

# Ensure we restrict to only Sample Sheet Selected Idats::
sdfs <- sdfs[sample_sheet_tib$Sentrix_Name]

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               Processing::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Loop Over Parameter_Vectors::
#   - Workflow
#   - Thresholds
#   - Sample_Sheets
# TBD:: Write Probe Level Summaries...
# TBD:: Write Summary Plots...
#

#  
# > sesame::prepSesameList()
# code                  func                                description
# 1     0              resetMask                    Reset mask to all FALSE
# 2     Q            qualityMask                 Mask probes of poor design
# 3     G        prefixMaskButCG                    Mask all but cg- probes
# 4     H         prefixMaskButC             Mask all but cg- and ch-probes
# 5     C* inferInfiniumIChannel        Infer channel for Infinium-I probes  [ switch_failed = TRUE ]
# 5     c* inferInfiniumIChannel        Infer channel for Infinium-I probes  [ switch_failed = FALSE ]
# 6     D*             dyeBiasNL           Dye bias correction (non-linear)  [ mask = TRUE ]
# 6     d*             dyeBiasNL           Dye bias correction (non-linear)  [ mask = FALSE ]
# 7     E               dyeBiasL               Dye bias correction (linear)  
# 8     P*                pOOBAH        Detection p-value masking using oob  [ return.pval = FALSE, combine.neg = TRUE,  pval.threshold = min_pval ]
# 8     p*                pOOBAH        Detection p-value masking using oob  [ return.pval = FALSE, combine.neg = FALSE, pval.threshold = min_pval ]
# 9     I            detectionIB Mask detection by intermediate beta values
# 10    B                   noob           Background subtraction using oob  [ combine.neg = TRUE,  offset = 15 ]
# 10    b                   noob           Background subtraction using oob  [ combine.neg = FALSE, offset = 15 ]
# 11    S           inferSpecies                  Set species-specific mask
# 12    T            inferStrain           Set strain-specific mask (mouse)
# 13    M            matchDesign        Match Inf-I/II in beta distribution
#
# 14    N*    detectionPnegEcdf2        Detection p-value masking using neg  [ return.pval = FALSE, pval.threshold = min_pval, use_type = FALSE ]
# 15    V*              getBetas                         Return Beta Values  [ mask = TRUE,  sum.TypeI = FALSE, collapseToPfx = FALSE ]
# 15    v*              getBetas                         Return Beta Values  [ mask = FALSE, sum.TypeI = FALSE, collapseToPfx = FALSE ]
# 15    W*              getBetas                         Return Beta Values  [ mask = TRUE,  sum.TypeI = TRUE,  collapseToPfx = FALSE ]
# 15    w*              getBetas                         Return Beta Values  [ mask = FALSE, sum.TypeI = TRUE,  collapseToPfx = FALSE ]
# 
#

# 
# Example Group Split::
# tmp_list <-  sample_sheet_tib %>% tib_group_split( vec = c("Chip_Name","Chip_Version","Sample_Base","Concentration","Sample_Group"), vb=vb+10,vt=vt,tc=tc )
#

opt$AvgHU_min <- 0.2
opt$MedHU_min <- 0.2

opt$AvgMH_min <- 0.1
opt$MedMH_min <- 0.1

opt$AvgMU_min <- 0.5
opt$MedMU_min <- 0.5


opt$PdP_Avg_Min <- 80
opt$PdP_Med_Min <- 80

opt$Sds_Avg_Min <- 0.10
opt$Sds_Med_Min <- 0.10

opt$Mad_Avg_Min <- 0.10
opt$Mad_Med_Min <- 0.10


opt$single <- 0
# opt$single <- 1
# opt$single <- 2
# opt$single <- 3
# opt$single <- 4

# Workflows:: [Cc][D][NPp][Bb][VW]
# All Combinations::
v1 <- c("C","c","")
v2 <- c("N","")
v3 <- c("D","")
v4 <- c("P","p","")
v5 <- c("B","b")
v6 <- c("V","W")

# Starting::
v1 <- c("C","")
v2 <- c("N","")
v3 <- c("D")
v4 <- c("P")
v5 <- c("B","b")
v6 <- c("V")

# Best Workflow:: CNDPBV
#
# In c++::
#   - Idat() 
#     - [Done]: Calculate DetP with all negative controls
#   - IdatPair()
#     - Duplicate each Infinium I probe for both colors
#       - Add original col column
#     - Mask poor DetP probes
#     - Duplicate failed DetP as Negative
#     - Return Sig_DF data.frame
# Sesame::
#   - Call detectionPneg (N)
#   - DyeBiasNL()
#   - PooBAH()
#   - Noob()
#   - Betas()
#

# work_vec <- rbind( expand.grid( v1,v2,v3,v4,v5,v6 ),
#                    expand.grid( v1,v3,v4,v2,v5,v6 ) )
# work_vec <- rbind( expand.grid( v1,v2,v3,v4,v5 ),
#                    expand.grid( v1,v3,v4,v2,v5 ) )

work_vec <- rbind( expand.grid( v1,v2,v3,v4,v5 ),
                   expand.grid( v1,v2,v4,v3,v5 ) )

work_tib <- NULL
work_tib <- work_vec %>% as.data.frame() %>% 
  tibble::as_tibble() %>% 
  tidyr::unite( work_str, Var1:Var5, sep="" )

poob_mins <- c( 0.05, 0.10, 1.00 )
# poob_mins <- c( 0.05, 1.00 )
negs_mins <- c( 0.05, 1.00 )

for ( negs_min in negs_mins ) {
  negs_str <- paste0( "negs-", stringr::str_pad( 
    string = negs_min * 100, 
    width = 2, side = "left", pad = "0" ) )
  
  for ( poob_min in poob_mins ) {
    poob_str <- paste0( "poob-", stringr::str_pad( 
      string = poob_min * 100, 
      width = 2, side = "left", pad = "0" ) )
    
    for ( work_idx in c(1:base::nrow(work_tib)) ) {
      work_str <- work_tib$work_str[work_idx]
      
      if ( p0 ) cat(glue::glue("{pmssg} {BRK}{RET}",
                               "{pmssg} {S25}     BEG{RET}",
                               "{pmssg} {S15}    Analyzing Work={work_str}, Pval={negs_min}/{poob_min}.{RET}"))
      
      parm_str <- paste( negs_str,poob_str,work_str, sep="." )
      data_dir <- safe_mkdir( file.path( opt$out_path, "params",negs_str,poob_str,work_str) )
      
      sdfs_rds <- file.path( data_dir, paste(parm_str,"sdfs.dat.rds", sep=".") )
      beta_rds <- file.path( data_dir, paste(parm_str,"beta.dat.rds", sep=".") )
      
      sdfs_list <- base::list()
      beta_list <- base::list()
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Load Beta Matrix if Available::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( file.exists(sdfs_rds) && file.exists(beta_rds) &&
           base::file.info( sdfs_rds, extra_cols = TRUE) %>% dplyr::pull(mtime) < 
           base::file.info( beta_rds, extra_cols = TRUE) %>% dplyr::pull(mtime) ) {
        beta_list <- readr::read_rds( file = beta_rds )
      } else if ( file.exists(sdfs_rds) ) {
        sdfs_list <- readr::read_rds( file = sdfs_rds )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #
      #                        Mutate SDFs with Workflow::
      #
      # TBD:: Understand why dopar doesn't work... validate_lapply1()
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # sentrix_vec <- names(sdfs)[ ! names(sdfs) %in% c("206712840012_R01C01","206712840015_R01C01","206712840015_R04C02","206712840015_R06C02") ]
      
      if ( beta_list %>% length() == 0 && sdfs_list %>% length() == 0 ) {
        mut_time <- base::system.time({
          sdfs_list <- foreach::foreach( sn = names(sdfs), .inorder=TRUE, .final = function(x) setNames(x, names(sdfs)) ) %dopar% {
            mutate_sdf_simple( sdf = sdfs[[sn]], 
                               steps = work_str, 
                               negs_min = negs_min,
                               poob_min = poob_min,
                               vb=vb,vt=vt+3,tc=tc,tt=tt )
          }
          if ( sdfs_list %>% length() != 0 ) 
            readr::write_rds( x = sdfs_list, file = sdfs_rds, compress = "gz" )
        })
        tt$addTime( mut_time, paste( "mut",parm_str, sep=".") )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #
      #                         Generate Beta Matrix::
      #
      # Call Betas separately for each manifest? Or pass all Probe_IDs in...
      #   manifest_tibs %>% dplyr::bind_rows() %>% dplyr::distinct( Probe_ID )
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( beta_list %>% length() == 0 ) {
        beta_time <- base::system.time({
          
          for ( man_key in names(manifest_sheets) ) {
            man_sdfs <- NULL
            man_sdfs <- sdfs_list[ manifest_sheets[[man_key]]$Sentrix_Name ]
            # man_sdfs <- sdfs_list[ sentrix_vec ]
            man_cnts <- man_sdfs %>% length()
            
            if ( p0 ) cat(glue::glue("{pmssg} {BRK}{RET}",
                                     "{pmssg} {S25}     BEG{RET}",
                                     "{pmssg} {S15}    Analyzing Work=Betas({man_cnts}), Pval={negs_min}/{poob_min}.{RET}"))
            
            beta_list[[man_key]] <- NULL
            beta_list[[man_key]] <- 
              foreach::foreach( 
                sn = names(man_sdfs), 
                .inorder=TRUE, 
                # .final = function(x) setNames(x, names(man_sdfs)), 
                .combine = cbind ) %dopar% {
                  mutate_sdf_simple( sdf = man_sdfs[[sn]], 
                                     steps = "V",
                                     vb=vb,vt=vt+3,tc=tc,tt=tt )
                }
            # Set Column Names::
            #   - TBD:: Verify that the ordering is correct...
            colnames( beta_list[[man_key]] ) <- names(man_sdfs)
            
            if ( p0 ) cat(glue::glue("{pmssg} {S15}    Analyzing Work=Betas({man_cnts}), Manifest={man_key}.{RET}",
                                     "{pmssg} {S25}     END{RET}",
                                     "{pmssg} {BRK}{RET2}"))
          }
          if ( beta_list %>% length() != 0 )
            readr::write_rds( x = beta_list, file = beta_rds, compress = "gz" )
        })
        tt$addTime( beta_time, paste( "beta",parm_str, sep=".") )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #
      #                           Calculate Stats::
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      opt$calc_stats <- FALSE
      
      if ( opt$calc_stats ) {
        # Need to use tib_group_split() and then just loop over each group and
        #  generate the same stats
        
        grp_col_vec <- NULL
        grp_col_vec <- c( "Sheet_Prep","Sheet_Proc", "Chip_Name","Chip_Version", "Sample_Group", "Sample_Base" )
        
        base_sheets <- NULL
        base_sheets <- sample_sheet_tib %>% tib_group_split( vec = grp_col_vec, sep="." )
        
        # Results::
        #  - order_pass_tib
        #  - Technical Reps: tech_stats_tibs[[sbase]]
        
        order_pass_tib   <- NULL
        titrate_pass_tib <- NULL
        tech_pass_tib    <- NULL
        tech_stats_tibs  <- NULL
        
        for ( grp_key in names(base_sheets) ) {
          # grp_key <- "unk.unk.Embark.3.Charlie"
          # grp_key <- "unk.unk.Embark.3.GD"
          
          cur_sheet <- base_sheets[[grp_key]]
          sheet_cnt <- cur_sheet %>% base::nrow()
          sampl_grp <- cur_sheet %>% dplyr::distinct( Sample_Group ) %>% dplyr::pull( Sample_Group )
          manif_key <- cur_sheet %>% dplyr::distinct( Manifest_Key ) %>% dplyr::pull( Manifest_Key )
          sbase_key <- cur_sheet %>% dplyr::distinct( Sample_Base ) %>% dplyr::pull( Sample_Base )
          
          # if ( p0 ) cat(glue::glue("{pmssg}{TAB} Sheet={grp_key} = {sheet_cnt}.{RET}"))
          
          if ( length(sampl_grp) == 1 && sampl_grp == "TechnicalReplicates" && sheet_cnt > 1 ) {
            #
            # Technical Replicate Analysis::
            #
            beta_mat <- NULL
            beta_mat <- beta_list[[manif_key]][ , cur_sheet$Sentrix_Name ]
            
            # tech_stats_tibs[[sbase_key]] <- NULL
            # tech_stats_tibs[[sbase_key]] <- tibble::tibble(
            tech_stats_tib <- NULL
            tech_stats_tib <- tibble::tibble(
              Probe_ID = beta_mat %>% rownames(),
              Tot = beta_mat %>% ncol(),
              Nan = beta_mat %>% matrixStats::rowCounts( value = NA ),
              PdP = round( 100*(Tot-Nan)/Tot, 2 ),
              Avg = beta_mat %>% matrixStats::rowMeans2( na.rm = TRUE ),
              Sds = beta_mat %>% matrixStats::rowSds( na.rm = TRUE ),
              Med = beta_mat %>% matrixStats::rowMedians( na.rm = TRUE ),
              Mad = beta_mat %>% matrixStats::rowMads( na.rm = TRUE ),
              Var = beta_mat %>% matrixStats::rowVars( na.rm = TRUE ),
            ) %>% 
              dplyr::filter( !Probe_ID %>% stringr::str_starts("^ctl_") ) %>%
              dplyr::select( Probe_ID, PdP,Sds,Mad,Var )
            
            tech_stats_tibs <- tech_stats_tibs %>%
              dplyr::bind_rows( tech_stats_tib )
            
            # if ( is.null(tech_stats_tibs) ) {
            #   tech_stats_tibs <- tech_stats_tib
            # } else {
            #   tech_stats_tibs <- tech_stats_tibs %>%
            #     dplyr::left_join(tech_stats_tib, by="Probe_ID")
            # }
            
          } else if ( length(sampl_grp) == 1 && sampl_grp == "MeTritration" && sheet_cnt > 1 ) {
            #
            # Methylation Titration Analysis::
            #   - Method 1: Get Avg/Med/PPP for each Concentration
            #
            
            # First Sort by Concentration::
            con_sheets <- cur_sheet %>% dplyr::arrange( Concentration ) %>%
              split(.$Concentration)
            
            stats <- NULL
            stats <- base::list()
            con_sheet_cnt <- length(con_sheets)
            for ( ii in c(1:con_sheet_cnt) ) {
              con_key  <- NULL
              con_ssh  <- NULL
              beta_mat <- NULL
              
              con_key  <- names(con_sheets)[ii]
              con_ssh  <- con_sheets[[con_key]]
              beta_mat <- beta_list[[manif_key]][ , con_ssh$Sentrix_Name ]
              
              stats[[con_key]] <- NULL
              stats[[con_key]] <- tibble::tibble(
                Probe_ID = beta_mat %>% rownames(),
                Tot = beta_mat %>% ncol(),
                Nan = beta_mat %>% matrixStats::rowCounts( value = NA ),
                PdP = round( 100*(Tot-Nan)/Tot, 2 ),
                Avg = beta_mat %>% matrixStats::rowMeans2( na.rm = TRUE ),
                Sds = beta_mat %>% matrixStats::rowSds( na.rm = TRUE ),
                Med = beta_mat %>% matrixStats::rowMedians( na.rm = TRUE ),
                Mad = beta_mat %>% matrixStats::rowMads( na.rm = TRUE ),
                Var = beta_mat %>% matrixStats::rowVars( na.rm = TRUE ),
              )
            }
            
            order_stats_tib <- NULL
            order_stats_tib <- tibble::tibble(
              Probe_ID = beta_mat %>% rownames(),
              
              dAvg12 = stats[[2]]$Avg - stats[[1]]$Avg,
              dMed12 = stats[[2]]$Med - stats[[1]]$Med,
              
              dAvg23 = stats[[3]]$Avg - stats[[2]]$Avg,
              dMed23 = stats[[3]]$Med - stats[[2]]$Med,
              
              dAvg34 = stats[[4]]$Avg - stats[[3]]$Avg,
              dMed34 = stats[[4]]$Med - stats[[3]]$Med,
              
              dAvg45 = stats[[5]]$Avg - stats[[4]]$Avg,
              dMed45 = stats[[5]]$Med - stats[[4]]$Med,
              
              dAvg56 = stats[[6]]$Avg - stats[[5]]$Avg,
              dMed56 = stats[[6]]$Med - stats[[5]]$Med,
              
              dAvg67 = stats[[7]]$Avg - stats[[6]]$Avg,
              dMed67 = stats[[7]]$Med - stats[[6]]$Med,
              
            )
            order_stats_mat <- order_stats_tib %>% 
              dplyr::filter( !Probe_ID %>% stringr::str_starts("^ctl_") ) %>% 
              tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()
            order_stats_mat[ which( order_stats_mat < 0 ) ] <- NA_real_
            
            #
            # Count Final Order Failures::
            #
            order_pass_tib <- tibble::tibble(
              Probe_ID = order_stats_mat %>% rownames(),
              
              Tot = order_stats_mat %>% ncol(),
              Nan = order_stats_mat %>% matrixStats::rowCounts( value = NA ),
              Pod = round( 100*(Tot-Nan)/Tot, 2 )
            )
            
            #
            # Build Old School 0%, 50% 100% U,H,M stats::
            #
            titrate_tib <- NULL
            titrate_tib <- tibble::tibble(
              Probe_ID = beta_mat %>% rownames(),
              
              dAvgHU = stats[["50"]]$Avg - stats[["0"]]$Avg,
              dMedHU = stats[["50"]]$Med - stats[["0"]]$Med,
              
              dAvgMH = stats[["100"]]$Avg - stats[["50"]]$Avg,
              dMedMH = stats[["100"]]$Med - stats[["50"]]$Med,
              
              dAvgMU = stats[["100"]]$Avg - stats[["0"]]$Avg,
              dMedMU = stats[["100"]]$Med - stats[["0"]]$Med,
            ) %>% dplyr::mutate( 
              pAvgHU = dAvgHU > opt$AvgHU_min,
              pMedHU = dMedHU > opt$MedHU_min,
              
              pAvgMH = dAvgMH > opt$AvgMH_min,
              pMedMH = dMedMH > opt$MedMH_min,
              
              pAvgMU = dAvgMU > opt$AvgMU_min,
              pMedMU = dMedMU > opt$MedMU_min
            )
            
            titrate_mat <- titrate_tib %>% 
              dplyr::select( - dplyr::starts_with("d") ) %>%
              dplyr::filter( !Probe_ID %>% stringr::str_starts("^ctl_") ) %>% 
              # dplyr::select( dplyr::starts_with("[Pp]") )
              tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()
            titrate_mat[ which( titrate_mat == 0 ) ] <- NA_real_
            
            #
            # Count Final Titration Failures::
            #
            titrate_pass_tib <- titrate_tib %>%
              dplyr::inner_join( 
                tibble::tibble(
                  Probe_ID = titrate_mat %>% rownames(),
                  
                  Tot = titrate_mat %>% ncol(),
                  Nan = titrate_mat %>% matrixStats::rowCounts( value = NA ),
                  Pod = round( 100*(Tot-Nan)/Tot, 2 ) ),
                by=c("Probe_ID") )
            
            
          }
          
          # base_sheets[["Lightning.Auto.EPIC.2.Epigen"]] %>% as.data.frame()
          # base_sheets[["Lightning.Auto.EPIC.2.Epigen"]]$Sentrix_Name
          
          # break
        }
        
        tech_pass_tib <- tech_stats_tibs %>% 
          dplyr::group_by( Probe_ID ) %>% 
          dplyr::summarise( PdP_Avg = mean( PdP, na.rm = TRUE ),
                            PdP_Med = median( PdP, na.rm = TRUE ),
                            
                            Sds_Avg = mean( Sds, na.rm = TRUE ),
                            Sds_Med = median( Sds, na.rm = TRUE ),
                            
                            Mad_Avg = mean( Mad, na.rm = TRUE ),
                            Mad_Med = median( Mad, na.rm = TRUE ),
                            
                            # Var_Avg = mean( Var, na.rm = TRUE ),
                            # Var_Med = median( Var, na.rm = TRUE ),
                            
                            .groups = "drop" ) %>% 
          dplyr::mutate( DetP_Pass = 
                           PdP_Avg > opt$PdP_Avg_Min & 
                           PdP_Med > opt$PdP_Med_Min,
                         
                         Sds_Pass = 
                           Sds_Avg <= opt$Sds_Avg_Min & 
                           Sds_Med <= opt$Sds_Med_Min,
                         
                         Mad_Pass = 
                           Mad_Avg <= opt$Mad_Avg_Min &
                           Mad_Med <= opt$Mad_Med_Min
          )
        
        ord_csv <- file.path( data_dir, paste(parm_str,"ord.sum.csv.gz", sep=".") )
        umh_csv <- file.path( data_dir, paste(parm_str,"umh.sum.csv.gz", sep=".") )
        rep_csv <- file.path( data_dir, paste(parm_str,"rep.sum.csv.gz", sep=".") )
        
        readr::write_csv( x = order_pass_tib,   file = ord_csv )
        readr::write_csv( x = titrate_pass_tib, file = umh_csv )
        readr::write_csv( x = tech_pass_tib,    file = rep_csv )
        
        # order_pass_tib   <- NULL
        # titrate_pass_tib <- NULL
        # tech_pass_tib    <- NULL
        
        # dplyr::select( Probe_ID,PdP_Avg,PdP_Med ) %>% 
        # tech_pass_mat <- tech_pass_tib %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()
        
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                            Workflow Complete::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( p0 ) cat(glue::glue("{pmssg} {S15}    Analyzing Work={work_str}, Pval={negs_min}/{poob_min}.{RET}",
                               "{pmssg} {S25}     END{RET}",
                               "{pmssg} {BRK}{RET2}"))
      
      if (opt$single > 2) break
    }
    if (opt$single > 1) break
  }
  if (opt$single > 0) break
}


opt$summarize <- FALSE
if ( opt$summarize ) {
  
  #
  # Ord::
  #
  ord_cut <- 95
  ord_tab <- file_list( path = file.path( opt$out_path, "params" ), 
                        suffix = "ord.sum.csv.gz$", 
                        pattern = "ord.sum.csv.gz$", 
                        recursive = TRUE ) %>% # head() %>%
    lapply( readr::read_csv, show_col_types = FALSE ) %>%
    dplyr::bind_rows( .id = "Params" )
  
  ord_sum <- ord_tab %>% dplyr::group_by( Params ) %>% 
    dplyr::summarise( Tot = n(),
                      Pas = sum( Pod > ord_cut, na.rm=TRUE ),
                      Per = round( 100 * Pas/Tot, 2 ),
                      .groups = "drop" ) %>%
    dplyr::arrange( -Per )
  
  #
  # umh::
  #
  umh_cut <- 70
  umh_tab <- file_list( path = file.path( opt$out_path, "params" ), 
                        suffix = "umh.sum.csv.gz$", 
                        pattern = "umh.sum.csv.gz$", 
                        recursive = TRUE ) %>% # head() %>%
    lapply( readr::read_csv, show_col_types = FALSE ) %>%
    dplyr::bind_rows( .id = "Params" )
  
  umh_sum <- umh_tab %>% dplyr::group_by( Params ) %>% 
    dplyr::summarise( Tot = n(),
                      Pas = sum( Pod > umh_cut, na.rm=TRUE ),
                      Per = round( 100 * Pas/Tot, 2 ),
                      .groups = "drop" ) %>%
    dplyr::arrange( -Per )
  
  #
  # Rep::
  #
  rep_cut <- 70
  rep_tab <- file_list( path = file.path( opt$out_path, "params" ), 
                        suffix = "rep.sum.csv.gz$", 
                        pattern = "rep.sum.csv.gz$", 
                        recursive = TRUE ) %>% # head() %>%
    lapply( readr::read_csv, show_col_types = FALSE ) %>%
    dplyr::bind_rows( .id = "Params" )
  
  rep_sum <- rep_tab %>% dplyr::group_by( Params ) %>% 
    dplyr::summarise( Tot = n(),
                      Pas = sum ( DetP_Pass & Sds_Pass & Mad_Pass, na.rm = TRUE ),
                      Per = round( 100 * Pas/Tot, 2 ),
                      .groups = "drop" ) %>%
    dplyr::arrange( -Per )
  
  ord_sum %>% head()
  umh_sum %>% head()
  rep_sum %>% head()
  
  # rep_sum %>% dplyr::filter( !Params %>% stringr::str_detect("poob-100") )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Code to be added back in above::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  #
  # Generate Stats::
  #  - Within Chip
  #    - Across All
  #    - Across Sample_Base
  #    - Across Sample_Name
  #  - Across Chips
  #    - r2/dB
  #
  
  # Need to use tib_group_split() and then just loop over each group and
  #  generate the same stats
  for ( chip_name in names(beta_list) ) {
    # TBD:: Calculate column-wise passing detection p-value and remove poor samples...
    #
    
    # Loop Over Prep/Proc
    
    chip_all_sum_tib <- NULL
    chip_all_sum_tib <- tibble::tibble(
      Chip_Name = chip_name,
      Sample_Name = "All",
      Workflow = work_str,
      Negs_Min = negs_min,
      Poob_Min = poob_min,
      
      Probe_ID = beta_list[[chip_name]] %>% rownames(),
      Sample_Tot = beta_list[[chip_name]] %>% nrow(),
      Probe_Tot  = beta_list[[chip_name]] %>% ncol(),
      Probe_Nan  = beta_list[[chip_name]] %>% matrixStats::rowCounts( value = NA ),
      Probe_Pas  = round( 100* (Probe_Tot - Probe_Nan) / Probe_Tot, 2 )
    )
    
    # Base Samples::
    base_sheets <- manifest_sheets[[chip_name]] %>% split(.$Sample_Base)
    
    
    # Exact Samples::
    
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                 OLD STATS::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# > for_time
# user  system elapsed 
# 54.922  11.689  79.725 

# > raw_time ( bettermc::mclapply )
# user  system elapsed 
# 165.571  17.870 124.471 

# > raw_time (linear lapply)
# user  system elapsed 
# 98.274   1.380  99.748 
#
# > raw_time (parallel mclapply)
# user  system elapsed 
# 1.414   2.049  62.909 

# > par_time (foreach %dopar% parallel n=168)
# user  system elapsed 
# 305.704 201.306 761.309
# 
# > par_time (foreach %dopar% parallel n=6)
# user  system elapsed 
# 10.666   4.008  18.988 

# sdfs_list <- NULL
# mut_time <- base::system.time({
#   sdfs_list <- foreach::foreach( sn = names(sdfs), .inorder=TRUE, .final = function(x) setNames(x, names(sdfs)) ) %do% {
#     mutate_sdf_simple( sdf = sdfs[[sn]], 
#                        steps = work_str, 
#                        min_pval = poob_min, 
#                        off_set = 15,
#                        vb=vb,vt=vt+3,tc=tc,tt=tt )
#   }
# })
# > mut_time (foreach %do% linear)
# user   system  elapsed 
# 1067.678   31.658 1098.768 

# Some extra testing::
if ( FALSE ) {
  sesame::controls( tmp_sdfs$`206203800149_R01C01` )
  
  tmp_dat <- sesame::controls( tmp_sdfs$`206203800149_R01C01` ) %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate( Type = Probe_ID %>% stringr::str_to_upper() ) %>%
    dplyr::filter( Type %>% stringr::str_detect("NEGATIVE") )
  
  poob_dat <- sesame::pOOBAH( sdf = tmp_sdfs$`206203800149_R01C01`, return.pval = TRUE )
  negs_dat <- detectionPnegEcdf2( sdf = tmp_sdfs$`206203800149_R01C01`, return.pval = TRUE )
  
  which( poob_dat > 0.05 ) %>% length()
  which( negs_dat > 0.05 ) %>% length()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                 OLD CODE::
#
#                    Pre-processing:: Manifests Comparison
#
# TBD:: Move this somewhere else, maybe once we're done with the full analysis
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  sdfs_list_lin <- sdfs_list
  # for ( ii in c(1:length(sdfs)) ) {
  for ( ii in c(84:length(sdfs)) ) {
    sentrix_name <- names(sdfs)[ii]
    
    if ( p0 ) cat(glue::glue("{pmssg} Mutating[{ii}]='{sentrix_name}'{RET}"))
    
    sdfs_list[[sentrix_name]] <- NULL
    
    if ( sentrix_name == "206712840012_R01C01" ) next
    if ( sentrix_name == "206712840015_R01C01" ) next
    if ( sentrix_name == "206712840015_R04C02" ) next
    if ( sentrix_name == "206712840015_R06C02" ) next
    
    sdfs_list[[sentrix_name]] <- mutate_sdf_simple( 
      sdf = sdfs[[sentrix_name]], 
      steps = work_str, 
      negs_min = negs_min,
      poob_min = poob_min,
      vb=vb,vt=vt+3,tc=tc,tt=tt )
    
    if ( ii %% 20 == 0 ) sdfs_list[[sentrix_name]] %>% head(n=3) %>% print()
  }
  # sesameQC_calcStats( sdf = sdfs[[sentrix_name]], "dyeBias" )@stat$RGdistort
  #
  # Error in if (sesameQC_calcStats(sdf, "dyeBias")@stat$RGdistort > 10) { : 
  #     missing value where TRUE/FALSE needed
  #   Called from: sesame::dyeBiasNL(sdf = sdf, mask = TRUE)
}

# > for_time
# user  system elapsed 
# 408.882  87.611 438.522

# OLD mclapply method
#
# sdfs <- sample_sheet_list %>% # head(n=3) %>% 
#   mclapply( function(x) { 
#     prefix_to_sdf( prefix   = x[["Prefix"]],
#                    platform = x[["Chip_Name"]], 
#                    manifest = manifest_tibs[[ x[["Manifest_Key"]] ]],
#                    out_dir  = file.path( opt$out_path,"sdf/raw" ),
#                    run_tag  = x[["Sentrix_Name"]],
#                    reload   = opt$reload, 
#                    reload_min = 0, 
#                    parallel = opt$parallel, 
#                    vb=vb,vt=vt+2,tc=tc,tt=tt )
#   } )

if ( FALSE ) {
  # A tibble: 865,918 × 47 = 865918
  v1_man_tib <- load_genome_studio_manifest( file = manifest_vec[1], 
                                             load_clean    = TRUE,
                                             load_controls = FALSE,
                                             cols_convert  = FALSE,
                                             write_clean   = FALSE,
                                             overwrite     = FALSE,
                                             ret_data      = FALSE,
                                             vb=vb,vt=vt,tc=tc )
  
  v1_ses_tib <- gs_to_sesame( tib = v1_man_tib, 
                              out_dir = opt$out_path, 
                              run_tag = "EPIC_v1", 
                              reload = opt$reload, 
                              reload_min = 10,
                              parallel = opt$par_csv,
                              vb=vb,vt=vt+4,tc=tc )
  
  # A tibble: 870,635 × 23 = 870635
  v2_man_tib <- load_genome_studio_manifest( file = manifest_vec[2], 
                                             load_clean    = TRUE,
                                             load_controls = FALSE,
                                             cols_convert  = FALSE,
                                             write_clean   = FALSE,
                                             overwrite     = FALSE,
                                             ret_data      = FALSE,
                                             vb=vb,vt=vt,tc=tc )
  
  v2_ses_tib <- gs_to_sesame( tib = v2_man_tib, 
                              out_dir = opt$out_path, 
                              run_tag = "EPIC_v2", 
                              reload = opt$reload, 
                              reload_min = 10,
                              parallel = opt$par_csv,
                              vb=vb,vt=vt+4,tc=tc )
  
  #
  # Compare Manifests::
  #
  # A tibble: 759,087 × 6
  int_man_prb_tib <- dplyr::inner_join(
    v1_man_tib %>% 
      dplyr::select( IlmnID,Name,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>% 
      dplyr::mutate(
        dplyr::across( AlleleB_ProbeSeq, ~tidyr::replace_na(.x, "" ) ),
        Probe_ID = IlmnID %>% stringr::str_remove("[0-9]$")
      ),
    v2_man_tib %>% 
      dplyr::select( IlmnID,Name,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>% 
      dplyr::mutate(
        dplyr::across( AlleleB_ProbeSeq, ~tidyr::replace_na(.x, "" ) ),
        Probe_ID = IlmnID %>% stringr::str_remove("[0-9]$")
      ),
    by=c("AlleleA_ProbeSeq","AlleleB_ProbeSeq"),
    suffix=c("_A","_B") )
  
  # A tibble: 18 × 8
  int_man_prb_tib %>% dplyr::filter( Probe_ID_A != Probe_ID_B )
  
  # A tibble: 759,084 × 9
  int_man_cgn_tib <- dplyr::inner_join(
    v1_man_tib %>% 
      dplyr::select( IlmnID,Name,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>% 
      dplyr::mutate(
        dplyr::across( AlleleB_ProbeSeq, ~tidyr::replace_na(.x, "" ) ),
        Probe_ID = IlmnID %>% stringr::str_remove("[0-9]$")
      ),
    v2_man_tib %>% 
      dplyr::select( IlmnID,Name,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>% 
      dplyr::mutate(
        dplyr::across( AlleleB_ProbeSeq, ~tidyr::replace_na(.x, "" ) ),
        Probe_ID = IlmnID %>% stringr::str_remove("[0-9]$")
      ),
    by=c("Probe_ID"),
    suffix=c("_A","_B") )
  
  # A tibble: 15 × 9
  int_man_cgn_tib %>% dplyr::filter( AlleleA_ProbeSeq_A != AlleleA_ProbeSeq_B )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt,tc=tc,tt=tt )

# sysTime <- Sys.time()
# cat(glue::glue("{pmssg} Finished(time={sysTime}){RET2}"))

# End of file

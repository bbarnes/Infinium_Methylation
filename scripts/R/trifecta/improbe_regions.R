
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Test Scratch Script
#                      Design Probes + Alignment + SNPs::
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
  base::require("microbenchmark", quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("ggplot2", quietly = TRUE) ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Set Run Environment:: RStudio/Command-Line
#                            Source All Functions
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par  <- list()
args <- commandArgs(trailingOnly = FALSE)

par$src_path <- NULL
par$run_mode <- args[1]
par$date_str <- Sys.Date() %>% as.character()
par$prgm_dir <- 'trifecta'
par$prgm_tag <- 'improbe_main'
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
par <- source_functions( pars = par, rcpp = 0, vb = par$verbose )
par <- params_check( pars = par, args = args, 
                     prgm_aux_check = FALSE, vb = par$verbose )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$version <- 1
par$version <- 2
par$version <- 3
par$min_maf <- 0.0001

opt <- NULL
opt <- trifecta_options( pars = par, args = args, vb = par$verbose )
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
prgm_dat <- program_init( name = par$prgm_tag,
                          opts = opt, opt_reqs = opt_reqs,
                          pars = par, par_reqs = par_reqs, 
                          rcpp = 1,
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

strandFR_vec <- NULL
strandCO_vec <- NULL
cpgRank_vec  <- NULL
scrRank_vec  <- NULL
scrType_vec  <- NULL
scrName_vec <-  NULL

strandFR_vec <- opt$strandFR %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
strandCO_vec <- opt$strandCO %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
cpgRank_vec  <- opt$min_cpg_rank %>% str_split(pattern=',', simplify=TRUE) %>% as.integer() %>% as.vector()
scrRank_vec  <- opt$min_scr_rank %>% str_split(pattern=',', simplify=TRUE) %>% as.double() %>% as.vector()
scrType_vec  <- opt$min_prb_score %>% str_split(pattern=',', simplify=TRUE) %>% as.double() %>% as.vector()
scrName_vec  <- opt$score_key %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         0.0.0 Load Pre-Defined::
#                                imGenomes
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

work_list <- 
  file_list( path = opt$workflow,
             dir_only = TRUE,
             paths_exists = FALSE,
             vb = vb )

gen_path <- 
  file_list( path = opt$ref_path, 
             subs = opt$ref_species,
             unique = FALSE,
             dir_only = TRUE, 
             ret_type = COM,
             subs_exists = FALSE,
             vb = vb ) %>%
  file_list( subs = "NCBI", # opt$ref_source,
             unique = FALSE,
             dir_only = TRUE, 
             ret_type = COM,
             subs_exists = FALSE,
             vb = vb ) %>%
  file_list( subs = opt$ref_build,
             unique = FALSE,
             dir_only = TRUE, 
             ret_type = COM,
             subs_exists = FALSE,
             paths_exists = FALSE,
             vb = vb )

dir_add_str <- NULL
# if ( opt$ref_source == "NCBI") dir_add_str <- "Fasta"

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
             vb = vb )

chr_dirs <- 
  file_list( path = gen_path, 
             subs = "Sequence/Chromosomes",
             names = names(ref_seqs),
             unique = FALSE,
             dir_only = TRUE, 
             ret_type = "list",
             subs_exists = FALSE,
             paths_exists = FALSE,
             vb = vb )

opt$ref_builds <- safe_mkdir( file.path( opt$out_path, "Genome_Builds" ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               Load Designs::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

design_key_sym <- rlang::sym( opt$design_key )
design_seq_sym <- rlang::sym( opt$design_seq )
design_prb_sym <- rlang::sym( opt$design_prb )

imp_tib <- safe_read( opt$improbe_path, vb=vb,vt=vt,tc=tc,tt=tt )
imp_cols <- names( imp_tib )

if ( !opt$design_key %in% imp_cols ) {
  fail_mssg <- glue::glue("Column key'{opt$design_key}' does NOT exist in data")
  stop(glue::glue("{perrs} {fail_mssg}!{perrs} Exiting...{RET2}"))
}
if ( !opt$design_seq %in% imp_cols ) {
  fail_mssg <- glue::glue("Column key'{opt$design_seq}' does NOT exist in data")
  stop(glue::glue("{perrs} {fail_mssg}!{perrs} Exiting...{RET2}"))
}

if ( !opt$design_prb %in% imp_cols ) {
  if ( p1 )
    cat( glue::glue("{pmssg} {opt$design_prb} not found in data will use the first.",
                    "two letters of {opt$design_key} {RET2}") )
  imp_tib <- imp_tib %>% 
    dplyr::mutate( !!design_prb_sym := stringr::str_sub( !!design_key_sym, 1,2 ) )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#              Filter Designs and Assign Infinium Design Type::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

min_scr_sym  <- rlang::sym( "Min_Score" )
min_scr1_sym <- rlang::sym( scrName_vec[1] )
if ( length(scrName_vec) == 1 ) {
  imp_tib <- imp_tib %>% dplyr::mutate(
    !!min_scr_sym := !!min_scr1_sym
  )
} else if ( length(scrName_vec) == 2 ) {
  min_scr2_sym <- rlang::sym( scrName_vec[2] )
  
  imp_tib <- imp_tib %>% dplyr::mutate(
    !!min_scr_sym := pmin( !!min_scr1_sym,!!min_scr2_sym )
  )
} else {
  fail_mssg <- glue::glue("Only one or two score columns are allowed: {opt$min_prb_score}!")
  stop(glue::glue("{perrs} {fail_mssg}!{perrs} Exiting...{RET2}"))
}

sel_tib <- NULL
srd_col_sym <- rlang::sym( opt$strand_co_key )
cpg_cnt_sym <- rlang::sym( opt$cpg_cnt_key )
for ( srd_idx in c(1:length(strandCO_vec)) ) {
  min_scr <- scrType_vec[srd_idx]
  srd_key <- strandCO_vec[srd_idx]

  sel_tib <- sel_tib %>% dplyr::bind_rows( 
    dplyr::filter( imp_tib, !!srd_col_sym == srd_key & !!min_scr_sym >= min_scr ) )
}

# sel_tib %>% dplyr::arrange( Seq_ID ) %>% dplyr::select( Seq_ID, opt$strand_co_key, !!min_scr1_sym, min_scr2_sym, min_scr_sym )
sel_tib <- sel_tib %>%
  dplyr::filter( !!srd_col_sym %in% strandCO_vec) %>% 
  dplyr::mutate( Design_Type=dplyr::case_when(
    !!cpg_cnt_sym==cpgRank_vec[1] & !!min_scr_sym>=scrRank_vec[1] ~ 2, # CpgCnt == 3 & Prb_Scr >= 0.6
    !!cpg_cnt_sym==cpgRank_vec[2] & !!min_scr_sym>=scrRank_vec[2] ~ 2, # CpgCnt == 2 & Prb_Scr >= 0.5
    !!cpg_cnt_sym==cpgRank_vec[3] & !!min_scr_sym>=scrRank_vec[3] ~ 2, # CpgCnt == 1 & Prb_Scr >= 0.4
    !!cpg_cnt_sym==cpgRank_vec[4] & !!min_scr_sym>=scrRank_vec[4] ~ 2, # CpgCnt == 0 & Prb_Scr >= 0.3
    TRUE ~ 1 ) %>% as.integer()
  ) %>% dplyr::arrange( Seq_ID )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               Build Designs::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_cpp_tib <- NULL
if ( !is.null(opt$strand_fr_key) && stringr::str_length(opt$strand_fr_key) != 0 &&
     !is.null(opt$strand_co_key) && stringr::str_length(opt$strand_co_key) != 0 &&
     !is.null(opt$strand_tb_key) && stringr::str_length(opt$strand_tb_key) != 0 ) {
  if ( p1 ) cat( glue::glue("{pmssg} Will use strand columns.{RET}",
                            "{pmssg}{TAB} opt$strand_fr_key='{opt$strand_fr_key}{RET}",
                            "{pmssg}{TAB} opt$strand_co_key='{opt$strand_co_key}{RET}",
                            "{pmssg}{TAB} opt$strand_tb_key='{opt$strand_tb_key}{RET}") )

  design_frs_sym <- rlang::sym( opt$strand_fr_key )
  design_cos_sym <- rlang::sym( opt$strand_co_key )
  design_tbs_sym <- rlang::sym( opt$strand_tb_key )
  design_ids_sym <- rlang::sym( "Ilmn_ID" )
  
  sel_tib <- sel_tib %>% 
    dplyr::mutate( 
      !!design_tbs_sym := stringr::str_sub( !!design_tbs_sym, 1,1 ),
      !!design_ids_sym := paste0( !!design_key_sym,"_",
                                  !!design_tbs_sym,
                                  !!design_cos_sym,
                                  Design_Type ) )
  
  if ( opt$unique_cpgs ) sel_tib <- sel_tib %>% 
    dplyr::group_by( Ilmn_ID ) %>%
    dplyr::add_count( Ilmn_ID, name="Rep_Num" ) %>%
    dplyr::mutate( Rep_Idx = dplyr::row_number() ) %>%
    dplyr::ungroup() %>%
    dplyr::filter( Rep_Num == 1 )

  imp_cpp_tib <- 
    improbe_seqs_cpp( 
      fwd_vec_r  = dplyr::select( sel_tib, !!design_seq_sym ) %>% pull( !!design_seq_sym ) %>% as.vector() %>% shear_brac(),
      din_vec_r  = dplyr::select( sel_tib, !!design_prb_sym ) %>% pull( !!design_prb_sym ) %>% as.vector(),
      ids_vec_r  = dplyr::select( sel_tib, !!design_ids_sym ) %>% pull( !!design_ids_sym ) %>% as.vector(),
      frs_vec_r_ = dplyr::select( sel_tib, !!design_frs_sym ) %>% pull( !!design_frs_sym ) %>% as.vector(),
      cos_vec_r_ = dplyr::select( sel_tib, !!design_cos_sym ) %>% pull( !!design_cos_sym ) %>% as.vector(),
      return_source = FALSE, uc = TRUE, 
      vb=vb, vt=vt+10 ) %>% dplyr::bind_rows() %>%
    dplyr::rename( Strand_FR = Strand_SR )
  
} else {
  if ( p1 ) cat( glue::glue("{pmssg} Will build all strand columns.{RET}",
                            "{pmssg}{TAB} opt$strandFR='{opt$strandFR}{RET}",
                            "{pmssg}{TAB} opt$strandCO='{opt$strandCO}{RET}" ) )
  
  imp_cpp_tib <- 
    improbe_seqs_cpp( 
      fwd_vec_r  = dplyr::select( sel_tib, !!design_seq_sym ) %>% pull( !!design_seq_sym ) %>% as.vector() %>% shear_brac(),
      din_vec_r  = dplyr::select( sel_tib, !!design_prb_sym ) %>% pull( !!design_prb_sym ) %>% as.vector(),
      ids_vec_r  = dplyr::select( sel_tib, !!design_key_sym ) %>% pull( !!design_key_sym ) %>% as.vector(),
      frs_str_r  = opt$strandFR %>% stringr::str_remove_all(","),
      cos_str_r  = opt$strandCO %>% stringr::str_remove_all(","),
      return_source = FALSE, uc = TRUE, 
      vb=vb, vt=vt+10 ) %>% dplyr::bind_rows() %>%
    dplyr::rename( Strand_FR = Strand_SR )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   Rejoin improbe_design to improbe_cpp::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_ord_man_csv <- file.path( opt$out_path, paste(opt$run_name, "improbe.order.manifest.csv.gz", sep='.') )
imp_ord_tab_csv <- file.path( opt$out_path, paste(opt$run_name, "improbe.order.stack.csv.gz", sep='.') )

if ( p2 )
  cat( glue::glue("{pmssg} Writing Improbe Order Manifest = {imp_ord_man_csv}...{RET}") )

# readr::write_csv( imp_ord_tib, imp_ord_man_csv )

if ( p2 )
  cat( glue::glue("{pmssg} Done. Writing Improbe Order Manifest = {imp_ord_man_csv}.{RET2}") )

imp_des_tib <- sel_tib %>% 
  dplyr::mutate(
    Coordinate = as.integer(Coordinate),
    Methyl_Allele_TB_Strand = Methyl_Allele_TB_Strand %>% stringr::str_sub(1,1),
    Improbe_48U = UnMethyl_Probe_Sequence %>% stringr::str_sub(2,49) ) %>%
  dplyr::rename(
    Improbe_SeqM = Methyl_Probe_Sequence,
    # Design_Seq = Methyl_Probe_Covered_Top_Sequence,
    Improbe_FR = Methyl_Allele_FR_Strand,
    Improbe_TB = Methyl_Allele_TB_Strand,
    Improbe_CO = Methyl_Allele_CO_Strand,
    Improbe_ScoreM = Methyl_Final_Score,
    Improbe_Cpg_Count = Methyl_Underlying_CpG_Count,
    # Improbe_Next_Base = Methyl_Next_Base,
    Improbe_SeqU = UnMethyl_Probe_Sequence,
    Improbe_ScoreU = UnMethyl_Final_Score
  ) %>% dplyr::inner_join( imp_cpp_tib, by=c("Ilmn_ID"="Probe_ID") ) %>%
  dplyr::mutate(
    Align_Seq_A = dplyr::case_when(
      Design_Type == 1 ~ PRB_U,
      Design_Type == 2 ~ PRB_D,
      TRUE ~ NA_character_
    ),
    Align_Seq_B = dplyr::case_when(
      Design_Type == 1 ~ PRB_M,
      Design_Type == 2 ~ NA_character_,
      TRUE ~ NA_character_
    )
  )

imp_ord_tab <- imp_des_tib %>%
  dplyr::distinct( Ilmn_ID, Chromosome, Coordinate, 
                   Strand_FR, Improbe_TB, Strand_CO,
                   Align_Seq_A, Align_Seq_B ) %>% 
  mutate_all( list(~na_if(.,"") ) ) %>%
  tidyr::pivot_longer( cols = c(Align_Seq_A, Align_Seq_B),
                       names_to = "Probe_Allele",
                       names_prefix = "Align_Seq_",
                       values_to = "Align_Fwd_Seq",
                       values_drop_na = TRUE ) %>%
  dplyr::mutate( 
    Align_Rev_Seq = mutate_seqs_cpp( seqs = Align_Fwd_Seq, m='RC', uc = TRUE )
  ) %>%
  dplyr::arrange( Ilmn_ID, Probe_Allele, Chromosome, Coordinate ) %>%
  dplyr::group_by( Ilmn_ID, Probe_Allele, Align_Fwd_Seq ) %>% 
  dplyr::mutate( Improbe_Map_Count = dplyr::row_number() ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select( Ilmn_ID, Improbe_Map_Count, Probe_Allele, 
                 Chromosome, Coordinate, 
                 Align_Fwd_Seq, Align_Rev_Seq )

if ( p2 )
  cat( glue::glue("{pmssg} Writing Improbe Order Stack = {imp_ord_tab_csv}...{RET}") )

# readr::write_csv( imp_ord_tab, imp_ord_tab_csv )

if ( p2 )
  cat( glue::glue("{pmssg} Done. Writing Improbe Order Stack = {imp_ord_tab_csv}.{RET2}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Write Passing Probe Fasta Files::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$bsp_path <- safe_mkdir( file.path( opt$out_path, "bsmap" ) )
opt$ord_fas  <- file.path( opt$bsp_path, paste(opt$run_name,"prb.fa.gz", sep='.') )

ord_fas_vec <- imp_ord_tab %>% 
  dplyr::filter(!is.na( Align_Fwd_Seq ) ) %>%
  dplyr::distinct( Ilmn_ID, Align_Fwd_Seq, Probe_Allele ) %>%
  dplyr::mutate( Fasta_Str = paste0(">",Ilmn_ID,"_",Probe_Allele,"\n",Align_Fwd_Seq) ) %>% 
  dplyr::arrange( Fasta_Str ) %>%
  dplyr::pull( Fasta_Str )

# ord_fas_vec <- dplyr::bind_rows(
#   imp_ord_tab %>% 
#     dplyr::filter(!is.na( PRB_U ) ) %>%
#     dplyr::distinct( Probe_ID, PRB_U ) %>%
#     dplyr::mutate( Fasta_Str = paste0(">",Probe_ID,"_A\n",PRB_U) ),
#   
#   imp_ord_tab %>% 
#     dplyr::filter(!is.na( PRB_M ) ) %>%
#     dplyr::distinct( Probe_ID, PRB_M ) %>%
#     dplyr::mutate( Fasta_Str = paste0(">",Probe_ID,"_B\n",PRB_M) )
# ) %>% 
#   dplyr::arrange( Fasta_Str ) %>%
#   dplyr::pull( Fasta_Str )

readr::write_lines( x = ord_fas_vec, file = opt$ord_fas )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Launch BSMAP on Fasta Files::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Chicken  = ref_seqs$galGal5
# CrawFish = ref_seqs$pvir_pacbio_dovetail
# Cho      = ref_seqs$criGriChoV1
opt$bsp_tsv <- run_bsmap( ref_fas = ref_seqs$pvir_pacbio_dovetail,
                          can_fas = opt$ord_fas, 
                          bsp_exe = opt$bsmap_exe, 
                          bsp_dir = opt$bsmap_dir, 
                          slim = TRUE, 
                          out_dir = opt$bsp_path,
                          run_tag = opt$run_name,
                          reload = opt$reload,
                          reload_min = 10,
                          vb=vb,vt=vt,tc=tc,tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Load BSMAP Alignment Results::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ord_bsp_tib <- load_bsmap( file = opt$bsp_tsv,
                           sort = FALSE,
                           add_cnt = TRUE,
                           parse_id = TRUE,
                           out_dir = opt$bsp_path,
                           run_tag = opt$run_name,
                           reload = opt$reload,
                           reload_min = 10,
                           vb=vb,vt=vt,tc=tc,tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Consolidate Alignment/Improbe Results:: UM
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Report on all alignments:: Not just Unique Matches
#
tags <- c("UM", "MA", "OF")
ord_bsp_list <- ord_bsp_tib %>% split(f = .$Tag)
bsp_imp_list <- list()
for ( tag in tags ) {
  if ( vb >= 3 )
    cat( glue::glue("{pmssg} Joining tag='{tag}'{RET}") )
  
  bsp_imp_list[[tag]] <- dplyr::bind_rows(
    ord_bsp_list[[tag]] %>% 
      dplyr::filter( Strand_BS=="++" | Strand_BS=="--" ) %>%
      dplyr::left_join( imp_ord_tab %>% dplyr::select( -Align_Rev_Seq ), 
                        by=c("Ilmn_ID","Allele"="Probe_Allele",
                             "Aln_Seq"="Align_Fwd_Seq" ) ),
    
    ord_bsp_list[[tag]] %>% 
      dplyr::filter( Strand_BS=="+-" | Strand_BS=="-+" ) %>%
      dplyr::left_join( imp_ord_tab %>% dplyr::select( -Align_Fwd_Seq ),
                        by=c("Ilmn_ID","Allele"="Probe_Allele",
                             "Aln_Seq"="Align_Rev_Seq" ) )
  ) %>% 
    dplyr::mutate(
      Bsp_Imp_Delta = dplyr::case_when(
        # Chromosome != Chr ~ 99999999.0,
        TRUE ~ as.double( CpG_Pos - Coordinate )
      ) %>% as.integer()
    ) %>%
    dplyr::arrange(Ilmn_ID, Allele, Bsp_Imp_Delta) %>%
    dplyr::distinct( Ilmn_ID, Allele, Map_Count, .keep_all = TRUE ) %>%
    dplyr::filter( !is.na(Bsp_Imp_Delta) ) %>%
    dplyr::filter( !is.na(Mis_Cnt) ) %>%
    dplyr::filter( !is.na(Ins_Len) ) %>%
    dplyr::filter( Mis_Cnt == 0 ) %>%
    dplyr::filter( base::abs(Bsp_Imp_Delta) <= 1 )
  
  cur_cnt <- bsp_imp_list[[tag]] %>% base::nrow()
  
  if ( vb >= 3 )
    cat( glue::glue("{pmssg} Done. Joining tag='{tag}'; Count={cur_cnt}.{RET2}") )
  
  break
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#               Consolidate Infinium I U/M Back Together::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

unq_col_vec <- c("Ilmn_ID", "Probe_ID", "Probe_Type", "Infinium_Design", 
                 "Strand_TB", "Strand_CO", "Improbe_Map_Count", # "Design_Seq",
                 "Chr", "Chromosome", "Coordinate", "Beg", "Ref_Seq" )

inf1_tib <- dplyr::inner_join(
  bsp_imp_list[["UM"]] %>% dplyr::filter(Infinium_Design == 1 & Allele == 'A'),
  bsp_imp_list[["UM"]] %>% dplyr::filter(Infinium_Design == 1 & Allele == 'B'),
  by = unq_col_vec,
  suffix = c("_A", "_B")
) %>% dplyr::distinct(Ilmn_ID, .keep_all = TRUE)

inf2_tib <- 
  bsp_imp_list[["UM"]] %>% 
  dplyr::filter(Infinium_Design == 2 & Allele == 'A') %>% 
  dplyr::distinct(Ilmn_ID, .keep_all = TRUE)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Recombine Infinium I/II into Manifest::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Add back alignment stats
#
sel_ids_vec <- c( dplyr::bind_rows(
  inf1_tib %>% dplyr::select(Ilmn_ID), 
  inf2_tib %>% dplyr::select(Ilmn_ID) ) %>% 
    dplyr::distinct() %>% dplyr::pull(Ilmn_ID) )

unq_sel_tib <- imp_des_tib %>% 
  dplyr::mutate(Strand_FR = dplyr::case_when(
    Improbe_FR == "F" ~ '+',
    Improbe_FR == "R" ~ '-',
    TRUE ~ NA_character_ )
  ) %>%
  dplyr::filter( Ilmn_ID %in% sel_ids_vec ) %>%
  dplyr::add_count( Ilmn_ID, name="Rep_Count" ) %>%
  dplyr::group_by( Ilmn_ID ) %>%
  dplyr::mutate( 
    Rep_Rank = dplyr::row_number(),
    Ilmn_ID = paste0( Ilmn_ID, Rep_Rank )
  ) %>% dplyr::ungroup()

unq_sel_sum <- unq_sel_tib %>% 
  dplyr::group_by( Design_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=100)

unq_sel_grs <- 
  GenomicRanges::GRanges(
    seqnames = Rle( paste0("chr",unq_sel_tib$Chromosome ) ),
    strand=Rle(unq_sel_tib$Strand_FR),
    
    Ilmn_ID  = unq_sel_tib$Ilmn_ID,
    
    Probe_Type = unq_sel_tib$Probe_Type,
    Improbe_TB = unq_sel_tib$Improbe_TB,
    Improbe_CO = unq_sel_tib$Improbe_CO,
    Infinium_Design = unq_sel_tib$Design_Type,
    
    Rep_Count = unq_sel_tib$Rep_Count,
    Rep_Rank = unq_sel_tib$Rep_Rank,
    
    # Target_Chip = unq_sel_tib$Platform,
    # Bead_Pool = unq_sel_tib$Bead_Pool,
    
    # Gen_Dup_Count = unq_sel_tib$Gen_Dup_Count,
    # Seq_Dup_Count = unq_sel_tib$Seq_Dup_Count,
    # Ids_Dup_Count = unq_sel_tib$Ids_Dup_Count,
    
    # MASK_general=unq_sel_tib$MASK_general,
    # MASK_mapping=unq_sel_tib$MASK_mapping,
    # MASK_typeINextBaseSwitch=unq_sel_tib$MASK_typeINextBaseSwitch,
    
    IRanges(start = unq_sel_tib$Coordinate,
            width = 2,
            names=paste0( unq_sel_tib$Ilmn_ID ) )
  )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Format Order File::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

raw_out_csv <- file.path( opt$out_path, "raw_manifest.csv.gz")
raw_out_tib <- unq_sel_tib %>% 
  dplyr::distinct(Ilmn_ID, .keep_all = TRUE ) %>% 
  dplyr::rename( Infinium_Design=Design_Type,
                 Probe_Seq_A=Align_Seq_A,
                 Probe_Seq_B=Align_Seq_B,
                 Improbe_Next_Base=NXB_D) %>%
  dplyr::select(
    Ilmn_ID, # Loci_ID, 
    Probe_Type, Infinium_Design, Probe_Seq_A, Probe_Seq_B,
    Chromosome, Coordinate, Improbe_FR, Improbe_TB, Improbe_CO, Improbe_Next_Base,
    Improbe_ScoreM, Improbe_ScoreU, Improbe_Cpg_Count, 
    Forward_Sequence, Min_Score, Improbe_Cpg_Count)

if ( opt$pick_best ) raw_out_tib <- raw_out_tib %>%
  dplyr::arrange( Improbe_CO, -Min_Score, Improbe_Cpg_Count ) %>% 
  dplyr::distinct( Ilmn_ID, .keep_all=TRUE ) %>%
  dplyr::arrange( Chromosome, Coordinate )

safe_write( x = raw_out_tib, file = raw_out_csv, vb=vb )

# raw_out_tib <- unq_sel_tib %>% 
#   dplyr::distinct(Ilmn_ID, .keep_all = TRUE ) %>% 
#   dplyr::select(
#     Ilmn_ID, Loci_ID, 
#     Probe_Type, Infinium_Design, Probe_Seq_A, Probe_Seq_B,
#     Chromosome, Coordinate, Improbe_FR, Improbe_TB, Improbe_CO, Improbe_Next_Base,
#     Improbe_ScoreM, Improbe_ScoreU, Improbe_Cpg_Count, 
#     Platform, Bead_Pool, 
#     Forward_Sequence, Top_Sequence, Design_Seq, 
#     Rep_Count, Rep_Rank,
#     Gen_Dup_Count, Seq_Dup_Count, Ids_Dup_Count )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file

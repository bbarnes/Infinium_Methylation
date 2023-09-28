
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                              triplecrown:: 
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
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

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
par$prgm_tag <- 'trifecta_genomic'
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
par <- source_functions( pars = par, rcpp = FALSE, vb = 10 )
par <- params_check( pars = par, args = args, prgm_aux_check = FALSE, vb = 10 )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt <- NULL
opt <- trifecta_genomic_options( pars = par, args = args, vb = 10 )
vb  <- opt$verbose

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
                          rcpp = 0,
                          vb = opt$verbose, vt=3, tc=0 )

opt <- prgm_dat$opt
par <- prgm_dat$par
opt_tib <- prgm_dat$opt_tib
par_tib <- prgm_dat$par_tib

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Initialize Run Objects
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new()

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

# Load cgn <=> top truth map::
tru_tib <- safe_read( file = opt$imap_csv, vb = opt$verbose, tt = pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         0.1.1 Load Pre-Defined::
#                                imGenomes
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

work_list <-
  file_list( path = opt$workflow,
             dir_only = TRUE,
             paths_exists = FALSE,
             vb = opt$verbose )

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

dir_add_str <- NULL
if ( opt$ref_source == "NCBI") dir_add_str <- "Fasta"

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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         0.1.2 Load User Inputs::::
#                                  HM450
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# opt$tar_manifest <- "/Users/bretbarnes/Documents/data/manifests/methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         0.1.2 Load User Inputs::::
#                                Chicago
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Scratch for Chicago::
#
prb_csv <- "/Users/bretbarnes/Documents/data/CustomContent/Chicago-Ober-Custom/Content_Order_Two/selected/top-selected-supplementary-pool-chicago.csv.gz"
sel_csv <- "/Users/bretbarnes/Documents/data/CustomContent/Chicago-Ober-Custom/Content_Order_Two/selected/Custom_Methylation_Array_V2_brief_Selected_from_Chicago.csv.gz"

sel_tib <- safe_read(sel_csv, vb=3)
prb_tib <- safe_read(prb_csv, vb=3)

sel_ids_tib <- prb_tib %>% 
  # dplyr::inner_join(sel_tib) %>% 
  dplyr::rename(Ilmn_ID = Probe_ID) %>%
  tidyr::separate( Ilmn_ID, 
                   into = c("Probe_ID", "SRD_Str"), 
                   sep="_", remove = FALSE ) %>%
  tidyr::separate( SRD_Str, 
                   into = c("Strand_TB", "Strand_CO", "Infinium_Design"), 
                   sep = c(1,2,3), convert = TRUE ) %>% 
  dplyr::mutate(
    Bead_Count = dplyr::case_when(
      Infinium_Design == 1 ~ 2,
      Infinium_Design == 2 ~ 1,
      TRUE ~ NA_real_
    ) %>% as.integer()
  )

sel_ids_sum <- sel_ids_tib %>% 
  dplyr::group_by( Strand_TB, Strand_CO, Infinium_Design ) %>% 
  dplyr::summarise( Count=n(), 
                    Bead_Sum = sum(Bead_Count ), 
                    .groups = "drop" )

sel_ids_sum %>% 
  dplyr::group_by( Infinium_Design ) %>%
  dplyr::summarise( Inf1_Cnt = sum(Count), 
                    Inf2_Cnt = sum(Count),
                    .groups = "drop" )

sel_cnt <- sel_ids_sum$Bead_Sum %>% sum()

tar_des_tsv <- "/Users/bretbarnes/Documents/data/CustomContent/Chicago-Ober-Custom/Content_Order_Two/2021-08-23_Custom_Methylation_Array_Additions_v2_all-designs.tsv.gz"
all_des_tib <- safe_read(tar_des_tsv, vb = opt$verbose, tt = pTracker) %>% 
  dplyr::rename( 
    Probe_ID = Seq_ID,
    Strand_FR = Methyl_Allele_FR_Strand,
    Strand_TB = Methyl_Allele_TB_Strand,
    Strand_CO = Methyl_Allele_CO_Strand ) %>%
  dplyr::mutate(
    Min_Score = base::round(10 * pmin(Methyl_Final_Score, UnMethyl_Final_Score), 0 ),
    Infinium_Design = dplyr::case_when(
      Min_Score < 3 & Strand_CO == "C" ~ 0,
      Min_Score < 2 & Strand_CO == "O" ~ 0,
      
      Methyl_Underlying_CpG_Count > 3 ~ 1,
      
      Min_Score < 4 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 1 ~ 1,
      Min_Score < 5 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 2 ~ 1,
      Min_Score < 6 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 3 ~ 1,
      
      Min_Score < 4 & Strand_CO == "O" & Methyl_Underlying_CpG_Count <= 1 ~ 1,
      Min_Score < 5 & Strand_CO == "O" & Methyl_Underlying_CpG_Count <= 2 ~ 1,
      Min_Score < 6 & Strand_CO == "O" & Methyl_Underlying_CpG_Count <= 3 ~ 1,
      
      Min_Score >= 3 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 0 ~ 2,
      Min_Score >= 4 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 1 ~ 2,
      Min_Score >= 5 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 2 ~ 2,
      Min_Score >= 6 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 3 ~ 2,
      
      Min_Score >= 2 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 0 ~ 2,
      Min_Score >= 4 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 1 ~ 2,
      Min_Score >= 5 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 2 ~ 2,
      Min_Score >= 6 & Strand_CO == "C" & Methyl_Underlying_CpG_Count <= 3 ~ 2,
      
      TRUE ~ -1
    ) %>% as.integer(),
    Strand_TB = stringr::str_sub(Strand_TB, 1,1),
    Probe_Type = stringr::str_sub(Probe_ID, 1,2),
    Ilmn_ID = paste0( Probe_ID,"_",
                      Strand_TB,
                      Strand_CO,
                      Infinium_Design )
  )

tar_des_tib <- all_des_tib %>% 
  dplyr::inner_join(sel_ids_tib, 
                    byc=c("Ilmn_ID", "Probe_ID", "Chromosome", "Coordinate", 
                          "Strand_TB", "Strand_CO", "Infinium_Design") ) %>%
  dplyr::select( Ilmn_ID, Probe_ID, Probe_Type, Strand_FR, Strand_TB, Strand_CO, 
                 Infinium_Design, Chromosome,Coordinate, Bead_Count,
                 dplyr::everything() )

r_improbe_dir <- file.path( opt$out_path, "r_improbe")
safe_mkdir( r_improbe_dir )

# r_imp_tib <- r_improbe( tib = tar_des_tib, 
#                         ids_key = "Ilmn_ID", 
#                         seq_key = "Forward_Sequence", 
#                         din_key = "Probe_Type", 
#                         out_dir = r_improbe_dir, 
#                         run_tag = opt$run_name,
#                         add_matseq = TRUE, 
#                         vb = vb )
# 
# r_sel_tib <- r_imp_tib %>% 
#   dplyr::inner_join(tar_des_tib, 
#                     by=c( # "Ilmn_ID", 
#                          "Prb_1M"="Methyl_Probe_Sequence",
#                          "Prb_1U"="UnMethyl_Probe_Sequence") )

#
# This version doesn't work::
#

# source("/Users/bretbarnes/Documents/tools/backup/Infinium_Methylation_Workhorse.scripts/25072021/R/probe_design/functions/r_improbe_functions.R")
# 
# r2_imp_tib <- r_improbe( tib = tar_des_tib, 
#                          ids_key = "Ilmn_ID", 
#                          seq_key = "Forward_Sequence", 
#                          din_key = "Probe_Type", 
#                          # out_dir = r_improbe_dir, 
#                          # run_tag = opt$run_name,
#                          add_matseq = TRUE, 
#                          verbose = vb )

#
# This version is the old school version that works!!!
#
source("/Users/bretbarnes/Documents/tools/backup/Infinium_Methylation_Workhorse.04012020/scripts/R/probe_design/functions/improbe_functions.R")

r3_imp_tib <- desSeq_to_prbs( tib = tar_des_tib,
                              idsKey = "Ilmn_ID", 
                              seqKey = "Forward_Sequence",
                              prbKey = "Probe_Type", 
                              # out_dir = r_improbe_dir, 
                              # run_tag = opt$run_name,
                              addMatSeq = TRUE, 
                              verbose = vb )

dat_csv <- file.path( opt$out_path, "Custom_Methylation_Array_Chicago_v1.1.dat.csv.gz" )
ord_csv <- file.path( opt$out_path, "Custom_Methylation_Array_Chicago_v1.1.ord.csv.gz" )
dat_tib <- r3_imp_tib %>% 
  dplyr::inner_join( tar_des_tib, 
                     by=c( "Ilmn_ID", 
                           "PRB1_M_MAT"="Methyl_Probe_Sequence",
                           "PRB1_U_MAT"="UnMethyl_Probe_Sequence",
                           "Forward_Sequence") ) %>%
  dplyr::mutate( 
    Assay_Design_Id=Ilmn_ID,
    AlleleA_Probe_Id = paste(Assay_Design_Id,"A", sep="_"),
    AlleleB_Probe_Id = dplyr::case_when(
      Infinium_Design == 1 ~ paste(Assay_Design_Id,"B", sep="_"),
      Infinium_Design == 2 ~ "",
      TRUE ~ NA_character_ ),
    AlleleA_Probe_Sequence = dplyr::case_when(
      Infinium_Design == 1 ~ PRB1_U_MAT,
      Infinium_Design == 2 ~ PRB2_D_MAT,
      TRUE ~ NA_character_ ),
    AlleleB_Probe_Sequence = dplyr::case_when(
      Infinium_Design == 1 ~ PRB1_M_MAT,
      Infinium_Design == 2 ~ "",
      TRUE ~ NA_character_ ),
    Next_Base = stringr::str_to_upper(NXB_M),
    Normalization_Bin = dplyr::case_when(
      Infinium_Design == 1 & Next_Base == "A" ~ "A",
      Infinium_Design == 1 & Next_Base == "T" ~ "A",
      Infinium_Design == 1 & Next_Base == "C" ~ "B",
      Infinium_Design == 1 & Next_Base == "G" ~ "B",
      Infinium_Design == 2 ~ "C",
      TRUE ~ NA_character_
    )
  )

readr::write_csv( dat_tib, dat_csv )

ord_tib <- dat_tib %>%
  dplyr::select( Assay_Design_Id,
                 AlleleA_Probe_Id,AlleleA_Probe_Sequence,
                 AlleleB_Probe_Id,AlleleB_Probe_Sequence,
                 Normalization_Bin )

readr::write_csv( ord_tib, ord_csv )


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Fasta for Methyl Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tar_bed_tib <- tar_des_tib %>% dplyr::distinct(Seq_ID,Chromosome,Coordinate)

outM_fas <- file.path(opt$out_path, "methyl.infinium-1.fas.gz")
fasM_vec <- tar_des_tib %>% 
  dplyr::filter( Infinium_Design > 0 ) %>%
  dplyr::mutate( 
    Fas_Str = paste0( ">",Probe_ID,"_",Min_Score,"_",Methyl_Underlying_CpG_Count,"\n",
                      Methyl_Probe_Sequence )
  ) %>%
  dplyr::pull(Fas_Str)
readr::write_lines(fasM_vec, outM_fas)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Fasta for UnMethyl Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

outU_fas <- file.path(opt$out_path, "unmethyl.infinium-1.fas.gz")
fasU_vec <- tar_des_tib %>% 
  dplyr::filter( Infinium_Design > 0 ) %>%
  dplyr::mutate( 
    Fas_Str = paste0( ">",Probe_ID,"_",Min_Score,"_",Methyl_Underlying_CpG_Count,"\n",
                      UnMethyl_Probe_Sequence )
  ) %>%
  dplyr::pull(Fas_Str)
readr::write_lines(fasU_vec, outU_fas)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Score Summary Stats::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tar_des_tib %>% 
  dplyr::filter( Infinium_Design > 0 ) %>%
  dplyr::select( Probe_ID, Min_Score,Methyl_Underlying_CpG_Count, 
                 Methyl_Probe_Sequence, UnMethyl_Probe_Sequence )

tar_sum_tib <- tar_des_tib %>%
  dplyr::group_by(Min_Score, Methyl_Underlying_CpG_Count, Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups = "drop") %>% print(n=1000)

tar_des_tib %>% 
  dplyr::filter(Infinium_Design > 0) %>%
  dplyr::group_by(Min_Score, Methyl_Underlying_CpG_Count, Infinium_Design, Strand_CO) %>% 
  dplyr::summarise(Count=n(), .groups = "drop") %>% print(n=1000)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Build Command Script::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

slim <- TRUE

if (slim) {
  bsp_cols <- cols(
    Probe_Id = col_character(),
    Aln_Seq  = col_character(),
    
    Tag  = col_character(),
    Chr  = col_character(),
    Beg  = col_integer(),
    Strand_BS  = col_character(),
    
    Mis_Cnt  = col_integer(),
    Ref_Seq  = col_character(),
    Gap_Cnt  = col_integer(),
    Mis_Str  = col_character()
  )
  chr_key <- names(bsp_cols$cols)[4]
  beg_key <- names(bsp_cols$cols)[5]
  
} else {
  bsp_cols <- cols(
    Probe_Id = col_character(),
    Aln_Seq  = col_character(),
    Qual     = col_character(),
    
    Tag  = col_character(),
    Chr  = col_character(),
    Beg  = col_integer(),
    Strand_BS  = col_character(),
    
    Mis_Cnt  = col_integer(),
    Ref_Seq  = col_character(),
    Gap_Cnt  = col_integer(),
    Mis_Str  = col_character()
  )
  chr_key <- names(bsp_cols$cols)[5]
  beg_key <- names(bsp_cols$cols)[6]
  
}

out_key <- "unmethyl"
ord_fas <- outU_fas

out_key <- "methyl"
ord_fas <- outM_fas

bsp_opt <- "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
bsp_ssh <- file.path(opt$out_path, paste(out_key,"bsmap.sh", sep='.') )
beg_tsv <- file.path(opt$out_path, paste(out_key,"beg.txt", sep='.') )
end_tsv <- file.path(opt$out_path, paste(out_key,"end.txt", sep='.') )
out_bsp <- file.path(opt$out_path, paste(out_key,"bsmap.bsp", sep='.') )
out_tsv <- file.path(opt$out_path, paste(out_key,"bsmap.tsv.gz", sep='.') )
ref_fas <- "/Users/bretbarnes/Documents/data/imGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta/GRCh37.genome.fa.gz"
bsp_dir <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz"
bsp_exe <- "bsmapz"

add_cmd <- ""
if (slim) add_cmd <- "cut -f 1,2,4-11 | "

if (!is.null(bsp_dir)) {
  bsp_exe_test <- file.path( bsp_dir, base::basename(bsp_exe) )
  if (file.exists(bsp_exe_test)) bsp_exe_test
}

# TBD:: Clean up this missing bsp_exe for docker...
#
if ( ( is.null(bsp_exe) || !file.exists(bsp_exe ) ) &&
     (!is.null(bsp_dir) &&  dir.exists(bsp_dir  ) ) )
  bsp_exe <- file.path(bsp_dir, bsp_exe)

if ( is.null(bsp_exe) || !file.exists(bsp_exe) ) {
  cat(glue::glue("{pmssg} Warning: Unable to locate bsp_exe={bsp_exe}. ",
                 "Will try docker version: ",
                 "2.90 '/repo/bsmap-2.90/bsmap'.{RET2}"))
  bsp_exe <- '/repo/bsmap-2.90/bsmap'
}

if (!file.exists(bsp_exe)) {
  fail_mssg <- glue::glue("bsp_exe='{bsp_exe}' does NOT exist")
  stop(glue::glue("{perrs} {fail_mssg}!{perrs} Exiting...{RET2}"))
}

bsp_cmd <- glue::glue(
  "{bsp_exe} -a {ord_fas} -d {ref_fas} {bsp_opt} -o {out_bsp}{RET}",
  "cat {out_bsp} | {add_cmd} gzip -c -> {out_tsv}{RET}",
  "rm -f {out_bsp}{RET}",
  "touch {end_tsv}{RET}" )

if ( opt$verbose >= 1 )
  cat(glue::glue("{pmssg} Writing BSMAP shell = {bsp_ssh}...{RET2}"))

out_cnt <- safe_write( x = bsp_cmd, type = "line", file = bsp_ssh, 
                       done = TRUE, write_spec = FALSE, append = FALSE, 
                       permissions = "0777", 
                       vb=opt$verbose, tt=pTracker )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Execute BSMAP::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (vb >= vt)
  cat(glue::glue("{pmssg} Running; CMD = '{bsp_cmd}'...{RET}"))

sys_ret <- 0
sys_ret <- base::system(bsp_ssh)

if ( sys_ret != 0 ) {
  fail_mssg <- glue::glue("sys_ret({sys_ret}) != 0!")
  stop(glue::glue("{perrs} {fail_mssg}!{perrs} Exiting...{RET2}"))
} else {
  if ( opt$verbose >= 1 )
    cat(glue::glue("{pmssg} BSMAP Completed Succesfully!{RET2}"))
  
  # sys_ret <- base::system( glue::glue("touch {end_txt}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Load Data::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

outU_tsv <- file.path(opt$out_path, "unmethyl.bsmap.tsv.gz")
bspU_tib <- safe_read( file = outU_tsv, type = "tsv",
                       clean = TRUE, 
                       use_spec = TRUE, spec_col = bsp_cols, 
                       has_head = FALSE, write_spec = TRUE, 
                       # fun_tag = fun_tag, 
                       vb=opt$verbose, tt=pTracker ) %>%
  dplyr::select( -dplyr::any_of( c("Bsp_Qual") ) )

if ( opt$verbose >= 1 )
  cat(glue::glue("{pmssg} Loaded BSP = {outU_tsv}.{RET}"))

outM_tsv <- file.path(opt$out_path, "methyl.bsmap.tsv.gz")
bspM_tib <- safe_read( file = outM_tsv, type = "tsv",
                       clean = TRUE, 
                       use_spec = TRUE, spec_col = bsp_cols, 
                       has_head = FALSE, write_spec = TRUE, 
                       # fun_tag = fun_tag, 
                       vb=opt$verbose, tt=pTracker ) %>%
  dplyr::select( -dplyr::any_of( c("Bsp_Qual") ) )

if ( opt$verbose >= 1 )
  cat(glue::glue("{pmssg} Loaded BSP = {outM_tsv}.{RET}"))


bsp_join_tib <- bspU_tib %>% 
  dplyr::filter(Tag == "UM") %>% 
  dplyr::inner_join( bspM_tib, 
                     by=c("Probe_Id", "Tag"), suffix=c("_U", "_M") ) %>% 
  tidyr::separate( Probe_Id, 
                   into=c("Seq_ID","Probe_Tag","Score_Min","Min_Cpg"), 
                   sep="_", convert = TRUE) %>% 
  tidyr::unite( Probe_ID, Seq_ID,Probe_Tag, sep="_", remove = FALSE )


bsp_top_tib <- bsp_join_tib %>%
  dplyr::arrange(-Score_Min) %>% 
  dplyr::distinct(Seq_ID, .keep_all = TRUE) %>% 
  dplyr::select(Probe_ID, Seq_ID)

pos_top_csv <- file.path(opt$out_path, "top-selected-supplementary-pool-chicago.csv.gz")
pos_top_tib <- tar_bed_tib %>%
  dplyr::inner_join(bsp_top_tib, by="Seq_ID") %>%
  dplyr::distinct(Chromosome,Coordinate,Probe_ID)

safe_write( x = pos_top_tib, file = pos_top_csv, vb = opt$verbose )


org_man_csv <- "/Users/bretbarnes/Documents/data/manifests/methylation/Chicago-Ober-Custom/Chicago-S40.manifest.sesame-base.cpg-sorted.csv.gz"
org_man_tib <- readr::read_csv(org_man_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Original Order::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# tar_csv <- "/Users/bretbarnes/Documents/data/CustomContent/Chicago-Ober-Custom/Round1_AddOn/Custom_Methylation_Array_Additions.csv.gz"
tar_csv <- "/Users/bretbarnes/Documents/data/CustomContent/Chicago-Ober-Custom/Round1_AddOn/2021-08-23_Custom_Methylation_Array_Additions_v2.csv.gz"
tar_tib <- readr::read_csv(tar_csv) %>% 
  purrr::set_names( c("Chromosome", "Coordinate") ) %>%
  dplyr::mutate(
    Chromosome = Chromosome %>%
      stringr::str_remove("^chr") %>%
      paste0("chr",.)) %>% 
  split( f = .$Chromosome )

chr_list <- file_list( 
  path = chr_dirs[[opt$ref_build]], 
  pattern = "\\.fa\\.gz$",
  suffix = c("\\.gz$", "\\.fa$") ) %>% 
  purrr::set_names( stringr::str_remove(names(.), "^chr") %>% paste0("chr", .) )






if (FALSE) {
  for ( chr in names(tar_tib) ) {
    
    chr_fas <- chr_list[[chr]]
    if ( file.exists( chr_fas ) ) {
      cat(glue::glue( "Success found fasta file: '{chr_fas}'{RET}"))
    } else {
      
    }  
    
    cur_tar_tib <- tar_tib[[chr]] %>%
      dplyr::mutate( 
        # Coordinate_m2 = Coordinate - 62, 
        # Coordinate_m1 = Coordinate - 61,
        Coordinate_p0 = Coordinate - 60
        # Coordinate_p1 = Coordinate - 59,
        # Coordinate_p2 = Coordinate - 58
      )
    
    chr_seq = Biostrings::readDNAStringSet(filepath = chr_fas, format = "fasta")
    
    chr_grs <- parse_short_seq( chr_seq = chr_seq[[1]],
                                sub_seq = "CG",
                                flank_len = 60,
                                remove_ns = TRUE,
                                vb = 10, vt = 1 )
    
    # tar_m2_idx <- which( chr_grs@ranges@start %in% cur_tar_tib$Coordinate_m2 )
    # tar_m1_idx <- which( chr_grs@ranges@start %in% cur_tar_tib$Coordinate_m1 )
    tar_p0_idx <- which( chr_grs@ranges@start %in% cur_tar_tib$Coordinate_p0 )
    # tar_p1_idx <- which( chr_grs@ranges@start %in% cur_tar_tib$Coordinate_p1 )
    # tar_p2_idx <- which( chr_grs@ranges@start %in% cur_tar_tib$Coordinate_p2 )
    
    tar_all_len <- tar_tib[[chr]] %>% 
      dplyr::distinct(Coordinate) %>% 
      dplyr::pull(Coordinate) %>% length()
    # tar_m2_len  <- tar_m2_idx %>% length()
    # tar_m1_len  <- tar_m1_idx %>% length()
    tar_p0_len  <- tar_p0_idx %>% length()
    # tar_p1_len  <- tar_p1_idx %>% length()
    # tar_p2_len  <- tar_p2_idx %>% length()
    
    # chr_grs[tar_idx, ]
    
    if ( opt$verbose >= 0 ) {
      
      cat(glue::glue("{pmssg} opt$ref_build({opt$ref_build}) tar_all_len = {tar_all_len}.{RET}"))
      cat(glue::glue("{pmssg} opt$ref_build({opt$ref_build})  tar_m2_len = {tar_m2_len}.{RET}"))
      cat(glue::glue("{pmssg} opt$ref_build({opt$ref_build})  tar_m1_len = {tar_m1_len}.{RET}"))
      cat(glue::glue("{pmssg} opt$ref_build({opt$ref_build})  tar_p0_len = {tar_p0_len}.{RET}"))
      cat(glue::glue("{pmssg} opt$ref_build({opt$ref_build})  tar_p1_len = {tar_p1_len}.{RET}"))
      cat(glue::glue("{pmssg} opt$ref_build({opt$ref_build})  tar_p2_len = {tar_p2_len}.{RET}"))
      
    }
    
    imp_dir <- file.path( opt$out_path, "improbe" )
    imp_tsv <- file.path( imp_dir, paste0(chr,".improbe-input.tsv.gz") )
    imp_tib <- 
      tibble::tibble( 
        Sequence = as.character( chr_grs ), 
        Chromosome = chr, 
        Coordinate = chr_grs@ranges@start + 60, 
        Genome_Build = opt$ref_build, 
        CpG_Island = "FALSE" ) %>% 
      dplyr::mutate( Seq_ID = dplyr::row_number() ) %>% 
      dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island )
    
    safe_write( x = imp_tib, file = imp_tsv, type = "tsv", vb = opt$verbose )
    
    # next
    break
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=pTracker)

# End of file

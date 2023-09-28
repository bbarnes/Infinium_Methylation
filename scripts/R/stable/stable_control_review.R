
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
par$prgm_tag <- 'stable_control_review'
par$verbose  <- 3
local_paths  <- c( 
  "/Users/bbarnes/Documents/tools/imSuite/scripts/R"
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
par <- source_functions( pars = par, rcpp = FALSE, vb = par$verbose )
par <- params_check( pars = par, args = args, 
                     prgm_aux_check = FALSE, vb = par$verbose )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$version <- 0

par$run_name <- "GRBta1"

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

# Two Bit Genomes::
# lapply( ref_seqs, stringr::str_replace,".gz$",".2bit")
ref_seqs <- NULL
ref_seqs <- file_list( path = gen_path, 
                       subs = paste0("Sequence/WholeGenome",dir_add_str),
                       suffix = c("\\.gz", "\\.fa$", "\\.genome$"), 
                       pattern = "\\.genome\\.fa\\.gz$",
                       ret_type = "list", 
                       recursive = TRUE,
                       subs_exists = FALSE,
                       paths_exists = FALSE,
                       files_exists = FALSE,
                       vb=vb,vt=vt,tc=tc, tt=tt )

# Two Bit Genomes::
# lapply( ref_seqs, stringr::str_replace,".gz$",".2bit")
bit_seqs <- NULL
bit_seqs <- file_list( path = gen_path, 
                       subs = paste0("Sequence/WholeGenome",dir_add_str),
                       suffix = c("\\.2bit", "\\.fa$", "\\.genome$"), 
                       pattern = ".2bit$",
                       ret_type = "list", 
                       recursive = TRUE,
                       subs_exists = FALSE,
                       paths_exists = FALSE,
                       files_exists = FALSE,
                       vb=vb,vt=vt,tc=tc, tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Load Design Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$des_path <- file.path( opt$top_path, "Projects/Controls/data/GRBta1" )

des_csv_vec <- c( file.path(opt$des_path, "GRBta1-bsconversion122mers-filtered.tsv"),
                  file.path(opt$des_path, "GRBta1-nonpolymorphic122mers-filtered.tsv"),
                  file.path(opt$des_path, "GRBta1-specificity122mers-filtered.tsv") )

des_file_list <- NULL
des_file_list <- file_list( path = opt$des_path, 
                            prefix = "GRBta1-",
                            suffix = "122mers-filtered.tsv",
                            pattern = "122mers-filtered.tsv",
                            vb=vb,vt=vt+1,tc=tc, tt=tt )

des_data_list <- NULL
des_data_list <- des_file_list %>% 
  lapply( readr::read_tsv, show_col_types = FALSE )

# des_data_list[["GRBta1-bsconversion122mers-filtered.tsv"]] %>% dplyr::mutate()

des_tib <- NULL
des_tib <- des_data_list %>% 
  dplyr::bind_rows( .id = "Design_Type" ) %>% 
  dplyr::mutate( 
    Probe_Type = Assay_Design_Id %>% stringr::str_sub(1,2),
    Probe_Type = dplyr::case_when(
      Probe_Type == "ch" ~ "ch",
      Probe_Type != "ch" ~ "cg",
      TRUE ~ NA_character_ ),
    SRD_Str = Assay_Design_Id %>% stringr::str_remove("^.*_") %>% stringr::str_remove("[0-9]+$"),
    Strand_TB = SRD_Str %>% stringr::str_sub(1,1),
    Strand_CO = SRD_Str %>% stringr::str_sub(2,2),
    Strand_FR = dplyr::case_when(
      Top_Sequence == Sequence ~ "F",
      Top_Sequence != Sequence ~ "R",
      TRUE ~ NA_character_ ),
    Infinium_Design = dplyr::case_when(
      !is.na(AlleleA_Probe_Sequence) &  is.na(AlleleB_Probe_Sequence) ~ 1.0,
      !is.na(AlleleA_Probe_Sequence) & !is.na(AlleleB_Probe_Sequence) ~ 2.0,
      TRUE ~ NA_real_ ) %>% as.integer(),
    Probe_ID = Assay_Design_Id,
    Forward_Sequence = Sequence
  )

des_sum <- NULL
des_sum <- des_tib %>% 
  print_sum( vec = c("Probe_Type","Design_Type","Infinium_Design",
                     "Strand_TB","Strand_CO","Strand_FR"),
             vb=vb,vt=vt+1,tc=tc, tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Design All Probes
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

design_ids_sym <- rlang::sym( "Probe_ID" )
design_seq_sym <- rlang::sym( "Forward_Sequence" )
design_prb_sym <- rlang::sym( "Probe_Type" )

design_frs_sym <- rlang::sym( "Strand_FR" )
design_cos_sym <- rlang::sym( "Strand_CO" )
design_tbs_sym <- rlang::sym( "Strand_TB" )

imp_cpp_tib <- NULL
imp_cpp_tib <- improbe_seqs_cpp( 
  fwd_vec_r  = dplyr::select( des_tib, !!design_seq_sym ) %>% pull( !!design_seq_sym ) %>% as.vector() %>% shear_brac(),
  din_vec_r  = dplyr::select( des_tib, !!design_prb_sym ) %>% pull( !!design_prb_sym ) %>% as.vector(),
  ids_vec_r  = dplyr::select( des_tib, !!design_ids_sym ) %>% pull( !!design_ids_sym ) %>% as.vector(),
  # Skip these inputs to design all probes...
  # frs_vec_r_ = dplyr::select( des_tib, !!design_frs_sym ) %>% pull( !!design_frs_sym ) %>% as.vector(),
  # cos_vec_r_ = dplyr::select( des_tib, !!design_cos_sym ) %>% pull( !!design_cos_sym ) %>% as.vector(),
  return_source = FALSE, uc = TRUE, 
  vb=vb, vt=vt+10 ) %>% dplyr::bind_rows() %>%
  dplyr::rename( Strand_FR = Strand_SR ) %>%
  dplyr::mutate(
    PRB_U = stringr::str_to_upper(PRB_U),
    PRB_M = stringr::str_to_upper(PRB_M),
    PRB_D = stringr::str_to_upper(PRB_D)
  )

des_prbAU_sum <- NULL
des_prbAU_sum <- des_tib %>% 
  dplyr::inner_join( imp_cpp_tib, 
                     by=c("AlleleA_Probe_Sequence"="PRB_U"),
                     suffix=c("_des","_imp"),
                     multiple = "all" ) %>%
  dplyr::distinct( Assay_Design_Id, .keep_all = TRUE ) %>%
  print_sum( vec = c("Design_Type","Infinium_Design","Probe_Type_des","Probe_Type_imp") )

des_prbAM_sum <- NULL
des_prbAM_sum <- des_tib %>% 
  dplyr::inner_join( imp_cpp_tib, 
                     by=c("AlleleA_Probe_Sequence"="PRB_M"),
                     suffix=c("_des","_imp"),
                     multiple = "all" ) %>%
  dplyr::distinct( Assay_Design_Id, .keep_all = TRUE ) %>%
  print_sum( vec = c("Design_Type","Infinium_Design","Probe_Type_des","Probe_Type_imp") )

des_prbAD_sum <- NULL
des_prbAD_sum <- des_tib %>% 
  dplyr::inner_join( imp_cpp_tib, 
                     by=c("AlleleA_Probe_Sequence"="PRB_D"),
                     suffix=c("_des","_imp"),
                     multiple = "all" ) %>%
  dplyr::distinct( Assay_Design_Id, .keep_all = TRUE ) %>%
  print_sum( vec = c("Design_Type","Infinium_Design","Probe_Type_des","Probe_Type_imp") )

#
# Basically no matches::
#
des_prbBU_sum <- NULL
des_prbBU_sum <- des_tib %>% 
  dplyr::inner_join( imp_cpp_tib, 
                     by=c("AlleleB_Probe_Sequence"="PRB_U"),
                     suffix=c("_des","_imp"),
                     multiple = "all" ) %>%
  dplyr::distinct( Assay_Design_Id, .keep_all = TRUE ) %>%
  print_sum( vec = c("Probe_Type_des","Probe_Type_imp") )

des_prbBM_sum <- NULL
des_prbBM_sum <- des_tib %>% 
  dplyr::inner_join( imp_cpp_tib, 
                     by=c("AlleleB_Probe_Sequence"="PRB_M"),
                     suffix=c("_des","_imp"),
                     multiple = "all" ) %>%
  dplyr::distinct( Assay_Design_Id, .keep_all = TRUE ) %>%
  print_sum( vec = c("Probe_Type_des","Probe_Type_imp") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Compare Probe Designs
#
# Note sure if the comparisons below are really needed...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

des_int_tib <- NULL
des_int_tib <- des_tib %>% 
  dplyr::full_join( imp_cpp_tib, 
                    by=c("Probe_ID"),
                    suffix=c("_des","_imp"),
                    multiple = "all" )

des_cmp_tib <- NULL
des_cmp_tib <- des_int_tib %>% 
  dplyr::mutate(
    Probe_Type = Assay_Design_Id %>% stringr::str_sub(1,2),
    Probe_Type = dplyr::case_when(
      Probe_Type == "ch" ~ "ch",
      Probe_Type != "ch" ~ "cg",
      TRUE ~ NA_character_ )
  ) %>%
  dplyr::select( Probe_ID,Probe_Type,
                 AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, 
                 PRB_U,PRB_M,PRB_D )

#
# TBD:: Need to substring the actual 122mer from the genome...
#
des_cmp_tmp <- NULL
des_cmp_tmp <- des_cmp_tib %>%
  dplyr::filter( Probe_Type != "ch" ) %>% 
  dplyr::select(-Probe_ID,-Probe_Type) %>% 
  dplyr::select( AlleleA_Probe_Sequence,PRB_U,
                 AlleleB_Probe_Sequence,PRB_M) %>%
  head(n=32) %>% 
  as.data.frame()

# Only Matches Found::
# des_cmp_tib %>% dplyr::filter( AlleleA_Probe_Sequence == PRB_D ) %>% dplyr::distinct( Probe_ID, .keep_all = TRUE ) 

# des_cmp_tmp %>% dplyr::filter( AlleleA_Probe_Sequence == PRB_U )
# ATAACACATACTACTATCAAATACTATCAACACACTTTCAAATTCTATTG
# GTAATATATATTATTGTTAAATATTATTAATATATTTTTAGGTTTTGTTG

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Manifest to 122mer::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epic_man_tib <- NULL
epic_man_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/EPIC-8v2-0_A1.body.c1-23.csv.gz" )
epic_man_tib <- readr::read_csv( file = epic_man_csv, show_col_types = FALSE ) %>%
  clean_tib()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Write Passing Probe Fasta Files::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$bsp_path <- safe_mkdir( file.path( opt$out_path, "bsmap" ) )
opt$epic_fas  <- file.path( opt$bsp_path, paste(opt$run_name,"prb.fa.gz", sep='.') )

epic_fas_vec <- NULL
epic_fas_vec <- epic_man_tib %>% 
  dplyr::filter(!is.na( AlleleA_ProbeSeq ) ) %>%
  dplyr::distinct( IlmnID, AlleleA_ProbeSeq ) %>%
  dplyr::mutate( Fasta_Str = paste0(">",IlmnID,"\n",AlleleA_ProbeSeq) ) %>% 
  dplyr::arrange( Fasta_Str ) %>%
  dplyr::pull( Fasta_Str )
readr::write_lines( x = epic_fas_vec, file = opt$epic_fas )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Launch BSMAP on Fasta Files::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Add below to imProbeQC_options()
opt$bsmap_dir <- "/Users/bbarnes/Documents/tools/programs/BSMAPz"
opt$bsmap_exe <- "bsmapz"
opt$bsmap_ref <- ref_seqs$GRCh37

# Two Bit Conversion Command:
# ./tools/ucsc/faToTwoBit /Users/bbarnes/Documents/data/imGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta/GRCh37.genome.fa.gz /Users/bbarnes/Documents/data/imGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta/GRCh37.genome.fa.2bit

epic_bsp_tsv <- NULL
epic_bsp_tsv <- run_bsmap( ref_fas = opt$bsmap_ref,
                           can_fas = opt$epic_fas, 
                           bsp_exe = opt$bsmap_exe, 
                           bsp_dir = opt$bsmap_dir, 
                           slim = TRUE, 
                           out_dir = file.path( opt$out_path ),
                           run_tag = opt$run_name,
                           reload  = opt$reload,
                           reload_min = 10,
                           vb=vb,vt=vt,tc=tc, tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Load BSMAP Alignment Results::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epic_bsp_tib <- NULL
epic_bsp_tib <- load_bsmap( file = epic_bsp_tsv,
                            sort = FALSE, 
                            slim = TRUE,
                            
                            add_cnt  = TRUE,
                            parse_id = FALSE,
                            out_dir  = file.path( opt$out_path ),
                            run_tag  = opt$run_name,
                            reload   = opt$reload,
                            reload_min = 10,
                            vb=vb,vt=vt,tc=tc, tt=tt )

epic_bsp_sum <- NULL
epic_bsp_sum <- epic_bsp_tib %>% 
  print_sum( vec = c("Tag"), vb=vb,vt=vt,tc=tc, tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Extract 122mer Template::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epic_bed_tib <- NULL
epic_bed_tib <- epic_bsp_tib %>% 
  dplyr::mutate(
    prb_len = Aln_Seq %>% stringr::str_length(),
    gen_beg = dplyr::case_when(
      Strand_BS == "+-" ~ Beg - 60,
      Strand_BS == "--" ~ Beg - 60,
      Strand_BS == "++" ~ Beg - 62,
      Strand_BS == "-+" ~ Beg - 62,
      TRUE ~ NA_real_ ),
    gen_end = gen_beg + 122,
    gen_len = gen_end - gen_beg + 1,
    
    Strand_CO = Strand_BS %>% stringr::str_sub(1,1),
    Strand_FR = Strand_BS %>% stringr::str_sub(2,2)
    # Probe_ID  = Probe_ID %>% stringr::str_replace_all("-","_")
  ) %>%
  dplyr::filter( gen_beg >= 0 ) %>%
  dplyr::filter( gen_end >= 0 ) %>%
  dplyr::arrange( Chr,gen_beg,gen_end )

# Dumb QC Summary::
# tmp_tib %>% dplyr::group_by( Strand_CO,Strand_FR,gen_len ) %>% dplyr::summarise( Count=n(), .groups = "drop" )

#
# TBD:: Filter for [CG], write BED and run ./tools/ucsc/faToTwoBit 
#
epic_bed_tsv <- file.path( opt$out_path, "mer122.bed.tsv.gz" )
epic_bed_tib %>% 
  # dplyr::select( Chr,gen_beg,gen_end,Probe_ID,Strand_BS ) %>%
  dplyr::select( Chr,gen_beg,gen_end,Probe_ID ) %>%
  readr::write_tsv( file = epic_bed_tsv, col_names = FALSE )

bsp_bed_tsv <- file.path( opt$out_path, "mer122.bed.tsv.gz" )
epic_bed_tib %>% 
  dplyr::mutate( Beg = Beg - 4 - 1, 
                 End = Beg + 54 ) %>% 
  dplyr::select( Chr,Beg,End,Probe_ID ) %>%
  dplyr::arrange( Chr,Beg,End ) %>%
  readr::write_tsv( file = bsp_bed_tsv, col_names = FALSE )


# epic_man_tib %>% dplyr::filter( IlmnID == "cg17145361_TC21" )
# epic_bsp_tib %>% dplyr::filter( Probe_ID == "cg17145361_TC21" )
#
#                                                           tgGTTCTGCCATTGCTGCTGTGTGGAAGTTCACTCCTGCCTTTTCCTTTCCCta
#                                                          TTGGTTCTGCCATTGCTGCTGTGTGGAAGTTCACTCCTGCCTTTTCCTTTCCCT
# TGTCCTGGACACGCTGTTGGCCTGGATCTGAGCCCTGGTGGAGGTCAAAGCCACCT TTGGTTCTGCCATTGCTGCTGTGTGGAAGTTCACTCCTGCCTTTTCCTTTCCCTAGAGCCTCCACC


#
# TBD: probe_to_122mer()
#      - run_bsmap()
#      - ref_seq -> bowtie()
#      - extract 122mer via 2bit
#
# Add IUPAC codes::
#   - Convert VCF to IUPAC
#   - biostrings::replaceLetterAt() per chromosome
#
# > GCF_000001405_39_all_cnt    <- 1112554629
# > GCF_000001405_39_com_cnt    <- 23194287
# > common_all_20180418_all_cnt <- 37303035
# > GCF_000001405_40_sub_cnt    <- 21506976
#

# Try using BioStrings instead because its dynamic...

# ./tools/ucsc/faToTwoBit 
# ./tools/ucsc/twoBitToFa

# Need to load results as a biostring and compare against original 122mer
#
# NOTE: Use this function dbsnp_table_to_IUPAC()

opt$tar_ref <- "GRCh37"
twoBit_out_fas <- file.path( opt$out_path, paste0(opt$tar_ref,".122mer.twoBit_out.fa") )
twoBitToFa_exe <- file.path( opt$top_path, "tools/ucsc/twoBitToFa" )
twoBitToFa_cmd <- paste0( twoBitToFa_exe," ",
                          "-bed=",epic_bed_tsv," ",
                          bit_seqs[[opt$tar_ref]]," ",
                          twoBit_out_fas )

opt$tar_ref <- "GRCh37"
twoBit_out_fas <- file.path( opt$out_path, paste0(opt$tar_ref,".prbmer.twoBit_out.fa") )
twoBitToFa_exe <- file.path( opt$top_path, "tools/ucsc/twoBitToFa" )
twoBitToFa_cmd <- paste0( twoBitToFa_exe," ",
                          "-bed=",bsp_bed_tsv," ",
                          bit_seqs[[opt$tar_ref]]," ",
                          twoBit_out_fas )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file

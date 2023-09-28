
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
par$prgm_tag <- 'stable_scratch_EPICv2_docker_snps'
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
par$version <- 1
par$version <- 2
par$version <- 3
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
#           SesameData Search for update EPIC v2 AnnoS/I files...
#                    Not Really Useful at this point...
#
# Use this link instead: http://zwdzwd.github.io/InfiniumAnnotation
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sesameData::sesameDataCacheAll()

ses_file_tib <- NULL
ses_file_tib <- sesameData::sesameDataList() %>% 
  as.data.frame() %>% tibble::as_tibble()

ses_epic_dat <- NULL
ses_epic_dat <- sesameData::sesameDataGet( title = "EPIC.probeInfo" )
ses_epic_tib <- NULL
ses_epic_tib <- ses_epic_dat$mapped.probes.hg38 %>% 
  as.data.frame() %>%
  tibble::rownames_to_column( var = "Loci_ID" ) %>% 
  tibble::as_tibble() %>%
  dplyr::rename(
    CHR = seqnames,
    MAPINFO = start,
  ) %>%
  dplyr::mutate(
    CHR = CHR %>% stringr::str_remove("^chr")
  ) %>%
  dplyr::select( Loci_ID,CHR,MAPINFO ) %>% 
  dplyr::distinct()

ses_epic_sum <- NULL
ses_epic_sum <- ses_epic_tib %>% 
  dplyr::mutate( Probe_Type = Loci_ID %>% stringr::str_sub(1,2) ) %>% 
  dplyr::group_by(Probe_Type) %>% 
  dplyr::summarise(Count=n(), .groups = "drop")
if ( p2 ) print( ses_epic_sum, n=base::nrow(ses_epic_sum) ) 

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        dbSNP 151 Sesame Files::
#                                 Anno S/I
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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

opt$probe_bed_dir <- safe_mkdir( file.path( opt$out_path, paste0("probe/bed") ) )
opt$probe_vcf_dir <- safe_mkdir( file.path( opt$out_path, paste0("probe/vcf") ) )
opt$dbSNP_vcf_dir <- safe_mkdir( file.path( opt$out_path, paste0("dbSNP/vcf") ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Format dbSNP_151_vcf
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$dbSNP_151_vcf <- file.path( opt$top_path, "data/dbSNP/NCBI/common_all_20180418.vcf.gz" )

# TBD:: Better done inside of tool...

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Direct Manifest Extraction::
#                                 EPICv2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$extract_manifest <- FALSE

if ( opt$extract_manifest ) {
  # Building Useful Manifest Partitions::
  #
  # Tots 937699
  # Head 8
  # Ctls 937064
  #
  # cat data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv | cut -d, -f 1-23 | tail -n 937692 > tmp/EPIC-8v2-0_A1.no-header.csv
  # Tots 937692
  # Ctls 937057
  # Last ch.21.43742285F_BC21
  #
  # cat data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv | cut -d, -f 1-23 | head -n 7 > tmp/EPIC-8v2-0_A1.head-1-23.csv 
  # cat tmp/EPIC-8v2-0_A1.no-header.csv | head -n 937056 > tmp/EPIC-8v2-0_A1.body-1-23.csv 
  # cat tmp/EPIC-8v2-0_A1.no-header.csv | tail -n 636 > tmp/EPIC-8v2-0_A1.ctls-1-23.csv
  #
  # cp tmp/EPIC-8v2-0_A1.head-1-23.csv tmp/EPIC-8v2-0_A1.body-1-23.csv tmp/EPIC-8v2-0_A1.ctls-1-23.csv data/manifests/methylation/MethylationEPIC_v2-2/split/
  #
  # Build Clean RDS::
  man_a1_path <- file.path( opt$top_path, "data/manifests/methylation/MethylationEPIC_v2-2/split" )
  man_a1_body_csv <- file.path( man_a1_path, "EPIC-8v2-0_A1.body-1-23.csv" )
  man_a1_body_rds <- file.path( man_a1_path, "EPIC-8v2-0_A1.body-1-23.rds" )
  man_a1_body_tib <- NULL
  man_a1_body_tib <- readr::read_csv( file = man_a1_body_csv, show_col_types = FALSE ) %>% clean_tib()
  # man_a1_body_tib %>% readr::write_rds( file = man_a1_body_rds, compress = "gz" )
  
  man_a1_ctls_tib <- NULL
  man_a1_ctls_csv <- file.path( opt$top_path, "data/manifests/methylation/bgz/epic_ctls.csv.gz" )
  man_a1_ctls_tib <- readr::read_csv( file = man_a1_ctls_csv, show_col_types = FALSE ) %>% clean_tib()
  
  # 450K
  # Tots: 486436
  # Head: 8
  # Ctls: 485586
  # cat data/manifests/methylation/GenomeStudio/humanmethylation450_15017482_v1-2.csv| cut -d, -f 1-17 | tail -n 486429 | head -n 485578 > data/manifests/methylation/GenomeStudio/humanmethylation450_15017482_v1-2.body-1-17.csv
  man_h4_body_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/humanmethylation450_15017482_v1-2.body-1-17.csv" )
  man_h4_body_rds <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/humanmethylation450_15017482_v1-2.body-1-17.rds" )
  man_h4_body_tib <- NULL
  man_h4_body_tib <- readr::read_csv( file = man_h4_body_csv, show_col_types = FALSE ) %>% clean_tib()
  # man_h4_body_tib %>% readr::write_rds( file = man_h4_body_rds, compress = "gz" )
  
  # EPICv1-B2 Body Manifest::
  # gzip -dc data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz | cut -d, -f 1-15 | tail -n 867532 | head -n 866896 > data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.body-1-15.csv
  man_B2_body_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.body-1-15.csv" )
  man_B2_body_rds <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.body-1-15.rds" )
  man_B2_body_tib <- NULL
  man_B2_body_tib <- readr::read_csv( file = man_B2_body_csv, show_col_types = FALSE ) %>% clean_tib()
  # man_B2_body_tib %>% readr::write_rds( file = man_B2_body_rds, compress = "gz" )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Direct Manifest Loading::
#                               EPICv1v2450
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Basic Full Manifests::
man_v1_base_csv <- file.path( opt$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full/EPIC_v1.gs_to_sesame.csv.gz" )
man_v2_base_csv <- file.path( opt$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full/EPIC_v2.gs_to_sesame.csv.gz" )
man_v3_base_csv <- file.path( opt$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full/EPIC_v1_v2.gs_to_sesame.csv.gz" )

man_v1_base_tib <- NULL
man_v1_base_tib <- readr::read_csv( file = man_v1_base_csv, show_col_types = FALSE ) %>% clean_tib() # %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
man_v2_base_tib <- NULL
man_v2_base_tib <- readr::read_csv( file = man_v2_base_csv, show_col_types = FALSE ) %>% clean_tib() # %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
man_v3_base_tib <- NULL
man_v3_base_tib <- readr::read_csv( file = man_v3_base_csv, show_col_types = FALSE ) %>% clean_tib() # %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )

man_a1_path <- file.path( opt$top_path, "data/manifests/methylation/MethylationEPIC_v2-2/split" )
man_a1_body_tib <- NULL
man_a1_body_rds <- file.path( man_a1_path, "EPIC-8v2-0_A1.body-1-23.rds" )
man_a1_body_tib <- readr::read_rds( file = man_a1_body_rds )

srd_a1_snps_tib <- NULL
srd_a1_snps_tib <- man_a1_body_tib %>% 
  dplyr::filter( Probe_Type=="rs" ) %>% 
  dplyr::select( Name,Strand_FR ) %>%
  dplyr::rename( Strand_rs = Strand_FR ) %>%
  dplyr::distinct()

man_h4_body_tib <- NULL
man_h4_body_rds <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/humanmethylation450_15017482_v1-2.body-1-17.rds" )
man_h4_body_tib <- readr::read_rds( file = man_h4_body_rds )

# [DONE]: LEFT OFF HERE: We need Chromosome_hg38 EPICv1 coordinates...
#.  - Pull them from Sesame... [ man_B2_body_tib ]
# Expectation = [ A tibble: 838,881 × 21 ]
man_B2_body_tib <- NULL
man_B2_body_rds <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.body-1-15.rds" )
man_B2_body_tib <- readr::read_rds( file = man_B2_body_rds ) %>%
  dplyr::mutate( Rep_Num = as.integer(1),
                 Strand_FR = Strand,
                 Strand_TB = "T",
                 Strand_CO = "C",
                 Infinium_Design = dplyr::case_when(
                   Infinium_Design_Type == "I"  ~ 1.0,
                   Infinium_Design_Type == "II" ~ 2.0,
                   TRUE ~ NA_real_ ) %>% as.integer(),
                 Probe_Type = IlmnID %>% stringr::str_sub(1,2) ) %>% 
  dplyr::left_join( srd_a1_snps_tib, by=c("Name") ) %>%
  dplyr::mutate(
    Strand_FR = dplyr::case_when(
      Probe_Type == "rs" & is.na(Strand_FR) ~ Strand_rs,
      TRUE ~ Strand_FR ) ) %>%
  dplyr::select( -CHR, -MAPINFO ) %>% 
  dplyr::inner_join( ses_epic_tib, 
                     by=c("IlmnID"="Loci_ID") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Rank Historic Probes::
#                              EPICv1-450K
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Expectation = [ A tibble: 838,881 × 29 ]
rnk_B2_body_tib <- NULL
rnk_B2_body_tib <- rank_historic_probes( 
  tib = man_B2_body_tib, 
  prb_tibs = list( man_v1_base_tib,man_h4_body_tib ), 
  add_build = FALSE,
  man_col_vec = c( 
    "Probe_Type",
    "AddressA_ID","AddressB_ID",
    "AlleleA_ProbeSeq","AlleleB_ProbeSeq",
    "Strand_FR","Strand_TB","Strand_CO","Next_Base",
    "Infinium_Design_Type","Infinium_Design","Rep_Num",
    "Canonical_Rank","History_Cnt","Loci_Cnt" ), 
  
  out_dir = file.path( opt$out_path, "EPICv1" ),
  run_tag = opt$run_name, 
  reload = opt$reload, 
  reload_min = 0,
  write_out = TRUE,
  vb=vb,vt=vt+1,tc=tc,tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Map Against dbSPN151::
#                                 EPICv1
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Expectation = [ A tibble: 133,183 × 34 ]
snp_B2_body_tib <- NULL
snp_B2_body_tib <- intersect_bed_tabix(
  tib = rnk_B2_body_tib, 
  file = opt$dbSNP_151_vcf,
  
  ref_vcf   = TRUE, 
  com_len   = 56, 
  rm_chr0   = TRUE,
  strip_chr = TRUE, 
  snps_only = TRUE,
  match_src = TRUE,
  parse_caf = TRUE,
  min_maf   = 0.0,
  
  ids_key = "Probe_ID", 
  chr_key = "Chromosome", 
  beg_key = "Nxb_Pos_Up", 
  end_key = "Nxb_Pos_Dn",
  beg_buf = 2, end_buf = 2, beg_off = 0,
  
  out_dir = file.path( opt$out_path, "EPICv1" ),
  run_tag = opt$run_name, 
  reload = opt$reload, 
  reload_min = 0,
  write_out = TRUE,
  vb=vb,vt=vt+1,tc=tc,tt=tt )

snp_B2_body_sum <- NULL
snp_B2_body_sum <- snp_B2_body_tib %>% 
  dplyr::group_by( Probe_Type,SNP_Match_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) print( snp_B2_body_sum, n=base::nrow(snp_B2_body_sum) ) 

# SNP Summary::
snp_B2_body_sum2 <- NULL
snp_B2_body_sum2 <- snp_B2_body_tib %>%
  dplyr::select( Chromosome:Loci_ID, MAF,ALT ) %>%
  dplyr::inner_join( annoI_tib, by=c("Loci_ID","SNP_ID") ) %>%
  dplyr::summarise(
    Min = min(MAF),
    Avg = mean(MAF),
    Med = median(MAF), 
    Max = max(MAF),
    Cnt = n()
  )
if ( p2 ) print( snp_B2_body_sum2, n=base::nrow(snp_B2_body_sum2) ) 

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Rank Historic Probes::
#                              EPICv1-450K
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Expectation = [ A tibble: 937,055 × 21 ]
# Expectation = [ A tibble: 937,055 × 29 ]
rnk_a1_body_tib <- NULL
rnk_a1_body_tib <- rank_historic_probes( 
  tib = man_a1_body_tib, 
  prb_tibs = list( man_v1_base_tib,man_h4_body_tib ), 
  add_build = FALSE,
  man_col_vec = c( 
    "Probe_Type",
    "AddressA_ID","AddressB_ID",
    "AlleleA_ProbeSeq","AlleleB_ProbeSeq",
    "Strand_FR","Strand_TB","Strand_CO","Next_Base",
    "Infinium_Design_Type","Infinium_Design","Rep_Num",
    "Canonical_Rank","History_Cnt","Loci_Cnt" ), 
  
  out_dir = file.path( opt$out_path, "EPICv2" ),
  run_tag = opt$run_name, 
  reload = opt$reload, 
  reload_min = 0,
  write_out = TRUE,
  vb=vb,vt=vt+1,tc=tc,tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Map Against dbSPN151::
#                                 EPICv2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Expectation(+1): [ A tibble: 132,307 × 8 ]
# Expectation(+2): [ A tibble: 148,103 × 8 ]
# Expectation(+2): [ A tibble: 156,299 × 5 ]
# Expectation = [ A tibble: 132,611 × 34 ]
# Expectation = [ A tibble: 132,611 × 34 ]
snp_a1_body_tib <- NULL
snp_a1_body_tib <- intersect_bed_tabix(
  tib = rnk_a1_body_tib, 
  file = opt$dbSNP_151_vcf,
  
  ref_vcf   = TRUE, 
  com_len   = 56, 
  rm_chr0   = TRUE,
  strip_chr = TRUE, 
  snps_only = TRUE,
  match_src = TRUE,
  parse_caf = TRUE,
  min_maf   = 0.0,

  ids_key = "Probe_ID", 
  chr_key = "Chromosome", 
  beg_key = "Nxb_Pos_Up", 
  end_key = "Nxb_Pos_Dn",
  beg_buf = 2, end_buf = 2, beg_off = 0,
  
  out_dir = file.path( opt$out_path, "EPICv2" ),
  run_tag = opt$run_name, 
  reload = opt$reload, 
  reload_min = 0,
  write_out = TRUE,
  vb=vb,vt=vt+1,tc=tc,tt=tt )

snp_a1_body_sum <- NULL
snp_a1_body_sum <- snp_a1_body_tib %>% 
  dplyr::group_by( Probe_Type,SNP_Match_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) print( snp_a1_body_sum, n=base::nrow(snp_a1_body_sum) ) 


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       [TBD]: Review Strand Issue::
#                     [TBD]: Review Duplicate Issue::
#
# CONCLUSION: This has to do with Probe_ID NOT Loci_ID
#
# [TBD]: Validate Probe_ID Unique  Conclusion
# [TBD]: Earlier strand swap into B2 manifest 
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  
  # annoS_a1_clean_tib %>% dplyr::filter( Probe_Type == "rs" )
  
  snp_B2_body_tib %>% 
    dplyr::group_by( Probe_Type,strand,Strand_FR ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  snp_a1_body_tib %>% 
    dplyr::group_by( Probe_Type,strand,Strand_FR ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  
  snp_B2_dupcnt_tib <- snp_B2_body_tib %>% dplyr::add_count( Chromosome,Loci_ID,SNP_ID,REF_SNP,ALT, name = "Dup_Cnt" ) %>% dplyr::filter( Dup_Cnt > 1 )
  snp_a1_dupcnt_tib <- snp_a1_body_tib %>% dplyr::add_count( Chromosome,Loci_ID,SNP_ID,REF_SNP,ALT, name = "Dup_Cnt" ) %>% dplyr::filter( Dup_Cnt > 1 )
  
  snp_a1_dupcnt_tib %>%
    dplyr::filter( Probe_Type == "rs" ) %>%
    dplyr::select( Chromosome,Coordinate_SNP,Loci_ID,SNP_ID,REF_SNP,ALT,Dup_Cnt ) %>%
    dplyr::arrange( Chromosome,Coordinate_SNP )
  
  annoI_B2_dupcnt_tib <- annoI_B2_clean_tib %>% dplyr::add_count( Chromosome,Loci_ID,SNP_ID,REF_SNP,ALT, name = "Dup_Cnt" ) %>% dplyr::filter( Dup_Cnt > 1 )
  annoS_B2_dupcnt_tib <- annoS_B2_clean_tib %>% dplyr::add_count( Chromosome,Loci_ID,SNP_ID,REF_SNP,ALT, name = "Dup_Cnt" ) %>% dplyr::filter( Dup_Cnt > 1 )
  annoI_a1_dupcnt_tib <- annoI_a1_clean_tib %>% dplyr::add_count( Chromosome,Loci_ID,SNP_ID,REF_SNP,ALT, name = "Dup_Cnt" ) %>% dplyr::filter( Dup_Cnt > 1 )
  annoS_a1_dupcnt_tib <- annoS_a1_clean_tib %>% dplyr::add_count( Chromosome,Loci_ID,SNP_ID,REF_SNP,ALT, name = "Dup_Cnt" ) %>% dplyr::filter( Dup_Cnt > 1 )
  
  annoI_a1_dupcnt_tib %>%
    # dplyr::filter( Probe_Type == "rs" ) %>%
    dplyr::select( Chromosome,Coordinate_SNP,Loci_ID,SNP_ID,REF_SNP,ALT,Dup_Cnt ) %>%
    dplyr::arrange( Chromosome,Coordinate_SNP )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Generate/Validate AnnoI GRS::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# [DONE]: write/validate snps_to_annoI()
#

# Expectation Before ALT Screening: [ A tibble: 1,919 × 48 ]
annoI_B2_clean_tib <- NULL
annoI_B2_clean_tib <- snps_to_annoI( tib = snp_B2_body_tib, vb=vb,vt,tc=tc )

annoI_B2_clean_sum <- NULL
annoI_B2_clean_sum <- annoI_B2_clean_tib %>%
  # Screen after historic joining...
  dplyr::inner_join( annoI_tib, by=c("Loci_ID","SNP_ID") ) %>%
  dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT.x ) %>% 
  #
  # Screen before historic joining...
  # dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) print( annoI_B2_clean_sum, n=base::nrow(annoI_B2_clean_sum) )

# Expectation Before ALT Screening: [ A tibble: 1,259 × 48 ]
annoI_a1_clean_tib <- NULL
annoI_a1_clean_tib <- snps_to_annoI( tib = snp_a1_body_tib, vb=vb,vt,tc=tc )

annoI_a1_clean_sum <- NULL
annoI_a1_clean_sum <- annoI_a1_clean_tib %>%
  # Screen after historic joining...
  dplyr::inner_join( annoI_tib, by=c("Loci_ID","SNP_ID") ) %>%
  dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT.x ) %>% 
  #
  # Screen before historic joining...
  # dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) print( annoI_a1_clean_sum, n=base::nrow(annoI_a1_clean_sum) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Generate/Validate AnnoS GRS::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# [DONE]: write/validate snps_to_annoI()
#

# Expectation Before ALT Screening: [ A tibble: 116 × 48 ]
annoS_B2_clean_tib <- NULL
annoS_B2_clean_tib <- snps_to_annoS( tib = snp_B2_body_tib, vb=vb,vt=vt,tc=tc )

annoS_B2_clean_sum <- NULL
annoS_B2_clean_sum <- annoS_B2_clean_tib %>%
  # Screen after historic joining...
  dplyr::inner_join( annoS_tib, by=c("Loci_ID","SNP_ID") ) %>%
  dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT.x ) %>% 
  #
  # Screen before historic joining...
  # dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) print( annoS_B2_clean_sum, n=base::nrow(annoS_B2_clean_sum) )


# Expectation Before ALT Screening: [ A tibble: 147 × 48 ]
annoS_a1_clean_tib <- NULL
annoS_a1_clean_tib <- snps_to_annoS( tib = snp_a1_body_tib, vb=vb,vt=vt,tc=tc )

annoS_a1_clean_sum <- NULL
annoS_a1_clean_sum <- annoS_a1_clean_tib %>%
  # Screen after historic joining...
  dplyr::inner_join( annoS_tib, by=c("Loci_ID","SNP_ID") ) %>%
  dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT.x ) %>% 
  #
  # Screen before historic joining...
  # dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) print( annoS_a1_clean_sum, n=base::nrow(annoS_a1_clean_sum) )


# [DONE]: Intersect all at once...
# [DONE]: Update intersect_bed_tabix()
#         - intersect_bed_tabix( tib = bed_tib, file = opt$dbSNP_151_vcf, ids_key = )
# [DONE]: Investigate MAF screening and removing low MAF
#         - Done before...
# [DONE]: Seems that we're off on our coordinates...
# [DONE]: Remove all the extra fields Up/Dn stuff
# [DONE]: Re-run everything but with EPICv1
# [DONE]: Add "EPICv1" to run names...
#
# [DONE]: Implement screening inside intersect_bed_tabix() with 'Match Types'
# [DONE]: Implement matching by REF
#
# [TBD]: Convert anno's to GRS/write output
# [TBD]: Add copy commands to Docker Image { Manifestv2, annoI, annoS }
#

annoI_B2_clean_grs <- NULL
annoI_B2_clean_grs <- annoI_to_grs( tib = annoI_B2_clean_tib, 
                                    add_chr = TRUE,
                                    vb=vb,vt=vt,tc=tc, tt=tt )
annoI_a1_clean_grs <- NULL
annoI_a1_clean_grs <- annoI_to_grs( tib = annoI_a1_clean_tib,
                                    add_chr = TRUE,
                                    vb=vb,vt=vt,tc=tc, tt=tt )

annoS_B2_clean_grs <- NULL
annoS_B2_clean_grs <- annoI_to_grs( tib = annoS_B2_clean_tib, 
                                    add_chr = TRUE,
                                    vb=vb,vt=vt,tc=tc, tt=tt )
annoS_a1_clean_grs <- NULL
annoS_a1_clean_grs <- annoI_to_grs( tib = annoS_a1_clean_tib,
                                    add_chr = TRUE,
                                    vb=vb,vt=vt,tc=tc, tt=tt )


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   Investigate/Summarize EPICv1 Results::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  
  #
  # Rebuild annoI and annoS:: EPICv1 (B2)
  #
  anno_B2_tib <- NULL
  anno_B2_tib <- snp_B2_body_tib %>% 
    dplyr::mutate( Chromosome = paste0("chr",Chromosome),
                   In.band = "REF" ) %>%
    dplyr::select( Chromosome,Coordinate_SNP,strand,Loci_ID,SNP_ID,
                   Infinium_Design_Type,In.band,REF_SNP,ALT_SNP,
                   SNP_Match_Type,Probe_Type )
  
  annoI_B2_mat_tib <- NULL
  annoI_B2_mat_tib <- anno_B2_tib %>% 
    dplyr::inner_join( annoI_tib, by=c("Loci_ID","SNP_ID") )
  
  annoS_B2_mat_tib <- NULL
  annoS_B2_mat_tib <- anno_B2_tib %>% 
    dplyr::inner_join( annoS_tib, by=c("Loci_ID","SNP_ID") )
  
  # Type Summary::
  #
  anno_B2_sum <- NULL
  anno_B2_sum <- anno_B2_tib %>% 
    dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) print( anno_B2_sum, n=base::nrow(anno_B2_sum) ) 
  
  # Full Summary::
  #
  annoI_B2_mat_sum1 <- NULL
  annoI_B2_mat_sum1 <- annoI_B2_mat_tib %>%
    # dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT_SNP ) %>% 
    dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) print( annoI_B2_mat_sum1, n=base::nrow(annoI_B2_mat_sum1) ) 
  
  annoS_B2_mat_sum1 <- NULL
  annoS_B2_mat_sum1 <- annoS_B2_mat_tib %>%
    dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, REF_SNP,ALT_SNP ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) print( annoS_B2_mat_sum1, n=base::nrow(annoS_B2_mat_sum1) ) 
  
  # Nucleotide Summary::
  #
  annoI_B2_mat_sum2 <- NULL
  annoI_B2_mat_sum2 <- annoI_B2_mat_tib %>%
    dplyr::group_by( REF_SNP,ALT_SNP ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) print( annoI_B2_mat_sum2, n=base::nrow(annoI_B2_mat_sum2) ) 
  
  annoS_B2_mat_sum2 <- NULL
  annoS_B2_mat_sum2 <- annoS_B2_mat_tib %>%
    dplyr::group_by( REF_SNP,ALT_SNP ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) print( annoS_B2_mat_sum2, n=base::nrow(annoS_B2_mat_sum2) ) 
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Investigation Code Below::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  
  tmp_tib <- NULL
  tmp_tib <- snp_a1_body_tib %>% dplyr::mutate( 
    REF_Match = dplyr::case_when( 
      REF_SNP == Nxb_Nuc_Up ~ "0",
      REF_SNP == Cpg_Nuc_Up ~ "2",
      REF_SNP == Cpg_Nuc_Dn ~ "3",
      REF_SNP == Nxb_Nuc_Dn ~ "1",
      TRUE ~ NA_character_ )
  )
  
  tmp_sum <- NULL
  tmp_sum <- tmp_tib %>% 
    dplyr::group_by( Probe_Type,SNP_Match_Type,REF_Match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) print( tmp_sum, n=base::nrow(tmp_sum) )
  
  annoI_join_tib <- NULL
  annoI_join_tib <- annoI_tib %>% 
    dplyr::inner_join( tmp_tib, by=c("Loci_ID"), multiple = "all" )
  
  annoI_join_sum <- NULL
  annoI_join_sum <- annoI_join_tib %>%
    dplyr::arrange( REF_Match ) %>%
    dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>%
    dplyr::group_by( Probe_Type,SNP_Match_Type,REF_Match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) print( annoI_join_sum, n=base::nrow(annoI_join_sum) )
  
  
  
  
  
  
  
  
  snp_a1_body_tib %>% head(n=3) %>% as.data.frame()
  snp_a1_body_tib %>% dplyr::select( Chromosome,Coordinate, Coordinate_SNP:Next_Base, Infinium_Design,Rep_Num,strand ) %>% head(n=3) %>% as.data.frame()
  
  tmp_tib <- NULL
  tmp_tib <- snp_a1_body_tib %>% dplyr::mutate( 
    REF_Match = dplyr::case_when( 
      REF_SNP == Nxb_Nuc_Up ~ "NU",
      REF_SNP == Cpg_Nuc_Up ~ "CU",
      REF_SNP == Cpg_Nuc_Dn ~ "CD",
      REF_SNP == Nxb_Nuc_Dn ~ "ND",
      TRUE ~ NA_character_ )
  )
  
  tmp_tib %>% head(n=3) %>% as.data.frame()
  
  tmp_tib %>% 
    dplyr::filter( Loci_ID %in% annoI_tib$Loci_ID ) %>%
    dplyr::group_by( Probe_Type,strand,SNP_Match_Type,REF_Match ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  
  
  tmp_tib %>% 
    dplyr::filter( Loci_ID %in% annoI_tib$Loci_ID ) %>%
    dplyr::filter( SNP_Match_Type == "Nxb_Dn" ) %>%
    dplyr::select( REF_SNP,Nxb_Nuc_Dn, 
                   Probe_Type,strand,SNP_Match_Type,REF_Match )
  
  tmp_tib %>% 
    dplyr::filter( Loci_ID %in% annoI_tib$Loci_ID ) %>%
    dplyr::filter( SNP_Match_Type == "Nxb_Dn" ) %>%
    dplyr::group_by( REF_SNP,Nxb_Nuc_Dn, Probe_Type,strand,SNP_Match_Type,REF_Match ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  
  tmp_tib %>% 
    dplyr::filter( Loci_ID %in% annoI_tib$Loci_ID ) %>%
    # dplyr::filter( SNP_Match_Type == "Nxb_Dn" ) %>%
    dplyr::group_by( REF_SNP,Nxb_Nuc_Dn, Probe_Type,strand,SNP_Match_Type,REF_Match ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% print( n=1000 )
  
  
  
  #
  # [Done]: Fix the naming below...
  #
  nxb_up_tib <- NULL
  nxb_up_tib <- snp_a1_body_tib %>%
    dplyr::inner_join( 
      dplyr::mutate( rnk_a1_body_tib, Chromosome = Chromosome %>% stringr::str_remove("^chr") ),
      by=c("Chromosome","Coordinate"="Nxb_Pos_Up"),
      multiple = "all"
    )
  cpg_up_tib <- NULL
  cpg_up_tib <- snp_a1_body_tib %>%
    dplyr::inner_join( 
      dplyr::mutate( rnk_a1_body_tib, Chromosome = Chromosome %>% stringr::str_remove("^chr") ),
      by=c("Chromosome","Coordinate"="Cpg_Pos_Up"),
      multiple = "all"
    )
  cpg_dn_tib <- NULL
  cpg_dn_tib <- snp_a1_body_tib %>%
    dplyr::inner_join( 
      dplyr::mutate( rnk_a1_body_tib, Chromosome = Chromosome %>% stringr::str_remove("^chr") ),
      by=c("Chromosome","Coordinate"="Cpg_Pos_Dn"),
      multiple = "all"
    )
  nxb_dn_tib <- NULL
  nxb_dn_tib <- snp_a1_body_tib %>%
    dplyr::inner_join( 
      dplyr::mutate( rnk_a1_body_tib, Chromosome = Chromosome %>% stringr::str_remove("^chr") ),
      by=c("Chromosome","Coordinate"="Nxb_Pos_Dn"),
      multiple = "all"
    )
  
  nxb_up_tib %>% base::nrow()
  cpg_up_tib %>% base::nrow()
  cpg_dn_tib %>% base::nrow()
  nxb_dn_tib %>% base::nrow()
  
  #
  # [TBD]:: Evaluate based on known anno matches
  #
  nxb_up_matI_tib <- nxb_up_tib %>% dplyr::inner_join( annoI_tib, by=c("Loci_ID"), multiple = "all" )
  cpg_up_matI_tib <- cpg_up_tib %>% dplyr::inner_join( annoI_tib, by=c("Loci_ID"), multiple = "all" )
  cpg_dn_matI_tib <- cpg_dn_tib %>% dplyr::inner_join( annoI_tib, by=c("Loci_ID"), multiple = "all" )
  nxb_dn_matI_tib <- nxb_dn_tib %>% dplyr::inner_join( annoI_tib, by=c("Loci_ID"), multiple = "all" )
  
  nxb_up_matI_tib %>% base::nrow()
  cpg_up_matI_tib %>% base::nrow()
  cpg_dn_matI_tib %>% base::nrow()
  nxb_dn_matI_tib %>% base::nrow()
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Addition One Off Copy Commands::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Also see end of: "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/notes/docker.EPICv2-VA.history-notes.29962023.txt"
#.  for more details...
#
# epicv1_annoS_cp_cmd <- "cp Projects.new/EPIC_v2/docker/dat/annoS.rds Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoS.rds"
# epicv1_annoI_cp_cmd <- "cp Projects.new/EPIC_v2/docker/dat/annoI.rds Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoI.rds"
#
# [Wrong]: docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/dat/EPICv2/EPICv2.annoS.rds dd00d8b92b93:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPICv2/
# [Wrong]: docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/dat/EPICv2/EPICv2.annoI.rds dd00d8b92b93:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPICv2/
#
# docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/dat/EPICv2/EPICv2.annoS.rds dd00d8b92b93:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-A1/EPIC-A1.annoS.rds
# docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/dat/EPICv2/EPICv2.annoI.rds dd00d8b92b93:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-A1/EPIC-A1.annoI.rds

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file

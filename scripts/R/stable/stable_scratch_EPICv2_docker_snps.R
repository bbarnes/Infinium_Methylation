
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
par$version <- 4
par$version <- 5

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

opt$probe_man_dir <- safe_mkdir( file.path( opt$out_path, paste0("probe/man") ) )
opt$probe_bed_dir <- safe_mkdir( file.path( opt$out_path, paste0("probe/bed") ) )
opt$probe_vcf_dir <- safe_mkdir( file.path( opt$out_path, paste0("probe/vcf") ) )
opt$dbSNP_vcf_dir <- safe_mkdir( file.path( opt$out_path, paste0("dbSNP/vcf") ) )

opt$docV1_man_dir <- safe_mkdir( file.path( opt$out_path, paste0("docker/man/EPICv1") ) )
opt$docV2_man_dir <- safe_mkdir( file.path( opt$out_path, paste0("docker/man/EPICv2") ) )

opt$docV1_grs_dir <- safe_mkdir( file.path( opt$out_path, paste0("docker/dat/EPICv1") ) )
opt$docV2_grs_dir <- safe_mkdir( file.path( opt$out_path, paste0("docker/dat/EPICv2") ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Format dbSNP_151_vcf
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$dbSNP_151_vcf <- file.path( opt$top_path, "data/dbSNP/NCBI/common_all_20180418.vcf.gz" )

# [TBD]:: Better done inside of tool...

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
ctl_v2_base_csv <- file.path( opt$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full/EPIC-A1.manifest.sesame-base.cpg-sorted.ctl.csv.gz" )

man_v1_base_tib <- NULL
man_v1_base_tib <- readr::read_csv( file = man_v1_base_csv, show_col_types = FALSE ) %>% clean_tib() # %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
man_v2_base_tib <- NULL
man_v2_base_tib <- readr::read_csv( file = man_v2_base_csv, show_col_types = FALSE ) %>% clean_tib() # %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
man_v3_base_tib <- NULL
man_v3_base_tib <- readr::read_csv( file = man_v3_base_csv, show_col_types = FALSE ) %>% clean_tib() # %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )

ctl_v2_base_tib <- NULL
ctl_v2_base_tib <- readr::read_csv( file = ctl_v2_base_csv, show_col_types = FALSE ) %>% clean_tib() %>%
  dplyr::mutate( Full_ID = Probe_ID ) %>%
  dplyr::select( -mask )

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
#.                         Manifest CpG Compare::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Just verify the CGN's didn't change, so check by probe sequence...
# Conclusion: Doesn't matter...
if ( FALSE ) {
  man_h4_body_tib %>% 
    dplyr::filter( Name %in% annoI_tib$Loci_ID ) %>%
    dplyr::distinct( Name ) %>% base::nrow()
  
  man_a1_body_tib %>% 
    dplyr::filter( Name %in% annoI_tib$Loci_ID ) %>%
    dplyr::distinct( Name ) %>% base::nrow()
  
  man_h4_body_tib %>% 
    dplyr::filter( Name %in% annoI_tib$Loci_ID ) %>%
    dplyr::filter( AlleleA_ProbeSeq %in% man_a1_body_tib$AlleleA_ProbeSeq ) %>%
    dplyr::distinct( Name ) %>% base::nrow()
}

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
    "Infinium_Design_Type","Infinium_Design",
    "Color_Channel","col","Rep_Num",
    "Canonical_Rank","History_Cnt","Loci_Cnt" ), 
  
  out_dir = file.path( opt$out_path, "EPICv2" ),
  run_tag = opt$run_name, 
  reload = opt$reload, 
  reload_min = 0,
  write_out = TRUE,
  vb=vb,vt=vt+1,tc=tc,tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Write EPICv2 A1 Manifest::
#                              Docker Version
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doc_a1_full_csv <- file.path( opt$docV2_man_dir, "EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
doc_a1_full_tib <- NULL
doc_a1_full_tib <- rnk_a1_body_tib %>%
  dplyr::mutate( Probe_Source = "EPIC-A1",
                 # Docker Sesame uses col=NA for Inf2
                 # col = dplyr::case_when(
                 #   is.na(col) ~ "2",
                 #   TRUE ~ col ),
                 # Docker Sesame doesn't need 'mask' column
                 # mask = FALSE,
                 Color_Channel = dplyr::case_when(
                   is.na(Color_Channel) ~ "Both",
                   TRUE ~ Color_Channel ) ) %>%
  dplyr::rename( U = AddressA_ID,
                 M = AddressB_ID,
                 DESIGN = Infinium_Design_Type,
                 Probe_Design = Infinium_Design,
                 COLOR_CHANNEL = Color_Channel
  ) %>%
  dplyr::select( Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,
                 Probe_Type,Probe_Source,Next_Base,Probe_Design,Full_ID ) %>%
  clean_tib() %>%
  dplyr::bind_rows( ctl_v2_base_tib )

readr::write_csv( x = doc_a1_full_tib, file = doc_a1_full_csv )

doc_a1_full_sum <- NULL
doc_a1_full_sum <- doc_a1_full_tib %>% 
  dplyr::group_by( Probe_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) print( doc_a1_full_sum, n=base::nrow(doc_a1_full_sum) ) 

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Map Against dbSPN151::
#                                 EPICv2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Expectation = [ A tibble: 132,615 × 48 ]
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

annoI_B2_clean_rds <- file.path( opt$docV1_grs_dir, "EPICv1.annoI.rds" )
annoI_B2_clean_grs <- NULL
annoI_B2_clean_grs <- anno_to_grs( tib = annoI_B2_clean_tib, 
                                   add_chr = TRUE, isAnnoS = FALSE,
                                   vb=vb,vt=vt,tc=tc, tt=tt )
readr::write_rds( x = annoI_B2_clean_grs, file = annoI_B2_clean_rds, compress = "gz" )

annoI_a1_clean_rds <- file.path( opt$docV2_grs_dir, "EPICv2.annoI.rds" )
annoI_a1_clean_grs <- NULL
annoI_a1_clean_grs <- anno_to_grs( tib = annoI_a1_clean_tib %>% dplyr::filter( SNP_ID!="rs782688735"),
                                   add_chr = TRUE, isAnnoS = FALSE,
                                   vb=vb,vt=vt,tc=tc, tt=tt )
readr::write_rds( x = annoI_a1_clean_grs, file = annoI_a1_clean_rds, compress = "gz" )


annoS_B2_clean_rds <- file.path( opt$docV1_grs_dir, "EPICv1.annoS.rds" )
annoS_B2_clean_grs <- NULL
annoS_B2_clean_grs <- anno_to_grs( tib = annoS_B2_clean_tib, 
                                   add_chr = TRUE, isAnnoS = TRUE,
                                   vb=vb,vt=vt,tc=tc, tt=tt )
readr::write_rds( x = annoS_B2_clean_grs, file = annoS_B2_clean_rds, compress = "gz" )

annoS_a1_clean_rds <- file.path( opt$docV2_grs_dir, "EPICv2.annoS.rds" )
annoS_a1_clean_grs <- NULL
annoS_a1_clean_grs <- anno_to_grs( tib = annoS_a1_clean_tib,
                                   add_chr = TRUE, isAnnoS = TRUE,
                                   vb=vb,vt=vt,tc=tc, tt=tt )
readr::write_rds( x = annoS_a1_clean_grs, file = annoS_a1_clean_rds, compress = "gz" )

# [TBD]: Multiple ALTS Still seems incorrect...
#   - Pretty sure you need to take into account the actual probe direction as
#.    well as the Nxb_Up/Nxb_Dn

#
# PREVIOUS ERROR CHECKING BELOW::
#
if ( FALSE ) {
  
  # [TBD]: NOT UNIQUE: annoI_a1_clean_grs by 1
  # [TBD]: NOT UNIQUE: annoS_a1_clean_grs by 147 - 81
  
  # 937690 = doc_a1_full_tib %>% base::nrow()
  # 937690 = doc_a1_full_tib %>% dplyr::distinct( Probe_ID ) %>% base::nrow()
  
  # 1258 = doc_a1_full_tib %>% dplyr::filter( Probe_ID %in% names(annoI_a1_clean_grs) ) %>% base::nrow()
  #   81 = doc_a1_full_tib %>% dplyr::filter( Probe_ID %in% names(annoS_a1_clean_grs) ) %>% base::nrow()
  # 1259 = annoI_a1_clean_tib %>% base::nrow()
  #  147 = annoS_a1_clean_tib %>% base::nrow()
  
  annoI_a1_clean_grs %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble() %>% dplyr::add_count( Probe_ID, name = "dupcnt") %>% dplyr::filter( dupcnt > 1 )
  
  
  annoI_mis_tib <- NULL
  annoI_mis_tib <- annoI_a1_clean_tib %>% dplyr::add_count( Probe_ID, name = "dupcnt" ) %>% dplyr::arrange(Probe_ID) %>% dplyr::filter( dupcnt > 1 )
  
  annoI_a1_clean_tib %>% dplyr::add_count( Probe_ID, name = "dupcnt" ) %>% dplyr::arrange(Probe_ID) %>% dplyr::filter( dupcnt > 1 ) %>% as.data.frame()
  # How many reps of this: cg24657313
  rnk_a1_body_tib %>% dplyr::filter( Loci_ID == "cg24657313" )
  # Solution Manifest Duplicate: cg24657313 to [ cg24657313, cg24657313_BC11 ]
  # Solution Remove SNP_ID == rs782688735
  
  annoS_mis_tib <- NULL
  annoS_mis_tib <- annoS_a1_clean_tib %>% dplyr::add_count( Probe_ID, name = "dupcnt" ) %>% dplyr::arrange(Probe_ID) %>% dplyr::filter( dupcnt > 1 )
  annoS_mis_sum <- NULL
  annoS_mis_sum <- annoS_mis_tib %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  # A tibble: 2 × 2
  # Probe_Type Count
  # <chr>      <int>
  # 1 nv             4
  # 2 rs           128
  annoS_mis_tib %>% dplyr::filter( Probe_Type == "nv" ) %>% as.data.frame()
  annoS_a1_clean_grs[ c("nv-GRCh38-chr17-7676040-7676040-C-G", "nv-GRCh38-chr17-7676040-7676040-C-T") ]
  
  annoS_a1_clean_grs[ annoS_mis_tib$Probe_ID ]
  
  annoS_a1_clean_grs[ annoS_mis_tib$Probe_ID ] 
  
  annoS_a1_clean_grs2 <- annoS_to_grs( tib = annoS_a1_clean_tib,
                                       add_chr = TRUE,
                                       vb=vb,vt=vt,tc=tc, tt=tt )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Addition One Off Copy Commands::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# LEFT OFF HERE:
# [TBD]: Modify A1 Manifest:
#.       - [Done]: Remove 'mask' field
#        - [Done]: Leave col=NA for Inf2
#        - [Postponed]: Update from A1 to E2
#.       - Update transfer directories...

# Also see end of: "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/notes/docker.EPICv2-VA.history-notes.29962023.txt"
#.  for more details...
#
# [OLD]: docker cp /Users/bbarnes/Documents/scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v3/probe/man/EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz 4eebd4d212fc:/repo/Infinium_Methylation_Workhorse/dat/manifest/core/
# docker cp /Users/bbarnes/Documents/scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v3/docker/man/EPICv2/EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz 4eebd4d212fc:/repo/Infinium_Methylation_Workhorse/dat/manifest/core/
#
# docker cp /Users/bbarnes/Documents/scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v3/docker/dat/EPICv2/EPICv2.annoS.rds 4eebd4d212fc:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-A1/EPIC-A1.annoS.rds
# docker cp /Users/bbarnes/Documents/scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v3/docker/dat/EPICv2/EPICv2.annoI.rds 4eebd4d212fc:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-A1/EPIC-A1.annoI.rds
#
# docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoS.rds 4eebd4d212fc:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-B4/EPIC-B4.annoS.rds
# docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoI.rds 4eebd4d212fc:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-B4/EPIC-B4.annoI.rds
#

# DOCKER TBD::
#
# [TBD]: Copy most recent code to Projects.new/EPIC_v2/docker/src/Infinium_Methylation_Workhorse/scripts.vXXXX
# [TBD]: Create new docker branch from v4
# [TBD]: Copy Data
# [TBD]: Copy Scripts
#

#
# OLD STUFF::
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

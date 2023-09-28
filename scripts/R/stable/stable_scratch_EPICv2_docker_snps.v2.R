
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
  tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()

annoI_tib <- NULL
annoI_tib <- annoI_dat %>% as.data.frame() %>% 
  tibble::rownames_to_column( var = "Loci_ID" ) %>% 
  tibble::as_tibble() %>%
  # magrittr::set_names( c("Loci_ID") )
  dplyr::rename(
    Chromosome_CpG_hg19 = seqnames,
    Beg_CpG_hg19 = start,
    End_CpG_hg19 = end,
    RS_ID = rs
  )
  

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Real Calls (Beta/Pval) Data From Docker::
#                                EPIC v1/v2
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$snp_dir <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/snps/docker-1.11.15.1.p.0.4.4" )
snps1_vcf <- file.path( opt$snp_dir, "EPICv1/chip-206203800149/swifthoof_main/206203800149_R01C01_EPIC_B4_ind.snps.vcf" )
beta1_csv <- file.path( opt$snp_dir, "EPICv1/chip-206203800149/swifthoof_main/206203800149_R01C01_EPIC_B4_ind.call.dat.csv.gz" )
beta2_csv <- file.path( opt$snp_dir, "EPICv2/chip-206891110001/swifthoof_main/206891110001_R01C01_EPIC_A1_ind.call.dat.csv.gz" )

beta1_tib <- NULL
beta1_tib <- readr::read_csv( file = beta1_csv, show_col_types = FALSE )
beta2_tib <- NULL
beta2_tib <- readr::read_csv( file = beta2_csv, show_col_types = FALSE ) %>% 
  dplyr::mutate( Loci_ID = Probe_ID %>% stringr::str_remove("_.*$") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Discrepency Analysis::
#                                EPIC v1/v2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# AnnoS Diff Analysis
#
betaS1_tib <- NULL
betaS1_tib <- beta1_tib %>% dplyr::filter( Probe_ID %in% names(annoS_dat) )
betaS2_tib <- NULL
betaS2_tib <- beta2_tib %>% dplyr::filter( Probe_ID %in% names(annoS_dat) )

#
# AnnoI Diff Analysis
#
betaI1_tib <- NULL
betaI1_tib <- beta1_tib %>% dplyr::filter( Probe_ID %in% names(annoI_dat) )
betaI2_tib <- NULL
betaI2_tib <- beta2_tib %>% dplyr::filter( Probe_ID %in% names(annoI_dat) )

#
# AnnoI EPICv2 Stripped Names Analysis::
#
lociS2_tib <- NULL
lociS2_tib <- beta2_tib %>% dplyr::inner_join( annoS_tib, by=c("Loci_ID"="Probe_ID") )
missS2_tib <- NULL
missS2_tib <- lociS2_tib %>% dplyr::filter( Probe_ID != Loci_ID )

# Missing Probes::
#
# lociS2_tib <- beta2_tib %>% dplyr::filter( Loci_ID %in% names(annoS_dat) )
# lociS2_tib %>% dplyr::distinct( Loci_ID, .keep_all = TRUE ) %>% dplyr::filter( Loci_ID %in% betaS1_tib$Probe_ID )
# betaS1_tib %>% dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>% dplyr::filter( !Probe_ID %in% lociS2_tib$Loci_ID )

#
# AnnoS EPICv2 Stripped Names Analysis::
#

lociI2_tib <- NULL
lociI2_tib <- beta2_tib %>% dplyr::filter( Loci_ID %in% names(annoI_dat) )
# joinI2_tib <- NULL
# joinI2_tib <- beta2_tib %>% dplyr::inner_join( annoI_tib, by=c("Probe_ID"="Loci_ID") )
# joinI2_tib <- beta2_tib %>% dplyr::inner_join( annoI_tib, by=c("Loci_ID") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Direct Manifest Investigation::
#                                 EPICv2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# cat data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv | head | grep IlmnID | perl -pe 's/,/\n/gi;' | perl -pe 'print $ii++."\t";'
# cat data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv | grep "^cg" | cut -d, -f 48 | perl -pe 's/;/\n/gi' | sort -n | uniq -c 
#. 10331 
#. 749942 0
#. 966702 1
#
# cat data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv | grep -n "Controls"
#. 937064:[Controls],,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
# cat data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv | wc -l             
#. 937699
#
# Last before controls: ch.21.43742285F_BC21
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Direct Manifest Extraction::
#                                 EPICv2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epicv2_man_csv <- file.path( opt$top_path, "data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv" )
epicv2_man_tib <- NULL
epicv2_man_tib <- readr::read_csv( file = epicv2_man_csv, 
                                   skip = 7, n_max = 937064 - 9, 
                                   show_col_types = FALSE )

epicv2_sub_tib <- epicv2_man_tib %>% 
  dplyr::filter( Infinium_Design == 1 ) %>%
  dplyr::filter( Probe_Type == "cg" ) %>%
  dplyr::select( IlmnID,Name, 
                 Strand_FR,Strand_TB,Strand_CO,Infinium_Design,Rep_Num,
                 CHR,MAPINFO,
                 SNP_ID,SNP_DISTANCE,SNP_MinorAlleleFrequency )

# Before Range Reduction::
# > epicv2_bed_RG_tib
# A tibble: 126,710 × 8
epicv2_bed_RG_tib <- NULL
epicv2_bed_RG_tsv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/bed/EPIC-8v2-0_A1.next-base-RG.bed" )
epicv2_bed_RG_tib <- epicv2_sub_tib %>% 
  dplyr::filter( CHR != "0" ) %>%
  dplyr::filter( MAPINFO != 0 ) %>%
  dplyr::mutate( 
    CHR = CHR %>% stringr::str_remove("^chr"), 
    Beg = dplyr::case_when( 
      Strand_FR == "F" ~ MAPINFO + 1, 
      Strand_FR == "R" ~ MAPINFO - 2, 
      TRUE ~ NA_real_ ),
    End = Beg + 1,
    
    # Previous Range Method::
    Beg_Rng_hg38 = MAPINFO - 2,
    End_Rng_hg38 = MAPINFO + 3
  ) %>% 
  dplyr::select( CHR,Beg,End,IlmnID, Strand_FR,Strand_TB,Strand_CO,
                 Beg_Rng_hg38,End_Rng_hg38 ) %>%
  dplyr::mutate( Loci_ID = IlmnID %>% stringr::str_remove("_.*$") ) %>%
  dplyr::arrange( CHR, Beg ) %>%
  dplyr::rename( 
    Chromosome_BED_hg38 = CHR,
    Beg_BED_hg38 = Beg,
    End_BED_hg38 = End,
    Probe_ID = IlmnID
  )
readr::write_tsv( x = epicv2_bed_RG_tib, file = epicv2_bed_RG_tsv, col_names = FALSE )

# Tabix Command::
# tabix -h -R Projects.new/EPIC_v2/docker/bed/EPIC-8v2-0_A1.next-base.bed data/dbSNP/NCBI/common_all_20180418.vcf.gz > Projects.new/EPIC_v2/docker/intersect/EPICv2.dbSNP151-common.next.base.vcf 
#
dbSNP_151_bgzip <- file.path( opt$top_path, "data/dbSNP/NCBI/common_all_20180418.vcf.gz" )
epicv2_int_dbSNP_RG_vcf <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/intersect/EPICv2.dbSNP151-common.next-base-RG.vcf" )
tabix_RG_cmd <- paste0("tabix -h -R ",epicv2_bed_RG_tsv," ",dbSNP_151_bgzip," > ",epicv2_int_dbSNP_RG_vcf)
base::system( tabix_RG_cmd )

# Before Range Reduction::
# > epicv2_int_dbSNP_RG_tib
# A tibble: 11,705 × 8
#
# > epicv2_int_dbSNP_RG_tib
# A tibble: 1,284 × 8
#
# > epicv2_int_dbSNP_RG_tib
# A tibble: 801 × 5
epicv2_int_dbSNP_RG_tib <- NULL
epicv2_int_dbSNP_RG_tib <- readr::read_tsv( file = epicv2_int_dbSNP_RG_vcf, 
                                         skip = 56, show_col_types = FALSE) %>% 
  magrittr::set_names( c("Chromosome_SNP_hg38","Coordinate_SNP_hg38","RS_ID",
                         "REF","ALT","QUAL","FILTER","INFO") ) %>%
  dplyr::filter( stringr::str_length(REF) == 1 ) %>%
  dplyr::filter( stringr::str_length(ALT) == 1 ) %>%
  dplyr::filter( REF != ALT ) %>%
  dplyr::filter( !RS_ID %>% stringr::str_detect(",") ) %>% 
  dplyr::distinct( Chromosome_SNP_hg38,Coordinate_SNP_hg38,RS_ID,REF,ALT )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Validation Using Existing Overlap::
#                                 EPICv2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

check_validation <- TRUE
check_validation <- FALSE
if ( check_validation ) {
  
  # Expectation Below for Overlapping EPICv2 Infinium I Trifecta Probes::
  #  > val_annoI_dbSNP_tib
  #  A tibble: 449 × 23
  val_annoI_dbSNP_tib <- NULL
  val_annoI_dbSNP_tib <- epicv2_int_dbSNP_RG_tib %>% 
    dplyr::select( Chromosome_SNP_hg38:ALT ) %>% 
    dplyr::inner_join( annoI_tib, by=c("RS_ID","REF","ALT") ) %>%
    dplyr::inner_join( epicv2_bed_RG_tib, 
                       by=c("Loci_ID","Chromosome_SNP_hg38"="Chromosome_BED_hg38"),
                       multiple = "all" ) %>% 
    dplyr::mutate( Beg_Dif = Coordinate_SNP_hg38 - Beg_BED_hg38,
                   End_Dif = Coordinate_SNP_hg38 - End_BED_hg38 )
  
  # Quick Print::
  val_annoI_dbSNP_tib %>% head() %>% as.data.frame()
  
  #
  # Investigate Distance by Strand::
  #.  Conclusion Use End_BED_hg38 as the Next Base Position!!!
  val_annoI_dbSNP_tib %>% 
    dplyr::group_by( Beg_Dif,Strand_FR,Strand_TB,Strand_CO ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  val_annoI_dbSNP_tib %>% 
    dplyr::group_by( End_Dif,Strand_FR,Strand_TB,Strand_CO ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  #
  # Validation: Unique Check
  #
  val_annoI_dbSNP_tib %>% dplyr::distinct( Loci_ID,RS_ID )
  val_annoI_dbSNP_tib %>% dplyr::distinct( Loci_ID,RS_ID,Probe_ID )
  val_annoI_dbSNP_tib %>% dplyr::distinct( Loci_ID,RS_ID,Probe_ID,Chromosome_SNP_hg38 )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Remapping Suffix Shorten Names::
#                             EPICv2 (annoI)
#
#
# LEFT OFF HERE::
#   [Done]: Validation Attempt:: Testing RG
#.  TBD:: Using lessons from below filter above to make a new AnnoI_v2
#.   - [Done]: Remove Indels, Multiple RS#'s...
#.   - Investigate MAF...
#.   - Also need to figure out the whole mapping back issue...
# IMPORTANT: Shouldn't run ind run idn and don't use i for SNP calling!
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Master Manifest to Monidfy::
#
man_v2_prbs_csv2 <- file.path( opt$top_path, "scratch/stable_docker_join_manifest/EPICv2-UCSC-v3/manifests/EPIC-8v2-0_A1.prb.csv.gz" )
epic_v2_all_tib2 <- NULL
epic_v2_all_tib2 <- readr::read_csv( file = man_v2_prbs_csv2, show_col_types = FALSE )

# Logic to the U field for annoS::
#
# > annoS_tib %>% dplyr::group_by( U, designType,REF, ALT ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
# A tibble: 4 × 4
# U     REF   ALT   Count
# <chr> <chr> <chr> <int>
# 1 ALT   C     T        10
# 2 ALT   G     A        10
# 3 REF   A     G        14
# 4 REF   T     C        10

epic_v2_cgn_tib2 <- NULL
epic_v2_cgn_tib2 <- epic_v2_all_tib2 %>%
  dplyr::mutate(
    Probe_ID = dplyr::case_when(
      Probe_ID == "rs10033147_BC11" ~ "rs10033147",
      Probe_ID == "rs6991394_BC11"  ~ "rs6991394",
      Probe_ID == "rs951295_BC11"   ~ "rs951295",
      # Probe_ID == "" ~ "",
      TRUE ~ Probe_ID
    )
  )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  Final Genomic Ranges (RDS) Generation::
#                             EPICv2 (annoS)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

write_annoS <- FALSE
write_annoS <- TRUE
if ( write_annoS ) {
  
  #
  # Get Manifest Designs::
  #
  epicv2_snp_tib <- NULL
  epicv2_snp_tib <- epic_v2_cgn_tib2 %>% 
    dplyr::filter( Probe_Type == "rs" ) %>% 
    dplyr::inner_join( 
      epicv2_man_tib %>% dplyr::select( IlmnID, Strand_FR, Infinium_Design,
                                        Infinium_Design_Type,
                                        CHR,MAPINFO,Probe_Type ) %>% 
        dplyr::rename( Full_ID = IlmnID ),
      by=c("Full_ID","Probe_Type")
    ) %>%
    dplyr::select( Probe_ID,Loci_ID,Full_ID, 
                   Strand_FR, # Strand_TB,Strand_CO,
                   Infinium_Design,
                   Infinium_Design_Type,
                   Rep_Num,
                   CHR,MAPINFO,Probe_Type ) %>%
    dplyr::rename(
      designType = Infinium_Design_Type
    ) %>%
    dplyr::mutate(
      strand = dplyr::case_when(
        Strand_FR == "F" ~ "+",
        Strand_FR == "R" ~ "-",
        TRUE ~ NA_character_
      )
    ) %>% clean_tib()

  # Need to load /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/snps/missing/missing.EPICv2_snps.19.csv
  #.  to get remaining ALT/REF Fields!!!
  #
  missAF_col <- NULL
  missAF_col <- readr::cols(
    Chromosome_SNP_hg38 = readr::col_character(),
    Coordinate_SNP_hg38 = readr::col_integer(),
    RS_ID               = readr::col_character(),
    REF                 = readr::col_character(),
    ALT                 = readr::col_character(),
    QUAL                = readr::col_character(),
    FILTER              = readr::col_character(),
    INFO                = readr::col_character()
  )
  missAF_tsv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/snps/missing/missing.EPICv2_snps.19.tsv.gz" )
  missAF_tib <- NULL
  missAF_tib <- readr::read_tsv( file = missAF_tsv, 
                                 col_names = names(missAF_col$cols), 
                                 col_types = missAF_col ) %>% 
    dplyr::select( RS_ID, Chromosome_SNP_hg38,Coordinate_SNP_hg38, REF,ALT ) %>%
    dplyr::rename( Loci_ID = RS_ID,
                   seqnames = Chromosome_SNP_hg38, 
                   start = Coordinate_SNP_hg38 ) %>% 
    dplyr::mutate( 
      end = start,
      seqnames = paste0("chr",seqnames)
    ) %>%
    dplyr::inner_join( epicv2_snp_tib, # %>% dplyr::select( Loci_ID,strand,designType ), 
                       by=c("Loci_ID",
                            "seqnames"="CHR",
                            "start"="MAPINFO") ) %>%
    dplyr::mutate( 
      # U = dplyr::case_when(
      #   REF == "C" & ALT == "T" ~ "ALT",
      #   REF == "G" & ALT == "A" ~ "ALT",
      #   REF == "A" & ALT == "G" ~ "REF",
      #   REF == "T" & ALT == "C" ~ "REF",
      #   TRUE ~ NA_character_ ),
      ALT2 = dplyr::case_when(
        REF == "C" ~ "T",
        REF == "G" ~ "A",
        REF == "A" ~ "G",
        REF == "T" ~ "C",
        TRUE ~ NA_character_
      ),
      U = dplyr::case_when(
        REF == "C" & ALT %>% stringr::str_detect("T") ~ "ALT",
        REF == "G" & ALT %>% stringr::str_detect("A") ~ "ALT",
        REF == "A" & ALT %>% stringr::str_detect("G") ~ "REF",
        REF == "T" & ALT %>% stringr::str_detect("C") ~ "REF",
        TRUE ~ NA_character_ ),
      rs = Loci_ID
    ) %>%
    dplyr::select( Loci_ID, seqnames, start, end, 
                   strand, rs, designType, U, REF, ALT2, ALT, Full_ID )
  
  #
  # Gather all the standard EPIC SNPs::
  #
  annoS_snp_all_tib <- NULL
  annoS_snp_all_tib <- annoS_tib %>% 
    dplyr::rename( Loci_ID = Probe_ID) %>% 
    dplyr::select( -width ) %>% clean_tib() %>%
    dplyr::inner_join( epicv2_snp_tib, 
                       by=c("Loci_ID",
                            "designType",
                            "strand",
                            # "end"="MAPINFO",
                            "seqnames"="CHR" ),
                       multiple = "all" ) %>% 
    dplyr::select( Probe_ID, seqnames, start, end, 
                   strand, rs, designType, U, REF, ALT, Full_ID ) %>%
    dplyr::rename( Loci_ID = Probe_ID ) %>%
    dplyr::mutate( ALT2 = ALT ) %>%
    dplyr::select( Loci_ID, seqnames, start, end, 
                   strand, rs, designType, U, REF, ALT2, ALT, Full_ID )
  
  epicv2_snp_all_tib <- NULL
  epicv2_snp_all_tib <- dplyr::bind_rows( missAF_tib, annoS_snp_all_tib ) %>% 
    dplyr::arrange( seqnames,start,end,Loci_ID )
  
  #
  # Write annoS::
  #
  #. TEST 1.0:: Attempting...
  #.   -  rs  = epicv2_snp_all_tib$Loci_ID %>% stringr::str_remove("_.*$")
  #.   -  ALT = epicv2_snp_all_tib$ALT2,
  #
  #. TEST 1.1:: NA
  #.   - Remove all probes with suffix...
  #
  epicv2_snp_all_tib <- epicv2_snp_all_tib %>% 
    dplyr::filter( !Loci_ID %>% stringr::str_detect("_") )

  #
  # IMPORANT:: THESE MATCH PEFECTLY::
  #
  match_v1v2_tib <- NULL
  match_v1v2_tib <- annoS_tib %>% 
    dplyr::inner_join( epicv2_snp_all_tib, 
                       by=c("Probe_ID"="Loci_ID", "seqnames", "start","end",
                            "strand","designType","rs","U","REF","ALT" ) )
  
  epicv2_annoS_dbSNP_rds <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/EPICv2/EPICv2.annoS.rds" )
  epicv2_annoS_dbSNP_grs <- NULL
  epicv2_annoS_dbSNP_grs <- GenomicRanges::GRanges(
    seqnames   = epicv2_snp_all_tib$seqnames, 
    strand     = epicv2_snp_all_tib$strand,
    
    # rs         = epicv2_snp_all_tib$Loci_ID,
    rs         = epicv2_snp_all_tib$Loci_ID %>% stringr::str_remove("_.*$"),
    designType = epicv2_snp_all_tib$designType,
    U          = epicv2_snp_all_tib$U,
    REF        = epicv2_snp_all_tib$REF,
    ALT        = epicv2_snp_all_tib$ALT2,
    # ALT        = epicv2_snp_all_tib$ALT,
    # Full_ID    = epicv2_snp_all_tib$Full_ID,
    
    IRanges::IRanges(
      start = epicv2_snp_all_tib$start,
      end   = epicv2_snp_all_tib$end,
      names = epicv2_snp_all_tib$Loci_ID
    )
  )
  readr::write_rds( x = epicv2_annoS_dbSNP_grs, 
                    file = epicv2_annoS_dbSNP_rds, 
                    compress = "gz" )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  Final Genomic Ranges (RDS) Generation::
#                             EPICv2 (annoI)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

write_manifest <- FALSE
write_manifest <- TRUE

if ( write_manifest ) {
  #
  # LEFT OFF HERE::
  #.  TBD:: Incorporate these changes: missS2_tib
  #.  TBD:: Remove this probe (rs1495031) from EPICv2.annoS.rds 
  #
  
  epic_v2_cgn_sum2 <- NULL
  epic_v2_cgn_sum2 <- epic_v2_cgn_tib2 %>%
    dplyr::select( Full_ID,Probe_ID, M,U,mask,Probe_Type ) %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) epic_v2_cgn_sum2 %>% print( n=base::nrow(epic_v2_cgn_sum2) )
  
  #
  # Write Mapping EPICv2 File::
  #
  #.  > epic_v2_cgn_map2
  #.  A tibble: 937,690 × 6
  epic_v2_cgn_map2_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/map/EPICv2.cgn-map.csv.gz" )
  epic_v2_cgn_map2 <- NULL
  epic_v2_cgn_map2 <- epic_v2_cgn_tib2 %>%
    dplyr::select( Full_ID,Probe_ID, M,U,mask,Probe_Type )
  readr::write_csv( x = epic_v2_cgn_map2, file = epic_v2_cgn_map2_csv )
  # docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/map/EPICv2.cgn-map.csv.gz dd00d8b92b93:/repo/Infinium_Methylation_Workhorse/dat/map/
  
  #
  # Previous Modification Step::
  #.  TBD:: Update rs#'s
  #
  epic_v2_cgn_csv3 <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/manifest/EPICv2/EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
  epic_v2_cgn_tib3 <- NULL
  epic_v2_cgn_tib3 <- epic_v2_cgn_tib2 %>% 
    dplyr::select( Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,Probe_Type,
                   Probe_Source,Next_Base,Probe_Design,mask,Full_ID ) %>%
    dplyr::mutate( Probe_Source = "EPIC-A1" )
  readr::write_csv( x = epic_v2_cgn_tib3, file = epic_v2_cgn_csv3 )
  # docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/manifest/EPICv2/EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz dd00d8b92b93:/repo/Infinium_Methylation_Workhorse/dat/manifest/core/
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  Final Genomic Ranges (RDS) Generation::
#                             EPICv2 (annoI)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

write_annoI <- FALSE
write_annoI <- TRUE
if ( write_annoI ) {
  
  epicv2_cpg_dbSNP_RG_tib <- NULL
  epicv2_cpg_dbSNP_RG_tib <- epicv2_int_dbSNP_RG_tib %>%
    dplyr::inner_join( epicv2_bed_RG_tib, 
                       by=c("Chromosome_SNP_hg38"="Chromosome_BED_hg38",
                            "Coordinate_SNP_hg38"="End_BED_hg38"),
                       multiple = "all" ) %>%
    dplyr::mutate(
      strand = dplyr::case_when(
        Strand_FR == "F" ~ "+",
        Strand_FR == "R" ~ "-",
        TRUE ~ NA_character_
      ),
      Chromosome_SNP_hg38 = paste0("chr",Chromosome_SNP_hg38)
    ) %>%
    dplyr::rename( 
      Full_ID = Probe_ID
    ) %>%
    dplyr::arrange( Chromosome_SNP_hg38,Coordinate_SNP_hg38,RS_ID )
  
  #
  # Set Manifest Used Names::
  #
  epicv2_annoI_dbSNP_tib <- NULL
  epicv2_annoI_dbSNP_tib <- dplyr::inner_join(
    dplyr::select( epicv2_cpg_dbSNP_RG_tib, Probe_ID, 
                   Chromosome_SNP_hg38,Coordinate_SNP_hg38,
                   strand,REF,ALT,RS_ID ) %>% dplyr::rename( Full_ID = Probe_ID ),
    dplyr::select( epic_v2_cgn_tib2, Full_ID,Probe_ID,Probe_Type ),
    by=c("Full_ID") #, multiple = "all"
  ) %>%
    dplyr::arrange( Chromosome_SNP_hg38,Coordinate_SNP_hg38,Probe_ID ) %>%
    dplyr::rename( )
  
  epicv2_cgn_dbSNP_RG_sum <- NULL
  epicv2_cgn_dbSNP_RG_sum <- epicv2_annoI_dbSNP_tib %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) epicv2_cgn_dbSNP_RG_sum %>% print( n=base::nrow(epicv2_cgn_dbSNP_RG_sum) )
  
  # Quick Sanity Check::
  #. epicv2_annoI_dbSNP_tib %>% dplyr::filter( Probe_ID %>% stringr::str_detect("_") )
  #
  # Quick Print::
  #  epicv2_cpg_dbSNP_RG_tib %>% head() %>% as.data.frame()
  #  epicv2_annoI_dbSNP_tib %>% head() %>% as.data.frame()
  
  #
  # Write annoI::
  #
  #. TEST 1.1:: NA
  #.   - Remove all probes with suffix...
  #
  epicv2_annoI_dbSNP_tib <- epicv2_annoI_dbSNP_tib %>% 
    dplyr::filter( !Probe_ID %>% stringr::str_detect("_") )
  
  #
  # IMPORANT:: THESE MATCH PEFECTLY:: NOT CHECKED YET!!!!
  #
  match_v1v2_annoI_tib <- NULL
  match_v1v2_annoI_tib <- annoI_tib %>% 
    dplyr::inner_join( epicv2_annoI_dbSNP_tib, 
                       by=c("Probe_ID"="Loci_ID", "seqnames", "start","end",
                            "strand","designType","rs","U","REF","ALT" ) )
  
  epicv2_annoI_dbSNP_rds <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/EPICv2/EPICv2.annoI.rds" )
  epicv2_annoI_dbSNP_grs <- NULL
  epicv2_annoI_dbSNP_grs <- GenomicRanges::GRanges( 
    seqnames   = paste0("chr",epicv2_annoI_dbSNP_tib$Chromosome_SNP_hg38), 
    strand     = epicv2_annoI_dbSNP_tib$strand,
    designType = "I",
    In.band    = "REF",
    REF        = epicv2_annoI_dbSNP_tib$REF,
    ALT        = epicv2_annoI_dbSNP_tib$ALT,
    rs         = epicv2_annoI_dbSNP_tib$RS_ID,
    
    IRanges::IRanges(
      start = epicv2_annoI_dbSNP_tib$Coordinate_SNP_hg38,
      end   = epicv2_annoI_dbSNP_tib$Coordinate_SNP_hg38,
      names = epicv2_annoI_dbSNP_tib$Probe_ID
    )
  )
  readr::write_rds( x = epicv2_annoI_dbSNP_grs, 
                    file = epicv2_annoI_dbSNP_rds, 
                    compress = "gz" )
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


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
par$prgm_tag <- 'stable_analysis_EPICv2_docker_snps'
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
par$version <- 12

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

# opt$pdf_path = safe_mkdir( dir = file.path( opt$out_path,"pdf") )
# opt$csv_path = safe_mkdir( dir = file.path( opt$out_path,"csv") )
# opt$rds_path = safe_mkdir( dir = file.path( opt$out_path,"rds") )
opt$snv_path = safe_mkdir( dir = file.path( opt$out_path,"snv") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Quick RDS Conversion::
#
# This is to provide Steven and team with a text version of the annoI/S
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {

  convert_rds_to_csv <- TRUE
  convert_rds_to_csv <- FALSE
  
  if ( convert_rds_to_csv ) {
    dat_dir <- file.path( "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/images/Infinium_Methylation_Workhorse_Centos.v.1.11.15.1.p.0.6.2/Infinium_Methylation_Workhorse/dat/dbSNP/b151" )
    
    annoI_b4_rds <- file.path( dat_dir, "EPIC-B4/EPIC-B4.annoI.rds" )
    annoS_b4_rds <- file.path( dat_dir, "EPIC-B4/EPIC-B4.annoS.rds" )
    annoI_b4_csv <- file.path( dat_dir, "EPIC-B4/EPIC-B4.annoI.csv" )
    annoS_b4_csv <- file.path( dat_dir, "EPIC-B4/EPIC-B4.annoS.csv" )
    
    annoI_b4_tib <- annoI_b4_rds %>% readr::read_rds() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
    annoS_b4_tib <- annoS_b4_rds %>% readr::read_rds() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
    readr::write_csv( x = annoI_b4_tib, file = annoI_b4_csv )
    readr::write_csv( x = annoS_b4_tib, file = annoS_b4_csv )
    
    annoI_a1_rds <- file.path( dat_dir, "EPIC-A1/EPIC-A1.annoI.rds" )
    annoS_a1_rds <- file.path( dat_dir, "EPIC-A1/EPIC-A1.annoS.rds" )
    annoI_a1_csv <- file.path( dat_dir, "EPIC-A1/EPIC-A1.annoI.csv" )
    annoS_a1_csv <- file.path( dat_dir, "EPIC-A1/EPIC-A1.annoS.csv" )
    
    annoI_a1_tib <- annoI_a1_rds %>% readr::read_rds() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
    annoS_a1_tib <- annoS_a1_rds %>% readr::read_rds() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
    readr::write_csv( x = annoI_a1_tib, file = annoI_a1_csv )
    readr::write_csv( x = annoS_a1_tib, file = annoS_a1_csv )
  }

  covnert_csv_to_rds <- FALSE
  covnert_csv_to_rds <- TRUE
  if ( covnert_csv_to_rds ) {
    #
    # Convert new annoI and annoS CSV to RDS
    #
    new_dat_dir <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/manifest/EPICv2/steven-test-swap" )
    
    annoI_a1_csv <- file.path( new_dat_dir, "EPIC-A1/EPIC-A1.annoI.csv.gz" )
    annoS_a1_csv <- file.path( new_dat_dir, "EPIC-A1/EPIC-A1.annoS.csv.gz" )
    
    annoI_a1_rds <- file.path( new_dat_dir, "EPIC-A1/EPIC-A1.annoI.rds" )
    annoS_a1_rds <- file.path( new_dat_dir, "EPIC-A1/EPIC-A1.annoS.rds" )

    # Load CSV::
    annoI_a1_tib <- NULL
    annoI_a1_tib <- readr::read_csv( file = annoI_a1_csv, show_col_types = FALSE )
    annoS_a1_tib <- NULL
    annoI_a1_tib <- readr::read_csv( file = annoS_a1_csv, show_col_types = FALSE )
    
    # Convert To Genomic Range::
    annoI_a1_grs <- NULL
    annoI_a1_grs <- GenomicRanges::GRanges( 
      seqnames = Rle( annoI_a1_tib$seqnames ),
      strand   = Rle( annoI_a1_tib$strand ),
      
      designType = annoI_a1_tib$designType,
      In.band    = annoI_a1_tib$In.band,
      REF = annoI_a1_tib$REF,
      ALT = annoI_a1_tib$ALT,
      rs  = annoI_a1_tib$rs,
      
      IRanges(start = annoI_a1_tib$start, 
              width = 1,
              # end   = annoI_a1_tib$end, 
              names = annoI_a1_tib$Probe_ID )
    )
    
    annoS_a1_grs <- NULL
    annoS_a1_grs <- GenomicRanges::GRanges( 
      seqnames = Rle( annoS_a1_tib$seqnames ),
      strand   = Rle( annoS_a1_tib$strand ),
      
      rs  = annoS_a1_tib$rs,
      designType = annoS_a1_tib$designType,
      U   = annoS_a1_tib$U,
      REF = annoS_a1_tib$REF,
      ALT = annoS_a1_tib$ALT,
      
      IRanges(start = annoS_a1_tib$start, 
              width = 1,
              # end   = annoS_a1_tib$end, 
              names = annoS_a1_tib$Probe_ID )
    )
    
    # Save RDS::
    readr::write_rds( x = annoI_a1_grs, file = annoI_a1_rds )
    readr::write_rds( x = annoS_a1_grs, file = annoS_a1_rds )
    
    # Copy to Docker::
    # docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/manifest/EPICv2/steven-test-swap/EPIC-A1/EPIC-A1.annoI.rds 04253273b6cf:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-A1/
    # docker cp /Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/manifest/EPICv2/steven-test-swap/EPIC-A1/EPIC-A1.annoS.rds 04253273b6cf:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-A1/

    # Test Docker::
    
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Final Manifest Swap::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  final_man_swap <- FALSE
  final_man_swap <- TRUE
  
  if ( final_man_swap ) {
    
    doc_dir <- "/Users/bbarnes/Documents/Projects/EPIC_v2/docker"

    new_man_dir <- file.path( doc_dir, "dat/EPICv2_VA_MVP-Data-2023_09_05" )
    new_man_csv <- file.path( new_man_dir, "EPICv2_VA_MVP_A1-Builder.body-sub.csv" )
    new_man_tib <- readr::read_csv( file = new_man_csv, show_col_types = FALSE )
    
    doc_man_dir <- file.path( doc_dir, "images/Infinium_Methylation_Workhorse_Centos.v.1.11.15.1.p.0.6.2/Infinium_Methylation_Workhorse/dat/manifest/core" )
    doc_man_csv <- file.path( doc_man_dir, "EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
    doc_man_tib <- readr::read_csv( file = doc_man_csv, show_col_types = FALSE )
    
    # new_man_tib: A tibble: 936,382 × 21
    # doc_man_tib: A tibble: 937,690 × 11
    # 937690 - 936382 = 1308
    
    dif_man_tib <- NULL
    dif_man_tib <- new_man_tib %>% 
      dplyr::anti_join( doc_man_tib, c("IlmnID"="Full_ID") ) 
    
    add_man_tib <- NULL
    add_man_tib <- dif_man_tib %>% 
      dplyr::rename( 
        Full_ID = IlmnID, 
        Probe_ID = Name, 
        U = AddressA_ID, 
        M = AddressB_ID,
        COLOR_CHANNEL = Color_Channel ) %>%
      dplyr::mutate( 
        DESIGN = dplyr::case_when(
          !is.na(U) &  is.na(M) ~ "II",
          !is.na(U) & !is.na(M) ~ "I",
          TRUE ~ NA_character_ ),
        Probe_Design = dplyr::case_when(
          DESIGN ==  "I" ~ 1.0,
          DESIGN == "II" ~ 2.0,
          TRUE ~ NA_real_ ) %>% as.integer(),
        COLOR_CHANNEL = dplyr::case_when(
          DESIGN == "II" ~ "Both",
          TRUE ~ COLOR_CHANNEL ),
        Probe_Source = "EPIC-A1"
      ) %>%
      dplyr::select( Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col,
                     Probe_Type, Probe_Source, Next_Base, Probe_Design, Full_ID )
    
    
    fin_man_csv <- file.path( doc_man_dir, "EPIC-A2.manifest.sesame-base.cpg-sorted.csv.gz" )
    fin_man_tib <- NULL
    fin_man_tib <- dplyr::bind_rows( doc_man_tib,add_man_tib ) %>% 
      dplyr::arrange( Probe_ID )
    readr::write_csv( x = fin_man_tib, file = fin_man_csv )
                                   
  }

}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Quick Math::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# 2^30 = 1,073,741,824
# 2^29 =   536,870,912
#
# 1.1^181 = 31,051,030
#  1.1^31 =         19.19434
#
# v1.1 = 1073741824 + 31051030 = 1104792854
# v2.1 =  536870912 + 19.19434 =  536870931
# v2.2 = 1073741824 + 19.19434 = 1073741843
#
# v2.1 / v1.1 =  536870931 / 1104792854 = 0.4859471
# v2.2 / v1.1 = 1073741843 / 1104792854 = 0.9718943
#

#
# New: rs = 64 - 59 = 5
# New: sv = 17 
# New: cg = 877
#


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Pre-processing::
#
#                    Compare/Merge All Allele Frequencies::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

db_snv_tib <- NULL
db_snv_csv <- NULL
db_snv_csv <- file.path( opt$snv_path, "dbSNP151_AF.csv.gz" )

if ( file.exists(db_snv_csv) ) {
  
  db_snv_tib <- readr::read_csv( file = db_snv_csv, show_col_types = FALSE )
  
} else {
  va_ids_tib <- NULL
  va_ids_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/VA_SNP_Selection/methGenotypeProbes.csv" )
  va_ids_tib <- readr::read_csv( file = va_ids_csv, show_col_types = FALSE ) %>%
    dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>%
    tidyr::separate( ID, into=c("Chr","Pos","Alt","Ref"), 
                     sep=":", remove = TRUE, convert = TRUE ) %>%
    dplyr::mutate( Source = "VA",
                   Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"),
                   Probe_Type = Probe_ID %>% stringr::str_sub(1,2)
    ) %>% dplyr::select( Source, Probe_ID, Loci_ID, Probe_Type, dplyr::everything() )
  
  #
  # [Done]: Re-run: stable_scratch_EPICv2_docker_snps.R
  #.       - Obtain new dbSNP.tabix files...
  # [Done]: Write a function that can parse a tabix correctly...
  # [TBD]: Distribution of REF AF: [ va,v1,v2,{v2!v1} ] x [ {rs,nv}, {cg,ch} ]
  #
  
  v1_body_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv1/rank_historic_probes/EPICv2-UCSC-v5.rank_historic_probes.csv.gz" )
  v2_body_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv2/rank_historic_probes/EPICv2-UCSC-v5.rank_historic_probes.csv.gz" )
  
  v1_dbSNP_vcf <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv1/intersect_bed_tabix/EPICv2-UCSC-v5.intersect_bed_tabix.tab.bed.gz" )
  v2_dbSNP_vcf <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv2/intersect_bed_tabix/EPICv2-UCSC-v5.intersect_bed_tabix.tab.bed.gz" )
  
  v1_body_tib <- NULL
  v1_body_tib <- readr::read_csv( file = v1_body_csv, show_col_types = FALSE )
  v2_body_tib <- NULL
  v2_body_tib <- readr::read_csv( file = v2_body_csv, show_col_types = FALSE )
  
  # [match_src = FALSE]: A tibble: 158,676 × 15
  # [match_src = TRUE]:  A tibble: 133,167 × 44
  v1_snv_tib <- NULL
  v1_snv_tib <- parse_dbSNP_vcf( 
    snp_vcf = v1_dbSNP_vcf, 
    src_tib = v1_body_tib %>% dplyr::mutate( Chromosome = stringr::str_remove(Chromosome, "^chr") ),
    
    ref_vcf   = TRUE, 
    com_len   = 56, 
    rm_chr0   = TRUE,
    strip_chr = TRUE,
    snps_only = TRUE,
    parse_caf = TRUE,
    # match_src = FALSE,
    match_src = TRUE,
    min_maf   = 0.0,
    
    out_dir = file.path( opt$out_path, "EPICv1" ),
    run_tag = opt$run_name, 
    reload = opt$reload, 
    reload_min = 10,
    write_out = FALSE,
    
    vb=vb,vt=vt+1,tc=tc,tt=tt )
  # v1_snv_tib %>% dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE )
  
  # [match_src = FALSE]: A tibble: 156,286 × 15
  # [match_src = TRUE]:  A tibble: 132,604 × 46
  v2_snv_tib <- NULL
  v2_snv_tib <- parse_dbSNP_vcf( 
    snp_vcf = v2_dbSNP_vcf, 
    src_tib = v2_body_tib %>% dplyr::mutate( Chromosome = stringr::str_remove(Chromosome, "^chr") ),
    
    ref_vcf   = TRUE, 
    com_len   = 56, 
    rm_chr0   = TRUE,
    strip_chr = TRUE,
    snps_only = TRUE,
    parse_caf = TRUE,
    # match_src = FALSE,
    match_src = TRUE,
    min_maf   = 0.0,
    
    out_dir = file.path( opt$out_path, "EPICv1" ),
    run_tag = opt$run_name, 
    reload = opt$reload, 
    reload_min = 10,
    write_out = FALSE,
    
    vb=vb,vt=vt+1,tc=tc,tt=tt )
  # v2_snv_tib %>% dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE )
  
  # snv_tib <- TRUE
  # snv_tib <- dplyr::bind_rows( v1_snv_tib,v2_snv_tib ) %>%
  #   dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE ) %>%
  #   dplyr::arrange( Chromosome,Coordinate_SNP,SNP_ID, -MAF )
  
  # A tibble: 175,479 × 10
  db_snv_tib <- TRUE
  db_snv_tib <- dplyr::bind_rows( 
    v1_snv_tib %>% dplyr::select( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, RAF,MAF, Probe_ID,Full_ID,Loci_ID),
    v2_snv_tib %>% dplyr::select( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, RAF,MAF, Probe_ID,Full_ID,Loci_ID)
  ) %>%
    dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE ) %>%
    dplyr::arrange( Chromosome,Coordinate_SNP,SNP_ID, -MAF )
  
  readr::write_csv( x = db_snv_tib, file = db_snv_csv )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Pre-processing::
#
#                        Load Truth Sample Sheets::
#                         Load Auto Sample Sheets::
#                             Merging VCFs::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# [TBD]:
#  - Truth Comparisons Sample Calling
#  - Comparisons Beta (r2/dB) & SNVs (gt) x EPICv1/v2
#

top_dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2" )

true_wd_ssh_tib <- NULL
true_wd_ssh_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/sampleSheets/EPICv2_GSE228820_Wanding.base.sampleSheet.v5.csv" )
true_wd_ssh_tib <- readr::read_csv( file = true_wd_ssh_csv, show_col_types = FALSE ) %>%
  clean_tib() %>%
  dplyr::mutate( Sample_Source = "GSE228820" )

true_ep_ssh_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/formatted/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz" )
true_ep_ssh_tib <- NULL
true_ep_ssh_tib <- readr::read_csv( file = true_ep_ssh_csv, show_col_types = FALSE ) %>%
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
                       Concentration,Sample_Group,Sample_Source )

true_ssh_tib <- NULL
true_ssh_tib <- dplyr::bind_rows( true_wd_ssh_tib,true_ep_ssh_tib ) %>% 
  dplyr::distinct( Sentrix_Name, .keep_all = TRUE )

true_ssh_sum <- NULL
true_ssh_sum <- true_ssh_tib %>% 
  dplyr::group_by( Sample_Source,Sample_Base,Sample_ID,Sample_Class ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) true_ssh_sum %>% print( n=base::nrow(true_ssh_sum) )

opt$min_call_rate <- 75

ssh_tib <- NULL
ssh_tib <- file_list( path    = top_dat_path, 
                      prefix  = top_dat_path,
                      suffix  = "_AutoSampleSheet.csv.gz", 
                      pattern = "_AutoSampleSheet.csv.gz$",
                      recursive = TRUE ) %>% 
  # head( n=50 ) %>%
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows() %>%
  dplyr::select( Sentrix_Name, 
                 detect_version,
                 cg_calls_pass_perc_1,
                 AutoSample_R2_Key_1,AutoSample_dB_Key_1,
                 AutoSample_R2_Val_1,AutoSample_dB_Val_1, 
                 AutoSample_dB_Cnt_1 ) %>%
  dplyr::full_join( true_ssh_tib, by=c("Sentrix_Name") ) %>%
  dplyr::mutate(
    Sample_Requeue = dplyr::case_when( 
      cg_calls_pass_perc_1 >= opt$min_call_rate ~ FALSE, 
      cg_calls_pass_perc_1 <  opt$min_call_rate ~ TRUE,
      TRUE ~ NA ),
    
    Sample_Base = dplyr::case_when(
      is.na(Sample_Base) ~ paste0( "X",dplyr::row_number() ), 
      TRUE ~ Sample_Base ),
    Sample_ID = dplyr::case_when(
      # is.na(Sample_Group) ~ "X",
      # is.na(Sample_Group) ~ "I",
      is.na(Sample_Group) ~ Sample_Base,
      TRUE ~ Sample_ID ),
    Sample_Class = dplyr::case_when(
      is.na(Sample_Group) ~ "GSE222131",
      TRUE ~ Sample_Group ),
    Concentration = dplyr::case_when(
      is.na(Sample_Group) ~ 0.0,
      TRUE ~ Concentration ) %>% as.integer(),
    Sample_Source = dplyr::case_when(
      is.na(Sample_Group) ~ "GSE222131",
      TRUE ~ Sample_Source ),
    Sample_Group = dplyr::case_when(
      is.na(Sample_Group) ~ "GSE222131",
      TRUE ~ Sample_Group )
  ) %>%
  dplyr::filter( !is.na(cg_calls_pass_perc_1) ) %>%
  dplyr::filter( !is.na(Sample_Requeue) ) %>%
  dplyr::group_by( Sample_Base,detect_version ) %>% 
  dplyr::mutate(
    Rep_Num = dplyr::row_number(),
    PPP_Int = cg_calls_pass_perc_1 %>% as.integer(),
    Plot_ID = paste( detect_version,Sample_ID,Concentration,Rep_Num,PPP_Int, sep="_")
    # Plot_ID = paste0( detect_version,"_",Sample_ID,Rep_Num,"_",PPP_Int)
  ) %>% dplyr::ungroup() %>%
  dplyr::select( Sentrix_Name,Plot_ID, dplyr::everything() )

ssh_sum <- NULL
ssh_sum <- ssh_tib %>% 
  dplyr::group_by( Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) ssh_sum %>% print( n=base::nrow(ssh_sum) )

# ssh_tib %>% dplyr::group_by( Sample_Class ) %>% dplyr::summarise( Count=n(), .groups = "drop" )

ssh_sum2 <- ssh_tib %>% 
  dplyr::group_by( detect_version,Sample_Source,Sample_Group,Sample_Class,Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)

#
# [TBD]: Build Experiment Groups
#  - Alpha, cg_calls_pass_perc_1 > 75 + GSE228820(Sample_Use=Correlation)
#. - Alpha, cg_calls_pass_perc_1 > 100
#  - 
#

ssh_tib0 <- NULL
ssh_tib0 <- dplyr::bind_rows(
  ssh_tib %>% dplyr::filter( Sample_Source == "Alpha" ),
  ssh_tib %>% dplyr::filter( Sample_Source == "GSE228820" & Sample_Use == "Correlation" ),
) %>% dplyr::filter( Sample_Base != "NIH3T3" )

ssh_tib1 <- NULL
ssh_tib1 <- ssh_tib0 %>% dplyr::filter( Sample_Requeue == FALSE, Sample_Group == "CellLine" )

exp_ssh_tibs <- NULL
# exp_ssh_tibs[["CellLine_Full"]] <- ssh_tib0 %>% dplyr::filter( Sample_Requeue == FALSE, Sample_Group == "CellLine" )
# exp_ssh_tibs[["CellLine_Pass"]] <- ssh_tib0 
# exp_ssh_tibs[["MeTritration"]] 

exp_ssh_tibs[["v1"]] <- ssh_tib1 %>% dplyr::filter( detect_version == "B4" ) %>%
  dplyr::add_count( Sample_Base, name="Sam_Cnt" ) %>% dplyr::filter( Sam_Cnt > 1 )
exp_ssh_tibs[["v2"]] <- ssh_tib1 %>% dplyr::filter( detect_version == "A1" ) %>%
  dplyr::add_count( Sample_Base, name="Sam_Cnt" ) %>% dplyr::filter( Sam_Cnt > 1 )
exp_ssh_tibs[["v12"]] <- ssh_tib1 %>%
  dplyr::add_count( Sample_Base, name="Sam_Cnt" ) %>% dplyr::filter( Sam_Cnt > 1 )

exp_ssh_tibs[["v1"]] %>%
  dplyr::group_by( Sample_Requeue, detect_version,Sample_Source,Sample_Group,Sample_Class,Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
exp_ssh_tibs[["v2"]] %>%
  dplyr::group_by( Sample_Requeue, detect_version,Sample_Source,Sample_Group,Sample_Class,Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
exp_ssh_tibs[["v12"]] %>%
  dplyr::group_by( Sample_Requeue, detect_version,Sample_Source,Sample_Group,Sample_Class,Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)

#
# Counts Calculation:
#
exp_ssh_tibs[["v1"]] %>%
  dplyr::group_by( Sample_Requeue, detect_version,Sample_Source,Sample_Group,Sample_Class,Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% dplyr::summarise( sum(Count) )
exp_ssh_tibs[["v2"]] %>%
  dplyr::group_by( Sample_Requeue, detect_version,Sample_Source,Sample_Group,Sample_Class,Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% dplyr::summarise( sum(Count) )
exp_ssh_tibs[["v12"]] %>%
  dplyr::group_by( Sample_Requeue, detect_version,Sample_Source,Sample_Group,Sample_Class,Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% dplyr::summarise( sum(Count) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                            Build VCF Lists::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

vcf_list <- NULL
vcf_list <- file_list( path    = top_dat_path, 
                       prefix  = top_dat_path, 
                       suffix  = "_EPIC_.*_raw.snps.vcf", 
                       pattern = "_raw.snps.vcf$",
                       recursive = TRUE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#.                         Compare Overall MAFs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

va_ids_tib <- NULL
va_ids_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/VA_SNP_Selection/methGenotypeProbes.csv" )
va_ids_tib <- readr::read_csv( file = va_ids_csv, show_col_types = FALSE ) %>%
  dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>%
  tidyr::separate( ID, into=c("Chr","Pos","Alt","Ref"), 
                   sep=":", remove = TRUE, convert = TRUE ) %>%
  dplyr::mutate( Source = "VA",
                 Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"),
                 Probe_Type = Probe_ID %>% stringr::str_sub(1,2)
  ) %>% dplyr::select( Source, Probe_ID, Loci_ID, Probe_Type, dplyr::everything() )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                         Fingerprint:: All-SNVs
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

gtc_vec <- NULL
gtc_vec <- c( 0, 20 )

if ( FALSE ) {
  
  opt$single <- TRUE
  opt$single <- FALSE
  
  opt$write_out <- TRUE
  opt$write_out <- FALSE
  
  opt$plot_heat <- FALSE
  opt$plot_heat <- TRUE
  
  sam_vec <- NULL
  sam_vec <- ssh_tib %>% 
    dplyr::filter( Sample_Base != "PROMEGAMOUSE" ) %>% 
    dplyr::pull( Sample_Base ) %>% unique()
  
  prd_vec <- NULL
  prd_vec <- unique(ssh_tib$detect_version)
  
  grp_vec <- NULL
  grp_vec <- ssh_tib %>% 
    dplyr::filter( Sample_Class != "NegativeMouse" ) %>% 
    dplyr::pull( Sample_Class ) %>% unique()
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Build Evaluation Sets::
#                           [ VA,V1,V2,v1,v2,V3 ]
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sub_lst <- NULL
sub_lst <- list()

sub_lst[["VA"]] <- NULL
sub_lst[["VA"]] <- va_ids_tib %>% dplyr::pull( Probe_ID )

sub_lst[["V1"]] <- NULL
sub_lst[["V1"]] <- c( 
  readr::read_rds( file = file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoS.rds" ) ) %>%
    as.data.frame() %>% 
    tibble::rownames_to_column( var = "Probe_ID" ) %>% 
    tibble::as_tibble() %>% dplyr::pull(Probe_ID),
  readr::read_rds( file = file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoI.rds" ) ) %>%
    as.data.frame() %>% 
    tibble::rownames_to_column( var = "Probe_ID" ) %>% 
    tibble::as_tibble() %>% dplyr::pull(Probe_ID)
) %>% unique()

sub_lst[["V2"]] <- NULL
sub_lst[["V2"]] <- c( 
  readr::read_rds( file = file.path(opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/docker/dat/EPICv2/EPICv2.annoI.rds") ) %>%
    as.data.frame() %>% 
    tibble::rownames_to_column( var = "Probe_ID" ) %>% 
    tibble::as_tibble() %>% dplyr::pull(Probe_ID),
  readr::read_rds( file = file.path(opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/docker/dat/EPICv2/EPICv2.annoS.rds") ) %>%
    as.data.frame() %>% 
    tibble::rownames_to_column( var = "Probe_ID" ) %>% 
    tibble::as_tibble() %>% dplyr::pull(Probe_ID)
) %>% unique()

sub_lst[["V3_all"]] <- NULL
sub_lst[["V3_all"]] <- c( sub_lst[["V1"]],sub_lst[["V2"]] ) %>% unique()

sub_lst[["v1_unq"]] <- NULL
sub_lst[["v1_unq"]] <- sub_lst[["V1"]] %>% setdiff( sub_lst[["V2"]] ) %>% unique()

sub_lst[["v2_unq"]] <- NULL
sub_lst[["v2_unq"]] <- sub_lst[["V2"]] %>% setdiff( sub_lst[["V1"]] ) %>% unique()

sub_lst[["v6"]] <- NULL
sub_lst[["v6"]] <- intersect( sub_lst[["VA"]], sub_lst[["V2"]] )

sub_lst %>% lapply( length )

vcf_list %>% length()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#.                          Evaluate:: by Params
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( TRUE ) {
  snvs_dats <- NULL
  for ( gts_min in gtc_vec ) {
    gts_key <- paste0("gts",gts_min)
    cur_out_path <- safe_mkdir( file.path( opt$out_path, "partitions", "byFullData" ) )

    snvs_gtc <- NULL
    snvs_gtc <- vcf_list %>% # head() %>%
      analyze_snvs( run_tag="All",
                    ssh_tib = ssh_tib,
                    sub_vec = NULL,
                    gts_min = gts_min, 
                    fval = 5, 
                    jcnt = 0, 
                    uval = TRUE,
                    outDir = cur_out_path,
                    write_out = FALSE, 
                    plot_heat = FALSE )
    
    snvs_dats[[gts_key]] <- NULL
    snvs_dats[[gts_key]] <- snvs_gtc
  }
  # Proof their different::
  # snvs_dats$gts0 - snvs_dats$gts20
  
  snvs_pers <- NULL
  for ( gts_key in names(snvs_dats) ) {
    
    snvs_gtc <- NULL
    snvs_gtc <- snvs_dats[[gts_key]]
    
    for ( exp_key in names(exp_ssh_tibs) ) {
      sam_ssh_tibs <- NULL
      sam_ssh_tibs <- exp_ssh_tibs[[exp_key]] %>% split(.$Sample_Base)
      
      for ( sam_key in names(sam_ssh_tibs) ) {
        if (p2) cat(glue::glue("{pmssg} gts={gts_key}, exp={exp_key}, sam={sam_key}...{RET}"))

        cur_ssh_tib <- NULL
        cur_ssh_tib <- sam_ssh_tibs[[sam_key]]
        
        sam_snvs_gtc <- NULL
        # sam_snvs_gtc <- snvs_gtc[ sub_vec, sam_vec ]
        sam_snvs_gtc <- snvs_gtc[ , cur_ssh_tib$Plot_ID ]
        
        per_snvs_tib <- NULL
        per_snvs_tib <- gtc_snv_to_performance( x = sam_snvs_gtc ) %>%
          dplyr::mutate( Gts_Min = gts_key, 
                         Exp_Set = exp_key, 
                         Sam_Set = sam_key )
        
        snvs_pers <- dplyr::bind_rows( snvs_pers,per_snvs_tib )
        
        # break
      }
      #cbreak
    }
    # break
  }
}

snvs_pers_tab <- NULL
snvs_pers_tab <- snvs_pers %>% 
  dplyr::filter( Gts_Min == "gts20" ) %>%
  dplyr::mutate( Mat_Per = round(Match_Cnt/Total_Cnt,3), 
                 Exp_Str = paste(Exp_Set,Gts_Min,Sam_Set, sep="_") ) %>% 
  dplyr::select( Probe_ID,Exp_Str,Mat_Per ) %>%
  tidyr::pivot_wider( id_cols = c(Probe_ID), 
                      names_from = c(Exp_Str), 
                      values_from = c(Mat_Per) )

snv_full_tab <- NULL
snv_full_tab <- snvs_pers_tab %>%
  dplyr::mutate(
    Probe_Source = dplyr::case_when(
      Probe_ID %in% sub_lst[["v6"]] ~ "VA_60",
      Probe_ID %in% sub_lst[["VA"]] ~ "VA_151",
      Probe_ID %in% sub_lst[["V1"]] & Probe_ID %in% sub_lst[["V2"]] ~ "V3_390",
      Probe_ID %in% sub_lst[["V2"]] ~ "V2_889",
      Probe_ID %in% sub_lst[["V1"]] ~ "V1_267",
      
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select( Probe_ID,Probe_Source, dplyr::everything() )

snv_full_sum <- NULL
snv_full_sum <- snv_full_tab %>% 
  dplyr::group_by( Probe_Source ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

snv_full_csv <- file.path( opt$out_path, "snv_cellLine_agreement.csv.gz" )
readr::write_csv( x = snv_full_tab, file = snv_full_csv )



# snv_full_tab %>% dplyr::filter( Probe_Source == "VA_60" ) %>% dplyr::select( 1,2,37:44 ) %>% dplyr::select(-Probe_Source) %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix() %>% matrixStats::rowMeans2()

#
# 
#
proof_path <- NULL
proof_path <- safe_mkdir( dir = file.path( opt$out_path, "proof") )

va60_proof_csv <- file.path( proof_path, "va60_proof.csv" )
va60_proof_tib <- NULL
va60_proof_tib <- snv_full_tab %>% 
  dplyr::filter( Probe_Source == "VA_60" ) %>% 
  dplyr::select( Probe_ID, Probe_Source, v12_gts20_GM12878 )
readr::write_csv( x = va60_proof_tib, file = va60_proof_csv )

v2_889_proof_csv <- file.path( proof_path, "v2_889_proof.csv" )
v2_889_proof_tib <- NULL
v2_889_proof_tib <- snv_full_tab %>% dplyr::filter( Probe_Source == "V2_889" ) %>% 
  dplyr::select( Probe_ID, Probe_Source, v2_gts20_GM12878)
readr::write_csv( x = v2_889_proof_tib, file = v2_889_proof_csv )

v2_889_proof_sum <- NULL
v2_889_proof_sum <- v2_889_proof_tib %>% 
  dplyr::group_by( v2_gts20_GM12878 ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )





snv_full_tab %>% dplyr::filter( Probe_Source == "V2_889" ) %>% 
  dplyr::select( Probe_ID, Probe_Source, v2_gts20_HELA) %>% 
  dplyr::group_by( v2_gts20_HELA ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

snv_full_tab %>% dplyr::filter( Probe_Source == "V2_889" ) %>% 
  dplyr::select( Probe_ID, Probe_Source, v2_gts20_GM12878) %>% 
  dplyr::group_by( v2_gts20_GM12878 ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )






# snvs_pers_tab <- NULL
# snvs_pers_tab <- snvs_pers %>% 
#   # dplyr::filter( Exp_Set == "" ) %>%
#   dplyr::mutate( Mat_Per = round(Match_Cnt/Total_Cnt,3), 
#                  Exp_Str = paste(Gts_Min,Exp_Set,Sam_Set, sep="_") ) %>% 
#   dplyr::select( Probe_ID,Exp_Str,Mat_Per ) %>%
#   tidyr::pivot_wider( id_cols = c(Probe_ID), 
#                       names_from = c(Exp_Str), 
#                       values_from = Mat_Per, 
#                       values_fill = 0 )
# 
# snvs_pers_tab <- NULL
# snvs_pers_tab <- snvs_pers %>% 
#   dplyr::mutate( Mat_Per = round(Match_Cnt/Total_Cnt,3), 
#                  Exp_Str = paste(Exp_Set,Gts_Min,Sam_Set, sep="_") ) %>% 
#   dplyr::select( Probe_ID,Exp_Str,Mat_Per ) %>%
#   split( .$Exp_Str ) %>% 
#   lapply( function( x ) {
#     x %>% tidyr::pivot_wider( 
#       id_cols = c(Probe_ID), 
#       names_from = c(Exp_Str), 
#       values_from = Mat_Per, 
#       values_fill = 0 ) 
#   })

snvs_pers_tab %>% dplyr::bind_rows( .id = "" )
  

hela_join_g0_tib <- NULL
hela_join_g0_tib <- dplyr::inner_join( 
  snvs_pers_tab$CellLine_Pass_v1_gts0_HELA, 
  snvs_pers_tab$CellLine_Pass_v2_gts0_HELA, by=c("Probe_ID") ) %>% 
  dplyr::inner_join( snvs_pers_tab$CellLine_Pass_v3_gts0_HELA, by=c("Probe_ID") )

hela_join_g0_tib$CellLine_Pass_v1_gts0_HELA %>% mean( na.rm = TRUE )
hela_join_g0_tib$CellLine_Pass_v2_gts0_HELA %>% mean( na.rm = TRUE )
hela_join_g0_tib$CellLine_Pass_v3_gts0_HELA %>% mean( na.rm = TRUE )

hela_join_g20_tib <- NULL
hela_join_g20_tib <- dplyr::inner_join( 
  snvs_pers_tab$CellLine_Pass_v1_gts20_HELA, 
  snvs_pers_tab$CellLine_Pass_v2_gts20_HELA, by=c("Probe_ID") ) %>% 
  dplyr::inner_join( snvs_pers_tab$CellLine_Pass_v3_gts20_HELA, by=c("Probe_ID") )

hela_join_g20_tib$CellLine_Pass_v1_gts20_HELA %>% mean( na.rm = TRUE )
hela_join_g20_tib$CellLine_Pass_v2_gts20_HELA %>% mean( na.rm = TRUE )
hela_join_g20_tib$CellLine_Pass_v3_gts20_HELA %>% mean( na.rm = TRUE )


if ( FALSE ) {
  src_ssh_tib <- NULL
  src_ssh_tib <- ssh_tib %>% split( .$Sample_Source )
  sam_ssh_tib <- NULL
  sam_ssh_tib <- ssh_tib %>% 
    dplyr::filter( Sample_Source == "Alpha" ) %>%
    split( .$Sample_Base )
  
  
  if ( FALSE ) {
    snvs_pers <- NULL
    for ( sub_key in names(sub_lst) ) {
      sub_vec <- NULL
      sub_vec <- sub_lst[[sub_key]]
      
      for ( sam_key in names(sam_ssh_tib) ) {
        # [TBD]: Skip Sample Count < 4
        
        if (p2) cat(glue::glue("[{pmssg}]: sub={sub_key}, sam={sam_key}...{RET}"))
        
        sam_vec <- NULL
        sam_vec <- ssh_tib %>% 
          dplyr::filter( Sample_Base == sam_key ) %>% dplyr::pull( Plot_ID )
        
        # Same Thing as Above...
        # sam_vec <- sam_ssh_tib[[sam_key]] %>% dplyr::pull( Plot_ID )
        
        sam_snvs_gtc <- NULL
        sam_snvs_gtc <- snvs_gtc[ sub_vec, sam_vec ]
        
        # [Done]: Functionalize Per Probe Performance [ rr x ii x jj ]
        per_snvs_tib <- NULL
        per_snvs_tib <- gtc_snv_to_performance( x = sam_snvs_gtc ) %>%
          dplyr::mutate( Sub_Set = sub_key, Sample_Set = sam_key )
        
        snvs_pers <- dplyr::bind_rows( snvs_pers,per_snvs_tib )
        
        # [TBD]: Add gts_min, call_rate_min to snvs_pers
        # [TBD]: Create new classes in ssh_tib to define experiments
        
        # By Source
        if ( FALSE ) {
          for ( src_key in names(src_ssh_tib) ) {
            src_vec <- NULL
            src_vec <- ssh_tib %>% 
              dplyr::filter( Sample_Base == sam_key ) %>%
              dplyr::filter( Sample_Source == src_key ) %>% dplyr::pull( Plot_ID )
            
            if ( length(src_vec) == 0 )
              
              src_snvs_gtc <- NULL
            src_snvs_gtc <- snvs_gtc[ sub_vec, src_vec ]
            
            break
          }
        }
        
        # break
      } # for ( sam_key in names(sam_ssh_tib) )
      
      
      # break
    } # for ( sub_key in names(sub_lst) )
    
    # break
  } # for ( gts_min in gtc_vec )
  
} # TRUE

snvs_pers %>% 
  dplyr::mutate( Mat_Per = round(Match_Cnt/Total_Cnt,3), 
                 Exp_Str = paste(Sub_Set,Sample_Set, sep="_") ) %>% 
  dplyr::select( Probe_ID,Exp_Str,Mat_Per ) %>%
  tidyr::pivot_wider( id_cols = c(Probe_ID), 
                      names_from = c(Exp_Str), 
                      values_from = Mat_Per, 
                      values_fill = 0 )
































if ( FALSE ) {
  
  ssh_tib %>% dplyr::group_by( Sample_Source ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  
  # [TBD]: Loop over: sub_lst x Sample_Source
  
  gtc_dat00 <- NULL
  gtc_dat00 <- tmp_tib[ sub_lst[["V2"]], ssh_tib %>% dplyr::filter( Sample_Source == "Alpha" ) %>% dplyr::pull( Plot_ID ) ] 
  
  # gtc_dat00 <- snvs_dats[[1]]$snv_dat$GTC[ , ssh_tib %>% dplyr::filter( Sample_Source == "Alpha" ) %>% dplyr::pull( Plot_ID ) ] 
  # 
  # gtc_dat00 <- snvs_dats[[1]]$snv_dat$GTC
  # gtc_dat00[ which( is.na(gtc_dat00) ) ] <- 6

  tot_cnt <- gtc_dat00 %>% base::nrow()
  out_cnt <- which( gtc_dat00[ , 1] >  5 | gtc_dat00[ , 2] >  5 ) %>% length()
  nan_cnt <- which( gtc_dat00[ , 1] == 5 | gtc_dat00[ , 2] == 5 ) %>% length()
  mat_cnt <- which( gtc_dat00[ , 1] <  5 & gtc_dat00[ , 2] <  5 & matrixStats::rowDiffs( gtc_dat00[ , c(1,2) ] ) == 0 ) %>% length()
  mis_cnt <- which( gtc_dat00[ , 1] <  5 & gtc_dat00[ , 2] <  5 & matrixStats::rowDiffs( gtc_dat00[ , c(1,2) ] ) != 0 ) %>% length()
  sum_cnt <- out_cnt + nan_cnt + mat_cnt + mis_cnt
  
  cat(glue::glue("{pmssg} out_cnt={out_cnt}.{RET}"))
  cat(glue::glue("{pmssg} nan_cnt={nan_cnt}.{RET}"))
  cat(glue::glue("{pmssg} mat_cnt={mat_cnt}.{RET}"))
  cat(glue::glue("{pmssg} mis_cnt={mis_cnt}.{RET}"))
  cat(glue::glue("{pmssg} tot_cnt={tot_cnt}.{RET}"))
  cat(glue::glue("{pmssg} sum_cnt={sum_cnt}.{RET}"))
  
  snv_cnt_tib <- NULL
  for ( rr in c(1:base::nrow(gtc_dat00) ) ) {
    tcnt <- 0
    mcnt <- 0
    for ( ii in c(1:base::ncol(gtc_dat00) ) ) {
      for ( jj in c(1:base::ncol(gtc_dat00) ) ) {
        if ( ii >= jj ) next
        tcnt = tcnt + 1
        if ( gtc_dat00[ rr,ii] < 5 && gtc_dat00[ rr,jj ] < 6 &&
             gtc_dat00[ rr,ii ] == gtc_dat00[ rr,jj] ) mcnt <- mcnt + 1
      }
    }
    snv_cnt_tib <- snv_cnt_tib %>% 
      dplyr::bind_rows(
        tibble::tibble(
          Probe_ID = rownames(gtc_dat00)[rr],
          Total_Cnt = tcnt,
          Match_Cnt = mcnt
        )
      )
  }

}

if ( FALSE ) {
  
  # Last Column has no name...
  # snvs_dats$snv_dat$GTC[ , 1:136 ]
  
  snv_dats <- list()
  for ( sub_key in names(sub_lst) ) {
    sub_vec <- NULL
    sub_vec <- sub_lst[[sub_key]]
    # if ( prb_key == "MVP" ) sub_vec <- va_ids_tib$Probe_ID
    
    for ( gts_min in gtc_vec ) {
      gts_key <- paste0("gts",gts_min)
      
      snv_key <- NULL
      snv_key <- paste( sub_key,gts_key, sep="_" )
      if ( p1 ) cat(glue::glue("{pmssg} Fingerprint Params: {snv_key}...{RET}"))
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #.                          Evaluate:: by Product/Sample
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      cur_out_path <- safe_mkdir( file.path( opt$out_path, "partitions", "byProductSample" ) )
      for ( grp_key in prd_vec ) {
        
        for ( sam_key in sam_vec ) {
          cur_key <- paste( snv_key,grp_key,sam_key, sep="_" )
          
          cur_ssh_tib <- NULL
          cur_ssh_tib <- ssh_tib %>% 
            dplyr::filter( detect_version == grp_key & Sample_Base == sam_key)
          
          if ( base::nrow( cur_ssh_tib ) < 2 ) next
          
          snv_dats[[cur_key]] <- NULL
          snv_dats[[cur_key]] <- analyze_snvs( 
            vcfs = vcf_list[cur_ssh_tib$Sentrix_Name], 
            run_tag = cur_key, 
            ssh_tib = cur_ssh_tib, 
            sub_vec = sub_vec, 
            gts_min = gts_min, 
            fval = NA_real_, 
            outDir = cur_out_path, 
            write_out = TRUE, plot_heat = TRUE )
        }
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #.                          Evaluate:: by Class
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      cur_out_path <- safe_mkdir( file.path( opt$out_path, "partitions", "byClass" ) )
      for ( grp_key in grp_vec ) {
        cur_key <- paste( snv_key,grp_key, sep="_" )
        
        cur_ssh_tib <- NULL
        cur_ssh_tib <- ssh_tib %>% dplyr::filter( Sample_Class == grp_key )
        
        if ( base::nrow( cur_ssh_tib ) < 2 ) next
        
        snv_dats[[cur_key]] <- NULL
        snv_dats[[cur_key]] <- analyze_snvs( 
          vcfs = vcf_list[cur_ssh_tib$Sentrix_Name], 
          run_tag = cur_key, 
          ssh_tib = cur_ssh_tib, 
          sub_vec = sub_vec, 
          gts_min = gts_min, 
          fval = NA_real_, 
          outDir = cur_out_path, 
          write_out = TRUE, plot_heat = TRUE )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #.                          Evaluate:: by Sample
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      cur_out_path <- safe_mkdir( file.path( opt$out_path, "partitions", "bySample" ) )
      for ( grp_key in sam_vec ) {
        cur_key <- paste( snv_key,grp_key, sep="_" )
        
        cur_ssh_tib <- NULL
        cur_ssh_tib <- ssh_tib %>% dplyr::filter( Sample_Base == grp_key )
        
        if ( base::nrow( cur_ssh_tib ) < 2 ) next
        if ( base::nrow( cur_ssh_tib ) < 5 ) next
        
        snv_dats[[cur_key]] <- NULL
        snv_dats[[cur_key]] <- analyze_snvs( 
          vcfs = vcf_list[cur_ssh_tib$Sentrix_Name], 
          run_tag = cur_key, 
          ssh_tib = cur_ssh_tib, 
          sub_vec = sub_vec, 
          gts_min = gts_min, 
          fval = NA_real_, 
          outDir = cur_out_path, 
          write_out = TRUE, plot_heat = TRUE )
      }
      
      if ( opt$single ) break
    }
    if ( opt$single ) break
  }
}

#
# tmp_dats <- analyze_snvs( vcfs = vcf_list, run_tag = "All", ssh_tib = cur_ssh_tib, sub_vec = sub_vec, gts_min = gts_min, fval = 5.0, outDir = opt$out_path, vb=vb,vt=vt+1,tc=tc, tt=NULL )
#

sig_cols <- NULL
sig_cols <- readr::cols(
  Sample_Name = readr::col_character(),
  Call_Rate   = readr::col_double(),
  Fingerprint = readr::col_character()
)

sig_vcf_path <- file.path( opt$top_path, "scratch/stable_analysis_EPICv2_docker_snps/EPICv2-UCSC-v5/fingerprint_vcfs" )
sig_vcf_list <- NULL
sig_vcf_list <- file_list( path    = sig_vcf_path, 
                       prefix  = sig_vcf_path, 
                       suffix  = ".fingerprint-signatures.csv", 
                       pattern = ".fingerprint-signatures.csv$",
                       recursive = FALSE )

sig_dat_list <- NULL
sig_dat_list <- sig_vcf_list %>%
  lapply( readr::read_csv, 
          col_names = names(sig_cols$cols), 
          col_types = sig_cols )

v2_tmp_tib <- NULL
v2_tmp_tib <- sig_dat_list$`EPICv2.gts_min-20`

# RS Names: v2_ids_vec

v2_sep_vec <- NULL
v2_sep_vec <- c( 1:(length(v2_ids_vec)-1) )
v2_snp_vec <- NULL
v2_snp_vec <- paste0( "v",c( 1:length(v2_ids_vec) ) )

v2_sig_mat <- NULL
v2_sig_mat <- v2_tmp_tib %>% 
  tidyr::separate( Fingerprint, into = v2_ids_vec, 
                   sep = v2_sep_vec, remove = TRUE,
                   convert = TRUE )


if ( FALSE ) {
  
  sel_ssh_tib <- NULL
  sel_ssh_tib <- ssh_tib %>% 
    dplyr::filter( Sample_Source == "Alpha" | Sample_Source == "GSE228820" ) %>%
    dplyr::select( Plot_ID,Sample_Base,detect_version,cg_calls_pass_perc_1,
                   Sample_Class, Concentration, Sample_Group ) %>%
    dplyr::rename( Product = detect_version, Call_Rate = cg_calls_pass_perc_1 ) %>%
    dplyr::filter( Sample_Base != "PROMEGAMOUSE" ) %>%
    dplyr::arrange( Call_Rate )

  sel_ssh_tib %>% 
    ggplot2::ggplot( aes(x=Call_Rate, fill=Sample_Base) ) + 
    ggplot2::geom_density( alpha=0.2 ) +
    ggplot2::facet_grid( rows = vars(Product) )
  
  # sel_ssh_tib %>% 
  #   ggplot2::ggplot( aes(x=Call_Rate) ) + 
  #   ggplot2::geom_density( alpha=0.2 ) +
  #   ggplot2::facet_grid( rows = vars(Product) )
  # 
  
}

# Split by Sample
# Produce Probe Level Genotype Agreement Rate: GAF (Full), GAR (RAF), GAA (MAF)
#

sel_ssh_tib <- NULL
sel_ssh_tib <- ssh_tib %>% 
  dplyr::filter( Sample_Source == "Alpha" | Sample_Source == "GSE228820" ) %>%
  dplyr::select( Plot_ID,Sample_Base,detect_version,cg_calls_pass_perc_1,
                 Sample_Class, Concentration, Sample_Group ) %>%
  dplyr::rename( Product = detect_version, Call_Rate = cg_calls_pass_perc_1 ) %>%
  dplyr::filter( Sample_Base != "PROMEGAMOUSE" ) %>%
  dplyr::arrange( Call_Rate ) %>% 
  dplyr::filter( Call_Rate >= 75 )

pcr_all_tab <- NULL

all_pcr_tab <- NULL
all_tot_tab <- NULL
for ( exp_key in names(sig_dat_list) ) {

  cur_fig_tib <- NULL
  cur_fig_tib <- sig_dat_list[[exp_key]] %>% 
    tidyr::separate( Fingerprint, into = v2_ids_vec, 
                     sep = v2_sep_vec, remove = TRUE,
                     convert = TRUE )
  
  cur_fig_lst <- NULL
  cur_fig_lst <- sel_ssh_tib %>%
    dplyr::inner_join( cur_fig_tib,
                       by=c("Sample_Base"="Sample_Name","Call_Rate") ) %>%
    dplyr::rename( Sample_Name = Sample_Base ) %>%
    split(.$Sample_Name)
  
  for ( sam_key in names(cur_fig_lst) ) {

    cur_mat <- NULL
    cur_mat <- cur_fig_lst[[sam_key]] %>% 
      dplyr::select( -Sample_Name, -Product, -Call_Rate, 
                     -Sample_Class, -Concentration, -Sample_Group ) %>% 
      tibble::column_to_rownames( var = "Plot_ID" ) %>% 
      as.matrix() %>% 
      as.data.frame() %>% as.matrix() %>% t()
    
    snv_tot_cnt <- base::ncol(cur_mat)
    
    hit_vec <- c()
    tot_vec <- c()
    for ( ii in c(1:base::nrow(cur_mat) ) ) {
      snv_hit_cnt <- 0
      snv_tot_cnt <- 0
      for ( jj in c(1:base::ncol(cur_mat) ) ) {
        for ( kk in c(1:base::ncol(cur_mat) ) ) {
          if ( jj >= kk ) next
          if ( !is.na(cur_mat[ii,jj]) && !is.na(cur_mat[ii,kk]) &&
               cur_mat[ii,jj] != 5 && cur_mat[ii,kk] != 5 &&
               cur_mat[ii,jj] == cur_mat[ii,kk] ) {
            snv_hit_cnt <- snv_hit_cnt + 1
          }
          snv_tot_cnt <- snv_tot_cnt + 1
        }
      }
      hit_vec <- c( hit_vec, snv_hit_cnt )
      tot_vec <- c( tot_vec, snv_tot_cnt )
      # if ( ii > 5 ) break
    }
    
    all_key <- NULL
    all_key <- paste( exp_key,sam_key, sep="_" )
    cur_pcr_tib <- NULL
    cur_pcr_tib <- tibble::tibble(
      Probe_ID   = rownames(cur_mat),
      !!all_key := hit_vec
    )
    cur_tot_tib <- NULL
    cur_tot_tib <- tibble::tibble(
      Probe_ID   = rownames(cur_mat),
      !!all_key := tot_vec
    )
    # Should divide by denominator...
    
    if ( !is.null(all_pcr_tab) ) {
      all_pcr_tab <- all_pcr_tab %>% 
        dplyr::full_join( cur_pcr_tib, by=c("Probe_ID") )
    } else {
      all_pcr_tab <- cur_pcr_tib
    }
    
    if ( !is.null(all_tot_tab) ) {
      all_tot_tab <- all_tot_tab %>% 
        dplyr::full_join( cur_tot_tib, by=c("Probe_ID") )
    } else {
      all_tot_tab <- cur_tot_tib
    }
    
    pcr_all_tab <- pcr_all_tab %>%
      dplyr::bind_rows(
        tibble::tibble(
          Exp_Key = exp_key,
          Sample_Key = sam_key,
          Probe_ID = rownames(cur_mat),
          Tot_Cnt = tot_vec,
          Mat_Cnt = hit_vec,
        )
      )
    
    # break
  }

  # dplyr::select( -Sample_Name, -Call_Rate ) %>% tidyr::pivot_longer( names_to = c(Probe_ID), values_to = c(Mat_Cnt) )
  # break
}

# As Matrix::
# all_pcr_tab / all_tot_tab

all_pcr_mat <- NULL
all_pcr_mat <- all_pcr_tab %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()
all_tot_mat <- NULL
all_tot_mat <- all_tot_tab %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()

all_per_mat <- NULL
all_per_mat <- all_pcr_mat / all_tot_mat

all_per_tib <- NULL
all_per_tib <- all_per_mat %>% 
  as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    Product = dplyr::case_when(
      Probe_ID %in% sub_lst[["VA"]] ~ "VA",
      Probe_ID %in% sub_lst[["V1"]] & !Probe_ID %in% sub_lst[["V2"]] ~ "v1_only",
      Probe_ID %in% sub_lst[["V2"]] & !Probe_ID %in% sub_lst[["V1"]] ~ "v2_only",
      Probe_ID %in% sub_lst[["V1"]] &  Probe_ID %in% sub_lst[["V2"]] ~ "v3_both",
      TRUE ~ NA_character_
    )
  ) %>% dplyr::select( Probe_ID,Product, dplyr::everything() )

# Fix Names to remove non-favorable characters...  
colnames(all_per_tib) <- all_per_tib %>% names() %>% 
  stringr::str_replace_all("-","_") %>% stringr::str_replace_all("\\.","_")

if ( FALSE ) {

  all_per_tib %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer( cols = c(Probe_ID,Product), 
                         names_to = "Experiment", 
                         values_to = "Match_Percent" )
  
  all_per_tib %>% 
    dplyr::select( -Probe_ID ) %>% 
    tidyr::pivot_longer( cols = c(Product), 
                         names_to = c(), 
                         values_to = "Values" )
  
}

# tidyr::pivot_longer( cols = c(!Probe_ID,!Product), names_to = "Metric", values_to = "Values" )
# all_per_tib %>% tidyr::pivot_longer()
# all_per_tib %>% dplyr::select( -Probe_ID ) %>% tidyr::gather( key = "Product" )


plt_tib <- NULL
plt_tib <- all_per_mat %>% 
  matrixStats::rowMedians( na.rm = TRUE, useNames = TRUE ) %>% as.data.frame() %>% 
  tibble::rownames_to_column( var = "Probe_ID" ) %>% 
  tibble::as_tibble() %>% magrittr::set_names(c("Probe_ID","Med")) %>% 
  dplyr::mutate(
    Product = dplyr::case_when(
      Probe_ID %in% sub_lst[["VA"]] ~ "VA",
      Probe_ID %in% sub_lst[["V1"]] & !Probe_ID %in% sub_lst[["V2"]] ~ "v1_only",
      Probe_ID %in% sub_lst[["V2"]] & !Probe_ID %in% sub_lst[["V1"]] ~ "v2_only",
      Probe_ID %in% sub_lst[["V1"]] &  Probe_ID %in% sub_lst[["V2"]] ~ "v3_both",
      TRUE ~ NA_character_
    )
  )

plt_sum <- NULL
plt_sum <- plt_tib %>% 
  dplyr::group_by( Product ) %>% 
  dplyr::summarise( Tot_Cnt = n(),
                    Avg = mean(Med, na.rm=TRUE ),
                    .groups = "drop" )

# plt_tib %>% ggplot2::ggplot( aes(x=Med, fill=Product) ) + geom_density( alpha=0.2 )

pcr_plt_tab <- NULL
pcr_plt_tab <- pcr_all_tab %>%
  dplyr::mutate(
    Product = dplyr::case_when(
      Probe_ID %in% sub_lst[["VA"]] ~ "VA",
      Probe_ID %in% sub_lst[["V1"]] & !Probe_ID %in% sub_lst[["V2"]] ~ "v1_only",
      Probe_ID %in% sub_lst[["V2"]] & !Probe_ID %in% sub_lst[["V1"]] ~ "v2_only",
      Probe_ID %in% sub_lst[["V1"]] &  Probe_ID %in% sub_lst[["V2"]] ~ "v3_both",
      TRUE ~ NA_character_
    ),
    Per = Mat_Cnt/Tot_Cnt
  ) %>% dplyr::select( Exp_Key, Product, Sample_Key, dplyr::everything() )

pcr_plt_ggg <- NULL
pcr_plt_ggg <- pcr_plt_tab %>% 
  ggplot2::ggplot( aes(x=Per, fill=Sample_Key) ) + 
  ggplot2::geom_density( alpha=0.2 ) +
  ggplot2::facet_grid( rows = vars(Exp_Key),
                       cols = vars(Product) )

pcr_plt_sum <- NULL
pcr_plt_sum <- pcr_plt_tab %>% 
  dplyr::group_by( Exp_Key,Product,Sample_Key ) %>%
  dplyr::summarise( Tot_Cnt = n(),
                    Avg_Per = mean(Per),
                    Sds_Per = sd(Per),
                    Med_Per = median(Per),
                    .groups = "drop" )
if ( p2 ) pcr_plt_sum %>% print( n=base::nrow(pcr_plt_sum) )

pcr_plt_csv <- file.path( opt$out_path, "snv-consistency-table.csv" )
readr::write_csv( x = pcr_plt_tab, file = pcr_plt_csv )

pcr_plt_sum2 <- NULL
pcr_plt_sum2 <- pcr_plt_tab %>% 
  dplyr::group_by( Probe_ID ) %>% 
  dplyr::filter( !is.na(Per) ) %>%
  dplyr::summarise( Tot = n(),
                    Avg = mean(Per, na.rm=TRUE ),
                    Med = median( Per, na.rm=TRUE ),
                    Sds = sd( Per, na.rm=TRUE ),
                    .groups = "drop" ) %>% 
  dplyr::arrange( -Avg )

pcr_plt_sum2 %>% 
  ggplot2::ggplot( aes(x=Avg) ) +
  ggplot2::geom_density( alpha=0.2 )

#
# This Plot is legit!!!
#
pcr_plt_ggg2 <- NULL
pcr_plt_ggg2 <- pcr_plt_tab %>% 
  dplyr::filter( Per != 0 ) %>%
  ggplot2::ggplot( aes(x=Per, fill=Sample_Key) ) +
  ggplot2::geom_density( alpha=0.2 ) +
  ggplot2::facet_grid( rows = vars(Product),
                       cols = vars(Exp_Key) )

pcr_plt_ggg3 <- NULL
pcr_plt_ggg3 <- pcr_plt_tab %>% 
  dplyr::filter( Per > 0.75 ) %>%
  ggplot2::ggplot( aes(x=Per, fill=Sample_Key) ) +
  ggplot2::geom_density( alpha=0.2 ) +
  ggplot2::facet_grid( rows = vars(Product),
                       cols = vars(Exp_Key) )

pcr_plt_tab %>% 
  dplyr::filter( Per > 0.75 ) %>%
  dplyr::group_by( Exp_Key,Product,Sample_Key ) %>%
  dplyr::summarise( Tot = n(),
                    Avg = mean( Per ),
                    Med = median( Per ),
                    Sds = sd( Per ),
                    .groups = "drop" ) %>% 
  dplyr::arrange( -Avg ) %>% print( n=1000 )

pcr_plt_tab %>% 
  dplyr::filter( Per > 0.75 ) %>%
  dplyr::filter( Exp_Key == "EPICv1-MVP.gts_min-20" ) %>% 
  dplyr::filter( Product == "v2_only" ) %>% 
  dplyr::filter( Sample_Key == "HELA" ) %>% 
  dplyr::summarise( mean( Per ) )















if ( FALSE ) {
  if ( FALSE ) {
    
    opt$pdf_path = safe_mkdir( dir = file.path( opt$out_path,"pdf") )
    opt$csv_path = safe_mkdir( dir = file.path( opt$out_path,"csv") )
    opt$rds_path = safe_mkdir( dir = file.path( opt$out_path,"rds") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             SNP Analysis::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    snv_rds <- file.path( opt$rds_path, paste0( snv_key,".rds") )
    hit_csv <- file.path( opt$csv_path, paste0( snv_key,".hit.csv.gz") )
    mis_csv <- file.path( opt$csv_path, paste0( snv_key,".mis.csv.gz") )
    nan_csv <- file.path( opt$csv_path, paste0( snv_key,".nan.csv.gz") )
    hit_pdf <- file.path( opt$pdf_path, paste0( snv_key,".hit.heatmap.pdf") )
    
    snv_dat <- NULL
    snv_dat <- vcf_list %>% # head(n=3) %>% 
      lapply( read_snv_vcf ) %>% 
      lapply( parse_snv_vcf, gts_min = gts_min, fval = NA_real_ ) %>%
      stack_snv_list( ssh = ssh_tib, sub = sub_vec ) %>%
      split_snv_stack( nstr = "Sample_Name", istr = "Target_ID" ) 
    
    snv_tab <- NULL
    snv_tab <- snv_stack_to_contingency_tab( x = snv_dat$GTC )
    
    snv_dats[[snv_key]] = NULL
    snv_dats[[snv_key]] = snv_dat
    
    #
    # Write/Plot Data::
    #
    if ( opt$write_out ) {
      snv_hit_tib <- NULL
      snv_hit_tib <- snv_tab$hit_mat %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
      snv_mis_tib <- NULL
      snv_mis_tib <- snv_tab$mis_mat %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
      snv_nan_tib <- NULL
      snv_nan_tib <- snv_tab$nan_mat %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
      
      readr::write_rds( x = snv_tab, file = snv_rds, compress = "gz" )
      readr::write_csv( x = snv_hit_tib, file = hit_csv )
      readr::write_csv( x = snv_mis_tib, file = mis_csv )
      readr::write_csv( x = snv_nan_tib, file = nan_csv )
    }
    
    if ( opt$plot_heat ) {
      pdf( file = hit_pdf, width = 10, height = 10 )
      heatmap( snv_tab$hit_mat )
      dev.off()
    }
    # tmp_dat <- snv_tab$hit_mat / ( snv_tab$hit_mat + snv_tab$mis_mat + snv_tab$nan_mat )
    
    prb_tibs[[snv_key]] <- NULL
    prb_tibs[[snv_key]] <- tibble::tibble(
      Target_ID = rownames(snv_dat$GTC),
      Probe_ID  = Target_ID %>% stringr::str_remove(":.*$"),
      Loci_ID   = Probe_ID %>% stringr::str_remove("_.*$"),
      Avg_Mat   = snv_dat$GTC %>% matrixStats::rowMeans2( na.rm = TRUE ),
      Med_Mat   = snv_dat$GTC %>% matrixStats::rowMedians( na.rm = TRUE ),
      Sds_Mat   = snv_dat$GTC %>% matrixStats::rowSds( na.rm = TRUE ),
      Mad_Mat   = snv_dat$GTC %>% matrixStats::rowMads( na.rm = TRUE ),
    )
    # gtc_prb_tib %>% dplyr::arrange( -Med_Mat )
    # gtc_prb_tib %>% dplyr::distinct( Loci_ID )
    # gtc_prb_tib %>% dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>% dplyr::arrange( -Med_Mat )
    
    # break
  }
  # break
}












snp_dats <- list()
opt$gts_mins <- c( 0, 20 )
opt$sub_strs <- c( "All", "MVP")

for ( sub_str in opt$sub_strs ) {
  sub_vec <- NULL
  if ( sub_str == "MVP" ) sub_vec <- va_ids_tib$Probe_ID
  sam_vec <- NULL
  sam_vec = unique(ssh_tib$Sample_Base)
  
  for ( gts_min in opt$gts_mins ) {
    gts_str <- paste0("gts",gts_min)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             SNP Analysis::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    snp_key <- paste0( gts_str,"_",sub_str  )
    if ( p1 ) cat(glue::glue("{pmssg} Fingerprint Params: {snp_key}...{RET}"))
    
    snp_dat <- NULL
    snp_dat <- fingerprint_vcfs( vcfs = vcf_list, 
                                 sam_tib = ssh_tib,
                                 sam_vec = sam_vec,
                                 sub_vec = sub_vec,
                                 run_sig = TRUE,
                                 gts_min = gts_min,
                                 out_dir = file.path( opt$out_path ), 
                                 run_tag = snp_key,
                                 reload = opt$reload, 
                                 reload_min = 10, 
                                 ret_data  = TRUE,
                                 write_out = FALSE, 
                                 write_sum = TRUE, 
                                 write_sig = TRUE,
                                 
                                 vb=vb+10,vt=vt+1,tc=tc+1, tt=tt )
    snp_dats[[snp_key]] = snp_dat
    
    
    break 
  }
  break
}
















ssh_tib %>% 
  dplyr::select( Sentrix_Name, Plot_ID ) %>% 
  dplyr::right_join( tmp_dat, by=c("Sentrix_Name"), multiple = "all" ) %>%
  dplyr::select( -Sentrix_Name ) %>%
  dplyr::rename( Sample_Name = Plot_ID )


tmp_dat <- NULL
tmp_dat <- vcf_list %>% head(n=3) %>% 
  lapply( readr::read_tsv, 
          col_names=names(vcf_cols$cols), 
          col_types=vcf_cols, skip=6 ) %>% 
  lapply( parse_snv_vcf ) %>%
  dplyr::bind_rows( .id = "Sentrix_Name" ) %>%
  dplyr::select( -GS ) %>%
  
  tidyr::pivot_longer( cols = c( Qual,Filter,PVF,GTS,RGT,AGT,GTC ), 
                       names_to  = c("Key"), 
                       values_to = c("Val") ) %>% split(.$Key) %>% 
  
  lapply( function(x) { 
    x %>% dplyr::select(-Key) %>% 
      tidyr::pivot_wider( id_cols = c(Target_ID), 
                          names_from = c(Sentrix_Name), 
                          values_from = c(Val) )
  }) %>%
  lapply( function(x) {
    x %>% tibble::column_to_rownames( var = "Target_ID") %>% 
      as.data.frame() %>% as.matrix()
  })

# Apply Filtering...

# Set Column Names



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                            Previous VCF Lists::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


v1_vcf_list <- NULL
v1_vcf_list <- file_list( path    = v1_dat_path, 
                          prefix  = v1_dat_path,
                          suffix  = "_EPIC_B4_raw.snps.vcf", 
                          pattern = "_EPIC_B4_raw.snps.vcf$",
                          recursive = TRUE )

v1_vcf_list <- NULL
v1_vcf_list <- file_list( path    = v1_dat_path, 
                          prefix  = v1_dat_path,
                          suffix  = "_EPIC_B4_raw.snps.vcf", 
                          pattern = "_EPIC_B4_raw.snps.vcf$",
                          recursive = TRUE )

v2_vcf_list <- NULL
v2_vcf_list <- file_list( path    = v2_dat_path, 
                          prefix  = v2_dat_path,
                          suffix  = "_EPIC_A1_raw.snps.vcf", 
                          pattern = "_EPIC_A1_raw.snps.vcf$",
                          recursive = TRUE )

wd_vcf_list <- NULL
wd_vcf_list <- file_list( path    = wd_dat_path, 
                          prefix  = wd_dat_path,
                          suffix  = "_EPIC_A1_raw.snps.vcf", 
                          pattern = "_EPIC_A1_raw.snps.vcf$",
                          recursive = TRUE )


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Load Sample Sheets:: Wanding
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  v1_dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/EPICv1" )
  v2_dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/EPICv2" )
  wd_dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/EPICv2_GSE228820" )
  
  true_wd_ssh_tib <- NULL
  true_wd_ssh_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/sampleSheets/EPICv2_GSE228820_Wanding.base.sampleSheet.v4.csv" )
  true_wd_ssh_tib <- readr::read_csv( file = true_wd_ssh_csv, show_col_types = FALSE ) %>%
    clean_tib()
  
  true_wd_ssh_tib %>% dplyr::group_by()
  
  true_wd_ssh_sum <- true_wd_ssh_tib %>% 
    dplyr::group_by( Sample_Base ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) true_wd_ssh_sum %>% print( n=base::nrow(true_wd_ssh_sum) )
  
  wd_ssh_tib <- NULL
  wd_ssh_tib <- file_list( path    = wd_dat_path, 
                           prefix  = wd_dat_path,
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
                   AutoSample_dB_Cnt_1 )
  
  dplyr::mutate( cg_calls_pass_perc_1 )
  
  wd_ssh_tib %>% dplyr::filter( !Sentrix_Name %in% true_wd_ssh_tib$Sentrix_Name )
  true_wd_ssh_tib %>% dplyr::filter( !Sentrix_Name %in% wd_ssh_tib$Sentrix_Name )
  
  # wd_ssh_tib$detect_manifest %>% unique()
  # wd_ssh_tib$Loci_Count_cg %>% unique()
  # wd_ssh_tib$Bead_Total_cg %>% unique()
  # wd_ssh_tib$Bead_Count_cg %>% unique()
  # wd_ssh_tib$Bead_Pool %>% unique()
  # wd_ssh_tib$detect_sample_cnt %>% unique()
  # wd_ssh_tib$Loci_Count_rs %>% unique()
  # wd_ssh_tib$detect_platform %>% unique()
  # wd_ssh_tib$detect_version %>% unique()
  # wd_ssh_tib$Bead_Pool %>% unique()
  # wd_ssh_tib$detect_match_cnt %>% unique()
  
  wd_sam_tib <- NULL
  wd_sam_tib <- dplyr::right_join(
    true_wd_ssh_tib %>% dplyr::select( Sentrix_Name, Sample_Base ) %>%
      dplyr::mutate( Sample_Base = Sample_Base %>% stringr::str_to_upper() ),
    wd_ssh_tib %>% 
      dplyr::select( Sentrix_Name, 
                     cg_calls_pass_perc_1,
                     AutoSample_R2_Key_1,AutoSample_dB_Key_1,
                     AutoSample_R2_Val_1,AutoSample_dB_Val_1, 
                     AutoSample_dB_Cnt_1 ),
    by=c("Sentrix_Name")
  ) %>%
    dplyr::group_by( Sample_Base ) %>% 
    dplyr::mutate( 
      Rep_Num = dplyr::row_number(),
      Sam_Suf = Sample_Base,
      PPP_Int = cg_calls_pass_perc_1 %>% as.integer(),
      Plot_ID = paste0( "v2_",Sam_Suf,Rep_Num,"_",PPP_Int)
    ) %>% dplyr::ungroup() %>%
    dplyr::select( Sentrix_Name,Plot_ID, dplyr::everything() )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Load Sample Sheets:: EPICv1/v2
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  true_ep_ssh_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/formatted/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz" )
  true_ep_ssh_tib <- NULL
  true_ep_ssh_tib <- readr::read_csv( file = true_ep_ssh_csv, show_col_types = FALSE )
  
  v1_ssh_list <- NULL
  v1_ssh_list <- file_list( path    = v1_dat_path, 
                            prefix  = v1_dat_path,
                            suffix  = "_EPIC_B4_AutoSampleSheet.csv.gz", 
                            pattern = "_EPIC_B4_AutoSampleSheet.csv.gz$",
                            recursive = TRUE )
  
  v1_ssh_tib <- NULL
  v1_ssh_tib <- v1_ssh_list %>% 
    lapply( readr::read_csv, show_col_types = FALSE ) %>% 
    dplyr::bind_rows() %>%
    dplyr::mutate( cg_calls_pass_perc_1 )
  
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
  
  
  
  
  
  
  
  
  v1_sam_tib <- NULL
  v1_sam_tib <- dplyr::right_join(
    true_ep_ssh_tib %>% dplyr::select( Sentrix_Name, Sample_Base ) %>%
      dplyr::mutate( 
        
        Sample_Base = Sample_Base %>% stringr::str_to_upper() ),
    v1_ssh_tib %>% 
      dplyr::select( Sentrix_Name, 
                     cg_calls_pass_perc_1,
                     AutoSample_R2_Key_1,AutoSample_dB_Key_1,
                     AutoSample_R2_Val_1,AutoSample_dB_Val_1, 
                     AutoSample_dB_Cnt_1 ),
    by=c("Sentrix_Name")
  )
  
  v1_sam_sum <- NULL
  v1_sam_sum <- v1_sam_tib %>% 
    dplyr::group_by( Sample_Base ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) v1_sam_sum %>% print( n=base::nrow(v1_sam_sum) )
  
  
  
  v2_sam_tib <- NULL
  v2_sam_tib <- dplyr::right_join(
    true_ep_ssh_tib %>% dplyr::select( Sentrix_Name, Sample_Base ) %>%
      dplyr::mutate( Sample_Base = Sample_Base %>% stringr::str_to_upper() ),
    v2_ssh_tib %>% 
      dplyr::select( Sentrix_Name, 
                     cg_calls_pass_perc_1,
                     AutoSample_R2_Key_1,AutoSample_dB_Key_1,
                     AutoSample_R2_Val_1,AutoSample_dB_Val_1, 
                     AutoSample_dB_Cnt_1 ),
    by=c("Sentrix_Name")
  )
  
  v2_sam_sum <- NULL
  v2_sam_sum <- v2_sam_tib %>% 
    dplyr::group_by( Sample_Base ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) v2_sam_sum %>% print( n=base::nrow(v2_sam_sum) )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#  .                         Merge Sample Sheetss::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Define Target Sample Names::
#
plot_sample_vec <- NULL
# plot_sample_vec <- v2_sam_tib$Sample_Base %>% unique() %>% as.vector()
# plot_sample_vec <- c("NA12873","HELA","JURKAT","MCF7","RAJI")
plot_sample_vec <- c("NA12873","HELA","JURKAT","RAJI")

# v2_sam_tib %>% dplyr::filter( Sample_Base == AutoSample_R2_Key_1 )

# dplyr::bind_rows( v1_sam_tib, v2_sam_tib, wd_sam_tib ) %>% 

# LEFT OFF HERE:
#
# Load Truth Sample Sheets [ Sentrix_Name, Sample_Name ]
#...
# Combine Sample Sheets
#. - [ Sentrix_Name, Plot_Name, Sample_Name, Sample_Rep, Call_Rate, ... ]
# Load all VCFs
# Load all Calls
#


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#.                         Compare Overall MAFs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

va_ids_tib <- NULL
va_ids_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/dat/VA_SNP_Selection/methGenotypeProbes.csv" )
va_ids_tib <- readr::read_csv( file = va_ids_csv, show_col_types = FALSE ) %>%
  dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>%
  tidyr::separate( ID, into=c("Chr","Pos","Alt","Ref"), 
                   sep=":", remove = TRUE, convert = TRUE ) %>%
  dplyr::mutate( Source = "VA",
                 Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"),
                 Probe_Type = Probe_ID %>% stringr::str_sub(1,2)
  ) %>% dplyr::select( Source, Probe_ID, Loci_ID, Probe_Type, dplyr::everything() )

#
# [Done]: Re-run: stable_scratch_EPICv2_docker_snps.R
#.       - Obtain new dbSNP.tabix files...
# [Done]: Write a function that can parse a tabix correctly...
# [TBD]: Distribution of REF AF: [ va,v1,v2,{v2!v1} ] x [ {rs,nv}, {cg,ch} ]
#

v1_body_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv1/rank_historic_probes/EPICv2-UCSC-v5.rank_historic_probes.csv.gz" )
v2_body_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv2/rank_historic_probes/EPICv2-UCSC-v5.rank_historic_probes.csv.gz" )

v1_dbSNP_vcf <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv1/intersect_bed_tabix/EPICv2-UCSC-v5.intersect_bed_tabix.tab.bed.gz" )
v2_dbSNP_vcf <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv2/intersect_bed_tabix/EPICv2-UCSC-v5.intersect_bed_tabix.tab.bed.gz" )

v1_body_tib <- NULL
v1_body_tib <- readr::read_csv( file = v1_body_csv, show_col_types = FALSE )
v2_body_tib <- NULL
v2_body_tib <- readr::read_csv( file = v2_body_csv, show_col_types = FALSE )

# [match_src = FALSE]: A tibble: 158,676 × 15
# [match_src = TRUE]:  A tibble: 133,167 × 44
v1_snv_tib <- NULL
v1_snv_tib <- parse_dbSNP_vcf( 
  snp_vcf = v1_dbSNP_vcf, 
  src_tib = v1_body_tib %>% dplyr::mutate( Chromosome = stringr::str_remove(Chromosome, "^chr") ),
  
  ref_vcf   = TRUE, 
  com_len   = 56, 
  rm_chr0   = TRUE,
  strip_chr = TRUE,
  snps_only = TRUE,
  parse_caf = TRUE,
  # match_src = FALSE,
  match_src = TRUE,
  min_maf   = 0.0,
  
  out_dir = file.path( opt$out_path, "EPICv1" ),
  run_tag = opt$run_name, 
  reload = opt$reload, 
  reload_min = 10,
  write_out = FALSE,
  
  vb=vb,vt=vt+1,tc=tc,tt=tt )
# v1_snv_tib %>% dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE )

# [match_src = FALSE]: A tibble: 156,286 × 15
# [match_src = TRUE]:  A tibble: 132,604 × 46
v2_snv_tib <- NULL
v2_snv_tib <- parse_dbSNP_vcf( 
  snp_vcf = v2_dbSNP_vcf, 
  src_tib = v2_body_tib %>% dplyr::mutate( Chromosome = stringr::str_remove(Chromosome, "^chr") ),
  
  ref_vcf   = TRUE, 
  com_len   = 56, 
  rm_chr0   = TRUE,
  strip_chr = TRUE,
  snps_only = TRUE,
  parse_caf = TRUE,
  # match_src = FALSE,
  match_src = TRUE,
  min_maf   = 0.0,
  
  out_dir = file.path( opt$out_path, "EPICv1" ),
  run_tag = opt$run_name, 
  reload = opt$reload, 
  reload_min = 10,
  write_out = FALSE,
  
  vb=vb,vt=vt+1,tc=tc,tt=tt )
# v2_snv_tib %>% dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE )

# snv_tib <- TRUE
# snv_tib <- dplyr::bind_rows( v1_snv_tib,v2_snv_tib ) %>%
#   dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE ) %>%
#   dplyr::arrange( Chromosome,Coordinate_SNP,SNP_ID, -MAF )

# A tibble: 175,479 × 10
snv_tib <- TRUE
snv_tib <- dplyr::bind_rows( 
  v1_snv_tib %>% dplyr::select( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, RAF,MAF, Probe_ID,Full_ID,Loci_ID),
  v2_snv_tib %>% dplyr::select( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, RAF,MAF, Probe_ID,Full_ID,Loci_ID)
) %>%
  dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE ) %>%
  dplyr::arrange( Chromosome,Coordinate_SNP,SNP_ID, -MAF )


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Plot Heatmaps:: Betas (v1/v2)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

v1_cgn_vec <- NULL
v1_cgn_vec <- v1_body_tib %>% 
  dplyr::filter( Probe_Type == "cg" ) %>% dplyr::pull( Loci_ID )

v2_calls_lst <- NULL
v2_calls_lst <- file_list( 
  path    = v2_dat_path, 
  prefix  = v2_dat_path,
  suffix  = "_EPIC_A1_ind.call.dat.csv.gz", 
  pattern = "_EPIC_A1_ind.call.dat.csv.gz$",
  recursive = TRUE )

v1_calls_lst <- NULL
v1_calls_lst <- file_list( 
  path    = v1_dat_path, 
  prefix  = v1_dat_path,
  suffix  = "_EPIC_B4_ind.call.dat.csv.gz", 
  pattern = "_EPIC_B4_ind.call.dat.csv.gz$",
  recursive = TRUE )

v2_sent_vec <- NULL
v2_sent_vec <- v2_sam_tib %>% 
  dplyr::filter( Sample_Base %in% plot_sample_vec ) %>%
  dplyr::pull( Sentrix_Name )

v1_sent_vec <- NULL
v1_sent_vec <- v1_sam_tib %>% 
  dplyr::filter( Sample_Base %in% plot_sample_vec ) %>%
  dplyr::pull( Sentrix_Name )

v2_sent_mat <- NULL
v2_sent_mat <- v2_sam_tib %>% 
  dplyr::filter( Sample_Base %in% plot_sample_vec ) %>%
  dplyr::group_by( Sample_Base ) %>% 
  dplyr::mutate( 
    Rep_Num = dplyr::row_number(),
    Sam_Suf = Sample_Base %>% stringr::str_sub(1,1),
    PPP_Int = cg_calls_pass_perc_1 %>% as.integer(),
    Plot_ID = paste0( "v2_",Sam_Suf,Rep_Num,"_",PPP_Int)
  ) %>% dplyr::ungroup() %>%
  dplyr::select( Sentrix_Name,Plot_ID ) %>%
  tibble::column_to_rownames( var = "Sentrix_Name" ) %>%
  as.matrix()

v1_sent_mat <- NULL
v1_sent_mat <- v1_sam_tib %>% 
  dplyr::filter( Sample_Base %in% plot_sample_vec ) %>%
  dplyr::group_by( Sample_Base ) %>% 
  dplyr::mutate( 
    Rep_Num = dplyr::row_number(),
    Sam_Suf = Sample_Base %>% stringr::str_sub(1,1),
    PPP_Int = cg_calls_pass_perc_1 %>% as.integer(),
    Plot_ID = paste0( "v1_",Sam_Suf,Rep_Num,"_",PPP_Int)
  ) %>% dplyr::ungroup() %>%
  dplyr::select( Sentrix_Name,Plot_ID ) %>%
  tibble::column_to_rownames( var = "Sentrix_Name" ) %>%
  as.matrix()

#
# Load EPICv2 Calls::
#
v2_data_list <- NULL
v2_data_list <- v2_calls_lst[v2_sent_vec] %>% # head() %>%
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows( .id = "Sentrix_Name" ) %>%
  tidyr::pivot_longer( cols = c( pvals_pOOBAH,pvals_PnegEcdf,betas ), 
                       names_to  = c("Key"), 
                       values_to = c("Val") ) %>% split(.$Key) %>% 
  lapply( function(x) { 
    x %>% dplyr::select(-Key) %>% 
      tidyr::pivot_wider( id_cols = c(Probe_ID), 
                          names_from = c(Sentrix_Name), 
                          values_from = c(Val) ) %>% 
      dplyr::filter( Probe_ID %in% v1_cgn_vec ) %>%
      tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()
  })

# Filter Matrix...
v2_data_list$betas[ which( v2_data_list$pvals_pOOBAH > 0.05) ] <- NA_real_
# Set Names::
colnames(v2_data_list$betas) <- v2_sent_mat[ colnames(v2_data_list$betas), ] %>% as.vector()

v2_cgn_vec <- NULL
v2_cgn_vec <- rownames(v2_data_list$betas)

#
#
# Load EPICv2 Calls::
#
#
v1_data_list <- NULL
v1_data_list <- v1_calls_lst[v1_sent_vec] %>% # head() %>%
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows( .id = "Sentrix_Name" ) %>%
  tidyr::pivot_longer( cols = c( pvals_pOOBAH,pvals_PnegEcdf,betas ), 
                       names_to  = c("Key"), 
                       values_to = c("Val") ) %>% split(.$Key) %>% 
  lapply( function(x) { 
    x %>% dplyr::select(-Key) %>% 
      tidyr::pivot_wider( id_cols = c(Probe_ID), 
                          names_from = c(Sentrix_Name), 
                          values_from = c(Val) ) %>% 
      dplyr::filter( Probe_ID %in% v2_cgn_vec ) %>%
      tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()
  })

# Filter Matrix...
v1_data_list$betas[ which( v1_data_list$pvals_pOOBAH > 0.05) ] <- NA_real_
# Set Names::
colnames(v1_data_list$betas) <- v1_sent_mat[ colnames(v1_data_list$betas), ] %>% as.vector()

# Ensure v1 is sorted before joining...
# Check dims()
cgn_ord_vec <- NULL
cgn_ord_vec <- intersect( rownames(v1_data_list$betas), rownames(v2_data_list$betas) ) %>% unique() %>% sort()

# v1_data_list$betas <- v1_data_list$betas[ cgn_ord_vec, ]
# v2_data_list$betas <- v2_data_list$betas[ cgn_ord_vec, ]
# betas <- cbind( v1_data_list$betas, v2_data_list$betas )

betas <- NULL
betas <- cbind( v1_data_list$betas[ cgn_ord_vec, ],
                v2_data_list$betas[ cgn_ord_vec, ] )
#
# Build Correlation::
#
cor_beta_mat <- NULL
cor_beta_mat <- cor( x = betas, use = "pairwise.complete.obs", method = "pearson" )

heatmap_pdf <- file.path( plot_dir, "heatmap_v1-v2.MVP.betas.pdf")
pdf( file = heatmap_pdf, width = 10, height = 10 )
heatmap( x = cor_beta_mat )
dev.off()


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Plot Heatmaps:: Betas (Wanding)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

wd_calls_lst <- NULL
wd_calls_lst <- file_list( 
  path    = wd_dat_path, 
  prefix  = wd_dat_path,
  suffix  = "_EPIC_A1_ind.call.dat.csv.gz", 
  pattern = "_EPIC_A1_ind.call.dat.csv.gz$",
  recursive = TRUE )

wd_sent_vec <- NULL
wd_sent_vec <- wd_sam_tib %>% 
  # dplyr::filter( Sample_Base %in% plot_sample_vec ) %>%
  dplyr::pull( Sentrix_Name )

wd_sent_mat <- NULL
wd_sent_mat <- wd_sam_tib %>% 
  # dplyr::filter( Sample_Base %in% plot_sample_vec ) %>%
  dplyr::group_by( Sample_Base ) %>% 
  dplyr::mutate( 
    Rep_Num = dplyr::row_number(),
    Sam_Suf = Sample_Base,
    PPP_Int = cg_calls_pass_perc_1 %>% as.integer(),
    Plot_ID = paste0( "v2_",Sam_Suf,Rep_Num,"_",PPP_Int)
  ) %>% dplyr::ungroup() %>%
  dplyr::select( Sentrix_Name,Plot_ID ) %>%
  tibble::column_to_rownames( var = "Sentrix_Name" ) %>%
  as.matrix()

#
# Load EPICv2 Calls::
#
wd_data_list <- NULL
wd_data_list <- wd_calls_lst[wd_sent_vec] %>% # head() %>%
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows( .id = "Sentrix_Name" ) %>%
  tidyr::pivot_longer( cols = c( pvals_pOOBAH,pvals_PnegEcdf,betas ), 
                       names_to  = c("Key"), 
                       values_to = c("Val") ) %>% split(.$Key) %>% 
  lapply( function(x) { 
    x %>% dplyr::select(-Key) %>% 
      tidyr::pivot_wider( id_cols = c(Probe_ID), 
                          names_from = c(Sentrix_Name), 
                          values_from = c(Val) ) %>% 
      dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
      tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()
  })

# Filter Matrix...
wd_data_list$betas[ which( wd_data_list$pvals_pOOBAH > 0.05) ] <- NA_real_
# Set Names::
colnames(wd_data_list$betas) <- wd_sent_mat[ colnames(wd_data_list$betas), ] %>% as.vector()

#
# Build Correlation::
#
wd_cor_beta_mat <- NULL
wd_cor_beta_mat <- cor( x = wd_data_list$betas, 
                        use = "pairwise.complete.obs", method = "pearson" )

wd_heatmap_pdf <- file.path( plot_dir, "heatmap_v2-Wanding.MVP.betas.pdf")
pdf( file = wd_heatmap_pdf, width = 10, height = 10 )
heatmap( x = wd_cor_beta_mat )
dev.off()



# call_path <- safe_mkdir( file.path( opt$out_path, "calls") )
#
# v2_calls_dat <- NULL
# v2_calls_dat <- call_list_to_matrix( 
#   calls = v2_calls_lst, 
#   out_path = file.path( call_path, "EPICv2"), 
#   out_name = "EPICv2"  )


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                              Header Notes::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

##fileformat=VCFv4.0
##fileDate=20230710
##reference=hg19
##INFO=<ID=PVF,Number=1,Type=Float,Description="Pseudo Variant Frequency">
##INFO=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=GS,Number=1,Type=Integer,Description="Genotyping score from 7 to 85">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

#
# Bret's Suggestions::
# [Done]: Implement coding below:
##INFO=<ID=PVF,Number=1,Type=Integer,Description="Pseudo Variant Frequency * 100">
##INFO=<ID=RGT,Number=1,Type=Integer,Description="Reference Genotype">
##INFO=<ID=AGT,Number=1,Type=Integer,Description="Allele Genotype">
##INFO=<ID=GTS,Number=1,Type=Integer,Description="Genotyping score from 7 to 85">
#
# [TBD]: Make weighted score (for plotting): [-1,1] * QUAL/PVF/GS | FILTER
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                            Build VCF Lists::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

v1_vcf_list <- NULL
v1_vcf_list <- file_list( path    = v1_dat_path, 
                          prefix  = v1_dat_path,
                          suffix  = "_EPIC_B4_raw.snps.vcf", 
                          pattern = "_EPIC_B4_raw.snps.vcf$",
                          recursive = TRUE )

v2_vcf_list <- NULL
v2_vcf_list <- file_list( path    = v2_dat_path, 
                          prefix  = v2_dat_path,
                          suffix  = "_EPIC_A1_raw.snps.vcf", 
                          pattern = "_EPIC_A1_raw.snps.vcf$",
                          recursive = TRUE )

wd_vcf_list <- NULL
wd_vcf_list <- file_list( path    = wd_dat_path, 
                          prefix  = wd_dat_path,
                          suffix  = "_EPIC_A1_raw.snps.vcf", 
                          pattern = "_EPIC_A1_raw.snps.vcf$",
                          recursive = TRUE )



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                         Fingerprint:: All-SNVs
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_dats <- list()
opt$gts_mins <- c( 0, 20 )
opt$sub_strs <- c( "All", "MVP")

for ( sub_str in opt$sub_strs ) {
  sub_vec <- NULL
  if ( sub_str == "MVP" ) sub_vec <- va_ids_tib$Probe_ID
  
  for ( gts_min in opt$gts_mins ) {
    gts_str <- paste0("gts",gts_min)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            EPICv1 SNP Analysis::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    v1_snp_dat <- NULL
    v1_snp_key <- paste0( "EPICv1","_",gts_str,"_",sub_str  )
    if ( p1 ) cat(glue::glue("{pmssg} Fingerprint Params: {v1_snp_key}...{RET}"))
    v1_snp_dat <- fingerprint_vcfs( vcfs = v1_vcf_list, 
                                    sam_tib = v1_sam_tib,
                                    sam_vec = plot_sample_vec,
                                    sub_vec = sub_vec,
                                    run_sig = TRUE,
                                    gts_min = gts_min,
                                    out_dir = file.path( opt$out_path ), 
                                    run_tag = v1_snp_key,
                                    reload = opt$reload, 
                                    reload_min = 10, 
                                    ret_data  = TRUE,
                                    write_out = FALSE, 
                                    write_sum = TRUE, 
                                    write_sig = TRUE,
                                    vb=vb,vt=vt+1,tc=tc+1, tt=tt )
    snp_dats[[v1_snp_key]] = v1_snp_dat
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            EPICv2 SNP Analysis::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    v2_snp_dat <- NULL
    v2_snp_key <- paste0( "EPICv2","_",gts_str,sub_str  )
    if ( p1 ) cat(glue::glue("{pmssg} Fingerprint Params: {v2_snp_key}...{RET}"))
    v2_snp_dat <- fingerprint_vcfs( vcfs = v2_vcf_list, 
                                    sam_tib = v2_sam_tib,
                                    sam_vec = plot_sample_vec,
                                    sub_vec = sub_vec,
                                    run_sig = TRUE,
                                    gts_min = gts_min,
                                    out_dir = file.path( opt$out_path ),
                                    run_tag = v2_snp_key,
                                    reload = opt$reload, 
                                    reload_min = 10, 
                                    ret_data  = TRUE,
                                    write_out = FALSE, 
                                    write_sum = TRUE, 
                                    write_sig = TRUE,
                                    vb=vb,vt=vt+1,tc=tc+1, tt=tt )
    snp_dats[[v2_snp_key]] = v2_snp_dat
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          EPICv2 Wanding SNP Analysis::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    wd_snp_dat <- NULL
    wd_snp_key <- paste0( "Wanding","_",gts_str,sub_str  )
    if ( p1 ) cat(glue::glue("{pmssg} Fingerprint Params: {wd_snp_key}...{RET}"))
    wd_snp_dat <- fingerprint_vcfs( vcfs = wd_vcf_list, 
                                    sam_tib = wd_sam_tib,
                                    sam_vec = wd_sam_tib$Sample_Base,
                                    sub_vec = sub_vec,
                                    run_sig = TRUE,
                                    gts_min = gts_min,
                                    out_dir = file.path( opt$out_path ),
                                    run_tag = wd_snp_key,
                                    reload = opt$reload, 
                                    reload_min = 10, 
                                    ret_data  = TRUE,
                                    write_out = FALSE, 
                                    write_sum = TRUE, 
                                    write_sig = TRUE,
                                    vb=vb,vt=vt+1,tc=tc+1, tt=tt )
    snp_dats[[wd_snp_key]] = wd_snp_dat
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Plot Heatmaps::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # snp_to_heatmap()
    # c( v1_snp_dat$snp_tab,v2_snp_dat$snp_tab,wd_snp_dat$snp_tab )
    
    
    
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Compare SNP Lists:: VA/EPICv1/v2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

v1_ids_tib <- NULL
v1_ids_tib <- v1_snp_dat$snp_tab %>% 
  dplyr::distinct( Target_ID ) %>%
  tidyr::separate( Target_ID, into=c("Probe_ID","Chr","Pos","Alt","Ref"), 
                   sep=":", remove = TRUE, convert = TRUE ) %>%
  dplyr::mutate( Source = "V1",
                 Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"),
                 Probe_Type = Probe_ID %>% stringr::str_sub(1,2)
  ) %>%
  dplyr::select( Source, Probe_ID, Loci_ID, Probe_Type, dplyr::everything() )

v2_ids_tib <- NULL
v2_ids_tib <- v2_snp_dat$snp_tab %>% 
  dplyr::distinct( Target_ID ) %>%
  tidyr::separate( Target_ID, into=c("Probe_ID","Chr","Pos","Alt","Ref"), 
                   sep=":", remove = TRUE, convert = TRUE ) %>%
  dplyr::mutate( Source = "V2",
                 Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"),
                 Probe_Type = Probe_ID %>% stringr::str_sub(1,2)
  ) %>%
  dplyr::select( Source, Probe_ID, Loci_ID, Probe_Type, dplyr::everything() )

ids_tib <- NULL
ids_tib <- dplyr::bind_rows( va_ids_tib, v1_ids_tib, v2_ids_tib ) %>% 
  dplyr::add_count( Probe_ID, name="Probe_Cnt" )

ids_sum <- NULL
ids_sum <- ids_tib %>% 
  dplyr::distinct( Source,Loci_ID, .keep_all = TRUE ) %>%
  dplyr::group_by( Source,Probe_Type ) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Overlap Comparisons::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# VA/V2 Overlap::
#
v2_mat_ids <- NULL
v2_mat_ids <- va_ids_tib %>% 
  dplyr::inner_join( v2_ids_tib %>% dplyr::distinct( Loci_ID ), 
                     by=c("Loci_ID"), 
                     suffix=c("_va","_v2") )

v2_mat_sum <- NULL
v2_mat_sum <- v2_mat_ids %>%
  dplyr::group_by( Probe_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) v2_mat_sum %>% print( n=base::nrow(v2_mat_sum) )

#
# V1 Replaecments (Old)::
#
v2_unq_ids <- NULL
v2_unq_ids <- v2_ids_tib %>% 
  dplyr::distinct( Loci_ID, .keep_all = TRUE ) %>% 
  dplyr::anti_join( va_ids_tib, 
                    by=c("Loci_ID") )

v2_unq_sum <- NULL
v2_unq_sum <- v2_unq_ids %>%
  dplyr::group_by( Probe_Type )%>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) v2_unq_sum %>% print( n=base::nrow(v2_unq_sum) )

#
# V2 Replacements (New)::
#
v2_new_ids <- NULL
v2_new_ids <- v2_unq_ids %>% 
  dplyr::filter( !Loci_ID %in% v1_ids_tib$Loci_ID )

v2_new_sum <- NULL
v2_new_sum <- v2_new_ids %>%
  dplyr::group_by( Probe_Type )%>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) v2_new_sum %>% print( n=base::nrow(v2_new_sum) )




# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Overlap Comparisons:: TRY2
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# move to new plotting page
# grid.newpage()

# create Venn diagram with three sets
# draw.triple.venn(area1=40, area2=15, area3=10, 
#                  n12=5, n23=12, n13=4, n123=2, 
#                  category=c("Science","Economics","English"),
#                  col="Red",fill=c("Green","Yellow","Blue"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

va_ids_tib %>% dplyr::distinct( Loci_ID, .keep_all = TRUE ) %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
v1_ids_tib %>% dplyr::distinct( Loci_ID, .keep_all = TRUE ) %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
v2_ids_tib %>% dplyr::distinct( Loci_ID, .keep_all = TRUE ) %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )

cg_venn_pdf <- file.path( plot_dir, "venndiagram.cg.pdf")
pdf( file = cg_venn_pdf, width = 10, height = 10 )
grid.newpage()
draw.triple.venn(area1=181, area2=824, area3=1234, 
                 n12=181, n23=407, n13=31, n123=31, 
                 category=c("VA(cg)","V1(cg)","V2(cg)"),
                 col="Red",fill=c("Green","Yellow","Blue"))
dev.off()


rs_venn_pdf <- file.path( plot_dir, "venndiagram.rs.pdf")
pdf( file = rs_venn_pdf, width = 10, height = 10 )
grid.newpage()
draw.triple.venn(area1=30, area2=44, area3=61, 
                 n12=30, n23=43, n13=29, n123=29, 
                 category=c("VA(rs)","V1(rs)","V2(rs)"),
                 col="Red",fill=c("Green","Yellow","Blue"))
dev.off()



A2_mat_ids <- NULL
A2_mat_ids <- va_ids_tib %>% 
  dplyr::inner_join( v2_ids_tib %>% dplyr::distinct( Loci_ID ), 
                     by=c("Loci_ID"), 
                     suffix=c("_va","_v2") )

A2_mat_sum <- NULL
A2_mat_sum <- A2_mat_ids %>%
  dplyr::group_by( Probe_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) A2_mat_sum %>% print( n=base::nrow(A2_mat_sum) )


O2_mat_ids <- NULL
O2_mat_ids <- v1_ids_tib %>% 
  dplyr::inner_join( v2_ids_tib %>% dplyr::distinct( Loci_ID ), 
                     by=c("Loci_ID"), 
                     suffix=c("_va","_v2") )

O2_mat_sum <- NULL
O2_mat_sum <- O2_mat_ids %>%
  dplyr::group_by( Probe_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) O2_mat_sum %>% print( n=base::nrow(O2_mat_sum) )

A1_mat_ids <- NULL
A1_mat_ids <- va_ids_tib %>% 
  dplyr::inner_join( v1_ids_tib %>% dplyr::distinct( Loci_ID ), 
                     by=c("Loci_ID"), 
                     suffix=c("_va","_v2") )

A1_mat_sum <- NULL
A1_mat_sum <- A1_mat_ids %>%
  dplyr::group_by( Probe_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) A1_mat_sum %>% print( n=base::nrow(A1_mat_sum) )


A12_mat_ids <- NULL
A12_mat_ids <- va_ids_tib %>% 
  dplyr::inner_join( v1_ids_tib %>% dplyr::distinct( Loci_ID ), 
                     by=c("Loci_ID"), 
                     suffix=c("_va","_v2") ) %>%
  dplyr::inner_join( v2_ids_tib %>% dplyr::distinct( Loci_ID ), 
                     by=c("Loci_ID"), 
                     suffix=c("_va","_v2") )

A12_mat_sum <- NULL
A12_mat_sum <- A12_mat_ids %>%
  dplyr::group_by( Probe_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) A12_mat_sum %>% print( n=base::nrow(A12_mat_sum) )



#
# Best Plot to Prove the Point::
# [Done]: Use RAF/MAF to plot distributions!!!
#
af_dat_tib <- NULL
af_dat_tib <- dplyr::bind_rows(
  va_ids_tib %>% dplyr::inner_join( snv_tib, by=c("Loci_ID"), multiple = "all" ) %>% dplyr::mutate( AF=RAF, Color="RAF"),
  va_ids_tib %>% dplyr::inner_join( snv_tib, by=c("Loci_ID"), multiple = "all" ) %>% dplyr::mutate( AF=MAF, Color="MAF"),
  v2_ids_tib %>% dplyr::inner_join( snv_tib, by=c("Loci_ID"), multiple = "all" ) %>% 
    dplyr::filter( Loci_ID %in% v2_mat_ids$Loci_ID ) %>% dplyr::mutate( Source="v2_VA") %>% dplyr::mutate( AF=RAF, Color="RAF"),
  v2_ids_tib %>% dplyr::inner_join( snv_tib, by=c("Loci_ID"), multiple = "all" ) %>% 
    dplyr::filter( Loci_ID %in% v2_mat_ids$Loci_ID ) %>% dplyr::mutate( Source="v2_VA") %>% dplyr::mutate( AF=MAF, Color="MAF")
)

raf_hisD_pdf <- file.path( plot_dir, "AF-VA-V2.density.pdf" )
raf_hisD_ggg <- af_dat_tib %>%
  # ggplot2::ggplot( aes(x=RAF, fill=Probe_Type) ) +
  ggplot2::ggplot( aes(x=AF, fill=Color ) ) +
  ggplot2::geom_density( alpha=0.2 ) +
  ggplot2::facet_grid( cols = vars(Source), 
                       rows = vars( Probe_Type),
                       scales = "free_y" )
ggplot2::ggsave( filename = raf_hisD_pdf, plot = raf_hisD_ggg, device = "pdf", width = 7, height = 7, dpi = 320 )


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                         Fingerprint:: VA-ONLY
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# snp_dats <- NULL
opt$gts_mins <- c( 0, 20 )
for ( gts_min in opt$gts_mins ) {
  # cur_dats <- NULL
  
  gts_min <- 20
  
  sub_vec <- NULL
  sub_str <- ""
  
  # sub_vec <- va_ids_tib$Probe_ID
  sub_vec <- v2_mat_ids$Probe_ID
  sub_str <- "-MVP"
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            EPICv1 SNP Analysis::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  v1_vcf_list <- NULL
  v1_vcf_list <- file_list( path    = v1_dat_path, 
                            prefix  = v1_dat_path,
                            suffix  = "_EPIC_B4_raw.snps.vcf", 
                            pattern = "_EPIC_B4_raw.snps.vcf$",
                            recursive = TRUE )
  v1_snp_dat <- NULL
  v1_snp_dat <- fingerprint_vcfs( vcfs = v1_vcf_list, 
                                  sam_tib = v1_sam_tib,
                                  sam_vec = plot_sample_vec,
                                  sub_vec = sub_vec,
                                  run_sig = TRUE,
                                  gts_min = gts_min,
                                  out_dir = file.path( opt$out_path), 
                                  run_tag = paste0( "EPICv1",sub_str ), 
                                  # run_tag = paste0( opt$run_name ), 
                                  reload = opt$reload, 
                                  reload_min = 10, 
                                  ret_data  = TRUE,
                                  write_out = FALSE, 
                                  write_sum = TRUE, 
                                  write_sig = TRUE,
                                  vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  # cur_dats[["v1"]] <- v1_snp_dat
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            EPICv2 SNP Analysis::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  v2_vcf_list <- NULL
  v2_vcf_list <- file_list( path    = v2_dat_path, 
                            prefix  = v2_dat_path,
                            suffix  = "_EPIC_A1_raw.snps.vcf", 
                            pattern = "_EPIC_A1_raw.snps.vcf$",
                            recursive = TRUE )
  
  v2_snp_dat <- NULL
  v2_snp_dat <- fingerprint_vcfs( vcfs = v2_vcf_list, 
                                  sam_tib = v2_sam_tib,
                                  sam_vec = plot_sample_vec,
                                  sub_vec = sub_vec,
                                  run_sig = TRUE,
                                  gts_min = gts_min,
                                  out_dir = file.path( opt$out_path), 
                                  run_tag = paste0( "EPICv2",sub_str ), 
                                  # run_tag = paste0( opt$run_name ), 
                                  reload = opt$reload, 
                                  reload_min = 10, 
                                  ret_data  = TRUE,
                                  write_out = FALSE, 
                                  write_sum = TRUE, 
                                  write_sig = TRUE,
                                  vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  # cur_dats[["v2"]] <- v2_snp_dat
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                          EPICv2 Wanding SNP Analysis::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  wd_vcf_list <- NULL
  wd_vcf_list <- file_list( path    = wd_dat_path, 
                            prefix  = wd_dat_path,
                            suffix  = "_EPIC_A1_raw.snps.vcf", 
                            pattern = "_EPIC_A1_raw.snps.vcf$",
                            recursive = TRUE )
  
  wd_snp_dat <- NULL
  wd_snp_dat <- fingerprint_vcfs( vcfs = wd_vcf_list, 
                                  sam_tib = wd_sam_tib,
                                  sam_vec = wd_sam_tib$Sample_Base,
                                  sub_vec = sub_vec,
                                  run_sig = TRUE,
                                  gts_min = gts_min,
                                  out_dir = file.path( opt$out_path), 
                                  run_tag = paste0( "Wanding",sub_str ), 
                                  # run_tag = paste0( opt$run_name ), 
                                  reload = opt$reload, 
                                  reload_min = 10, 
                                  ret_data  = TRUE,
                                  write_out = FALSE, 
                                  write_sum = TRUE, 
                                  write_sig = TRUE,
                                  vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  # cur_dats[["v2"]] <- v2_snp_dat
  
  break
}

# GREAT EXAMPLE::
#. v1_snp_dat$fig_gtc_lst$NA12873
#. v2_snp_dat$fig_gtc_lst$NA12873

#
# [TBD]: Split fingerprint_vcfs() so we can combine data sets...
# [TBD]: Create HeatPlot of best matching fingerprint...
#
# LEFT OFF HERE!!!
#
# [TBD]: Functionalize Heatmap...
# [TBD]: Move Sample Sheet Plot_Name code above!
# [TBD]: Write Plot_Name and Full Name map files.
# [TBD]: Write Max Genotype Funciton...
# [TBD]: Clean up code space.
# [TBD]: Convert to RMD (R-Mark-Down) format...
#

#
# Just for plotting
#
if ( FALSE ) {
  
  dplyr::bind_rows( 
    dplyr::bind_rows(v1_snp_dat$fig_gtc_lst) %>% dplyr::mutate( Name="v1" ),
    dplyr::bind_rows(v2_snp_dat$fig_gtc_lst) %>% dplyr::mutate( Name="v2" )
  ) %>% dplyr::select( Name, dplyr::everything() ) %>%
    dplyr::arrange( Sample_Base,PPP ) %>%
    dplyr::rename( Sample=Sample_Base ) %>%
    dplyr::filter( Sample=="HELA" | Sample=="NA12873") %>%
    dplyr::mutate( FingerPrint = paste0("s",FingerPrint) )
  
  
  dplyr::bind_rows( 
    dplyr::bind_rows(v1_snp_dat$fig_gtc_lst) %>% dplyr::mutate( Name="v1" ),
    dplyr::bind_rows(v2_snp_dat$fig_gtc_lst) %>% dplyr::mutate( Name="v2" )
  ) %>% dplyr::select( Name, dplyr::everything() ) %>%
    dplyr::arrange( Sample_Base,PPP ) %>%
    dplyr::rename( Sample=Sample_Base ) %>%
    dplyr::filter( Sample=="HELA" | Sample=="NA12873" | Sample=="RAJI") %>%
    dplyr::mutate( FingerPrint = paste0("s",FingerPrint) )
  
  
}

sep_vec <- NULL
sep_vec <- c( 1:(length(sub_vec)-1) )
snp_vec <- NULL
snp_vec <- paste0( "v",c( 1:length(sub_vec) ) )
# snp_vec <- paste0( "v",sep_vec )

va_fig_mat <- NULL
va_fig_mat <- dplyr::bind_rows( 
  dplyr::bind_rows(v1_snp_dat$fig_gtc_lst) %>% dplyr::mutate( Name="v1", PPP=as.integer(PPP) ),
  dplyr::bind_rows(v2_snp_dat$fig_gtc_lst) %>% dplyr::mutate( Name="v2", PPP=as.integer(PPP) )
) %>%
  dplyr::group_by( Sample_Base,Name ) %>%
  dplyr::mutate( Rep = dplyr::row_number() ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate( 
    Sample_Base = Sample_Base %>% stringr::str_sub(1,1),
    Name=paste0( Name,"_",Sample_Base,Rep,"_",PPP) ) %>%
  dplyr::select( Name,FingerPrint ) %>%
  tidyr::separate( FingerPrint, into = snp_vec, sep = sep_vec, remove = TRUE,
                   convert = TRUE ) %>%
  tibble::column_to_rownames( var = "Name" ) %>%
  as.matrix() %>% t()

# [1] 60 32
va_fig_mat %>% dim()

va_fig_mat[ which(va_fig_mat == 5 ) ] <- NA_real_

col_vec <- c(1:base::ncol(va_fig_mat) )
row_vec <- c(1:base::nrow(va_fig_mat) )

heat_mat <- NULL
heat_mat <- matrix( data=0, nrow = length(col_vec), ncol=length(col_vec) )
colnames(heat_mat) <- colnames(va_fig_mat)
rownames(heat_mat) <- colnames(va_fig_mat)

for ( ii in col_vec ) {
  for ( jj in col_vec ) {
    mat_cnt <- which( va_fig_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 ) %>% length()
    mis_cnt <- which( va_fig_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() != 0 ) %>% length()
    nan_cnt <- which( is.na(va_fig_mat[ ,ii ]) | is.na(va_fig_mat[ ,jj ]) ) %>% length()
    
    if ( mat_cnt + mis_cnt + nan_cnt != length(row_vec) ) {
      stop(glue::glue("{pmssg} Failed: ii={ii}, jj={jj}, mat={mat_cnt}, mis={mis_cnt}, nan={nan_cnt}.{RET2}"))
    }
    
    heat_mat[ii,jj] = mat_cnt
    
    # break
  }
  # break
}

heatmap_pdf <- file.path( plot_dir, "heatmap_v1-v2.MVP.pdf")
pdf( file = heatmap_pdf, width = 10, height = 10 )
heatmap( heat_mat )
dev.off()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                            Wanding Heatmap::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  # wd_snp_dat
  
  sep_vec <- NULL
  sep_vec <- c( 1:(length(sub_vec)-1) )
  snp_vec <- NULL
  snp_vec <- paste0( "v",c( 1:length(sub_vec) ) )
  # snp_vec <- paste0( "v",sep_vec )
  
  wd_fig_mat <- NULL
  wd_fig_mat <- dplyr::bind_rows( 
    # dplyr::bind_rows(wd_snp_dat$fig_gtc_lst) %>% dplyr::mutate( Name="v1", PPP=as.integer(PPP) ),
    dplyr::bind_rows(wd_snp_dat$fig_gtc_lst) %>% dplyr::mutate( Name="v2", PPP=as.integer(PPP) )
  ) %>%
    dplyr::group_by( Sample_Base,Name ) %>%
    dplyr::mutate( Rep = dplyr::row_number() ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate( 
      # Sample_Base = Sample_Base %>% stringr::str_sub(1,1),
      Name=paste0( Name,"_",Sample_Base,"r",Rep,"_",PPP) ) %>%
    dplyr::select( Name,FingerPrint ) %>%
    tidyr::separate( FingerPrint, into = snp_vec, sep = sep_vec, remove = TRUE,
                     convert = TRUE ) %>%
    tibble::column_to_rownames( var = "Name" ) %>%
    as.matrix() %>% t()
  
  # [1] 60 38
  wd_fig_mat %>% dim()
  
  wd_fig_mat[ which(wd_fig_mat == 5 ) ] <- NA_real_
  
  col_vec <- c(1:base::ncol(wd_fig_mat) )
  row_vec <- c(1:base::nrow(wd_fig_mat) )
  
  wd_heat_mat <- NULL
  wd_heat_mat <- matrix( data=0, nrow = length(col_vec), ncol=length(col_vec) )
  colnames(wd_heat_mat) <- colnames(wd_fig_mat)
  rownames(wd_heat_mat) <- colnames(wd_fig_mat)
  
  for ( ii in col_vec ) {
    for ( jj in col_vec ) {
      mat_cnt <- which( wd_fig_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 ) %>% length()
      mis_cnt <- which( wd_fig_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() != 0 ) %>% length()
      nan_cnt <- which( is.na(wd_fig_mat[ ,ii ]) | is.na(wd_fig_mat[ ,jj ]) ) %>% length()
      
      if ( mat_cnt + mis_cnt + nan_cnt != length(row_vec) ) {
        stop(glue::glue("{pmssg} Failed: ii={ii}, jj={jj}, mat={mat_cnt}, mis={mis_cnt}, nan={nan_cnt}.{RET2}"))
      }
      
      wd_heat_mat[ii,jj] = mat_cnt
      
      # break
    }
    # break
  }
  
  wd_heatmap_pdf <- file.path( plot_dir, "heatmap_v2-Wanding-Names.MVP.pdf")
  pdf( file = wd_heatmap_pdf, width = 10, height = 10 )
  heatmap( wd_heat_mat )
  dev.off()
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Beta Pairs Plotting:: raw/ind
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$load_calls <- TRUE
opt$load_calls <- FALSE

opt$plot_calls <- TRUE
opt$plot_calls <- FALSE

if ( opt$load_calls ) {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Load Manifests:: B4/A1
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
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Load Calls:: raw/ind
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # [TBD]: Functionalize all of this so we can pass a list of files to generate 
  #. a set of plots...
  #
  
  v2_calls_list <- NULL
  v2_calls_list[["raw"]] <- NULL
  v2_calls_list[["raw"]] <- file_list( 
    path    = v2_dat_path, 
    prefix  = v2_dat_path,
    suffix  = "_EPIC_A1_raw.call.dat.csv.gz", 
    pattern = "_EPIC_A1_raw.call.dat.csv.gz$",
    recursive = TRUE )
  
  v2_calls_list[["ind"]] <- NULL
  v2_calls_list[["ind"]] <- file_list( 
    path    = v2_dat_path, 
    prefix  = v2_dat_path,
    suffix  = "_EPIC_A1_ind.call.dat.csv.gz", 
    pattern = "_EPIC_A1_ind.call.dat.csv.gz$",
    recursive = TRUE )
  
  analysis_val <- "CellLine"
  samp_name <- "swap"
  detp_keys <- c("poob","negs")
  beta_keys <- c("ind","raw")
  prod_keys <- c("EPICv2")
  pval_mins <- c(1.0, 0.05 )
  # pval_mins <- c(1.0, 0.05, 0.01 )
  dB_mins   <- c( 0.2 )
  # dB_mins   <- c( 0.2, 0.1 )
  
  #
  # Simplest::
  #
  # analysis_val <- "CellLine"
  # samp_name <- "swap"
  # detp_keys <- c("poob")
  # beta_keys <- c("ind")
  # prod_keys <- c("EPICv2")
  # pval_mins <- c(0.05 )
  # dB_mins   <- c( 0.2 )
  
  # [TBD]: Investigate idn instead of ind...
  
  # [TBD]: To use this correctly we need to have EPICv1/EPICv2 together...
  # [TBD]: Load all samples in a single pass instead of v1/v2
  
  plot_single <- TRUE
  plot_single <- FALSE
  
  v3_man_sum <- NULL
  v3_man_sum <- v3_man_tib %>% 
    dplyr::group_by( Probe_Group_Mask,BP_Group,Probe_Group_BP ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) v3_man_sum %>% print( n=base::nrow(v3_man_sum) )
  
  for ( sample_base in plot_sample_vec ) {
    # Need to extract the correct Sentrix_Names
    sentrix_vec <- NULL
    sentrix_vec <- v2_sam_tib %>% 
      dplyr::filter( Sample_Base == sample_base ) %>% 
      dplyr::arrange( cg_calls_pass_perc_1 ) %>%
      dplyr::pull( Sentrix_Name )
    
    # if ( sentrix_vec %>% length() > 3 ) {
    #   sentrix_vec = c( sentrix_vec[1],
    #                    sentrix_vec[ as.integer(sentrix_vec %>% length() / 2) ],
    #                    sentrix_vec[ sentrix_vec %>% length() ]
    #   )
    # }
    if ( sentrix_vec %>% length() > 4 ) {
      sentrix_len <- sentrix_vec %>% length()
      sentrix_vec <- c( sentrix_vec[1],
                        sentrix_vec[2],
                        sentrix_vec[ sentrix_len - 1 ],
                        sentrix_vec[ sentrix_len ] )
    }
    
    for ( beta_key in beta_keys ) {
      beta_str <- paste0("beta-",beta_key)
      
      call_list <- NULL
      call_list <- v2_calls_list[[beta_key]]
      data_list <- NULL
      data_list <- call_list[sentrix_vec] %>%
        lapply( readr::read_csv, show_col_types = FALSE ) %>%
        dplyr::bind_rows( .id = "Sentrix_Name" ) %>%
        tidyr::pivot_longer( cols = c( pvals_pOOBAH,pvals_PnegEcdf,betas ), 
                             names_to  = c("Key"), 
                             values_to = c("Val") ) %>% split(.$Key) %>% 
        lapply( function(x) { 
          x %>% dplyr::select(-Key) %>% 
            tidyr::pivot_wider( id_cols = c(Probe_ID), 
                                names_from = c(Sentrix_Name), 
                                values_from = c(Val) )
        })
      
      betas <- NULL
      betas <- data_list$betas
      
      for ( detp_key in detp_keys ) {
        detp_str <- paste0("detp-",detp_key)
        
        detps <- NULL
        detps <- data_list$pvals_PnegEcdf
        if ( detp_key == "poob" ) detps <- data_list$pvals_pOOBAH
        
        for ( prod_key in prod_keys ) {
          prod_str <- paste0("prod-",prod_key)
          for ( pval_min in pval_mins ) {
            pval_str <- paste0("pval-",pval_min)
            for ( dB_min in dB_mins ) {
              dB_str <- paste0("dB-",dB_min)
              
              # break } break } break } break } break }
              cat(glue::glue("{pmssg} Params: {detp_key}, {beta_key}, {prod_key}, {pval_min}, {dB_min}...{RET}"))
              
              sample_red <- 0
              sample_red <- 3
              # betas_ncol <- base::ncol( betas %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) )
              # detps_ncol <- base::ncol( detps %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) )
              betas_ncol <- base::ncol( betas )
              detps_ncol <- base::ncol( detps )
              
              # if ( betas_ncol > sample_red + 1 ) {
              #   betas_ncol <- base::ncol( betas ) - sample_red
              #   detps_ncol <- base::ncol( detps ) - sample_red
              # }
              plot_dir <- safe_mkdir( file.path( opt$out_path, sample_base,beta_str,detp_str,prod_str,pval_str,dB_str ) )
              
              if ( opt$plot_calls ) {
                plot_ret <- NULL
                plot_ret <- plot_beta_gg( 
                  # betas = betas %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) %>% dplyr::select( dplyr::all_of( 1:betas_ncol) ),
                  # detps = detps %>% dplyr::select( dplyr::all_of( c("Probe_ID",sentrix_vec) ) ) %>% dplyr::select( dplyr::all_of( 1:detps_ncol) ),
                  betas = betas,
                  detps = detps,
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
                  
                  out_dir = file.path( plot_dir ),
                  run_tag = paste( opt$run_name,
                                   detp_str,beta_str,prod_str,pval_str,dB_str,sample_base,
                                   sep="_"),
                  
                  reload     = opt$reload,
                  reload_min = 10,
                  reload_pre = NULL,
                  
                  # ret_data   = FALSE,
                  ret_data   = TRUE,
                  parallel   = FALSE,
                  
                  vb=vb,vt=vt,tc=tc, tt=tt )
              }
              cat(glue::glue("{pmssg} Done.{RET2}"))
              
              if ( plot_single ) break
            }
            if ( plot_single ) break
          }
          if ( plot_single ) break
        }
        if ( plot_single ) break
      }
      if ( plot_single ) break
    }
    if ( plot_single ) break
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file

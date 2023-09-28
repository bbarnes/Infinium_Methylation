
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
annoS_rds <- file.path( "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/dat/annoS.rds" )
annoI_rds <- NULL
annoI_rds <- file.path( "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/dat/annoI.rds" )

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

top_dir <- "/Users/bbarnes/Documents/Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.4.4"
snps1_vcf <- file.path( top_dir, "EPICv1/chip-206203800149/swifthoof_main/206203800149_R01C01_EPIC_B4_ind.snps.vcf" )
beta1_csv <- file.path( top_dir, "EPICv1/chip-206203800149/swifthoof_main/206203800149_R01C01_EPIC_B4_ind.call.dat.csv.gz" )
beta2_csv <- file.path( top_dir, "EPICv2/chip-206891110001/swifthoof_main/206891110001_R01C01_EPIC_A1_ind.call.dat.csv.gz" )

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
lociS2_tib <- beta2_tib %>% dplyr::filter( Loci_ID %in% names(annoS_dat) )

# Missing Probes::
#
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

#
# Manifest Investigation::
#
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

epicv2_man_csv <- file.path( "/Users/bbarnes/Documents/data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv" )
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

#
# Attempt:: Range (RG)
#
epicv2_bed_RG_tib <- NULL
epicv2_bed_RG_tsv <- file.path( "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/bed/EPIC-8v2-0_A1.next-base-RG.bed" )
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

# Before Range Reduction::
# > epicv2_bed_RG_tib
# A tibble: 126,710 × 8

# Tabix Command::
# tabix -h -R Projects.new/EPIC_v2/docker/bed/EPIC-8v2-0_A1.next-base.bed data/dbSNP/NCBI/common_all_20180418.vcf.gz > Projects.new/EPIC_v2/docker/intersect/EPICv2.dbSNP151-common.next.base.vcf 
#
dbSNP_151_bgzip <- file.path( "/Users/bbarnes/Documents/data/dbSNP/NCBI/common_all_20180418.vcf.gz" )
epicv2_int_dbSNP_RG_vcf <- file.path( "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/intersect/EPICv2.dbSNP151-common.next-base-RG.vcf" )
tabix_RG_cmd <- paste0("tabix -h -R ",epicv2_bed_RG_tsv," ",dbSNP_151_bgzip," > ",epicv2_int_dbSNP_RG_vcf)
base::system( tabix_RG_cmd )

# Fix the names stuff below...

# Before Range Reduction::
# > epicv2_int_dbSNP_RG_tib
# A tibble: 11,705 × 8
epicv2_int_dbSNP_RG_tib <- NULL
epicv2_int_dbSNP_RG_tib <- readr::read_tsv( file = epicv2_int_dbSNP_RG_vcf, 
                                         skip = 56, show_col_types = FALSE) %>% 
  magrittr::set_names( c("Chromosome_SNP_hg38","Coordinate_SNP_hg38","RS_ID",
                         "REF","ALT","QUAL","FILTER","INFO") )

epicv2_sub_dbSNP_RG_tib <- NULL
epicv2_sub_dbSNP_RG_tib <- epicv2_int_dbSNP_RG_tib %>%
  dplyr::filter( stringr::str_length(REF) == 1 ) %>%
  dplyr::filter( stringr::str_length(ALT) == 1 ) %>%
  dplyr::filter( REF != ALT ) %>%
  dplyr::filter( !RS_ID %>% stringr::str_detect(",") )

epicv2_cpg_dbSNP_RG_tib <- NULL
epicv2_cpg_dbSNP_RG_tib <- epicv2_sub_dbSNP_RG_tib %>%
  dplyr::inner_join( epicv2_bed_RG_tib, 
                     by=c("Chromosome_SNP_hg38"="Chromosome_BED_hg38",
                          "Coordinate_SNP_hg38"="End_BED_hg38"),
                     multiple = "all" ) %>%
  dplyr::mutate(
    strand = dplyr::case_when(
      Strand_FR == "F" ~ "+",
      Strand_FR == "R" ~ "-",
      TRUE ~ NA_character_
    )
  )

# Quick Print::
epicv2_cpg_dbSNP_RG_tib %>% head() %>% as.data.frame()

# Delete this later, we can just build the GRS directly below...
#
# epicv2_cpg_dbSNP_RG_df <- NULL
# epicv2_cpg_dbSNP_RG_df <- epicv2_cpg_dbSNP_RG_tib %>% 
#   dplyr::mutate( 
#     seqnames = paste0("chr",Chromosome_SNP_hg38),
#     start = Coordinate_SNP_hg38 %>% as.integer(),
#     end = Coordinate_SNP_hg38 %>% as.integer(),
#     width = end - start + 1 %>% as.integer(),
#     strand = dplyr::case_when(
#       Strand_FR == "F" ~ "+",
#       Strand_FR == "R" ~ "-",
#       TRUE ~ NA_character_
#     ), 
#     rs = RS_ID,
#     designType = "I",
#     In.band = "REF"
#   ) %>%
#   dplyr::select( Probe_ID, # Loci_ID,
#                  seqnames,start,end,strand,
#                  designType,In.band,REF,ALT,rs ) %>%
#   dplyr::distinct() %>%
#   tibble::column_to_rownames( var = "Probe_ID" )

# LEFT OFF HERE::
#.  - Need to look up how to build a GRange again...
#.  - Also need to figure out the whole mapping back issue...
#
epicv2_cpg_dbSNP_RG_grs <- NULL
epicv2_cpg_dbSNP_RG_grs <- GenomicRanges::GRanges( 
  seqnames   = paste0("chr",epicv2_cpg_dbSNP_RG_tib$Chromosome_SNP_hg38), 
  strand     = epicv2_cpg_dbSNP_RG_tib$strand,
  designType = "I",
  In.band    = "REF",
  REF        = epicv2_cpg_dbSNP_RG_tib$REF,
  ALT        = epicv2_cpg_dbSNP_RG_tib$ALT,
  rs         = epicv2_cpg_dbSNP_RG_tib$RS_ID,
  
  IRanges::IRanges(
    start = epicv2_cpg_dbSNP_RG_tib$Coordinate_SNP_hg38,
    end = epicv2_cpg_dbSNP_RG_tib$Coordinate_SNP_hg38,
    names = epicv2_cpg_dbSNP_RG_tib$Probe_ID
  )
)

#
# LEFT OFF HERE::
#   [Done]: Validation Attempt:: Testing RG
#.  TBD:: Using lessons from below filter above to make a new AnnoI_v2
#.   - Remove Indels, Multiple RS#'s...
#.   - Investigate MAF...
#

#
# IMPORTANT: Shouldn't run ind run idn and don't use i for SNP calling!
#
check_validation <- TRUE
check_validation <- FALSE
if ( check_validation ) {
  
  # Expectation Below for Overlapping EPICv2 Infinium I Trifecta Probes::
  #  > val_annoI_dbSNP_tib
  #  A tibble: 449 × 23
  val_annoI_dbSNP_tib <- NULL
  val_annoI_dbSNP_tib <- 
    # epicv2_int_dbSNP_RG_tib %>% 
    epicv2_sub_dbSNP_RG_tib %>% 
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







#
# Poor Coordinate Guessing for Tabix Search::
#
if ( FALSE ) {
  
  epicv2_bed_FR_tib <- NULL
  epicv2_bed_FR_tsv <- file.path( "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/bed/EPIC-8v2-0_A1.next-base-FR.bed" )
  epicv2_bed_FR_tib <- epicv2_sub_tib %>% 
    dplyr::filter( CHR != "0" ) %>%
    dplyr::filter( MAPINFO != 0 ) %>%
    dplyr::mutate( 
      CHR = CHR %>% stringr::str_remove("^chr"), 
      Beg = dplyr::case_when( 
        Strand_FR == "F" ~ MAPINFO - 1, 
        Strand_FR == "R" ~ MAPINFO + 2, 
        TRUE ~ NA_real_ ),
      End = Beg + 1
    ) %>% 
    dplyr::select( CHR,Beg,End,IlmnID, Strand_FR,Strand_TB,Strand_CO ) %>%
    dplyr::arrange( CHR, Beg )
  readr::write_tsv( x = epicv2_bed_FR_tib, file = epicv2_bed_FR_tsv, col_names = FALSE )
  
  # Tabix Command::
  # tabix -h -R Projects.new/EPIC_v2/docker/bed/EPIC-8v2-0_A1.next-base.bed data/dbSNP/NCBI/common_all_20180418.vcf.gz > Projects.new/EPIC_v2/docker/intersect/EPICv2.dbSNP151-common.next.base.vcf 
  #
  dbSNP_151_bgzip <- file.path( "/Users/bbarnes/Documents/data/dbSNP/NCBI/common_all_20180418.vcf.gz" )
  epicv2_int_dbSNP_FR_vcf <- file.path( "/Users/bbarnes/Documents/Projects.new/EPIC_v2/docker/intersect/EPICv2.dbSNP151-common.next-base-FR.vcf" )
  tabix_FR_cmd <- paste0("tabix -h -R ",epicv2_bed_FR_tsv," ",dbSNP_151_bgzip," > ",epicv2_int_dbSNP_FR_vcf)
  base::system( tabix_FR_cmd )
  
  epicv2_int_dbSNP_tib <- NULL
  epicv2_int_dbSNP_tib <- readr::read_tsv( file = epicv2_int_dbSNP_FR_vcf, 
                                           skip = 56, show_col_types = FALSE) %>% 
    magrittr::set_names( c("Chromosome","Coordinate","RS_ID",
                           "REF","ALT","QUAL","FILTER","INFO") )
  
  
  #
  # Validation Attempt:: Failed only 9 matchin in each case::
  #
  if ( FALSE ) {
    #
    # Validation Attempt:: Looks Good!!!
    #
    epicv2_map_FR_dbSNP <- NULL
    epicv2_map_FR_dbSNP <- epicv2_bed_FR_tib %>% 
      dplyr::inner_join( epicv2_int_dbSNP_tib, 
                         by=c("CHR"="Chromosome", "End"="Coordinate"), 
                         multiple = "all" ) %>%
      dplyr::mutate( Loci_ID = IlmnID %>% stringr::str_remove("_.*$") )
    
    annoI_tib %>% 
      dplyr::rename( Loci_ID = Probe_ID,
                     RS_ID = rs ) %>%
      dplyr::inner_join( epicv2_map_FR_dbSNP, by=c("Loci_ID") )
    
    annoI_tib %>% 
      dplyr::rename( Loci_ID = Probe_ID,
                     RS_ID = rs ) %>%
      dplyr::inner_join( epicv2_map_FR_dbSNP, by=c("RS_ID") )
  }
  
  
}

#
# Previous Validation Attempts That Failed!!!
#
if ( FALSE ) {
  #
  # Validation Attempt:: Wrong Coordinate see above...
  #
  epicv2_bed_tib %>% 
    dplyr::inner_join( epicv2_int_dbSNP_tib, 
                       by=c("CHR"="Chromosome", "Beg"="Coordinate"), 
                       multiple = "all" )
  
  #
  # Validation Attempt:: NOT MUCH...
  #
  joinI2_tib %>% dplyr::rename( RS_ID = rs) %>% 
    dplyr::inner_join( epicv2_int_dbSNP_tib, by=c("RS_ID"), 
                       suffix=c("_I2", "_dbSNP") )
  
  #
  # Validation Attempt:: NOT MUCH...
  #
  epicv2_joinI2_tib <- NULL
  epicv2_joinI2_tib <- epicv2_man_tib %>% 
    dplyr::select( IlmnID,Name, 
                   Strand_FR,Strand_TB,Strand_CO,Infinium_Design,Rep_Num,
                   CHR,MAPINFO,
                   SNP_ID,SNP_DISTANCE,SNP_MinorAlleleFrequency ) %>% 
    dplyr::rename( Loci_ID = Name ) %>%
    dplyr::inner_join( joinI2_tib, by=c("Loci_ID") )
  
  # Simple Search for starting 0
  epicv2_joinI2_tib %>% dplyr::filter( SNP_DISTANCE %>% stringr::str_starts( "0;" ) ) %>% as.data.frame()
  
  # Should just search for rs#...
}


#
# NOTES::
#. - Build AnnoS_v2 - remove rs1495031 and trim other missing rs#'s
#. - Build AnnoI_v2 
#.   - TBD:: Search for additional SNPs in manifest
#.   - TBD:: Email Wanding
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                 Original Old School Sesame Function::
#                            Used in Docker Image
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  formatVCF_swifthoof = function (sset, vcf = NULL, refversion = "hg19", annoS = NULL, annoI = NULL) 
  {
    platform <- sset@platform
    if (is.null(annoS)) 
      annoS <- sesameDataPullVariantAnno_SNP(platform)
    betas <- sesame::getBetas(sset)[names(annoS)]
    vafs <- ifelse(annoS$U == "REF", betas, 1 - betas)
    gts <- lapply(vafs, genotype)
    GT <- vapply(gts, function(g) g$GT, character(1))
    GS <- vapply(gts, function(g) g$GS, numeric(1))
    vcflines_snp <- cbind(as.character(GenomicRanges::seqnames(annoS)), 
                          as.character(end(annoS)), names(annoS), annoS$REF, annoS$ALT, 
                          GS, ifelse(GS > 20, "PASS", "FAIL"), sprintf("PVF=%1.3f;GT=%s;GS=%d", 
                                                                       vafs, GT, GS))
    if (is.null(annoI)) 
      annoI <- sesameDataPullVariantAnno_InfiniumI(platform)
    af <- c(pmax(rowSums(oobR(sset)), 1)/(pmax(rowSums(oobR(sset)) + 
                                                 rowSums(IG(sset)), 2)), pmax(rowSums(oobG(sset)), 1)/(pmax(rowSums(oobG(sset)) + 
                                                                                                              rowSums(IR(sset)), 2)))
    af <- af[names(annoI)]
    vafs <- ifelse(annoI$In.band == "REF", af, 1 - af)
    gts <- lapply(vafs, genotype)
    GT <- vapply(gts, function(g) g$GT, character(1))
    GS <- vapply(gts, function(g) g$GS, numeric(1))
    vcflines_typeI <- cbind(as.character(GenomicRanges::seqnames(annoI)), 
                            as.character(end(annoI)), names(annoI), annoI$REF, annoI$ALT, 
                            GS, ifelse(GS > 20, "PASS", "FAIL"), sprintf("PVF=%1.3f;GT=%s;GS=%d", 
                                                                         vafs, GT, GS))
    header <- c("##fileformat=VCFv4.0", sprintf("##fileDate=%s", 
                                                format(Sys.time(), "%Y%m%d")), sprintf("##reference=%s", 
                                                                                       refversion), paste0("##INFO=<ID=PVF,Number=1,Type=Float,", 
                                                                                                           "Description=\"Pseudo Variant Frequency\">"), paste0("##INFO=<ID=GT,Number=1,Type=String,", 
                                                                                                                                                                "Description=\"Genotype\">"), paste0("##INFO=<ID=GS,Number=1,Type=Integer,", 
                                                                                                                                                                                                     "Description=\"Genotyping score from 7 to 85\">"))
    out <- data.frame(rbind(vcflines_snp, vcflines_typeI))
    colnames(out) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                       "FILTER", "INFO")
    rownames(out) <- out$ID
    out <- out[order(out[["#CHROM"]], as.numeric(out[["POS"]])), 
    ]
    if (is.null(vcf)) {
      out
    }
    else {
      writeLines(header, vcf)
      write.table(out, file = vcf, append = TRUE, sep = "\t", 
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
}


# End of file

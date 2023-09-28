
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           imProbeQC Functions:: 
#                             MSA Sample Sheets
#
# Short Summary: A bunch of annoying functions to format Sample Sheets into
#  a usauble formate...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Tidyverse Core Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("tidyverse",  quietly = TRUE) ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Common Human Abbreviations for String Variables::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
TAB2 <- "\t\t"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Idat Reader Functions::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

read_idat_pair_r = function( prefix,
                             
                             neg_vec,
                             man_tib,
                             ctl_tib = NULL,
                             
                             workflow_vec = "i",
                             sesame_work_str = NULL,
                             
                             min_pval = 0.05,
                             min_beta = 0.3,
                             max_beta = 0.7,
                             min_perO = 0.75,
                             min_perI = 0.05,
                             
                             rm_outliers = TRUE,
                             platform = "EPIC",
                             
                             out_dir,
                             run_tag,
                             
                             reload     = 0,
                             reload_min = 2,
                             reload_pre = NULL,
                             
                             ret_data   = FALSE,
                             parallel   = FALSE,
                             write_out  = FALSE,
                             
                             vb=0,vt=3,tc=1, tt=NULL,
                             fun_tag='read_idat_pair_r')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p4  <- vb > vt + 4
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_rds <- file.path( out_dir, paste(out_tag, 'rds', sep='.') )
  beg_txt <- paste(out_rds, 'start.txt', sep='.')
  end_txt <- paste(out_rds, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_rds, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_rds, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}       prefix = '{prefix}'.{RET}"))
    cat(glue::glue("{mssg}     min_pval = '{min_pval}'.{RET}"))
    cat(glue::glue("{mssg}     min_beta = '{min_beta}'.{RET}"))
    cat(glue::glue("{mssg}     max_beta = '{max_beta}'.{RET}"))
    cat(glue::glue("{mssg}     min_perO = '{min_perO}'.{RET}"))
    cat(glue::glue("{mssg}     min_perI = '{min_perI}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid  = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload  = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min  = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre  = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data  = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel  = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     write_out = '{write_out}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_tag = '{out_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_rds = '{out_rds}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, aux_csv, out_rds, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Outdir out='{out_dir}' does not exist")
  if ( !dir.exists( out_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ftime <- base::system.time({
    
    ret_tib <- read_idat_pair_rcpp( 
      prefix_path      = prefix,
      output_path      = out_dir,
      workflow_vec     = workflow_vec %>% as.vector(),
      
      pval_add_vec     = neg_vec %>% as.vector(),
      addU_man_vec     = man_tib$U %>% as.vector(),
      addM_man_vec     = man_tib$M %>% as.vector(),
      
      cgns_man_vec     = man_tib$Probe_ID %>% as.vector(), 
      cols_man_vec     = man_tib$col %>% as.vector(),
      keys_man_vec     = man_tib$Manifest %>% as.vector(),
      anns_man_vec     = man_tib$Annotation %>% as.vector(), 
      chrs_man_vec     = man_tib$Chromosome %>% as.vector(),

      min_pval         = min_pval,
      min_beta         = min_beta,
      max_beta         = max_beta,
      min_perO         = min_perO,
      min_perI         = min_perI,
      
      read_bgz         = FALSE,
      write_bgz        = FALSE, 
      
      rm_pval_outliers = rm_outliers,
      return_df        = ret_data,
      
      vb=vb,vt=vt,tc=tc ) %>%
      # clean_tib() %>%
      # dplyr::mutate(across(where(is.factor), as.character) ) %>%
      dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Run Sesame::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( FALSE ) {
      if ( !is.null(sesame_work_str) ) {
        ret_tib <- ret_tib %>%
          dplyr::mutate( M = dplyr::case_when( M==0 ~ NA_integer_, TRUE ~ M ) ) %>%
          # dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") ) %>%
          as.data.frame() %>%
          sesame::SigDF( platform = platform, 
                         ctl = ctl_tib
          ) %>%
          mutate_sdf_simple(
            steps = sesame_work_str,
            negs_min = 1.0, 
            poob_min = 1.0, 
            vb=vb,vt=vt+1,tc=tc, tt=tt )
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( write_out )
      out_cnt <- safe_write( x = ret_tib, file = out_rds, type = "rds", 
                             done = TRUE, write_spec = TRUE, append = FALSE, 
                             fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Genotyping Functions:: MSA.v.1.0
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

## very simple genotyper
genotyper2 = function( x, model_background=0.1, model_nbeads=40) {
  
  GL <- vapply(
    c(model_background, 0.5, 1-model_background),
    function(af) {
      dbinom(
        round(x*model_nbeads),
        size=model_nbeads, prob=af)}, numeric(1))
  
  ind <- which.max(GL)
  GT <- c('0/0','0/1','1/1')[ind]
  GS <- floor(-log10(1-GL[ind] / sum(GL))*10) # assuming equal prior
  list(GT=GT, GS=GS)
}

formatVCF2 = function( sdf, annoS, annoI, 
                       vcf=NULL, genome="hg19", 
                       verbose = FALSE)
{
  
  platform <- sdfPlatform(sdf, verbose = verbose)
  betas <- getBetas(sdf)[names(annoS)]
  vafs <- ifelse(annoS$U == 'REF', betas, 1-betas)
  gts <- lapply(vafs, genotyper2)
  GT <- vapply(gts, function(g) g$GT, character(1))
  GS <- vapply(gts, function(g) g$GS, numeric(1))
  vcflines_snp <- cbind(as.character(GenomicRanges::seqnames(annoS)),
                        as.character(GenomicRanges::end(annoS)),
                        names(annoS), annoS$REF, annoS$ALT, GS, ifelse(GS>20,'PASS','FAIL'),
                        sprintf("PVF=%1.3f;GT=%s;GS=%d", vafs, GT, GS))
  
  # Basically the result of getAFTypeIbySumAlleles()
  # c(setNames(pmax(dG$MR+dG$UR,1)/pmax(dG$MR+dG$UR+dG$MG+dG$UG,2), dG$Probe_ID),
  #   setNames(pmax(dR$MG+dR$UG,1)/pmax(dR$MR+dR$UR+dR$MG+dR$UG,2), dR$Probe_ID) )
  # red / sum and grn / sum
  
  af <- getAFTypeIbySumAlleles(sdf, known.ccs.only=FALSE)
  af <- af[names(annoI)]
  vafs <- ifelse(annoI$In.band == 'REF', af, 1-af)
  gts <- lapply(vafs, genotyper2)
  GT <- vapply(gts, function(g) g$GT, character(1))
  GS <- vapply(gts, function(g) g$GS, numeric(1))
  vcflines_typeI <- cbind(as.character(GenomicRanges::seqnames(annoI)),
                          as.character(GenomicRanges::end(annoI)),
                          annoI$rs, annoI$REF, annoI$ALT, GS, ifelse(GS>20,'PASS','FAIL'),
                          sprintf("PVF=%1.3f;GT=%s;GS=%d", vafs, GT, GS))
  
  header <- vcf_header(genome)
  out <- data.frame(rbind(vcflines_snp, vcflines_typeI))
  colnames(out) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  rownames(out) <- out$ID
  out <- out[order(out[['#CHROM']], as.numeric(out[['POS']])),]
  
  if(is.null(vcf)) { return(out);
  } else {
    writeLines(header, vcf)
    write.table(out, file=vcf, append=TRUE, sep='\t',
                row.names = FALSE, col.names = FALSE, quote = FALSE) }
}

load_annoS_EPICv1_local_tib <- function( top_path ) {
  
  annoS_rds <- NULL
  annoS_rds <- file.path( top_path, "Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoS.rds" )
  
  annoS_dat <- NULL
  annoS_dat <- readr::read_rds( file = annoS_rds )
  
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
  
  return( annoS_tib )
}

load_annoI_EPICv1_local_tib <- function( top_path ) {
  
  annoI_rds <- NULL
  annoI_rds <- file.path( top_path, "Projects.new/EPIC_v2/docker/dat/EPICv1/EPICv1.annoI.rds" )
  
  annoI_dat <- NULL
  annoI_dat <- readr::read_rds( file = annoI_rds )
  
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
  
  return( annoI_tib )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                 Manifest Formatting Functions:: MSA.v.1.0
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

build_core_msa_man <- function( neg_rds, m03_csv, m10_csv, out_csv,
                                vb=0,vt=3,tc=1, tt=NULL,
                                fun_tag='build_core_msa_man')
{
  org_neg_tib  <- NULL
  org_neg_tib  <- readr::read_rds( neg_rds ) %>%
    dplyr::mutate( Probe_ID = paste0(Probe_ID,"_0") )
  
  # A tibble: 77,838 × 15
  v03_man_tib <- NULL
  v03_man_tib <- readr::read_csv( file = m03_csv, show_col_types = FALSE ) %>%
    clean_tib() %>%
    dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("^ZZ", "ctl_"),
                   Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_neg", "ctl_Negative"),
                   # Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_neg", "ctl_Negative_"),
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
  
  # Extract Controls: All
  v03_ctl_tib <- NULL
  v03_ctl_tib <- v03_man_tib %>%
    dplyr::filter( !Probe_Type == "cg" ) %>%
    dplyr::filter( !Probe_Type == "ch" ) %>%
    dplyr::filter( !Probe_Type == "rs" )
  
  v03_ctl_sum <- NULL
  v03_ctl_sum <- print_sum( tib = v03_ctl_tib, vec = c("Probe_Type"), 
                            vb=vb,vt=vt,tc=tc )
  
  # Extract Controls: Negative Only
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
  
  # Combine All Negative Controls::
  # A tibble: 4,957 × 3
  all_neg_tib <- NULL
  all_neg_tib <- dplyr::bind_rows( org_neg_tib,v03_neg_tib ) %>%
    dplyr::arrange( Sequence) %>%
    dplyr::distinct( Address, .keep_all = TRUE )
  
  # MSA v.1.0 Manifest::
  # A tibble: 293,649 × 27
  v10_man_tib <- NULL
  v10_man_tib <- readr::read_csv( file = m10_csv, show_col_types = FALSE ) %>% clean_tib()
  
  # A tibble: 301,805 × 10
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
                   Annotation, Chromosome ) %>%
    clean_tib()
  readr::write_csv( x = msa_man_tib, file = out_csv )
  
  msa_man_tib
}

msa_ctls_to_negs = function( x ) {
  
  x %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") ) %>%
    dplyr::mutate(
      # Probe_ID = Probe_ID %>% stringr::str_remove("^ctl_"),
      # Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_Negative_", "neg"),
      Address = U,
      Sequence = AlleleA_ProbeSeq
    ) %>%
    dplyr::select( Address, Probe_ID, Sequence ) %>%
    dplyr::arrange( Probe_ID ) %>%
    dplyr::distinct( Address, .keep_all = TRUE )
}

#
# Need to split this into v03_man_csv controls and merge with v10_man_csv analytical...
#
format_msa_ctls = function( neg_man_csv = "/Users/bbarnes/Documents/data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_SS_BP123_A1.csv.gz",
                            neg_ctl_tib = NULL,
                            
                            out_dir,
                            run_tag,
                            
                            reload     = 0,
                            reload_min = 2,
                            reload_pre = NULL,
                            
                            ret_data   = FALSE,
                            parallel   = FALSE,
                            write_out  = FALSE,
                            
                            vb=0, vt=3, tc=1, tt=NULL,
                            fun_tag='format_msa_ctls')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p4  <- vb > vt + 4
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p4 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid  = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload  = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min  = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre  = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data  = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel  = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     write_out = '{write_out}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_tag = '{out_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, aux_csv, out_csv, end_txt) )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Outdir out='{out_dir}' does not exist")
  if ( !dir.exists( out_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Pre-processing:: Manifest (MSA10)
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # COMMAND LINE:: Quick to build Genome Studio converted manifest sections...
  #   - Plus some historical notes...
  #
  # pwd
  # /Users/bbarnes/Documents/data/pre-idats/MSA
  # gzip -dc MSA-48v1-0-Post-PQC_2B.csv.gz | wc -l                  
  # 296809
  # gzip -dc MSA-48v1-0-Post-PQC_2B.csv.gz | grep -n "^\[Controls\]"
  # 293658:[Controls],,,,,,,,,,,,,,,,,,,,,,,,,,
  #
  # 296809 - 293658 = 3151
  #
  # gzip -dc MSA-48v1-0-Post-PQC_2B.csv.gz | head -n 293657 | tail -n 293650 | gzip -c > MSA-48v1-0-Post-PQC_2B.body.csv.gz
  # gzip -dc MSA-48v1-0-Post-PQC_2B.csv.gz | tail -n 3151 | gzip -c > MSA-48v1-0-Post-PQC_2B.ctls-noheader.gz
  #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Pre-processing:: Most Recent Full Controls
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ctl_man_tib <- NULL
  ctl_man_tib <- readr::read_csv( file = neg_man_csv, show_col_types = FALSE ) %>% 
    clean_tib() %>%
    dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("^ZZ", "ctl_"),
                   Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_neg", "ctl_Negative"),
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
                   SNP_Probe_ID, AlleleA_ProbeSeq ) %>%
    dplyr::filter( !Probe_Type == "cg" ) %>% 
    dplyr::filter( !Probe_Type == "ch" ) %>% 
    dplyr::filter( !Probe_Type == "rs" ) %>%
    dplyr::filter( !Probe_Type == "nv" ) %>%
    dplyr::distinct( U, .keep_all = TRUE )
  
  ctl_man_sum <- NULL
  ctl_man_sum <- print_sum( tib = ctl_man_tib, vec = c("Probe_Type"),
                            vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  ret_tib <- NULL
  ret_tib <- ctl_man_tib
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Write Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( write_out )
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Print Summary::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

format_msa_mans = function( v10_man_csv = "/Users/bbarnes/Documents/data/pre-idats/MSA/MSA-48v1-0-Post-PQC_2B.body.csv.gz",
                            neg_ctl_tib = NULL,
                            
                            out_dir,
                            run_tag,
                            
                            reload     = 0,
                            reload_min = 2,
                            reload_pre = NULL,
                            
                            ret_data   = FALSE,
                            parallel   = FALSE,
                            write_out  = FALSE,
                            
                            vb=0, vt=3, tc=1, tt=NULL,
                            fun_tag='format_msa_mans')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p4  <- vb > vt + 4
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p4 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid  = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload  = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min  = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre  = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data  = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel  = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     write_out = '{write_out}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_tag = '{out_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, aux_csv, out_csv, end_txt) )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Outdir out='{out_dir}' does not exist")
  if ( !dir.exists( out_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Pre-processing:: Manifest (MSA10)
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # COMMAND LINE:: Quick to build Genome Studio converted manifest sections...
  #   - Plus some historical notes...
  #
  # pwd
  # /Users/bbarnes/Documents/data/pre-idats/MSA
  # gzip -dc MSA-48v1-0-Post-PQC_2B.csv.gz | wc -l                  
  # 296809
  # gzip -dc MSA-48v1-0-Post-PQC_2B.csv.gz | grep -n "^\[Controls\]"
  # 293658:[Controls],,,,,,,,,,,,,,,,,,,,,,,,,,
  #
  # 296809 - 293658 = 3151
  #
  # gzip -dc MSA-48v1-0-Post-PQC_2B.csv.gz | head -n 293657 | tail -n 293650 | gzip -c > MSA-48v1-0-Post-PQC_2B.body.csv.gz
  # gzip -dc MSA-48v1-0-Post-PQC_2B.csv.gz | tail -n 3151 | gzip -c > MSA-48v1-0-Post-PQC_2B.ctls-noheader.gz
  #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Pre-processing:: Most Recent Analytical Body
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  v10_man_tib <- NULL
  v10_man_tib <- readr::read_csv( file = v10_man_csv, show_col_types = FALSE ) %>% clean_tib()
  
  v10_man_sum <- NULL
  v10_man_sum <- print_sum( tib = v10_man_tib, vec = c("Probe_Type"),
                            vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Pre-processing:: Rcpp Manifest
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_tib <- NULL
  ret_tib <- # dplyr::bind_rows(
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
    dplyr::select( Probe_ID, U,M, col, Manifest, Annotation, Chromosome,
                   Probe_Type ) %>%
    #   v03_ctl_tib %>% 
    #     dplyr::mutate( 
    #       Annotation = Probe_Type,
    #       Chromosome = "chr0"
    #     ) %>%
    #     dplyr::select( Probe_ID, U,M, col, 
    #                    Species, 
    #                    Manifest, Manifest_Version,
    #                    Annotation, Chromosome )
    # ) %>%
    dplyr::mutate(
      Name = Probe_ID %>% stringr::str_remove("_.*$"),
      Species = "homo_sapiens",
      Manifest = "MSA",
      Manifest_Version = "V1"
    ) %>%
    dplyr::select( Probe_ID, U,M, col, Name,
                   Species, Manifest, Manifest_Version,
                   Annotation, Chromosome ) #, Probe_Type )
  
  # msa1_man_tib %>% dplyr::filter( is.na(col) )
  # msa_man_tib <- NULL
  # msa_man_tib <- dplyr::bind_rows( all_man_tib,msa1_man_tib )
  # 
  # opt$return_df <- 3
  #
  # ann_sum_tib <- readr::read_rds( file = file.path( opt$top_path, "data/manifests/methylation/bgz/all_manifests.sub8.rds" ) ) %>% dplyr::group_by( Annotation ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Write Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( write_out )
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Print Summary::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             Sample Sheet Formatting Functions:: EPICv1/v2 (Alpha)
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_epic_ssh_alpha = function( file = "/Users/bbarnes/Documents/Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/formatted/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz",
                                  
                                  user_format = "8x1",
                                  write_out  = FALSE,
                                  
                                  vb=0, vt=6, tc=1, tt=NULL,
                                  fun_tag='format_epic_ssh_alpha')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p4  <- vb > vt + 4
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p4 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}         path = '{path}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Input File='{file}' does not exist")
  if ( !file.exists( file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_tib <- NULL
  ret_tib <- readr::read_csv( file = file, show_col_types = FALSE ) %>%
    clean_tib() %>%
    dplyr::mutate(
      Sample_Base = Sample_Base %>% stringr::str_to_upper(),
      Sample_Base = dplyr::case_when(
        Sample_Base == "EPIGEN" ~ "EPIGENDX",
        TRUE ~Sample_Base ),
      Sample_ID = Sample_Base %>% stringr::str_sub(1,1),
      Sample_Source = "Alpha",
      Sample_Group = dplyr::case_when(
        Sample_Group == "MeTritration" ~ "MeTitration",
        TRUE ~ Sample_Group ),
      Sample_Class = dplyr::case_when(
        Sample_Group == "CellLine" ~ "r2",
        Sample_Group == "Coriell"  ~ "r2",
        Sample_Group == "MeTitration" ~ "me",
        Sample_Group == "NegativeMouse" ~ "nn",
        TRUE ~ NA_character_ ),
      Source_Key = paste0( "E",Chip_Version ),
      Source_Name = paste( Source_Key,Sheet_Prep,Sheet_Proc, sep="_" ),
      Sample_Input = dplyr::case_when(
        Sample_Group == "MeTitration" ~ 500,
        TRUE ~ Concentration ) %>% as.integer(),
      Sample_Titration = dplyr::case_when(
        Sample_Group == "MeTitration" ~ Concentration,
        TRUE ~ 0.0 ) %>% as.integer(),
      User_Format = user_format
    ) %>% 
    dplyr::select( Source_Key,Source_Name,Sentrix_Name,
                   Sample_Base,Sample_Input,Sample_Titration,
                   Sample_Group,User_Format )
  
  
  ret_sum <- NULL
  ret_sum <- ret_tib %>% 
    print_sum( vec = c("Source_Key","Source_Name","Sample_Base",
                       "Sample_Input","Sample_Titration","User_Format"),
               vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Write Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( write_out )
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Print Summary::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             Sample Sheet Formatting Functions:: MSA.v.1.0
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_msa_sshs = function( paths,
                            
                            out_dir,
                            run_tag,
                            
                            reload     = 0,
                            reload_min = 2,
                            reload_pre = NULL,
                            
                            ret_data   = FALSE,
                            parallel   = FALSE,
                            write_out  = FALSE,
                            
                            vb=0, vt=3, tc=1, tt=NULL,
                            fun_tag='format_msa_sshs')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p4  <- vb > vt + 4
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p4 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid  = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload  = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min  = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre  = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data  = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel  = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     write_out = '{write_out}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_tag = '{out_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, aux_csv, out_csv, end_txt) )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Outdir out='{out_dir}' does not exist")
  if ( !dir.exists( out_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Pre-processing:: Sample Sheet (Replicates v.0)
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ssh_rep0_tib <- NULL
  ssh_rep0_tib <- format_msa_ssh_rep0( path = paths[1], user_format = "48x1",
                                       vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Pre-processing:: Sample Sheet (Titration v.0)
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ssh_uhm0_tib <- NULL
  ssh_uhm0_tib <- format_msa_ssh_uhm0( path = paths[2], user_format = "48x1",
                                       vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Pre-processing:: Sample Sheet (Titration v.1)
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ssh_uhm1_tib <- NULL
  ssh_uhm1_tib <- format_msa_ssh_uhm1( path = paths[3],
                                       vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Pre-processing:: Join All Manifests
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_tib <- NULL
  ret_tib <- ssh_rep0_tib %>%
    dplyr::bind_rows( ssh_uhm0_tib ) %>%
    dplyr::bind_rows( ssh_uhm1_tib ) %>%
    # NOTE: Basic fix for User_Format::
    dplyr::mutate( Sentrix_ID = Sentrix_Name %>% stringr::str_remove("_.*$") ) %>%
    dplyr::add_count( Sentrix_ID, name = "User_Format" ) %>%
    dplyr::mutate( User_Format = paste0(User_Format,"x1") ) %>%
    dplyr::select( -Sentrix_ID )

  ret_sum0 <- ret_tib %>%
    print_sum( vec = c("Source_Key","Source_Name","Sample_Base",
                       "Sample_Input","Sample_Titration","User_Format"),
               vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  ret_sum1 <- ret_tib %>%
    print_sum( vec = c("Sample_Base",
                       "Sample_Input","Sample_Titration","User_Format"),
               vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Write Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( write_out )
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Print Summary::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

format_msa_ssh_rep0 = function( path = "/Users/bbarnes/Documents/data/pre-idats/MSA/sample_sheets/rep0",
                                
                                user_format = NA_character_, # "48x1",
                                write_out  = FALSE,
                                
                                vb=0, vt=6, tc=1, tt=NULL,
                                fun_tag='format_msa_ssh_rep0')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p4  <- vb > vt + 4
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p4 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}         path = '{path}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Input Path='{path}' does not exist")
  if ( !dir.exists( path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ssh_dat <- NULL
  ssh_dat <- file_list( path = path, 
                        prefix  = path,
                        suffix  = "_sampleSheet.csv.gz",
                        pattern = "_sampleSheet.csv.gz$", 
                        recursive = TRUE ) %>% 
    lapply( readr::read_csv, skip=10, show_col_types=FALSE  ) %>%
    dplyr::bind_rows( .id = "Source_Name" ) %>%
    dplyr::rename( Sentrix_Name = Sample_ID ) %>% 
    dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
    dplyr::select( Source_Name,Sentrix_Name, dplyr::everything() )
  
  ret_tib <- NULL
  ret_tib <- ssh_dat %>%
    set_names( ssh_dat %>% names() %>% stringr::str_replace_all(" ","_") ) %>%
    dplyr::mutate( 
      Sample_Base = dplyr::case_when(
        !is.na(Sample_Well) ~ Sample_Well, 
        TRUE ~ NA_character_ ) %>% 
        stringr::str_replace_all("-","") %>% 
        stringr::str_to_upper(),
      Sample_Fix = Sample_Name %>% 
        stringr::str_remove("^[^-]+-") %>% 
        # stringr::str_remove("^[^\\s]+\\\\s") %>%
        stringr::str_remove_all("[A-Z]") %>%
        stringr::str_remove_all("[a-z]") %>%
        stringr::str_remove_all( " " ) %>%
        stringr::str_remove( "ng$" ),
      Sample_Input = dplyr::case_when( 
        !is.na( Sample_Group ) ~ Sample_Group,
        is.na(Sample_Group) ~ Sample_Fix, 
        TRUE ~ NA_character_ ) %>%
        stringr::str_remove_all( " " ) %>%
        stringr::str_remove( "ng$" ) %>% as.integer(),
      Sample_Titration = 0.0 %>% as.integer(),
      User_Format = user_format,
      Sample_Group = dplyr::case_when( 
        Sample_Base == "CORIELL"  ~ "Coriell", 
        Sample_Base == "EPIGENDX" ~ "MeTitration",
        Sample_Base == "HELA"     ~ "CellLine",
        Sample_Base == "JURKAT"   ~ "CellLine",
        Sample_Base == "MCF7"     ~ "CellLine",
        Sample_Base == "K562"     ~ "CellLine",
        Sample_Base == "RAJI"     ~ "CellLine",
        Sample_Base == "HELA"     ~ "CellLine",
        TRUE ~ NA_character_ ),
      Source_Key = "R0"
    ) %>%
    dplyr::select( Source_Key,Source_Name,Sentrix_Name,
                   Sample_Base,Sample_Input,Sample_Titration,
                   Sample_Group,User_Format )
  
  ret_sum <- NULL
  ret_sum <- ret_tib %>% 
    print_sum( vec = c("Source_Key","Source_Name","Sample_Base",
                       "Sample_Input","Sample_Titration","User_Format"),
               vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Write Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( write_out )
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Print Summary::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

format_msa_ssh_uhm0 = function( path = "/Users/bbarnes/Documents/data/pre-idats/MSA/sample_sheets/uhm0",
                                
                                user_format = NA_character_, # "48x1",
                                write_out  = FALSE,
                                
                                vb=0, vt=6, tc=1, tt=NULL,
                                fun_tag='format_msa_ssh_uhm0')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p4  <- vb > vt + 4
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p4 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}         path = '{path}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Input Path='{path}' does not exist")
  if ( !dir.exists( path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ssh_dat <- NULL
  ssh_dat <- file_list( path    = path, 
                        prefix  = path,
                        suffix  = "_sampleSheet.csv.gz",
                        pattern = "_sampleSheet.csv.gz$", 
                        recursive = TRUE ) %>% 
    lapply( readr::read_csv, skip=10, show_col_types=FALSE  ) %>%
    dplyr::bind_rows( .id = "Source_Name" ) %>%
    dplyr::rename( Sentrix_Name = Sample_ID ) %>% 
    dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
    dplyr::select( Source_Name,Sentrix_Name, dplyr::everything() )
  
  # return(ssh_dat)
  # head(ssh_dat)
  # tail(ssh_dat)
  
  ret_tib <- NULL
  ret_tib <- ssh_dat %>%
    set_names( ssh_dat %>% names() %>% stringr::str_replace_all(" ","_") ) %>%
    dplyr::mutate( 
      Sample_Base = dplyr::case_when(
        !is.na(Sample_Well) ~ Sample_Well, 
        TRUE ~ NA_character_ ) %>% 
        stringr::str_to_upper(),
      Sample_Fix = Sample_Name %>% 
        stringr::str_remove("^[^-]+-") %>% 
        stringr::str_remove_all("[A-Z]") %>%
        stringr::str_remove_all("[a-z]") %>%
        stringr::str_remove_all( " " ) %>%
        stringr::str_remove( "ng$" ),
      
      Sample_Input = Sample_Name %>% 
        stringr::str_remove("^.*-") %>% 
        stringr::str_remove(" ng$") %>% 
        stringr::str_remove( "^.* " ) %>% as.integer(),
      
      Sample_Titration = Sample_Plate %>%
        stringr::str_remove_all( " " ) %>%
        stringr::str_remove( "%" ) %>%
        stringr::str_remove( "ng$" ),
      Sample_Titration = dplyr::case_when(
        is.na(Sample_Titration) ~ "0",
        TRUE ~ Sample_Titration ) %>% as.integer(),
      User_Format = user_format,
      Sample_Group = dplyr::case_when(
        Sample_Base == "RAJI"     ~ "CellLine",
        Sample_Base == "CORIELL"  ~ "Coriell",
        Sample_Base == "EPIGENDX" ~ "MeTitration",
        TRUE ~ NA_character_ ),
      Source_Key = "T0"
    ) %>%
    dplyr::select( Source_Key,Source_Name,Sentrix_Name,
                   Sample_Base,Sample_Input,Sample_Titration,
                   Sample_Group,User_Format )
  
  ret_sum <- NULL
  ret_sum <- ret_tib %>% 
    print_sum( vec = c("Source_Key","Source_Name","Sample_Base",
                       "Sample_Input","Sample_Titration","User_Format"),
               vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Write Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( write_out )
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Print Summary::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

format_msa_ssh_uhm1 = function( path = "/Users/bbarnes/Documents/data/pre-idats/MSA/sample_sheets/uhm1",
                                
                                user_format = NA_character_, # "48x1",
                                write_out  = FALSE,
                                
                                vb=0, vt=6, tc=1, tt=NULL,
                                fun_tag='format_msa_ssh_uhm1')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p4  <- vb > vt + 4
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p4 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}         path = '{path}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Input Path='{path}' does not exist")
  if ( !dir.exists( path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ssh_dat <- NULL
  ssh_dat <- file_list( path = path, 
                        prefix  = path,
                        suffix  = "_sampleSheet.csv.gz",
                        pattern = "_sampleSheet.csv.gz$", 
                        recursive = TRUE ) %>% 
    lapply( readr::read_csv, skip=10, show_col_types=FALSE  ) %>%
    dplyr::bind_rows( .id = "Source_Name" ) %>%
    dplyr::rename( Sentrix_Name = Sample_ID ) %>% 
    dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
    # dplyr::mutate( Sentrix_ID = Sentrix_Name %>% stringr::str_remove("_.*$") ) %>%
    # dplyr::add_count( Sentrix_ID, name = "User_Format" ) %>%
    dplyr::select( Source_Name,Sentrix_Name, dplyr::everything() )
  
  # return( ssh_dat )
  ret_tib <- NULL
  ret_tib <- ssh_dat %>%
    set_names( ssh_dat %>% names() %>% stringr::str_replace_all(" ","_") ) %>%
    dplyr::mutate( 
      Sample_Base = dplyr::case_when(
        Sample %>% stringr::str_detect("Epigen") ~ "EpigenDX",
        TRUE ~ Sample ) %>% 
        stringr::str_to_upper(),
      Sample_Input = DNA_Input %>%
        stringr::str_remove( "ng$" ) %>% as.integer(),
      Sample_Titration = dplyr::case_when(
        Sample_Base == "EPIGENDX" & DNA %>% stringr::str_starts("0%")   ~ 0.0,
        Sample_Base == "EPIGENDX" & DNA %>% stringr::str_starts("50%")  ~ 50.0,
        Sample_Base == "EPIGENDX" & DNA %>% stringr::str_starts("100%") ~ 100.0,
        TRUE ~ 0.0 ) %>% as.integer(),
      #  TRUE ~ NA_real_ ) %>% as.integer(),
      User_Format = user_format,
      Sample_Group = Sample_Type %>% stringr::str_remove_all(" "),
      Sample_Group = dplyr::case_when(
        Sample_Base == "R311874" ~ "Blood",
        Sample_Base == "R315156" ~ "Blood",
        Sample_Group == "Control" ~ "MeTitration",
        TRUE ~ Sample_Group ),
      Sample_Base = dplyr::case_when(
        Sample_Group == "FFPE" ~ paste0( Sample_Group,Sample_Base ),
        TRUE ~ Sample_Base ),
      Source_Key = "T1"
    ) %>%
    dplyr::select( Source_Key,Source_Name,Sentrix_Name,
                   Sample_Base,Sample_Input,Sample_Titration,
                   Sample_Group,User_Format )
  
  ret_sum <- NULL
  ret_sum <- ret_tib %>% 
    print_sum( vec = c("Source_Key","Source_Name","Sample_Base",
                       "Sample_Input","Sample_Titration","User_Format"),
               vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Write Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( write_out )
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Print Summary::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}


# End of file

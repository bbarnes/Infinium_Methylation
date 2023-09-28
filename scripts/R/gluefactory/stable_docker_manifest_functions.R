
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       EPICv2 Docker Manifest Functions::
#                                 VA/MVP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Command Line Options Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("optparse",   quietly = TRUE) ) )

# Tidyverse Core Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("tidyverse",  quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("stringr", quietly = TRUE) ))
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("readr",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("glue",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("furrr",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("purrr",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("magrittr", quietly = TRUE) ) )

# Parallel Processing Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("doParallel", quietly = TRUE) ) )

# Plotting Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("ggplot2", quietly = TRUE) ) )

# Matrix Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("matrixStats", quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("scales",      quietly = TRUE) ))

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
#                       Fingerprinting (SNV) Fuctions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

vcf_cols <- NULL
vcf_cols <- readr::cols(
  Chromosome = readr::col_character(),
  Coordinate = readr::col_integer(),
  SNP_ID     = readr::col_character(),
  REF        = readr::col_character(),
  ALT        = readr::col_character(),
  Qual       = readr::col_character(),
  Filter     = readr::col_character(),
  Info_Str   = readr::col_character()
)

fingerprint_vcfs = function( vcfs,
                             sam_tib,
                             run_sig = FALSE,
                             
                             out_dir,
                             run_tag,
                             
                             reload     = 0,
                             reload_min = 2,
                             reload_pre = NULL,
                             
                             ret_data   = FALSE,
                             parallel   = FALSE,
                             write_out  = FALSE,
                             write_sum  = FALSE,
                             write_sig  = FALSE,
                             
                             vb=0, vt=3, tc=1, tt=NULL,
                             fun_tag='fingerprint_vcfs')
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
  if ( p1 ) {
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
    cat(glue::glue("{mssg}     write_sum = '{write_sum}'.{RET}"))
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
  
  errs_mssg <- glue::glue("VCF List is emplty")
  if ( is.null( vcfs) || length(vcfs) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  snp_tab <- NULL
  snp_tab <- vcfs %>% 
    lapply( readr::read_tsv, 
            col_names=names(vcf_cols$cols), 
            col_types=vcf_cols, skip=6 ) %>% 
    dplyr::bind_rows( .id = "Sentrix_Name" ) %>%
    tidyr::separate( Info_Str, into = c("PVF","GT","GS"), sep=";", remove = TRUE ) %>%
    tidyr::unite( Target_ID, SNP_ID,Chromosome,Coordinate,REF,ALT, 
                  sep=":", remove = TRUE) %>%
    dplyr::mutate(
      Filter = dplyr::case_when(
        Filter == "PASS" ~ 1.0,
        Filter == "FAIL" ~ 0.0,
        TRUE ~ NA_real_ ) %>% as.integer(),
      PVF = PVF %>% stringr::str_remove("^PVF=") %>% as.double(),
      PVF = as.integer( PVF * 1000 ),
      GT  = GT  %>% stringr::str_remove( "^GT="),
      GTS = GS  %>% stringr::str_remove( "^GS=") %>% as.integer(),
      RGT = GT  %>% stringr::str_remove("/.*$") %>% as.integer(),
      AGT = GT  %>% stringr::str_remove("^.*/") %>% as.integer(),
      GTC = dplyr::case_when(
        RGT==0 & AGT==0 ~ 0.0,
        RGT==0 & AGT==1 ~ 1.0,
        RGT==1 & AGT==0 ~ 2.0,
        RGT==1 & AGT==1 ~ 3.0,
        TRUE ~ NA_real_ ) %>% as.integer()
    ) %>% dplyr::select( -GT ) %>% clean_tib()
  
  snp_sum <- NULL
  snp_sum <- snp_tab %>% 
    dplyr::group_by( RGT,AGT ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  print_tib( snp_tab, fun_tag = fun_tag, name = "SNP-Summary", 
             vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  ret_tib <- snp_tab
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Format SNP Table to Matrix List::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  snp_tab_lst <- NULL
  snp_tab_lst <- snp_tab %>%
    tidyr::pivot_longer( cols = c( Qual,Filter,PVF,GTS,RGT,AGT,GTC ), 
                         names_to  = c("Key"), 
                         values_to = c("Val") ) %>% split(.$Key)
  
  snp_mat_lst <- NULL
  snp_mat_lst <- snp_tab_lst %>% 
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
  
  vet_gtc_mat <- NULL
  vet_gtc_mat <- snp_mat_lst$GTC
  vet_agt_mat <- NULL
  vet_agt_mat <- snp_mat_lst$AGT
  vet_agt_mat <- NULL
  vet_agt_mat <- snp_mat_lst$RGT
  
  # vet_tib_list <- NULL
  # vet_tib_list <- snp_tab_lst %>% 
  #   lapply( function(x) { 
  #     x %>% dplyr::select(-Key) %>% 
  #       tidyr::pivot_wider( id_cols = c(Target_ID), 
  #                           names_from = c(Sentrix_Name), 
  #                           values_from = c(Val) )
  #   })
  # 
  # vet_gtc_mat <- NULL
  # vet_gtc_mat <- vet_tib_list$GTC %>% 
  #   tibble::column_to_rownames( var = "Target_ID") %>% 
  #   as.data.frame() %>% as.matrix()
  # 
  # vet_agt_mat <- NULL
  # vet_agt_mat <- vet_tib_list$AGT %>% 
  #   tibble::column_to_rownames( var = "Target_ID") %>% 
  #   as.data.frame() %>% as.matrix()
  # 
  # vet_agt_mat <- NULL
  # vet_agt_mat <- vet_tib_list$RGT %>% 
  #   tibble::column_to_rownames( var = "Target_ID") %>% 
  #   as.data.frame() %>% as.matrix()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Format SNP Martix Lists to Matrix::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  col_vec <- c(1:base::ncol(vet_gtc_mat) )
  row_vec <- c(1:base::nrow(vet_gtc_mat) )
  
  sam_mat <- NULL
  sam_mat <- sam_tib %>% 
    dplyr::select( Sentrix_Name,Sample_Base,cg_calls_pass_perc_1 ) %>% 
    tibble::column_to_rownames( var = "Sentrix_Name" ) %>% 
    as.data.frame() %>% as.matrix()
  
  # Expectation = [ A tibble: 1,600 × 11 ] - All Samples
  # Expectation = [ A tibble:   576 × 11 ] - Non-plot_sample_vec Samples
  stats_tab <- NULL
  for ( ii in col_vec ) {
    if ( ! (sam_mat[ colnames(vet_gtc_mat)[ii], 1 ] %in% plot_sample_vec) ) next
    
    for ( jj in col_vec ) {
      if ( ! (sam_mat[ colnames(vet_gtc_mat)[jj], 1 ] %in% plot_sample_vec) ) next
      
      tot_cnt <- vet_gtc_mat[ ,c(ii,jj) ] %>% base::nrow()
      mat_cnt <- which( vet_gtc_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 ) %>% length()
      mis_cnt <- which( vet_gtc_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() != 0 ) %>% length()
      agt_cnt <- which( vet_agt_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 ) %>% length()
      rgt_cnt <- which( vet_agt_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 ) %>% length()
      
      stats_tab <- stats_tab %>%
        dplyr::bind_rows(
          tibble::tibble(
            Sentrix_NameA = colnames(vet_gtc_mat)[ii],
            Sentrix_NameB = colnames(vet_gtc_mat)[jj],
            Sample_BaseA = sam_mat[ Sentrix_NameA, ][1],
            Sample_BaseB = sam_mat[ Sentrix_NameB, ][1],
            Idx = ii,
            Jdx = jj,
            Tot = tot_cnt,
            Mat = mat_cnt / tot_cnt,
            Mis = mis_cnt / tot_cnt,
            Agt = agt_cnt / tot_cnt,
            Rgt = rgt_cnt / tot_cnt )
        )
      
      # fig_gtc_tib %>% dplyr::filter( Sample_Base %in% plot_sample_vec )
      
      # cat(glue::glue("ii={ii}, jj={jj} = {tot_cnt}, {mat_cnt}, {mis_cnt}, {agt_cnt}, {rgt_cnt}.{RET}") )
    }
  }

  stats_sum_csv <- NULL
  stats_sum_csv <- file.path( out_dir, paste0( run_tag,".fingerprint-signatures.stats.csv.gz" ) )
  stats_sum <- NULL
  stats_sum <- stats_tab %>% 
    # dplyr::filter( Sample_BaseA == Sample_BaseB ) %>% 
    dplyr::filter( Sentrix_NameA != Sentrix_NameB ) %>% 
    dplyr::mutate(
      Sample_Map = dplyr::case_when(
        Sentrix_NameA == Sentrix_NameB ~ "ID",
        Sample_BaseA == Sample_BaseB ~ "==",
        TRUE ~ "!=" )
    ) %>%
    dplyr::group_by( Sample_Map,Sample_BaseA,Sample_BaseB ) %>% 
    dplyr::summarise( Mat_Sds = sd( Mat, na.rm = TRUE ),
                      Agt_Sds = sd( Agt, na.rm = TRUE ),
                      Rgt_Sds = sd( Rgt, na.rm = TRUE ),
                      
                      Mat_Mad = mad( Mat, na.rm = TRUE ),
                      Agt_Mad = mad( Agt, na.rm = TRUE ),
                      Rgt_Mad = mad( Rgt, na.rm = TRUE ),
                      
                      Mat_Avg = mean( Mat, na.rm=TRUE ),
                      Agt_Avg = mean( Agt, na.rm=TRUE ),
                      Rgt_Avg = mean( Rgt, na.rm=TRUE ),
                      
                      Mat_Med = median( Mat, na.rm=TRUE ),
                      Agt_Med = median( Agt, na.rm=TRUE ),
                      Rgt_Med = median( Rgt, na.rm=TRUE ),
                      
                      .groups = "drop" )
  if ( write_sum ) readr::write_csv( x = stats_sum, file = stats_sum_csv )

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Build Binary Signature::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( run_sig ) {
    
    fig_gtc_tib <- NULL
    fig_gtc_tib <- vet_gtc_mat %>% t() %>%
      as.data.frame() %>% 
      tibble::rownames_to_column( var = "Sentrix_Name" ) %>% 
      tibble::as_tibble() %>%
      # head() %>%
      tidyr::unite( FingerPrint, 2:base::ncol(.), sep="", remove = TRUE ) %>%
      dplyr::left_join( sam_tib, by=c("Sentrix_Name") ) %>%
      dplyr::filter( Sample_Base %in% plot_sample_vec )
    
    fig_gtc_lst <- NULL
    fig_gtc_lst <- fig_gtc_tib %>% 
      dplyr::rename( PPP=cg_calls_pass_perc_1 ) %>%
      dplyr::select( Sample_Base,PPP,FingerPrint) %>% 
      dplyr::arrange( Sample_Base,PPP ) %>% 
      split(.$Sample_Base)
    
    # Write Binary Signature...
    #
    if ( write_sig ) {
      fig_gtc_csv <- NULL
      fig_gtc_csv <- file.path( out_dir, paste0( run_tag,".fingerprint-signatures.csv" ) )
      if ( fig_gtc_csv %>% file.exists() ) unlink( x = fig_gtc_csv )
      for ( sample in names(fig_gtc_lst) ) {
        if ( p0 ) cat(glue::glue("{pmssg} Sample={sample}.{RET}"))
        readr::write_csv( x = fig_gtc_lst[[sample]], file = fig_gtc_csv, append = TRUE )
      }
    }
  }

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
  
  if ( ret_data ) {
    
    ret_dat$snp_tab <- snp_tab
    ret_dat$sam_mat <- sam_mat
    
    ret_dat$snp_tab_lst <- snp_tab_lst
    ret_dat$snp_mat_lst <- snp_mat_lst

    ret_dat$stats_tab <- stats_tab
    ret_dat$stats_sum <- stats_sum
    
    ret_dat$fig_gtc_tib <- fig_gtc_tib
    ret_dat$fig_gtc_lst <- fig_gtc_lst

    return( ret_dat )
  }
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       AnnoI/AnnoS SNP Fuctions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

anno_to_grs = function( tib,
                        add_chr = FALSE,
                        isAnnoS = TRUE,
                        
                        chr_key = "Chromosome",
                        pos_key = "Coordinate_SNP",
                        ids_key = "Probe_ID",
                        snp_key = "SNP_ID",
                        
                        vb=0, vt=6, tc=1, tt=NULL,
                        fun_tag='anno_to_grs')
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      chr_key = '{chr_key}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  chr_sym <- rlang::sym( chr_key )
  pos_sym <- rlang::sym( pos_key )
  ids_sym <- rlang::sym( ids_key )
  snp_sym <- rlang::sym( snp_key )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Ensure its sorted...
  tib <- tib %>% dplyr::arrange( !!chr_sym, !!pos_sym, !!ids_sym, !!snp_sym )
  
  if ( add_chr ) tib <- tib %>%
    dplyr::mutate( !!chr_sym := paste0("chr",!!chr_sym) )
  
  if ( isAnnoS ) {
    tib <- tib %>%
      dplyr::mutate(
        U = dplyr::case_when(
          REF_SNP == "C" ~ "ALT",
          REF_SNP == "G" ~ "ALT",
          REF_SNP == "A" ~ "REF",
          REF_SNP == "T" ~ "REF",
          TRUE ~ NA_character_ )
      ) %>% dplyr::filter( !is.na(U) ) %>% 
      dplyr::distinct( !!ids_sym, !!snp_sym, .keep_all = TRUE )
    
    ret_dat <- GenomicRanges::GRanges(
      seqnames = Rle( tib %>% dplyr::pull( chr_key ) ),
      strand = Rle( tib$strand ),
      
      rs         = tib$SNP_ID,
      designType = tib$Infinium_Design_Type,
      U          = tib$U,
      
      REF = tib$REF_SNP,
      ALT = tib$ALT,
      
      IRanges( start = tib$Coordinate_SNP,
               end   = tib$Coordinate_SNP,
               names = tib %>% dplyr::pull( ids_key ) )
    )
  } else {
    ret_dat <- GenomicRanges::GRanges(
      seqnames = Rle( tib %>% dplyr::pull( chr_key ) ),
      strand = Rle( tib$strand ),
      
      designType = tib$Infinium_Design_Type,
      In.band    = "REF",
      
      REF = tib$REF_SNP,
      ALT = tib$ALT,
      rs  = tib$SNP_ID,
      
      IRanges( start = tib$Coordinate_SNP,
               end   = tib$Coordinate_SNP,
               names = tib %>% dplyr::pull( ids_key ) )
    )
  }
  
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

snps_to_annoS = function( tib,
                          vb=0, vt=6, tc=1, tt=NULL,
                          fun_tag='snps_to_annoS')
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_tib <- tib %>%
    # dplyr::filter( Infinium_Design_Type == "I" ) %>%
    dplyr::filter( Probe_Type == "rs" | Probe_Type == "nv" ) %>%
    # dplyr::filter( SNP_Match_Type == "Cpg_Up" | SNP_Match_Type == "Cpg_Dn" ) %>%
    dplyr::filter( 
      (REF_SNP=="A" & ALT=="G") |
        (REF_SNP=="A" & ALT=="G") |
        (REF_SNP=="C" & ALT=="T") |
        (REF_SNP=="G" & ALT=="A") |
        (REF_SNP=="T" & ALT=="C")
    ) %>% clean_tib()
  
  ret_sum <- NULL
  ret_sum <- ret_tib %>%
    dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, 
                     REF_SNP,ALT ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p0 ) print( ret_sum, n=base::nrow(ret_sum) )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

snps_to_annoI = function( tib,
                          vb=0, vt=6, tc=1, tt=NULL,
                          fun_tag='snps_to_annoI')
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_tib <- tib %>%
    dplyr::filter( Infinium_Design_Type == "I" ) %>%
    dplyr::filter( Probe_Type == "cg" | Probe_Type == "ch" ) %>%
    # dplyr::filter( SNP_Match_Type == "Nxb_Dn" | SNP_Match_Type == "Nxb_Up" ) %>%
    dplyr::filter( 
      (SNP_Match_Type == "Nxb_Dn" & REF_SNP=="A" & ALT=="C") |
        (SNP_Match_Type == "Nxb_Dn" & REF_SNP=="A" & ALT=="C") |
        (SNP_Match_Type == "Nxb_Dn" & REF_SNP=="C" & ALT=="A") |
        (SNP_Match_Type == "Nxb_Dn" & REF_SNP=="C" & ALT=="G") |
        (SNP_Match_Type == "Nxb_Dn" & REF_SNP=="C" & ALT=="T") |
        (SNP_Match_Type == "Nxb_Dn" & REF_SNP=="G" & ALT=="C") |
        (SNP_Match_Type == "Nxb_Dn" & REF_SNP=="T" & ALT=="C") |
        
        (SNP_Match_Type == "Nxb_Up" & REF_SNP=="A" & ALT=="G") |
        (SNP_Match_Type == "Nxb_Up" & REF_SNP=="C" & ALT=="G") |
        (SNP_Match_Type == "Nxb_Up" & REF_SNP=="G" & ALT=="A") |
        (SNP_Match_Type == "Nxb_Up" & REF_SNP=="G" & ALT=="C") |
        (SNP_Match_Type == "Nxb_Up" & REF_SNP=="G" & ALT=="T") |
        (SNP_Match_Type == "Nxb_Up" & REF_SNP=="T" & ALT=="G")
    ) %>% clean_tib()
  
  ret_sum <- NULL
  ret_sum <- ret_tib %>%
    dplyr::group_by( Probe_Type,Infinium_Design_Type,SNP_Match_Type, 
                     REF_SNP,ALT ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p0 ) print( ret_sum, n=base::nrow(ret_sum) )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       EPICv2 Docker Manifest Fuctions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rank_historic_probes = function( tib,
                                 prb_tibs = NULL,
                                 
                                 add_build = FALSE,
                                 man_col_vec = c( 
                                   "Probe_ID","Full_ID","Loci_ID","Probe_Type",
                                   "AddressA_ID","AddressB_ID",
                                   "AlleleA_ProbeSeq","AlleleB_ProbeSeq",
                                   "Strand_FR","Strand_TB","Strand_CO","Next_Base",
                                   "Infinium_Design_Type","Infinium_Design","Rep_Num",
                                   # "Chromosome_hg38","Coordinate_hg38","strand",
                                   "Canonical_Rank","History_Cnt","Loci_Cnt" ),
                                 
                                 out_dir,
                                 run_tag,
                                 
                                 reload     = 0,
                                 reload_min = 2,
                                 reload_pre = NULL,
                                 
                                 ret_data   = FALSE,
                                 parallel   = FALSE,
                                 write_out  = FALSE,
                                 
                                 vb=0, vt=3, tc=1, tt=NULL,
                                 fun_tag='rank_historic_probes')
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
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}    add_build = '{add_build}'.{RET}"))
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
  
  errs_mssg <- glue::glue("Input manifest (tib) is empty")
  if ( is.null(tib) || base::nrow(tib) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return( NULL )
  
  errs_mssg <- glue::glue("Probes List (prb_tibs) is empty")
  if ( is.null(prb_tibs) || length(prb_tibs) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return( NULL )
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Calculate Historic Counts::
  #
  hist_tib <- NULL
  hist_tib <- lapply( prb_tibs, function(x) { 
    dplyr::select( x, AlleleA_ProbeSeq,AlleleB_ProbeSeq ) } ) %>%
    dplyr::bind_rows() %>% 
    dplyr::arrange( AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>% 
    dplyr::group_by( AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>% 
    dplyr::summarise( History_Cnt=n(), .groups = "drop" )
  
  cnts_tib <- NULL
  cnts_tib <- tib %>% 
    dplyr::rename(
      Probe_ID = IlmnID,
      Loci_ID  = Name
    ) %>%
    dplyr::add_count( Loci_ID, name="Loci_Cnt" ) %>%
    dplyr::left_join( hist_tib, by=c("AlleleA_ProbeSeq","AlleleB_ProbeSeq") ) %>%
    dplyr::mutate(
      History_Cnt = dplyr::case_when(
        is.na(History_Cnt) ~ 0,
        TRUE ~ History_Cnt
      ) %>% as.integer()
    ) %>%
    dplyr::arrange( -History_Cnt,Loci_Cnt,Rep_Num )
  
  rank_tib <- NULL
  rank_tib <- cnts_tib %>%
    dplyr::group_by( Loci_ID ) %>%
    dplyr::mutate( Canonical_Rank = dplyr::row_number() ) %>%
    dplyr::ungroup()
  
  # Update man_col_vec to include new fields::
  #
  chr_key <- "Chromosome"
  pos_key <- "Coordinate"
  if ( add_build ) {
    chr_key <- paste( "Chromosome",tib$Genome_Build %>% head(n=1), sep="_" )
    pos_key <- paste( "Coordinate",tib$Genome_Build %>% head(n=1), sep="_" )
  }
  chr_sym <- rlang::sym( chr_key )
  pos_sym <- rlang::sym( pos_key )
  
  man_col_vec <- c(
    "Probe_ID","Full_ID","Loci_ID",man_col_vec,
    chr_key,pos_key,"strand",
    "Nxb_Nuc_Up","Cpg_Nuc_Up","Cpg_Nuc_Dn","Nxb_Nuc_Dn",
    "Nxb_Pos_Up","Cpg_Pos_Up","Cpg_Pos_Dn","Nxb_Pos_Dn"
    # "Nxb_Pos_U2","Nxb_Pos_U1","Nxb_Pos_D1","Nxb_Pos_D2"
  ) %>% unique()
  
  ret_tib <- NULL
  ret_tib <- rank_tib %>%
    dplyr::mutate(
      Full_ID  = Probe_ID,
      Loci_ID  = Full_ID %>% stringr::str_remove("_.*$"),
      Probe_ID = dplyr::case_when(
        Canonical_Rank == 1 ~ Loci_ID,
        TRUE ~ Probe_ID
      ),
      strand = dplyr::case_when(
        Strand_FR == "F" ~ "+",
        Strand_FR == "R" ~ "-",
        TRUE ~ "*" ),
      
      # Update Nucleotides: Next_Base (Nxb) and Target (Cpg) Info::
      Fwd_Seq = Forward_Sequence %>% stringr::str_remove("/"),
      Nxb_Nuc_Up = Fwd_Seq %>% stringr::str_sub(60,60),
      Cpg_Nuc_Up = Fwd_Seq %>% stringr::str_sub(62,62),
      Cpg_Nuc_Dn = Fwd_Seq %>% stringr::str_sub(63,63),
      Nxb_Nuc_Dn = Fwd_Seq %>% stringr::str_sub(65,65),
      
      # Update Coordinates: Next_Base (Nxb) and Target (Cpg) Info::
      # Nxb_Pos_U2 = MAPINFO - 3,
      # Nxb_Pos_U1 = MAPINFO - 2,
      Nxb_Pos_Up = MAPINFO - 1,
      Cpg_Pos_Up = MAPINFO - 0,
      Cpg_Pos_Dn = dplyr::case_when(
        Probe_Type == "cg" ~ MAPINFO + 1,
        Probe_Type == "ch" ~ MAPINFO + 1,
        Probe_Type == "nv" ~ MAPINFO + 0,
        Probe_Type == "rs" ~ MAPINFO + 0,
        TRUE ~ NA_real_ ) %>% as.integer(),
      Nxb_Pos_Dn = Cpg_Pos_Dn + 1,
      # Nxb_Pos_D1 = Cpg_Pos_Dn + 2,
      # Nxb_Pos_D2 = Cpg_Pos_Dn + 3
    ) %>%
    dplyr::rename(
      !!chr_sym := CHR,
      !!pos_sym := MAPINFO
    )  %>% 
    dplyr::select( dplyr::all_of(man_col_vec) ) %>%
    dplyr::arrange( !!chr_sym,!!pos_sym,Canonical_Rank,
                    Loci_Cnt,Rep_Num )
  
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
  
  #
  # [Done]: Calculate Stats
  #
  ret_sum <- NULL
  ret_sum <- ret_tib %>% 
    dplyr::group_by( Canonical_Rank, History_Cnt ) %>%
    dplyr::summarise( 
      Total = n(),
      Min = min(Loci_Cnt, na.rm = TRUE),
      Avg = mean(Loci_Cnt, na.rm = TRUE),
      Med = median(Loci_Cnt, na.rm = TRUE),
      Sds = sd(Loci_Cnt, na.rm = TRUE),
      Mad = mad(Loci_Cnt, na.rm = TRUE),
      Max = max(Loci_Cnt, na.rm = TRUE),
      .groups = "drop"
    )
  if ( p2 ) ret_sum %>% print( n=base::nrow(ret_sum) )
  
  # Basically Proof that most Probes are Canonical...
  # [Done]: Print baic summary stats::
  tot_id_cnt <- ret_tib %>% base::nrow()
  add_id_cnt <- ret_tib %>% dplyr::distinct( AddressA_ID,AddressB_ID ) %>% base::nrow()
  seq_id_cnt <- ret_tib %>% dplyr::distinct( AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>% base::nrow()
  prb_id_cnt <- ret_tib %>% dplyr::distinct( Probe_ID ) %>% base::nrow()
  ful_id_cnt <- ret_tib %>% dplyr::distinct( Full_ID ) %>% base::nrow()
  loc_id_cnt <- ret_tib %>% dplyr::distinct( Loci_ID ) %>% base::nrow()
  
  if ( p2 ) cat(glue::glue("{pmssg} tot_id_cnt = '{tot_id_cnt}'{RET}"))
  if ( p2 ) cat(glue::glue("{pmssg} add_id_cnt = '{add_id_cnt}'{RET}"))
  if ( p2 ) cat(glue::glue("{pmssg} seq_id_cnt = '{seq_id_cnt}'{RET}"))
  if ( p2 ) cat(glue::glue("{pmssg} prb_id_cnt = '{prb_id_cnt}'{RET}"))
  if ( p2 ) cat(glue::glue("{pmssg} ful_id_cnt = '{ful_id_cnt}'{RET}"))
  if ( p2 ) cat(glue::glue("{pmssg} loc_id_cnt = '{loc_id_cnt}'{RET}"))
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}


# End of file

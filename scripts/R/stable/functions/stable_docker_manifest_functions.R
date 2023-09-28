
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
#                       Fingerprinting (SNV) Functions::
#                                  BASIC
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

read_snv_vcf = function( x ) {
  readr::read_tsv( file = x, 
                   col_names=names(vcf_cols$cols), 
                   col_types=vcf_cols, skip=6 )
}

parse_snv_vcf = function( x, gts_min = 0, fval = 5.0, jcnt = 0, uval = TRUE )
{
  x <- x %>% 
    tidyr::separate( Info_Str, into = c("PVF","GT","GS"), sep=";", remove = TRUE )
  
  if ( jcnt == 4 ) x <- x %>% 
      tidyr::unite( Target_ID, SNP_ID,Chromosome,Coordinate,REF,ALT, sep=":", remove = TRUE)
  if ( jcnt == 3 ) x <- x %>% 
      tidyr::unite( Target_ID, SNP_ID,Chromosome,Coordinate,REF, sep=":", remove = TRUE) %>%
      dplyr::select( -ALT )
  if ( jcnt == 2 ) x <- x %>% 
      tidyr::unite( Target_ID, SNP_ID,Chromosome,Coordinate, sep=":", remove = TRUE) %>%
      dplyr::select( -REF,-ALT )
  if ( jcnt == 1 ) x <- x %>% 
      tidyr::unite( Target_ID, SNP_ID,Chromosome, sep=":", remove = TRUE) %>%
      dplyr::select( -Coordinate,-REF,-ALT )
  if ( jcnt == 0 ) x <- x %>% 
      tidyr::unite( Target_ID, SNP_ID, sep=":", remove = TRUE) %>%
      dplyr::select( -Chromosome,-Coordinate,-REF,-ALT )
  
  if ( uval ) x <- x %>% dplyr::distinct( Target_ID, .keep_all = TRUE )
  
  x <- x %>%
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
      RGT = dplyr::case_when(
        GTS <= gts_min ~ NA_real_,
        is.na(GTS) ~ NA_real_,
        TRUE ~ RGT
      ) %>% as.integer(),
      AGT = dplyr::case_when(
        GTS <= gts_min ~ NA_real_,
        is.na(GTS) ~ NA_real_,
        TRUE ~ AGT
      ) %>% as.integer(),
      GTC = dplyr::case_when(
        # GTS <= gts_min ~ 5.0,
        # is.na(RGT) ~ 5.0,
        # is.na(AGT) ~ 5.0,
        GTS <= gts_min ~ fval,
        is.na(RGT) ~ fval,
        is.na(AGT) ~ fval,
        
        RGT==0 & AGT==0 ~ 0.0,
        RGT==0 & AGT==1 ~ 1.0,
        RGT==1 & AGT==0 ~ 2.0,
        RGT==1 & AGT==1 ~ 3.0,
        TRUE ~ NA_real_ ) %>% as.integer()
    ) %>% dplyr::select( -GT, -GS ) %>% clean_tib()
}

stack_snv_list = function( x, ssh = NULL, sub = NULL )
{
  x <- x %>% dplyr::bind_rows( .id = "Sentrix_Name" )
  if ( !is.null(ssh) ) x <- ssh %>% 
      dplyr::select( Sentrix_Name, Plot_ID ) %>% 
      dplyr::right_join( x, by=c("Sentrix_Name"), multiple = "all" ) %>%
      dplyr::select( -Sentrix_Name ) %>%
      dplyr::rename( Sample_Name = Plot_ID )
  
  if ( !is.null(sub) ) x <- x %>% 
      dplyr::mutate( Probe_ID = Target_ID %>% stringr::str_remove(":.*$") ) %>% 
      dplyr::filter( Probe_ID %in% sub ) %>% 
      dplyr::select( -Probe_ID )
  
  x
}

split_snv_stack = function( x, nstr="Sentrix_ID", istr="Target_ID" ) {
  x %>%
    tidyr::pivot_longer( cols = c( Qual,Filter,PVF,GTS,RGT,AGT,GTC ), 
                         names_to  = c("Key"), 
                         values_to = c("Val") ) %>% split(.$Key) %>% 
    
    lapply( function(x) { 
      x %>% dplyr::select(-Key) %>% 
        tidyr::pivot_wider( id_cols = dplyr::all_of(istr), 
                            names_from = dplyr::all_of(nstr), 
                            values_from = c(Val) )
    }) %>%
    lapply( function(x) {
      x %>% tibble::column_to_rownames( var = istr ) %>% 
        as.data.frame() %>% as.matrix()
    })
}

gtc_snv_to_performance = function( x ) {
  
  rcnt <- base::nrow( x )
  ccnt <- base::ncol( x )
  
  ret_tib <- NULL
  for ( rr in c(1:rcnt) ) {
    tcnt <- 0
    mcnt <- 0
    
    for ( ii in c(1:ccnt) ) {
      for ( jj in c(1:ccnt) ) {
        if ( ii >= jj ) next
        
        if ( x[ rr,ii] < 6 && x[ rr,jj] < 6 ) tcnt = tcnt + 1
        if ( x[ rr,ii] < 5 && x[ rr,jj ] < 5 &&
             x[ rr,ii ] == x[ rr,jj] ) mcnt <- mcnt + 1
      }
    }
    
    ret_tib <- ret_tib %>% 
      dplyr::bind_rows(
        tibble::tibble(
          Probe_ID = rownames(x)[rr],
          Total_Cnt = tcnt,
          Match_Cnt = mcnt
        )
      )
    
  }
  
  ret_tib
}

snv_stack_to_contingency_tab = function( x, single=FALSE ) {
  
  col_vec <- c(1:base::ncol(x) )
  row_vec <- c(1:base::nrow(x) )
  
  hit_mat <- NULL
  hit_mat <- matrix( data=0, nrow = length(col_vec), ncol=length(col_vec) )
  colnames(hit_mat) <- colnames(x)
  rownames(hit_mat) <- colnames(x)
  
  mis_mat <- NULL
  mis_mat <- matrix( data=0, nrow = length(col_vec), ncol=length(col_vec) )
  colnames(mis_mat) <- colnames(x)
  rownames(mis_mat) <- colnames(x)
  
  nan_mat <- NULL
  nan_mat <- matrix( data=0, nrow = length(col_vec), ncol=length(col_vec) )
  colnames(nan_mat) <- colnames(x)
  rownames(nan_mat) <- colnames(x)
  
  for ( ii in col_vec ) {
    for ( jj in col_vec ) {
      hit_cnt <- which( x[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 ) %>% length()
      mis_cnt <- which( x[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() != 0 ) %>% length()
      nan_cnt <- which( is.na(x[ ,ii ]) | is.na(x[ ,jj ]) ) %>% length()
      
      if ( hit_cnt + mis_cnt + nan_cnt != length(row_vec) ) {
        stop(glue::glue("{pmssg} Failed: ii={ii}, jj={jj}, mat={hit_cnt}, mis={mis_cnt}, nan={nan_cnt}.{RET2}"))
      }
      
      hit_mat[ii,jj] = hit_cnt
      mis_mat[ii,jj] = mis_cnt
      nan_mat[ii,jj] = nan_cnt
      
      if ( single ) break
    }
    if ( single ) break
  }
  
  ret_dat <- NULL
  ret_dat$hit_mat <- hit_mat
  ret_dat$mis_mat <- mis_mat
  ret_dat$nan_mat <- nan_mat
  
  ret_dat
}

# filter_snv_stack( x, dat_key, min_key, min_val = 0 ) {
#   x[ which( x[[min_key]] < min_val ) ] <- NA_real_
# }



analyze_snvs = function( vcfs, run_tag, ssh_tib, 
                         sub_vec = NULL,
                         gts_min = 0, 
                         fval = 5.0, # fval=NA_real_, 
                         jcnt = 0, 
                         uval = TRUE,
                         outDir, write_out=FALSE, plot_heat=FALSE ) 
{
  run_tag <- paste0( run_tag, "_gts",gts_min)
  
  outDir <- safe_mkdir( dir = file.path( outDir, run_tag ) )
  rdsDir <- safe_mkdir( dir = file.path( outDir, "rds") )
  csvDir <- safe_mkdir( dir = file.path( outDir, "csv") )
  pdfDir <- safe_mkdir( dir = file.path( outDir, "pdf") )
  
  snv_rds <- file.path( rdsDir, paste0( run_tag,".rds") )
  hit_csv <- file.path( csvDir, paste0( run_tag,".hit.csv.gz") )
  mis_csv <- file.path( csvDir, paste0( run_tag,".mis.csv.gz") )
  nan_csv <- file.path( csvDir, paste0( run_tag,".nan.csv.gz") )
  hit_pdf <- file.path( pdfDir, paste0( run_tag,".hit.heatmap.pdf") )
  
  snv_dat <- NULL
  snv_dat <- vcfs %>% # head(n=3) %>% 
    lapply( read_snv_vcf ) %>% 
    lapply( parse_snv_vcf, gts_min = gts_min, fval = fval, 
            jcnt = jcnt, uval = uval ) %>%
    stack_snv_list( ssh = ssh_tib, sub = sub_vec ) %>%
    split_snv_stack( nstr = "Sample_Name", istr = "Target_ID" ) 
  
  snv_tab <- NULL
  if ( write_out || plot_heat ) snv_tab <- snv_stack_to_contingency_tab( x = snv_dat$GTC )
  
  if ( write_out ) {
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
  
  if ( plot_heat ) {
    pdf( file = hit_pdf, width = 10, height = 10 )
    heatmap( snv_tab$hit_mat )
    dev.off()
  }
  
  if ( FALSE ) {
    
    prb_tib <- NULL
    prb_tib <- tibble::tibble(
      Exp_Tag   = run_tag,
      Target_ID = rownames(snv_dat$GTC),
      Probe_ID  = Target_ID %>% stringr::str_remove(":.*$"),
      Loci_ID   = Probe_ID %>% stringr::str_remove("_.*$"),
      Avg_Mat   = snv_dat$GTC %>% matrixStats::rowMeans2( na.rm = TRUE ),
      Med_Mat   = snv_dat$GTC %>% matrixStats::rowMedians( na.rm = TRUE ),
      Sds_Mat   = snv_dat$GTC %>% matrixStats::rowSds( na.rm = TRUE ),
      Mad_Mat   = snv_dat$GTC %>% matrixStats::rowMads( na.rm = TRUE ) )
    
    ret_dat <- NULL
    ret_dat$snv_dat = snv_dat
    if ( !is.null(snv_tab) ) ret_dat$snv_tab = snv_tab
    ret_dat$prb_tib = prb_tib
    
    return( ret_dat )
  }
 
  snv_dat$GTC[ which( is.na(snv_dat$GTC) ) ] <- 6
  
  ret_dat <- NULL
  ret_dat <- snv_dat$GTC
  
  ret_dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Plot SNV Heatmap Fuctions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_to_heatmap = function( gtc_list,
                           sub_vec,
                           
                           pre_key = "Exp_Name",
                           ids_key = "Sample_Name",
                           suf_key = "Call_Rate",
                           
                           out_dir,
                           run_tag,
                           
                           reload     = 0,
                           reload_min = 2,
                           reload_pre = NULL,
                           
                           ret_data   = FALSE,
                           parallel   = FALSE,
                           write_out  = FALSE,
                           
                           vb=0, vt=3, tc=1, tt=NULL,
                           fun_tag='snp_to_heatmap')
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
  
  ids_sym <- rlang::sym( ids_key )
  
  sep_vec <- NULL
  sep_vec <- c( 1:(length(sub_vec)-1) )
  snp_vec <- NULL
  snp_vec <- paste0( "v",c( 1:length(sub_vec) ) )
  # snp_vec <- paste0( "v",sep_vec )
  
  gtc_tib <- NULL
  for ( gtc_key in names(snp_list) ) {
    gtc_tib <- gtc_tib %>%
      dplyr::bind_rows( gtc_list[[gtc_key]] )
  }
  
  gtc_tib %>% 
    dplyr::group_by( ids_sym ) %>%
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
  
  
  # Set Name Vectors
  
  # Convert to Matrix
  
  # Build Output Heatmap Matrix
  
  # Write Heatamp PDF, data & mapping files...

  if ( FALSE ) {
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
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Parse dbSNP Tabix Fuctions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

parse_dbSNP_vcf = function( snp_vcf,
                            src_tib = NULL,
                            
                            ref_vcf   = TRUE, 
                            com_len   = 56, 
                            rm_chr0   = TRUE,
                            strip_chr = TRUE, 
                            snps_only = TRUE,
                            match_src = TRUE,
                            parse_caf = TRUE,
                            min_maf   = 0.0,
                            
                            out_dir,
                            run_tag,
                            
                            reload     = 0,
                            reload_min = 2,
                            reload_pre = NULL,
                            
                            ret_data   = FALSE,
                            parallel   = FALSE,
                            write_out  = FALSE,
                            
                            vb=0, vt=3, tc=1, tt=NULL,
                            fun_tag='parse_dbSNP_vcf')
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
    cat(glue::glue("{mssg}      snp_vcf = '{snp_vcf}'.{RET}"))
    cat(glue::glue("{mssg}      ref_vcf = '{ref_vcf}'.{RET}"))
    cat(glue::glue("{mssg}      com_len = '{com_len}'.{RET}"))
    cat(glue::glue("{mssg}      rm_chr0 = '{rm_chr0}'.{RET}"))
    cat(glue::glue("{mssg}    strip_chr = '{strip_chr}'.{RET}"))
    cat(glue::glue("{mssg}    snps_only = '{snps_only}'.{RET}"))
    cat(glue::glue("{mssg}    match_src = '{match_src}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
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
  
  errs_mssg <- glue::glue("File(snp_vcf)='{snp_vcf}' does not exist")
  if ( !file.exists( snp_vcf) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  tab_tib <- NULL
  if ( ref_vcf ) {
    tab_tib <- readr::read_tsv( snp_vcf,
                                skip = com_len, show_col_types = FALSE  ) %>% 
      magrittr::set_names( c("Chromosome","Coordinate_SNP","SNP_ID",
                             "REF_SNP","ALT_STR","QUAL","FILTER","INFO") )
    
    ret_cnt <- print_tib( tab_tib, fun_tag = fun_tag, name = "SNP_VCF.tib", 
                          vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    if ( snps_only ) tab_tib <- tab_tib %>%
      dplyr::filter( stringr::str_length(REF_SNP) == 1 ) %>%
      # This removes most large indels::
      dplyr::filter( ALT_STR %>% stringr::str_length() <= 5 ) %>%
      dplyr::filter( REF_SNP != ALT_STR ) %>%
      dplyr::filter( !SNP_ID %>% stringr::str_detect(",") ) %>% 
      dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_STR, .keep_all = TRUE ) %>%
      dplyr::select( Chromosome, Coordinate_SNP, SNP_ID,REF_SNP,ALT_STR, INFO ) %>%
      clean_tib() %>% 
      dplyr::mutate( ALT_LenT = ALT_STR %>% stringr::str_length(), 
                     ALT_LenC = ALT_STR %>% stringr::str_remove_all(",") %>% stringr::str_length(),
                     ALT_LenD = ALT_LenT - ALT_LenC ) %>% 
      dplyr::filter( !(ALT_LenT > 1 & ALT_LenD == 0) ) %>%
      dplyr::distinct( SNP_ID, .keep_all = TRUE ) %>%
      dplyr::select( Chromosome, Coordinate_SNP, SNP_ID, REF_SNP, ALT_STR, INFO )
    
    if ( parse_caf ) {
      
      caf_tib <- NULL
      caf_tib <- tab_tib %>%
        dplyr::mutate(
          CAF_STR = INFO %>% stringr::str_remove("^.*CAF=") %>% stringr::str_remove(";.*$")
        ) %>%
        dplyr::select( -INFO ) %>%
        tidyr::separate( ALT_STR, into = c("ALT1","ALT2","ALT3"), sep = ",", 
                         convert = FALSE, fill = "right", remove = FALSE, extra = "drop" ) %>%
        tidyr::separate( CAF_STR, into=c("RAF","MAF1","MAF2","MAF3"), sep = ",", 
                         convert = FALSE, fill = "right", remove = FALSE, extra = "drop" ) %>%
        dplyr::mutate(
          RAF = dplyr::case_when( 
            is.na(RAF) ~ "0.0",
            RAF == "." ~ "0.0",
            TRUE ~ RAF ) %>% as.double(),
          MAF1 = dplyr::case_when( 
            is.na(MAF1) ~ "0.0",
            MAF1 == "." ~ "0.0",
            TRUE ~ MAF1 ) %>% as.double(),
          MAF2 = dplyr::case_when( 
            is.na(MAF2) ~ "0.0",
            MAF2 == "." ~ "0.0",
            TRUE ~ MAF2 ) %>% as.double(),
          MAF3 = dplyr::case_when( 
            is.na(MAF3) ~ "0.0",
            MAF3 == "." ~ "0.0",
            TRUE ~ MAF3 ) %>% as.double()
        )
      
      # [TBD]: This should be done with pivots...
      # Split and Rejoin::
      sel1_tib <- NULL
      sel1_tib <- caf_tib %>% 
        dplyr::mutate( MAF = MAF1, ALT_SNP = as.character(ALT1) ) %>% 
        dplyr::filter( MAF > min_maf & !is.na(ALT_SNP) & stringr::str_length(ALT_SNP) == 1 )
      
      sel2_tib <- NULL
      sel2_tib <- caf_tib %>% 
        dplyr::mutate( MAF = MAF2, ALT_SNP = as.character(ALT2) ) %>% 
        dplyr::filter( MAF > min_maf & !is.na(ALT_SNP) & stringr::str_length(ALT_SNP) == 1 )
      
      sel3_tib <- NULL
      sel3_tib <- caf_tib %>% 
        dplyr::mutate( MAF = MAF3, ALT_SNP = as.character(ALT3) ) %>% 
        dplyr::filter( MAF > min_maf & !is.na(ALT_SNP) & stringr::str_length(ALT_SNP) == 1 )
      
      # [TBD]: Report sel1,sel2_sel3 lengths after filter...
      
      tab_tib <- dplyr::bind_rows( sel1_tib,sel2_tib,sel3_tib ) %>%
        dplyr::arrange( Chromosome,Coordinate_SNP,SNP_ID, -MAF ) %>%
        dplyr::select( Chromosome,Coordinate_SNP,SNP_ID,
                       REF_SNP,ALT_SNP,RAF,MAF, dplyr::everything() )
    }
    
  } else {
    tab_tib <- readr::read_tsv( snp_vcf )
  }
  
  # Default Return Type: Tabix Results
  ret_tib <- tab_tib
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Analyze Match Types::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( match_src && !is.null(src_tib) ) {
    ret_tib <- dplyr::bind_rows(
      dplyr::inner_join( 
        tab_tib, src_tib %>% dplyr::mutate( Coordinate_SNP=Nxb_Pos_Up ),
        by=c("Chromosome","Coordinate_SNP"), suffix=c("_SNP","_CPG"),
        multiple = "all" ) %>% 
        dplyr::mutate( SNP_Match_Type="Nxb_Up"),
      
      dplyr::inner_join( 
        tab_tib, src_tib %>% dplyr::mutate( Coordinate_SNP=Cpg_Pos_Up ),
        by=c("Chromosome","Coordinate_SNP"), suffix=c("_SNP","_CPG"),
        multiple = "all"  ) %>% 
        dplyr::mutate( SNP_Match_Type="Cpg_Up"),
      
      dplyr::inner_join( 
        tab_tib, src_tib %>% dplyr::mutate( Coordinate_SNP=Cpg_Pos_Dn ),
        by=c("Chromosome","Coordinate_SNP"), suffix=c("_SNP","_CPG"),
        multiple = "all"  ) %>% 
        dplyr::mutate( SNP_Match_Type="Cpg_Dn"),
      
      dplyr::inner_join( 
        tab_tib, src_tib %>% dplyr::mutate( Coordinate_SNP=Nxb_Pos_Dn ),
        by=c("Chromosome","Coordinate_SNP"), suffix=c("_SNP","_CPG"),
        multiple = "all"  ) %>% 
        dplyr::mutate( SNP_Match_Type="Nxb_Dn")
    )
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
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Fingerprinting (SNV) Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fingerprint_vcfs = function( vcfs,
                             sam_tib,
                             sam_vec,
                             sub_vec = NULL,
                             run_sig = FALSE,
                             gts_min = 0,
                             
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
  p4  <- vb > vt + 4
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  run_tag <- paste0( run_tag,".gts_min-",gts_min )
  
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
    cat(glue::glue("{mssg}      gts_min = '{gts_min}'.{RET}"))
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
      RGT = dplyr::case_when(
        GTS <= gts_min ~ NA_real_,
        is.na(GTS) ~ NA_real_,
        TRUE ~ RGT
      ) %>% as.integer(),
      AGT = dplyr::case_when(
        GTS <= gts_min ~ NA_real_,
        is.na(GTS) ~ NA_real_,
        TRUE ~ AGT
      ) %>% as.integer(),
      GTC = dplyr::case_when(
        GTS <= gts_min ~ 5.0,
        is.na(RGT) ~ 5.0,
        is.na(AGT) ~ 5.0,
        
        RGT==0 & AGT==0 ~ 0.0,
        RGT==0 & AGT==1 ~ 1.0,
        RGT==1 & AGT==0 ~ 2.0,
        RGT==1 & AGT==1 ~ 3.0,
        TRUE ~ NA_real_ ) %>% as.integer()
    ) %>% dplyr::select( -GT ) %>% clean_tib()
  
  # Subset to a set targets
  if ( !is.null(sub_vec) && length(sub_vec) != 0 ) {
    snp_tab <- snp_tab %>% 
      dplyr::mutate( Probe_ID = Target_ID %>% stringr::str_remove(":.*$") ) %>% 
      dplyr::filter( Probe_ID %in% sub_vec ) %>%
      dplyr::select( -Probe_ID )
  }
  
  snp_sum <- NULL
  snp_sum <- snp_tab %>% 
    dplyr::group_by( RGT,AGT ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  print_tib( snp_tab, fun_tag = fun_tag, name = "SNP-Summary", 
             vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  ret_tib <- snp_tab
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Format Table & Matrix Lists::
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
  vet_rgt_mat <- NULL
  vet_rgt_mat <- snp_mat_lst$RGT
  
  print_tib( vet_gtc_mat, fun_tag = fun_tag, name = "GTC-Matrix", 
             vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Format SNP Martix Lists to Matrix::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  col_vec <- c(1:base::ncol(vet_gtc_mat) )
  row_vec <- c(1:base::nrow(vet_gtc_mat) )
  col_len <- col_vec %>% length()
  row_len <- row_vec %>% length()
  if ( p0 ) cat(glue::glue("{mssg} rows={row_len}, cols={col_len}.{RET}"))
  
  print_tib( sam_tib, fun_tag = fun_tag, name = "SAM_TIB", 
             vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  sam_mat <- NULL
  sam_mat <- sam_tib %>% 
    dplyr::select( Sentrix_Name,Sample_Base,cg_calls_pass_perc_1 ) %>% 
    tibble::column_to_rownames( var = "Sentrix_Name" ) %>% 
    # as.data.frame() %>% 
    as.matrix()
  
  cat(glue::glue("{mssg} col_vec:{RET}"))
  print(col_vec)
  cat(glue::glue("{mssg} col_vec.{RET2}{RET2}{RET2}{tabs}{BRK}{RET2}"))
  
  if ( ret_data ) {
    
    ret_dat$snp_tab <- snp_tab
    ret_dat$sam_mat <- sam_mat
    
    ret_dat$snp_tab_lst <- snp_tab_lst
    ret_dat$snp_mat_lst <- snp_mat_lst
    
    ret_dat$row_vec <- row_vec
    ret_dat$col_vec <- col_vec
    
    # ret_dat$stats_tab <- stats_tab
    # ret_dat$stats_sum <- stats_sum
    
    # ret_dat$fig_gtc_tib <- fig_gtc_tib
    # ret_dat$fig_gtc_lst <- fig_gtc_lst
    
    # return( ret_dat )
  }
  
  # print_tib( sam_mat, fun_tag = fun_tag, name = "SAM_MAT", 
  #            vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  
  # Expectation = [ A tibble: 1,600 × 11 ] - All Samples
  # Expectation = [ A tibble:   576 × 11 ] - Non-sam_vec Samples
  stats_tab <- NULL
  for ( ii in col_vec ) {
    ii_set_key <- colnames(vet_gtc_mat)[ii]
    ii_sam_key <- sam_mat[ ii_set_key, 1 ]
    
    # if ( p0 ) cat(glue::glue("{mssg} colnames(vet_gtc_mat)[ii={ii}/{col_len}] = '{ii_set_key}, {ii_sam_key}'.{RET}"))
    if ( !(ii_sam_key %in% sam_vec) ) next

    for ( jj in col_vec ) {
      jj_set_key <- colnames(vet_gtc_mat)[jj]
      jj_sam_key <- sam_mat[ ii_set_key, 1 ]
      
      # if ( p0 ) cat(glue::glue("{mssg}{TAB} colnames(vet_gtc_mat)[ii={ii},{jj}/{col_len}] = '{jj_set_key}, {jj_sam_key}'.{RET}"))
      if ( !(jj_sam_key %in% sam_vec) ) next
      
      tot_cnt <- vet_gtc_mat[ ,c(ii,jj) ] %>% base::nrow()
      rgc_cnt <- which( !vet_rgt_mat[ ,c(ii,jj) ] %>% matrixStats::rowAnyNAs() ) %>% length()
      agc_cnt <- which( !vet_agt_mat[ ,c(ii,jj) ] %>% matrixStats::rowAnyNAs() ) %>% length()
      
      if ( rgc_cnt != agc_cnt ) {
        errs_mssg <- glue::glue("RGT != AGT Counts ({rgc_cnt} != {agc_cnt})")
        if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
        if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
        if ( eflag ) return(NULL)
      }
      pas_cnt <- rgc_cnt
      
      nan_cnt <- which( vet_gtc_mat[ ,ii ] == 5 | vet_gtc_mat[ ,jj ] == 5 ) %>% length()
      mat_cnt <- which( vet_gtc_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 &
                          vet_gtc_mat[ ,ii ] != 5 &
                          vet_gtc_mat[ ,jj ] != 5 ) %>% length()
      agt_cnt <- which( vet_agt_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 ) %>% length()
      rgt_cnt <- which( vet_rgt_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() == 0 ) %>% length()
      
      stats_tab <- stats_tab %>%
        dplyr::bind_rows(
          tibble::tibble(
            Sentrix_NameA = ii_set_key,
            Sentrix_NameB = jj_set_key,
            Sample_BaseA  = ii_sam_key,
            Sample_BaseB  = jj_sam_key,
            # Sample_BaseA = sam_mat[ Sentrix_NameA, ][1],
            # Sample_BaseB = sam_mat[ Sentrix_NameB, ][1],
            Idx = ii,
            Jdx = jj,
            
            Tot_Cnt = tot_cnt,
            Pas_Cnt = pas_cnt,
            
            Nan_Cnt = nan_cnt,
            Mat_Cnt = mat_cnt,
            RGT_Cnt = rgt_cnt,
            AGT_Cnt = agt_cnt,
            
            Nan_Per = nan_cnt / tot_cnt,
            Mat_Per = mat_cnt / pas_cnt,
            RGT_Per = rgt_cnt / pas_cnt,
            AGT_Per = agt_cnt / pas_cnt
          )
        )
      
      # print( stats_tab )
      # if ( p2 ) cat(glue::glue("{mssg} ii={ii}, jj={jj} = {tot_cnt}, ",
      #                          "mat={mat_cnt}, nan={nan_cnt}, agt={agt_cnt}, ",
      #                          "rgt={rgt_cnt}.{RET}") )
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
    dplyr::summarise( Mat_Sds = sd( Mat_Per, na.rm = TRUE ),
                      AGT_Sds = sd( AGT_Per, na.rm = TRUE ),
                      RGT_Sds = sd( RGT_Per, na.rm = TRUE ),
                      
                      Mat_Mad = mad( Mat_Per, na.rm = TRUE ),
                      AGT_Mad = mad( AGT_Per, na.rm = TRUE ),
                      RGT_Mad = mad( RGT_Per, na.rm = TRUE ),
                      
                      Nan_Avg = mean( Nan_Per, na.rm=TRUE ),
                      Mat_Avg = mean( Mat_Per, na.rm=TRUE ),
                      AGT_Avg = mean( AGT_Per, na.rm=TRUE ),
                      RGT_Avg = mean( RGT_Per, na.rm=TRUE ),
                      
                      Nan_Med = median( Nan_Per, na.rm=TRUE ),
                      Mat_Med = median( Mat_Per, na.rm=TRUE ),
                      AGT_Med = median( AGT_Per, na.rm=TRUE ),
                      RGT_Med = median( RGT_Per, na.rm=TRUE ),
                      
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
      dplyr::filter( Sample_Base %in% sam_vec )
    
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
    # ret_dat$snp_mat_lst <- snp_mat_lst
    
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

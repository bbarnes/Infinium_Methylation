
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               TripleCrown::
#                         CG# Database Functions::
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
#                        Replicate Screening Functions:
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

consolidate_tibs_to_matrix = function( a, b, keys = NULL,
                                       vb=0, vt=3, tc=1, tt=NULL,
                                       fun_tag='consolidate_tibs_to_matrix')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  eflag <- FALSE
  
  matA <- NULL
  matA <- a %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% 
    as.matrix()
  
  matB <- NULL
  matB <- b %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% 
    as.matrix()
  
  # Set rownames and consolidate matricies
  if ( is.null(keys) ) keys <- base::intersect( rownames(matA), rownames(matB) )
  row_cnts <- length(keys)
  
  errs_mssg <- glue::glue("There is zero overlap between datasets")
  if ( row_cnts == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  matA <- matA[ keys, ]
  matB <- matB[ keys, ]
  
  errs_mssg <- glue::glue("Failed to match rows")
  if ( base::nrow(matA) != base::nrow(matB) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  cbind( matA,matB )
}

calc_dBs = function( beta_tibA,
                     pval_tibA = NULL,
                     name_strA = "UnkA",

                     beta_tibB = NULL,
                     pval_tibB = NULL,
                     name_strB = NULL, # "UnkB",
                     
                     min_pval = 1.0,
                     min_dB   = 0.2,
                     cmp_str  = "lt", # c("lt","gt")
                     is_abs   = TRUE,
                     dBe_max  = 1.0,
                     dPe_max  = 1.0,
                     dBe_fac  = 1.0,
                     dPe_fac  = 1.0,
                     
                     out_dir,
                     run_tag,
                     
                     reload     = 0,
                     reload_min = 2,
                     reload_pre = NULL,
                     
                     ret_data   = FALSE,
                     parallel   = FALSE,
                     
                     vb=0, vt=3, tc=1, tt=NULL,
                     fun_tag='calc_dBs')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt
  p1  <- vb > vt + 1
  p1  <- vb > vt + 2
  
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
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      name_strA = '{name_strA}'.{RET}"))
    cat(glue::glue("{mssg}      name_strAB = '{name_strAB}'.{RET}"))
    cat(glue::glue("{mssg}     min_pval = '{min_pval}'.{RET}"))
    cat(glue::glue("{mssg}       min_dB = '{min_dB}'.{RET}"))
    cat(glue::glue("{mssg}      cmp_str = '{cmp_str}'.{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # errs_mssg <- glue::glue("Outdir out='{out_dir}' does not exist")
  # if ( !dir.exists( out_dir) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Pre-processing:: Matricies
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  is_single = FALSE
  if ( is.null(beta_tibB) ) is_single = TRUE
  if ( is.null(beta_tibB) ) beta_tibB <- beta_tibA
  if ( is.null(pval_tibB) ) pval_tibB <- pval_tibA
  if ( is.null(name_strB) ) name_strB <- name_strA

  ftime <- base::system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Format Matricies:: Beta
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ncolA <- base::ncol( pval_tibA ) - 1
    ncolB <- base::ncol( pval_tibB ) - 1
    idx_vecA <- c( 1:ncolA )
    idx_vecB <- c( 1:ncolA ) + ncolB
    if ( p1 ) cat(glue::glue("{mssg} ncolA={ncolA} x noclB={noclB}.{RET}"))

    beta_mat <- NULL
    beta_mat <- consolidate_tibs_to_matrix( beta_tibA, beta_tibB )
    row_keys <- rownames( beta_mat )
    
    # ret_cnt <- print_tib( beta_mat, fun_tag = fun_tag, name="beta_mat", 
    #                       vb=vb,vt=vt,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Format Matricies:: Pval
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

    pval_mat <- NULL
    if ( !is.null( pval_tibA) ) pval_mat <- 
      consolidate_tibs_to_matrix( pval_tibA, pval_tibB, keys = row_keys )
      
    #
    # [TBD]: This is where we need to split/improve the logic to get both the
    #    weighted dB score and the % passing threshold scores by masking the 
    #    beta_mat with pval_mat
    #  Will require an input flag to decide if we should mask for titration exps
    #
    # e.g.
    # beta_mat[ which( pval_mat > pval_min ) ] <- NA_real_
    #
    beta_matW <- NULL
    beta_matW <- beta_mat
    if ( min_pval < 1.0 ) beta_mat[ which( pval_mat > min_pval ) ] <- NA_real_
    # if ( min_pval < 1.0 ) beta_mat[ which( pval_mat > min_pval ) ] <- NA
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Calculate dBs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    dB_mat <- NULL # Delta Betas
    mP_mat <- NULL # Min Pvals
    wS_mat <- NULL # Weight Score ( Delta Betas * Min_Pvals )
    
    for ( ii in idx_vecA ) {
      for ( jj in idx_vecB ) {
        if ( is_single && ii+ncolA >= jj ) next
        
        if ( p1 ) cat(glue::glue("{mssg} {name_strA} x {name_strB}: {ii} x {jj}.{RET}"))
        
        dP_vec <- NULL
        dB_vec <- NULL
        dB_vec <- cbind( beta_mat[ ,ii ], beta_mat[ ,jj ] ) %>% matrixStats::rowDiffs()
        # dB_vec[ which(!is.finite(dB_vec)) ] <- NA_real_
        if ( is_abs ) dB_vec <- base::abs( dB_vec )
        # dB_vec[ which(!is.finite(dB_vec)) ] <- NA_real_
        
        if ( is.null(dB_mat) ) {
          dB_mat <- base::matrix( data = dB_vec )
          rownames( dB_mat ) <- rownames( beta_mat )
        } else {
          dB_mat <- cbind( dB_mat, dB_vec )
        }

        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Calculate Weighted dB Score::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if ( !is.null(pval_mat) ) {
          
          dB_vecW <- NULL
          dB_vecW <- dB_vec
          if ( min_pval < 1.0 ) {
            dB_vecW <- cbind( beta_matW[ ,ii ], beta_matW[ ,jj ] ) %>% matrixStats::rowDiffs()
            if ( is_abs ) dB_vecW <- base::abs( dB_vecW )
          }
          
          # dB_vecW <- 1 - base::abs( dB_vecW )
          dB_vecW <- ( dBe_max - dB_vecW ) * dBe_fac
          dP_vec <- ( dPe_max - matrixStats::rowMaxs( cbind( pval_mat[ ,ii ], pval_mat[ ,jj ] ) ) ) * dPe_fac
          
          wS_mat <- cbind( wS_mat, cbind( dB_vecW, dP_vec ) %>% matrixStats::rowProds() )
          mP_mat <- cbind( mP_mat, dP_vec )
        }
      }
    }
    # dB_vec
    
    ret_tib <- tibble::tibble( Probe_ID = row_keys )
    
    ret_tib <- ret_tib %>% 
      dplyr::mutate(
        dB_min = matrixStats::rowMins( dB_mat, na.rm = TRUE ),
        dB_avg = matrixStats::rowMeans2( dB_mat, na.rm = TRUE ),
        dB_med = matrixStats::rowMedians( dB_mat, na.rm = TRUE ),
        dB_sds = matrixStats::rowSds( dB_mat, na.rm = TRUE ),
        dB_mad = matrixStats::rowMads( dB_mat, na.rm = TRUE ),
        dB_max = matrixStats::rowMaxs( dB_mat, na.rm = TRUE )
      )
    
    ret_tib <- ret_tib %>% 
      dplyr::mutate(
        dB_tot_cnt = base::ncol( dB_mat ),
        dB_nan_cnt = matrixStats::rowCounts( is.na(dB_mat) ) )
    
    if ( cmp_str == "lte" ) {
      ret_tib <- ret_tib %>% 
        dplyr::mutate(
          dB_pas_cnt = matrixStats::rowCounts( dB_mat <= opt$min_dB, na.rm = TRUE ),
          db_pas_per = base::round( 100 * dB_pas_cnt / dB_tot_cnt, 3 ),
          # db_pas_per = base::round( 100 * dB_pas_cnt / (dB_tot_cnt - dB_nan_cnt), 3 ),
          dB_pas_call = db_pas_per == 100
        )
    } else if ( cmp_str == "gte" ) {
      ret_tib <- ret_tib %>% 
        dplyr::mutate(
          dB_pas_cnt = matrixStats::rowCounts( dB_mat >= opt$min_dB, na.rm = TRUE ),
          db_pas_per = base::round( 100 * dB_pas_cnt / dB_tot_cnt, 3 ),
          dB_pas_call = db_pas_per == 100
        )
    } else if ( cmp_str == "lt" ) {
      ret_tib <- ret_tib %>% 
        dplyr::mutate(
          dB_pas_cnt = matrixStats::rowCounts( dB_mat < opt$min_dB, na.rm = TRUE ),
          db_pas_per = base::round( 100 * dB_pas_cnt / dB_tot_cnt, 3 ),
          dB_pas_call = db_pas_per == 100
        )
    } else if ( cmp_str == "gt" ) {
      ret_tib <- ret_tib %>% 
        dplyr::mutate(
          dB_pas_cnt = matrixStats::rowCounts( dB_mat > opt$min_dB, na.rm = TRUE ),
          db_pas_per = base::round( 100 * dB_pas_cnt / dB_tot_cnt, 3 ),
          dB_pas_call = db_pas_per == 100
        )
    } else {
      errs_mssg <- glue::glue("Compare String='{cmp_str} not supported'")
      eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    
    # dB_key <- paste( "dB_Pass",name_strA,name_strAB, sep="_" )
    # dB_sym <- rlang::sym( dB_key )
    # 
    # ret_tib <- dB_pass_mat %>% as.data.frame() %>% 
    #   tibble::as_tibble() %>% dplyr::rename( !!dB_sym := V1 )
    
    if ( !is.null(pval_mat) ) {
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          mP_min = matrixStats::rowMins( mP_mat, na.rm = TRUE ),
          mP_avg = matrixStats::rowMeans2( mP_mat, na.rm = TRUE ),
          mP_med = matrixStats::rowMedians( mP_mat, na.rm = TRUE ),
          mP_sds = matrixStats::rowSds( mP_mat, na.rm = TRUE ),
          mP_mad = matrixStats::rowMads( mP_mat, na.rm = TRUE ),
          mP_max = matrixStats::rowMaxs( mP_mat, na.rm = TRUE )
        )
      
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          wS_min = matrixStats::rowMins( wS_mat, na.rm = TRUE ),
          wS_avg = matrixStats::rowMeans2( wS_mat, na.rm = TRUE ),
          wS_med = matrixStats::rowMedians( wS_mat, na.rm = TRUE ),
          wS_sds = matrixStats::rowSds( wS_mat, na.rm = TRUE ),
          wS_mad = matrixStats::rowMads( wS_mat, na.rm = TRUE ),
          wS_max = matrixStats::rowMaxs( wS_mat, na.rm = TRUE )
        )
      
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          mP_sum = matrixStats::rowSums2( mP_mat, na.rm = TRUE ),
          wS_sum = matrixStats::rowSums2( wS_mat, na.rm = TRUE ),
          wS_per = base::round( 100 * wS_sum / mP_sum, 3 )
        )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
    #                        done = TRUE, write_spec = TRUE, append = FALSE, 
    #                        fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if ( ret_data ) {
    ret_dat$ret_tib <- ret_tib
    ret_dat$dB_mat  <- dB_mat

    return( ret_dat )
  }
  
  ret_tib
}

probe_replicate_screen = function( tib,
                                   pval_tib = NULL,
                                   pval_min = 0.05,
                                   
                                   sample = "Unk",
                                   min_dB = 0.2,
                                   
                                   out_dir,
                                   run_tag,
                                   
                                   reload     = 0,
                                   reload_min = 2,
                                   reload_pre = NULL,
                                   
                                   ret_data   = FALSE,
                                   parallel   = FALSE,
                                   
                                   vb=0, vt=3, tc=1, tt=NULL,
                                   fun_tag='probe_replicate_screen')
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
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    # cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       sample = '{sample}'.{RET}"))
    cat(glue::glue("{mssg}       min_dB = '{min_dB}'.{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # errs_mssg <- glue::glue("File file='{file}' does not exist")
  # if ( !file.exists( file) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    beta_mat <- NULL
    beta_mat <- tib %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% 
      as.matrix()
    # if ( p4 ) beta_mat %>% head() %>% print()
    ret_cnt <- print_tib( beta_mat, fun_tag = fun_tag, name = "Pre-Beta", 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    if ( !is.null(pval_tib) && pval_min < 1.0 ) {
      
      pval_mat <- NULL
      pval_mat <- tib %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% 
        as.matrix()
      # if ( p4 ) pval_mat %>% head() %>% print()
      ret_cnt <- print_tib( pval_mat, fun_tag = fun_tag, name = "Pval", 
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )
      
      beta_mat[ which( pval_mat > pval_min ) ] <- NA_real_
      # if ( p4 ) beta_mat %>% head() %>% print()
      ret_cnt <- print_tib( beta_mat, fun_tag = fun_tag, name = "Post-Beta", 
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )
      
    }
    
    #
    # Calculate dB Failures::
    #
    dB_pass_mat <- NULL
    dB_pass_mat <- beta_mat[ ,c(1,2) ] %>% matrixStats::rowDiffs() < min_dB
    # if ( p4 ) dB_pass_mat %>% head() %>% print()
    ret_cnt <- print_tib( beta_mat, fun_tag = fun_tag, name = "dB_pass_mat", 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    ncol_cnt <- ncol( beta_mat )
    for ( ii in c(1:ncol_cnt) ) {
      for ( jj in c(1:ncol_cnt) ) {
        if ( ii >= jj ) next
        
        if ( p1 ) cat(glue::glue("{mssg} {sample}: {ii} x {jj}.{RET}"))
        
        dB_pass_mat <- dB_pass_mat &
          beta_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() < min_dB
      }
    }
    
    dB_key <- paste( "dB_Pass",sample, sep="_" )
    dB_sym <- rlang::sym( dB_key )
    
    ret_tib <- dB_pass_mat %>% as.data.frame() %>% 
      tibble::as_tibble() %>% dplyr::rename( !!dB_sym := V1 )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
    #                        done = TRUE, write_spec = TRUE, append = FALSE, 
    #                        fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

probe_replicate_screen2 = function( tibA,
                                    tibB,
                                    sampleA = "UnkA",
                                    sampleB = "UnkB",
                                    min_dB = 0.2,
                                    less_than = TRUE,
                                    
                                    out_dir,
                                    run_tag,
                                    
                                    reload     = 0,
                                    reload_min = 2,
                                    reload_pre = NULL,
                                    
                                    ret_data   = FALSE,
                                    parallel   = FALSE,
                                    
                                    vb=0, vt=3, tc=1, tt=NULL,
                                    fun_tag='probe_replicate_screen2')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt
  p1  <- vb > vt + 1
  p1  <- vb > vt + 2
  
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
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    # cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      sampleA = '{sampleA}'.{RET}"))
    cat(glue::glue("{mssg}      sampleB = '{sampleB}'.{RET}"))
    cat(glue::glue("{mssg}       min_dB = '{min_dB}'.{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # errs_mssg <- glue::glue("File file='{file}' does not exist")
  # if ( !file.exists( file) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    beta_matA <- NULL
    beta_matA <- tibA %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% 
      as.matrix()
    
    beta_matB <- NULL
    beta_matB <- tibB %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% 
      as.matrix()
    
    #
    # Calculate dB Failures::
    #
    dB_pass_mat <- NULL
    dB_pass_mat <- beta_matA[ ,c(1,2) ] %>% matrixStats::rowDiffs() < min_dB
    
    dB_pass_matA <- NULL
    dB_pass_matA <- beta_matA[ ,c(1,2) ] %>% matrixStats::rowDiffs() < min_dB
    
    dB_pass_matB <- NULL
    dB_pass_matB <- beta_matB[ ,c(1,2) ] %>% matrixStats::rowDiffs() < min_dB
    
    ncol_cntA <- ncol( beta_matA )
    ncol_cntB <- ncol( beta_matB )
    for ( ii in c(1:ncol_cntA) ) {
      for ( jj in c(1:ncol_cntB) ) {
        if ( ii >= jj ) next
        
        if ( p1 ) cat(glue::glue("{mssg} {sampleA} x {sampleB}: {ii} x {jj}.{RET}"))
        
        if ( less_than ) {
          dB_pass_mat <- dB_pass_mat &
            cbind( beta_matA[ ,ii ], beta_matB[ , jj ] )  %>% 
            matrixStats::rowDiffs() < min_dB
        } else {
          dB_pass_mat <- dB_pass_mat &
            cbind( beta_matA[ ,ii ], beta_matB[ , jj ] )  %>% 
            matrixStats::rowDiffs() > min_dB
        }
        
        #  beta_matA[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() < min_dB
      }
    }
    
    dB_key <- paste( "dB_Pass",sampleA,sampleB, sep="_" )
    dB_sym <- rlang::sym( dB_key )
    
    ret_tib <- dB_pass_mat %>% as.data.frame() %>% 
      tibble::as_tibble() %>% dplyr::rename( !!dB_sym := V1 )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
    #                        done = TRUE, write_spec = TRUE, append = FALSE, 
    #                        fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Local Run Time Defaults::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imProbeQC_options = function( pars, args,
                              
                              vb=0, vt=4, tc=1, tt=NULL,
                              fun_tag='imProbeQC_options') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  eflag <- FALSE
  wflag <- FALSE
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb>=vt+2) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   args={args}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Parse Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts <- base::list()
  
  # opts$run_name     <- NULL
  opts$out_path     <- NULL
  
  opts$sample_max   <- 0
  opts$platform     <- "EPIC"
  
  opts$ref_path     <- NULL
  opts$ref_file     <- NULL
  opts$ref_build    <- NULL
  opts$ref_source   <- NULL
  opts$ref_species  <- NULL
  
  opts$canonical_csv <- file.path( pars$aux_path, "canonical_cgn_top_grp.csv.gz")
  opts$temp_off_csv  <- file.path( pars$prgm_aux_path, "template_offset.csv.gz")
  
  opts$Rscript      <- "Rscript"
  opts$single       <- FALSE
  opts$cluster      <- FALSE
  opts$parallel     <- FALSE
  opts$track_time   <- TRUE
  opts$clean        <- FALSE
  opts$reload       <- 0
  opts$verbose      <- 3
  
  if (opts$verbose > 0)
    cat(glue::glue("[{pars$prgm_tag}]: Starting; {pars$prgm_tag}.{RET2}"))
  
  if (pars$run_mode == 'RStudio') {
    
    pars$local_run_type <- NULL
    
    pars$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R") %>%
      stringr::str_remove("/tools/imProbeQC/scripts/R") %>%
      stringr::str_remove("/tools/imSuite/scripts/R")
    
    opts$top_path <- pars$top_path
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Params Defaults::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( is.null(pars$version) ) pars$version <- '1'
    if ( is.null(pars$version_key) ) pars$version_key <- "v"
    opts$version  <- paste0(pars$version_key,pars$version)
    
    opts$out_path  <- file.path( pars$top_path, 'scratch' )
    opts$imp_path  <- file.path( pars$top_path, 'data/improbe' )
    opts$ann_path  <- file.path( pars$top_path, 'data/annotation' )
    opts$man_path  <- file.path( pars$top_path, 'data/manifests' )
    opts$idat_path <- file.path( pars$top_path, 'data/idats' )
    
    if (  is.null(pars$run_name) ) opts$run_name = "SNP"
    if ( !is.null(pars$run_name) ) opts$run_name = pars$run_name
    
    opts$ref_source <- paste( "NCBI", sep = ',' )
    # opts$ref_source <- paste( "UCSC", sep = ',' )
    
    opts$run_name <- paste(opts$run_name,opts$ref_source,opts$version, sep='-')
    if ( opts$ref_source == "UCSC" ) {
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$chr_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "hg38.fa.gz",
                              "hg19.fa.gz",
                              "hg18.fa.gz",
                              "mm10.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 sep = ',' )
      
      opts$ref_build <- paste( "hg38",
                               "hg19",
                               "hg18",
                               "mm10",
                               sep = ',' )
      
      
    } else if ( opts$ref_source == "NCBI" ) {
      
      # opts$ref_path <- paste(
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta' ),
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta' ),
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh36/Sequence/WholeGenomeFasta' ),
      #   file.path( pars$top_path, 'data/imGenomes/Mus_musculus/NCBI/GRCm10/Sequence/WholeGenomeFasta' ),
      #   sep = ',' )
      # 
      # opts$chr_path <- paste(
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes' ),
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/Chromosomes' ),
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh36/Sequence/Chromosomes' ),
      #   file.path( pars$top_path, 'data/imGenomes/Mus_musculus/NCBI/GRCm10/Sequence/Chromosomes' ),
      #   sep = ',' )
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$chr_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "GRCh38.genome.fa.gz",
                              "GRCh37.genome.fa.gz",
                              "GRCh36.genome.fa.gz",
                              "GRCm38.genome.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 sep = ',' )
      
      opts$ref_build <- paste( "GRCh38",
                               "GRCh37",
                               "GRCh36",
                               "GRCm10",
                               sep = ',' )
      
    } else {
      eflag <- TRUE
      errs_mssg <- glue::glue("Unsupported default ref_source = {ref_source}")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    # Define tru-cgn-map csv::
    opts$cgn_tru_csv <- file.path( 
      pars$top_path, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-top/bins-100/map.csv.gz" )
    
    # opts$single       = TRUE
    opts$single       = FALSE
    
    opts$cluster      = FALSE
    opts$parallel     = TRUE
    opts$track_time   = TRUE
    opts$clean        = FALSE
    opts$reload       = 1
    opts$verbose      = 3
    
    opts$workflow = "parse_genome,parse_chromosome,update_database"
    
    #
    # ewas_loci_variation::
    #
    if ( pars$prgm_tag == "ewas_loci_variation" ||
         pars$prgm_tag == "ewas_swifthoof" ||
         pars$prgm_tag == "ewas_loci_variation_fileBased" ||
         pars$prgm_tag == "stable_loci_variation_anaysis_scratch" ||
         pars$prgm_tag == "stable_registration_score_scratch" ) {
      
      opts$sample_max <- 0
      opts$platform <- "EPIC"
      opts$manifest <- file.path( opts$top_path, "Projects/EWAS/data/manifests/EWAS_PQC122021-NA-NA-GRCh37_sesame.beta.historic-beadPool.csv.gz" )
      opts$sample_sheet <- file.path( opts$top_path, "Projects/EWAS/data/sample_sheets/EWAS-alpha-cross-product-testing.07012022.sampleSheet.csv.gz" )
    }
    if ( pars$prgm_tag == "ewas_loci_variation"  ||
         pars$prgm_tag == "ewas_loci_variation_fileBased" ||
         pars$prgm_tag == "stable_loci_variation_anaysis_scratch" ) {
      opts$sdf_path  <- file.path( opts$top_path, "scratch/stable_ewas_swifthoof_scratch/sesame-UCSC-v3/prefix_to_sdf" )
      # opts$sdf_path  <- file.path( opts$top_path, "scratch/stable_loci_variation_R_scratch/CEPH/delta_beta/CEPH" )
    }
    
  } else if (pars$run_mode == 'Command_Line') {
    
    options_list <- list(
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Run Time Version Options:: 
      #                       Platform, Genome Build, etc
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--run_name"), type="character", default=opts$run_name, 
        help=glue::glue(
          "Build run name.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--workflow"), type="character", default=opts$workflow, 
        help=glue::glue(
          "Workflow(s) to be executed.",
          "{RET}{TAB2} e.g. parse_genome,parse_chromosome,update_database, etc.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--top_path"), type="character", default=opts$top_path,
        help=glue::glue(
          "Top directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--out_path"), type="character", default=opts$out_path,
        help=glue::glue(
          "Build output directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                 Options:: Loci Variation/EWAS-Swifthoof
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--sdf_path"), type="character", default=opts$sdf_path,
        help=glue::glue(
          "Pre-built SDF path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--idat_path"), type="character", default=opts$idat_path,
        help=glue::glue(
          "Idats path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--sample_sheet"), type="character", default=opts$sample_sheet,
        help=glue::glue(
          "Sample Sheet path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--manifest"), type="character", default=opts$manifest,
        help=glue::glue(
          "Manifest path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--platform"), type="character", default=opts$platform,
        help=glue::glue(
          "Infinium Methylation Platform path (i.e. EPIC, HM450, MM285, etc.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--sample_max"), type="integer", default=opts$sample_max, 
        help=glue::glue(
          "Maximum Number of Samples to Process.{RET}",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Program Specific Options:
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--ref_path"), type="character", default=opts$ref_path, 
        help=glue::glue(
          "Reference Genome directory path(s).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--chr_path"), type="character", default=opts$chr_path, 
        help=glue::glue(
          "Reference Genome directory path(s) for individual chromosomes.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_file"), type="character", default=opts$ref_file, 
        help=glue::glue(
          "Reference Genome fasta file name(s).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_build"), type="character", default=opts$ref_build, 
        help=glue::glue(
          "Reference Genome build names(s).",
          "{RET}{TAB2} e.g. GRch38,GRCh37,GRCh36,GRCm37.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_source"), type="character", default=opts$ref_source, 
        help=glue::glue(
          "Reference Source(s).",
          "{RET}{TAB2} e.g. UCSC, NCBI.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_species"), type="character", default=opts$ref_species, 
        help=glue::glue(
          "Reference Specie(s).",
          "{RET}{TAB2} e.g. Homo_sapiens, Mus_musculus, Rattus_norvegicus, ",
          "SARS-CoV-2.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--canonical_csv"), type="character", default=opts$canonical_csv, 
        help=glue::glue(
          "CSV file(s) containg canonical (already defined) CG#/Top Sequence ",
          "Template definitions.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--temp_off_csv"), type="character", default=opts$temp_off_csv, 
        help=glue::glue(
          "CSV file containg Top Sequence Template start and end offsets from ",
          "CG# genomic postion.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Run Time Mode Options::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--Rscript"), type="character", default=opts$Rscript,
        help=glue::glue(
          "Rscript executable path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # Process Parallel/Cluster Parameters::
      optparse::make_option(
        c("--single"), action="store_true", default=opts$single, 
        help=glue::glue(
          "Boolean variable to run a single sample on a single-core.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--parallel"), action="store_true", default=opts$parallel, 
        help=glue::glue(
          "Boolean variable to run parallel on multi-core.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--cluster"), action="store_true", default=opts$cluster,
        help=glue::glue(
          "Boolean variable to run jobs on cluster.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      # Run=time Options::
      optparse::make_option(
        c("--track_time"), action="store_true", default=opts$track_time,
        help=glue::glue(
          "Boolean variable tack run times.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--clean"), action="store_true", default=opts$clean, 
        help=glue::glue(
          "Boolean variable to run a clean build.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--reload"), type="integer", default=opts$reload, 
        help=glue::glue(
          "Integer value to reload intermediate files (for testing).{RET}",
          "{TAB2} Zero indicates no-reloading, higher numbers more reloads.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer"),
      
      # Verbosity level::
      optparse::make_option(
        c("-v", "--verbose"), type="integer", default=opts$verbose, 
        help=glue::glue(
          "Verbosity level: 0-5 (5 is very verbose).",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer")
    )
    
    opt_parser = optparse::OptionParser( option_list = options_list )
    opts = optparse::parse_args( opt_parser )
    
    if ( is.null(opts$top_path) && !is.null(pars$top_path) )
      opts$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R")
    
  } else {
    stop( glue::glue("{RET}[{pars$prgm_tag}]: ERROR: Unrecognized run_mode = ",
                     "'{pars$run_mode}'!{RET2}") )
  }
  
  #
  # Parameter Validation:: ewas_loci_variation
  #
  if ( pars$prgm_tag == "ewas_loci_variation" ||
       pars$prgm_tag == "ewas_swifthoof"  ||
       pars$prgm_tag == "ewas_loci_variation_fileBased" ||
       pars$prgm_tag == "stable_loci_variation_anaysis_scratch" ) {
    
    if ( !file.exists(opts$manifest) ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("Manifest does not exist: '{opts$manifest}'")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    if ( !file.exists(opts$sample_sheet) ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("Sample Sheet does not exist: '{opts$sample_sheet}'")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
  }
  if ( pars$prgm_tag == "ewas_loci_variation" ) {
    if ( !dir.exists(opts$sdf_path) ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("SDF Directory does not exist: '{opts$sdf_path}'")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
  }
  
  if ( opts$clean ) opts$reload <- -1
  
  ret_cnt <- length(opts)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

# End of file

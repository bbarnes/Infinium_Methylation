
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               TripleCrown::
#                         CG# Database Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Command Line Options Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )

# Tidyverse Core Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("stringr", quietly = TRUE) ))
suppressWarnings(suppressPackageStartupMessages( 
  base::require("readr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("glue",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("furrr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("purrr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("magrittr", quietly = TRUE) ) )

# Parallel Processing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

# Plotting Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("ggplot2", quietly = TRUE) ) )

# Matrix Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("matrixStats", quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("scales",      quietly = TRUE) ))

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
#                        Swifthoof IO Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


idat_manifests_overlap = function( x,
                                   mans,
                                   
                                   field = "Quants",
                                   
                                   add_key  = "ordering",
                                   add_keys = c("U","M"),
                                   
                                   ctl_key  = "controls",
                                   ctl_keys = "Address",
                                   
                                   col = "Red",
                                   zip = FALSE,
                                   
                                   vb=0, vt=6, tc=1, tt=NULL,
                                   fun_tag='idat_manifests_overlap')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    # cat(glue::glue("{mssg}    prefix(x) = '{x}'.{RET}"))
    cat(glue::glue("{mssg}          col = '{col}'.{RET}"))
    cat(glue::glue("{mssg}          zip = '{zip}'.{RET}"))
    cat(glue::glue("{mssg}        field = '{field}'.{RET}"))
    cat(glue::glue("{mssg}      add_key = '{add_key}'.{RET}"))
    cat(glue::glue("{mssg}      ctl_key = '{ctl_key}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Validate Inputs:: Idat
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  can_vec  <- NULL
  can_type <- NA_character_
  
  if ( purrr::is_character( x ) && base::length( x ) ==  1 ) {
    idat_path <- NULL
    if ( file.exists( x ) ) {
      can_type <- "Idat_File"
      idat_path <- x
    } else if ( dir.exists( base::dirname( x ) ) ) {
      can_type <- "Idat_Prefix"
      idat_path <- prefix_to_idat_path( x = x, col = col, zip = zip, 
                                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    } else {
      eflag <- TRUE
      errs_mssg <- glue::glue("Input is neither file or prefix='{x}'")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    can_vec <- illuminaio::readIDAT( idat_path, what="IlluminaID" ) %>%
      base::as.integer()
    
  } else if ( purrr::is_vector( x ) && !purrr::is_list( x ) ) {
    can_type <- "Vector"
    can_vec <- x %>% base::as.integer()
    
  } else if ( tibble::is_tibble( x ) || 
              base::is.data.frame( x ) ) {
    
    if ( field %in% names(x) ) {
      can_type <- paste("Idat_Tibble",field, sep="_")
      can_vec <- dplyr::pull(x, field ) %>% purrr::as_vector() %>%
        as.integer()
    } else {
      eflag <- TRUE
      errs_mssg <- glue::glue("Idat Tibble Field={field} does not exist.")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    
  } else if ( purrr::is_list( x ) && base::length(x) > 1 ) {
    if ( field %in% names(x) ) {
      can_type <- paste("Idat_DataList",field, sep="_")
      can_vec <- x[[field]] %>% rownames() %>% as.integer()
    } else {
      eflag <- TRUE
      errs_mssg <- glue::glue("Idat List Field={field} does not exist.")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    
  } else {
    eflag <- TRUE
    errs_mssg <- glue::glue("Unssuported idat type")
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
  }
  can_vec <- unique( can_vec )
  can_len <- base::length( can_vec )
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Input Type = '{can_type}', n = {can_len}.{RET}"))
  
  if ( !purrr::is_vector( can_vec) || 
       !purrr::is_integer( can_vec ) ) {
    eflag <- TRUE
    errs_mssg <- glue::glue("Input data did not convert to integer vector")
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Validate Inputs:: Manifests
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ref_keys <- mans %>% names()
  
  for ( ref_key in ref_keys ) {
    if ( vb >= vt+3 ) cat(glue::glue("{mssg} Manifest='{ref_key}'{RET}"))
    
    cur_list <- mans[[ref_key]]
    
    add_tib <- cur_list[[add_key]]
    ctl_tib <- cur_list[[ctl_key]]
    
    sum_cnt <- 0
    add_cnt <- 0
    ctl_cnt <- 0
    
    if ( !is.null(add_tib) ) add_cnt <- add_tib %>% base::nrow()
    if ( !is.null(ctl_tib) ) ctl_cnt <- ctl_tib %>% base::nrow()
    
    if ( add_cnt != 0 ) {
      for ( cur_key in add_keys ) {
        cur_sym <- rlang::sym( cur_key )
        cur_vec <- add_tib %>% 
          dplyr::filter( !is.na( !!cur_sym ) ) %>%
          dplyr::distinct( !!cur_sym ) %>%
          dplyr::pull( !!cur_sym ) %>% 
          purrr::as_vector() %>%
          as.integer() %>%
          unique()
        cur_len <- cur_vec %>% base::length()
        
        int_len <- intersect( can_vec, cur_vec ) %>% unique() %>% base::length()
        sum_cnt <- sum_cnt + int_len
        int_per <- round( 100 * int_len / can_len, 2 )
        
        if ( vb >= vt+4 ) 
          cat(glue::glue("{mssg}{TAB} cur_key='{cur_key}'{RET}",
                         "{mssg}{TAB2} cur_cnt={cur_len}.{RET}",
                         "{mssg}{TAB2} int_len={int_len}.{RET}",
                         "{mssg}{TAB2} can_len={can_len}.{RET}",
                         "{mssg}{TAB2} int_per={int_per}.{RET}" ) )
        
        ret_dat <- ret_dat %>% dplyr::bind_rows(
          tibble::tibble( Manifest = ref_key,
                          Group_Name = cur_key,
                          Group_Match = int_len,
                          Group_Total = cur_len ) )
        
      }
    }
    
    if ( ctl_cnt != 0 ) {
      for ( cur_key in ctl_keys ) {
        cur_sym <- rlang::sym( cur_key )
        cur_vec <- ctl_tib %>% 
          dplyr::filter( !is.na( !!cur_sym ) ) %>%
          dplyr::distinct( !!cur_sym ) %>%
          dplyr::pull( !!cur_sym ) %>% 
          purrr::as_vector() %>%
          as.integer() %>%
          unique()
        cur_len <- cur_vec %>% base::length()
        
        int_len <- intersect( can_vec, cur_vec ) %>% unique() %>% base::length()
        sum_cnt <- sum_cnt + int_len
        int_per <- round( 100 * int_len / can_len, 2 )
        
        if ( vb >= vt+4 ) 
          cat(glue::glue("{mssg}{TAB} cur_key='{cur_key}'{RET}",
                         "{mssg}{TAB2} cur_cnt={cur_len}.{RET}",
                         "{mssg}{TAB2} int_len={int_len}.{RET}",
                         "{mssg}{TAB2} can_len={can_len}.{RET}",
                         "{mssg}{TAB2} int_per={int_per}.{RET}" ) )
        
        #
        # NOTE:: Massive Hack For Control Sub Interesection Nameing::
        #
        if ( cur_key == "Address" ) cur_key <- "C"
        
        ret_dat <- ret_dat %>% dplyr::bind_rows(
          tibble::tibble( Manifest = ref_key,
                          Group_Name = cur_key,
                          Group_Match = int_len,
                          Group_Total = cur_len ) )
        
      }
    }
    
    ret_dat <- ret_dat %>% dplyr::bind_rows(
      tibble::tibble( Manifest = ref_key,
                      Group_Name = "Sum",
                      Group_Match = sum_cnt,
                      Group_Total = can_len ) )
    
  }
  
  ret_tib <- ret_dat %>% 
    tidyr::pivot_wider( names_from  = c( Group_Name ),
                        values_from = c( Group_Match, Group_Total ), values_fill = 0 ) %>% 
    dplyr::mutate( Total_Percent = base::round( 100 * Group_Match_Sum / Group_Total_Sum, 2 ) ) %>% 
    dplyr::arrange( -Group_Match_Sum )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

prefix_to_idat_path = function( x,
                                col,
                                zip = FALSE,
                                
                                vb=0, vt=6, tc=1, tt=NULL,
                                fun_tag='prefix_to_idat_path')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}    prefix(x) = '{x}'.{RET}"))
    cat(glue::glue("{mssg}          col = '{col}'.{RET}"))
    cat(glue::glue("{mssg}          zip = '{zip}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Invalid color='{col}'")
  if ( col != "Grn" && col != "Red" ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  idat_path = paste0( x,"_",col,".idat" )
  errs_mssg <- glue::glue("Idat path is null")
  if ( is.null(idat_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  prefix_dir <- base::dirname( idat_path )
  errs_mssg <- glue::glue("Failed to parse prefix directory")
  if ( is.null(prefix_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("Prefix directory='{prefix_dir}' does not exist")
  if ( !dir.exists( prefix_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( file.exists(idat_path) && zip ) {
    cmd <- paste("gzip",idat_path, sep=' ')
    if ( vb >= vt+1 ) cat(glue::glue("{mssg} Gzipping: cmd='{cmd}'...{RET2}"))
    cmd_ret <- 0
    cmd_ret <- base::system( command = cmd, 
                             ignore.stdout = TRUE, ignore.stderr = TRUE, 
                             wait = TRUE, intern = FALSE )
    
    errs_mssg <- glue::glue("Command return='{cmd_ret}' failed on '{cmd}'")
    if ( is.null(cmd_ret) || cmd_ret != 0 ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    idat_path <- paste( idat_path,"gz", sep='.' )
  } else if ( !file.exists(idat_path) ) {
    idat_path <- paste( idat_path,"gz", sep='.' )
  }
  
  errs_mssg <- glue::glue("Idat path='{idat_path}' does not exist")
  if ( !file.exists(idat_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  # ret_key <- glue::glue("final-ret-tib")
  # ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
  #                       vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  idat_path
}

prefix_to_SigDF = function( x,
                            
                            zip = FALSE,
                            
                            out_dir,
                            run_tag,
                            
                            reload     = 0,
                            reload_min = 2,
                            reload_pre = NULL,
                            
                            ret_data   = FALSE,
                            parallel   = FALSE,
                            
                            vb=0, vt=3, tc=1, tt=NULL,
                            fun_tag='prefix_to_SigDF') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
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
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}    prefix(x) = '{x}'.{RET}"))
    cat(glue::glue("{mssg}          zip = '{zip}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
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
  
  grn_idat_path <- prefix_to_idat_path( x = x, col = "Grn", zip = zip, 
                                        vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  red_idat_path <- prefix_to_idat_path( x = x, col = "Red", zip = zip, 
                                        vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    
    
    warn_mssg <- glue::glue("WARN_MESSAGE")
    if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    wflag <- FALSE
    
    errs_mssg <- glue::glue("ERROR_MESSAGE")
    if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

load_aux_manifests = function(path,
                              vb=0, vt=3, tc=1, tt=NULL,
                              fun_tag='load_aux_manifests') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}      path = '{path}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    
    
    
    warn_mssg <- glue::glue("WARN_MESSAGE")
    if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    wflag <- FALSE
    
    errs_mssg <- glue::glue("ERROR_MESSAGE")
    if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Local Run Time Defaults::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

swifthoof_local_defaults = function( pars, args,
                                     
                                     vb=0, vt=4, tc=1, tt=NULL,
                                     fun_tag='swifthoof_local_defaults')
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
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}{TAB} Local run_mode='{pars$run_mode}'{RET}"))
    cat(glue::glue("{mssg}{TAB} Local top_path='{pars$top_path}'{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  if ( p2 ) cat(glue::glue("{mssg}   args={args}.{RET}"))
  if ( p0 ) cat(glue::glue("{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Parse Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts <- base::list()
  
  #
  # OPTSION FROM ORIGIN::
  #
  
  opts$auto_detect     <- FALSE
  opts$auto_detect_csv <- file.path( pars$prgm_aux_path, "cell_line_auto_detect.csv.gz" )
  
  opts$workflow   <- NULL
  
  opts$pval     <- "pOOBAH,PnegEcdf"
  opts$pval_min <- "0.1,0.02"
  opts$pval_per <- "90,98"
  opts$delta_beta_min <- 0.2
  
  # Output Options::
  opt$load_idat    <- FALSE
  opt$save_idat    <- TRUE
  
  opt$load_sset   <- FALSE
  opt$save_sset   <- FALSE
  
  #
  # TBD: Add new variable names::
  #
  opt$write_beta  <- FALSE
  opt$write_bsum  <- FALSE
  
  opt$write_pval  <- FALSE
  opt$write_psum  <- FALSE
  
  opt$write_sigs  <- FALSE
  opt$write_ssum  <- FALSE
  
  opt$write_call  <- TRUE
  opt$write_csum  <- FALSE
  
  opt$write_snps  <- TRUE
  opt$write_auto  <- FALSE
  
  opt$mask_general <- FALSE
  
  # Threshold Options::
  
  opt$percision_sigs <- 1
  opt$percision_beta <- 4
  opt$percision_pval <- 6
  
  # Parallel/Cluster Options::
  opt$single   <- FALSE
  opt$parallel <- FALSE
  opt$cluster  <- FALSE
  
  # Plotting Options::
  opt$plotSset  <- FALSE
  opt$plotCalls <- FALSE
  opt$plotAuto  <- FALSE
  
  opt$make_pred <- TRUE
  
  opt$plotFormat <- 'pdf'
  opt$plotFormat <- 'png'
  
  opt$dpi <- 72
  opt$dpi <- 120
  
  opt$plotMax <- 10000
  opt$plotSub <- 5000
  
  opt$time_org_txt <- NULL
  opt$trackTime    <- FALSE
  
  # verbose Options::
  opt$verbose <- 3  
  
  
  #
  # COPIED CODE BELOW
  #
  
  opts$run_name     <- NULL
  opts$out_path     <- NULL
  opt$manifest      <- NULL
  opt$platform      <- NULL
  
  opts$aux_man_path  <- file.path( pars$aux_path, "manifests" )
  opts$canonical_csv <- file.path( pars$aux_path, "canonical_cgn_top_grp.csv.gz")
  opts$temp_off_csv  <- file.path( pars$prgm_aux_path, "template_offset.csv.gz" )
  
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
    
    pars$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R")
    
    opts$out_path  <- file.path( pars$top_path, 'scratch' )
    opts$idat_path <- file.path( pars$top_path, 'data/idats' )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Params Local Defaults::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( is.null(pars$version) ) pars$version <- '1'
    if ( is.null(pars$version_key) ) pars$version_key <- "v"
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Params to Options::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( is.null(pars$verbose) ) opts$verbose <- pars$verbose
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                              Target Project::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    opts$platform  <- 'EPIC'
    opts$version   <- paste0(pars$version_key,pars$version)
    opts$run_name  <- paste( par$local_run_type, opts$version, sep="-" )
    opts$idat_path <- 
      file.path( opts$idat_path, paste("idats",par$local_run_type, sep="_" ) )
    
    if ( is.null(par$local_run_type) ) {
      eflag <- TRUE
      errs_mssg <- 
        glue::glue("Parameter local_run_type = '{pars$local_run_type}' is NULL")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
      
    } else if ( par$local_run_type == "NA12878" ) {
      
    } else if ( par$local_run_type == "EX_CAM1" ) {
      
    } else if ( par$local_run_type == "AKE-MVP-Failed-v1" ) {
      
    } else {
      eflag <- TRUE
      errs_mssg <- 
        glue::glue("Unsupported local_run_type = '{pars$local_run_type}'")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    # Ensure idats directory exists!!!
    errs_mssg <- glue::glue("Idats Directory does NOT exist: {opts$idat_path}")
    if ( !dir.exists( opts$idat_path ) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
    
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
        c("--out_path"), type="character", default=opts$out_path,
        help=glue::glue(
          "Build output directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--manifest"), type="character", default=opts$manifest, 
        help=glue::glue(
          "Directory with multiple manifests to auto-detect or full path(s) ",
          "{RET}{TAB2} to manifes(s).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--platform"), type="character", default=opts$platform, 
        help=glue::glue(
          "Forced Sesame platform to use for calculations.",
          "{RET}{TAB2} e.g. EPIC, HM450, HM27, MM285",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      
      
      
      #
      # LEFT OFF HERE
      #
      
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
    
    # Read Options from pre-defined function::
    # option_list <- program_options(pars = par, opts = opt, vb = opts$verbose)
    
    opt_parser = optparse::OptionParser( option_list = options_list )
    opts = optparse::parse_args( opt_parser )
    
  } else {
    stop( glue::glue("{RET}[{pars$prgm_tag}]: ERROR: Unrecognized run_mode = ",
                     "'{pars$run_mode}'!{RET2}") )
  }
  
  ret_cnt <- length(opts)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

# End of file

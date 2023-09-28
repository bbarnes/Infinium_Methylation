
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Common Functions/Methods Used in Many Programs:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )

# Parallel Processing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

# Command Line Options Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("optparse",   quietly = TRUE) ) )
#

# [TBD]: REMOVE OLD STUFF:
#
# Additional Tidy Practices Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("dplyr",  quietly = TRUE) ))
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("stringr",  quietly = TRUE) ))
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("glue",     quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("furrr",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("magrittr", quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("stringr",  quietly = TRUE) ) )

# Matrix Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("matrixStats", quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("scales",      quietly = TRUE) ))

# Rcpp (c++) Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("Rcpp", quietly = TRUE) ) )

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
S10  <- paste0("  ",
               paste(rep("     ",2),collapse=" "))
S15  <- paste0("  ",
               paste(rep("     ",3),collapse=" "))
S20  <- paste0("  ",
               paste(rep("     ",4),collapse=" "))
S25  <- paste0("  ",
               paste(rep("     ",5),collapse=" "))
S30  <- paste0("  ",
               paste(rep("     ",6),collapse=" "))
S35  <- paste0("  ",
               paste(rep("     ",7),collapse=" "))
S40  <- paste0("  ",
               paste(rep("     ",8),collapse=" "))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Standard Template Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func0 = function( tib,
                           vb=0, vt=6, tc=1, tt=NULL,
                           fun_tag='template_func0')
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
  
  warn_mssg <- glue::glue("WARN_MESSAGE")
  if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
  if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
  wflag <- FALSE
  
  errs_mssg <- glue::glue("ERROR_MESSAGE")
  if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

template_func1 = function( tib,
                           vb=0, vt=3, tc=1, tt=NULL,
                           fun_tag='template_func1')
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
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

template_func2 = function( tib,

                           out_dir,
                           run_tag,
                           
                           reload     = 0,
                           reload_min = 2,
                           reload_pre = NULL,
                           
                           ret_data   = FALSE,
                           parallel   = FALSE,
                           write_out  = FALSE,
                           
                           vb=0, vt=3, tc=1, tt=NULL,
                           fun_tag='template_func2')
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
  
  if ( write_out )
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
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

template_func3 = function( tib,
                           file,
                           
                           out_dir,
                           run_tag,
                           
                           reload     = 0,
                           reload_min = 2,
                           reload_pre = NULL,
                           
                           ret_data   = FALSE,
                           parallel   = FALSE,
                           write_out  = FALSE,

                           vb=0,vt=3,tc=1, tt=NULL,
                           fun_tag='template_func2')
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
    cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("File file='{file}' does not exist")
  if ( !file.exists( file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
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
    
    if ( write_out )
      out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
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
#               Silly lapply/mclapply Initialize Functions::
#
#  validate_lapply()
#  validate_lapply1()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

validate_lapply = function( sdfs )
{
  tmp_beta_mat <- NULL
  tmp_beta_mat <- BiocGenerics::do.call(
    BiocGenerics::cbind,
    lapply( sdfs %>% head(n=1), base::nrow )
  )
  # tmp_beta_mat <- BiocGenerics::do.call(
  #   BiocGenerics::cbind,
  #   mclapply( sdfs %>% head(n=1), base::nrow )
  # )
  # tmp_beta_mat %>% dim() %>% print()
  return( TRUE )
}

validate_lapply1 = function()
{
  BiocGenerics::do.call(
    BiocGenerics::cbind,
    lapply( sesameData::sesameDataGet("EPIC.5.SigDF.normal"), base::nrow ) )
  return( TRUE )
}

validate_mclapply1 = function()
{
  validate_lapply1()
  BiocGenerics::do.call(
    BiocGenerics::cbind,
    mclapply( sesameData::sesameDataGet("EPIC.5.SigDF.normal"), base::nrow ) )
  return( TRUE )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Standard File IO Functions::
#
#  safe_mkdir()
#  guess_file_extension()
#  guess_file_delimiter()
#  safe_write()
#  safe_read()
#  valid_time_stamp()
#  sort_files_by_date()
#  get_file_list()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

number_as_commaK <- function(value, cutoff=10000) {
  if (!is.numeric(value))
    cat(glue::glue("[ERROR]: Not-Numeric={value}.{RET}{RET}"))
  stopifnot(is.numeric(value))
  
  s <- ''
  if (value>=cutoff) {
    s <- 'k'
    value <- as.integer(value/1000)
  }
  paste0(format(value, big.mark = ",", scientific = FALSE, digits = 22),s)
}

struct_to_str = function( x, empty = "none", del = "_" ) {
  
  x <- x %>% unlist() %>% paste( collapse = del )
  if ( length(x) == 0 || x == "" ) x <- empty
  
  x
}

safe_mkdir = function( dir, recursive = TRUE, vb=0, vt=6, tc=1, tt=NULL) {
  if ( file.exists(dir) && !dir.exists(dir) ) dir <- base::dirname(dir)
  if ( !dir.exists(dir) ) dir.create(dir, recursive = recursive)
  
  dir
}

guess_file_extension = function( file,
                                 
                                 vb=0, vt=6, tc=1, tt=NULL,
                                 fun_tag='guess_file_extension') {
  
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
    cat(glue::glue("{mssg}   file = '{file}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  errs_mssg <- glue::glue("File is null or length zero")
  if ( is.null(file) || length(file) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  file <- stringr::str_remove(file, ".gz$")
  if ( stringr::str_ends(file,".tsv") ) {
    type <- "tsv"
  } else if ( stringr::str_ends(file,".txt") ) {
    type <- "tsv"
  } else if ( stringr::str_ends(file,".csv") ) {
    type <- "csv"
  } else if ( stringr::str_ends(file,".rds") ) {
    type <- "rds"
  } else {
    type <- "line"
  }
  
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} type='{type}'{RET}"))
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  type
}

guess_file_delimiter = function( file,
                                 n_max = 100,
                                 
                                 vb=0, vt=6, tc=1, tt=NULL,
                                 fun_tag='guess_file_delimiter') {
  
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
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}    file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}   n_max = '{n_max}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  errs_mssg <- glue::glue("File is null or length zero")
  if ( is.null(file) || length(file) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} Reading lines; file='{file}'{RET}"))
  
  del_tib <- NULL
  del_tib <- readr::read_lines( file = file, n_max = n_max )
  
  if ( vb >= vt+6 ) {
    cat(glue::glue("{mssg} del_tib defined 1...{RET2}"))
    del_tib %>% print()
    
    cat(glue::glue("{RET2}{RET2}{mssg} Starting NO vb/vt...{RET2}{RET2}"))
    cat(glue::glue("{RET2}{RET2}{mssg} Starting vb='{vb}'...{RET2}{RET2}"))
    cat(glue::glue("{RET2}{RET2}{mssg} Starting vt='{vt}'...{RET2}{RET2}"))
  }
  
  del_tib <- del_tib %>% 
    tibble::as_tibble() %>%
    dplyr::mutate( Row_Num = dplyr::row_number(),
                   com_cnt = stringr::str_count(value, COM), 
                   tab_cnt = stringr::str_count(value, TAB),
                   spc_cnt = stringr::str_count(value, " ") )
  
  if ( vb >= vt+6 ) cat(glue::glue("{mssg} del_tib defined 3...{RET2}"))
  
  ret_key <- glue::glue("guess-file-delimter-del")
  ret_cnt <- print_tib( del_tib, fun_tag = fun_tag, name = ret_key,
                        vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  med_tib <- NULL
  med_tib <- del_tib %>%
    dplyr::summarise( com = median(com_cnt,na.rm=TRUE),
                      tab = median(tab_cnt,na.rm=TRUE),
                      spc = median(spc_cnt,na.rm=TRUE) )
  
  ret_key <- glue::glue("guess-file-delimter-med")
  ret_cnt <- print_tib( med_tib, fun_tag = fun_tag, name = ret_key,
                        vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  med_key <- NULL
  med_key <- med_tib %>% which.max() %>% names()
  
  sum_tib <- NULL
  sum_tib <- del_tib %>%
    dplyr::summarise( com = sum(com_cnt,na.rm=TRUE),
                      tab = sum(tab_cnt,na.rm=TRUE),
                      spc = sum(spc_cnt,na.rm=TRUE) )
  
  sum_key <- NULL
  sum_key <- sum_tib %>% which.max() %>% names()
  
  ret_key <- glue::glue("guess-file-delimter-sum")
  ret_cnt <- print_tib( sum_tib, fun_tag = fun_tag, name = ret_key,
                        vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  ret_val <- NULL
  if ( med_key == sum_key ) {
    if ( med_key == "com" ) ret_val=COM
    if ( med_key == "tab" ) ret_val=TAB
    if ( med_key == "spc" ) ret_val=" "
  } else {
    eflag <- TRUE
  }
  errs_mssg <- glue::glue("Failed to match file delimiter! ",
                          "med_key={med_key} != sum_key={sum_key}")
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+1 )
    cat(glue::glue("{mssg} med_key='{med_key}'; ret_val='{ret_val}'{RET}"))
  
  ret_tib <- dplyr::bind_rows(med_tib, sum_tib) %>%
    dplyr::mutate(ret_val = ret_val)
  
  ret_key <- glue::glue("guess-file-delimter-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                        vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_val
}

safe_write = function( x, 
                       file,
                       
                       max  = 0,
                       type = NULL,
                       done = TRUE, 
                       write_spec = TRUE, 
                       append = FALSE, 
                       compress = "gz",
                       permissions = NULL,
                       
                       sub_cols = NULL,  # Used to subset data if tibble
                       out_cols = NULL,  # Used to rename columns if tibble
                       
                       vb=0, vt=3, tc=1, tt=NULL,
                       fun_tag='safe_write') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  spec_rds <- stringr::str_remove(file, '.[a-zA-Z]+(.gz)$') %>%
    paste('spec.rds', sep = '.')
  end_txt <- paste(file, 'done.txt', sep='.')
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}   end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{mssg}  spec_rds = '{spec_rds}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      type = '{type}'.{RET}"))
    cat(glue::glue("{mssg}       max = '{max}'.{RET}"))
    cat(glue::glue("{mssg}      done = '{done}'.{RET}"))
    cat(glue::glue("{mssg}    append = '{append}'.{RET}"))
    cat(glue::glue("{mssg}  compress = '{compress}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    out_dir <- safe_mkdir( base::dirname( file ), vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Guess File Type::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    file_type <- type
    if (is.null(file_type)) 
      file_type <- guess_file_extension( file, vb=vb,vt=vt+4,tc=tc+1 )
    
    if (vb >= vt+1)
      cat(glue::glue("{mssg} Will use file_type = '{file_type}'.{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Get Data Type::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    data_type <- "rds"
    if ( purrr::is_vector(x) ) data_type <- "vec"
    if ( purrr::is_list(x) ) data_type <- "list"
    if ( base::is.data.frame(x) || tibble::is_tibble(x) ) data_type <- "tib"
    
    if (vb >= vt+1)
      cat(glue::glue("{mssg} Will use data_type = '{data_type}'.{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Check Data-Type/File-Type Pairing::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    warn <- FALSE
    if (file_type == "csv" && data_type != "tib") warn <- TRUE
    if (file_type == "tsv" && data_type != "tib") warn <- TRUE
    if (warn && vb > 0)
      cat(glue::glue("{warn} Data and File Type are a poor combination:: ",
                     "file_type={file_type}, data_type={data_type}.{RET}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                Modify Data:: Subset-Columns/Rename-Columns
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (!is.null(sub_cols) && !is.null(out_cols)) {
      sub_len <- length(sub_cols)
      out_len <- length(out_cols)
      dat_len <- names(x) %>% length()
      
      if (sub_len != out_len) {
        fail_mssg <- glue::glue("sub_len={sub_len} != out_len={out_len}")
        stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
        return(NULL)
      }
    }
    if (!is.null(sub_cols) && data_type=="tib") {
      sub_len <- length(sub_cols)
      dat_len <- names(x) %>% length()
      
      if (sub_len > dat_len) {
        fail_mssg <- glue::glue("sub_len={sub_len} > dat_len={dat_len}")
        stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
        return(NULL)
      }
      x <- x %>% dplyr::select(dplyr::all_of(sub_cols))
    }
    if (!is.null(out_cols)) {
      out_len <- length(out_cols)
      dat_len <- names(x) %>% length()
      
      if (sub_len > dat_len) {
        fail_mssg <- glue::glue("sub_len={sub_len} > out_len={out_len}")
        stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
        return(NULL)
      }
      x <- x %>% purrr::set_names(out_cols)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Get/Write Spec Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (write_spec && data_type == "tib") {
      spec_str <- x %>% 
        summarise_all(class) %>% 
        gather(col_name, col_type) %>% 
        dplyr::pull(col_type) %>% stringr::str_sub(1,1) %>% 
        paste0(collapse = '')
      
      file_spec <- 
        readr::read_csv(
          readr::format_csv(x, col_names = TRUE), 
          col_types = spec_str ) %>%
        readr::spec()
      
      readr::write_rds( file_spec, spec_rds )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (file_type == "line") {
      readr::write_lines( x = x, file = file, append = append )
    } else if (file_type == "csv") {
      readr::write_csv( x = x, file = file, append = append )
    } else if (file_type == "tsv") {
      readr::write_tsv( x = x, file = file, append = append )
    } else if (file_type == "rds") {
      readr::write_rds( x = x, file = file, compress = compress )
    } else {
      stop(cat(glue::glue(
        "{error} Invalid file_type={file_type}!","{error} Exiting...{RET2}") ) )
      return(NULL)
    }
    
    if (vb >= vt+1)
      cat(glue::glue("{mssg} Done writing data file = '{file}'.{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Set Permissions::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # permissions = "0777"
    if (!is.null(permissions)) {
      Sys.chmod(paths = file, mode = permissions)
      if (vb >= vt+1)
        cat(glue::glue("{mssg} Set permissions = '{permissions}'.{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Write Done File::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (done) {
      cmd <- glue::glue("touch {end_txt}")
      if (vb >= vt) cat(glue::glue("{mssg} Running cmd='{cmd}'...{RET}"))
      cmd_ret <- base::system(cmd)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (purrr::is_vector(x)) ret_cnt <- length(x)
    if (purrr::is_list(x)) ret_cnt <- length(x)
    if (base::is.data.frame(x) || tibble::is_tibble(x)) {
      ret_key <- glue::glue("final-ret-tib")
      ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                            vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    }
  })
  etime <- as.double(ftime[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(ftime,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_cnt
}

safe_read = function( file,
                      
                      skip  = 0,
                      head  = NULL,
                      type  = NULL,
                      clean = TRUE,
                      fresh = FALSE,
                      guess_max = 1000,
                      
                      use_spec = TRUE,
                      spec_rds = NULL,
                      spec_col = NULL,
                      has_head = TRUE,
                      write_spec = TRUE,
                      over_write = FALSE,
                      spec_prefernce = "rds",
                      
                      vb=0, vt=4, tc=1, tt=NULL,
                      fun_tag='safe_read') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  # First make sure an tibble wasn't passed in::
  if ( !purrr::is_character( file ) ) return( file )
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( is.null(spec_rds) )
    spec_rds <- stringr::str_remove(file, '.[a-zA-Z]+(.gz)$') %>%
    paste('spec.rds', sep = '.')
  
  if ( is.null(head) ) head <- "Starting"
  if ( vb>=vt ) cat(glue::glue("{mssg} {head}...{RET}"))
  if ( vb>=vt+2 ) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}  Verbose Parameters::{RET}"))
    cat(glue::glue("{mssg}               vb = '{vb}'.{RET}"))
    cat(glue::glue("{mssg}               vt = '{vt}'.{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}             skip = '{skip}'.{RET}"))
    cat(glue::glue("{mssg}             type = '{type}'.{RET}"))
    cat(glue::glue("{mssg}            clean = '{clean}'.{RET}"))
    cat(glue::glue("{mssg}            fresh = '{fresh}'.{RET}"))
    cat(glue::glue("{mssg}        guess_max = '{guess_max}'.{RET}"))
    cat(glue::glue("{mssg}         has_head = '{has_head}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}             file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}         use_spec = '{use_spec}'.{RET}"))
    cat(glue::glue("{mssg}         spec_rds = '{spec_rds}'.{RET}"))
    cat(glue::glue("{mssg}         spec_col = '{spec_col}'.{RET}"))
    cat(glue::glue("{mssg}       write_spec = '{write_spec}'.{RET}"))
    cat(glue::glue("{mssg}   spec_prefernce = '{spec_prefernce}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  if ( fresh && file.exists(spec_rds) ) unlink( spec_rds )
  
  errs_mssg <- glue::glue("Input file is null or length zero")
  if ( is.null(file) || length(file) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return( NULL )
  
  errs_mssg <- glue::glue("File: {file} does NOT exist")
  if ( !file.exists(file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return( NULL )
  
  if (use_spec && is.null(spec_col) && !file.exists(spec_rds) ) {
    warn_mssg <- glue::glue("Failed to find spec_rds={spec_rds}! ",
                            "Will not attempt to use spec columns")
    if ( vb >= vt ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    use_spec <- FALSE
  }
  if ( use_spec && !is.null(spec_col) && !purrr::is_list(spec_col) ){
    warn_mssg <- glue::glue("spec_col is not of type list! ",
                            "Will not use spec columns (caught in purrr ",
                            "is_list() )")
    if ( vb >= vt ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    use_spec <- FALSE
  }
  if ( use_spec && !is.null(spec_col) &&  base::typeof(spec_col) != "list" ) {
    warn_mssg <- glue::glue("spec_col is not of type list! ",
                            "Will not use spec columns (caught in base ",
                            "typeof() )")
    if ( vb >= vt ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    use_spec <- FALSE
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Guess File Type::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (vb >= vt+6)
      cat(glue::glue("{mssg} Guessing Stage; initial type = '{type}'{RET}"))
    
    file_type <- type
    if (is.null(file_type)) 
      file_type <- guess_file_extension( file, vb=vb,vt=vt+6,tc=tc+1 )
    
    if (vb >= vt+6)
      cat(glue::glue("{mssg} Guessing Stage; file type = '{file_type}'{RET}"))
    
    if (file_type != "rds" && is.null(type) )
      file_type <- guess_file_delimiter( file, vb=vb,vt=vt+6,tc=tc+1 )
    
    if (vb >= vt+6)
      cat(glue::glue("{mssg} Guessing Stage; final type = '{file_type}'{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Load File Specifications::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( ( is.null(spec_col) || spec_prefernce == "rds" ) && 
         ( !is.null(spec_rds) && file.exists(spec_rds) ) ) {
      
      if (vb >= vt+6)
        cat(glue::glue("{mssg} Reading spec_rds = {spec_rds}'...{RET}"))
      spec_col <- readr::read_rds(spec_rds)
      
      if (vb >= vt+8) print(spec_col)
    }
    
    if ( use_spec ) {
      
      errs_mssg <- glue::glue("For data file = {file}; Provided ",
                              "Specifications File: {spec_rds} does NOT ",
                              "have required $cols sub fields")
      if ( is.null(spec_col$cols) ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
      
      if (vb >= vt+4)
        cat(glue::glue("{mssg} Will use column specs...{RET}"))
      if (vb >= vt+8) print(spec_col)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Read Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (file_type == "line") {
      ret_tib <- suppressMessages(suppressWarnings( 
        readr::read_lines(file = file, skip = skip ) ) )
      
    } else if (file_type == "csv" || file_type == COM) {
      
      if ( !is.null(spec_col) && has_head ) {
        if (vb >= vt+6)
          cat(glue::glue("{mssg} Will use column specs for only ",
                         "variable types...{RET}"))
        if (vb >= vt+8) print(spec_col)
        
        ret_tib <- suppressMessages(suppressWarnings( 
          readr::read_csv( file=file, col_types = spec_col, skip = skip ) ) )
      } else if ( !is.null(spec_col) && !has_head ) {
        
        if (vb >= vt+6)
          cat(glue::glue("{mssg} Will use column specs for both names and ",
                         "variable types...{RET}"))
        if (vb >= vt+8) print(spec_col)
        
        ret_tib <- suppressMessages(suppressWarnings(
          readr::read_csv( file=file, 
                           col_names = names(spec_col$cols), 
                           col_types = spec_col,
                           skip = skip ) ) )
      } else {
        
        if (vb >= vt+6)
          cat(glue::glue("{mssg} Will NOT use column specs...{RET}"))
        
        ret_tib <- suppressMessages(suppressWarnings( 
          readr::read_csv( file=file, guess_max=guess_max, skip = skip ) ) )
      }
      
    } else if (file_type == "tsv" || file_type == TAB || file_type == " ") {
      
      if ( !is.null(spec_col) && has_head ) {
        ret_tib <- suppressMessages(suppressWarnings(
          readr::read_tsv( file=file, col_types = spec_col, skip = skip ) ) )
      } else if ( !is.null(spec_col) && !has_head ) {
        ret_tib <- suppressMessages(suppressWarnings(
          readr::read_tsv( file=file, 
                           col_names = names(spec_col$cols), 
                           col_types = spec_col, 
                           skip = skip ) ) )
      } else {
        ret_tib <- suppressMessages(suppressWarnings(
          readr::read_tsv( file=file, guess_max=guess_max, skip = skip ) ) )
      }
      
    } else if (file_type == "rds") {
      
      ret_tib <- readr::read_rds(file)
      
    } else {
      fail_mssg <- glue::glue("Invalid file_type='{file_type}'")
      stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
      return(NULL)
    }
    clean_key <- glue::glue("cleaning-ret-tib")
    if (clean) 
      ret_tib <- clean_tib( ret_tib, fun_tag = fun_tag, name = clean_key,
                            vb=vb,vt=vt+6,tc=tc+1 )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Get/Write Spec Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    data_type = "unk"
    if ( tibble::is_tibble(ret_tib) ) data_type = "tib"
    if (write_spec && is.null(spec_col) && data_type == "tib") {
      spec_str <- ret_tib %>% 
        summarise_all(class) %>% 
        gather(col_name, col_type) %>% 
        dplyr::pull(col_type) %>% 
        stringr::str_sub(1,1) %>% 
        paste0(collapse = '')
      
      spec_col <- 
        readr::read_csv(
          readr::format_csv(ret_tib, col_names = TRUE), 
          col_types = spec_str, 
          skip = skip ) %>%
        readr::spec()
      
      if ( vb >= vt+4 ) print(spec_col)
    }
    
    if ( write_spec && !is.null(spec_col) && 
         ( over_write || !file.exists(spec_rds) ) )
      readr::write_rds( spec_col, spec_rds )
    
    # file.mtime(file) > file.mtime(spec_rds) ) ) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  })
  etime <- as.double(ftime[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(ftime,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

valid_time_stamp = function( files,
                             out_dir = NULL,
                             vb=0, vt=6, tc=1, tt=NULL,
                             fun_tag='valid_time_stamp') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( vb >= vt+3 ) cat(glue::glue("{mssg} Starting...{RET}"))
  if ( vb >= vt+4 ) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   fun_tag={fun_tag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  ret_val <- TRUE
  
  ftime <- base::system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Check if Files Exist::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    file_cnt <- length(files)
    for (ii in c(1:file_cnt)) {
      cur_val <- TRUE
      if (!file.exists(files[ii])) cur_val <- FALSE
      if (vb>=vt+4)
        cat(glue::glue("{mssg} File({ii})={files[ii]} exist={cur_val}.{RET}"))
      
      ret_tib <- ret_tib %>% dplyr::bind_rows(
        tibble::tibble(Index = ii,
                       File = files[ii],
                       Valid = cur_val)
      )
      if (cur_val == FALSE) ret_val = FALSE
    }
    if (vb >= vt+5)
      cat(glue::glue("{mssg} Files exist={ret_val}.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Compare File Time Stamps::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (ret_val) {
      file_cnt <- file_cnt - 1
      for (ii in c(1:file_cnt)) {
        f1 <- files[ii]
        f2 <- files[ii+1]
        
        # t1 <- file.mtime(f1)
        # t2 <- file.mtime(f2)
        
        # Improved lookup time file.info() vs. file.mtime()
        #   extra_cols=TRUE forces R to use native C back-end code
        t1 <- base::file.info( f1, extra_cols = TRUE) %>%
          dplyr::pull(mtime)
        t2 <- base::file.info( f2, extra_cols = TRUE) %>%
          dplyr::pull(mtime)
        
        if ( vb >= vt+5 )
          cat(glue::glue("{mssg} ii={ii}, t1={t1}; f1={f1}{RET}",
                         "{mssg} ii={ii}, t2={t2}; f2={f2}{RET}"))
        
        cur_val <- TRUE
        if ( t1 > t2 ) cur_val <- FALSE
        ret_tib <- ret_tib %>% dplyr::bind_rows(
          tibble::tibble(Index = ii,
                         File = files[ii],
                         Valid = cur_val)
        )
        if ( cur_val == FALSE ) ret_val = FALSE
        
        if ( vb >= vt+5 )
          cat(glue::glue("{mssg} ii={ii}, Valid Order={ret_val}.{RET2}"))
      }
    }
    ret_cnt <- files %>% length()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( !is.null(out_dir) ) {
      safe_mkdir( out_dir )
      sum_csv <- file.path( out_dir, "valid_time_stamp.csv.gz" )
      safe_write( x = ret_tib, file = sum_csv, done = FALSE, write_spec = FALSE, 
                  append = FALSE, fun_tag = fun_tag, 
                  vb=vb,vt=vt+7,tc=tc+1,tt=tt )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, max = 0,
                          vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  })
  etime <- as.double(ftime[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(ftime,fun_tag)
  if ( vb >= vt+3 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_val
}

sort_files_by_date = function( path, 
                               pattern = NULL, 
                               by = 'mtime',
                               ret = "tib",
                               desc = FALSE,
                               
                               vb=0, vt=6, tc=1, tt=NULL,
                               fun_tag='sort_files_by_date') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}        by = '{by}'.{RET}"))
    cat(glue::glue("{mssg}       ret = '{ret}'.{RET}"))
    cat(glue::glue("{mssg}      path = '{path}'.{RET}"))
    cat(glue::glue("{mssg}   pattern = '{pattern}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if (!by %>% stringr::str_detect('^(m|a|c)time$')) {
    fail_mssg <- glue::glue("Argument `by` must be one of [ 'mtime', ",
                            "'atime', or 'ctime' ]")
    stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    return(NULL)
  }
  
  if ( !ret %in% c("tib", "vec", "list") ) {
    fail_mssg <- glue::glue("Argument `ret` must be one of [ 'tib', ",
                            "'vec', or 'list' ]")
    stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    return(NULL)
  }
  
  file_names <-
    
    # Gets all the file paths for a given path and pattern
    base::list.files( path = path, pattern = pattern, full.names = TRUE ) %>%
    
    # Turns into a one-column tibble (see below)
    tibble::tibble( file_names = . )
  
  ret_tib <-
    
    base::suppressWarnings(
      furrr::future_map_dfr(
        .x = file_names,
        .f = file.info,
        .progress = TRUE,
        extra_cols = FALSE # passed to file.info for performance boost
      )
    ) %>%
    
    # gets expanded file info, then select the last-modified time and last-accessed time
    select( mtime, atime, ctime ) %>%
    
    # reintroduce our original 'file_names'
    dplyr::bind_cols( file_names )
  
  ret_key <- glue::glue("Pre-sort-ret-tib(desc={desc})")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                        vb=vb,vt=vt+5,tc=tc+1,tt=tt )
  
  if (desc)
    ret_tib <- ret_tib %>%
    # arrange by descending time (depending on the users choice)
    dplyr::arrange(
      dplyr::desc(
        dplyr::case_when(
          (by == 'mtime') ~ mtime,
          (by == 'atime') ~ atime,
          (by == 'ctime') ~ ctime
        )
      )
    )
  else
    ret_tib <- ret_tib %>%
    # arrange by descending time (depending on the users choice)
    dplyr::arrange(
      dplyr::case_when(
        (by == 'mtime') ~ mtime,
        (by == 'atime') ~ atime,
        (by == 'ctime') ~ ctime
      )
    )
  
  
  if ( ret == "tib" ) {
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  } else if ( ret == "vec" ) {
    ret_tib <- ret_tib %>% dplyr::pull( file_names )
    
    ret_cnt <- length(ret_tib)
  } else if ( ret == "list" ) {
    ret_tib <- ret_tib %>% dplyr::pull( file_names )
    
    ret_fns <- ret_tib %>% 
      stringr::str_remove("\\.gz$") %>%
      stringr::str_remove("\\.[ct]sv$") %>%
      base::basename()
    
    ret_tib <- as.list(ret_tib)
    names(ret_tib) <- ret_fns
    
    ret_cnt <- length(ret_tib)
  }
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

get_file_list = function( x,
                          
                          clean     = TRUE,
                          suffix    = c("\\.gz", "\\.[a-z]+"),
                          pattern   = NULL,
                          validate  = TRUE,
                          recursive = FALSE, 
                          
                          vb=0, vt=6, tc=1, tt=NULL,
                          fun_tag='get_file_list') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p3  <- vb > vt + 3
  p4  <- vb > vt + 4
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET}"))
  if ( p2 ) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       clean = '{clean}'.{RET}"))
    cat(glue::glue("{mssg}      suffix = '{suffix}'.{RET}"))
    cat(glue::glue("{mssg}     pattern = '{pattern}'.{RET}"))
    cat(glue::glue("{mssg}    validate = '{validate}'.{RET}"))
    cat(glue::glue("{mssg}   recursive = '{recursive}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_dat <- NULL
  
  if ( purrr::is_vector(x) ) {
    
    file_vec <- x
    
  } else if ( purrr::is_character(x) &&
              dir.exists(x) ) {
    if ( is.null(pattern) ) {
      fail_mssg <- glue::glue( "When x is a directory you must provide a ",
                               "search pattern")
      stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
      return(NULL)
    }
    file_vec <- list.files( dir, 
                            pattern    = pattern, 
                            full.names = TRUE, 
                            recursive  = recursive)
    if ( length(file_vec) == 0 ) {
      fail_mssg <- glue::glue( "Failed to find {pattern} in direcotry = {dir}")
      stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
      return(NULL)
    }
    
  } else if ( purrr::is_character(x) &&
              file.exists(x) ) {
    
    file_vec <- c(x)
    
  } else {
    fail_mssg <- glue::glue( "x must be a vector, existing file or directory" )
    stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
    return(NULL)
  }
  
  if (validate)
    # Validate all file exist!
    for ( file in file_vec )
      if ( !file.exists(file) ) {
        fail_mssg <- glue::glue( "Failed to find {file}")
        stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
        return(NULL)
      }
  
  # Get Names::
  name_vec <- file_vec %>% base::basename()
  
  # Trim Names::
  if ( !is.null(suffix) && length(suffix) != 0 )
    for (s in suffix) name_vec <- stringr::str_remove( name_vec, paste0(s,"$") )
  
  # Clean Names::
  if (clean) name_vec <- name_vec %>% 
    stringr::str_replace_all("[^[:alnum:]]", "_")
  
  # Build List
  ret_dat <- stats::setNames( as.list(file_vec), name_vec )
  ret_cnt <- ret_dat %>% names() %>% length()
  
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

file_list = function( path = NULL,
                      subs = NULL,
                      file = NULL,
                      
                      prefix = NULL,
                      suffix = NULL,
                      
                      names = NULL,
                      dir_only = FALSE,
                      
                      unique = TRUE,
                      pattern = NULL,
                      recursive = FALSE,
                      ret_type  = "list",
                      subs_exists  = TRUE,
                      paths_exists = TRUE,
                      files_exists = TRUE,
                      
                      vb=0, vt=6, tc=1, tt=NULL,
                      fun_tag='file_list') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p3  <- vb > vt + 3
  p4  <- vb > vt + 4
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p3 ) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   path = '{path}'.{RET}"))
    cat(glue::glue("{mssg}   subs = '{subs}'.{RET}"))
    cat(glue::glue("{mssg}   file = '{file}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  paths <- NULL
  files <- NULL
  
  ret_cnt <- 0
  
  errs_mssg <- glue::glue("Both path and file can't be NULL")
  if ( is.null(path) && is.null(file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( purrr::is_character( path ) ) path <- split_str_to_vec( path, unique = unique )
  if ( purrr::is_character( subs ) ) subs <- split_str_to_vec( subs, unique = unique )
  if ( purrr::is_character( file ) ) file <- split_str_to_vec( file, unique = unique )
  
  # Check if directories exist::
  errs_mssg <- glue::glue("Some directory paths don't exist = {path}")
  if ( !is.null(path) ) flag <- FALSE %in% lapply(path, dir.exists) %>% unlist()
  if ( eflag ) cat(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  path_len <- length(path)
  subs_len <- length(subs)
  file_len <- length(file)
  
  if ( !is.null(path) && paths_exists ) {
    # Check if paths exist::
    errs_mssg <- glue::glue("Some file paths don't exist = {path}")
    if ( !is.null(path) ) eflag <- FALSE %in% lapply(path,dir.exists) %>% unlist()
    if ( eflag ) cat(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
  }
  
  if ( p4 ) {
    cat(glue::glue("{mssg} path = {path_len}.{RET}"))
    cat(glue::glue("{mssg} subs = {subs_len}.{RET}"))
    cat(glue::glue("{mssg} file = {file_len}.{RET}"))
  }
  
  # Merge paths and subs
  if ( !is.null(path) && !is.null(subs) ) {
    if ( path_len == subs_len )
      for (ii in c(1:path_len) ) paths[ii] <- file.path(path[ii], subs[ii])
    else if ( path_len == 1 && subs_len != 0 )
      for (ii in c(1:subs_len) ) paths[ii] <- file.path(path, subs[ii])
    else if ( path_len != 0 && subs_len == 1 )
      for (ii in c(1:path_len) ) paths[ii] <- file.path(path[ii], subs)
    
    path <- paths
    
    if ( subs_exists ) {
      # Check if paths exist::
      errs_mssg <- glue::glue("Some sub-directory paths don't exist = {path}")
      if ( !is.null(path) ) eflag <- FALSE %in% lapply(path,dir.exists) %>% unlist()
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
  }
  
  # Merge paths and files
  if ( !is.null(path) && !is.null(file) ) {
    if ( path_len == file_len )
      for (ii in c(1:path_len) ) files[ii] <- file.path(path[ii], file[ii])
    else if ( path_len == 1 && file_len != 0 )
      for (ii in c(1:file_len) ) files[ii] <- file.path(path, file[ii])
    
  } else if ( !is.null(path) && is.null(file) && !dir_only ) {
    files <- list.files( path = path, pattern = pattern, full.names = TRUE, 
                         recursive = recursive )
  } else if ( is.null(path) && !is.null(file) ) {
    files <- file
  } else {
    # Don't think this should ever come up...
    #
    # errs_mssg <- glue::glue("Both path and file can't be NULL")
    # if ( is.null(path) && is.null(pattern) ) eflag <- TRUE
    # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    # if ( eflag ) return(NULL)
  }
  
  if ( dir_only ) files <- path
  
  errs_mssg <- glue::glue("Failed to find any files")
  if ( is.null(files) || length(files) == 0) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( !dir_only && files_exists ) {
    # Check if files exist::
    errs_mssg <- glue::glue("Some file paths don't exist = {files}")
    if ( !is.null(files) ) eflag <- FALSE %in% lapply(files,file.exists) %>% unlist()
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
  }
  
  if ( ret_type == "list" ) {
    
    if ( is.null(names) ) names <- base::basename( files )
    
    # Clean up names::
    if ( !is.null(prefix) )
      for (p in prefix) names <- stringr::str_remove(names, p)
    
    if ( !is.null(suffix) )
      for (s in suffix) names <- stringr::str_remove(names, s)
    
    # Create List::
    files <- as.list(files)
    names(files) <- names
    
    # Validate List::
    errs_mssg <- glue::glue("Failed to build file LIST")
    if ( !purrr::is_list(files) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
  } else if ( ret_type == "vec") {
    # Do nothing...
  } else {
    files <- files %>% paste( collapse = ret_type )
  }
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( files, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  files
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Standard Tibble Functions::
#
#  print_tib()
#  clean_tib()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tib_group_split = function( tib,
                            vec,
                            sep = ",",
                            vb=0, vt=6, tc=1, tt=NULL,
                            fun_tag='tib_group_split')
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
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( vec %>% length() == 0 ) {
    warn_mssg <- glue::glue("Group Vector is size zero returning original tib")
    if ( p2 ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    return( tib )
  }
  
  ret_dat <- tib %>% 
    tidyr::unite( NEW_GROUP, dplyr::all_of(vec), sep = sep, remove = FALSE ) %>%
    # dplyr::mutate( NEW_GROUP = paste( , sep = sep ) ) %>%
    dplyr::select( NEW_GROUP, dplyr::everything() ) %>%
    split( .$NEW_GROUP )
  
  # errs_mssg <- glue::glue("ERROR_MESSAGE")
  # if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_dat, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

print_sum = function( tib, vec,
                      
                      vb=0, vt=3, tc=1, tt=NULL,
                      fun_tag='print_sum') {
  
  p0  <- vb > vt + 0
  
  sum_tib <- NULL
  sum_tib <- tib %>% 
    dplyr::group_by( dplyr::across( dplyr::all_of(vec) ) ) %>%
    dplyr::summarise( Count=n(), .groups = "drop")
  
  if ( p0 ) print( sum_tib, n=base::nrow(sum_tib) )
  
  sum_tib
}

print_tib = function( x, 
                      name = "",
                      max  = 3,
                      
                      vb=0, vt=4, tc=1, tt=NULL,
                      fun_tag='print_tib') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  p0  <- vb > vt + 0
  
  ret_dat <- NULL
  
  # Get number of elements::
  ret_cnt <- 0
  col_cnt <- 0
  
  dat_type <- "unk"
  if ( base::typeof(x) == "S4" ) {
    ret_cnt  <- x %>% names() %>% length()
    dat_type <- "S4"
  } else if ( base::is.data.frame(x) || tibble::is_tibble(x) ) {
    ret_cnt <- base::nrow( x )
    col_cnt  <- base::ncol(x)
    dat_type <- "tib"
    
    if ( max == 0 )
      ret_dat <- x
    else 
      ret_dat <- head( x, n = max )
    
  } else if ( purrr::is_list(x) ) {
    ret_cnt <- length( names( x ) )
    dat_type <- "list"
    
    if ( max == 0 )
      ret_dat <- names( x )
    else
      ret_dat <- names( x ) %>% head( n = max )
    
  } else if ( base::is.matrix(x) ) {
    ret_cnt  <- base::nrow(x)
    col_cnt  <- base::ncol(x)
    dat_type <- "mat"
    
    if ( max == 0 )
      ret_dat <- x
    else
      ret_dat <- x[ 1:min(max,col_cnt), 1:min(max,ret_cnt) ]
    
    # ret_dat <- x[ 1:min(max,ret_cnt), 1:min(max,col_cnt) ]
    
  } else if ( purrr::is_vector(x) ) {
    ret_cnt <- length( x )
    dat_type <- "vec"
    
    if ( max == 0 )
      ret_dat <- x
    else
      ret_dat <- head( x, n = max )
  }
  
  if ( p0 ) {
    cat(glue::glue("{mssg} {name}; row count({dat_type}) = [{ret_cnt}] ",
                   "row count({dat_type}) = [{col_cnt}] ",
                   "print = {max}::{RET}"))
    if ( dat_type == "S4" ) {
      print( x )
    } else if ( dat_type == "tib" || dat_type == "mat" ) {
      print( ret_dat )
    } else {
      cat( glue::glue("{mssg} '{ret_dat}'.{RET}") )
    }
    cat(glue::glue("{RET}{mssg}{BRK}{RET2}") )
  }
  
  ret_cnt
}

clean_tib = function(tib,
                     name = NULL,
                     vb=0, vt=6, tc=1, tt=NULL,
                     fun_tag='clean_tib') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ret_tib <- tib %>%
    dplyr::select(where(~sum(!is.na(.x)) > 0)) %>%
    utils::type.convert(as.is=TRUE) %>%
    dplyr::mutate(across(where(is.factor), as.character) )
  
  # NOTE:: Don't use the readr version. It does something different than 
  #  the utils version...
  #  readr::type_convert() %>% 
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Standard String Functions::
#
#  split_str_to_vec()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

split_str_to_vec = function(x, 
                            del = COM, 
                            unique = TRUE,
                            vb=0, vt=6, tc=1, tt=NULL,
                            fun_tag='split_str_to_vec') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p3  <- vb > vt + 3
  p4  <- vb > vt + 4
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET}"))
  if ( p2 ) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}        x={x}.{RET}"))
    cat(glue::glue("{mssg}      del={del}.{RET}"))
    cat(glue::glue("{mssg}   unique={unique}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_vec <- NULL
  
  if ( !is.null(x) ) ret_vec <- 
    stringr::str_split(x, pattern=del, simplify=TRUE) %>% as.vector()
  
  if (unique) ret_vec <- ret_vec %>% unique()
  
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_vec
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                General Program Initialization Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_init = function( name,
                         opts, opt_reqs = NULL,
                         pars, par_reqs = NULL,
                         rcpp = 0,
                         long_out_dir = TRUE,
                         vb=0, vt=3, tc=1, tt=NULL,
                         fun_tag='program_init')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0 <- vb > vt + 0
  p1 <- vb > vt + 1
  p2 <- vb > vt + 2
  p3 <- vb > vt + 3
  p4 <- vb > vt + 4
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p2 ) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   fun_tag={fun_tag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Validate Required Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  if ( !is.null(par_reqs) ) {
    par_names   <- names(pars)
    par_int_vec <- intersect( par_names, par_reqs )
    
    if ( !setequal( par_int_vec, par_reqs ) ) {
      mis_vec <- setdiff( par_int_vec, par_reqs )
      mis_cnt <- length( mis_vec )
      
      cat(glue::glue("{errs}: Missing required parameters = {mis_vec}.{RET}"))
      for ( cur_par in par_reqs ) {
        if ( !cur_par %in% par_names ) cat(glue::glue(
          "{errs}: {TAB}Missing required parameter = '{cur_par}'."))
      }
      cat("\n\n")
      
      fail_mssg <- glue::glue("Failed to find all required parameters! ",
                              "Missing {mis_cnt} required parameters")
      stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
      return(NULL)
    }
  }
  if ( !is.null(opt_reqs) ) {
    opt_names   <- names(opts)
    opt_int_vec <- intersect( opt_names, opt_reqs )
    
    if ( !setequal( opt_int_vec, opt_reqs ) ) {
      mis_vec <- setdiff( opt_int_vec, opt_reqs )
      mis_cnt <- length( mis_vec )
      
      cat(glue::glue("{errs}: Missing required options = {mis_vec}.{RET}"))
      for ( cur_opt in opt_reqs ) {
        if ( !cur_opt %in% opt_names ) cat(glue::glue(
          "{errs}: {TAB}Missing required option = '{cur_opt}'."))
      }
      cat("\n\n")
      
      fail_mssg <- glue::glue("Failed to find all required options! ",
                              "Missing {mis_cnt} required options")
      stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
      return(NULL)
    }
  }
  
  # TBD:: Redundant:: Will remove after testing!!!
  # stopifnot(!is.null(opts[['out_path']]))
  stopifnot(!is.null(pars[['prgm_tag']]))
  stopifnot(!is.null(pars[['exe_path']]))
  
  # TBD:: Redundant:: Will remove after testing!!!
  # if (!is.null(opt[['verbose']])) verbose <- opts$verbose
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                              Load Libraries::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # NOTE: This is already done now in the program executable script::
  #
  #  src_cnt <- load_libraries(pars=pars, rcpp=rcpp, vb,vt=vt+4,tc=tc+1 )
  if ( rcpp > 0 ) {
    opts$rcpp <- TRUE
    
    if ( rcpp > 0 ) {
      # Formerly bisulfite_functions.cpp = assoc_stack_functions.cpp
      # pars$src_1_rcpp_path <- file.path(pars$src_path, 'Rcpp/assoc_stack_functions.cpp')
      pars$src_1_rcpp_path <- file.path(pars$src_path, 'Rcpp/bisulfite_functions.cpp')
      # pars$src_1_rcpp_path <- file.path(pars$src_path, 'Rcpp/idat_functions.cpp')
      if ( p1 )
        cat(glue::glue("{mssg} Will load source file = {pars$src_1_rcpp_path}.{RET}"))
      ret_cnt <- ret_cnt + 1
      Rcpp::sourceCpp( pars$src_1_rcpp_path )
    }
    
    if ( rcpp > 1 ) {
      pars$src_1_rcpp_path <- file.path(pars$src_path, 'Rcpp/idat_functions.cpp')
      if ( p1 )
        cat(glue::glue("{mssg} Will load source file = {pars$src_1_rcpp_path}.{RET}"))
      ret_cnt <- ret_cnt + 1
      Rcpp::sourceCpp( pars$src_1_rcpp_path )
    }

    if ( rcpp > 2 ) {
      pars$src_1_rcpp_path <- file.path(pars$src_path, 'Rcpp/idat_pair_functions.cpp')
      if ( p1 )
        cat(glue::glue("{mssg} Will load source file = {pars$src_1_rcpp_path}.{RET}"))
      ret_cnt <- ret_cnt + 1
      Rcpp::sourceCpp( pars$src_1_rcpp_path )
    }
    
    if ( FALSE && rcpp > 3 ) {
      pars$src_2_rcpp_path <- file.path(pars$src_path, 'Rcpp/infinium_bisulfite_functions.cpp')
      if ( p1 )
        cat(glue::glue("{mssg} Will load source file = {pars$src_1_rcpp_path}.{RET}"))
      ret_cnt <- ret_cnt + 1
      Rcpp::sourceCpp( pars$src_2_rcpp_path )
    }
    
    if ( FALSE && rcpp > 4 ) {
      pars$src_3_rcpp_path <- file.path(pars$src_path, 'Rcpp/cpg_loci_variation.cpp')
      if ( p1 )
        cat(glue::glue("{mssg} Will load source file = {pars$src_2_rcpp_path}.{RET}"))
      ret_cnt <- ret_cnt + 1
      Rcpp::sourceCpp( pars$src_3_rcpp_path )
    }
  } else {
    opts$rcpp <- FALSE
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Set Common Defaults:: docker, etc.
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # TBD:: Update these values after building new docker branch::
  #   NOTE:: Should actually be done in program_options() NOT HERE
  #    Although we should check they're not null here...
  #
  if ( is.null(opts[["image_key"]]) )
    opts$image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
  
  if ( is.null(opts[["image_ver"]]) )
    opts$image_ver <- "v.1.25"
  
  if ( is.null(opts[["doc_shell"]]) )
    opts$doc_shell <- "run_improbe.sh"
  
  if ( is.null(opts[["doc_image"]]) )
    opts$doc_image <- glue::glue("{opts$image_key}.{opts$image_ver}")
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Build Directories::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( p0 ) cat(glue::glue("{mssg} Output Directory = {opts$out_path}...{RET}"))
  
  if (long_out_dir) {
    opts$out_path <- file.path( opts$out_path, name )
    if (!is.null( opts$run_name ) ) 
      opts$out_path <- file.path( opts$out_path, opts$run_name )
  }
  
  if ( !is.null(opts[['clean']]) && opts$clean && dir.exists(opts$out_path ) ) {
    rm_files <- list.files(opts$out_path, full.names = TRUE, recursive = TRUE, 
                           include.dirs = FALSE)
    rm_dirs <- list.dirs(opts$out_path, full.names = TRUE, recursive = TRUE)
    
    if ( p0 ) cat(glue::glue("{mssg} Removing file(s) = '{rm_files}'...{RET}"))
    unlink(rm_files)
    
    if ( p0 )
      cat(glue::glue("{mssg} Removing Directories(s) = '{rm_dirs}'...{RET}"))
    unlink(rm_dirs, recursive = TRUE)
  }
  safe_mkdir( opts$out_path )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Program Start Time Stamp::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts$time_org_txt <- 
    file.path(opts$out_path,paste(name,'time_stamp-org.txt', sep='.'))
  
  if (!file.exists(opts$time_org_txt) || 
      ( !is.null(opts[['clean']]) && opts$clean ) )
    readr::write_lines(x=date(), file=opts$time_org_txt,sep='\n', append=FALSE)
  
  opts$opt_csv  <- file.path(
    opts$out_path, paste(pars$prgm_tag,'program-options.csv', sep='.') )
  opts$par_csv  <- file.path(
    opts$out_path, paste(pars$prgm_tag,'program-parameters.csv', sep='.') )
  opts$time_csv <- file.path(
    opts$out_path, paste(pars$prgm_tag,'time-tracker.csv.gz', sep='.') )
  
  if (file.exists(opts$opt_csv))  unlink(opts$opt_csv)
  if (file.exists(opts$par_csv))  unlink(opts$par_csv)
  if (file.exists(opts$time_csv)) unlink(opts$time_csv)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                          Program Command Shell::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (pars$prgm_tag==name) cmd_shell_name <- name
  
  opts_tib <- opts %>%
    unlist(recursive = TRUE) %>%
    dplyr::bind_rows() %>% 
    tidyr::gather("Option", "Value") %>%
    dplyr::distinct()
  
  if ( p4 ) opts_tib %>% base::print(n=base::nrow(opts_tib) )
  
  pars_tib <- pars %>%
    unlist(recursive = TRUE) %>%
    dplyr::bind_rows() %>% 
    tidyr::gather("Params", "Value") %>%
    dplyr::distinct()
  
  if ( p4 ) pars_tib %>% base::print(n=base::nrow(pars_tib) )
  
  pars$cmd_shell <- 
    file.path(opts$out_path,paste(cmd_shell_name,'command.sh', sep='.'))
  pars$cmd_str <- 
    opts_to_command( opts = opts_tib,
                     exe  = pars$exe_path,
                     pre  = opts$Rscript,
                     file = pars$cmd_shell, 
                     vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  ret_cnt <- names(opts) %>% length()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                                  Done::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat$opt <- opts
  ret_dat$par <- pars
  ret_dat$opt_tib <- opts_tib
  ret_dat$par_tib <- pars_tib
  
  ret_dat
}

program_done = function(opts, pars, precision=3,
                        
                        vb=0, vt=6, tc=1, tt=NULL,
                        fun_tag='program_done') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p3  <- vb > vt + 3
  p4  <- vb > vt + 4
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  
  ret_cnt <- 0
  
  opts_tib <- opts %>%
    unlist(recursive = TRUE) %>%
    dplyr::bind_rows() %>% 
    tidyr::gather("Option", "Value") %>%
    dplyr::distinct()
  
  pars_tib <- pars %>%
    unlist(recursive = TRUE) %>%
    dplyr::bind_rows() %>% 
    tidyr::gather("Params", "Value") %>%
    dplyr::distinct()
  
  time_tib <- tt$time %>% 
    dplyr::mutate_if(is.numeric, list(round), 4)
  
  if ( !is.null(tt) ) {
    stamp_str <- Sys.time() %>% 
      stringr::str_replace_all(" ","_") %>% 
      stringr::str_replace_all(":","_") %>% 
      stringr::str_replace_all("-","_")
    
    rtime_csv <- file.path( opts$out_path, paste0(opts$run_name,".run_time_summary.",stamp_str,".csv.gz") )
    rtime_sum <- tt$time %>% 
      dplyr::group_by(Method) %>% 
      dplyr::summarise( Cnt=n(), 
                        Avg=mean(elapsed), 
                        Med=median(elapsed), 
                        Sum=sum(elapsed) ) %>% 
      dplyr::mutate( dplyr::across( where(is.numeric), ~round(.x, precision ) ) ) %>%
      dplyr::arrange( Sum )
    if ( p1 ) rtime_sum %>% print( n=base::nrow(rtime_sum) )
    readr::write_csv( x = rtime_sum, file = rtime_csv )
  }

  readr::write_csv(opts_tib, opts$opt_csv)
  readr::write_csv(pars_tib, opts$par_csv)
  readr::write_csv(time_tib, opts$time_csv)
  
  ret_key <- glue::glue("final-time-tib")
  ret_cnt <- print_tib( time_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  sysTime <- Sys.time()
  cat(glue::glue("{RET}[{pars$prgm_tag}]: Finished! (time={sysTime}){RET2}"))
  
  return(0)
}

opts_to_command = function(opts, 
                           exe,
                           pre  = NULL, 
                           rm   = NULL,
                           add  = NULL, 
                           file = NULL,
                           key  = "Option",
                           val  = "Value",
                           cluster=TRUE,
                           
                           vb=0, vt=6, tc=1, tt=NULL,
                           fun_tag='opts_to_command') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p3  <- vb > vt + 3
  p4  <- vb > vt + 4
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET}"))
  if ( p2 ) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   fun_tag={fun_tag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ret_val  <- NULL
  add_str  <- ''
  opt_str  <- ''
  bool_str <- ''
  key_sym <- key %>% as.character() %>% rlang::sym()
  val_sym <- val %>% as.character() %>% rlang::sym()
  
  # Options:: Exclude
  if (!is.null(rm)) opts <- opts %>% dplyr::filter(! (!!key_sym %in% rm) )
  
  # Handel Boolean Variables::
  bool <- NULL
  
  # TBD:: Test the below to make sure it works::
  bool <- opts %>% dplyr::filter(!!val == "TRUE")
  opts <- opts %>% dplyr::filter(!!val != "TRUE")
  opts <- opts %>% dplyr::filter(!!val != "FALSE")
  
  if (!is.null(bool)) {
    bool_str <- bool %>% 
      dplyr::mutate(!!key_sym := stringr::str_c('--',!!key_sym)) %>%
      dplyr::pull(!!key_sym) %>% 
      stringr::str_c(collapse=" ")
    if ( p4 ) cat(glue::glue("{mssg} bool_str='{bool_str}'.{RET2}"))
  }
  
  # Merge Options::
  if (!is.null(opts)) {
    opt_str <- opts %>% dplyr::arrange(!!key_sym) %>%
      dplyr::mutate(!!val_sym := paste0('"',!!val_sym,'"')) %>%
      tidyr::unite(Param, !!key_sym, !!val_sym, sep='=') %>%
      dplyr::mutate(Param=stringr::str_c('--',Param)) %>% 
      dplyr::pull(Param) %>%
      stringr::str_c(collapse=" ")
    if ( p4 ) cat(glue::glue("{mssg} opt_str='{opt_str}'.{RET2}"))
    
    # Second Removal Attempt of rm fields::
    if (!cluster) {
      opt_str <- opt_str %>% stringr::str_remove('--cluster')
      if ( p4 ) 
        cat(glue::glue("{mssg} opt_str(-cluster)='{opt_str}'.{RET2}"))
    }
  }
  
  # Options:: Add
  if (!is.null(add)) {
    add_str <- add %>% tidyr::unite(Param, !!key_sym, !!val_sym, sep='=') %>%
      dplyr::mutate(Param=stringr::str_c('--',Param)) %>% 
      dplyr::pull(Param) %>%
      stringr::str_c(collapse=" ")
    if ( p4 ) cat(glue::glue("{mssg} add_str='{add_str}'.{RET2}"))
  }
  
  # Add Executable and Join Options::
  cmd <- ''
  if (!is.null(pre) && length(pre)!=0) cmd <- pre
  cmd <- stringr::str_c(cmd,exe,opt_str,add_str,bool_str, sep=' ')
  if ( p4 ) cat(glue::glue("{mssg} cmd='{cmd}'.{RET2}"))
  
  ret_val <- cmd
  if (!is.null(file)) {
    dir <- base::dirname(file)
    if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    base::unlink(file)
    
    if ( p1 ) 
      cat(glue::glue("{mssg} Writing program shell={file}...{RET}"))
    readr::write_lines(cmd, file=file)
    Sys.chmod(file, mode="0777")
    
    ret_val <- file
  }
  
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_val
}

source_functions = function(pars,
                            rcpp = FALSE,
                            vb=0, vt=6, tc=1, tt=NULL,
                            fun_tag='source_functions') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p3  <- vb > vt + 3
  p4  <- vb > vt + 4
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   fun_tag={fun_tag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  
  source_files <- 
    base::list.files(pars$src_path, full.names = TRUE) %>% 
    file.path("functions") %>% 
    base::list.files(pattern='\\.R$', full.names=TRUE)
  ret_cnt <- length(source_files)
  
  if ( p1 )
    cat(glue::glue("{mssg} Will load source file = {source_files}.{RET}"))
  
  # lapply(source_files, base::source)
  
  for ( src in source_files ) {
    if ( p1 )
      cat(glue::glue("{mssg} {TAB} Sourcing = {src}.{RET}"))
    base::source( src )
  }
  pars$src_file_cnt <- ret_cnt
  
  if ( p1 ) cat(glue::glue("{mssg} Sourced {pars$src_file_cnt} files.{RET}"))
  
  if ( rcpp ) {
    pars$src_1_rcpp_path <- file.path(pars$src_path, 'Rcpp/infinium_bisulfite_functions.cpp')
    if ( p1 )
      cat(glue::glue("{mssg} Will load source file = {pars$src_1_rcpp_path}.{RET}"))
    ret_cnt <- ret_cnt + 1
    Rcpp::sourceCpp( pars$src_1_rcpp_path )
    
    pars$src_2_rcpp_path <- file.path(pars$src_path, 'Rcpp/cpg_loci_variation.cpp')
    if ( p1 )
      cat(glue::glue("{mssg} Will load source file = {pars$src_2_rcpp_path}.{RET}"))
    ret_cnt <- ret_cnt + 1
    Rcpp::sourceCpp( pars$src_2_rcpp_path )
  }
  pars$src_file_cnt <- ret_cnt
  
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  pars
}

params_check = function( pars, args, prgm_aux_check = TRUE,
                         
                         vb=0, vt=4, tc=1, tt=NULL,
                         fun_tag='params_check') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  p3  <- vb > vt + 3
  p4  <- vb > vt + 4
  p6  <- vb > vt + 6
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb>=vt+2) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   args={args}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Validate Input Parameters::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Input: args is empty or invalid")
  if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("Input: pars is empty or invalid")
  if ( is.null(pars) || length(pars) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("Input: run_mode is empty or invalid")
  if ( is.null(pars$run_mode) || length(pars$run_mode) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  if ( p6 ) cat(glue::glue("{mssg} Params: run_mode='{pars$run_mode}'{RET}"))
  
  errs_mssg <- 
    glue::glue("Input: source directory (src_path) is empty or invalid")
  if ( is.null(pars$src_path) || length(pars$src_path) == 0 || 
       !dir.exists(pars$src_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  if ( p6 ) cat(glue::glue("{mssg} Params: src_path='{pars$src_path}'{RET}"))
  
  pars$dat_path <- 
    file.path( base::dirname( pars$src_path ) %>% base::dirname(), "dat" )
  errs_mssg <- glue::glue("Input: data directory (dat_path) is empty or ",
                          "invalid or does not exist")
  if ( is.null(pars$dat_path) || length(pars$dat_path) == 0 || 
       !dir.exists(pars$dat_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  if ( p6 ) cat(glue::glue("{mssg} Params: dat_path='{pars$dat_path}'{RET}"))
  
  pars$aux_path <- file.path( pars$dat_path, "auxiliary" )
  errs_mssg <- 
    glue::glue("Input: auxiliary directory (aux_path) is empty or invalid")
  if ( is.null(pars$aux_path) || length(pars$aux_path) == 0 || 
       !dir.exists(pars$aux_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  if ( p6 ) cat(glue::glue("{mssg} Params: aux_path='{pars$aux_path}'{RET}"))
  
  if (pars$run_mode == 'RStudio') {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Validate Rstudio Parameters::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    errs_mssg <- 
      glue::glue("Input: scripts directory (prgm_dir) is empty or invalid")
    if ( is.null(pars$prgm_dir) || length(pars$prgm_dir) == 0 ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    errs_mssg <- 
      glue::glue("Input: progarm name (prgm_tag) is empty or invalid")
    if ( is.null(pars$prgm_tag) || length(pars$prgm_tag) == 0 ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    pars$scr_path <- file.path( pars$src_path, pars$prgm_dir )
    pars$exe_path <- file.path( pars$scr_path, paste0( pars$prgm_tag,".R" ) )
    
    errs_mssg <- 
      glue::glue("Input: scripts direcotry (scr_path) is empty or invalid")
    if ( is.null(pars$scr_path) || !dir.exists(pars$scr_path) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    errs_mssg <- 
      glue::glue("Input: program executable (exe_path) is empty or invalid")
    if ( is.null(pars$scr_path) || !dir.exists(pars$scr_path) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
  } else if ( pars$run_mode == "Command_Line" ) {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Validate Command Line Parameters::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    pars$prgm_dir <- base::dirname( pars$exe_path )
    pars$prgm_tag <- base::basename( pars$exe_path ) %>% 
      stringr::str_remove("\\.R$")
    
  } else {
    errs_mssg <- 
      glue::glue("Input: Unrecognized run_mode = '{pars$run_mode}'")
    if ( is.null(pars$scr_path) || !dir.exists(pars$scr_path) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Validate Data Directory Parameters::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( prgm_aux_check ) {
    pars$prgm_aux_path <- file.path( pars$aux_path, pars$prgm_tag )
    errs_mssg <- glue::glue( "Input: program data directory (prgm_aux_path) is ",
                             "empty or invalid or does not exist" )
    if ( is.null(pars$prgm_aux_path) || length(pars$prgm_aux_path) == 0 ||
         !dir.exists( pars$prgm_aux_path) ) eflag <- TRUE
    if ( eflag ) stop( glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}") )
    if ( eflag ) return(NULL)
    if ( p6 ) 
      cat(glue::glue("{mssg} Params: prgm_aux_path='{pars$prgm_aux_path}'{RET}"))
  }
  
  ret_cnt <- length(pars)
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  pars
}

# End of file

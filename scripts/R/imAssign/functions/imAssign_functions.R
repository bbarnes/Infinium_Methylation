
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               imAssign::
#                         CG# Database Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Command Line Options Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )

# Tidyverse Core Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
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
#                            File IO Functions:: Tru Imap
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_tru_imap = function( file,
                          spec = NULL,
                          
                          vb=0, vt=3, tc=1, tt=NULL,
                          fun_tag='load_tru_imap') {
  
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
  
  errs_mssg <- glue::glue("File = '{file} does not exist!'")
  if ( is.null(file) || !file.exists(file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if (is.null(spec)) {
    warn_mssg <- glue::glue("Using default spec...")
    if ( vb >= vt+1 ) wflag <- TRUE
    if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    wflag <- FALSE
    
    spec = cols(
      cgn = col_integer(),
      top = col_character(),
    )
    if ( vb >= vt+3 ) print(spec)
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    ret_tib <- safe_read( file = file, use_spec = TRUE, has_head = FALSE, 
                          spec_col = spec, vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    
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
#                            File IO Functions::
#
#  tops_to_cgn()
#  chrom_to_top()
#  chroms_to_tops()
#  genome_to_chrom()
#  valid_bgzipped()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# LEFT OFF HERE!!!!
#
tops_to_cgn = function(file,
                       chr_tag,
                       tru_csv = NULL,
                       cgn_vec = NULL,
                       cgn_max = 0,
                       
                       cols = cols(
                         Fwd_Sequence = col_character(),
                         Exp_Sequence = col_character(),
                         Top_Sequence = col_character(),
                         Cpg_Position = col_integer(),
                         Strand_TB    = col_character()
                       ),
                       
                       out_dir,
                       run_tag,
                       pre_tag = NULL,
                       
                       reload     = 0,
                       reload_min = 2,
                       ret_data   = FALSE,
                       parallel   = FALSE,
                       
                       vb=0, vt=3, tc=1, tt=NULL,
                       fun_tag='tops_to_cgn') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, chr_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  idx_csv <- file.path( out_dir, paste(out_tag, 'idx.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}      chr_tag = '{chr_tag}'.{RET}"))
    cat(glue::glue("{mssg}      tru_csv = '{tru_csv}'.{RET}"))
    cat(glue::glue("{mssg}      cgn_max = '{cgn_max}'.{RET}"))
    cat(glue::glue("{mssg}     func_tag = '{func_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
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
  
  errs_mssg <- glue::glue("File '{file}' does not exist")
  if ( !file.exists(top) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Load Top Sequences::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    top_tib <- readr::read_csv( file = file,
                                # col_names = names(cols),
                                col_types = cols )
    
    ret_dat = push_tops( top_vec = top_tib$Top_Sequence,
                         cgn_vec = cgn_vec,
                         path = out_dir,
                         name = out_tag,
                         write_idx = TRUE,
                         write_seq = FALSE,
                         sep = COM,
                         cgn_max = cgn_max,
                         trim = TRUE,
                         vb = opt$verbose, vt = 2 )
    
    if ( !is.null(tru_csv) ) {
      
    }
    
    # warn_mssg <- glue::glue("WARN_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # wflag <- FALSE
    
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
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

DNAStringSet_to_tops = function( set,
                                 di_nuc = "CG",
                                 chr_str,
                                 
                                 off_tib,
                                 remove_ns = TRUE,

                                 out_dir,
                                 run_tag,
                                 
                                 reload     = 0,
                                 reload_min = 2,
                                 reload_pre = NULL,
                                 
                                 ret_data   = FALSE,
                                 parallel   = FALSE,
                                 
                                 vb=0, vt=3, tc=1, tt=NULL,
                                 fun_tag='DNAStringSet_to_tops')
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
    cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
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
  
  errs_mssg <- glue::glue("File file='{file}' does not exist")
  if ( !file.exists( file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
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
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

chrom_to_top = function( fas,
                         di_nuc = "CG",
                         chr_str,
                         
                         off_tib,
                         remove_ns = TRUE,
                         
                         out_dir,
                         run_tag,
                         pre_tag = NULL,
                         
                         reload     = 0,
                         reload_min = 2,
                         ret_data   = FALSE,
                         parallel   = FALSE,
                         ret_type   = "file",
                         
                         vb=0, vt=3, tc=1, tt=NULL,
                         fun_tag='chrom_to_top') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, chr_str, sep='.' )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid ) {
    if ( ret_type == "file") {
      if ( vb >= vt+1 ) cat(glue::glue("{mssg} All files already exist! {RET}"))
      return( out_csv )
    }
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  }
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     ret_type = '{ret_type}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    if ( purrr::is_character(fas) && file.exists(fas) )
      cat(glue::glue("{mssg}          fas = '{fas}'.{RET}"))
    cat(glue::glue("{mssg}       di_nuc = '{di_nuc}'.{RET}"))
    cat(glue::glue("{mssg}    remove_ns = '{remove_ns}'.{RET}"))
    cat(glue::glue("{mssg}      chr_str = '{chr_str}'.{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      pre_tag = '{pre_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c( out_csv, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    din_dns_pos <- base::abs( off_tib$dns_len[1] ) + 1
    din_ups_pos <- din_dns_pos + off_tib$din_len[1] - 1
    
    flank_len <- off_tib$ups_len[1]
    chr_seq_set <- fas_to_dnaset( x = fas, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    chr_seq_len <- chr_seq_set[[1]]@length - flank_len
    
    # TBD:: This is where we load chromosome partitioned dbSNP VCFs
    #
    #   - all_commnon_YYYYMMDD.vcf.gz
    #        [ https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/ ]
    #        [ rsync -a -P rsync://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153Common.bb ]
    #   - clinvar_YYYYMMDD.vcf.gz 
    #        [ https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/ ]
    #        [ rsync -a -P rsync://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153ClinVar.bb ]
    #
    
    # Skip Chromosomes that are too short
    if ( chr_seq_len >= off_tib$tmp_len[1] ) {
      
      fwd_seqs <- parse_short_seq( chr_seq = chr_seq_set[[1]], 
                                   sub_seq = di_nuc,
                                   
                                   flank_len = flank_len,
                                   remove_ns = remove_ns,
                                   
                                   vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      fwd_seq_vec <- base::as.character(fwd_seqs)
      fwd_pos_vec <- fwd_seqs@ranges@start + flank_len
      
      # Extract test sequences only::
      # tmp_sel_idx_vec <- which( fwd_pos_vec >= 248754851 & fwd_pos_vec <= 248755270 )
      # print(tmp_sel_idx_vec)
      # fwd_seq_vec <- fwd_seq_vec[tmp_sel_idx_vec]
      # fwd_pos_vec <- fwd_pos_vec[tmp_sel_idx_vec]
      #
      
      # TBD:: Should perform expansion here::
      #  - expand top sequences
      #  - write genomic coordinates
      #  - expand all probe designs
      #  - write fasta file and launch bsmap
      #
      
      ret_tib <- fwd2tops_cpp(
        seq_vec_r = fwd_seq_vec,
        pos_vec_r = fwd_pos_vec,
        var_len = 2,
        return_source = TRUE,
        uc = TRUE, 
        vb = vb,
        vt = vt+1 ) %>%
        #
        # NOTE:: Need to sort by position and sequence to ensure reproducible 
        #   order during cg-number assignment!!!
        #
        dplyr::arrange(Cpg_Position, Top_Sequence) %>%
        #
        # NOTE:: Need to filter out non-CG sites after expansion!!! 
        #   e.g. [YG] => [CG]/[TG] TG is not a TG number template sequence!!!
        #
        dplyr::filter( stringr::str_to_upper( 
          stringr::str_sub( Top_Sequence, din_dns_pos,din_ups_pos ) ) == "CG" )
      
      # ret_tib <- 
      #   top_cpp( seqs = fwd_seq_vec, return_source = FALSE, vb = 0 ) %>%
      #   dplyr::mutate(
      #     # Chromosome = chr, # Chromosome is redundant
      #     Coordinate = fwd_pos_vec,
      #     Strand_FR = "F",
      #   )
    }
    
    #
    # TBD:: This case of zero hits needs better handling...
    #
    warn_mssg <- glue::glue("Found zero {di_nuc} sites in chromosome {chr_str}")
    if ( is.null(ret_tib) || base::nrow(ret_tib) == 0 ) wflag <- TRUE
    if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    wflag <- FALSE
    
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
  
  if ( ret_type == "file" ) return( out_csv )
  
  ret_tib
}

chroms_to_tops = function( x,
                           
                           di_nuc  = "CG",
                           chr_key = "Chr_Key",
                           fas_key = "Chr_Bgz",
                           top_key = "Top_Fas",
                           
                           off_tib,
                           remove_ns = TRUE,
                           
                           build,
                           out_dir,
                           run_tag,
                           pre_tag = NULL,
                           
                           max = 0,
                           reload     = TRUE,
                           reload_min = 1,
                           ret_data   = FALSE,
                           parallel   = FALSE,
                           ret_type   = "file",
                           
                           vb=0, vt=3, tc=1, tt=NULL,
                           fun_tag='chroms_to_tops') {
  
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
  out_csv <- file.path( out_dir, "chromosome-order-list.csv.gz")
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          max = '{max}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     ret_type = '{ret_type}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       di_nuc = '{di_nuc}'.{RET}"))
    cat(glue::glue("{mssg}      chr_key = '{chr_key}'.{RET}"))
    cat(glue::glue("{mssg}      fas_key = '{fas_key}'.{RET}"))
    cat(glue::glue("{mssg}      top_key = '{top_key}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}    remove_ns = '{remove_ns}'.{RET}"))
    cat(glue::glue("{mssg}        build = '{build}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      pre_tag = '{pre_tag}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    chr_sym <- rlang::sym( chr_key )
    top_sym <- rlang::sym( top_key )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       LOAD:: Chromosome Order List
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    chr_tib <- NULL
    chr_tib <- safe_read( file = x, spec_prefernce = "rds", type = "csv",
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    ret_key <- glue::glue("chr_files_tib")
    ret_cnt <- print_tib( chr_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    if ( !tibble::is_tibble( chr_tib ) ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("Input: x is neither tibble or file")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #           Split Chromosomes and Set Order For Sorting Later::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( vb >= vt+1 )
      cat(glue::glue("{mssg} Processing Reference: '{build}'...{RET}"))
    
    chr_names <- chr_tib %>% dplyr::pull( !!chr_sym )
    chr_split <- split( chr_tib, f = chr_names )
    
    ret_key <- glue::glue("chr_split")
    ret_cnt <- print_tib( chr_split, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # Record the order of the chromosomes
    chr_count <- chr_names %>% length()
    unq_count <- chr_names %>% unique() %>% length()
    if ( max != 0 ) chr_names <- head( chr_names, n = max )
    tar_count <- chr_names %>% unique() %>% length()
    
    if ( vb >= vt+1 ) 
      cat(glue::glue("{mssg}  Total Ref Chromosomes = {chr_count}.{RET}",
                     "{mssg} Unique Ref Chromosomes = {unq_count}.{RET}",
                     "{mssg} Target Ref Chromosomes = {tar_count}.{RET2}") )
    if ( vb >= vt+4 )
      cat(glue::glue("{mssg}{TAB}Will target Chromosomes: '{chr_names}'{RET}") )
    
    all_csvs <- list()
    
    if ( parallel ) {
      if ( vb >= vt+1 ) cat(glue::glue(
        "{mssg} Processing Reference: '{build}' in parallele mode.{RET}") )
      
      all_csvs <- 
        foreach( chr_str = chr_names, .inorder = TRUE, 
                 .final = function(x) setNames(x,chr_names)) %dopar% {
                   
                   chr_dir <- file.path( out_dir, chr_str )
                   chr_fas <- chr_split[[chr_str]] %>% head(n=1) %>% 
                     dplyr::pull( fas_key )
                   
                   if ( vb >= vt+1 ) 
                     cat(glue::glue("{mssg}{TAB}Current: '{build}' chromosome ",
                                    "'{chr_str}'{RET}") )
                   
                   cur_csv <- chrom_to_top( fas = chr_fas,
                                            di_nuc  = di_nuc,
                                            chr_str = chr_str, 
                                            off_tib = off_tib,
                                            remove_ns = remove_ns,
                                            
                                            out_dir = chr_dir,
                                            run_tag = run_tag,
                                            pre_tag = pre_tag,
                                            
                                            reload     = reload, 
                                            reload_min = reload_min,
                                            ret_type   = ret_type,
                                            
                                            vb=vb,vt=vt+4,tc=tc+1,tt=tt )
                   if ( vb >= vt+1 )
                     cat(glue::glue("{mssg} Done. chr: {chr_str} = ",
                                    "'{cur_csv}'{RET2}") )
                   
                   cur_csv
                 }
      
    } else {
      if ( vb >= vt+1 ) cat(glue::glue(
        "{mssg} Processing Reference: '{build}' in linear mode.{RET}") )
      
      for (chr_str in chr_names ) {
        chr_dir <- file.path( out_dir, chr_str )
        chr_fas <- chr_split[[chr_str]] %>% head(n=1) %>% dplyr::pull(fas_key)
        
        if ( vb >= vt+1 ) cat(glue::glue(
          "{mssg}{TAB}Current: '{build}' chromosome '{chr_str}'{RET}") )
        
        
        cur_csv <- chrom_to_top( fas = chr_fas,
                                 di_nuc  = di_nuc,
                                 chr_str = chr_str, 
                                 off_tib = off_tib,
                                 remove_ns = remove_ns,
                                 
                                 out_dir = chr_dir,
                                 run_tag = run_tag,
                                 pre_tag = pre_tag,
                                 
                                 reload     = reload, 
                                 reload_min = reload_min,
                                 ret_type   = ret_type,
                                 
                                 vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        all_csvs[[chr_str]] <- cur_csv
        
        if ( vb >= vt+1 )
          cat(glue::glue("{mssg} Done. chr: {chr_str} = '{cur_csv}'{RET2}") )
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #           Force Original Order and Write Chromosome Order File::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- all_csvs %>% unlist() %>% 
      as.data.frame() %>% 
      as_tibble(rownames = chr_key ) %>%
      purrr::set_names(c(chr_key,top_key)) %>%
      dplyr::arrange(factor(!!chr_sym, levels = chr_names)) %>%
      dplyr::right_join(chr_tib, by = chr_key ) %>%
      dplyr::filter( !is.na(!!top_sym) )
    
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
  
  if ( ret_type == "file") return( out_csv )
  
  ret_tib
}

genome_to_chrom = function(file,
                           
                           skip_vec = NULL,
                           
                           chr_key  = "Chr_Key",
                           fas_key  = "Chr_Fas",
                           bgz_key  = "Chr_Bgz",
                           gzi_key  = "Chr_Gzi",
                           fai_key  = "Chr_Fai",
                           head_key = "Fasta_Header",
                           
                           build,
                           source,
                           valid_source = c("UCSC"),
                           validate_only = FALSE,
                           ucsc_ftp = 'ftp://hgdownload.cse.ucsc.edu/goldenPath',
                           
                           out_dir,
                           run_tag,
                           pre_tag = NULL,
                           
                           max = 0,
                           reload     = 0,
                           reload_min = 1,
                           ret_data   = FALSE,
                           parallel   = FALSE,
                           ret_type   = "file",
                           
                           vb=0, vt=3, tc=1, tt=NULL,
                           fun_tag='genome_to_chrom') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  eflag <- FALSE
  wflag <- FALSE
  
  out_tag <- file %>% stringr::str_remove("\\.gz$")
  chr_csv <- paste( out_tag, "chrs-info.csv.gz", sep = '.' )
  out_csv <- file.path( out_dir, "chromosome-order-list.csv.gz")
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  # Fastest way to check chromosomes from original fasta, large number of
  #   individual time stamps and account for all chromosomes, plus corrrect
  #   file output order...
  if ( reload >= reload_min && file.exists(out_csv) && file.exists(chr_csv) ) {
    chr_sym <- rlang::sym(chr_key)
    
    ret_tib <- NULL
    ret_tib <- safe_read( file = out_csv, use_spec = TRUE, type = "csv",
                          spec_prefernce = "rds", vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    ret_cnt <- ret_tib %>% dplyr::distinct(!!chr_sym) %>% base::nrow()
    
    # Load chr_tib from Referece_File_Tag.chrs-info.csv.gz
    chr_tib <- safe_read( file = chr_csv, type = "csv", use_spec = TRUE, 
                          write_spec = TRUE, spec_prefernce = "rds",
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    fas_tib <- chr_tib %>%
      dplyr::left_join(
        sort_files_by_date(path = out_dir, 
                           pattern = "\\.fa\\.[fagzi]+$", by = "mtime") %>% 
          tibble::as_tibble() %>% 
          dplyr::select(file_names) %>% 
          dplyr::mutate(file_names = base::basename(file_names)) %>% 
          tidyr::separate(file_names, into=c(chr_key, "suffix"), sep="\\.fa\\.") %>% 
          dplyr::group_by( !!chr_sym ) %>% 
          dplyr::mutate(Rank = dplyr::row_number() ) %>% 
          dplyr::ungroup() %>% 
          dplyr::arrange(!!chr_sym, Rank),
        by=( chr_key )
      )
    
    fas_cnt <- as.integer(fas_tib %>% base::nrow() / 3)
    fas_val <- dplyr::filter( fas_tib, is.na(suffix) ) %>% base::nrow() == 0
    ret_val <- ret_tib %>% 
      dplyr::filter(!!chr_sym %in% dplyr::pull( fas_tib, !!chr_sym ) ) %>%
      dplyr::distinct(!!chr_sym) %>%
      base::nrow() == ret_cnt
    cnt_val <- ret_cnt == fas_cnt
    
    idx1_vec <- seq(1, fas_cnt, by = 3)
    idx2_vec <- seq(2, fas_cnt, by = 3)
    idx3_vec <- seq(3, fas_cnt, by = 3)
    
    bgz_val <- fas_tib[ idx1_vec, ]$suffix %>% unique() == "gz"
    gzi_val <- fas_tib[ idx2_vec, ]$suffix %>% unique() == "gzi"
    fai_val <- fas_tib[ idx3_vec, ]$suffix %>% unique() == "fai"
    
    if ( vb >= vt+8 ) {
      cat(glue::glue("{mssg} {TAB}ret_cnt = '{ret_cnt}'{RET}"))
      cat(glue::glue("{mssg} {TAB}fas_cnt = '{fas_cnt}'{RET}"))
      cat(glue::glue("{mssg} {TAB}ret_val = '{ret_val}'{RET}"))
      cat(glue::glue("{mssg} {TAB}cnt_val = '{cnt_val}'{RET}"))
      cat(glue::glue("{mssg} {TAB}fas_val = '{fas_val}'{RET}"))
      cat(glue::glue("{mssg} {TAB}bgz_val = '{bgz_val}'{RET}"))
      cat(glue::glue("{mssg} {TAB}gzi_val = '{gzi_val}'{RET}"))
      cat(glue::glue("{mssg} {TAB}fai_val = '{fai_val}'{RET}"))
    }
    
    if ( fas_val && ret_val && cnt_val && bgz_val && gzi_val && fai_val ) {
      if ( vb >= vt+1 ) cat(glue::glue("{mssg} All files already exist!{RET}",
                                       "{mssg}{TAB}     reload='{reload}'{RET}",
                                       "{mssg}{TAB} reload_min='{reload_min}'{RET2}"))
      if ( ret_type == "file") return( out_csv )
      return( ret_tib )
    }
  }
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          max = '{max}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     ret_type = '{ret_type}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   chr_key = '{chr_key}'.{RET}"))
    cat(glue::glue("{mssg}   fas_key = '{fas_key}'.{RET}"))
    cat(glue::glue("{mssg}  head_key = '{head_key}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}      file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}     build = '{build}'.{RET}"))
    cat(glue::glue("{mssg}   out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  
  errs_mssg <- glue::glue("Reference genome fasta {file} does not exist")
  if ( !file.exists(file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  safe_mkdir( out_dir )
  
  # TBD:: make the bio_tib first!!!
  is_valid_bgzip <- valid_bgzipped( file = file,
                                    
                                    chr_key = chr_key,
                                    
                                    build  = build,
                                    source = source,
                                    valid_source = valid_source,
                                    ucsc_ftp = ucsc_ftp,
                                    
                                    reload = reload,
                                    reload_min = reload_min,
                                    
                                    vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  if ( validate_only ) {
    if ( vb >= vt+1 ) {
      cat(glue::glue("{mssg} Validation done. Exiting...{RET2}"))
    }
    return( out_csv )
  }
  
  ftime <- base::system.time({
    
    chr_sym  <- rlang::sym( chr_key )
    fas_sym  <- rlang::sym( fas_key )
    bgz_sym  <- rlang::sym( bgz_key )
    gzi_sym  <- rlang::sym( gzi_key )
    fai_sym  <- rlang::sym( fai_key )
    head_sym <- rlang::sym( head_key )
    
    if ( vb >= vt+1 ) {
      cat(glue::glue("{mssg} Processing Reference: '{build}'...{RET}"))
      cat(glue::glue("{mssg}   Fasta = '{file}'...{RET}"))
    }
    
    # Load Genome File:: Note this will maintain order when we re-load 
    #  everything after parallelization...
    #
    chr_seqs <- Biostrings::readDNAStringSet(filepath = file, format = "fasta" )
    if ( vb >= vt+1 )
      cat(glue::glue("{mssg}   Done loading Fasta.{RET2}"))
    
    # if ( vb >= vt+3 ) print(chr_seqs)
    
    chr_names <- names( chr_seqs )
    chr_count <- chr_names %>% length()
    unq_count <- chr_names %>% unique() %>% length()
    if ( max != 0 ) chr_names <- head( chr_names, n = max )
    max_count <- chr_names %>% length()
    if ( vb >= vt+1 ) cat(glue::glue(
      "{mssg} Total Reference Chromosomes: {chr_count}; ",
      "unique = {unq_count}, max = {max_count}.{RET}") )
    
    if ( !is.null(skip_vec) && length(skip_vec) > 0 ) {
      chr_names <- chr_names[ which( !chr_names %in% skip_vec ) ]
      chr_count <- chr_names %>% length()
      unq_count <- chr_names %>% unique() %>% length()
      
      if ( vb >= vt+1 ) cat(glue::glue(
        "{mssg} Total Reference Chromosomes: {chr_count}; ",
        "unique = {unq_count}, max = {max_count}.{RET}") )
      # if ( vb >= vt+1 ) cat(glue::glue(
      #   "{mssg} '{build}' Chromosome: {chr_str} exists in skip list...{RET}",
      #   "{mssg}{TAB}Fasta = '{chr_bgz}'{RET}") )
      # next
      
      # return( chr_names )
    }
    
    # Split Chromosomes from Genome Build::
    split_seqs <- base::split( x = chr_seqs, f = names(chr_seqs) )
    ret_key <- glue::glue("split_seqs-list")
    ret_cnt <- print_tib( split_seqs, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    if ( vb >= vt+4 )
      cat(glue::glue("{mssg}{TAB}Will target Chromosomes: '{chr_names}'{RET}") )
    
    chr_fastas <- list()
    if ( parallel ) {
      if ( vb >= vt+1 ) cat(glue::glue(
        "{mssg} Processing Reference: '{build}' in parallele mode.{RET}") )
      
      chr_fastas <- foreach( 
        chr_tag = chr_names, .inorder = TRUE, 
        .final = function(x) setNames(x, chr_names ) ) %dopar% {
          
          chr_str <- chr_tag %>% stringr::str_remove("\\s+.*$")
          chr_fas <- file.path( out_dir, paste(chr_str,"fa", sep = '.') )
          chr_fai <- paste(chr_fas,"fai", sep=".")
          chr_gzi <- paste(chr_fas,"gzi", sep=".")
          chr_bgz <- paste(chr_fas,"gz",  sep=".")
          
          is_valid <- valid_time_stamp( c(pre_tag, file, 
                                          chr_bgz, chr_gzi, chr_fai ), 
                                        out_dir = out_dir,
                                        vb=vb,vt=vt+4,tc=tc+1,tt=tt )
          
          if ( reload >= reload_min && is_valid ) {
            if ( vb >= vt+1 ) cat(glue::glue(
              "{mssg} '{build}' Chromosome: {chr_str} already exists...{RET}",
              "{mssg}{TAB}Fasta = '{chr_bgz}'{RET}") )
            
          } else {
            if ( vb >= vt+1 ) cat(glue::glue(
              "{mssg} Building '{build}' Chromosome: {chr_str}...{RET}",
              "{mssg}{TAB}Fasta = '{chr_fas}'{RET}",
              "{mssg}{TAB}Fasta = '{chr_bgz}'{RET}",
              "{mssg}{TAB}Fasta = '{chr_gzi}'{RET}",
              "{mssg}{TAB}Fasta = '{chr_fai}'{RET}" ) )
            
            Biostrings::writeXStringSet(x = split_seqs[[ chr_tag ]], 
                                        filepath = chr_fas, 
                                        append = FALSE, 
                                        compress = FALSE, 
                                        format = "fasta" )
            Biostrings::fasta.index( filepath = chr_fas, seqtype = "DNA" )
            
            cmd <- glue::glue("bgzip -f {chr_fas}")
            if ( vb >= vt+1 ) 
              cat(glue::glue("{mssg} Bgzip Compressing:: '{cmd}'...{RET}"))
            system( cmd )
            
            cmd <- glue::glue("samtools faidx --fai-idx {chr_fai} ",
                              "--gzi-idx {chr_gzi} {chr_bgz}")
            if ( vb >= vt+1 ) 
              cat(glue::glue("{mssg} Creating Faidx:: '{cmd}'...{RET}"))
            system( cmd )
            
          }
          
          chr_bgz
        }
      
    } else {
      if ( vb >= vt+1 ) cat(glue::glue(
        "{mssg} Processing Reference: '{build}' in linear mode.{RET}") )
      
      for (chr_tag in chr_names ) {
        chr_str <- chr_tag %>% stringr::str_remove("\\s+.*$")
        chr_fas <- file.path( out_dir, paste(chr_str,"fa", sep = '.') )
        chr_fai <- paste(chr_fas,"fai", sep=".")
        chr_gzi <- paste(chr_fas,"gzi", sep=".")
        chr_bgz <- paste(chr_fas,"gz",  sep=".")
        
        is_valid <- valid_time_stamp( c(pre_tag, file, 
                                        chr_bgz, chr_gzi, chr_fai ), 
                                      out_dir = out_dir,
                                      vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        
        if ( reload >= reload_min && is_valid ) {
          if ( vb >= vt+1 ) cat(glue::glue(
            "{mssg} '{build}' Chromosome: {chr_str} already exists...{RET}",
            "{mssg}{TAB}Fasta = '{chr_bgz}'{RET}") )
          
        } else {
          if ( vb >= vt+1 ) cat(glue::glue(
            "{mssg} Building '{build}' Chromosome: {chr_str}...{RET}",
            "{mssg}{TAB}Fasta = '{chr_fas}'{RET}",
            "{mssg}{TAB}Fasta = '{chr_bgz}'{RET}",
            "{mssg}{TAB}Fasta = '{chr_gzi}'{RET}",
            "{mssg}{TAB}Fasta = '{chr_fai}'{RET}" ) )
          
          # unlink( chr_fas, chr_bgz, chr_fai, chr_gzi )
          
          Biostrings::writeXStringSet(x = split_seqs[[ chr_tag ]], 
                                      filepath = chr_fas, 
                                      append = FALSE, 
                                      compress = FALSE, 
                                      format = "fasta" )
          Biostrings::fasta.index( filepath = chr_fas, seqtype = "DNA" )
          
          cmd <- glue::glue("bgzip -f {chr_fas}")
          if ( vb >= vt+1 ) 
            cat(glue::glue("{mssg} Bgzip Compressing:: '{cmd}'...{RET}"))
          system( cmd )
          
          cmd <- glue::glue("samtools faidx --fai-idx {chr_fai} ",
                            "--gzi-idx {chr_gzi} {chr_bgz}")
          if ( vb >= vt+1 ) 
            cat(glue::glue("{mssg} Creating Faidx:: '{cmd}'...{RET}"))
          system( cmd )
          
        }
        chr_fastas[[ chr_tag ]] <- chr_bgz
      }
      
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #           Force Original Order and Write Chromosome Order File::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Force the correct order::
    #
    ret_tib <- chr_fastas %>% unlist() %>% 
      as.data.frame() %>% as_tibble(rownames = head_key ) %>% 
      purrr::set_names( c(head_key,bgz_key) ) %>%
      dplyr::arrange( factor(!!head_sym, levels = chr_names) ) %>%
      dplyr::mutate( !!chr_sym := !!head_sym %>% 
                       stringr::str_remove(" .*$"),
                     !!fas_sym := !!bgz_sym %>% 
                       stringr::str_remove("\\.gz$"),
                     !!gzi_sym := paste0(!!fas_sym,".gzi"),
                     !!fai_sym := paste0(!!fas_sym,".fai") ) %>%
      dplyr::select( !!chr_sym, 
                     !!fas_sym, !!bgz_sym, !!gzi_sym, !!fai_sym,
                     !!head_sym )
    
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
  
  if ( ret_type == "file") return( out_csv )
  
  ret_tib
}

valid_bgzipped = function(file,
                          
                          chr_key = "Chr_Key",
                          
                          build  = NULL,
                          source = NULL,
                          valid_source = c("UCSC"),
                          ucsc_ftp = 'ftp://hgdownload.cse.ucsc.edu/goldenPath',
                          
                          reload     = 0,
                          reload_min = 2,
                          
                          vb=0, vt=3, tc=1, tt=NULL,
                          fun_tag='valid_bgzipped') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  eflag <- FALSE
  wflag <- FALSE
  
  out_dir <- base::dirname(file)
  out_tag <- file %>% stringr::str_remove("\\.gz$")
  zip_txt <- paste( out_tag, "file-info.txt", sep = '.' )
  chr_csv <- paste( out_tag, "chrs-info.csv.gz", sep = '.' )
  out_gzi <- paste( out_tag, "gzi", sep=".")
  out_fai <- paste( out_tag, "fai", sep=".")
  end_txt <- paste( out_tag, 'gz.done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c( chr_csv, out_gzi, out_fai, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid ) {
    if ( vb >= vt+1 ) cat(glue::glue("{mssg} All files already exist! {RET}"))
    return( TRUE )
  }
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}  chr_key = '{chr_key}'.{RET}"))
    cat(glue::glue("{mssg}    build = '{build}'.{RET}"))
    cat(glue::glue("{mssg}   source = '{source}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  chr_sym <- rlang::sym( chr_key )
  ret_cnt <- 0
  
  # NOTE:: If file doesn't exist we'll have to get it from UCSC/NCBI.
  #  TBD:: Right now I only pull from UCSC. Need to add NCBI!!!
  #
  if ( !file.exists(file) ) {
    errs_mssg <- glue::glue("File OR (source AND build) must exist")
    if ( is.null(source) || is.null(build) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    errs_mssg <- 
      glue::glue("Unsupported source = '{source}'; valid = '{valid_source}'")
    if ( !source %in% valid_source ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    warn_mssg <- glue::glue("File '{file}' does not exist. Will download ",
                            "source from '{source}'")
    if ( !file.exists(file) == 0 ) wflag <- TRUE
    if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    wflag <- FALSE
    
    if ( source == "UCSC" ) {
      
      out_fas <- paste0(build,'.fa.gz')
      ftp_str <- file.path( ucsc_ftp, build, "bigZips", out_fas )
      get_cmd <- glue::glue("wget --timestamping {ftp_str} -O {file}")
      
      if ( vb >= vt+1 ) cat(glue::glue("{mssg} FTP: '{get_cmd}'...{RET}"))
      get_ret <- system(get_cmd)
      if ( vb >= vt+1 ) cat(glue::glue("{mssg} RET: '{get_ret}'.{RET2}"))
    }
  }
  
  # Sanity Check:: Validate File Really Exists Now!!!
  #
  errs_mssg <- glue::glue("Somehow file still doesn't exist even after source")
  if ( !file.exists(file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  # Write out chromosome list for look up later. This will probably require
  #  unzipping and re-zipping, but its done once and then everything is 
  #  reliable moving forward (chromsome count, order, etc.)
  if ( !file.exists( chr_csv ) || 
       file.info(chr_csv, extra_cols = FALSE)[["mtime"]] >
       file.info(file, extra_cols = FALSE)[["mtime"]] ) {
    
    chr_tib <- 
      Biostrings::fasta.index( filepath = file, seqtype = "DNA" ) %>% 
      tibble::as_tibble() %>% 
      dplyr::select( desc ) %>% 
      dplyr::mutate( !!chr_sym := desc %>% stringr::str_remove(" .*$") ) %>%
      dplyr::select( !!chr_sym, dplyr::everything() )
    
    safe_write( x = chr_tib, file = chr_csv, type = "csv", write_spec = TRUE,
                done = TRUE, append = FALSE, vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    tap_cmd <- glue::glue("touch {file}")
    tap_ret <- system(tap_cmd)
    
    if ( file.exists(out_fai) ) unlink( out_fai )
    if ( file.exists(out_gzi) ) unlink( out_gzi )
    if ( file.exists(end_txt) ) unlink( end_txt )
    
    # if ( file.exists(out_gzi) ) {
    #   tap_cmd <- glue::glue("touch {out_gzi}")
    #   tap_ret <- system(tap_cmd)
    # }
    # if ( file.exists(out_fai) ) {
    #   tap_cmd <- glue::glue("touch {out_fai}")
    #   tap_ret <- system(tap_cmd)
    # }
  }
  
  # Check if file was bgzipped!!!
  #
  cmd <- glue::glue("file {file} > {zip_txt}")
  system(cmd)
  lines = readr::read_lines(zip_txt)
  unlink( zip_txt )
  
  is_gzipped <- FALSE
  is_bgziped <- FALSE
  is_gzipped <- stringr::str_detect( lines, "compressed" )
  if ( is_gzipped )
    is_bgziped <- stringr::str_detect( lines, "extra field" )
  
  # is_faidx <- FALSE
  # if ( file.exists(out_fai) ) is_faidx <- TRUE
  # is_gzidx <- FALSE
  # if ( file.exists(out_gzi) ) is_gzidx <- TRUE
  
  is_valid <- valid_time_stamp( c( chr_csv, out_gzi, out_fai, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( is_bgziped && is_valid ) {
    if ( vb >= vt+1 ) cat(glue::glue("{mssg} All files up to date! {RET}"))
    return( TRUE )
  }
  
  warn_mssg <- glue::glue("File {file} is not bgzip compressed")
  if ( !is_bgziped ) wflag <- TRUE
  if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
  wflag <- FALSE
  
  ftime <- base::system.time({
    
    if ( file.exists(out_fai) ) unlink( out_fai )
    if ( file.exists(out_gzi) ) unlink( out_gzi )
    if ( file.exists(end_txt) ) unlink( end_txt )
    
    if ( !is_bgziped ) {
      if ( is_gzipped ) {
        cmd <- glue::glue("gzip -d -f {file} ")
        if ( vb >= vt+1 ) cat(glue::glue("{mssg} Uncompressing:: '{cmd}'...{RET}"))
        system( cmd )
        file <- file %>% stringr::str_remove("\\.gz$")
      }
      cmd <- glue::glue("bgzip -f {file}")
      if ( vb >= vt+1 ) cat(glue::glue("{mssg} Bgzip Compressing:: '{cmd}'...{RET}"))
      system( cmd )
      file = paste(file,"gz", sep=".")
    }
    
    cmd <- glue::glue("samtools faidx --fai-idx {out_fai} --gzi-idx {out_gzi} {file}")
    if ( vb >= vt+1 ) cat(glue::glue("{mssg} Creating Faidx:: '{cmd}'...{RET}"))
    system( cmd )
    
    tap_cmd <- glue::glue("touch {end_txt}")
    tap_ret <- system(tap_cmd)
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  TRUE
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Local Run Time Defaults::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imAssign_options_old = function( pars, args,
                                    
                                    vb=0, vt=4, tc=1, tt=NULL,
                                    fun_tag='imAssign_options_old') {
  
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
  
  opts$run_name     <- NULL
  opts$out_path     <- NULL
  
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
    
    # Read Options and Test Cases from pre-defined function
    # loc_dat <- lighthoof_local_defaults( opts = opt,
    #                                      pars = par,
    #                                      vb = opts$verbose )
    # opt <- loc_dat$opt
    # par <- loc_dat$par
    
    pars$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R")
    
    opts$out_path  <- file.path( pars$top_path, 'scratch' )
    opts$imp_path  <- file.path( pars$top_path, 'data/improbe' )
    opts$ann_path  <- file.path( pars$top_path, 'data/annotation' )
    opts$man_path  <- file.path( pars$top_path, 'data/manifests' )
    opts$idat_path <- file.path( pars$top_path, 'data/idats' )
    
    opts$run_name = "test"
    
    opts$ref_source <- paste( "NCBI", sep = ',' )
    # opts$ref_source <- paste( "UCSC", sep = ',' )
    
    opts$run_name <- paste(opts$run_name,opts$ref_source, sep='-')
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
        c("--out_path"), type="character", default=opts$out_path,
        help=glue::glue(
          "Build output directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
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
    
  } else {
    stop( glue::glue("{RET}[{pars$prgm_tag}]: ERROR: Unrecognized run_mode = ",
                     "'{pars$run_mode}'!{RET2}") )
  }
  
  ret_cnt <- length(opts)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

imAssign_options = function( pars, args,
                                
                                vb=0, vt=4, tc=1, tt=NULL,
                                fun_tag='imAssign_options') {
  
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
  
  opts$run_name     <- NULL
  opts$out_path     <- NULL
  
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
    
    # Read Options and Test Cases from pre-defined function
    # loc_dat <- lighthoof_local_defaults( opts = opt,
    #                                      pars = par,
    #                                      vb = opts$verbose )
    # opt <- loc_dat$opt
    # par <- loc_dat$par
    
    # pars$top_path <- pars$src_path %>% 
    #   stringr::str_remove("/tools/Workhorse-Unstained/scripts/R")
    pars$top_path <- pars$src_path %>% 
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
    
    # opts$ref_source <- paste( "NCBI", sep = ',' )
    # opts$ref_source <- paste( "UCSC", sep = ',' )
    opts$ref_source <- paste( "Both", sep = ',' )

    # Trying to make this work for imAssign::    
    opts$ref_source <- paste( "NCBI", sep = ',' )
    
    
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
                              "galGal5.fa.gz",
                              "criGriChoV1.fa.gz",
                              "pvirPacbioDovetail.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 "Gallus_gallus",
                                 "ChineseHamster",
                                 "MarbledCrayfish",
                                 sep = ',' )
      
      opts$ref_build <- paste( "hg38",
                               "hg19",
                               "hg18",
                               "mm10",
                               "galgal5",
                               "criGrichov1",
                               "GRCc00",
                               sep = ',' )
      
    } else if ( opts$ref_source == "NCBI" ) {
      
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
      
    } else if ( opts$ref_source == "Both" ) {
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "GRCh38.genome.fa.gz",
                              "GRCh37.genome.fa.gz",
                              "GRCh36.genome.fa.gz",
                              "GRCm38.genome.fa.gz",
                              
                              "galGal5.fa.gz",
                              "criGriChoV1.fa.gz",
                              "pvirPacbioDovetail.fa.gz",
                              
                              "canFam3.genome.fa.gz",
                              "canFam6.genome.fa.gz",

                              # "hg38.fa.gz",
                              # "hg19.fa.gz",
                              # "hg18.fa.gz",
                              # "mm10.fa.gz",
                              
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 
                                 "Gallus_gallus",
                                 "ChineseHamster",
                                 "MarbledCrayfish",
                                 
                                 "Canis_lupus_familiaris",
                                 "Canis_lupus_familiaris",
                                 
                                 # "Homo_sapiens",
                                 # "Homo_sapiens",
                                 # "Homo_sapiens",
                                 # "Mus_musculus",
                                 
                                 sep = ',' )
      
      opts$ref_build <- paste( "GRCh38",
                               "GRCh37",
                               "GRCh36",
                               "GRCm10",
                               
                               "galgal5",
                               "criGrichov1",
                               "GRCc00",
                               
                               "canFam3",
                               "canFam6",
                               
                               # "hg38",
                               # "hg19",
                               # "hg18",
                               # "mm10",
                               
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
  
  ret_cnt <- length(opts)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

# End of file

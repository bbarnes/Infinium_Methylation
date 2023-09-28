
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               imGenomes::
#                          BioStrings Functions::
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
#                            File IO Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

extract_genomic_templates = function( tib,
                                      path,
                                      
                                      ids_key = "Genomic_ID",
                                      chr_key = "Chromosome",
                                      pos_key = "Coordinate",
                                      fwd_key = "Forward_Sequence",
                                      din_key = "Din_Seq",
                                      
                                      buf_len = 60,
                                      var_len = 2,
                                      
                                      out_dir,
                                      run_tag,
                                      pre_tag = NULL,
                                      
                                      reload     = 0,
                                      reload_min = 2,
                                      ret_data   = FALSE,
                                      parallel   = FALSE,
                                      
                                      vb=0, vt=3, tc=1, tt=NULL,
                                      fun_tag='extract_genomic_templates' )
{
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
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}         path = '{path}'.{RET}"))
    cat(glue::glue("{mssg}      ids_key = '{ids_key}'.{RET}"))
    cat(glue::glue("{mssg}      chr_key = '{chr_key}'.{RET}"))
    cat(glue::glue("{mssg}      pos_key = '{pos_key}'.{RET}"))
    cat(glue::glue("{mssg}      fwd_key = '{fwd_key}'.{RET}"))
    cat(glue::glue("{mssg}      din_key = '{din_key}'.{RET}"))
    cat(glue::glue("{mssg}      buf_len = '{buf_len}'.{RET}"))
    cat(glue::glue("{mssg}      var_len = '{var_len}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      pre_tag = '{pre_tag}'.{RET}"))
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
  
  errs_mssg <- glue::glue("Directory path path='{path}' does not exist")
  if ( !dir.exists( path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ids_sym = rlang::sym( ids_key )
  chr_sym = rlang::sym( chr_key )
  pos_sym = rlang::sym( pos_key )
  fwd_sym = rlang::sym( fwd_key )
  din_sym = rlang::sym( din_key )
  
  seq_len = var_len + ( buf_len * 2 )
  DN0_idx = buf_len + 1
  DN1_idx = buf_len + 2
  
  chrom_names = tib %>% 
    dplyr::distinct( !!chr_sym ) %>% 
    dplyr::pull( !!chr_sym )

  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    for ( chr in chrom_names ) {
      chr_fas <- file.path( path, paste0(chr,".fa.gz") )
      
      if ( !file.exists(chr_fas) ) {
        if ( vb >= vt+1 )
          cat( glue::glue("{mssg} Failed to find chromosome fasta = '{chr_fas}'{RET}") )
        next;
      }
      chr_tib <- tib %>% 
        dplyr::filter( !!chr_sym == chr ) %>%
        dplyr::distinct( !!pos_sym, .keep_all=TRUE ) %>%
        dplyr::mutate( !!ids_sym := paste( !!chr_sym, !!pos_sym, sep="_") ) %>%
        dplyr::select( dplyr::all_of( c(ids_key,chr_key,pos_key) ), 
                       dplyr::everything() )
        # dplyr::group_by( !!ids_sym ) %>% 
        # dplyr::mutate( Rep_Idx = dplyr::row_number() ) %>%
        # dplyr::ungroup()
      
      chr_cnt <- chr_tib %>% base::nrow()
      
      if ( chr_cnt == 0 ) {
        if ( vb >= vt+1 )
          cat( glue::glue("{mssg} Skipping... chr_cnt = '{chr_cnt}'{RET}") )
        next;
      }
      ret_key <- glue::glue("{chr}-tib")
      ret_cnt <- print_tib( chr_tib, fun_tag = fun_tag, name = ret_key, 
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )

      cat( glue::glue("{mssg} Loading chromosome fasta '{chr_fas}'{RET}") )
      chr_seq <- Biostrings::readDNAStringSet(filepath = chr_fas, format = "fasta")
      
      pos_vec <- chr_tib %>% dplyr::pull( !!pos_sym )
      ids_vec <- chr_tib %>% dplyr::pull( !!ids_sym )

      cpg_iRange <- IRanges::IRanges( 
        start = pos_vec - buf_len, 
        width = seq_len,
        names = ids_vec )
      
      fwd_seqs <- Biostrings::extractAt( chr_seq[[1]], at = cpg_iRange )
      print( fwd_seqs )
      
      fwd_seq_vec <- base::as.character( fwd_seqs )
      fwd_pos_vec <- fwd_seqs@ranges@start + buf_len
      
      dif_pos_tib <- tibble::tibble(
        pos1 = pos_vec,
        pos2 = fwd_pos_vec,
      ) %>% dplyr::mutate(
        pos_delta = pos1 - pos2,
      )
      dif_pos_cnt <- dif_pos_tib %>% 
        dplyr::filter( pos_delta != 0 ) %>% 
        base::nrow()
      
      if ( vb >= vt+1 )
        cat( glue::glue("{mssg} dif_pos_cnt='{dif_pos_cnt}'.{RET}") )

      cur_tib <- fwd2tops_cpp(
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
        dplyr::arrange( 
          Cpg_Position, Top_Sequence ) %>%
        dplyr::mutate( 
          !!chr_sym := chr,
          Din_Fwd = Fwd_Sequence %>%
            stringr::str_sub( DN0_idx, DN1_idx ),
          Din_Top = Top_Sequence %>%
            stringr::str_sub( DN0_idx, DN1_idx ),
          SNP_Range_Beg = Cpg_Position - buf_len,
          SNP_Range_End = Cpg_Position + seq_len,
          Probe_Type = dplyr::case_when(
            Din_Fwd == "CG" & Din_Fwd == Din_Top ~ "cg",
            Din_Fwd %>% stringr::str_starts("C") | Din_Top %>% stringr::str_starts("C") ~ "ch",
            TRUE ~ "rs"
          )
        )
      
      cur_bed <- file.path( out_dir, paste0(chr,"-snp-range-buf-",buf_len,".bed") )
      bed_tib <- cur_tib %>% 
        dplyr::distinct( !!chr_sym, SNP_Range_Beg, SNP_Range_End ) %>%
        dplyr::arrange( !!chr_sym, SNP_Range_Beg )
      readr::write_tsv( bed_tib, cur_bed, col_names = FALSE )
      
      ret_tib <- ret_tib %>% dplyr::bind_rows( cur_tib )
    }
    
    non_cpg_fwd_cnt <- ret_tib %>% 
      dplyr::filter(Din_Fwd != "CG") %>% base::nrow()
    non_cpg_top_cnt <- ret_tib %>% 
      dplyr::filter(Din_Top != "CG") %>% base::nrow()
    
    if ( vb >= vt+1 )
      cat( glue::glue("{mssg} non_cpg_fwd_cnt='{non_cpg_fwd_cnt}'.{RET}") )

    # This calculation makes assumptions on existence of Forward Sequence...
    #
    # mis_fwd_cnt <- prb_dat_tib %>% 
    #   dplyr::inner_join( gen_fwd_tib, by=c("Ilmn_ID") ) %>% 
    #   dplyr::filter( Forward_Sequence != Forward_Sequence_Genome ) %>%
    #   base::nrow()
    # 
    # if ( vb >= 3 )
    #   cat( glue::glue("{mssg} non_cpg_cnt='{non_cpg_cnt}'.{RET}",
    #                   "{mssg} mis_fwd_cnt='{mis_fwd_cnt}'.{RET2}") )

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

fas_to_dnaset = function(x,
                         vb=0, vt=3, tc=1, tt=NULL,
                         fun_tag='fas_to_dnaset') {
  
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
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  f_time <- base::system.time({
    
    seq_set <- NULL
    
    if ( purrr::is_character(x) && file.exists(x) ) {
      if ( vb >= vt+1 )
        cat(glue::glue("{mssg} Loading Fasta File: '{x}'...{RET}"))
      
      errs_mssg <- glue::glue("File '{x} does not exist'")
      if ( !file.exists(x) ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
      
      seq_set <- Biostrings::readDNAStringSet(filepath = x, format = "fasta" )
      
    } else {
      if ( vb >= vt+1 )
        cat(glue::glue("{mssg} DNA String Set Already Loaded!{RET}"))
      
      seq_set <- x
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-dna-seq-set")
    ret_cnt <- print_tib( seq_set, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  seq_set
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  Filtering, Subsets and Modification 
#                               Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

parse_short_seq = function(chr_seq,
                           sub_seq = "CG",
                           
                           flank_len = 0,
                           remove_ns = TRUE,
                           
                           vb=0, vt=3, tc=1, tt=NULL,
                           fun_tag='parse_short_seq') {
  
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
    cat(glue::glue("{mssg}    sub_seq = '{sub_seq}'.{RET}"))
    cat(glue::glue("{mssg}  flank_len = '{flank_len}'.{RET}"))
    cat(glue::glue("{mssg}  remove_ns = '{remove_ns}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  f_time <- base::system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      PARSE:: Coordinates of Sub-Seq 
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    pos_seqs <- NULL
    sub_seq_len <- sub_seq %>% length()
    if ( vb >= vt ) cat(glue::glue("{mssg} sub_seq_len = '{sub_seq_len}'{RET}"))
    
    for (sseq in sub_seq) {
      if ( vb >= vt+3 ) cat(glue::glue("{mssg} sseq = '{sseq}' = "))
      
      cur_seqs <- Biostrings::matchPattern( pattern = sseq, subject = chr_seq )
      cur_slen <- length( cur_seqs@ranges@start )
      
      if ( vb >= vt+3 ) cat( glue::glue("'{cur_slen}'{RET}") )
      
      if ( cur_slen > 0 ) {
        if ( is.null(pos_seqs) ) {
          pos_seqs = cur_seqs
        } else {
          pos_seqs <- c( pos_seqs, cur_seqs )
        }
        if ( vb >= vt+4 ) print( pos_seqs )
      }
    }
    if ( vb >= vt+4 ) print(pos_seqs)
    if ( vb >= vt+4 ) cat( glue::glue("pos_seq above!!!{RET2}{RET2}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      FILTER:: Templates Overhanging 
    #                Beginning and End of Chromosome Sequence
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if ( flank_len != 0 ) {
      
      chr_len <- chr_seq@length - flank_len
      if ( vb >= vt+4 ) cat( glue::glue("Chromosome Length = '{chr_len}'{RET}") )
      
      if (purrr::is_list(pos_seqs)) pos_seqs <- pos_seqs %>% unlist()
      
      ups_idx <- which( pos_seqs@ranges@start < flank_len + 1 )
      ups_cnt <- length(ups_idx)
      if ( ups_cnt != 0 ) pos_seqs <- pos_seqs[ -ups_idx ]
      if ( vb >= vt+4 ) cat( glue::glue("Up Counts = '{ups_cnt}'{RET}") )
      
      dns_idx <- which( pos_seqs@ranges@start+pos_seqs@ranges@width > chr_len )
      dns_cnt <- length(dns_idx)
      if ( dns_cnt != 0 ) pos_seqs <- pos_seqs[ -dns_idx ]
      if ( vb >= vt+4 ) cat( glue::glue("Down Counts = '{dns_cnt}'{RET}") )
      
      if ( vb >= vt+3 ) {
        cat(glue::glue(
          "{mssg} Created IRanges(pos_seqs), flank_len={flank_len}.{RET}",
          "{mssg}{TAB} Removed {ups_cnt} upstream short flanks.{RET}",
          "{mssg}{TAB} Removed {dns_cnt} downstream short flanks.{RET}") )
        ret_key <- glue::glue("filtered-pos_seqs")
        ret_cnt <- print_tib( pos_seqs, fun_tag = fun_tag, name = ret_key, 
                              vb=vb,vt=vt+3,tc=tc+1,tt=tt )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #           EXPAND:: Extend Flanking Length on Both Up/DownStream 
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      fwd_seqs <- pos_seqs + flank_len
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                 FILTER:: Remove Any Sequences with N's 
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( remove_ns ) {
        pass_idx <- base::as.character(fwd_seqs) %>% 
          stringr::str_detect("N", negate = TRUE)
        pass_cnt <- length( pass_idx )
        fail_cnt <- pos_seqs@subject@length - pass_cnt
        fwd_seqs <- fwd_seqs[ pass_idx, ]
      }
      
      warn_mssg <- glue::glue("Found zero acceptable {sub_seq} templates")
      if ( pass_cnt == 0 ) wflag <- TRUE
      if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
      wflag <- FALSE
      
      if ( vb >= vt+3 ) {
        cat(glue::glue(
          "{mssg} Expanding IRanges(pos_seqs) by {flank_len} = {RET}",
          "{mssg}{TAB} Removed {fail_cnt} Templates with N's.{RET}",
          "{mssg}{TAB} Final Passing Count = {pass_cnt}!{RET2}") )
        print( fwd_seqs )
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-dna-seq-set")
    ret_cnt <- print_tib( fwd_seqs, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  fwd_seqs
}



# End of file
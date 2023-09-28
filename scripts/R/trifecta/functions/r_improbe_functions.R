
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   improbe (Infinium Methylation) Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Infinium Methylation Probe Design Methods::
#                        r-improbe re-implemented
#                    Allows All Probe Designs:: cg,ch,rs 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Infinium Methylation Probe Design Methods::
#                           r-improbe re-implemented
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

r_improbe = function(tib,
                     
                     ids_key, 
                     seq_key, 
                     din_key,
                     
                     top_col = NULL,
                     top_key = NULL,
                     
                     srsplit=FALSE,
                     srd_key=NULL,
                     srd_str='FR',
                     
                     cosplit=FALSE,
                     cos_key=NULL,
                     cos_str='CO',
                     
                     ups_len = 60, 
                     seq_len = 122,
                     
                     join     = FALSE,
                     subset   = FALSE,
                     sub_cols = NULL,
                     
                     add_matseq = TRUE,
                     add_topseq = FALSE,
                     
                     del='_',
                     max=0,
                     
                     out_dir,
                     run_tag,
                     pre_tag = NULL,
                     
                     reload     = 0,
                     reload_min = 2,
                     ret_data   = FALSE,
                     parallel   = FALSE,
                     
                     vb=0, vt=3, tc=1, tt=NULL,
                     fun_tag='r_improbe') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
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
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      ids_key = '{ids_key}'.{RET}"))
    cat(glue::glue("{mssg}      seq_key = '{seq_key}'.{RET}"))
    cat(glue::glue("{mssg}      din_key = '{din_key}'.{RET}"))
    cat(glue::glue("{mssg}      top_col={top_col}.{RET}"))
    cat(glue::glue("{mssg}      top_key={top_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      srsplit={srsplit}.{RET}"))
    cat(glue::glue("{mssg}      srd_key={srd_key}.{RET}"))
    cat(glue::glue("{mssg}      srd_str={srd_str}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      cosplit={cosplit}.{RET}"))
    cat(glue::glue("{mssg}      cos_key={cos_key}.{RET}"))
    cat(glue::glue("{mssg}      cos_str={cos_str}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      ups_len={ups_len}.{RET}"))
    cat(glue::glue("{mssg}      seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}         join={join}.{RET}"))
    cat(glue::glue("{mssg}       subset={subset}.{RET}"))
    cat(glue::glue("{mssg}     sub_cols={sub_cols}.{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, aux_csv, out_csv, end_txt) )
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if ( ( srsplit && !is.null(srd_key) ) || 
       ( cosplit && !is.null(cos_key) ) ) {
    
    if ( vb >= vt ) 
      cat(glue::glue("{mssg} Will split by strand...{RET}"))
    
    if (srsplit &&!is.null(srd_key)) {
      tib_list <- split(tib, f=dplyr::pull(tib, srd_key))
      
      for (srd in names(tib_list)) {
        if ( vb >= vt ) 
          cat(glue::glue("{RET}{BRK}{RET}{mssg} Starting::split={tc}/{srd}{RET2}"))
        
        if ( vb >= vt+4 )
          tib_list[[srd]] %>% dplyr::select(
            dplyr::all_of(c(!!ids_key,!!seq_key,!!din_key,
                            !!srd_key,!!cos_key))) %>% print()
        
        cur_tib <- NULL
        cur_tib <- r_improbe(tib = tib_list[[srd]], 
                             
                             srd_str = srd_str, 
                             cos_str = cos_str,
                             
                             ids_key = ids_key,
                             seq_key = seq_key,
                             din_key = din_key,
                             
                             top_col = NULL,
                             top_key = NULL,
                             
                             srsplit = FALSE,
                             srd_key = srd_key,
                             cosplit = cosplit,
                             cos_key = cos_key,
                             
                             ups_len = ups_len, 
                             seq_len = seq_len, 
                             
                             join = join,
                             subset = subset,
                             sub_cols = sub_cols,
                             
                             add_matseq = add_matseq,
                             add_topseq = add_topseq,
                             
                             del = del,
                             max = max,
                             
                             out_dir = out_dir,
                             run_tag = run_tag,
                             pre_tag = pre_tag,
                             
                             reload     = reload,
                             reload_min = reload_min,
                             ret_data   = ret_data,
                             parallel   = parallel,
                             
                             vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
        
        if ( vb >= vt ) 
          cat(glue::glue("{RET}{mssg} DONE::split={srd}{RET}{BRK}{RET}"))
      }
      
    } else if (cosplit && !is.null(cos_key)) {
      tib_list <- split(tib, f=dplyr::pull(tib, cos_key))
      
      for (srd in names(tib_list)) {
        if ( vb >= vt ) 
          cat(glue::glue("{RET}{BRK}{RET}{mssg} ",
                         "Starting::split={tc}/{srd}{RET2}"))
        if ( vb >= vt+4 )
          tib_list[[srd]] %>% dplyr::select(
            dplyr::all_of(c(!!ids_key,!!seq_key,!!din_key,
                            !!srd_key,!!cos_key))) %>% print()
        
        cur_tib <- NULL
        cur_tib <- r_improbe(tib = tib_list[[srd]], 
                             srd_str = srd_str, 
                             cos_str = cos_str,
                             
                             ids_key = ids_key, 
                             seq_key = seq_key, 
                             din_key = din_key,
                             
                             top_col = NULL,
                             top_key = NULL,
                             
                             srsplit = srsplit,
                             srd_key = srd_key,
                             cosplit = FALSE,
                             cos_key = cos_key,
                             
                             ups_len = ups_len, 
                             seq_len = seq_len, 
                             
                             join = join,
                             subset = subset,
                             sub_cols = sub_cols,
                             
                             add_matseq = add_matseq,
                             add_topseq = add_topseq,
                             
                             del = del,
                             max = max,
                             
                             out_dir = out_dir,
                             run_tag = run_tag,
                             pre_tag = pre_tag,
                             
                             reload     = reload,
                             reload_min = reload_min,
                             ret_data   = ret_data,
                             parallel   = parallel,
                             
                             vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
        
        if ( vb >= vt ) 
          cat(glue::glue("{RET}{mssg} DONE::split={srd}{RET}{BRK}{RET}"))
      }
    }
  } else {
    
    if ( vb >= vt ) 
      cat(glue::glue("{mssg} Running r-improbe...{RET2}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    
    e_time <- 0
    f_time <- 0
    
    ftime <- system.time({
      
      # if (TRUE) {
      # } else if (FALSE) {
      #   warn_mssg <- glue::glue("WARN_MESSAGE")
      #   cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
      # } else {
      # fail_mssg <- glue::glue("FAIL_MESSAGE")
      # stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
      # return(NULL)
      # }
      
      valid_srd_vec <- c("FR","TB")
      if (!srd_str %in% valid_srd_vec) {
        fail_mssg <- glue::glue("Invalid SR string={srd_str}")
        stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
        return(NULL)
      }
      valid_cos_vec <- c("CO")
      if (!cos_str %in% valid_cos_vec) {
        fail_mssg <- glue::glue("Invalid CO string={cos_str}")
        stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
        return(NULL)
      }
      
      if (purrr::is_character(tib))
        tib <- safe_read(file = tib, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      if (max>0) {
        if ( vb >= vt ) 
          cat(glue::glue("{mssg} Will subset input to max={max}.{RET}"))
        tib <- tib %>% head(n=max)
      }
      
      # Ambiguous Source Design Sequence Strand
      ids_sym <- rlang::sym(ids_key)
      seq_sym <- rlang::sym(seq_key)
      din_sym <- rlang::sym(din_key)
      
      srd_sym  <- rlang::sym(srd_str)
      cos_sym  <- rlang::sym(cos_str)
      
      srd_vec  <- stringr::str_split(srd_str, '', simplify=TRUE) %>% as.vector()
      cos_vec  <- stringr::str_split(cos_str, '', simplify=TRUE) %>% as.vector()
      
      tar_srd_vec <- srd_vec
      if (srsplit && !is.null(srd_key))
        tar_srd_vec <- tib %>% dplyr::pull(!!srd_key) %>% unique() %>% dplyr::pull()
      
      tar_cos_vec <- cos_vec
      if (cosplit && !is.null(cos_key))
        tar_cos_vec <- tib %>% dplyr::pull(!!cos_key) %>% unique() %>% dplyr::pull()
      
      tar_srd_vec <- expand.grid(tar_srd_vec, tar_cos_vec) %>% 
        tibble::as_tibble() %>% 
        tidyr::unite(srd, Var1,Var2, sep='', remove=TRUE)
      tar_srd_vec <- tar_srd_vec %>% dplyr::pull()
      
      #
      #
      #
      # Do we still need the sill Seq_ID and PRB_DES???
      #
      #
      #
      src_man_tib <- tib %>% 
        dplyr::select(!!ids_sym, !!din_sym, !!seq_sym) %>%
        dplyr::mutate(Seq_ID:=!!ids_sym, PRB_DES:=!!din_sym)
      
      ret_key <- glue::glue("src_man_tib-1({fun_tag})")
      ret_cnt <- print_tib(src_man_tib, fun_tag = fun_tag, name = ret_key,
                           vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      # Ensure we have 122 mer format 60[NN]60
      src_man_tib <- validate_templates(tib = src_man_tib,
                                        seq_key = seq_key,
                                        ups_len = ups_len,
                                        seq_len = seq_len,
                                        
                                        vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #
      #          Build Ref/BSC Template Sequences for Target Strands::
      #                              bsc_templates()
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( vb >= vt ) 
        cat(glue::glue("{mssg} Building forward & reverse ",
                       "reference template sequneces...{RET}"))
      
      bsc_tibs <- NULL
      bsc_tibs <- c(
        bsc_templates(tib=src_man_tib, 
                      srd_chr="F",
                      srd_vec=tar_srd_vec, 
                      seq_key=seq_key,
                      
                      srd_bol=TRUE,
                      srd_str=srd_str,  # Previously "SR"
                      cos_bol=TRUE,
                      cos_str=cos_str,  # Previously "CO"
                      
                      vb=vb,vt=vt+4,tc=tc+1,tt=tt ),
        
        bsc_templates(tib=src_man_tib, 
                      srd_chr="R",
                      srd_vec=tar_srd_vec, 
                      seq_key=seq_key, 
                      
                      srd_bol=FALSE,
                      srd_str=srd_str,  # Previously "SR"
                      cos_bol=TRUE,
                      cos_str=cos_str,  # Previously "CO"
                      
                      vb=vb,vt=vt+4,tc=tc+1,tt=tt ) )
      
      srd_names <- names(bsc_tibs)
      srd_count <- srd_names %>% length()
      
      if ( vb >= vt )
        cat(glue::glue("{mssg} Strand Names({srd_names})={RET}"))
      if ( vb >= vt ) 
        cat(glue::glue("{mssg} Done. Building forward & reverse ",
                       "reference template sequneces(strands={srd_count})!",
                       "{RET2}"))
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                    Build all Probes on Each Strand::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (parallel) {
        if ( vb >= vt ) 
          cat(glue::glue("{mssg}{TAB} Building probes for each ",
                         "strand (Parallel)...{RET}"))
        
        ret_tib <- foreach (srd=srd_names, .combine=rbind) %dopar% {
          lapply(split(bsc_tibs[[srd]], dplyr::pull(bsc_tibs[[srd]],din_key)), 
                 des_all_prbs, srd_str=srd_str, cos_str=cos_str, din_key=din_key,
                 vb=vb,vt=vt+4,tc=tc+1,tt=tt ) %>% dplyr::bind_rows()
        }
      } else {
        if ( vb >= vt ) 
          cat(glue::glue("{mssg}{TAB} Building probes for each ",
                         "strand (Linear)...{RET}"))
        
        for (srd in srd_names) {
          ret_tib <- ret_tib %>% dplyr::bind_rows(
            lapply(split(bsc_tibs[[srd]], dplyr::pull(bsc_tibs[[srd]],din_key)), 
                   des_all_prbs, srd_str=srd_str, cos_str=cos_str, din_key=din_key,
                   vb=vb,vt=vt+4,tc=tc+1,tt=tt ) %>% 
              dplyr::bind_rows() )
        }
      }
      if ( vb >= vt ) 
        cat(glue::glue("{mssg}{TAB} Done. Building probes for ",
                       "each strand.{RET2}"))
      ret_key <- glue::glue("prbs-on-all-srds:ret-tib({fun_tag})")
      ret_cnt <- print_tib(ret_tib,fun_tag = fun_tag, name = ret_key,
                           vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      # Update Keys::
      if ( vb >= vt ) 
        cat(glue::glue("{mssg} Updating keys and strands.{RET}"))
      
      ids_unq_key <- paste(ids_key,"srd", sep=del)
      ids_unq_sym <- rlang::sym(ids_unq_key)
      
      sr_out_str <- paste("Strand",srd_str, sep=del)
      co_out_str <- paste("Strand",cos_str, sep=del)
      sr_out_sym <- rlang::sym(sr_out_str)
      co_out_sym <- rlang::sym(co_out_str)
      
      prb1U_str  <- "Prb_1U"
      prb1M_str  <- "Prb_1M"
      prb2D_str  <- "Prb_2D"
      
      prb1U_sym  <- rlang::sym(prb1U_str)
      prb1M_sym  <- rlang::sym(prb1M_str)
      prb2D_sym  <- rlang::sym(prb2D_str)
      
      ret_tib <- ret_tib %>% 
        dplyr::mutate(
          !!sr_out_sym:=case_when(!!srd_sym ~ srd_vec[1], 
                                  !(!!srd_sym) ~ srd_vec[2], 
                                  TRUE ~ NA_character_),
          !!co_out_sym:=case_when(!!cos_sym ~ cos_vec[1], 
                                  !(!!cos_sym) ~ cos_vec[2], 
                                  TRUE ~ NA_character_),
          !!ids_unq_sym:=paste(!!ids_sym,
                               paste0(!!sr_out_sym,!!co_out_sym), sep=del)
        ) %>% dplyr::arrange(ids_unq_key)
      ret_key <- glue::glue("Updated-keys/srds:ret-tib({fun_tag})")
      ret_cnt <- print_tib(ret_tib,fun_tag = fun_tag, name = ret_key,
                           vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      sub_cols <- c(sub_cols,ids_unq_key,ids_key,din_key,
                    sr_out_str,co_out_str,seq_key)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                    Add Formatted Match Probe Sequences::
      #                       Upper Case, Non-cryptic Names....
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (add_matseq) {
        if ( vb >= vt ) 
          cat(glue::glue("{mssg}{TAB} Formatting Match Probe Sequences...{RET}"))
        
        ret_tib <- ret_tib %>% 
          dplyr::mutate(!!prb1U_sym:=stringr::str_to_upper(PRB1_U),
                        !!prb1M_sym:=stringr::str_to_upper(PRB1_M),
                        !!prb2D_sym:=stringr::str_to_upper(PRB2_D) )
        
        sub_cols <- c(sub_cols, prb1U_str,prb1M_str,prb2D_str)
      } else {
        sub_cols <- c(sub_cols, "PRB1_U","PRB1_M","PRB2_D")
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                       Add Top-Sequence Designation::
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (add_topseq && !is.null(top_col)) {
        if ( vb >= vt ) 
          cat(glue::glue("{mssg}{TAB} Adding Top-Sequences Designation...{RET}"))
        
        ret_tib <- ret_tib %>% 
          set_topbot_tib(seq_key=seq_key, 
                         srd_key=top_col,
                         top_key=top_key, 
                         
                         vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      }
      
      # Order and subset output::
      #  TBD:: Clean the sub_cols/subset up. Pretty messy right now...
      #
      if (subset) {
        old_cols <- ret_tib %>% dplyr::select(
          dplyr::starts_with("NXB_"), dplyr::starts_with("TAR_") ) %>%
          names() %>% as.vector()
        new_cols <- old_cols %>%
          stringr::str_replace("^NXB_","Nxb_") %>%
          stringr::str_replace("^TAR_","Ext_") %>% as.vector()
        sub_cols <- c(sub_cols, new_cols) %>% unique()
        
        ret_tib <- ret_tib %>% 
          dplyr::rename_with( ~ new_cols, dplyr::all_of(old_cols) ) %>%
          dplyr::select( dplyr::any_of(sub_cols) )
        
      } else {
        sub_cols <- sub_cols %>% c("NXB_U","NXB_M","NXB_D",
                                   "CPN_U","CPN_M","CPN_D",
                                   "TAR_U","TAR_M","TAR_D") %>% unique()
        
        ret_tib <- ret_tib %>% 
          dplyr::select(dplyr::any_of(sub_cols), dplyr::everything())
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                               Write Data::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # safe_write( sum_csv, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      # safe_write( aux_csv, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      # safe_write( out_csv, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
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
    e_time <- as.double(f_time[3]) %>% round(2)
    if (!is.null(tt)) tt$addTime(f_time,fun_tag)
    if (vb >= vt) cat(glue::glue(
      "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
    
    ret_tib
  }
}

validate_templates = function(tib,
                              seq_key,
                              ups_len=60, 
                              seq_len=122, 
                              del="_",
                              pad="N",
                              
                          vb=0, vt=6, tc=1, tt=NULL,
                          fun_tag='validate_templates') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}   seq_key={seq_key}.{RET}"))
    cat(glue::glue("{mssg}   ups_len={ups_len}.{RET}"))
    cat(glue::glue("{mssg}   seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  e_time  <- 0
  f_time  <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  seq_sym <- rlang::sym(seq_key)
  
  if ( vb >= vt ) 
    cat(glue::glue("{mssg} Validating design sequneces...{RET}"))
  
  ret_tib <- tib %>%
    dplyr::mutate(
      !!seq_sym := stringr::str_replace(!!seq_sym, '\\[','_') %>% 
        stringr::str_replace('\\]','_')) %>%
    tidyr::separate(!!seq_sym, 
                    into=c("PreSeqN", "MidSeqN", "PosSeqN"), sep=del) %>%
    dplyr::mutate(
      PreSeqN=stringr::str_sub(PreSeqN,   -ups_len),
      PosSeqN=stringr::str_sub(PosSeqN, 1, ups_len),
      PreSeqN=stringr::str_pad(string=PreSeqN, 
                               width=ups_len, side='left', pad=pad),
      PosSeqN=stringr::str_pad(string=PosSeqN, 
                               width=ups_len, side='right', pad=pad),
      DesNucA=stringr::str_sub(MidSeqN, 1,1), 
      DesNucB=stringr::str_sub(MidSeqN, 2,2),
      !!seq_sym :=paste0(PreSeqN,'[',MidSeqN,']',PosSeqN) )
  
  if ( vb >= vt ) 
    cat(glue::glue("{mssg} Done. Validating design sequences.{RET2}"))
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

bsc_templates = function(tib,
                         srd_chr,
                         srd_vec,
                         seq_key,
                         
                         srd_bol,
                         srd_str="SR",
                         cos_bol,
                         cos_str="CO",
                         
                         vb=0, vt=6, tc=1, tt=NULL,
                         fun_tag='bsc_templates') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}   srd_chr={srd_chr}.{RET}"))
    cat(glue::glue("{mssg}   srd_vec={srd_vec}.{RET}"))
    cat(glue::glue("{mssg}   seq_key={seq_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}   srd_bol={srd_bol}.{RET}"))
    cat(glue::glue("{mssg}   srd_str={srd_str}.{RET}"))
    cat(glue::glue("{mssg}   cos_bol={cos_bol}.{RET}"))
    cat(glue::glue("{mssg}   cos_str={cos_str}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  e_time  <- 0
  f_time  <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ret_key <- glue::glue("Initial-tib")
  ret_cnt <- print_tib( tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  srdC <- paste0(srd_chr,"C")
  srdO <- paste0(srd_chr,"O")
  
  seq_sym <- rlang::sym(seq_key)
  srd_sym <- rlang::sym(srd_str)
  cos_sym <- rlang::sym(cos_str)
  
  des_seq_tib <- NULL
  if (srdC %in% srd_vec || srdO %in% srd_vec)
    des_seq_tib <- tib %>% dplyr::mutate(
      !!srd_sym:=srd_bol,
      !!cos_sym:=cos_bol, 
      DesSeqN=shearBrac(!!seq_sym) )
  
  ret_key <- glue::glue("des_seq_tib")
  ret_cnt <- print_tib( des_seq_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  # BSC-[Forward|Reverse]-Converted::
  #
  if (!is.null(des_seq_tib)) {
    if ( vb >= vt ) 
      cat(glue::glue("{mssg} Building bisulfite convertering ",
                     "forward strand template sequence(s)...{RET}"))
    
    ret_dat[[srdC]] <- des_seq_tib %>% 
      dplyr::mutate(
        DesBscU = bscUs(DesSeqN),
        DesBscM = bscMs(DesSeqN),
        DesBscD = bscDs(DesSeqN) )
    
    if ( vb >= vt+1 ) 
      cat(glue::glue("{mssg} Done. Building bisulfite ",
                     "template sequences ({srdC}).{RET}"))
    
    ret_key <- glue::glue("bsc_tibs-{srdC}")
    ret_cnt <- print_tib( ret_dat[[srdC]], fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # BSC-[Forward|Reverse]-Opposite::
    #
    if (srdO %in% srd_vec) {
      ret_dat[[srdO]] <- ret_dat[[srdC]] %>% 
        dplyr::mutate(
          !!srd_sym:=srd_bol,
          !!cos_sym:=!cos_bol,
          DesBscU=revCmp(DesBscU),
          DesBscM=revCmp(DesBscM),
          DesBscD=revCmp(DesBscD) )
      
      if ( vb >= vt+1 )
        cat(glue::glue("{mssg} Done. Building bisulfite ",
                       "template sequences ({srdO}).{RET}"))
      
      ret_key <- glue::glue("bsc_tibs-{srdO}")
      ret_cnt <- print_tib( ret_dat[[srdO]], fun_tag = fun_tag, name = ret_key, 
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )
      
    }
    if (!srdC %in% srd_vec) ret_dat[[srdC]] <- NULL
    if ( vb >= vt+1 ) cat(glue::glue("{RET}"))
  }
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

des_all_prbs = function(tib,
                        srd_str="SR",
                        cos_str="CO",
                        din_key,
                        
                        vb=0, vt=6, tc=1, tt=NULL,
                        fun_tag='des_all_prbs') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}   srd_str={srd_str}.{RET}"))
    cat(glue::glue("{mssg}   cos_str={cos_str}.{RET}"))
    cat(glue::glue("{mssg}   din_key={din_key}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  e_time  <- 0
  f_time  <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  srd_sym <- rlang::sym(srd_str)
  cos_sym <- rlang::sym(cos_str)
  din_sym <- rlang::sym(din_key)
  
  fr <- tib %>% dplyr::distinct(!!srd_sym) %>% base::as.logical()
  co <- tib %>% dplyr::distinct(!!cos_sym) %>% base::as.logical()
  pr <- tib %>% dplyr::distinct(!!din_sym) %>% base::as.character()
  
  ret_tib <- dplyr::bind_rows(
    tib %>% 
      des_to_prbs(fwd=fr, con=co, pr=pr, mu='N', des_seq='DesSeqN', 
                  vb=vb,vt=vt+3,tc=tc+1,tt=tt ) %>%
      des_to_prbs(fwd=fr, con=co, pr=pr, mu='U', des_seq='DesBscU', 
                  vb=vb,vt=vt+3,tc=tc+1,tt=tt ) %>%
      des_to_prbs(fwd=fr, con=co, pr=pr, mu='M', des_seq='DesBscM', 
                  vb=vb,vt=vt+3,tc=tc+1,tt=tt ) %>%
      des_to_prbs(fwd=fr, con=co, pr=pr, mu='D', des_seq='DesBscD', 
                  vb=vb,vt=vt+3,tc=tc+1,tt=tt ) )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# Function to design probes from single design strand and orientations::
#   TBD:: Previous default was 'QC_CPN=TRUE' Not sure if that is needed...
#
des_to_prbs = function(tib, 
                       fwd,
                       con, 
                       pr, 
                       mu, 
                       des_seq='DesSeqN', 
                       len=48, 
                       del='_',
                       QC_CPN=FALSE,
                       
                       vb=0, vt=6, tc=1, tt=NULL,
                       fun_tag='des_to_prbs') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}       pr={pr}.{RET}"))
    cat(glue::glue("{mssg}       mu={mu}.{RET}"))
    cat(glue::glue("{mssg}      fwd={fwd}.{RET}"))
    cat(glue::glue("{mssg}      con={con}.{RET}"))
    cat(glue::glue("{mssg}  des_seq={des_seq}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  e_time  <- 0
  f_time  <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  stopifnot(is.logical(fwd))
  stopifnot(is.logical(con))
  
  if (!is.logical(fwd)) {
    fail_mssg <- glue::glue("fwd={fwd} Must be logical")
    stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    return(NULL)
  }
  if (!is.logical(con)) {
    fail_mssg <- glue::glue("con={con} Must be logical")
    stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    return(NULL)
  }
  
  if (mu!='N' && mu!='U' && mu!='M' && mu!='D') {
    fail_mssg <- glue::glue("mu={mu} Only Supported=[N,U,M,D]")
    stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    return(NULL)
  }
  if (pr!='cg' && pr!='ch' && pr!='rs' && pr!='rp' && pr!='mu' && pr!='bc') {
    fail_mssg <- glue::glue("pr={pr} Only Supported=[cg,ch,rp,rs]")
    stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    return(NULL)
    
    
    stop(glue::glue("{RET}{mssg} ERROR: pr={pr} Only ",
                    "Supported=[cg,ch,rp,rs]!{RET2}"))
    return(ret_tib)
  }
  
  des_seq_sym <- rlang::sym(des_seq)
  if (pr=='rs') {
    if      ( fwd &&  con) nxb_pos <- 60
    else if (!fwd &&  con) nxb_pos <- 61
    else if ( fwd && !con) nxb_pos <- 61
    else if (!fwd && !con) nxb_pos <- 60
    else {
      fail_mssg <- glue::glue("Unsupported combination fwd={fwd}, con={con}")
      stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
      return(NULL)
    }
  } else if (pr=='ch') {
    # Originally this was identical to rs format above, but for forward 
    #  sequences needs to be shifted upstream for converted and downstream 
    #  for opposite::
    if      ( fwd &&  con) nxb_pos <- 61 # Previously = 60
    else if (!fwd &&  con) nxb_pos <- 61
    else if ( fwd && !con) nxb_pos <- 60 # Previously = 61
    else if (!fwd && !con) nxb_pos <- 60
    else {
      fail_mssg <- glue::glue("Unsupported combination fwd={fwd}, con={con}")
      stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
      return(NULL)
    }
    
    # NEW:: CpH Code::
    #
    if      ( fwd &&  con) nxb_pos <- 60 # Previously = 60
    else if (!fwd &&  con) nxb_pos <- 61
    else if ( fwd && !con) nxb_pos <- 61 # Previously = 61
    else if (!fwd && !con) nxb_pos <- 60
    else {
      fail_mssg <- glue::glue("Unsupported combination fwd={fwd}, con={con}")
      stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
      return(NULL)
    }
    
  } else if (pr=='cg' || pr=='rp' || pr=='mu' || pr=='bc' || 
             stringr::str_starts(pr,'ct')) {
    # $prb_F_C_I  = revCmp(substr($des_F_C, 60, 50));
    # $prb_F_O_I  = revCmp(substr($des_F_O, 61, 50));
    # $prb_R_C_I  = revCmp(substr($des_R_C, 60, 50));
    # $prb_R_O_I  = revCmp(substr($des_R_O, 61, 50));
    # 
    # $prb_F_C_II  = revCmp(substr($des_F_C, 61, 50));
    # $prb_F_O_II  = revCmp(substr($des_F_O, 62, 50));
    # $prb_R_C_II  = revCmp(substr($des_R_C, 61, 50));
    # $prb_R_O_II  = revCmp(substr($des_R_O, 62, 50));
    nxb_pos <- 60
    if (!con) nxb_pos <- 61
    
  } else {
    fail_mssg <- glue::glue("Probe_Type={pr} is currently NOT supported")
    stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    return(NULL)
  }
  cpg_pos <- nxb_pos + 1
  sec_pos <- cpg_pos + 1
  bod_pos <- sec_pos + 1
  end_pos <- bod_pos + len
  
  # Special consideration is needed for U/M strands at the query site. 
  #  For CN (i.e. cg or ch) this is actually done naturally in U/M conversion
  #  However, for  non-CN probes (i.e. rs) this needs to be forced to U/M
  #
  # This is handled by the TAR (Target/Query Nucleotide). This should only 
  #  change for U/M (QMAP_U/QMAP_M) for D its just itself.
  #
  ret_tib <- tib %>% dplyr::mutate(
    NXB=stringr::str_sub(!!des_seq_sym, nxb_pos, nxb_pos),
    CPN=stringr::str_sub(!!des_seq_sym, cpg_pos, cpg_pos),
    TAR=qmaps(CPN, mu=mu),
    SEC=stringr::str_sub(!!des_seq_sym, sec_pos, sec_pos),
    BOD=stringr::str_sub(!!des_seq_sym, bod_pos, end_pos-1),
    END=stringr::str_sub(!!des_seq_sym, end_pos, end_pos)
  )
  
  #  QC TEST:: for CpN (cg or ch) verify that the probes are equal. Well call
  #   this PRB0 (CGN) and PRB1 (TAR). After testing remove PRB0
  #
  if (QC_CPN && (pr=='cg')) {
    ret_tib <- ret_tib %>%
      tidyr::unite(PRB0, CPN,SEC,BOD, sep='', remove=FALSE) %>%
      tidyr::unite(PRB1, TAR,SEC,BOD, sep='', remove=FALSE) %>%
      tidyr::unite(PRB2, SEC,BOD,END, sep='', remove=FALSE) %>%
      dplyr::mutate(PRB0=revCmp(PRB0), PRB1=revCmp(PRB1), PRB2=revCmp(PRB2))
    
    qc_tib <- ret_tib %>% filter(PRB0!=PRB1)
    qc_len <- qc_tib %>% base::nrow()
    if (qc_len != 0) {
      qc_tib %>% dplyr::select(1,PRB0,PRB1) %>% print()
      fail_mssg <- glue::glue("pr={pr}, qc_len={qc_len} != 0")
      stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
      return(NULL)
    }
  } else {
    ret_tib <- ret_tib %>%
      tidyr::unite(PRB1, TAR,SEC,BOD, sep='', remove=FALSE) %>%
      tidyr::unite(PRB2, SEC,BOD,END, sep='', remove=FALSE) %>%
      dplyr::mutate(PRB1=revCmp(PRB1), PRB2=revCmp(PRB2))
  }
  
  # Add suffix to sequences for merging later
  ret_tib <- ret_tib %>%
    dplyr::select(PRB1,PRB2, NXB,CPN,TAR,BOD,END, everything()) %>%
    dplyr::rename(!!paste('PRB1',mu, sep=del):=PRB1,
                  !!paste('PRB2',mu, sep=del):=PRB2,
                  !!paste('NXB', mu, sep=del):=NXB,
                  !!paste('CPN', mu, sep=del):=CPN,
                  !!paste('TAR', mu, sep=del):=TAR,
                  !!paste('SEC', mu, sep=del):=SEC,
                  !!paste('BOD', mu, sep=del):=BOD,
                  !!paste('END', mu, sep=del):=END)
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

des2prbsNOTES = function(srd, desSeq,
                         verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'des2prbs'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} Starting...{RET}"))
  
  # bscNewU <- bscNewU(desSeq)
  # bscNewM <- bscNewM(desSeq)
  # bscNewD <- bscNewD(desSeq)
  
  # my @desSetU = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewU);
  # my @desSetM = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewM);
  # my @desSetD = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewD);
  
  #    my @desSetU=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewU($desSeq,$retUC));                                                                                                             
  #    my @desSetM=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewM($desSeq,$retUC));                                                                                                             
  
  # my ($desNxbU, $desCpgU, $desSeqU, $desEndU);
  # my ($desNxbM, $desCpgM, $desSeqM, $desEndM);
  # my ($desNxbD, $desCpgD, $desSeqD, $desEndD);
  
  #         (  $desNxbU,    $desCpgU,               $desSeqU,                $desEndU) =                                                                                                       
  #   return( $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7],                                                                                                      
  #           $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7]) if ($desCO eq $C);                                                                                   
  #          ( $desNxbM,    $desCpgM,               $desSeqM,                $desEndM) =                                                                                                       
  
  #                  ( $desNxbU,            $desCpgU,                    $desSeqU,                  $desEndU) =                                                                                
  #   return( revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]),                                                                              
  #           revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0])) if ($desCO eq $O);                                                           
  #                  ( $desNxbM,            $desCpgM,                    $desSeqM,                  $desEndM) =                                                                                
  
  
  # $$prbRef[$srd][$iU]=[ $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7] ];
  # $$prbRef[$srd][$iM]=[ $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7] ];
  # $$prbRef[$srd][$iD]=[ $desSetD[2], $desSetD[3], $desSetD[4].$desSetD[5].$desSetD[6], $desSetD[7] ];
  # 
  # $srd++;
  # $$prbRef[$srd][$iU]=[ revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]) ];
  # $$prbRef[$srd][$iM]=[ revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0]) ];
  # $$prbRef[$srd][$iD]=[ revCmpl($desSetD[4]), revCmpl($desSetD[3]), revCmpl($desSetD[1].$desSetD[2]), revCmpl($desSetD[0]) ];
  
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} Done.{RET2}"))
  
  NULL
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#          Infinium Methylation Probe toString/printing Methodss::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Replace Seq_ID and PRB_DES with symbolic links
#  Seq_ID::  seq_sym -> rlang::sym(seq_key)
#  PRB_DES:: din_sym -> rlang::sym(din_key)
#
print_prbs = function(tib, 
                      tar_des='cg',
                      ids_key,
                      prb_key,
                      des_key,
                      din_key,
                      org=NULL,
                      outDir, 
                      plotName, 
                      max=0,
                      verbose=0,vt=5,tc=1,tt=NULL, funcTag='print_prbs') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} Starting...{RET}"))
  if ( vb >= vt+1 +4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}   tar_des={tar_des}.{RET}"))
    cat(glue::glue("{mssg}   prb_key={prb_key}.{RET}"))
    cat(glue::glue("{mssg}   tar_des={tar_des}.{RET}"))
    cat(glue::glue("{mssg}   des_key={des_key}.{RET}"))
    cat(glue::glue("{mssg}   din_key={din_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}     outDir={outDir}.{RET}"))
    cat(glue::glue("{mssg}   plotName={plotName}.{RET2}"))
  }
  
  ids_sym <- rlang::sym(ids_key)
  prb_sym <- rlang::sym(prb_key)
  des_sym <- rlang::sym(des_key)
  din_sym <- rlang::sym(din_key)
  
  prb_mat_tibs <- tib %>% 
    dplyr::filter(des_key==tar_des) %>% 
    dplyr::distinct()
  prb_mat_cnt  <- prb_mat_tibs %>% base::nrow()
  
  plot_ord_tib <- prb_mat_tibs %>% 
    dplyr::distinct(!!ids_sym, .keep_all=TRUE)
  plot_ord_cnt <- plot_ord_tib %>% base::nrow()
  if ( vb >= vt+1 ) 
    cat(glue::glue("{mssg} plot_ord_cnt={plot_ord_cnt}, ",
                   "prb_mat_cnt={prb_mat_cnt}.{RET}"))
  
  if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
  
  passFile <- file.path(outDir, paste(plotName,tar_des,'pass.txt', sep='_') )
  failFile <- file.path(outDir, paste(plotName,tar_des,'fail.txt', sep='_') )
  unlink(passFile)
  unlink(failFile)
  
  if (!is.null(org)) org <- org %>% dplyr::distinct()
  
  tibs <- NULL
  for (ii in seq(1,plot_ord_cnt)) {
    tag_tib <- plot_ord_tib[ii,]
    cur_SeqID <- tag_tib %>% head(n=1) %>% dplyr::pull(ids_key)
    if ( vb >= vt+1 ) 
      cat(glue::glue("{mssg}{TAB} ii={ii}; cur_SeqID={cur_SeqID}.{RET}"))
    
    cur_org <- NULL
    org_str <- ''
    if (!is.null(org)) {
      cur_org <- org %>% 
        dplyr::filter(ids_key==cur_SeqID) %>% head(n=1)
      
      if (FALSE) {
        for (jj in c(1:length(cur_org))) {
          org_str[jj] <- stringr::str_c(cur_org[jj], collapse=',')
          break
        }
      }
      org_str <- cur_org %>% paste(collapse=', ')
    }
    
    cur_tib <- prb_mat_tibs %>% dplyr::filter(ids_key==cur_SeqID)
    tibs <- cur_tib %>% srds_to_brac()
    print(tibs)
    
    # Build Printable Strings for reach strand::
    strs <- NULL
    strs$F <- prbsToStr(tibs$F, pr=tag_tib[[des_key]][1], verbose=verbose, vt=6)
    strs$R <- prbsToStr(tibs$R, pr=tag_tib[[des_key]][1], verbose=verbose, vt=6)
    
    # Determine which file to write to::
    #   - Failed if both F/R strands have equal PRB1_U==PRB1_M
    #   - Passed otherwise
    fwd_cnt <- tibs$F %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
    rev_cnt <- tibs$R %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
    
    space_str <- paste(stringr::str_pad("# ", width=150, side="right", pad="-"), "\n", sep='')
    seq_id_str <- ''
    for (jj in c(1:length(org_str))) {
      
      seq_id_str <- 
        glue::glue("Probe_Type={tar_des}: ids_key: ",
                   paste(unique(cur_tib[[ids_key]]), collapse='\t'),
                   "; Original: {org_str[jj]}{RET}")
    }
    out_lines <- NULL
    out_lines[1] <- space_str
    out_lines[2] <- seq_id_str
    out_lines[3] <- stringr::str_c(strs$F)
    out_lines[4] <- stringr::str_c(strs$R)
    
    if ( vb >= vt+1 ) cat(out_lines)
    
    if (fwd_cnt>0 || rev_cnt>0) {
      readr::write_lines(out_lines, passFile, append=TRUE)
    } else {
      readr::write_lines(out_lines, failFile, append=TRUE)
    }
    
    # This was done for known comparisons::
    #
    if (FALSE) {
      if (FALSE) { # && cur_tib$Infinium_Design_Type[1]=='II') {
        cur_tib %>% dplyr::select(FR,CO, PRB2_D_IUP, PRB2_D_IMP,Man_MisMatch,Man_TarMatch,Bad_Design) %>% print()
      } else {
        prb1U_iup <- cur_tib %>% dplyr::select(FR,CO, PRB1_U_IUP, PRB1_U_IMP,Man_MisMatch,Man_TarMatch,Bad_Design)
        prb1M_iup <- cur_tib %>% dplyr::select(FR,CO, PRB1_M_IUP, PRB1_M_IMP,Man_MisMatch,Man_TarMatch,Bad_Design)
        # readr::write_file(str_c(prb1U_iup), passFile, append=TRUE)
        # readr::write_file(str_c(prb1M_iup), passFile, append=TRUE)
      }
    }
    
    if (max!=0 && ii>=max) break
  }
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} Done.{RET2}"))
  
  tibs
}

srds_to_brac = function(tib, 
                        beg1=1, end1=60, mid1=61,
                        beg2=63,end2=122,mid2=62) {
  # TBD:: Calculate all data points based on sequence length
  
  tib %>% 
    dplyr::mutate(StrandFR=case_when(FR ~ 'F', !FR ~ 'R', TRUE ~ NA_character_),
                  StrandCO=case_when(CO ~ 'C', !CO ~ 'O', TRUE ~ NA_character_),
                  DesSeqN=paste0(stringr::str_sub(DesSeqN,beg1,end1),{BNG},
                                 stringr::str_sub(DesSeqN,mid1,mid2),{BNG},
                                 stringr::str_sub(DesSeqN,beg2,end2)),
                  
                  DesBscU=paste0(stringr::str_sub(DesBscU,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscU,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscU,beg2,end2)),
                  
                  DesBscM=paste0(stringr::str_sub(DesBscM,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscM,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscM,beg2,end2)),
                  
                  DesBscD=paste0(stringr::str_sub(DesBscD,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscD,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscD,beg2,end2))) %>%
    dplyr::arrange(StrandFR, StrandCO) %>% split(.$StrandFR)
  
}

prbsToStr = function(tib,
                     pr, 
                     verbose=0, vt=5,tc=1,tt=NULL, funcTag='prbsToStr') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} Starting; pr={pr}.{RET}"))
  
  # fr1Key <- tib %>% dplyr::select(!!frKey) %>% head(n=1) %>% pull()
  # fr2Key <- tib$StrandFR[2]
  # co1Key <- tib$StrandCO[1]
  # co2Key <- tib$StrandCO[2]
  
  fr1Key <- tib$StrandFR[1]
  fr2Key <- tib$StrandFR[2]
  co1Key <- tib$StrandCO[1]
  co2Key <- tib$StrandCO[2]
  
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} fr1Key={fr1Key}, fr2Key={fr2Key}, co1Key={co1Key}, co2Key={co2Key}.{RET}"))
  
  mud <- list('U'='U','M'='M','D'='D')
  
  desSeq <- 'DesSeqN'
  bscKey <- lapply(mud, function(x) { paste('DesBsc',x,sep='')} )
  nxbKey <- lapply(mud, function(x) { paste('NXB',x,sep='_')} )
  cpnKey <- lapply(mud, function(x) { paste('CPN',x,sep='_')} )
  tarKey <- lapply(mud, function(x) { paste('TAR',x,sep='_')} )
  secKey <- lapply(mud, function(x) { paste('SEC',x,sep='_')} )
  bodKey <- lapply(mud, function(x) { paste('BOD',x,sep='_')} )
  endKey <- lapply(mud, function(x) { paste('END',x,sep='_')} )
  
  # Dertermine if the probes on this strand were designable::
  pass_cnt <- tib %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
  pass_str <- 'PASS: '
  if (pass_cnt==0) pass_str <- 'FAIL: '
  
  # TBD:: Note on the Opposite Strand we should reverse all framents, but currently fragLen==1 are left alone for effiecntcy...
  # Sketch Output::
  #
  # F_C_N    DesSeqN[CG]DesSeqN
  #
  #                  D2 22222
  #                N M1 1111
  #                N U1 1111
  # F_C_U    DesBscU[tG]DesBscU
  # F_O_U    DesBscU[Ca]DesBscU
  #             1111 1U N
  #             1111 1M N
  #            22222 2D
  #
  if (pr=='rs'||pr=='ch') {
    if (fr1Key=='F' && fr2Key=='F') {
      bufC <- 0
      bufO <- 0
      str <- glue::glue(
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 61+bufC), collapse=''),"{tib[[tarKey$D]][1]}{tib[[secKey$D]][1]}{BNG}{tib[[bodKey$D]][1]}{tib[[endKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 59+bufC), collapse=''),"{tib[[nxbKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[secKey$M]][1]}{BNG}{tib[[bodKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 59+bufC), collapse=''),"{tib[[nxbKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[secKey$U]][1]}{BNG}{tib[[bodKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_U {TAB}{tib[[bscKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_M {TAB}{tib[[bscKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{tib[[bscKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{tib[[desSeq]][1]}{RET}",
        # "FwdSeq{TAB}{tib$Forward_Sequence[1]}{RET}",
        "{RET}",
        "{pass_str}{fr1Key}_{co2Key}_N {TAB}{cmpl(tib[[desSeq]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_M {TAB}{Biostrings::reverse(tib[[bscKey$M]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_U {TAB}{Biostrings::reverse(tib[[bscKey$U]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 11-bufO), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][2])}{tib[[secKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[nxbKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 11-bufO), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][2])}{tib[[secKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[nxbKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 10-bufO), collapse=''),"{tib[[endKey$D]][2]}{Biostrings::reverse(tib[[bodKey$D]][2])}{tib[[secKey$D]][2]}{BNG}{tib[[tarKey$D]][2]}{RET}",
        "{RET}")
    } else if (fr1Key=='R' && fr2Key=='R') {
      bufO <- 0
      bufC <- 0
      str <- glue::glue(
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 61+bufO), collapse=''),"{tib[[tarKey$D]][2]}{tib[[secKey$D]][2]}{BNG}{tib[[bodKey$D]][2]}{tib[[endKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 59+bufO), collapse=''),"{tib[[nxbKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[secKey$M]][2]}{BNG}{tib[[bodKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 59+bufO), collapse=''),"{tib[[nxbKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[secKey$U]][2]}{BNG}{tib[[bodKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{tib[[bscKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_N {TAB}{revCmp(tib[[desSeq]][2])}{RET}",
        # "{fr2Key}_{co2Key}_N {TAB}{tib[[desSeq]][2]}{RET}",
        "{RET}",
        # "{fr1Key}_{co1Key}_N {TAB}{cmpl(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{Biostrings::reverse(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 11-bufC), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][1])}{tib[[secKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[nxbKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 11-bufC), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][1])}{tib[[secKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[nxbKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 10-bufC), collapse=''),"{tib[[endKey$D]][1]}{Biostrings::reverse(tib[[bodKey$D]][1])}{tib[[secKey$D]][1]}{BNG}{tib[[tarKey$D]][1]}{RET}",
        "{RET}")
    } else {
      stop(glue::glue("{RET}{mssg} ERROR: fr1Key={fr1Key}, fr2Key={fr2Key}, Allowed Values=[F,R]!{RET2}"))
    }
  } else {
    if (fr1Key=='F' && fr2Key=='F') {
      buf <- 0
      str <- glue::glue(
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 61+buf), collapse=''),"{tib[[tarKey$D]][1]}{tib[[secKey$D]][1]}{BNG}{tib[[bodKey$D]][1]}{tib[[endKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 59+buf), collapse=''),"{tib[[nxbKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[secKey$M]][1]}{BNG}{tib[[bodKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 59+buf), collapse=''),"{tib[[nxbKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[secKey$U]][1]}{BNG}{tib[[bodKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{tib[[bscKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{tib[[desSeq]][1]}{RET}",
        "{RET}",
        "{pass_str}{fr1Key}_{co2Key}_N {TAB}{cmpl(tib[[desSeq]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 11-buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][2])}{tib[[secKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[nxbKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 11-buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][2])}{tib[[secKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[nxbKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 10-buf), collapse=''),"{tib[[endKey$D]][2]}{Biostrings::reverse(tib[[bodKey$D]][2])}{tib[[secKey$D]][2]}{BNG}{tib[[tarKey$D]][2]}{RET}",
        "{RET}")
    } else if (fr1Key=='R' && fr2Key=='R') {
      buf <- 0
      buf <- 1
      str <- glue::glue(
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 61+buf), collapse=''),"{tib[[tarKey$D]][2]}{BNG}{tib[[secKey$D]][2]}{tib[[bodKey$D]][2]}{tib[[endKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 60+buf), collapse=''),"{tib[[nxbKey$M]][2]}{tib[[tarKey$M]][2]}{BNG}{tib[[secKey$M]][2]}{tib[[bodKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 60+buf), collapse=''),"{tib[[nxbKey$U]][2]}{tib[[tarKey$U]][2]}{BNG}{tib[[secKey$U]][2]}{tib[[bodKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{tib[[bscKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_N {TAB}{revCmp(tib[[desSeq]][2])}{RET}",
        # "{fr2Key}_{co2Key}_N {TAB}{tib[[desSeq]][2]}{RET}",
        "{RET}",
        # "{fr1Key}_{co1Key}_N {TAB}{cmpl(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{Biostrings::reverse(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 11+buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][1])}{BNG}{tib[[secKey$U]][1]}{tib[[tarKey$U]][1]}{BNG}{tib[[nxbKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 11+buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][1])}{BNG}{tib[[secKey$M]][1]}{tib[[tarKey$M]][1]}{BNG}{tib[[nxbKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 10+buf), collapse=''),"{tib[[endKey$D]][1]}{Biostrings::reverse(tib[[bodKey$D]][1])}{BNG}{tib[[secKey$D]][1]}{tib[[tarKey$D]][1]}{RET}",
        "{RET}")
    } else {
      stop(glue::glue("{RET}{mssg} ERROR: fr1Key={fr1Key}, fr2Key={fr2Key}, Allowed Values=[F,R]!{RET2}"))
    }
  }
  if ( vb >= vt+1 ) cat(str)
  if ( vb >= vt+1 ) cat(glue::glue("{mssg} Done.{RET2}"))
  
  str
}

prbsToStrMUD = function(tib, 
                        mu='U', 
                        verbose=0,vt=5,tc=1,tt=NULL, funcTag='prbsToStr') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  fr1Key <- tib$StrandFR[1]
  fr2Key <- tib$StrandFR[2]
  
  co1Key <- tib$StrandCO[1]
  co2Key <- tib$StrandCO[2]
  
  bscKey <- paste('DesBsc',mu, sep='')
  
  nxbKey <- paste('NXB',mu, sep='_')
  cpnKey <- paste('CPN',mu, sep='_')
  tarKey <- paste('TAR',mu, sep='_')
  secKey <- paste('SEC',mu, sep='_')
  bodKey <- paste('BOD',mu, sep='_')
  endKey <- paste('END',mu, sep='_')
  
  str <- glue::glue(
    paste0(rep(" ", 59+10), collapse=''),"{tib[[cpnKey]][1]}{tib[[secKey]][1]}{BNG}{tib[[bodKey]][1]}{tib[[endKey]][1]}{RET}",
    paste0(rep(" ", 59+8),  collapse=''),"{tib[[nxbKey]][1]}{BNG}{tib[[cpnKey]][1]}{tib[[secKey]][1]}{BNG}{tib[[bodKey]][1]}{RET}",
    "{fr1Key}_{co1Key}_{mu}{TAB}{tib[[bscKey]][1]}{RET}",
    "{fr2Key}_{co2Key}_{mu}{TAB}{Biostrings::reverse(tib[[bscKey]][2])}{RET}",
    paste0(rep(" ", 9+10), collapse=''),"{Biostrings::reverse(tib[[bodKey]][2])}{tib[[secKey]][2]}{BNG}{tib[[cpnKey]][2]}{tib[[nxbKey]][2]}{RET}",
    paste0(rep(" ", 9+9),  collapse=''),"{tib[[endKey]][2]}{Biostrings::reverse(tib[[bodKey]][2])}{tib[[secKey]][2]}{BNG}{tib[[cpnKey]][2]}{RET}",
    "{RET}")
  
  if ( vb >= vt+1 ) cat(str)
  
  str
}

# End of file

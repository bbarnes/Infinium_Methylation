
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       4.0 Analyze and Assign Cgn:: 
#                      CGN-Map/BSMAP/dbGCGN look-up
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   funcTag={funcTag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       4.0 Analyze and Assign Cgn:: 
#                      CGN-Map/BSMAP/dbGCGN look-up
#                           Main Workflow Driver
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_mapping_workflow = function(ord_tib,
                                bsp_tib,
                                seq_tib,
                                
                                ids_key = "Prb_Key",
                                des_key = "Prb_Des", 
                                din_key = "Prb_Din",
                                map_key = "Prb_Map",
                                
                                Cgn_Int = "Cgn_Int",
                                Can_Cgn = "Can_Cgn",
                                Ord_Cgn = "Ord_Cgn",
                                Bsp_Cgn = "Bsp_Cgn",
                                Imp_Cgn = "Imp_Cgn",

                                can_csv,
                                
                                join    = "inner",
                                merge   = FALSE,
                                retData = FALSE,
                                
                                out_csv = NULL,
                                out_dir,
                                run_tag, 
                                re_load = FALSE,
                                pre_tag = NULL,
                                end_str = 'csv.gz',
                                sep_chr = '.',
                                out_col = c("Prb_Key","Address","Ord_Des",
                                            "Ord_Din","Ord_Map","Ord_Prb"),
                                unq_col = c("Ord_Din","Ord_Map","Cgn_Int"),

                                verbose=0,vt=3,tc=1,tt=NULL,
                                funcTag='cgn_mapping_workflow') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  out_csv <- redata(out_dir, run_tag, funcTag, re_load, 
                    pre_tag, end_str=end_str, sep=sep_chr,
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  if (tibble::is_tibble(out_csv)) return(out_csv)
  if (is.null(out_csv)) {
    stop(glue::glue("{RET}{mssg} ERROR: out_csv is NULL!{RET2}"))
    return(out_csv) }
  out_dir <- base::dirname(out_csv)
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   funcTag={funcTag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  stime <- base::system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Consolidate/Assign CGNs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # TBD:: This function needs serious re-writing!!!!
    #
    ret_tib <- assign_cgn(ord_tib = ord_tib,
                          bsp_tib = bsp_tib,
                          seq_tib = seq_tib,
                          
                          ids_key = ids_key,
                          des_key = des_key,
                          din_key = din_key,
                          map_key = map_key,
                          
                          Cgn_Int = Cgn_Int,
                          Can_Cgn = Can_Cgn,
                          Ord_Cgn = Ord_Cgn,
                          Bsp_Cgn = Bsp_Cgn,
                          Imp_Cgn = Imp_Cgn,
                          
                          join    = join,
                          merge   = merge,
                          retData = retData,

                          can_csv = can_csv,
                          out_csv = out_csv,
                          out_col = out_col,
                          unq_col = unq_col,
                          
                          end_str = end_str,
                          sep_chr = sep_chr,
                          
                          verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)

    # Format Final Output::
    out_cols <- c("Prb_Key","Address","Ord_Des","Ord_Din","Ord_Map","Ord_Prb")
    ret_tib <- ret_tib %>% 
      dplyr::select(dplyr::all_of(out_cols), dplyr::everything())

    # Format Final Output::
    ret_tib <- ret_tib %>% 
      dplyr::select(dplyr::all_of(out_col), dplyr::everything())
    
    out_cnt <- safe_write(x=ret_tib, file=out_csv, funcTag=funcTag, done=TRUE,
                          verbose=verbose, vt=vt+2,tc=tc+1,tt=tt)
    
    tt$addFile(out_csv)
    
    if (retData) ret_dat$cgn_assign <- ret_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Assign Best CGN from:: BSP & SEQ
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This function needs to be cleaned up! Logic probably needs to be fixed
#   as well. 
#
assign_cgn = function(ord_tib,
                      bsp_tib,
                      seq_tib,

                      can_csv,
                      can_tib = NULL,
                      out_csv = NULL,
                      out_col = NULL,
                      unq_col = NULL,
                      
                      ids_key = "Prb_Key",
                      des_key = "Prb_Des", 
                      din_key = "Prb_Din",
                      map_key = "Prb_Map",
                      
                      Cgn_Int = "Cgn_Int",
                      Can_Cgn = "Can_Cgn",
                      Ord_Cgn = "Ord_Cgn",
                      Bsp_Cgn = "Bsp_Cgn",
                      Imp_Cgn = "Imp_Cgn",

                      join    = "inner",
                      merge   = FALSE,
                      retData = FALSE,
                      
                      end_str = 'csv.gz', 
                      sep_chr = '.',
                      
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='assign_cgn') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  sum_csv <- out_csv %>% 
    stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
    paste("cgn-counts-table.csv.gz", sep=sep_chr)
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}    can_csv={can_csv}.{RET}"))
    cat(glue::glue("{mssg}    sum_csv={sum_csv}.{RET}"))
    cat(glue::glue("{mssg}    out_csv={out_csv}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      ids_key={ids_key}.{RET}"))
    cat(glue::glue("{mssg}      join={join}.{RET}"))
    cat(glue::glue("{mssg}     merge={merge}.{RET}"))
    cat(glue::glue("{mssg}   retData={retData}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if (is.null(can_tib) && 
      !is.null(can_csv) &&
      !file.exists(can_csv)) {
    stop(glue::glue("{RET}{mssg} ERROR: File does NOT exist: ",
                    "canonical_csv={can_csv}!{RET2}"))
    return(ret_tib)
  }
  
  stime <- base::system.time({
    
    ids_sym <- rlang::sym(ids_key)
    cgn_sym <- rlang::sym(Cgn_Int)
    
    # Set up column variables::
    base_cols <- out_col
    if (is.null(out_col) || length(out_col)==0)
      base_cols <- c(ids_key,des_key,din_key,map_key)
    
    grps_cols <- c(base_cols, Cgn_Int)
    keys_cols <- c(base_cols, Cgn_Int)
    cnts_cols <- c("Ord_Cnt", "Can_Cnt", "Bsp_Cnt", "Seq_Cnt")
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Load Canonical CGNs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (is.null(can_tib) && 
        !is.null(can_csv) &&
        file.exists(can_csv)) {
      
      can_tib <- safe_read(file=can_csv, funcTag = funcTag,
                           verbose=verbose, vt=vt+1,tc=tc+1,tt=tt) %>% 
        dplyr::select(dplyr::all_of(Can_Cgn)) %>% 
        dplyr::rename(!!cgn_sym:=Can_Cgn) %>% 
        dplyr::mutate(Can_Cnt=1) %>%
        clean_tibble()
    }
    can_cnt <- base::nrow(can_tib)
    if (is.null(can_tib) || can_cnt==0) {
      stop(glue::glue("{RET}{mssg} ERROR: Canonical Tib is null or zero ",
                      "length; can_cnt={can_cnt}!{RET2}"))
      return(ret_tib)
    }
    can_key <- glue::glue("can_tib({funcTag})")
    can_cnt <- print_tib(can_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=can_key)

    if (retData) ret_dat$can_tib <- can_tib
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Defined Order tib to ord original cgn::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ord_tib <- ord_tib %>% 
      dplyr::select(dplyr::all_of(c(base_cols,Ord_Cgn)) ) %>%
      dplyr::rename(!!cgn_sym:=dplyr::all_of(Ord_Cgn)) %>% 
      dplyr::mutate(Ord_Cnt=1) %>%
      dplyr::distinct() %>%
      clean_tibble()
    
    ord_key <- glue::glue("ord_tib({funcTag})")
    ord_cnt <- print_tib(ord_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ord_key)

    if (retData) ret_dat$ord_tib <- ord_tib
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Format BSP::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    bsp_tib <- bsp_tib %>% 
      dplyr::filter(!is.na(Bsp_Cgn)) %>% 
      dplyr::select(dplyr::all_of(c(base_cols,Bsp_Cgn)) ) %>% 
      dplyr::rename(!!cgn_sym:=Bsp_Cgn) %>%
      dplyr::distinct() %>%
      # dplyr::arrange(ids_key, Cgn) %>%
      dplyr::group_by_at(grps_cols) %>%
      dplyr::summarise(Bsp_Cnt=n(), .groups = "drop") %>%
      clean_tibble()
    
    bsp_key <- glue::glue("bsp_tib({funcTag})")
    bsp_cnt <- print_tib(bsp_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=bsp_key)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Format Seq::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    seq_tib <- seq_tib %>% 
      dplyr::filter(!is.na(Imp_Cgn)) %>% 
      dplyr::select(dplyr::all_of(c(base_cols,Imp_Cgn)) ) %>%
      dplyr::rename(!!cgn_sym:=Imp_Cgn) %>%
      dplyr::distinct() %>%
      # dplyr::arrange(ids_key, Cgn) %>% 
      dplyr::group_by_at(grps_cols) %>% 
      dplyr::summarise(Seq_Cnt=n(), .groups = "drop") %>%
      clean_tibble()
    
    seq_key <- glue::glue("seq_tib({funcTag})")
    seq_cnt <- print_tib(seq_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=seq_key)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #             Critical Step:: Build and Sort Counts Tables
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # All in-one command::
    cnt_tib <- ord_tib %>%
      dplyr::left_join(can_tib, by=Cgn_Int) %>% dplyr::distinct() %>%
      dplyr::left_join(bsp_tib, by=c(keys_cols) ) %>% dplyr::distinct() %>%
      dplyr::left_join(seq_tib, by=c(keys_cols) ) %>% dplyr::distinct() %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(cnts_cols), tidyr::replace_na,0),
                    Sum_Cnt=Bsp_Cnt+Seq_Cnt,
                    Prd_Cnt=Bsp_Cnt*Seq_Cnt) %>% 
      dplyr::add_count(!!ids_sym, name="Cgn_Cnt") %>% 
      dplyr::arrange(-Can_Cnt,-Prd_Cnt,-Sum_Cnt,-Ord_Cnt) %>%
      dplyr::mutate(Rank=dplyr::row_number()) %>%
      clean_tibble()

    sum_cnt <- safe_write(cnt_tib, file = sum_csv, done = TRUE, 
                          verbose=verbose, vt=vt+2,tc=tc+1,tt=tt)
    
    cnt_key <- glue::glue("cnt_tib({funcTag})")
    cnt_cnt <- print_tib(cnt_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=cnt_key)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                Infinium II::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cnt_list <- cnt_tib %>% split(.$Ord_Des)
    
    inf2_tib <- cnt_list[["2"]] %>%
      # dplyr::arrange(!!ids_sym, Rank) %>%
      dplyr::arrange(Rank) %>%
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    
    inf2_key <- glue::glue("inf2_tib({funcTag})")
    inf2_cnt <- print_tib(inf2_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=inf2_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                              Infinium I:: Full Join
    #
    #  TBD:: The joining should really be done by sequence: Ord_Prb
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (join=="full") {
      inf1_tib <- dplyr::full_join(
        cnt_list[["U"]], 
        cnt_list[["M"]], 
        # by=c("Ord_Din","Ord_Map","Cgn"),
        by=c(unq_col),
        suffix=c("_U","_M")
      )
    } else if (join=="inner") {
      inf1_tib <- dplyr::inner_join(
        cnt_list[["U"]], 
        cnt_list[["M"]], 
        # by=c("Ord_Din","Ord_Map","Cgn"), 
        by=c(unq_col),
        suffix=c("_U","_M")
      )
    } else {
      stop(glue::glue("{mssg} Unsupported join type={join}.{RET}"))
      return(NULL)
    }
    inf1_tib <- inf1_tib %>%
      dplyr::mutate(Rank=pmin(Rank_U,Rank_M)) %>%
      # dplyr::arrange(Ord_Map, Rank) %>%
      dplyr::arrange(Rank) %>%
      dplyr::distinct(Ord_Map,Prb_Key_U,Prb_Key_M, .keep_all = TRUE)
    
    inf1_key <- glue::glue("inf1_tib({funcTag})")
    inf1_cnt <- print_tib(inf1_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=inf1_key)

    if (retData) ret_dat$inf1_tib <- inf1_tib
    if (retData) ret_dat$inf2_tib <- inf2_tib
    
    selU_cols <- c(out_col[!out_col %in% unq_col] %>% paste("U",sep="_"),
                   cnts_cols %>% paste("U",sep="_"),
                   unq_col, "Rank")
    selM_cols <- c(out_col[!out_col %in% unq_col] %>% paste("M",sep="_"),
                   cnts_cols %>% paste("M",sep="_"),
                   unq_col, "Rank")
    sel2_cols <- c(out_col[!out_col %in% unq_col], 
                   cnts_cols, 
                   unq_col, "Rank")
    name_cols <- c(out_col[!out_col %in% unq_col], cnts_cols, unq_col, "Rank")
    
    # Re-Bind all three: U/M/2::
    #
    ret_tib <- dplyr::bind_rows(
      dplyr::select(inf1_tib, dplyr::all_of(selU_cols) ) %>% 
        purrr::set_names(name_cols),
      
      dplyr::select(inf1_tib, dplyr::all_of(selM_cols) ) %>% 
        purrr::set_names(name_cols),
      
      dplyr::select(inf2_tib, dplyr::all_of(sel2_cols) ) %>%
        purrr::set_names(name_cols)
    ) %>% dplyr::distinct()

    ret_key <- glue::glue("ret_tib({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
    
    #
    # Need to fix these missing fields!!!!
    #
    mis_tib <- dplyr::anti_join(ord_tib, ret_tib, by=c(ids_key))
    
    sig_tib <- dplyr::filter(cnt_tib, ids_key %in% mis_tib[[ids_key]]) %>%
      dplyr::arrange(Ord_Map,Rank) %>%
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    mis_key <- glue::glue("mis_tib({funcTag})")
    mis_cnt <- print_tib(mis_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=mis_key)
    
    sig_key <- glue::glue("sig_tib({funcTag})")
    sig_cnt <- print_tib(sig_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=sig_key)
    
    mis_cnt <- ret_tib %>% dplyr::filter(is.na(!!ids_sym)) %>% base::nrow()
    mul_cnt <- ret_tib %>% dplyr::add_count(!!ids_sym,Cgn_Int, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()

    if (verbose>=vt) {
      cat(glue::glue("{mssg}   Miss Count={mis_cnt}.{RET}"))
      cat(glue::glue("{mssg}  Multi Count={mul_cnt}.{RET}"))
      cat(glue::glue("{mssg} Single Count={sig_cnt}.{RET2}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Merge all data together::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_cnt <- ret_tib %>% base::nrow()
    ret_tib <- dplyr::bind_rows(
      #
      # Add Formatted cg#'s
      #
      ret_tib %>% 
        dplyr::mutate(
          Cgn_Tag=dplyr::case_when(
            Ord_Din=="rs" ~ Ord_Din,
            Ord_Din=="ch" ~ Ord_Din,
            TRUE ~ "cg"
          ),
          Cgn_Str=dplyr::case_when(
            Ord_Din=="rs" ~ stringr::str_remove(Ord_Map, "[-_:].*$"),
            Ord_Din=="ch" ~ stringr::str_remove(Ord_Map, "[-_:].*$"),
            TRUE ~ paste0("cg",stringr::str_pad(Cgn_Int,width=8,side="left",pad="0"))
          )
        ),
      mis_tib %>% 
        # dplyr::select(!!ids_sym,Cgn,Ord_Des,Ord_Din,Ord_Map) %>% 
        dplyr::mutate(
          Can_Cnt=0, 
          Rank=dplyr::row_number() + ret_cnt,
          Cgn_Tag="uk",
          Cgn_Str=paste0(Cgn_Tag,stringr::str_pad(Cgn_Int,width=8,side="left",pad="0"))
        )
    ) %>%
      # TBD:: Capture other CGN's in separate column:: actual CGN's not Count!!
      dplyr::add_count(!!ids_sym, name="Alt_Cgn_Cnt") %>%
      # One Final Clean Up To Ensure Uniqueness::
      dplyr::arrange(Rank) %>% 
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    
    mul_cnt <- ret_tib %>% 
      dplyr::add_count(!!ids_sym,Cgn_Int, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
    
    if (verbose>=vt)
      cat(glue::glue("{mssg}  Multi Count Final={mul_cnt}.{RET2}"))
    
    if (mul_cnt!=0) {
      stop(glue::glue("{RET}{mssg} Multi-Count Final={mul_cnt} ",
                      "not equal to zero!!!{RET2}"))
      return(NULL)
    }
    
    if (merge) {
      ret_tib <- bsp_tib %>%
        dplyr::left_join(
          ret_tib, 
          by=c(ids_key,"Ord_Des","Ord_Din","Ord_Map"),
          suffix=c("_bsp","_cgn"))
      
    }
    ret_tib <- ret_tib %>% clean_tibble()
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file

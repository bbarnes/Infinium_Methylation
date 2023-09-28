
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#         3.0 Intersect Sequences Address and improbe:: U49/M49
#                         CGN Mapping Workflow()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

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
#                  Convert Annotation::Genomic Range Methods
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ann_to_grs = function( x,
                       vb=0,vt=3,tc=1,tt=NULL,
                       fun_tag='ann_to_grs' ) {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (vb>=vt) cat(glue::glue("[{fun_tag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (tibble::is_tibble(x)) {
      if (vb>=vt+4) {
        cat(glue::glue("[{fun_tag}]:{tabsStr} Usung tibble={RET}"))
        print(x)
      }
      tib <- x
    } else if (file.exists(x) && !dir.exists(x)) {
      if (vb>=vt+4)
        cat(glue::glue("[{fun_tag}]:{tabsStr} Loading from file={x}.{RET}"))
      
      tib <- suppressMessages(suppressWarnings(readr::read_tsv(x)))
    } else {
      stop(glue::glue("[{fun_tag}]:{tabsStr} Unknown type={RET}"))
      print(x)
      return(ret_tib)
    }
    ret_key <- glue::glue("tib-1")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    if (vb>=vt)
      cat(glue::glue("[{fun_tag}]:{tabsStr} Building GRanges...{RET}"))
    
    tib <- tib %>% dplyr::mutate(unq_key=paste(unq_key,rank,sep="."))
    
    # Double check unique key::
    tib_tot_cnt <- tib %>% base::nrow()
    tib_unq_cnt <- tib$unq_key %>% unique() %>% length()
    
    if (tib_tot_cnt != tib_unq_cnt) {
      cat(glue::glue("{RET}[{fun_tag}]:{tabsStr} ERROR: Total and Unique ",
                     "counts done match: tib_tot_cnt={tib_tot_cnt}, ",
                     "tib_unq_cnt={tib_unq_cnt}.{RET}"))
      tib <- tib %>% dplyr::add_count(unq_key, name="Unq_Key_Cnt")
      return(tib)
    }
    
    ret_tib =
      GenomicRanges::GRanges(
        seqnames=Rle(tib$chr),
        strand=Rle(tib$srd),
        
        name=tib$name,
        name2=tib$name2,
        class=tib$class,
        source=tib$source,
        tissue=tib$tissue,
        rank=tib$rank,
        evidence=tib$evidence,
        
        IRanges(start=tib$beg,
                end=tib$end,
                names=paste(tib$source,tib$unq_key, sep="_") )
      )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb>=vt)
    cat(glue::glue("[{fun_tag}]:{tabsStr} Done; ",
                   "Return Count={ret_cnt}; ",
                   "elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ",
                   "----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Quick Tib to GRS Converter::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tib_to_grs = function(tib, ids_key="Probe_ID", chr_key="Chromosome",
                      vb=0,vt=3,tc=1,tt=NULL,
                      fun_tag='tib_to_grs') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   fun_tag={fun_tag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ret_tib <- tib %>% dplyr::filter(! Prb_Des=="M")
    # ret_tib <- tib %>% dplyr::filter(! Prb_Des=="B")
    
    # Double check unique key::
    tib_tot_cnt <- ret_tib %>% base::nrow()
    tib_unq_cnt <- ret_tib %>% dplyr::pull(ids_key) %>% unique() %>% length()
    
    if (tib_tot_cnt != tib_unq_cnt) {
      cat(glue::glue("{RET}[{fun_tag}]:{tabs} ERROR: Total and Unique ",
                     "counts done match: tib_tot_cnt={tib_tot_cnt}, ",
                     "tib_unq_cnt={tib_unq_cnt}.{RET}"))
      ret_tib <- ret_tib %>% dplyr::add_count(unq_key, name="Unq_Key_Cnt")
      return(ret_tib)
    }
    
    man_grs <-
      GenomicRanges::GRanges(
        seqnames = Rle(ret_tib %>% dplyr::pull( chr_key ) ),
        # seqnames = Rle(ret_tib$Chromosome),
        # strand=Rle(ret_tib$Strand_FR),
        
        # Strand_FR=ret_tib$Strand_FR,
        Prb_Des=ret_tib$Prb_Des,
        # Prb_Din=ret_tib$Prb_Din,
        # Prb_Seq=ret_tib$Prb_Seq,
        Address=ret_tib$Address,
        # Prb_Mapq=ret_tib$Prb_Mapq,
        
        MASK_general=ret_tib$MASK_general,
        MASK_mapping=ret_tib$MASK_mapping,
        MASK_typeINextBaseSwitch=ret_tib$MASK_typeINextBaseSwitch,
        
        IRanges(start = ret_tib$Coordinate,
                width = 2,
                names=paste(ret_tib %>% dplyr::pull(ids_key)) )
      )
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  man_grs
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#         3.0 Intersect Sequences Address and improbe:: U49/M49
#                           Main Workflow Driver
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

seq_mapping_workflow = function( ord_tib,
                                 
                                 seq_dir, 
                                 pattern_u,
                                 pattern_m,
                                 
                                 prb_key = "Prb_Seq",
                                 add_key = "Address", 
                                 des_key = "Prb_Des", 
                                 din_key = "Prb_Din",
                                 ids_key = "Prb_Key",
                                 aln_key = "Aln_P49",
                                 
                                 prefix = NULL,
                                 suffix = NULL,
                                 
                                 idxA = 1,
                                 idxB = 1,
                                 
                                 out_dir,
                                 run_tag,
                                 pre_tag = NULL,
                                 
                                 reload     = 0,
                                 reload_min = 2,
                                 ret_data   = FALSE,
                                 parallel   = FALSE,
                                 ret_type   = "file",
                                 
                                 out_col = c("Prb_Key","Address","Ord_Des",
                                             "Ord_Din","Ord_Map","Ord_Prb"),
                                 
                                 vb=0, vt=3, tc=1, tt=NULL,
                                 fun_tag='seq_mapping_workflow') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+10,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+10,tc=tc+1,tt=tt ) )
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   pattern_u={pattern_u}.{RET}"))
    cat(glue::glue("{mssg}   pattern_m={pattern_m}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}     prb_key={prb_key}.{RET}"))
    cat(glue::glue("{mssg}     add_key={add_key}.{RET}"))
    cat(glue::glue("{mssg}     des_key={des_key}.{RET}"))
    cat(glue::glue("{mssg}     din_key={din_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}     out_dir={out_dir}.{RET}"))
    cat(glue::glue("{mssg}     out_csv={out_csv}.{RET}"))
    cat(glue::glue("{mssg}      prefix={prefix}.{RET}"))
    cat(glue::glue("{mssg}      suffix={suffix}.{RET}"))
    cat(glue::glue("{mssg}        idxA={idxA}.{RET}"))
    cat(glue::glue("{mssg}        idxB={idxB}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      reload={reload}.{RET}"))
    cat(glue::glue("{mssg}    parallel={parallel}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ids_sym <- rlang::sym(ids_key)
    prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    din_sym <- rlang::sym(din_key)
    aln_sym <- rlang::sym(aln_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Build Candidate Sub-String Probes:: U49/M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # NOTE: Build U49/M49 Sub-sequences and split by 5' two nucelotide prefix
    #  into sorted output files. Splitting by prefix makes the join later much
    #  faster...
    #
    
    ret_tib <- ord_tib %>%
      dplyr::mutate(
        Aln_Prb = deMs(!!prb_sym, uc=TRUE),
        !!aln_sym := dplyr::case_when(
          !!des_sym == '2' ~ stringr::str_sub(Aln_Prb, 2),
          !!des_sym == 'U' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          !!des_sym == 'M' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::distinct( !!ids_sym,Aln_Prb, .keep_all=TRUE ) %>%
      clean_tib()
    
    ret_key <- glue::glue("add-P49")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: U49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    u49_tib <-
      intersect_seq_workflow( tib = ret_tib,
                              
                              seq_dir = seq_dir,
                              tar_bsc = "U49",
                              tar_des = c("2","U"),
                              pattern = pattern_u,
                              
                              prefix  = prefix,
                              suffix  = suffix,
                              
                              prb_key = prb_key,
                              des_key = des_key,
                              din_key = din_key,
                              ids_key = ids_key,
                              aln_key = aln_key,
                              
                              idxA = idxA,
                              idxB = idxB,
                              
                              out_dir = out_dir,
                              run_tag = run_tag,
                              pre_tag = pre_tag,
                              
                              reload   = reload,
                              parallel = parallel,
                              
                              vb=vb, vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    m49_tib <-
      intersect_seq_workflow( tib = ret_tib,
                              
                              seq_dir = seq_dir,
                              tar_bsc = "M49",
                              tar_des = c("M"),
                              pattern = pattern_m,
                              
                              prefix  = prefix,
                              suffix  = suffix,
                              
                              prb_key = prb_key,
                              des_key = des_key,
                              din_key = din_key,
                              ids_key = ids_key,
                              aln_key = aln_key,
                              
                              idxA = idxA,
                              idxB = idxB,
                              
                              out_dir = out_dir,
                              run_tag = run_tag,
                              pre_tag = pre_tag,
                              
                              reload   = reload,
                              parallel = parallel,
                              
                              vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- join_seq_intersect( tib = ord_tib,
                                   u49 = u49_tib,
                                   m49 = m49_tib,
                                   
                                   ids_key = ids_key,
                                   out_col = out_col,
                                   
                                   vb=vb, vt=vt+1,tc=tc+1,tt=tt )
    
    # Format Final Output::
    if ( !is.null(out_col) )
    {
      ret_tib <- ret_tib %>% 
        dplyr::select( dplyr::all_of(out_col), dplyr::everything() )
      
      # Format Final Output::
      ret_tib <- ret_tib %>% 
        dplyr::select(dplyr::all_of(out_col), dplyr::everything())
    }
    
    out_cnt <- safe_write( ret_tib, file=out_csv, vb=vb, vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

intersect_seq_workflow = function(tib,
                                  
                                  seq_dir,
                                  pattern,
                                  
                                  tar_bsc,
                                  tar_des,
                                  
                                  prefix,
                                  suffix,
                                  
                                  prb_key = "Prb_Seq",
                                  des_key = "Prb_Des",
                                  din_key = "Prb_Din",
                                  ids_key = "Prb_Key",
                                  aln_key = "Aln_P49",
                                  
                                  idxA=1,idxB=1,
                                  
                                  out_dir,
                                  run_tag,
                                  pre_tag = NULL,
                                  
                                  reload     = 0,
                                  reload_min = 2,
                                  ret_data   = FALSE,
                                  parallel   = FALSE,
                                  
                                  vb=0,vt=3,tc=1,tt=NULL,
                                  fun_tag='intersect_seq_workflow') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+10,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+10,tc=tc+1,tt=tt ) )
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb>=vt+2) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   seq_dir={seq_dir}.{RET}"))
    cat(glue::glue("{mssg}   pattern={pattern}.{RET}"))
    cat(glue::glue("{mssg}   prefix={prefix}.{RET}"))
    cat(glue::glue("{mssg}   suffix={suffix}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}     tar_bsc={tar_bsc}.{RET}"))
    cat(glue::glue("{mssg}     tar_des={tar_des}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}     prb_key={prb_key}.{RET}"))
    cat(glue::glue("{mssg}     des_key={des_key}.{RET}"))
    cat(glue::glue("{mssg}     din_key={din_key}.{RET}"))
    cat(glue::glue("{mssg}     ids_key={ids_key}.{RET}"))
    cat(glue::glue("{mssg}     aln_key={aln_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}        idxA={idxA}.{RET}"))
    cat(glue::glue("{mssg}        idxB={idxB}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      reload={reload}.{RET}"))
    cat(glue::glue("{mssg}    parallel={parallel}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    din_sym <- rlang::sym(din_key)
    aln_sym <- rlang::sym(aln_key)
    
    tar_des_str <- paste(tar_des, collapse="-")
    
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #      Find All Pre-Built Reference Prefix-Partition Files:: U49/M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cgn_tsvs <- file_list( path = seq_dir,
                           pattern = pattern,
                           suffix  = pattern,
                           vb = vb)
    
    tsv_cnt <- cgn_tsvs %>% names() %>% length()
    pre_len <- cgn_tsvs %>% names() %>% stringr::str_length() %>% max()
    
    if (vb>=vt+6) {
      cat(glue::glue("{mssg} Found {tsv_cnt} {tar_bsc}-files, ",
                     "prefix length={pre_len}, tar_des={tar_des_str}.{RET}"))
      cgn_tsvs %>% head(n=3) %>% print()
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: U49/M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_tibs <- tib %>%
      dplyr::filter(!!des_sym %in% tar_des) %>%
      dplyr::filter(!is.na(!!aln_sym)) %>%
      dplyr::select(dplyr::all_of(c(aln_key,ids_key)) ) %>%
      dplyr::mutate(
        Pre_Nuc=stringr::str_sub(!!aln_sym, 1, pre_len)
      ) %>%
      dplyr::arrange(!!aln_sym) %>%
      split(f=.$Pre_Nuc, drop = TRUE)
    
    ret_key <- glue::glue("prb_tibs")
    ret_cnt <- print_tib( prb_tibs, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    if (vb>=vt+6) {
      nuc_cnt <- prb_tibs %>% names() %>% length()
      cat(glue::glue("{mssg} Found {nuc_cnt} {tar_bsc}-tibs, ",
                     "prefix length={pre_len}, tar_des={tar_des_str}.{RET}"))
      
      ret_key <- glue::glue("{tar_bsc}/{tar_des_str}-{tar_des}-pre1-int-tib")
      ret_cnt <- print_tib( prb_tibs[[2]], fun_tag = fun_tag, name = ret_key, 
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    }
    
    ret_tib <- NULL
    if (parallel) {
      if (vb>=vt)
        cat(glue::glue("{mssg}{TAB} Intersecting probe sequences ",
                       "{tar_bsc}-{tar_des_str} (Parallel)...{RET}"))
      
      ret_tib <- foreach (pre_nuc=names(prb_tibs), .combine=rbind) %dopar% {
        cur_tib <- prb_tibs[[pre_nuc]] %>% dplyr::select(-Pre_Nuc)
        
        intersect_seq_strand(
          can = cur_tib,
          ref = cgn_tsvs[[pre_nuc]],
          bsc_str = tar_bsc,
          pre_nuc = pre_nuc,
          idxA = idxA, idxB = idxB,
          
          out_dir = out_dir,
          run_tag = run_tag,
          pre_tag = pre_tag,
          
          prefix = prefix,
          suffix = suffix,
          reload = reload,
          
          vb=vb,vt=vt+1,tc=tc+1,tt=tt)
      }
    } else {
      if (vb>=vt)
        cat(glue::glue("{mssg}{TAB} Intersecting probe sequences ",
                       "{tar_bsc}-{tar_des_str} (Linear)...{RET}"))
      
      for (pre_nuc in names(prb_tibs)) {
        cur_tib <- prb_tibs[[pre_nuc]] %>% dplyr::select(-Pre_Nuc)
        
        prb_tib <- intersect_seq_strand(
          can = cur_tib,
          ref = cgn_tsvs[[pre_nuc]],
          bsc_str = tar_bsc,
          pre_nuc = pre_nuc,
          idxA = idxA, idxB = idxB,
          
          out_dir = out_dir, 
          run_tag = run_tag,
          pre_tag = pre_tag,
          
          prefix = prefix,
          suffix = suffix,
          reload = reload,
          
          vb=vb,vt=vt+1,tc=tc+1,tt=tt)
        ret_tib <- dplyr::bind_rows(ret_tib, prb_tib)
        
        if (vb>=vt)
          cat(glue::glue("{mssg}{TAB}{TAB} Done. Intersecting probe ",
                         "sequences {tar_bsc}/{tar_des_str} nuc={pre_nuc}{RET2}"))
      }
    }
    if (vb>=vt)
      cat(glue::glue("{mssg}{TAB} Done. Intersecting probe sequences ",
                     "{tar_bsc}/{tar_des_str}!{RET2}"))
    
    ret_key <- glue::glue("{tar_bsc}/{tar_des_str}-{tar_des}-int-tib({fun_tag})")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

intersect_seq_strand = function( can, 
                                 ref,
                                 
                                 bsc_str,
                                 pre_nuc,
                                 
                                 idxA, 
                                 idxB,
                                 
                                 prefix, 
                                 suffix,
                                 
                                 out_dir,
                                 run_tag,
                                 pre_tag = NULL,
                                 
                                 reload     = 0,
                                 reload_min = 2,
                                 ret_data   = FALSE,
                                 parallel   = FALSE,
                                 
                                 vb=0,vt=3,tc=1,tt=NULL,
                                 fun_tag='intersect_seq_strand') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+10,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+10,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}     fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          ref = '{ref}'.{RET}"))
    cat(glue::glue("{mssg}         idxA = '{idxA}'.{RET}"))
    cat(glue::glue("{mssg}         idxB = '{idxB}'.{RET}"))
    cat(glue::glue("{mssg}      bsc_str = '{bsc_str}'.{RET}"))
    cat(glue::glue("{mssg}      pre_nuc = '{pre_nuc}'.{RET}"))
    cat(glue::glue("{mssg}       prefix = '{prefix}'.{RET}"))
    cat(glue::glue("{mssg}       suffix = '{suffix}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
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
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }

  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    file_base = ""
    if ( !is.null(prefix) ) file_base <- paste(file_base, prefix, sep='.')
    if ( !is.null(suffix) ) file_base <- paste(file_base, suffix, sep='.')
    file_base <- paste(file_base, bsc_str, sep='.')
    
    ref_tsv <- ref
    can_tsv <- file.path( out_dir, paste( file_base, "tsv.gz", sep='.') )
    out_tsv <- file.path( out_dir, paste( file_base, "intersect.tsv.gz", sep='.') )
    safe_write( x = can, file = can_tsv, vb=vb, vt=vt+4,tc=tc+1,tt=tt )
    
    ret_tib <- intersect_seq( ref = ref_tsv,
                              can = can_tsv,
                              out_tsv = out_tsv,
                              
                              out_dir = out_dir,
                              run_tag = run_tag,
                              
                              idxA = idxA,
                              idxB = idxB,
                              reload = reload,
                              
                              vb=vb, vt=vt+1,tc=tc+1,tt=tt )
    
    ret_key <- glue::glue("{pre_nuc}-{bsc_str}-int-tib({fun_tag})")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

intersect_seq = function( ref, 
                          can,
                          
                          out_tsv,
                          
                          idxA = 1,
                          idxB = 1,
                          
                          out_dir,
                          run_tag,
                          pre_tag = NULL,
                          
                          reload     = 0,
                          reload_min = 2,
                          ret_data   = FALSE,
                          parallel   = FALSE,
                          
                          vb=0,vt=3,tc=1,tt=NULL,
                          fun_tag='intersect_seq') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  etime <- 0
  ftime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+10,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+10,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters:: {RET}"))
    cat(glue::glue("{mssg}     fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          ref = '{ref}'.{RET}"))
    cat(glue::glue("{mssg}          can = '{can}'.{RET}"))
    cat(glue::glue("{mssg}         idxA = '{idxA}'.{RET}"))
    cat(glue::glue("{mssg}         idxB = '{idxB}'.{RET}"))
    cat(glue::glue("{mssg}      out_tsv = '{out_tsv}'.{RET}"))
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
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  int_seq_cols <-
    cols(
      Imp_Seq  = col_character(),
      Imp_Nuc  = col_character(),
      
      Imp_SrdI = col_integer(),
      Imp_Srd3 = col_character(),
      
      Imp_Key  = col_character(),
      Imp_Scr  = col_character(),
      
      Imp_Cnt  = col_integer(),
      Prb_Key  = col_character()
    )
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (reload && file.exists(out_tsv)) {
      if ( vb >= vt+1) cat(glue::glue("{mssg} Reloading!{RET}"))
    } else {
      clean <- FALSE
      if (stringr::str_ends(can, '.gz')) {
        cmd_str <- glue::glue("gzip -f -k -d {can}")
        if ( vb >= vt )
          cat(glue::glue("[{fun_tag}]: Running cmd={cmd_str}...{RET}"))
        
        cmd_ret <- base::system(cmd_str)
        if (cmd_ret!=0) {
          stop(glue::glue("{RET}[{fun_tag}]: ERROR: Failed(cmd_ret={cmd_ret}) ",
                          "cmd={cmd_str}!{RET2}"))
          return(ret_tib)
        }
        can <- stringr::str_remove(can, ".gz$")
        clean <- TRUE
      }
      
      cmd_str = glue::glue("gzip -dc {ref} | join -t $'\t' -1{idxA} ",
                           "-2{idxB} - {can} | gzip -c - > {out_tsv}")
      if ( vb >= vt )
        cat(glue::glue("[{fun_tag}]: Running cmd='{cmd_str}'...{RET}"))
      cmd_ret <- system(cmd_str)
      
      if (clean && !stringr::str_ends(can, '.gz')) {
        cmd_str <- glue::glue("rm {can}")
        if ( vb >= vt )
          cat(glue::glue("[{fun_tag}]: Running cmd='{cmd_str}'...{RET}"))
        
        cmd_ret <- system(cmd_str)
        if (cmd_ret!=0) {
          stop(glue::glue("{RET}[{fun_tag}]: ERROR: Failed(cmd_ret='{cmd_ret}') ",
                          "cmd={cmd_str}!{RET2}"))
          return(ret_tib)
        }
        can <- paste(can,'gz', sep='.')
      }
    }
    
    if ( vb >= vt )
      cat(glue::glue("[{fun_tag}]: Loading intersection out_tsv={out_tsv}...{RET}"))
    
    ret_tib <- suppressMessages(suppressWarnings(
      readr::read_tsv(out_tsv, col_names=names(int_seq_cols$cols),
                      col_types=int_seq_cols) )) # %>% clean_tib()
    
    ret_key <- glue::glue("ret-fin({fun_tag})")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

join_seq_intersect = function(tib, u49, m49,
                              
                              ids_key = "Prb_Key",
                              out_col,
                              
                              vb=0,vt=3,tc=1,tt=NULL,
                              fun_tag='join_seq_intersect') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ret_tib <- tib %>% 
      dplyr::select(dplyr::all_of(out_col)) %>%
      dplyr::right_join(
        dplyr::bind_rows(u49,m49),
        by=ids_key,
      ) %>%
      dplyr::select(-Imp_SrdI,-Imp_Scr) %>% 
      dplyr::rename(Aln_Prb=Imp_Seq, Aln_Nuc=Imp_Nuc) %>%
      dplyr::mutate(Imp_Key=stringr::str_split(Imp_Key, pattern=",") ) %>% 
      tidyr::unnest(Imp_Key) %>%
      tidyr::separate(Imp_Key, 
                      into=c("Imp_Cgn",
                             "Imp_Hit_hg38", "Imp_Hit_hg37", 
                             "Imp_Hit_hg36", "Imp_Hit_mm10"), 
                      sep="_", remove=TRUE) %>%
      tidyr::separate(Imp_Srd3, into=c("Imp_TB","Imp_CO", "Imp_Nxb"),
                      sep=c(1,2)) %>%
      dplyr::select(dplyr::all_of(out_col), Imp_Cgn, Imp_TB, Imp_CO, Imp_Nxb, 
                    Aln_Prb, Aln_Nuc, dplyr::everything()) %>%
      clean_tib()
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file


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
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Infinium Methylation Probe Design Methods::
#               s-improbe = improbe by BSC genome sub-string
#                    Allows All Probe Designs:: cg,ch,rs 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

s_improbe_workflow = function(prb_tib, 
                              gen_tib,

                              # srd_key,  # F/R -- Not sure if this is needed
                              pos_key,  # "Bsp_Pos/Coordinate",
                              chr_key,  # "Bsp_Chr/Chromosome",
                              
                              nrec = 0,
                              
                              # Probe Info::
                              ids_vec = NULL,  # c(ids_key="Aln_Key_Unq", des_key="Ord_Des", din_key="Ord_Din")
                              prb_srd = NULL,  # "Strand_FR/Bsp_FR",
                              prb_cos = NULL,  # "Strand_CO/Bsp_CO",
                              prb_tbs = NULL,  # "Strand_TB"
                              
                              query_n= NULL, # IUPAC insertion nucleotides
                              
                              tri_din = NULL, # "rs" tri-fecta search type
                              
                              build=c("Prb_1C","Prb_2C","Prb_1O","Prb_2O"),
                              
                              
                              # imp_seq = "Template_Sequence",
                              # ext_seq = "Ext_Template_Seq",
                              # iup_seq = "Iupac_Forward_Sequence",
                              # iupac = NULL,
                              
                              # ref_col = "Ref",
                              # alt_col = "Alt",
                              # iup_col = "Iupac",
                              
                              
                              ups_len = 60, 
                              seq_len = 122, 

                              reload  = FALSE,
                              retData = FALSE,
                              
                              parallel   = FALSE,
                              add_flanks = FALSE, 
                              add_probes = FALSE,
                              
                              out_csv = NULL,
                              out_dir,
                              run_tag, 
                              re_load = FALSE,
                              pre_tag = NULL,
                              end_str = 'csv.gz',
                              sep_chr = '.',
                              
                              verbose=0,vt=4,tc=1,tt=NULL,
                              funcTag='s_improbe_workflow') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  out_csv <- redata(out_dir, run_tag, funcTag, re_load, 
                    pre_tag, end_str=end_str, sep=sep_chr,
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  if (tibble::is_tibble(out_csv)) return(out_csv)
  if (is.null(out_csv)) {
    stop(glue::glue("{RET}{mssg} ERROR: out_csv is NULL!{RET2}"))
    return(out_csv)
  }
  out_dir <- base::dirname(out_csv)
  
  sum_csv <- out_csv %>% 
    stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
    paste("summary.csv.gz", sep=sep_chr)
  
  prb_cnt <- prb_tib %>% base::nrow()
  gen_cnt <- gen_tib %>% base::nrow()

  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Genome Parameters::{RET}"))
    cat(glue::glue("{mssg}          nrec={nrec}.{RET}"))
    cat(glue::glue("{mssg}       prb_cnt={prb_cnt}.{RET}"))
    cat(glue::glue("{mssg}       gen_cnt={gen_cnt}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Output File Parameters::{RET}"))
    cat(glue::glue("{mssg}       out_csv={out_csv}{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Field Parameters::{RET}"))
    # cat(glue::glue("{mssg}       ids_key={ids_key}{RET}"))
    # cat(glue::glue("{mssg}       des_key={des_key}{RET}"))
    # cat(glue::glue("{mssg}       din_key={din_key}{RET}"))
    # cat(glue::glue("{mssg}       srd_key={srd_key}{RET}"))
    # cat(glue::glue("{mssg}       cos_key={cos_key}{RET}"))
    cat(glue::glue("{mssg}       tar_din={tar_din}{RET}"))
    cat(glue::glue("{RET}"))
    # cat(glue::glue("{mssg}       ext_seq={ext_seq}.{RET}"))
    # cat(glue::glue("{mssg}       iup_seq={iup_seq}.{RET}"))
    # cat(glue::glue("{mssg}       imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
    # cat(glue::glue("{mssg}       srd_str={srd_str}.{RET}"))
    cat(glue::glue("{mssg}       pos_key={pos_key}.{RET}"))
    cat(glue::glue("{mssg}       chr_key={chr_key}.{RET}"))
    cat(glue::glue("{RET}"))
    # cat(glue::glue("{mssg}       ref_col={ref_col}.{RET}"))
    # cat(glue::glue("{mssg}       alt_col={alt_col}.{RET}"))
    # cat(glue::glue("{mssg}       iup_col={iup_col}.{RET}"))
    cat(glue::glue("{mssg}         query_n={query_n}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    # cat(glue::glue("{mssg}           del={del}.{RET}"))
    cat(glue::glue("{mssg}       ups_len={ups_len}.{RET}"))
    cat(glue::glue("{mssg}       seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
    # cat(glue::glue("{mssg}      subset={subset}.{RET}"))
    # cat(glue::glue("{mssg}    sub_cols={sub_cols}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      reload={reload}.{RET}"))
    cat(glue::glue("{mssg}     retData={retData}.{RET}"))
    cat(glue::glue("{mssg}    parallel={parallel}.{RET}"))
    cat(glue::glue("{mssg}  add_flanks={add_flanks}.{RET}"))
    cat(glue::glue("{mssg}  add_probes={add_probes}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  etime   <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  stime <- base::system.time({
    if (verbose>=vt) 
      cat(glue::glue("{mssg} Building fresh...{RET}"))
    
    # Define symbolic variables::
    #
    # ids_sym  <- rlang::sym(ids_key)
    # ext_sym  <- rlang::sym(ext_seq)
    # iup_sym  <- rlang::sym(iup_seq)
    # imp_sym  <- rlang::sym(imp_seq)
    
    pos_sym  <- rlang::sym(pos_key)
    chr_sym  <- rlang::sym(chr_key)
    
    # ref_col_sym  <- rlang::sym(ref_col)
    # alt_col_sym  <- rlang::sym(alt_col)
    # iup_col_sym  <- rlang::sym(iup_col)
    
    # Load Genome::
    #
    gen_chr_dat <- load_genome(file=gen_fas,
                           nrec=nrec,
                           chr_key=chr_key,
                           ret_map=TRUE,
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (retData) ret_dat$gen <- seq_dat
    
    # Split Data by Chromosome::
    prb_chroms_list <- prb_tib %>% 
      dplyr::arrange(!!chr_sym,!!pos_sym) %>%
      split(.[[chr_key]])
    
    gen_chroms_list <- gen_chr_dat$maps %>%
      split(.[[chr_key]])
    
    chrom_names <- intersect(names(prb_chroms_list), names(gen_chroms_list))
    
    if (FALSE && parallel) {
      if (verbose>=vt) 
        cat(glue::glue("{mssg}{TAB} Extacting sequence templates ",
                       "from genome (Parallel)...{RET}"))
      
      ret_tib <- foreach (chrom=chrom_names, .combine=rbind) %dopar% {
        chr_idx <- chr_maps[[chr_str]] %>% head(n=1) %>% pull(Idx) %>% as.integer()
        s_improbe_template_workflow(
          tib=chr_list[[chr_str]], seq=seq_dat$seqs[[chr_idx]],
          srd_str=srd_str, pos_key=pos_key, chr_key=chr_key, chr_str=chr_str,
          ext_seq=ext_seq, iup_seq=iup_seq, imp_seq=imp_seq,
          ups_len=ups_len, seq_len=seq_len, iupac=iupac, del=del,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      }
      
    } else {
      if (verbose>=vt) 
        cat(glue::glue("{mssg}{TAB} Extacting sequence templates ",
                       "from genome (Linear)...{RET}"))
      
      for (chrom in chrom_names) {
        chr_idx <- chr_maps[[chr_str]] %>% head(n=1) %>% pull(Idx) %>% as.integer()
        
        
        cur_tib <- s_improbe_template_workflow(
          tib=chr_list[[chr_str]], seq=seq_dat$seqs[[chr_idx]], 
          srd_str=srd_str, pos_key=pos_key, chr_key=chr_key, chr_str=chr_str,
          ext_seq=ext_seq, iup_seq=iup_seq, imp_seq=imp_seq,
          ups_len=ups_len, seq_len=seq_len, iupac=iupac, del=del,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                          Merge Probe Designs::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
        cur_key <- glue::glue("cur-tib({funcTag}-{chr_str})")
        cur_cnt <- print_tib(cur_tib %>% dplyr::select(!!ids_sym,!!imp_sym),
                             funcTag, verbose,vt=vt+4,tc=tc+1, n=cur_key)
        
        if (verbose>=vt+1)
          cat(glue::glue("{mssg}{TAB} Done. Substring ",
                         "chr_str={chr_str}.{RET2}"))
      }
    }
    cat(glue::glue("{mssg}{TAB} Done. Extacting sequence templates ",
                   "from genome.{RET2}"))
    
    if (retData) ret_dat$seq <- ret_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc, n=ret_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #      Forward Template Sequence Generation:: Tri-fecta s-improbe
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (add_flanks) {
      
      tri_tib <- NULL
      tri_tib <- s_improbe_trifecta(tib=ret_tib,
                                    
                                    tar_din=tar_din, 
                                    ids_key=ids_key, 
                                    din_key=din_key,
                                    
                                    pos_key=pos_key,
                                    
                                    ext_seq=ext_seq, 
                                    iup_seq=iup_seq, 
                                    imp_seq=imp_seq, 
                                    
                                    ref_col=ref_col,
                                    alt_col=alt_col,
                                    iup_col=iup_col,
                                    
                                    verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      
      if (retData) ret_dat$tri <- tri_tib
      ret_tib <- dplyr::bind_rows(ret_tib, tri_tib)
    }
    
    if (add_probes) {
      
      prb_tib <- NULL
      prb_tib <- s_improbe(tib = ret_tib,
                           ids_key = ids_key,
                           
                           des_key = des_key,
                           din_key = din_key,
                           srd_key = srd_key,
                           cos_key = cos_key,
                           
                           pos_key = pos_key,
                           chr_key = chr_key,
                           build   = build, 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      if (retData) ret_dat$prb <- prb_tib
      
      prb_cols <- c(ids_key,din_key,pos_key,chr_key)
      ret_tib <- ret_tib %>% 
        dplyr::inner_join(prb_tib,
                          by=c(prb_cols),
                          suffix=c("","_prb"))
    }
    
    ret_tib <- ret_tib %>% dplyr::arrange(!!ids_sym)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Rename Strands To Match Standards::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- ret_tib %>% 
      dplyr::mutate(improbe_type = "s") %>%
      dplyr::rename(Strand_FR = dplyr::all_of(srd_key),
                    Strand_CO = dplyr::all_of(cos_key) )
    
    if (!is.null(gen_tib)) ret_tib <- cbind(gen_tib, ret_tib) %>% 
      tibble::as_tibble()

    ret_key <- glue::glue("After-gen-tib-cbind({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Calculate Data Summary:: CSV
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt)
      cat(glue::glue("{mssg}{TAB} Calculating Summary...{RET}"))
    
    sum_tib <- ret_tib %>% 
      dplyr::group_by(up61,dn61,Des_Din) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    sum_key <- glue::glue("ret-summary({funcTag})")
    sum_cnt <- print_tib(sum_tib,funcTag, verbose,vt+4,tc, n=sum_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Write Full Data Set:: CSV
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt+1)
      cat(glue::glue("{mssg} Writing seq table CSV={out_csv}...{RET}"))
    
    sum_cnt <- safe_write(x=sum_tib, file=sum_csv, funcTag=funcTag, 
                          verbose=verbose,vt=vt,tc=tc,append=FALSE)
    seq_cnt <- safe_write(x=ret_tib, file=out_csv, funcTag=funcTag, done=TRUE,
                          verbose=verbose,vt=vt,tc=tc,append=FALSE)
    
    # tt$addFile(out_csv)
    
    if (retData) ret_dat$sum <- sum_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

s_improbe_workflow_old = function(tib,
                              
                              build=c("Prb_1C","Prb_2C","Prb_1O","Prb_2O"),
                              
                              out_csv=NULL, out_dir, run_tag, 
                              re_load=FALSE, pre_tag=NULL,
                              end_str='csv.gz', sep_chr='.',

                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='s_improbe_workflow_old') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  out_csv <- redata(out_dir, run_tag, funcTag, re_load, 
                    pre_tag, end_str=end_str, sep=sep_chr,
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  if (tibble::is_tibble(out_csv)) return(out_csv)
  if (is.null(out_csv)) {
    stop(glue::glue("{RET}{mssg} ERROR: out_csv is NULL!{RET2}"))
    return(out_csv)
  }
  out_dir <- base::dirname(out_csv)
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_dir={out_dir}.{RET}"))
    cat(glue::glue("{mssg}      out_csv={out_csv}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       build={build}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ret_tib <- s_improbe(tib = tib, build = build, 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    out_cnt <- safe_write(x=ret_tib, file=out_csv, funcTag=funcTag, done=TRUE,
                          verbose=verbose, vt=vt,tc=tc,append=FALSE)

    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Infinium Methylation Probe Design Methods::
#               s-improbe = improbe by BSC genome sub-string
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

s_improbe = function(tib, 
                     
                     ids_key,
                     
                     des_key,
                     din_key,
                     srd_key,
                     cos_key,
                     
                     pos_key,
                     chr_key,
                     
                     build=c("Prb1C","Prb2C","Prb1O","Prb2O"),
                     
                     verbose=0,vt=3,tc=1,tt=NULL,
                     funcTag='s_improbe') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2-2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Field Parameters::{RET}"))
    cat(glue::glue("{mssg}       ids_key={ids_key}{RET}"))
    cat(glue::glue("{mssg}       des_key={des_key}{RET}"))
    cat(glue::glue("{mssg}       din_key={din_key}{RET}"))
    cat(glue::glue("{mssg}       srd_key={srd_key}{RET}"))
    cat(glue::glue("{mssg}       cos_key={cos_key}{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}       pos_key={pos_key}{RET}"))
    cat(glue::glue("{mssg}       chr_key={chr_key}{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}         build={build}{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    #
    # Build individual parts of the template sequence:: 
    #
    #                                            iupac
    #     up01.up02...up11.up12.up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn12.dn11.dn02.dn01
    #                                       Nxb [  C   G  ] Nxb
    #
    # Probe Design Formulas::
    #
    # Inf1C                               Nxb60* up61------------------dn58
    # Inf2C                                     ext61* dn61-------------------dn12
    #
    # Inf1O                up12------------------up61 Nxb61*
    # Inf2O           up11------------------up60 ext61*
    #
    #
    #  TBD:: Add reverse/complements::
    #
    # Build requested sub-string probes::
    #
    sel_cols <- c(ids_key,din_key,pos_key,chr_key,
                  "Prb_1C","Nxb_1C","Len_1C",
                  "Prb_2C","Nxb_2C","Len_2C",
                  "Prb_1O","Nxb_1O","Len_1O",
                  "Prb_2O","Nxb_2O","Len_2O",
                  "Tmp_Seq","Tmp_Len","Tmp_Pad")
    
    ret_tib <- tib
    if ("Prb_1C" %in% build)
      ret_tib <- ret_tib %>% dplyr::mutate(
        Prb_1C=paste0(iupac,dn61,dn60,dn59,dn58) %>% revCmp(),
        Nxb_1C=paste0(up60) %>% cmpl(),
        Len_1C=stringr::str_length(Prb_1C) )
    
    if ("Prb_2C" %in% build)
      ret_tib <- ret_tib %>% dplyr::mutate(
        Prb_2C=paste0(dn61,dn60,dn59,dn58,dn12) %>% revCmp(),
        Nxb_2C=paste0(iupac) %>% cmpl(),
        Len_2C=stringr::str_length(Prb_2C))
    
    if ("Prb_1O" %in% build)
      ret_tib <- ret_tib %>% dplyr::mutate(
        Prb_1O=paste0(up12,up58,up59,up60,iupac) %>% revCmp(),
        Nxb_1O=paste0(dn61) %>% cmpl(),
        Len_1O=stringr::str_length(Prb_1O) )
    
    if ("Prb_2O" %in% build)
      ret_tib <- ret_tib %>% dplyr::mutate(
        Prb_2O=paste0(up11,up12,up58,up59,up60) %>% revCmp(),
        Nxb_2O=paste0(iupac) %>% cmpl(),
        Len_2O=stringr::str_length(Prb_2O))

      # NOTE:: I think you can ignore FR strand since the template sequence
      #   should already be 5' -> 3' for the strand of interest
      #
      # if ("Prb2_RC" %in% build)
      #   ret_tib <- ret_tib %>% dplyr::mutate(
      #     Prb2_RC=paste0(up12,up58,up59,up60,up61), 
      #     Prb2_RC_Len=stringr::str_length(Prb2_RC))
      # 
      # if ("Prb1_RC" %in% build)
      #   ret_tib <- ret_tib %>% dplyr::mutate(
      #     Prb1_RC=paste0(up58,up59,up60,up61,iupac), 
      #     Prb1_RC_Len=stringr::str_length(Prb1_RC))

    # NOTE::This is just for graphical sanity checks and should be removed
    #   once its validated...
    #
    ret_tib <- ret_tib %>% dplyr::mutate(
      Tmp_Seq=paste(up01,up02,up10,up11,up12,up58,up59,up60,up61,
                     dn61,dn60,dn59,dn58,dn12,dn11,dn10,dn02,dn01, sep=""),
      Tmp_Len=stringr::str_length(Tmp_Seq) - 2,
      Tmp_Rvc=revCmp(Tmp_Seq) %>% addBrac(),
      Tmp_Seq=addBrac(Tmp_Seq),
      Tmp_Pad=paste(up01,up02,up10,up11,up12,up58,up59,up60,"[",up61,
                     dn61,"]",dn60,dn59,dn58,dn12,dn11,dn10,dn02,dn01, sep=" ") )
    
    ret_tib <- ret_tib %>% dplyr::select(dplyr::all_of(sel_cols))
    
    ret_key <- glue::glue("s-improbe-tib({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
}

s_improbe_template_workflow = function(tib, 
                                       seq,
                                       
                                       srd_str="F",
                                       pos_key="Coordinate",
                                       chr_key="Chromosome",
                                       chr_str,
                                       
                                       ext_seq="Ext_Forward_Seq",
                                       iup_seq="Iupac_Forward_Sequence",
                                       imp_seq="Forward_Sequence",
                                       
                                       ups_len=60, 
                                       seq_len=122, 
                                       iupac=NULL,
                                       del="_",
                                       
                                       verbose=0,vt=3,tc=1,tt=NULL,
                                       funcTag='s_improbe_template_workflow') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Field Parameters::{RET}"))
    cat(glue::glue("{mssg}   srd_str={srd_str}.{RET}"))
    cat(glue::glue("{mssg}   pos_key={pos_key}.{RET}"))
    cat(glue::glue("{mssg}   chr_key={chr_key}.{RET}"))
    cat(glue::glue("{mssg}   chr_str={chr_str}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}   ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("{mssg}   iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("{mssg}   imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}     iupac={iupac}.{RET}"))
    cat(glue::glue("{mssg}   ups_len={ups_len}.{RET}"))
    cat(glue::glue("{mssg}   seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # map_key <- glue::glue("map-tib({funcTag})")
    # map_cnt <- print_tib(map,funcTag, verbose,vt=0,tc, n=map_key)
    
    
    chr_sym <- rlang::sym(chr_key)
    cur_tib <- tib
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Parse Forward Template Sequence from Genome::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ref_seqs <- parse_genomic_seqs(tib=cur_tib,
                                   seq=as.character(seq),
                                   
                                   srd_str=srd_str,
                                   pos_key=pos_key,
                                   chr_str=chr_str,
                                   
                                   ups_len=ups_len, 
                                   seq_len=seq_len,
                                   verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                Forward Template Sequence Generation:: s-improbe
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- s_improbe_template(tib=cur_tib,
                                  seqs=ref_seqs,
                                  
                                  chr_str=chr_str, 
                                  
                                  ext_seq=ext_seq,
                                  iup_seq=iup_seq,
                                  imp_seq=imp_seq,
                                  
                                  iupac=iupac, 
                                  ups_len=ups_len,
                                  seq_len=seq_len,
                                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
}

# s_improbe_parse(tib=cur_tib, seqs=ref_seqs, chr_str=chr_str, iupac=iupac, ext_seq=ext_seq, iup_seq=iup_seq, imp_seq=imp_seq)
s_improbe_template = function(tib, seqs, 
                              chr_str="chr0", 
                              
                              ext_seq="Ext_Forward_Seq",
                              iup_seq="Iupac_Forward_Sequence",
                              imp_seq="Forward_Sequence",
                              
                              iupac=NULL,
                              ups_len=60,
                              seq_len=122,
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='s_improbe_template') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Field Parameters::{RET}"))
    cat(glue::glue("{mssg}   chr_str={chr_str}{RET}"))
    cat(glue::glue("{mssg}   ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("{mssg}   iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("{mssg}   imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}       iupac={iupac}.{RET}"))
    cat(glue::glue("{mssg}     ups_len={ups_len}.{RET}"))
    cat(glue::glue("{mssg}     seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ext_sym  <- rlang::sym(ext_seq)
    iup_sym  <- rlang::sym(iup_seq)
    imp_sym  <- rlang::sym(imp_seq)
    
    #
    # Build individual parts of the template sequence:: 
    #
    #                                            iupac
    #     up01.up02...up11.up12.up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn12.dn11.dn02.dn01
    #                                       Nxb [  C   G  ] Nxb
    #
    
    ref_up01 <- stringr::str_sub(seqs, 1,1)
    ref_up02 <- stringr::str_sub(seqs, 2,2)
    
    ref_up10 <- stringr::str_sub(seqs,  3,12)
    ref_up11 <- stringr::str_sub(seqs, 13,13)
    ref_up12 <- stringr::str_sub(seqs, 14,14)
    
    ref_up58 <- stringr::str_sub(seqs, 15,ups_len-2+2)
    
    ref_up59 <- stringr::str_sub(seqs, ups_len-2+3,ups_len-2+3)
    ref_up60 <- stringr::str_sub(seqs, ups_len-2+4,ups_len-2+4)
    ref_up61 <- stringr::str_sub(seqs, ups_len-2+5,ups_len-2+5)
    
    iupac_vec <- ref_up61
    if (!is.null(iupac)) iupac_vec <- tib %>% dplyr::pull(!!iupac)
    
    ref_dn61 <- stringr::str_sub(seqs, ups_len-2+6,ups_len-2+6)
    ref_dn60 <- stringr::str_sub(seqs, ups_len-2+7,ups_len-2+7)
    ref_dn59 <- stringr::str_sub(seqs, ups_len-2+8,ups_len-2+8)
    
    ref_dn58 <- stringr::str_sub(seqs, ups_len-2+9,ups_len-2+ups_len-2+8-12)
    
    ref_dn12 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+8-12+1,ups_len-2+ups_len-2+8-12+1)
    # ref_dn11 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+8-12+2,ups_len-2+ups_len-2+8)
    ref_dn11 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+8-12+2,ups_len-2+ups_len-2+8-12+2)
    ref_dn10 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+8-12+3,ups_len-2+ups_len-2+8)
    
    ref_dn02 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
    ref_dn01 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
    
    if (verbose>=vt+4) {
      cat(glue::glue("{mssg}{TAB} seqs({chr_str})={RET}"))
      seqs %>% head(n=2) %>% print()
    }
    
    # Add all fields to current tibble::
    #
    ret_tib <- tib %>%
      dplyr::mutate(
        !!ext_sym := seqs %>% addBrac(),
        
        up01=ref_up01,
        up02=ref_up02,
        
        up10=ref_up10,
        up11=ref_up11,
        up12=ref_up12,
        
        up58=ref_up58,
        
        up59=ref_up59,
        up60=ref_up60,
        up61=ref_up61,
        
        iupac=iupac_vec,
        
        dn61=ref_dn61,
        dn60=ref_dn60,
        dn59=ref_dn59,
        
        dn58=ref_dn58,
        
        dn12=ref_dn12,
        dn11=ref_dn11,
        dn10=ref_dn10,
        
        dn02=ref_dn02,
        dn01=ref_dn01,
        
        # Now we can assemble to optimal template::
        #
        !!iup_sym := paste0(up10,up11,up12, up58, up59,up60, iupac,dn61, dn60,dn59, dn58, dn12,dn11,dn10) %>% addBrac(),
        !!imp_sym := paste0(up10,up11,up12, up58, up59,up60,     "CG",   dn60,dn59, dn58, dn12,dn11,dn10) %>% addBrac(),
        
        # TBD:: Pretty sure this right, but really needs more testing...
        #
        Din_Str=stringr::str_to_lower(paste0(iupac,dn61)),
        Des_Din=dplyr::case_when(
          up61=='C' & dn61=='G' ~ "cg",
          up61=='C' & dn61!='G' ~ "ch",
          up61!='C' ~ "rs",
          TRUE ~ NA_character_
        )
      )
    ret_key <- glue::glue("ret_tib_2=chr_str={chr_str}")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
}

s_improbe_trifecta = function(tib, 
                              tar_din="rs",
                              ids_key="Aln_Key_Unq",
                              din_key="Ord_Din",
                              
                              pos_key="Coordinate",
                              chr_str="chr0", 
                              
                              ext_seq="Ext_Forward_Seq",
                              iup_seq="Iupac_Forward_Sequence",
                              imp_seq="Forward_Sequence",
                              
                              ref_col="Ref",
                              alt_col="Alt",
                              iup_col="Iupac",
                              
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='s_improbe_trifecta') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Field Parameters::{RET}"))
    cat(glue::glue("{mssg}   tar_din={tar_din}{RET}"))
    cat(glue::glue("{mssg}   ids_key={ids_key}{RET}"))
    cat(glue::glue("{mssg}   din_key={din_key}{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}   pos_key={pos_key}.{RET}"))
    cat(glue::glue("{mssg}   chr_str={chr_str}{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}   ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("{mssg}   iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("{mssg}   imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}   ref_col={ref_col}.{RET}"))
    cat(glue::glue("{mssg}   alt_col={alt_col}.{RET}"))
    cat(glue::glue("{mssg}   iup_col={iup_col}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ids_sym  <- rlang::sym(ids_key)
    din_sym  <- rlang::sym(din_key)
    ext_sym  <- rlang::sym(ext_seq)
    iup_sym  <- rlang::sym(iup_seq)
    imp_sym  <- rlang::sym(imp_seq)
    
    ref_col_sym  <- rlang::sym(ref_col)
    alt_col_sym  <- rlang::sym(alt_col)
    iup_col_sym  <- rlang::sym(iup_col)
    
    pos_sym  <- rlang::sym(pos_key)
    
    if (verbose>=vt+1)
      cat(glue::glue("{mssg}{TAB} Adding Upstream CG ",
                     "Flank Seuqnces({chr_str})...{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Check for Trifecta Probes:: Upstream of SNP
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ups_tib <- NULL
    ups_tib <- tib %>% 
      dplyr::filter(!!din_sym==tar_din) %>%
      dplyr::filter(up59=="C" & up60=="G") %>%
      dplyr::mutate(
        !!din_sym := "cg",
        !!ext_sym := paste0(up01,up02,up10,up11,up12,up58,up59,up60,  up61,dn61, dn60,dn59,dn58,dn12,dn11,dn10) %>% addBrac(), 
        !!iup_sym := paste0(up01,up02,up10,up11,up12,up58,up59,up60, iupac,dn61, dn60,dn59,dn58,dn12,dn11,dn10) %>% addBrac(), 
        !!imp_sym := paste0(up01,up02,up10,up11,up12,up58,               "CG",   up61,dn61,dn60,dn59,dn58,dn12,dn11,dn10) %>% addBrac(), 
        !!ref_col_sym := "C",
        !!alt_col_sym := "T",
        !!iup_col_sym := "Y",
        Din_Str=stringr::str_to_lower(paste0(up59,up60)),
        Des_Din="cg",
        !!ids_sym:=paste(!!ids_sym,"CG-UP", sep='-'),
        !!pos_sym:=as.integer(!!pos_sym - 2)    # Previously:: Coordinate=as.integer(Coordinate+2)
      )
    ups_key <- glue::glue("ups_tib({funcTag})")
    ups_cnt <- print_tib(ups_tib,funcTag, verbose,vt+4,tc, n=ups_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #               Check for Trifecta Probes:: Downstream of SNP
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt+1)
      cat(glue::glue("{mssg}{TAB} Adding Downstream CG ",
                     "Flank Seuqnces({chr_str})...{RET}"))
    
    dns_tib <- NULL
    dns_tib <- tib %>% 
      dplyr::filter(!!din_sym==tar_din) %>%
      dplyr::filter(dn61=="C" & dn60=="G") %>%
      dplyr::mutate(
        !!din_sym := "cg",
        !!ext_sym := paste0(up10,up11,up12,up58,up59,up60,up61,  dn61,dn60, dn59,dn58,dn12,dn11,dn10,dn02) %>% addBrac(),
        !!iup_sym := paste0(up10,up11,up12,up58,up59,up60,iupac, dn61,dn60, dn59,dn58,dn12,dn11,dn10,dn02) %>% addBrac(),
        !!imp_sym := paste0(up10,up11,up12,up58,up59,up60,up61,     "CG",   dn59,dn58,dn12,dn11,dn10,dn02) %>% addBrac(),
        !!ref_col_sym := "C",
        !!alt_col_sym := "T",
        !!iup_col_sym := "Y",
        Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
        Des_Din="cg",
        !!ids_sym:=paste(!!ids_sym,"CG-DN", sep='-'),
        !!pos_sym:=as.integer(!!pos_sym - 2)    # Previously:: Coordinate=as.integer(Coordinate+2)
      )
    dns_key <- glue::glue("dns_tib({funcTag})")
    dns_cnt <- print_tib(dns_tib,funcTag, verbose,vt+4,tc, n=dns_key)
    
    ret_tib <- dplyr::bind_rows(ups_tib,dns_tib)
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Genome FASTA Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_genome = function(chr_fas,
                       chr_max,
                       ret_map = TRUE,
                       chr_key = "Chrom",
                       verbose=0 ,vt=3,tc=1,tt=NULL,
                       funcTag='load_genome') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   chr_fas = {chr_fas}.{RET}"))
    cat(glue::glue("{mssg}   chr_max = {chr_max}.{RET}"))
    cat(glue::glue("{mssg}   chr_key = {chr_key}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  stime <- base::system.time({
    
    if (chr_max==0) {
      seqs <- Biostrings::readDNAStringSet(filepath=chr_fas, format="fasta")
    } else {
      seqs <-  Biostrings::readDNAStringSet(filepath=chr_fas, format="fasta", 
                                            nrec=nrec)
    }
    
    ret_cnt <- seqs %>% length()
    if (verbose>=vt) 
      cat(glue::glue("{mssg} Loaded {ret_cnt} chromosomes.{RET}"))
    
    if (verbose>=vt+2) print(seqs)
    
    if (ret_map) {
      
      chr_sym <- rlang::sym(chr_key)
      
      ret_tib <- seqs %>% 
        names() %>% 
        stringr::str_remove(" .*$") %>% 
        stringr::str_remove("^chr") %>%
        tibble::tibble() %>% 
        purrr::set_names(chr_key) %>% 
        dplyr::mutate(Idx=dplyr::row_number(),
                      !!chr_sym := paste0("chr",!!chr_sym) )
      
      ret_dat$seqs <- seqs
      ret_dat$maps <- ret_tib
      
      ret_key <- glue::glue("ret-FIN({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  if (ret_map) return(ret_dat)
  
  ret_tib
}

parse_genomic_seqs = function(tib, seq,
                              srd_str="F", 
                              pos_key="Coordinate",
                              chr_str="chr0",
                              ups_len=60,
                              seq_len=122,
                              pad_len=2,
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='parse_genomic_seqs') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Field Parameters::{RET}"))
    cat(glue::glue("{mssg}   srd_str={srd_str}{RET}"))
    cat(glue::glue("{mssg}   pos_key={pos_key}.{RET}"))
    cat(glue::glue("{mssg}   chr_str={chr_str}{RET}"))
    cat(glue::glue("{mssg}   ups_len={ups_len}{RET}"))
    cat(glue::glue("{mssg}   seq_len={seq_len}{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    pos_sym  <- rlang::sym(pos_key)
    
    cur_key <- glue::glue("cur_tib_1=chr_str={chr_str}")
    cur_cnt <- print_tib(tib,funcTag, verbose,vt+4,tc, n=cur_key)
    
    # Format Coordinates::
    #   Subtract two for cg upstream 122mer formation and add two+two as well
    begs <- NULL
    ends <- NULL
    begs <- tib %>% dplyr::pull(!!pos_sym) - ups_len - pad_len
    ends <- begs + seq_len - 1 + (2 * pad_len)
    
    if (verbose>=vt+6) {
      cat(glue::glue("{mssg}{TAB} begs={RET}"))
      begs %>% head() %>% print()
      cat(glue::glue("{mssg}{TAB} ends={RET}"))
      ends %>% head() %>% print()
    }
    
    # Sub-string ref sequences::
    seq_vec <- NULL
    seq_vec <- stringr::str_sub( as.character(seq), begs, ends ) %>%
      stringr::str_to_upper()
    ret_cnt <- seq_vec %>% length()
    
    # *** NOTE:: DO NOT use this if loading BSC genomes ***
    #
    #         [ they're already strand specific ]
    #
    if (srd_str=="R") seq_vec <- seq_vec %>% cmpl()
    
    if (verbose>=vt+4) {
      cat(glue::glue("{mssg} Parsed Sequences({ret_cnt}):{RET}"))
      seq_vec %>% head(n=2) %>% print()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  seq_vec
}

# End of file


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
#             Basic Functions to Loop Over Probe Design Results::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

compare_seqs = function(i,j, t, d) {
  tabA_tib <- t[i, ]
  tabB_tib <- t[j, ]
  
  keyA <- tabA_tib$Key
  keyB <- tabB_tib$Key
  
  cat(glue::glue("[{par$prgmTag}]: Current i={i}, keyA={keyA}.{RET}"))
  cat(glue::glue("[{par$prgmTag}]: Current j={j}, keyB={keyB}.{RET}"))
  cat(glue::glue("[{par$prgmTag}]:{RET}"))
  
  datA_tib <- d[[keyA]]
  datB_tib <- d[[keyB]]
  
  # rnA_cnt <- base::nrow(datA_tib[[keyA]])
  # rnB_cnt <- base::nrow(datA_tib[[keyB]])
  rnA_cnt <- base::nrow(datA_tib)
  rnB_cnt <- base::nrow(datB_tib)
  max_cnt <- max(rnA_cnt, rnB_cnt)
  
  # cat(glue::glue("[{par$prgmTag}]: datA_tib={RET}"))
  # datA_tib %>% print(n=3)
  # cat(glue::glue("[{par$prgmTag}]:{RET}"))
  # 
  # cat(glue::glue("[{par$prgmTag}]: datB_tib={RET}"))
  # datB_tib %>% print(n=3)
  # cat(glue::glue("[{par$prgmTag}]:{RET2}"))
  
  seqA_keys <- NULL
  if (tabA_tib$improbe_type=="c") {
    seqA_keys <- c("Probe_Seq_U","Probe_Seq_M")
  }
  if (tabA_tib$improbe_type=="r") {
    seqA_keys <- c("Prb_1U","Prb_1M","Prb_2D")
  }
  if (tabA_tib$improbe_type=="s") {
    seqA_keys <- c("Prb1C","Prb2C","Prb1O","Prb2O")
  }
  
  seqB_keys <- NULL
  if (tabA_tib$improbe_type=="c") {
    seqB_keys <- c("Probe_Seq_U","Probe_Seq_M")
  }
  if (tabB_tib$improbe_type=="r") {
    seqB_keys <- c("Prb_1U","Prb_1M","Prb_2D")
  }
  if (tabB_tib$improbe_type=="s") {
    seqB_keys <- c("Prb1C","Prb2C","Prb1O","Prb2O")
  }
  
  keysA_len <- seqA_keys %>% length()
  keysB_len <- seqB_keys %>% length()
  
  for (seqA in seqA_keys) {
    cat(glue::glue("[{par$prgmTag}]: keyA={keyA}({keysA_len}) vs. keyB={keyB}({keysB_len}).{RET}"))
    
    for (seqB in seqB_keys) {
      # cat(glue::glue("[{par$prgmTag}]:{TAB}SeqA={seqA}={RET}"))
      # datA_tib[[seqA]] %>% head(n=3) %>% print()
      
      # cat(glue::glue("[{par$prgmTag}]:{TAB}SeqB={seqB}={RET}"))
      # datB_tib[[seqB]] %>% head(n=3) %>% print()
      # cat("\t\tIntersection=\n")
      
      per_mat_cnt <-  -1
      per_mat_cnt <- which(datA_tib[[seqA]] %in% datB_tib[[seqB]] == TRUE) %>% length()
      per_mat_per <- round(100*per_mat_cnt/max_cnt)
      
      cat(glue::glue("[{par$prgmTag}]: ",
                     "match_cnt={per_mat_cnt}{TAB}",
                     "match_per={per_mat_per}{TAB}",
                     "rnA_cnt={rnA_cnt}{TAB}",
                     "rnB_cnt={rnB_cnt}{TAB}",
                     "max_cnt={max_cnt}{TAB}"))
      cat(glue::glue("{seqA}{TAB}{seqB}{RET}"))
    }
    cat(glue::glue("[{par$prgmTag}]:{RET}"))
  }
  
  cat(glue::glue("[{par$prgmTag}]:{RET}"))
}

files_to_tib = function(dir, pat) {
  
  list.files(dir, pattern = pat, full.names = TRUE) %>% 
    tibble::as_tibble() %>% purrr::set_names("Path") %>%  
    dplyr::mutate(Base_Name=base::basename(Path) %>% 
                    stringr::str_remove(".[ct]sv.gz$")) %>%
    tidyr::separate(Base_Name, into=c("Run_Tag","improbe_type"), sep="\\.",
                    remove = FALSE) %>% 
    dplyr::mutate(Alphabet=stringr::str_remove(Run_Tag, "^.*_"),
                  Bsc_Tag=stringr::str_remove(Run_Tag, "^.*-") %>%
                    stringr::str_remove("_.*$")) %>% 
    tidyr::separate(Bsc_Tag, into=c("Strand_FR","Strand_CO", "Strand_BSC"),
                    sep=c(1,2), remove=FALSE) %>% 
    dplyr::mutate(improbe_type=improbe_type %>% 
                    stringr::str_remove("_improbe_workflow$"),
                  Key=paste(improbe_type, Alphabet, Bsc_Tag, sep='_')) %>%
    dplyr::select(Key, improbe_type, Alphabet, Bsc_Tag, 
                  Strand_FR, Strand_CO, Strand_BSC, Path)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Loop Over Genomes Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This function should be completely written after individual modules
#    are cleaned up and re-written::
#
#    - s_improbe_workflow()
#    - r_improbe_workflow()
#    - c_improbe_workflow()
#

prb_designs_workflow = function(
  tib,
  
  # Output Parameters::
  out_dir,
  run_tag,
  pre_tag,
  re_load = TRUE,
  end_str = 'csv.gz',
  sep_chr='.',
  
  add_inf    = TRUE,
  imp_level  = 3,
  
  # Genomes Parameters::
  gen_bld  = "UNK",
  gen_nrec = 0,
  gen_key  = "Genome_Key",
  gen_tib  = NULL,
  
  # Field Parameters:: general
  ids_key = "Prb_Key",
  des_key = "Ord_Des",
  din_key = "Ord_Din", 
  
  # bsp_srd = "Bsp_FR",
  # bsp_cos = "Bsp_CO",
  pos_key = "Bsp_Pos", 
  chr_key = "Bsp_Chr", 
  
  # Field Parameters:: s-improbe
  ext_seq = "Ext_Forward_Seq",
  iup_seq = "Iupac_Forward_Sequence",
  imp_seq = "Forward_Sequence",
  iupac=NULL,
  
  # Field Parameters:: r-improbe
  srsplit = FALSE,
  srd_key = NULL,
  srd_str = "FR",
  
  cosplit = FALSE,
  cos_key = NULL,
  cos_str = "CO",
  
  # Docker Parameters::
  doc_image,
  doc_shell = "run_improbe.sh",
  
  join     = FALSE,
  join_new = c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos","Bsp_FR","Bsp_CO"),
  join_old = c("Seq_ID","Chromosome","Coordinate","Strand_FR","Strand_CO"),
  
  # Run Parameters::
  max = 0,
  del = "_",
  ups_len  = 60,
  seq_len  = 122,
  subset   = FALSE,
  sub_cols = NULL,
  
  reload   = FALSE,
  retData  = FALSE, 
  parallel = FALSE,
  
  r_improbe = FALSE,
  s_improbe = FALSE,
  c_improbe = FALSE,
  
  add_flanks = FALSE,
  add_probes = FALSE,
  add_matseq = TRUE,
  
  verbose=0,vt=3,tc=1,tt=NULL, funcTag='prb_designs_workflow') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} improbe Parameters::{RET}"))
    cat(glue::glue("{mssg}        out_dir={out_dir}.{RET}"))
    cat(glue::glue("{mssg}       run_tag={run_tag}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}        add_inf={add_inf}.{RET}"))
    cat(glue::glue("{mssg}      imp_level={imp_level}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Genome Parameters::{RET}"))
    cat(glue::glue("{mssg}        gen_bld={gen_bld}.{RET}"))
    cat(glue::glue("{mssg}       gen_nrec={gen_nrec}.{RET}"))
    cat(glue::glue("{mssg}        gen_key={gen_nrec}.{RET}"))
    print_tib(gen_tib, funcTag, verbose,vt=vt+8, n="gen_tib")
    # cat(glue::glue("{RET}"))
    # cat(glue::glue("{mssg}   Tmp_Strand_FR={Tmp_Strand_FR}.{RET}"))
    # cat(glue::glue("{mssg}   Tmp_Strand_CO={Tmp_Strand_CO}.{RET}"))
    # cat(glue::glue("{mssg}   Tmp_Strand_BC={Tmp_Strand_BC}.{RET}"))
    # cat(glue::glue("{mssg}   Tmp_Char_Code={Tmp_Char_Code}.{RET}"))
    # cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Field Parameters::{RET}"))
    cat(glue::glue("{mssg}        ids_key={ids_key}.{RET}"))
    cat(glue::glue("{mssg}        des_key={des_key}.{RET}"))
    cat(glue::glue("{mssg}        din_key={din_key}.{RET}"))
    # cat(glue::glue("{mssg}        bsp_srd={bsp_srd}.{RET}"))
    # cat(glue::glue("{mssg}        bsp_cos={bsp_cos}.{RET}"))
    cat(glue::glue("{mssg}        pos_key={pos_key}.{RET}"))
    cat(glue::glue("{mssg}        chr_key={chr_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}        ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("{mssg}        iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("{mssg}        imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{mssg}          iupac={iupac}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}        srsplit={srsplit}.{RET}"))
    cat(glue::glue("{mssg}        srd_key={srd_key}.{RET}"))
    cat(glue::glue("{mssg}        srd_str={srd_str}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}        cosplit={cosplit}.{RET}"))
    cat(glue::glue("{mssg}        cos_key={cos_key}.{RET}"))
    cat(glue::glue("{mssg}        cos_str={cos_str}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Docker Parameters::{RET}"))
    cat(glue::glue("{mssg}      doc_shell={doc_shell}.{RET}"))
    cat(glue::glue("{mssg}      doc_image={doc_image}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Joining Parameters::{RET}"))
    cat(glue::glue("{mssg}           join={join}.{RET}"))
    cat(glue::glue("{mssg}       join_new={join_new}.{RET}"))
    cat(glue::glue("{mssg}       join_old={join_old}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}            max={max}.{RET}"))
    cat(glue::glue("{mssg}            del={del}.{RET}"))
    cat(glue::glue("{mssg}        ups_len={ups_len}.{RET}"))
    cat(glue::glue("{mssg}        seq_len={seq_len}.{RET}"))
    cat(glue::glue("{mssg}         subset={subset}.{RET}"))
    cat(glue::glue("{mssg}       sub_cols={sub_cols}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}         reload={reload}.{RET}"))
    cat(glue::glue("{mssg}        retData={retData}.{RET}"))
    cat(glue::glue("{mssg}       parallel={parallel}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      r_improbe={r_improbe}.{RET}"))
    cat(glue::glue("{mssg}      s_improbe={s_improbe}.{RET}"))
    cat(glue::glue("{mssg}      c_improbe={c_improbe}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}     add_flanks={add_flanks}.{RET}"))
    cat(glue::glue("{mssg}     add_probes={add_probes}.{RET}"))
    cat(glue::glue("{mssg}     add_matseq={add_matseq}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ref_fas_list <- gen_tib %>%
    split(f=.[[gen_key]])
  
  ret_cnt <- names(ref_fas_list) %>% length()
  ret_tib <- NULL
  ret_dat <- NULL
  
  for (cur_gen_key in names(ref_fas_list)) {
    
    cur_gen_tib <- ref_fas_list[[cur_gen_key]] %>% head(n=1)
    
    cur_gen_fas <- cur_gen_tib$Path
    cur_gen_srdFR <- cur_gen_tib$Strand_FR  %>% as.character()
    cur_gen_cosCO <- cur_gen_tib$Strand_CO  %>% as.character()
    cur_gen_bscBC <- cur_gen_tib$Strand_BSC %>% as.character()
    cur_gen_Alpha <- cur_gen_tib$Alphabet   %>% as.character()
    
    cur_gen_tib <- cur_gen_tib %>% 
      dplyr::select(-Path) %>%
      dplyr::rename(Gen_Strand_FR=Strand_FR,
                    Gen_Strand_CO=Strand_CO,
                    Gen_Strand_BSC=Strand_BSC,
                    Gen_Alphabet=Alphabet,
                    Gen_Source=Source,
                    Gen_Genome_Build=Genome_Build)
    
    if (verbose>=vt) {
      cat(glue::glue("{mssg} cur_gen_srdFR={cur_gen_srdFR}.{RET}"))
      cat(glue::glue("{mssg} cur_gen_cosCO={cur_gen_cosCO}.{RET}"))
      cat(glue::glue("{mssg} cur_gen_bscBC={cur_gen_bscBC}.{RET}"))
      cat(glue::glue("{mssg} cur_gen_Alpha={cur_gen_Alpha}.{RET}"))
      cat(glue::glue("{mssg}   cur_gen_key={cur_gen_key}.{RET}"))
      cat(glue::glue("{mssg}   cur_gen_fas={cur_gen_fas}.{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    if (is.null(cur_gen_srdFR)) {
      stop(glue::glue("{RET}{mssg} cur_gen_srdFR is NULL!{RET2}"))
      return(ret_tib)
    }
    
    s_imp_tib <- NULL
    s_imp_tib <- 
      s_improbe_workflow(tib = tib,
                         
                         # Add genome source row::
                         gen_tib = cur_gen_tib,
                         
                         gen_bld = gen_bld,
                         gen_fas = cur_gen_fas,
                         
                         ids_key = ids_key,
                         des_key = des_key,
                         din_key = din_key,
                         
                         srd_key = srd_key,
                         cos_key = cos_key,
                         
                         tar_din = "rs",
                         
                         # srd_str = "F",
                         srd_str = cur_gen_srdFR,
                         pos_key = pos_key,
                         chr_key = chr_key,
                         
                         ext_seq = "Ext_Forward_Seq",
                         iup_seq = "Iupac_Forward_Sequence",
                         imp_seq = "Forward_Sequence",
                         iupac = NULL,
                         
                         ups_len = 60,
                         seq_len = 122,
                         
                         # del = "_",
                         # subset   = subset,
                         # sub_cols = sub_cols,
                         
                         reload  = reload,
                         # retData = retData,
                         # reload  = FALSE,
                         retData = FALSE,
                         
                         parallel   = parallel,
                         add_flanks = add_flanks,
                         add_probes = add_probes,
                         
                         out_dir = out_dir,
                         run_tag = paste(run_tag,cur_gen_key,sep='-'),
                         re_load = TRUE,
                         pre_tag = tt$file_vec,
                         
                         verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    if (cur_gen_tib$Gen_Strand_BSC=="N") {
      
      if (r_improbe) {
        #
        # When pulling from only forward strand we must fix the Forward FR Strand!!!
        #
        
        # If dna calculation Top Sequence
        # If snp do NOT calculation Top Sequence
        
        top_col <- "Strand_TB"
        top_key <- "Top_Sequence"
        add_topseq <- FALSE
        if (cur_gen_tib$Gen_Alphabet=="dna") 
          add_topseq <- TRUE
        
        r_imp_tib <- NULL
        r_imp_tib <- 
          r_improbe_workflow(tib = s_imp_tib,
                             
                             # Add genome source row::
                             gen_tib = cur_gen_tib,
                             
                             ids_key = ids_key,
                             seq_key = imp_seq,
                             din_key = din_key,
                             
                             top_col = top_col,
                             top_key = top_key,
                             
                             srsplit = srsplit,
                             srd_key = "Strand_FR",
                             srd_str = srd_str,
                             
                             cosplit = cosplit,
                             cos_key = "Strand_CO",
                             cos_str = cos_str,
                             
                             ups_len = ups_len,
                             seq_len = seq_len,
                             
                             # subset   = subset,
                             # sub_cols = sub_cols,
                             
                             reload     = reload,
                             parallel   = parallel,
                             add_matseq = add_matseq,
                             add_topseq = add_topseq,
                             
                             out_dir = out_dir,
                             run_tag = paste(run_tag,cur_gen_key,sep='-'),
                             re_load = TRUE,
                             pre_tag = tt$file_vec,
                             
                             verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      }
      
      if (c_improbe && cur_gen_tib$Gen_Alphabet == "dna") {
        
        #
        # TBD:: Update s_imp_tib with setTopBot_tib()
        # TBD:: rename setTopBot_tib() -> set_TopBot_tib() and all params
        #       This should be done in s_improbe_workflow() with a flag
        #        for set_TB=TRUE/FALSE and targets imp_seq key.
        #
        # TBD:: Validate s_imp_tib.Top_Sequence == c_imp_tib.Top_Sequence
        #        inside of c_improbe_workflow() if set_TB==TRUE
        #
        
        c_imp_tib <- NULL
        c_imp_tib <- c_improbe_workflow(imp_tib = s_imp_tib, 
                                        
                                        # Add genome source row::
                                        gen_tib = cur_gen_tib,
                                        
                                        ids_key = ids_key,
                                        imp_seq = imp_seq,                              
                                        pos_key = pos_key,
                                        chr_key = chr_key,
                                        gen_bld = gen_bld,
                                        
                                        doc_image = doc_image,
                                        doc_shell = doc_shell,
                                        
                                        level   = imp_level,
                                        add_inf = add_inf, 
                                        
                                        join     = join,
                                        join_new = join_new,
                                        join_old = join_old,
                                        
                                        reload = reload,
                                        
                                        out_dir = out_dir,
                                        run_tag = paste(run_tag,cur_gen_key,sep='-'),
                                        re_load = TRUE,
                                        pre_tag = tt$file_vec,
                                        end_str = 'tsv.gz',
                                        
                                        verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
        
        if (retData) ret_dat$c_imp <- c_imp_tib
      }
      
      # ret_key <- glue::glue("ret-FIN({funcTag})")
      # ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    }
    
  }
  
  etime <- 0
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
  
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

compare_probes = function(ref, can, ref_keys, can_keys,
                          by=c("Seq_ID","Strand_FR","Strand_CO",
                               "Ord_Des","Ord_Din"),
                          grp=c("Ord_Des","Ord_Din","Strand_FR","Strand_CO"),
                          retData=FALSE, pivot=FALSE,
                          verbose=0,vt=4,tc=1,tt=NULL,
                          funcTag='compare_probes') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  ret_tab <- NULL
  stime <- base::system.time({
    
    # Inner Join::
    inn_tib <- dplyr::inner_join(ref,can, by=by,suffix=c("_ref","_can"))
    ret_key <- glue::glue("inn-tib({funcTag})")
    ret_cnt <- print_tib(inn_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    for (ii in c(1:length(ref_keys))) {
      fieldA <- ref_keys[ii]
      fieldB <- can_keys[ii]
      
      # Fix any names that were duplicated::
      if (!fieldA %in% names(inn_tib))
        fieldA <- paste0(fieldA,"_ref")
      
      if (!fieldB %in% names(inn_tib))
        fieldB <- paste0(fieldB,"_can")
      
      if (verbose>=vt)
        cat(glue::glue("{mssg} ii={ii}, A={fieldA}, B={fieldB}{RET}"))
      
      sel_vec <- c(by, fieldA, fieldB)
      new_vec <- c(by, "Probe_A", "Probe_B", 
                   "Man_MisMatch", "DI_NUC_AB", "Man_TarMatch")
      
      ret_tib <- inn_tib %>% 
        dplyr::select(dplyr::all_of(sel_vec)) %>%
        purrr::set_names(sel_vec) %>%
        cmpInfII_MisMatch(fieldA=fieldA,fieldB=fieldB, 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
        purrr::set_names(new_vec) %>%
        dplyr::bind_rows(ret_tib)
      
      # Exact Match::
      #
      # inn_tib %>% 
      #   dplyr::select(dplyr::all_of(sel_vec)) %>%
      #   purrr::set_names(sel_vec) %>%
      #   cmpInfII(fieldA=fieldA,fieldB=fieldB, 
      #            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% print()
      
      if (verbose>=vt)
        cat(glue::glue("{mssg}{BRK}{RET2}"))
    }
    
    if (pivot) ret_tab <- ret_tib %>% tidyr::pivot_longer(
      cols=c("Probe_A","Probe_B"),names_to="Prb_Source", values_to="Probe_Seq")
    
    if (retData) ret_dat$cmp <- ret_tib
    if (retData) ret_dat$tab <- ret_tab
    
    sum_key_vec <- c("Man_MisMatch", "Man_TarMatch", "DI_NUC_AB")
    for (key in sum_key_vec) {
      if (verbose>=vt)
        cat(glue::glue("{mssg} Summarizing {key}...{RET}"))
      
      grp_vec <- c(key,grp)
      sum_tib <- ret_tib %>% 
        dplyr::group_by(dplyr::across(dplyr::all_of(grp_vec)) ) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      
      if (retData) ret_dat[[key]] <- sum_tib
      
      if (verbose>=vt)
        cat(glue::glue("{mssg} Done.Summarizing {key}.{RET2}"))
    }
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      MisMatch Probe Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cmpInfIMU_MisMatch = function(tib, fieldA, fieldB, mu, del='_',
                              verbose=0,vt=5,tc=1,tt=NULL,
                              funcTag='cmpInfIMU_MisMatch') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}:      mu={mu}.{RET}"))
    cat(glue::glue("{mssg}:  fieldA={fieldA}.{RET}"))
    cat(glue::glue("{mssg}:  fieldB={fieldB}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(
      BOD_NumMM=mapply(
        adist,
        stringr::str_sub(!!fieldA,1,stringr::str_length(!!fieldA)-1),
        stringr::str_sub(!!fieldB,1,stringr::str_length(!!fieldB)-1),
        ignore.case=TRUE ),
      
      DI_NUC_AB=paste0(
        stringr::str_to_upper(
          stringr::str_sub(!!fieldA,
                           stringr::str_length(!!fieldA),
                           stringr::str_length(!!fieldA)) ),
        stringr::str_to_upper(
          stringr::str_sub(!!fieldB,
                           stringr::str_length(!!fieldB),
                           stringr::str_length(!!fieldB)) )
      ),
      TAR_EQU=cmpIUPACs(DI_NUC_AB)
    ) %>%
    dplyr::rename(!!paste('BOD_NumMM',mu, sep=del):=BOD_NumMM,
                  !!paste('TAR_EQU',  mu, sep=del):=TAR_EQU)
  
  if (verbose>=vt) cat(glue::glue("{mssg} Done.{RET2}"))
  
  tib
}

cmpInfI_MisMatch = function(tib, fieldAU, fieldBU, fieldAM, fieldBM, del='_',
                            verbose=0,vt=5,tc=1,tt=NULL,
                            funcTag='cmpInfI_MisMatch') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg} fieldAU={fieldAU}.{RET}"))
    cat(glue::glue("{mssg} fieldBU={fieldBU}.{RET}"))
    cat(glue::glue("{mssg} fieldAM={fieldAM}.{RET}"))
    cat(glue::glue("{mssg} fieldBM={fieldBM}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  tib <- tib %>% 
    cmpInfIMU_MisMatch(fieldAU, fieldBU, mu='U', del=del,
                       verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% 
    cmpInfIMU_MisMatch(fieldAM, fieldBM, mu='M', del=del,
                       verbose=verbose, vt=vt+1,tc=tc+1)
  
  tib <- tib %>% dplyr::mutate(
    Man_MisMatch=(BOD_NumMM_U+BOD_NumMM_M)/2, #, na.rm=TRUE),
    Man_TarMatch=case_when(TAR_EQU_U & TAR_EQU_M ~ TRUE, TRUE ~ FALSE) ) %>%
    dplyr::select(-c(BOD_NumMM_U,BOD_NumMM_M,TAR_EQU_U,TAR_EQU_M))
  
  if (verbose>=vt) cat(glue::glue("{mssg} Done.{RET2}"))
  
  tib
}

cmpInfII_MisMatch = function(tib, fieldA, fieldB, mu='D', del='_',
                             verbose=0,vt=4,tc=1,tt=NULL,
                             funcTag='cmpInfII_MisMatch') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{mssg} fieldA={fieldA}.{RET}"))
    cat(glue::glue("{mssg} fieldB={fieldB}.{RET2}"))
  }
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    cmpInfIMU_MisMatch(fieldA, fieldB, mu=mu, del=del,
                       verbose=verbose, vt=vt+1,tc=tc+1) %>%
    dplyr::rename(
      Man_MisMatch=BOD_NumMM_D,
      Man_TarMatch=TAR_EQU_D)
  # dplyr::select(-c(BOD_NumMM_U,BOD_NumMM_M,TAR_EQU_U,TAR_EQU_M))
  
  if (verbose>=vt) cat(glue::glue("{mssg} Done.{RET2}"))
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Exact Probe Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cmpInfIMU= function(tib, fieldA, fieldB, mu, del='_',
                    verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfMU'
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} fieldA={fieldA}, fieldB={fieldB} mu={mu}{RET}"))
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(
      SUB_SEQ_A=stringr::str_to_upper(stringr::str_sub(!!fieldA,1,stringr::str_length(!!fieldA)-1)),
      SUB_SEQ_B=stringr::str_to_upper(stringr::str_sub(!!fieldB,1,stringr::str_length(!!fieldB)-1)),
      DI_NUC_AB=paste0(
        stringr::str_to_upper(stringr::str_sub(!!fieldA,stringr::str_length(!!fieldA),stringr::str_length(!!fieldA)) ),
        stringr::str_to_upper(stringr::str_sub(!!fieldB,stringr::str_length(!!fieldB),stringr::str_length(!!fieldB)) )
      ),
      
      BOD_EQU=case_when(SUB_SEQ_A==SUB_SEQ_B ~ TRUE, TRUE ~ FALSE),
      TAR_EQU=cmpIUPACs(DI_NUC_AB),
      Inf1_Match=case_when(BOD_EQU & BOD_EQU==TAR_EQU ~ TRUE, TRUE ~ FALSE)
    ) %>%
    dplyr::select(-c(SUB_SEQ_A,SUB_SEQ_B,BOD_EQU,DI_NUC_AB,TAR_EQU)) %>%
    dplyr::rename(!!paste('Inf1_Match',mu, sep=del):=Inf1_Match)
  
  if (verbose>=vt+1) print(tib)
  
  tib
}

cmpInfI = function(tib, fieldAU, fieldBU, fieldAM, fieldBM, del='_',
                   verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfI'
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{tabs} fieldAU={fieldAU}, fieldBU={fieldBU}{RET}"))
  if (verbose>=vt) cat(glue::glue("{tabs} fieldAM={fieldAM}, fieldBM={fieldBM}{RET}"))
  
  tib <- tib %>% cmpInfIMU(fieldAU, fieldBU, mu='U', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% cmpInfIMU(fieldAM, fieldBM, mu='M', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>%
    dplyr::mutate(Man_Match=case_when(Inf1_Match_U & Inf1_Match_M ~ TRUE, TRUE ~ FALSE) ) %>%
    dplyr::select(-c(Inf1_Match_U,Inf1_Match_M))
  
  tib
}

cmpInfII = function(tib, fieldA, fieldB, mu='D', del='_',
                    verbose=0,vt=4,tc=1,tt=NULL, funcTag='cmpInfI') {
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(Man_Match=stringr::str_to_upper(!!fieldA)==stringr::str_to_upper(!!fieldB) )
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Illumina Strand Methods:: TOP/BOT
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Formally called:: setTopBot_tib()
set_topbot_tib = function(tib, 
                          seq_key, 
                          top_col, 
                          top_key = NULL, 
                          max = 0,
                          
                          verbose=0, vt=4,tc=1,tt=NULL, 
                          funcTag='set_topbot_tib') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}      seq_key={seq_key}.{RET}"))
    cat(glue::glue("{mssg}      top_col={top_col}.{RET}"))
    cat(glue::glue("{mssg}      top_key={top_key}.{RET}"))
    cat(glue::glue("{mssg}          max={max}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  seq_key_sym <- rlang::sym(seq_key)
  top_col_sym <- rlang::sym(top_col)
  if (!is.null(top_key)) top_key_sym <- rlang::sym(top_key)
  
  # TBD:: Should remove any seq_key that is na or not the correct length...
  etime <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    if (max!=0) tib <- tib %>% head(n=max)
    
    # Make template sequences unique to reduce run time()
    unq_tib <- tib %>% 
      dplyr::distinct(dplyr::across(dplyr::all_of( seq_key)) )
    
    bit_tib <- unq_tib %>%
      dplyr::mutate(
        pre_seq = !!seq_key_sym %>% 
          stringr::str_remove('\\[.*$') %>% 
          Biostrings::reverse() %>%
          stringr::str_to_upper() %>% 
          stringr::str_replace_all('A', '1') %>% 
          stringr::str_replace_all('T', '1') %>% 
          stringr::str_replace_all('C', '2') %>% 
          stringr::str_replace_all('G', '2'),
        
        pos_seq = !!seq_key_sym %>% 
          stringr::str_remove('^.*\\]') %>%
          stringr::str_to_upper() %>% 
          stringr::str_replace_all('A', '1') %>% 
          stringr::str_replace_all('T', '1') %>% 
          stringr::str_replace_all('C', '2') %>% 
          stringr::str_replace_all('G', '2')
      )
    
    bit_key <- glue::glue("bit-tib({funcTag})")
    bit_cnt <- print_tib(bit_tib,funcTag, verbose,vt+4,tc, n=bit_key)
    
    pre_bit <- bit_tib$pre_seq %>% 
      stringr::str_split('', simplify = TRUE)
    pre_key <- glue::glue("pre-bit({funcTag})")
    pre_cnt <- print_tib(tibble::as_tibble(pre_bit), funcTag, 
                         verbose,vt+4,tc, n=pre_key)
    
    pos_bit <- bit_tib$pos_seq %>% stringr::str_split('', simplify = TRUE)
    pos_key <- glue::glue("pos-bit({funcTag})")
    pos_cnt <- print_tib(tibble::as_tibble(pos_bit), funcTag, 
                         verbose,vt+4,tc, n=pos_key)
    
    pre_mat <- pre_bit %>% as.data.frame() %>% data.matrix()
    pos_mat <- pos_bit %>% as.data.frame() %>% data.matrix()
    
    if (verbose>=vt+4) {
      cat(glue::glue("{mssg} dim(pre_mat) = {dim(pre_mat)}{RET}"))
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} dim(pos_mat) = {dim(pos_mat)}{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    dif_mat <- pre_mat-pos_mat
    dif_mat <- dif_mat %>% cbind(rep(2, dim(dif_mat)[1]))
    
    if (verbose>=vt+4) {
      cat(glue::glue("{mssg} dim(dif_mat) = {dim(dif_mat)}{RET}"))
      cat(glue::glue("{mssg} dif_mat = "))
      dif_mat %>% head(n=2) %>% print()
      cat(glue::glue("{RET}"))
    }
    
    # This is just for a sanity check...
    #  dif_mat[1,1:60] = 0
    bit_vec <- apply(dif_mat,1, function(x) head(x[x!=0],1)) + 1
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}{TAB} bit_vec = {head(bit_vec,n=3)}{RET}"))

    # Convert numeric values back to characters (T=TOP,B=BOT,U=Unknown)
    tb_vec <- bit_vec %>% 
      as.character() %>% 
      stringr::str_replace('0', 'T') %>% 
      stringr::str_replace('2', 'B') %>%
      stringr::str_replace('3', 'U')
    
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}{TAB} tb_vec = {head(tb_vec,n=3)}{RET}"))
    
    unq_tib <- unq_tib %>% dplyr::mutate(!!top_col := tb_vec)
      
    if (!is.null(top_key)) {
      unq_tib <- unq_tib %>% 
        dplyr::mutate(
          rev_seq = shearBrac(!!seq_key_sym) %>%
            revCmp() %>% addBrac(),
          !!top_key_sym := dplyr::case_when(
            !!top_col_sym=='T' ~ !!seq_key_sym,
            !!top_col_sym=='B' ~ rev_seq,
            TRUE ~ NA_character_)
        ) %>% dplyr::select(-rev_seq)
    }
    ret_tib <- dplyr::left_join(tib, unq_tib, by=c(seq_key))
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# Single variable method that is OUT-OF-DATE; Should use the method above...
isTopBot_single = function(seq, verbose=0, vt=4) {
  seq <- seq %>% stringr::str_to_upper() %>% 
    stringr::str_replace_all('T', 'A') %>% 
    stringr::str_replace_all('G', 'C')
  
  seq_vec <- stringr::str_split(seq,'\\[', simplify=TRUE) %>% as.vector()
  pre_seq <- seq_vec[1] %>% Biostrings::reverse()
  seq_vec <- stringr::str_split(seq,'\\]', simplify=TRUE) %>% as.vector()
  pos_seq <- seq_vec[2]
  
  pre_vec <- pre_seq %>% stringr::str_split('', simplify=TRUE) %>% as.vector()
  pos_vec <- pos_seq %>% stringr::str_split('', simplify=TRUE) %>% as.vector()
  
  if (verbose>=vt) {
    print(seq)
    print(pre_seq)
    print(pos_seq)
  }
  
  min_len <- min(stringr::str_length(pre_seq),stringr::str_length(pos_seq))
  for (ii in c(1:min_len)) {
    if (pre_vec[ii]==pos_vec[ii]) next
    if (pre_vec[ii]!=pos_vec[ii]) {
      if (pre_vec[ii]=='A') return('T')
      return('B')
    }
  }
  return('U')
}

isTopBots = function(x, verbose=0, vt=4) {
  x <- lapply(x,isTopBot_single)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Basic Bisulfite Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

test_cpp = function(seqs = c("CCCATCGCGTGGGAACttttCC", "CCCATCGCGTGGGAACttttCCC"),
                    uc = FALSE, vb = 0) {
  
  # Verbose Check
  if (vb > 1) {
    bscU_cpp( seqs = seqs, uc = uc, vb = vb )
    bscM_cpp( seqs = seqs, uc = uc, vb = vb )
    bscD_cpp( seqs = seqs, uc = uc, vb = vb )
    
    rev_cpp( seqs = seqs, uc = uc, vb = vb )
    cmp_cpp( seqs = seqs, uc = uc, vb = vb )
    revCmp_cpp( seqs = seqs, uc = uc, vb = vb )
    
  } else {
    # Quick check upper case = TRUE/FALSE

    bscU_cpp( seqs = seqs, uc = TRUE,  vb = vb ) %>% print()
    bscU_cpp( seqs = seqs, uc = FALSE, vb = vb ) %>% print()
    bscM_cpp( seqs = seqs, uc = TRUE,  vb = vb ) %>% print()
    bscM_cpp( seqs = seqs, uc = FALSE, vb = vb ) %>% print()
    bscD_cpp( seqs = seqs, uc = TRUE,  vb = vb ) %>% print()
    bscD_cpp( seqs = seqs, uc = FALSE, vb = vb ) %>% print()
    
    rev_cpp( seqs = seqs, uc = TRUE,  vb = vb ) %>% print()
    rev_cpp( seqs = seqs, uc = FALSE, vb = vb ) %>% print()
    cmp_cpp( seqs = seqs, uc = TRUE,  vb = vb ) %>% print()
    cmp_cpp( seqs = seqs, uc = FALSE, vb = vb ) %>% print()
    revCmp_cpp( seqs = seqs, uc = TRUE,  vb = vb ) %>% print()
    revCmp_cpp( seqs = seqs, uc = FALSE, vb = vb ) %>% print()
    
  }
  
}

MAPDi = function(x) {
  if (length(MAP_DI[[x]])==0) return(NA)
  MAP_DI[[x]]
}
mapDIs = function(x) {
  x <- lapply(x, MAPDi) %>% BiocGenerics::unlist()
}

# De-Methylate Probe Sequence for 4-base aligners
#
deM = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  
  # aln_seq=prb_seq %>%
  #   stringr::str_replace_all("R","A") %>% # A/G
  #   stringr::str_replace_all("Y","T") %>% # C/T
  #   
  #   stringr::str_replace_all("S","C") %>% # G/C
  #   stringr::str_replace_all("W","A") %>% # A/T
  #   stringr::str_replace_all("K","T") %>% # G/T
  #   stringr::str_replace_all("M","A") %>% # A/C
  #   
  #   stringr::str_replace_all("B","T") %>% # C/G/T
  #   stringr::str_replace_all("D","A") %>% # A/G/T
  #   stringr::str_replace_all("H","A") %>% # A/C/T
  #   stringr::str_replace_all("V","A") %>% # A/C/G
  #   
  #   stringr::str_replace_all("N","A"), # A/C/T/G
  
  if (uc) x <- tr(x, 'RYSWKMBDHV', 'ATCATATAAA')
  else    x <- tr(x, 'RYSWKMBDHV', 'atcatataaa')
  x
}
deMs = function(x, uc=FALSE) { deM(x, uc) }

bscU = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  if (uc) x <- tr(x, 'CYSMBHV', 'TTKWKWD')
  else    x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  x
}
bscUs = function(x, uc=FALSE) { bscU(x, uc) }

MAPM = function(x) {
  if (length(MAP_M[[x]])==0) return(x)
  MAP_M[[x]]
}
bscM = function(x) { stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPM) }
bscMs = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  x <- lapply(x, bscM) %>% BiocGenerics::unlist()
  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  if (uc) x <- stringr::str_to_upper(x)
  x
}

MAPD = function(x) {
  if (length(MAP_D[[x]])==0) return(x)
  MAP_D[[x]]
}
bscD = function(x) { stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPD) }
bscDs = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  x <- lapply(x, bscD) %>% BiocGenerics::unlist()
  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  if (uc) x <- stringr::str_to_upper(x)
  x
}


bscXs = function(x, umd_key, uc = TRUE,
                 vb=0, vt=3, tc=1, tt=NULL,
                 fun_tag='bscXs') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  tar_umd = x[[umd_key]] %>% unique()
  # cat("\n\tumd =",umd_key,"\n")
  
  ret_tib <- NULL
  if (tar_umd == "U") {
    ret_tib <- x %>% 
      dplyr::bind_cols(
        as.data.frame( bscUs( x$Bsp_Ref, uc = uc ) ) %>% 
          purrr::set_names(paste0("Bsp_Ref_Bsc") ) )
    
  } else if (tar_umd == "M") {
    ret_tib <- x %>% 
      dplyr::bind_cols(
        as.data.frame( bscMs( x$Bsp_Ref, uc = uc ) ) %>% 
          purrr::set_names(paste0("Bsp_Ref_Bsc") ) )
    
    
  } else if (tar_umd == "D") {
    ret_tib <- x %>% 
      dplyr::bind_cols(
        as.data.frame( bscDs( x$Bsp_Ref, uc = uc ) ) %>% 
          purrr::set_names(paste0("Bsp_Ref_Bsc") ) )
    
  } else {
    fail_mssg <- glue::glue("Unsupported Target umd = {tar_umd}")
    stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    return(NULL)
  }
}

QMAP = function(x, mu) {
  if (mu=='U') {
    return(QMAP_U[[x]])
  } else if (mu=='M') {
    return(QMAP_M[[x]])
  }
  x
}
qmaps = function(x, mu) {
  x <- lapply(x, QMAP, mu) %>% BiocGenerics::unlist()
}

cmpIUPAC = function(x) {
  if (base::is.element(x, names(IUPAC_EQ) )) return(IUPAC_EQ[[x]])
  # if (is.null(IUPAC_EQ[[x]])) return(FALSE)
  # if (length(IUPAC_EQ[[x]])==0) return(FALSE)
  FALSE
}

cmpIUPACs = function(x) {
  x <- lapply(x, cmpIUPAC) %>% BiocGenerics::unlist()
}

shearBrac = function(x) {
  x %>% stringr::str_remove('\\[') %>% stringr::str_remove('\\]')
}

shear_brac = function(x) {
  x %>% stringr::str_remove('\\[') %>% stringr::str_remove('\\]')
}

add_brac = function(seq) {
  seq_len <- stringr::str_length(seq)
  mid_len <- 2
  snp_idx <- as.integer(seq_len / 2)
  
  pre_beg <- 1
  pre_end <- snp_idx-1
  
  mid_beg <- snp_idx
  mid_end <- mid_beg+mid_len-1
  
  pos_beg <- mid_end+1
  pos_end <- seq_len
  
  pre_seq <- seq %>% stringr::str_sub(pre_beg,pre_end)
  mid_seq <- seq %>% stringr::str_sub(mid_beg,mid_end)
  pos_seq <- seq %>% stringr::str_sub(pos_beg,pos_end)
  
  new_seq <- paste0(pre_seq,'[',mid_seq,']',pos_seq)
  
  new_seq
}

addBrac = function(seq) {
  seq_len <- stringr::str_length(seq)
  mid_len <- 2
  snp_idx <- as.integer(seq_len / 2)
  
  pre_beg <- 1
  pre_end <- snp_idx-1
  
  mid_beg <- snp_idx
  mid_end <- mid_beg+mid_len-1
  
  pos_beg <- mid_end+1
  pos_end <- seq_len
  
  pre_seq <- seq %>% stringr::str_sub(pre_beg,pre_end)
  mid_seq <- seq %>% stringr::str_sub(mid_beg,mid_end)
  pos_seq <- seq %>% stringr::str_sub(pos_beg,pos_end)
  
  new_seq <- paste0(pre_seq,'[',mid_seq,']',pos_seq)
  
  new_seq
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Basic Reverse/Complement Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
revCmp = function(x) {
  Biostrings::reverse(x) %>% cmpl()
}

revSeq = function(x) {
  Biostrings::reverse(x)
}

cmpl_srd = function(x) {
  if (x=='F' || x == 'f') return('R')
  if (x=='R' || x == 'r') return('F')
  if (x=='T' || x == 't') return('B')
  if (x=='B' || x == 'b') return('T')
  if (x=='C' || x == 'c') return('O')
  if (x=='O' || x == 'o') return('C')
  
  x
}

cmpl_FR = function(x) {
  if (x=='F' || x == 'f') return('R')
  if (x=='R' || x == 'r') return('F')
  
  x
}

cmpl_TB = function(x) {
  if (x=='T' || x == 't') return('B')
  if (x=='B' || x == 'b') return('T')
  
  x
}

cmpl_CO = function(x) {
  if (x=='C' || x == 'c') return('O')
  if (x=='O' || x == 'o') return('C')
  
  x
}

cmpl = function(x) {
  tr(x, 'ACTGRYSWKMBDHVactgryswkmbdhv[]','TGACYRSWMKVHDBtgacyrswmkvhdb][')
  #
  #  ACTGRYSWKMBDHVactgryswkmbdhv[]
  #  TGACYRSWMKVHDBtgacyrswmkvhdb][
  #
  # x <- tr(x, 'ACTGRYSWKMBDHV','TGACYRSWMKVHDB')
  # x <- tr(x, 'actgryswkmbdhv','tgacyrswmkvhdb')
  # x
}

tr = function(x, old, new) {
  Biostrings::chartr(old, new, x)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Improbe Header Mapping Structure::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

INIT_IMP_HEADER = function() {
  IMP_TYPES = cols(
    Seq_ID = col_character(),
    Forward_Sequence = col_character(),
    Genome_Build = col_character(),
    Chromosome = col_character(),
    Coordinate = col_double(),
    Design_State = col_character(),
    Seq_Length = col_double(),
    Forward_CpG_Coord = col_double(),
    TB_Strand = col_character(),
    Top_Sequence = col_character(),
    Top_CpG_Coord = col_double(),
    Probe_Type = col_character(),
    Probeset_ID = col_character(),
    Probeset_Score = col_double(),
    Methyl_Probe_ID = col_character(),
    Methyl_Probe_Sequence = col_character(),
    Methyl_Probe_Length = col_double(),
    Methyl_Start_Coord = col_double(),
    Methyl_End_Coord = col_double(),
    Methyl_Probe_Covered_Top_Sequence = col_character(),
    Methyl_Allele_FR_Strand = col_character(),
    Methyl_Allele_TB_Strand = col_character(),
    Methyl_Allele_CO_Strand = col_character(),
    Methyl_Allele_Type = col_character(),
    Methyl_Final_Score = col_double(),
    Methyl_Tm = col_double(),
    Methyl_Tm_Score = col_double(),
    Methyl_GC_Percent = col_double(),
    Methyl_GC_Score = col_double(),
    Methyl_13mer_Count = col_double(),
    Methyl_13mer_Score = col_double(),
    Methyl_Address_Count = col_double(),
    Methyl_Address_Score = col_double(),
    Methyl_Self_Complementarity = col_double(),
    Methyl_Self_Complementarity_Score = col_double(),
    Methyl_Mono_Run = col_double(),
    Methyl_Mono_Run_Score = col_double(),
    Methyl_Ectopic_Count = col_double(),
    Methyl_Ectopic_Score = col_double(),
    Methyl_Underlying_CpG_Count = col_double(),
    Methyl_Underlying_CpG_Min_Dist = col_double(),
    Methyl_Underlying_CpG_Score = col_double(),
    Methyl_In_CpG_Island_Relaxed = col_logical(),
    Methyl_CpG_Island_Score = col_double(),
    Methyl_Next_Base = col_character(),
    Methyl_Next_Base_Score = col_double(),
    UnMethyl_Probe_ID = col_character(),
    UnMethyl_Probe_Sequence = col_character(),
    UnMethyl_Probe_Length = col_double(),
    UnMethyl_Start_Coord = col_double(),
    UnMethyl_End_Coord = col_double(),
    UnMethyl_Probe_Covered_Top_Sequence = col_character(),
    UnMethyl_Allele_FR_Strand = col_character(),
    UnMethyl_Allele_TB_Strand = col_character(),
    UnMethyl_Allele_CO_Strand = col_character(),
    UnMethyl_Allele_Type = col_character(),
    UnMethyl_Final_Score = col_double(),
    UnMethyl_Tm = col_double(),
    UnMethyl_Tm_Score = col_double(),
    UnMethyl_GC_Percent = col_double(),
    UnMethyl_GC_Score = col_double(),
    UnMethyl_13mer_Count = col_double(),
    UnMethyl_13mer_Score = col_double(),
    UnMethyl_Address_Count = col_double(),
    UnMethyl_Address_Score = col_double(),
    UnMethyl_Self_Complementarity = col_double(),
    UnMethyl_Self_Complementarity_Score = col_double(),
    UnMethyl_Mono_Run = col_double(),
    UnMethyl_Mono_Run_Score = col_double(),
    UnMethyl_Ectopic_Count = col_double(),
    UnMethyl_Ectopic_Score = col_double(),
    UnMethyl_Underlying_CpG_Count = col_double(),
    UnMethyl_Underlying_CpG_Min_Dist = col_double(),
    UnMethyl_Underlying_CpG_Score = col_double(),
    UnMethyl_In_CpG_Island_Relaxed = col_logical(),
    UnMethyl_CpG_Island_Score = col_double(),
    UnMethyl_Next_Base = col_character(),
    UnMethyl_Next_Base_Score = col_double()
  )
  
  IMP_TYPES
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Nucelotide Mapping Structures::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# mapA = function(x) {
#   if (length(MAP_A[[x]])==0) return(x)
#   MAP_A[[x]]
# }
# INIT_MAP_A = function() {
#   MAP_A <- NULL
#   MAP_A[['AA']] <- 'aa'
#   MAP_A[['Aa']] <- 'aa'
#   MAP_A[['aA']] <- 'aa'
#   MAP_A[['aa']] <- 'aa'
#   
#   MAP_A
# }
# MAP_A <- INIT_MAP_A()

# Di-Nuc to IUPAC::
INIT_MAP_DI = function() {
  MAP <- NULL
  
  # Mono-Nuc Map::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'
  
  # Di-Nuc Map::
  MAP[['AC']] <- 'M'
  MAP[['AG']] <- 'R'
  MAP[['AT']] <- 'W'
  MAP[['AA']] <- 'A'
  
  MAP[['CA']] <- 'M'
  MAP[['CT']] <- 'Y'
  MAP[['CG']] <- 'S'
  MAP[['CC']] <- 'C'
  
  MAP[['GA']] <- 'R'
  MAP[['GT']] <- 'K'
  MAP[['GC']] <- 'S'
  MAP[['GG']] <- 'G'
  
  MAP[['TC']] <- 'Y'
  MAP[['TG']] <- 'K'
  MAP[['TA']] <- 'W'
  MAP[['TT']] <- 'T'
  
  #
  # Tri-Nuc Map:: A
  #
  MAP[['AAC']] <- 'M'
  MAP[['AAG']] <- 'R'
  MAP[['AAT']] <- 'W'
  MAP[['AAA']] <- 'A'
  
  MAP[['ACA']] <- 'M'
  MAP[['ACT']] <- 'H' # NEW-3
  MAP[['ACG']] <- 'V' # NEW-3
  MAP[['ACC']] <- 'M' # NEW-2
  
  MAP[['AGA']] <- 'R'
  MAP[['AGT']] <- 'D' # NEW-3
  MAP[['AGC']] <- 'V' # NEW-3
  MAP[['AGG']] <- 'R' # NEW-2
  
  MAP[['ATC']] <- 'H' # NEW-3
  MAP[['ATG']] <- 'D' # NEW-3
  MAP[['ATA']] <- 'W'
  MAP[['ATT']] <- 'W' # NEW-2
  
  #
  # Tri-Nuc Map:: C
  #
  MAP[['CAC']] <- 'M'
  MAP[['CAG']] <- 'V' # NEW-3
  MAP[['CAT']] <- 'H' # NEW-3
  MAP[['CAA']] <- 'M' # NEW-2
  
  MAP[['CCA']] <- 'M'
  MAP[['CCT']] <- 'Y'
  MAP[['CCG']] <- 'S'
  MAP[['CCC']] <- 'C'
  
  MAP[['CGA']] <- 'V' # NEW-3
  MAP[['CGT']] <- 'B' # NEW-3
  MAP[['CGC']] <- 'S'
  MAP[['CGG']] <- 'S' # NEW-2
  
  MAP[['CTC']] <- 'Y'
  MAP[['CTG']] <- 'B' # NEW-3
  MAP[['CTA']] <- 'H' # NEW-3
  MAP[['CTT']] <- 'Y' # NEW-2
  
  #
  # Tri-Nuc Map:: G
  #
  MAP[['GAC']] <- 'V' # NEW-3
  MAP[['GAG']] <- 'R'
  MAP[['GAT']] <- 'D' # NEW-3
  MAP[['GAA']] <- 'R' # NEW-2
  
  MAP[['GCA']] <- 'V' # NEW-3
  MAP[['GCT']] <- 'B' # NEW-3
  MAP[['GCG']] <- 'S'
  MAP[['GCC']] <- 'S' # NEW-2
  
  MAP[['GGA']] <- 'R'
  MAP[['GGT']] <- 'K'
  MAP[['GGC']] <- 'S'
  MAP[['GGG']] <- 'G'
  
  MAP[['GTC']] <- 'B' # NEW-3
  MAP[['GTG']] <- 'K'
  MAP[['GTA']] <- 'D' # NEW-3
  MAP[['GTT']] <- 'K' # NEW-2
  
  #
  # Tri-Nuc Map:: T
  #
  MAP[['TAC']] <- 'H' # NEW-3
  MAP[['TAG']] <- 'D' # NEW-3
  MAP[['TAT']] <- 'W'
  MAP[['TAA']] <- 'W' # NEW-2
  
  MAP[['TCA']] <- 'H' # NEW-3
  MAP[['TCT']] <- 'Y'
  MAP[['TCG']] <- 'B' # NEW-3
  MAP[['TCC']] <- 'Y' # NEW-2
  
  MAP[['TGA']] <- 'D' # NEW-3
  MAP[['TGT']] <- 'K'
  MAP[['TGC']] <- 'B' # NEW-3
  MAP[['TGG']] <- 'K' # NEW-2
  
  MAP[['TTC']] <- 'Y'
  MAP[['TTG']] <- 'K'
  MAP[['TTA']] <- 'W'
  MAP[['TTT']] <- 'T'
  
  acgt <- c('A','C','G','T')
  for ( n1 in acgt ) {
    for ( n2 in acgt ) {
      for ( n3 in acgt ) {
        for ( n4 in acgt ) {
          n <- paste0(n1,n2,n3,n4)
          MAP[[n]] <- 'N'
        }
      }
    }
  }
  
  MAP
}
MAP_DI <- INIT_MAP_DI()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Initialize M/D Maps for Bisulfite Conversion::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Data generated with:: /Users/bbarnes/Documents/Projects/scripts/mapMD.pl
INIT_MAP_M = function() {
  MAP <- NULL
  MAP[['CG']] <- 'cG'
  MAP[['CR']] <- 'cR'
  MAP[['CS']] <- 'cS'
  MAP[['CK']] <- 'cK'
  MAP[['CB']] <- 'cB'
  MAP[['CD']] <- 'cD'
  MAP[['CV']] <- 'cV'
  MAP[['YG']] <- 'yG'
  MAP[['YR']] <- 'yR'
  MAP[['YS']] <- 'yS'
  MAP[['YK']] <- 'yK'
  MAP[['YB']] <- 'yB'
  MAP[['YD']] <- 'yD'
  MAP[['YV']] <- 'yV'
  MAP[['SG']] <- 'sG'
  MAP[['SR']] <- 'sR'
  MAP[['SS']] <- 'sS'
  MAP[['SK']] <- 'sK'
  MAP[['SB']] <- 'sB'
  MAP[['SD']] <- 'sD'
  MAP[['SV']] <- 'sV'
  MAP[['MG']] <- 'mG'
  MAP[['MR']] <- 'mR'
  MAP[['MS']] <- 'mS'
  MAP[['MK']] <- 'mK'
  MAP[['MB']] <- 'mB'
  MAP[['MD']] <- 'mD'
  MAP[['MV']] <- 'mV'
  MAP[['BG']] <- 'bG'
  MAP[['BR']] <- 'bR'
  MAP[['BS']] <- 'bS'
  MAP[['BK']] <- 'bK'
  MAP[['BB']] <- 'bB'
  MAP[['BD']] <- 'bD'
  MAP[['BV']] <- 'bV'
  MAP[['HG']] <- 'hG'
  MAP[['HR']] <- 'hR'
  MAP[['HS']] <- 'hS'
  MAP[['HK']] <- 'hK'
  MAP[['HB']] <- 'hB'
  MAP[['HD']] <- 'hD'
  MAP[['HV']] <- 'hV'
  MAP[['VG']] <- 'vG'
  MAP[['VR']] <- 'vR'
  MAP[['VS']] <- 'vS'
  MAP[['VK']] <- 'vK'
  MAP[['VB']] <- 'vB'
  MAP[['VD']] <- 'vD'
  MAP[['VV']] <- 'vV'
  
  MAP
}
MAP_M <- INIT_MAP_M()

INIT_MAP_D = function() {
  MAP <- NULL
  
  MAP[['CG']] <- 'yG'
  MAP[['CR']] <- 'yR'
  MAP[['CS']] <- 'yS'
  MAP[['CK']] <- 'yK'
  MAP[['CB']] <- 'yB'
  MAP[['CD']] <- 'yD'
  MAP[['CV']] <- 'yV'
  MAP[['YG']] <- 'yG'
  MAP[['YR']] <- 'yR'
  MAP[['YS']] <- 'yS'
  MAP[['YK']] <- 'yK'
  MAP[['YB']] <- 'yB'
  MAP[['YD']] <- 'yD'
  MAP[['YV']] <- 'yV'
  MAP[['SG']] <- 'bG'
  MAP[['SR']] <- 'bR'
  MAP[['SS']] <- 'bS'
  MAP[['SK']] <- 'bK'
  MAP[['SB']] <- 'bB'
  MAP[['SD']] <- 'bD'
  MAP[['SV']] <- 'bV'
  MAP[['MG']] <- 'hG'
  MAP[['MR']] <- 'hR'
  MAP[['MS']] <- 'hS'
  MAP[['MK']] <- 'hK'
  MAP[['MB']] <- 'hB'
  MAP[['MD']] <- 'hD'
  MAP[['MV']] <- 'hV'
  MAP[['BG']] <- 'bG'
  MAP[['BR']] <- 'bR'
  MAP[['BS']] <- 'bS'
  MAP[['BK']] <- 'bK'
  MAP[['BB']] <- 'bB'
  MAP[['BD']] <- 'bD'
  MAP[['BV']] <- 'bV'
  MAP[['HG']] <- 'hG'
  MAP[['HR']] <- 'hR'
  MAP[['HS']] <- 'hS'
  MAP[['HK']] <- 'hK'
  MAP[['HB']] <- 'hB'
  MAP[['HD']] <- 'hD'
  MAP[['HV']] <- 'hV'
  MAP[['VG']] <- 'nG'
  MAP[['VR']] <- 'nR'
  MAP[['VS']] <- 'nS'
  MAP[['VK']] <- 'nK'
  MAP[['VB']] <- 'nB'
  MAP[['VD']] <- 'nD'
  MAP[['VV']] <- 'nV'
  
  MAP
}
MAP_D <- INIT_MAP_D()



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Initialize Query Maps::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# These are for mapping Query positions to their possible outcomes based on
#  the assumption of unmethylated and methylated
#
# OLD_INIT_QMAP_* are previous versions that did not account for Infinium I
#  probes with some degenerate bases like W={A,T} that would still work
#

# TBD:: Address the tri-nucelotide version!

# Unmethylated will always use the letter earliest in the alphabet
#  e.g. 
#   A over T, 
#   A over C, 
#   A over G, 
#   C over G, 
#   C over T, 
#   G over T,
#
INIT_QMAP_U = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'  # '-'
  MAP[['G']] <- 'G'  # '-'
  MAP[['T']] <- 'T'
  
  MAP[['R']] <- 'A'
  MAP[['Y']] <- 'T'
  MAP[['S']] <- 'C'  # C over G 
  MAP[['W']] <- 'A'  # A over T
  # MAP[['S']] <- 'S'  # '-' Old Version
  # MAP[['W']] <- 'W'  # '-' Old Version
  
  MAP[['K']] <- 'T'
  MAP[['M']] <- 'A'
  MAP[['B']] <- 'T'
  MAP[['D']] <- 'W'
  MAP[['H']] <- 'W'
  MAP[['V']] <- 'A'
  
  MAP[['N']] <- 'W'
  
  # Lower Case::
  MAP[['a']] <- 'a'
  MAP[['c']] <- 'c'  # '-'
  MAP[['g']] <- 'g'  # '-'
  MAP[['t']] <- 't'
  
  MAP[['r']] <- 'a'
  MAP[['y']] <- 't'
  MAP[['s']] <- 'c'  # c over g
  MAP[['w']] <- 'a'  # a over t
  # MAP[['s']] <- 's'  # '-' Old Version
  # MAP[['w']] <- 'w'  # '-' Old Version
  
  MAP[['k']] <- 't'
  MAP[['m']] <- 'a'
  MAP[['b']] <- 't'
  MAP[['d']] <- 'w'
  MAP[['h']] <- 'w'
  MAP[['v']] <- 'a'
  
  MAP[['n']] <- 'w'
  
  MAP
}
QMAP_U <- INIT_QMAP_U()

# Methylated will always use the letter latest in the alphabet
#  e.g. 
#   T over A, 
#   T over C, 
#   T over G,
#   G over A, 
#   G over C, 
#   C over A, 
#
INIT_QMAP_M = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'  # '-'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'  # '-'
  
  MAP[['R']] <- 'G'
  MAP[['Y']] <- 'C'
  MAP[['S']] <- 'G'  # G over C
  MAP[['W']] <- 'T'  # T over A
  # MAP[['S']] <- 'S'  # '-' Old Version
  # MAP[['W']] <- 'W'  # '-' Old Version
  
  MAP[['K']] <- 'G'
  MAP[['M']] <- 'C'
  MAP[['B']] <- 'S'
  MAP[['D']] <- 'G'
  MAP[['H']] <- 'C'
  MAP[['V']] <- 'S'
  
  MAP[['N']] <- 'S'
  
  # Lower Case::
  MAP[['a']] <- 'a'  # '-'
  MAP[['c']] <- 'c'
  MAP[['g']] <- 'g'
  MAP[['t']] <- 't'  # '-'
  
  MAP[['r']] <- 'g'
  MAP[['y']] <- 'c'
  MAP[['s']] <- 'g'  # g over c
  MAP[['w']] <- 't'  # t over a
  # MAP[['s']] <- 's'  # '-' Old Version
  # MAP[['w']] <- 'w'  # '-' Old Version
  
  MAP[['k']] <- 'g'
  MAP[['m']] <- 'c'
  MAP[['b']] <- 's'
  MAP[['d']] <- 'g'
  MAP[['h']] <- 'c'
  MAP[['v']] <- 's'
  
  MAP[['n']] <- 's'
  
  MAP
}
QMAP_M <- INIT_QMAP_M()

# OLD Versions to be used now::
#
INIT_QMAP_U = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'  # '-'
  MAP[['G']] <- 'G'  # '-'
  MAP[['T']] <- 'T'
  
  MAP[['R']] <- 'A'
  MAP[['Y']] <- 'T'
  MAP[['S']] <- 'S'  # '-'
  MAP[['W']] <- 'W'  # '-'
  
  MAP[['K']] <- 'T'
  MAP[['M']] <- 'A'
  MAP[['B']] <- 'T'
  MAP[['D']] <- 'W'
  MAP[['H']] <- 'W'
  MAP[['V']] <- 'A'
  
  MAP[['N']] <- 'W'
  
  # Lower Case::
  MAP[['a']] <- 'a'
  MAP[['c']] <- 'c'  # '-'
  MAP[['g']] <- 'g'  # '-'
  MAP[['t']] <- 't'
  
  MAP[['r']] <- 'a'
  MAP[['y']] <- 't'
  MAP[['s']] <- 's'  # '-'
  MAP[['w']] <- 'w'  # '-'
  
  MAP[['k']] <- 't'
  MAP[['m']] <- 'a'
  MAP[['b']] <- 't'
  MAP[['d']] <- 'w'
  MAP[['h']] <- 'w'
  MAP[['v']] <- 'a'
  
  MAP[['n']] <- 'w'
  
  MAP
}
# OLD_QMAP_U <- OLD_INIT_QMAP_U()

OLD_INIT_QMAP_M = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'  # '-'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'  # '-'
  
  MAP[['R']] <- 'G'
  MAP[['Y']] <- 'C'
  MAP[['S']] <- 'S'  # '-'
  MAP[['W']] <- 'W'  # '-'
  
  MAP[['K']] <- 'G'
  MAP[['M']] <- 'C'
  MAP[['B']] <- 'S'
  MAP[['D']] <- 'G'
  MAP[['H']] <- 'C'
  MAP[['V']] <- 'S'
  
  MAP[['N']] <- 'S'
  
  # Lower Case::
  MAP[['a']] <- 'a'  # '-'
  MAP[['c']] <- 'c'
  MAP[['g']] <- 'g'
  MAP[['t']] <- 't'  # '-'
  
  MAP[['r']] <- 'g'
  MAP[['y']] <- 'c'
  MAP[['s']] <- 's'  # '-'
  MAP[['w']] <- 'w'  # '-'
  
  MAP[['k']] <- 'g'
  MAP[['m']] <- 'c'
  MAP[['b']] <- 's'
  MAP[['d']] <- 'g'
  MAP[['h']] <- 'c'
  MAP[['v']] <- 's'
  
  MAP[['n']] <- 's'
  
  MAP
}
# QMAP_M <- OLD_INIT_QMAP_M()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Initialize IUPAC Equivelency Tables::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# This is really for convience...
#
INIT_IUPAC_EQ = function() {
  MAP <- NULL
  
  MAP[['AA']] <- TRUE
  MAP[['AD']] <- TRUE
  MAP[['AH']] <- TRUE
  MAP[['AM']] <- TRUE
  MAP[['AN']] <- TRUE
  MAP[['AR']] <- TRUE
  MAP[['AV']] <- TRUE
  MAP[['AW']] <- TRUE
  MAP[['BB']] <- TRUE
  MAP[['BC']] <- TRUE
  MAP[['BG']] <- TRUE
  MAP[['BK']] <- TRUE
  MAP[['BN']] <- TRUE
  MAP[['BS']] <- TRUE
  MAP[['BT']] <- TRUE
  MAP[['BY']] <- TRUE
  MAP[['CB']] <- TRUE
  MAP[['CC']] <- TRUE
  MAP[['CH']] <- TRUE
  MAP[['CM']] <- TRUE
  MAP[['CN']] <- TRUE
  MAP[['CS']] <- TRUE
  MAP[['CV']] <- TRUE
  MAP[['CY']] <- TRUE
  MAP[['DA']] <- TRUE
  MAP[['DD']] <- TRUE
  MAP[['DG']] <- TRUE
  MAP[['DK']] <- TRUE
  MAP[['DN']] <- TRUE
  MAP[['DR']] <- TRUE
  MAP[['DT']] <- TRUE
  MAP[['DW']] <- TRUE
  MAP[['GB']] <- TRUE
  MAP[['GD']] <- TRUE
  MAP[['GG']] <- TRUE
  MAP[['GK']] <- TRUE
  MAP[['GN']] <- TRUE
  MAP[['GR']] <- TRUE
  MAP[['GS']] <- TRUE
  MAP[['GV']] <- TRUE
  MAP[['HA']] <- TRUE
  MAP[['HC']] <- TRUE
  MAP[['HH']] <- TRUE
  MAP[['HM']] <- TRUE
  MAP[['HN']] <- TRUE
  MAP[['HT']] <- TRUE
  MAP[['HW']] <- TRUE
  MAP[['HY']] <- TRUE
  MAP[['KB']] <- TRUE
  MAP[['KD']] <- TRUE
  MAP[['KG']] <- TRUE
  MAP[['KK']] <- TRUE
  MAP[['KN']] <- TRUE
  MAP[['KT']] <- TRUE
  MAP[['MA']] <- TRUE
  MAP[['MC']] <- TRUE
  MAP[['MH']] <- TRUE
  MAP[['MM']] <- TRUE
  MAP[['MN']] <- TRUE
  MAP[['MV']] <- TRUE
  MAP[['NA']] <- TRUE
  MAP[['NB']] <- TRUE
  MAP[['NC']] <- TRUE
  MAP[['ND']] <- TRUE
  MAP[['NG']] <- TRUE
  MAP[['NH']] <- TRUE
  MAP[['NK']] <- TRUE
  MAP[['NM']] <- TRUE
  MAP[['NN']] <- TRUE
  MAP[['NR']] <- TRUE
  MAP[['NS']] <- TRUE
  MAP[['NT']] <- TRUE
  MAP[['NV']] <- TRUE
  MAP[['NW']] <- TRUE
  MAP[['NY']] <- TRUE
  MAP[['RA']] <- TRUE
  MAP[['RD']] <- TRUE
  MAP[['RG']] <- TRUE
  MAP[['RN']] <- TRUE
  MAP[['RR']] <- TRUE
  MAP[['RV']] <- TRUE
  MAP[['SB']] <- TRUE
  MAP[['SC']] <- TRUE
  MAP[['SG']] <- TRUE
  MAP[['SN']] <- TRUE
  MAP[['SS']] <- TRUE
  MAP[['SV']] <- TRUE
  MAP[['TB']] <- TRUE
  MAP[['TD']] <- TRUE
  MAP[['TH']] <- TRUE
  MAP[['TK']] <- TRUE
  MAP[['TN']] <- TRUE
  MAP[['TT']] <- TRUE
  MAP[['TW']] <- TRUE
  MAP[['TY']] <- TRUE
  MAP[['VA']] <- TRUE
  MAP[['VC']] <- TRUE
  MAP[['VG']] <- TRUE
  MAP[['VM']] <- TRUE
  MAP[['VN']] <- TRUE
  MAP[['VR']] <- TRUE
  MAP[['VS']] <- TRUE
  MAP[['VV']] <- TRUE
  MAP[['WA']] <- TRUE
  MAP[['WD']] <- TRUE
  MAP[['WH']] <- TRUE
  MAP[['WN']] <- TRUE
  MAP[['WT']] <- TRUE
  MAP[['WW']] <- TRUE
  MAP[['YB']] <- TRUE
  MAP[['YC']] <- TRUE
  MAP[['YH']] <- TRUE
  MAP[['YN']] <- TRUE
  MAP[['YT']] <- TRUE
  MAP[['YY']] <- TRUE
  
  MAP
}
IUPAC_EQ <- INIT_IUPAC_EQ()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Temp Code to be deleted::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  ret_tib <- readr::read_tsv(file, # guess_max=1000000)
                             col_types = cols(
                               Seq_ID = col_character(),
                               Forward_Sequence = col_character(),
                               Genome_Build = col_character(),
                               Chromosome = col_character(),
                               Coordinate = col_double(),
                               Design_State = col_character(),
                               Seq_Length = col_double(),
                               Forward_CpG_Coord = col_double(),
                               TB_Strand = col_character(),
                               Top_Sequence = col_character(),
                               Top_CpG_Coord = col_double(),
                               Probe_Type = col_character(),
                               Probeset_ID = col_character(),
                               Probeset_Score = col_double(),
                               Methyl_Probe_ID = col_character(),
                               Methyl_Probe_Sequence = col_character(),
                               Methyl_Probe_Length = col_double(),
                               Methyl_Start_Coord = col_double(),
                               Methyl_End_Coord = col_double(),
                               Methyl_Probe_Covered_Top_Sequence = col_character(),
                               Methyl_Allele_FR_Strand = col_character(),
                               Methyl_Allele_TB_Strand = col_character(),
                               Methyl_Allele_CO_Strand = col_character(),
                               Methyl_Allele_Type = col_character(),
                               Methyl_Final_Score = col_double(),
                               Methyl_Tm = col_double(),
                               Methyl_Tm_Score = col_double(),
                               Methyl_GC_Percent = col_double(),
                               Methyl_GC_Score = col_double(),
                               Methyl_13mer_Count = col_double(),
                               Methyl_13mer_Score = col_double(),
                               Methyl_Address_Count = col_double(),
                               Methyl_Address_Score = col_double(),
                               Methyl_Self_Complementarity = col_double(),
                               Methyl_Self_Complementarity_Score = col_double(),
                               Methyl_Mono_Run = col_double(),
                               Methyl_Mono_Run_Score = col_double(),
                               Methyl_Ectopic_Count = col_double(),
                               Methyl_Ectopic_Score = col_double(),
                               Methyl_Underlying_CpG_Count = col_double(),
                               Methyl_Underlying_CpG_Min_Dist = col_double(),
                               Methyl_Underlying_CpG_Score = col_double(),
                               Methyl_In_CpG_Island_Relaxed = col_logical(),
                               Methyl_CpG_Island_Score = col_double(),
                               Methyl_Next_Base = col_character(),
                               Methyl_Next_Base_Score = col_double(),
                               UnMethyl_Probe_ID = col_character(),
                               UnMethyl_Probe_Sequence = col_character(),
                               UnMethyl_Probe_Length = col_double(),
                               UnMethyl_Start_Coord = col_double(),
                               UnMethyl_End_Coord = col_double(),
                               UnMethyl_Probe_Covered_Top_Sequence = col_character(),
                               UnMethyl_Allele_FR_Strand = col_character(),
                               UnMethyl_Allele_TB_Strand = col_character(),
                               UnMethyl_Allele_CO_Strand = col_character(),
                               UnMethyl_Allele_Type = col_character(),
                               UnMethyl_Final_Score = col_double(),
                               UnMethyl_Tm = col_double(),
                               UnMethyl_Tm_Score = col_double(),
                               UnMethyl_GC_Percent = col_double(),
                               UnMethyl_GC_Score = col_double(),
                               UnMethyl_13mer_Count = col_double(),
                               UnMethyl_13mer_Score = col_double(),
                               UnMethyl_Address_Count = col_double(),
                               UnMethyl_Address_Score = col_double(),
                               UnMethyl_Self_Complementarity = col_double(),
                               UnMethyl_Self_Complementarity_Score = col_double(),
                               UnMethyl_Mono_Run = col_double(),
                               UnMethyl_Mono_Run_Score = col_double(),
                               UnMethyl_Ectopic_Count = col_double(),
                               UnMethyl_Ectopic_Score = col_double(),
                               UnMethyl_Underlying_CpG_Count = col_double(),
                               UnMethyl_Underlying_CpG_Min_Dist = col_double(),
                               UnMethyl_Underlying_CpG_Score = col_double(),
                               UnMethyl_In_CpG_Island_Relaxed = col_logical(),
                               UnMethyl_CpG_Island_Score = col_double(),
                               UnMethyl_Next_Base = col_character(),
                               UnMethyl_Next_Base_Score = col_double()
                             )
  )
}


# End of file

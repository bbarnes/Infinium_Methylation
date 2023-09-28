
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
#             c-improbe = c++ improbe (traditional) via docker image
#                    Includes Thermodynamic Calculations
#                   Only Designs Infinium I U/M Probes 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

c_improbe_workflow = function(template_tib,
                              
                              # Add genome source row::
                              imGenome = NULL,
                              
                              ids_key = "Prb_Key_Unq",
                              fwd_seq = "Forward_Sequence",
                              pos_key = "Bsp_Pos",
                              chr_key = "Bsp_Chr",
                              gen_key = NULL,
                              
                              doc_image,
                              doc_shell,
                              
                              outlevel = 0,
                              call_inf = FALSE,
                              
                              re_join  = FALSE,
                              new_join = NULL,
                              old_join = NULL,
                              
                              reload   = FALSE,
                              
                              out_dir,
                              run_tag,
                              re_load = FALSE,
                              pre_tag = NULL,
                              end_str = 'tsv.gz',
                              sep_chr = '.',
                              
                              verbose=0, vt=3,tc=1,tt=NULL,
                              funcTag='c_improbe_workflow') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  out_csv <- NULL
  out_csv <- redata(out_dir, run_tag, funcTag, re_load, 
                    pre_tag, end_str=end_str, sep=sep_chr,
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  if (tibble::is_tibble(out_csv)) return(out_csv)
  if (is.null(out_csv)) {
    stop(glue::glue("{RET}{mssg} ERROR: out_csv is NULL!{RET2}"))
    return(out_csv)
  }
  out_dir <- base::dirname(out_csv)
  
  imp_tsv <- out_csv %>% 
    stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
    paste("improbe-inputs.tsv.gz", sep=sep_chr)
  
  imp_des_tsv <- out_csv %>% 
    stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
    paste("improbe-designOutput.tsv.gz", sep=sep_chr)
  
  if (is.null(gen_key))
    if (!is.null(imGenome_dat) && !is.null(imGenome_dat$Genome_Version))
      gen_key <- imGenome_dat$Genome_Version
  else gen_key <- "x"
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}       out_csv={out_csv}.{RET}"))
    cat(glue::glue("{mssg}       imp_tsv={imp_tsv}.{RET}"))
    cat(glue::glue("{mssg}   imp_des_tsv={imp_des_tsv}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}      ids_key={ids_key}.{RET}"))
    cat(glue::glue("{mssg}      fwd_seq={fwd_seq}.{RET}"))
    cat(glue::glue("{mssg}      pos_key={pos_key}.{RET}"))
    cat(glue::glue("{mssg}      chr_key={chr_key}.{RET}"))
    cat(glue::glue("{mssg}      gen_key={gen_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}       outlevel={outlevel}.{RET}"))
    cat(glue::glue("{mssg}     call_inf={call_inf}.{RET}"))
    cat(glue::glue("{mssg}        re_join={re_join}.{RET}"))
    cat(glue::glue("{mssg}      reload={reload}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Docker Parameters::{RET}"))
    cat(glue::glue("{mssg}     doc_image={doc_image}.{RET}"))
    cat(glue::glue("{mssg}     doc_shell={doc_shell}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #               1.0 Build Input:: Forward Template Sequence
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (!is.null(imp_tsv) && file.exists(imp_tsv)) {
      if (is.null(out_dir) || !dir.exists(out_dir)) {
        out_dir = base::dirname(imp_tsv)
        dir.create(out_dir, recursive = TRUE)
      }
    } else {
      if (!dir.exists(out_dir)) dir.create(out_dir)
      
      imp_col <- c("Seq_ID","Sequence","Genome_Build",
                   "Chromosome","Coordinate","CpG_Island")
      
      # imp_tsv <- file.path(out_dir, paste(prefix,inp_suffix,"tsv.gz", sep='.'))
      template_tib <- template_tib %>% 
        # dplyr::select(-Genome_Build, -Strand_CO, -Strand_FR) %>%
        # dplyr::select(-dplyr::any_of("Genome_Build","Strand_CO","Strand_FR")) %>%
        dplyr::mutate(Genome_Build=!!gen_key, CpG_Island="FALSE") %>%
        dplyr::select(dplyr::all_of(c(!!ids_key, !!fwd_seq, "Genome_Build", 
                                      !!chr_key, !!pos_key, "CpG_Island"))) %>%
        purrr::set_names(imp_col)
      
      cat(glue::glue("{mssg} imp_tsv={imp_tsv}; template_tib::{RET}"))
      template_tib %>% head(n=3) %>% print()
      
      out_cnt <- safe_write(x=template_tib, file=imp_tsv, funcTag=funcTag, 
                            verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                 1.1 Run improbe designs:: via docker
    #                     (improbe = original c++ version)
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # TBD:: Make Docker Options more descriptive...
    ret_val <- run_improbe_docker(file   = imp_tsv, 
                                  name   = paste(run_tag,funcTag, sep=sep_chr),
                                  image  = doc_image, 
                                  shell  = doc_shell,
                                  reload = reload,
                                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     1.2 Load improbe designs:: filter
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- load_improbe_design(des_tsv  = imp_des_tsv, 
                                   out_tsv  = out_csv,
                                   imGenome = imGenome,
                                   
                                   outlevel = outlevel,
                                   call_inf  = call_inf, 
                                   
                                   re_join  = re_join,
                                   new_join = new_join,
                                   old_join = old_join,
                                   
                                   verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Docker improbe Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_improbe_docker = function(file,
                              out  = NULL, 
                              name = "unk", 
                              image = "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos.v.1.15",
                              shell = "run_improbe.sh",
                              suffix = 'improbe-design',
                              reload = FALSE,
                              
                              verbose=1, vt=3,tc=1,tt=NULL,
                              funcTag = 'run_improbe_docker') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}      out={out}.{RET}"))
    cat(glue::glue("{mssg}     file={file}.{RET}"))
    cat(glue::glue("{mssg}     name={name}.{RET}"))
    cat(glue::glue("{mssg}    image={image}.{RET}"))
    cat(glue::glue("{mssg}    shell={shell}.{RET}"))
    cat(glue::glue("{mssg}   reload={reload}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  out_dir <- out
  if (is.null(out)) out_dir <- base::dirname(file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
  
  cat(glue::glue("{mssg}      out={out}.{RET}"))
  
  ret_log <- file.path(out_dir, paste(name,suffix,'log', sep='.'))
  ret_tsv <- file.path(out_dir, paste(name,suffix,'tsv.gz', sep='.'))
  base_file <- base::basename( file )
  
  cat(glue::glue("{mssg}         file={file}.{RET}"))
  cat(glue::glue("{mssg}    base_file={base_file}.{RET}"))
  cat(glue::glue("{mssg}      ret_log={ret_log}.{RET}"))
  cat(glue::glue("{mssg}      ret_tsv={ret_tsv}.{RET}"))
  
  if (FALSE) {
    if (reload && 
        file.exists(file) &&
        file.exists(ret_log) &&
        file.exists(ret_tsv) &&
        file.mtime(file) <= file.mtime(ret_log) &&
        file.mtime(file) <= file.mtime(ret_tsv)) {
      if (verbose>=vt) 
        cat(glue::glue("{mssg} Reloading...{RET}"))
      
      ret_cnt <- 0
      etime   <- 0
      if (!is.null(tt)) tt$addTime(stime,funcTag)
      if (verbose>=vt) cat(glue::glue(
        "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
        "{RET}{tabs}{BRK}{RET2}"))
      
      return(ret_cnt)
    }
  }
  
  ret_cnt <- 1
  stime <- base::system.time({
    
    if (file.exists(ret_log)) unlink(ret_log)
    if (file.exists(ret_tsv)) unlink(ret_tsv)
    
    base::system(glue::glue("touch {ret_log}"))
    base::system(glue::glue("touch {ret_tsv}"))
    
    # TBD:: Update Docker with more options::
    imp_doc_cmd <- glue::glue("docker run -i --rm ",
                              "-v {out_dir}:/input ",
                              "-v {out_dir}:/output ",
                              "-w /work ",
                              "{image} {shell} {base_file} {name}")
    
    if (verbose>=vt)
      cat(glue::glue("{mssg} Running improbe cmd='{imp_doc_cmd}'...{RET}"))
    ret_cnt <- base::system(imp_doc_cmd)
    
    if (ret_cnt != 0) {
      stop(glue::glue("{RET}{mssg} ERROR: cmd return={ret_cnt} ",
                      "cmd='{imp_doc_cmd}'!{RET2}"))
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_cnt
}

load_improbe_design = function(des_tsv, 
                               out_tsv, 
                               imGenome = NULL,
                               
                               outlevel = 0,
                               re_join = FALSE, 
                               new_join = c("Aln_Key_Unq","Bsp_Chr",
                                            "Bsp_Pos","Bsp_FR","Bsp_CO"),
                               old_join = c("Seq_ID","Chromosome","Coordinate",
                                            "Strand_FR","Strand_CO"),
                               call_inf  = TRUE,
                               
                               verbose=0, vt=3,tc=1,tt=NULL, 
                               funcTag='load_improbe_design') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}    des_tsv={des_tsv}.{RET}"))
    cat(glue::glue("{mssg}    out_tsv={out_tsv}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Run Parameters::{RET}"))
    cat(glue::glue("{mssg}    re_join={re_join}.{RET}"))
    cat(glue::glue("{mssg}    call_inf={call_inf}.{RET}"))
    cat(glue::glue("{mssg}   outlevel={outlevel}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  # improbe output original fields::
  #
  org_cols <- NULL
  new_cols <- NULL
  if (outlevel>0) {
    # outlevel 1::
    ids_org_cols <- c("Seq_ID", "Genome_Build", "Chromosome", "Coordinate")
    seq_org_cols <- c("Forward_Sequence", "Top_Sequence")
    prb_org_cols <- c("Methyl_Probe_Covered_Top_Sequence", "UnMethyl_Probe_Sequence", "Methyl_Probe_Sequence")
    srd_org_cols <- c("Methyl_Allele_FR_Strand", "Methyl_Allele_TB_Strand", "Methyl_Allele_CO_Strand", "Methyl_Next_Base")
    
    # outlevel 2:: Final Scores::
    scrF_org_cols <- c("UnMethyl_Final_Score", "Methyl_Final_Score")
    
    # outlevel 3:: Both Metric/Score::
    rawB_org_cols <- c("Methyl_Underlying_CpG_Count", "Methyl_Underlying_CpG_Min_Dist")
    scrB_org_cols <- c("Methyl_Underlying_CpG_Score", "Methyl_Next_Base_Score")
    
    # outlevel 4:: Metrics
    rawM_org_cols <- c("Methyl_Tm",    "Methyl_GC_Percent", "Methyl_13mer_Count", "Methyl_Address_Count",
                       "Methyl_Self_Complementarity",       "Methyl_Mono_Run",    "Methyl_Ectopic_Count")
    rawU_org_cols <- paste0("Un",rawM_org_cols)
    
    # outlevel 5:: Scores
    scrM_org_cols <- c("Methyl_Tm_Score","Methyl_GC_Score",   "Methyl_13mer_Score",  "Methyl_Address_Score",
                       "Methyl_Self_Complementarity_Score", "Methyl_Mono_Run_Score", "Methyl_Ectopic_Score")
    scrU_org_cols <- paste0("Un",scrM_org_cols)
    
    # improbe output new fields::
    #
    ids_new_cols <- c("Seq_ID", "Genome_Build", "Chromosome", "Coordinate")
    seq_new_cols <- c("Forward_Sequence", "Top_Sequence")
    prb_new_cols <- c("Probe_Seq_T", "Probe_Seq_U", "Probe_Seq_M")
    srd_new_cols <- c("Strand_FR", "Strand_TB", "Strand_CO", "Next_Base")
    
    # outlevel 2:: Final Scores
    scrF_new_cols <- c("Scr_U", "Scr_M")
    
    # outlevel 3: Both Metric/Score::
    rawB_new_cols <- c("Cpg_Cnt", "Cpg_Dis")
    scrB_new_cols <- c("Cpg_Scr", "Nxb_Scr")
    
    # outlevel 4:: New Metrics::
    met_new_cols  <- c("Tm", "GC", "Mer13", "Address", "SelfCmpl", "Mono", "Ectopic")
    rawU_new_cols <- paste(met_new_cols, "Raw_U", sep="_")
    rawM_new_cols <- paste(met_new_cols, "Raw_M", sep="_")
    
    # outlevel 5:: New Metrics::
    scrU_new_cols <- paste(met_new_cols, "Scr_U", sep="_")
    scrM_new_cols <- paste(met_new_cols, "Scr_M", sep="_")
    
    #
    # Old/New Field Names to Select/Return::
    #
    if (outlevel>=1) {
      org_cols <- c(org_cols, ids_org_cols,seq_org_cols,prb_org_cols,srd_org_cols)
      new_cols <- c(new_cols, ids_new_cols,seq_new_cols,prb_new_cols,srd_new_cols)
    }
    if (outlevel>=2) {
      org_cols <- c(org_cols, scrF_org_cols)
      new_cols <- c(new_cols, scrF_new_cols)
    }
    if (outlevel>=3) {
      org_cols <- c(org_cols, rawB_org_cols,scrB_org_cols)
      new_cols <- c(new_cols, rawB_new_cols,scrB_new_cols)
    }
    if (outlevel>=4) {
      org_cols <- c(org_cols, rawU_org_cols,rawM_org_cols)
      new_cols <- c(new_cols, rawU_new_cols,rawM_new_cols)
    }
    if (outlevel>=5) {
      org_cols <- c(org_cols, scrU_org_cols,scrM_org_cols)
      new_cols <- c(new_cols, scrU_new_cols,scrM_new_cols)
    }
    
    if (verbose>=vt+6) {
      org_len <- org_cols %>% length()
      cat(glue::glue("{mssg} org_cols(org_len)={RET}"))
      print(org_cols)
      
      new_len <- new_cols %>% length()
      cat(glue::glue("{mssg} new_cols(new_len)={RET}"))
      print(new_cols)
    }
  } else {
    cat(glue::glue("{mssg} Returning full data (outlevel={outlevel})...{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (verbose>=vt) 
      cat(glue::glue("{mssg} Loading designs={des_tsv}...{RET}"))
    
    ret_tib <- 
      suppressMessages(suppressWarnings( readr::read_tsv(des_tsv) ))
    ret_key <- glue::glue("s-improbe-Raw-Improbe-Data-tib({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
    
    # Now subset and rename based on input out-level::
    #
    if (!is.null(org_cols) && !is.null(new_cols)) {
      if (verbose>=vt+2)
        cat(glue::glue("{mssg} Subseting and renaming based on outlevel={outlevel} ",
                       "input...{RET}"))
      
      ret_tib <- ret_tib %>% 
        dplyr::select(dplyr::all_of(org_cols)) %>%
        dplyr::mutate(Methyl_Allele_TB_Strand=
                        stringr::str_sub(Methyl_Allele_TB_Strand, 1,1)) %>%
        purrr::set_names(new_cols)
    }
    ret_tib <- ret_tib %>% clean_tibble()
    ret_key <- glue::glue("s-improbe-subset/renaming-tib({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
    
    if (call_inf && outlevel>=3) {
      if (verbose>=vt+2)
        cat(glue::glue("{mssg} Adding Infinium Preference outlevel={outlevel}.{RET}"))
      
      # Adding some extra required fields for matching later::
      #
      ret_tib <- ret_tib %>%
        dplyr::rename(Strand_Ref_FR=Strand_FR) %>%
        dplyr::mutate(
          Scr_Min=pmin(Scr_U,Scr_M),
          Inf_Type=dplyr::case_when(
            Scr_Min<0.2                  ~ 0,
            Scr_Min<0.3 & Strand_CO=="C" ~ 0,
            
            Scr_Min< 0.3 & Strand_CO=="O" & Cpg_Cnt==0 ~ 1,
            Scr_Min< 0.4 & Strand_CO=="O" & Cpg_Cnt==1 ~ 1,
            Scr_Min< 0.5 & Strand_CO=="O" & Cpg_Cnt==2 ~ 1,
            Scr_Min< 0.6 & Strand_CO=="O" & Cpg_Cnt==3 ~ 1,
            
            Scr_Min< 0.4 & Strand_CO=="C" & Cpg_Cnt==0 ~ 1,
            Scr_Min< 0.5 & Strand_CO=="C" & Cpg_Cnt==1 ~ 1,
            Scr_Min< 0.6 & Strand_CO=="C" & Cpg_Cnt==2 ~ 1,
            Scr_Min< 0.7 & Strand_CO=="C" & Cpg_Cnt==3 ~ 1,
            
            Scr_Min>=0.7 ~ 1,
            
            TRUE ~ 2
          ),
          Imp_U49=dplyr::case_when(
            Strand_CO=="C" ~ stringr::str_sub(Probe_Seq_U,1,49),
            Strand_CO=="O" ~ stringr::str_sub(Probe_Seq_U,2,50),
            TRUE ~ NA_character_
          ),
          Imp_M49=dplyr::case_when(
            Strand_CO=="C" ~ stringr::str_sub(Probe_Seq_M,1,49),
            Strand_CO=="O" ~ stringr::str_sub(Probe_Seq_M,2,50),
            TRUE ~ NA_character_
          ),
          Strand_FR=dplyr::case_when(
            Strand_Ref_FR=="F" ~ "R",
            Strand_Ref_FR=="R" ~ "F",
            TRUE ~ NA_character_
          )
        )
      ret_key <- glue::glue("s-improbe-after-Infinium-Preference({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
    }
    
    # NOTE: The more I think of it we should avoid this process of rejoining
    #   all the data together after each step...
    #
    # Merge data back with BSP results::
    # if (re_join) {
    #   if (verbose>=vt+2)
    #     cat(glue::glue("{mssg} Joining data with improbe data...{RET}"))
    #   
    #   ret_tib <- dplyr::left_join(
    #     join, dplyr::rename_with(ret_tib, ~ new_join, dplyr::all_of(old_join) ),
    #     by=new_join, suffix=c("_bsp","_imp")
    #   )
    #   ret_key <- glue::glue("s-improbe-merging-back-with-input-tib({funcTag})")
    #   ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
    # }
    
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Re-extracting Ord_Des/Ord_Din from Seq_ID...{RET}"))
    
    ret_tib <- ret_tib %>% dplyr::mutate(
      Info=Seq_ID %>% 
        stringr::str_remove("^[^_]+_") %>% 
        stringr::str_remove("_.*$"),
      improbe_type="c") %>% 
      tidyr::separate(Info, into=c("Ord_Des","Ord_Din"), sep=c(1), remove=TRUE)
    
    ret_key <- glue::glue("s-improbe-After-Des/Din-Extraction({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
    
    if (!is.null(imGenome)) ret_tib <- imGenome %>%
      dplyr::select(-Path, -dplyr::ends_with("_Cnt"), 
                    -dplyr::ends_with("_Int")) %>% 
      cbind(ret_tib) %>%
      tibble::as_tibble()
    
    out_cnt <- safe_write(x=ret_tib, file=out_tsv, funcTag=funcTag,
                          verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)

    ret_key <- glue::glue("Clean-Improbe-Data({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file

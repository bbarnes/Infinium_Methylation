
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Load Core Packages::
suppressWarnings(suppressPackageStartupMessages( base::require("sesame",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("minfi",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse",quietly=TRUE) ))

# Load Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel",quietly=TRUE) ))

COM <- ","
TAB <- "\t"
RET <- "\n"
BNG <- "|"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Calls Summary Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

callToPassPerc = function(tib=NULL, file=NULL, key, name=NULL, idx=NULL,
                          min=0.1, type='cg',
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'callToPassPerc'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting;  key={key}..{RET}"))

  ret_cnt <- 0
  ret_tib <- NULL
  if (is.null(tib) && is.null(file)) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Must provide tib or file!{RET}{RET}"))
    return(ret_tib)
  }
  
  stime <- system.time({
    
    if (is.null(tib)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} file={file}{RET}"))
      
      tib <- suppressMessages(suppressWarnings( readr::read_csv(file)))
    }
    
    min_int <- as.numeric(min)
    key_sys <- rlang::sym(key)
    tot_cnt <- tib %>% dplyr::filter(stringr::str_starts(Probe_ID, !!type)) %>% base::nrow()
    min_cnt <- tib %>% dplyr::filter(stringr::str_starts(Probe_ID, !!type)) %>% 
      dplyr::filter(!!key_sys <= min) %>% base::nrow()
    pas_cnt <- tib %>% dplyr::filter(stringr::str_starts(Probe_ID, !!type)) %>% 
      dplyr::filter(!!key_sys <= min_int) %>% base::nrow()
    pas_per  <- round(100*pas_cnt/tot_cnt, 3)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} per={pas_per}, pas={pas_cnt}, tot={tot_cnt}{RET}"))
    
    ret_tib <- tibble(
      pass_perc=pas_per,
      mins_count=min_cnt,
      total_count=tot_cnt,
      min_cutoff=min,
      metric=key
    )
    if (!is.null(name)) ret_tib <- ret_tib %>% dplyr::mutate(name=!!name)
    if (!is.null(idx)) {
      ret_tib <- ret_tib %>% purrr::set_names(paste(type,'calls',names(.),idx, sep='_'))
    } else {
      ret_tib <- ret_tib %>% purrr::set_names(paste(type,'calls',names(.), sep='_'))
    }
      
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_tib={RET}"))
      print(ret_tib)
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame SSET Access Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToPassPercSsheet = function(sset, man=NULL, min, per, idx=0, type='cg',
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPassPercSsheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}  min={min}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}  per={per}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}  idx={idx}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} type={type}{RET}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_tab <- NULL
  stime <- system.time({
    
    min <- as.numeric(min)
    
    # Open Sesame Comparison::
    pval_tib <- sset@pval %>% 
      tibble::enframe(name="Probe_ID", value="pvals")
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval_tib={RET}"))
      print(pval_tib)
    }
    
    # TBD:: Manifest loading function needs to retain original IDs
    #
    if (is.null(man)) {
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Manifest present (type={type}); pval_tib={RET}"))
        pval_tib %>% 
          dplyr::filter(stringr::str_starts(Probe_ID,type)) %>% print()
      }

      ret_tab <- pval_tib %>% 
        dplyr::filter(stringr::str_starts(Probe_ID,type)) %>% 
        dplyr::summarise(total_cnt=n(),
                         pass_cnt=base::sum(pvals<=min, na.rm=TRUE),
                         pass_perc=round(100*pass_cnt/total_cnt, 3),
                         .groups='drop')
      
      req_cnt <- ret_tab %>% 
        dplyr::filter(pass_perc<per) %>%
        base::nrow()
      
      ret_tib <- ret_tab %>%
        dplyr::mutate(Failed_QC_cnt=req_cnt,
                      Failed_QC=dplyr::case_when(
                        Failed_QC_cnt>0 ~ TRUE, TRUE ~ FALSE)) %>%
        purrr::set_names(paste(type,names(.), sep='_'))

    } else {
      ret_tab <- pval_tib %>% 
        dplyr::left_join(man, by="Probe_ID") %>% 
        dplyr::filter(Probe_Type==type) %>%
        dplyr::group_by(Probe_Type,Probe_Design) %>% 
        dplyr::summarise(total_cnt=n(),
                         pass_cnt=base::sum(pvals<=min, na.rm=TRUE),
                         pass_perc=round(100*pass_cnt/total_cnt, 3),
                         .groups='drop')

      tot_cnt_tib <- ret_tab %>% 
        sigsSumToSSheet(metric="total_cnt")
      
      pas_cnt_tib <- ret_tab %>% 
        sigsSumToSSheet(metric="pass_cnt")
      
      pas_per_tib <- ret_tab %>% 
        sigsSumToSSheet(metric="pass_perc")
      
      ret_tib <- dplyr::bind_cols(tot_cnt_tib,pas_cnt_tib,pas_per_tib)
    }
    
    ret_tib <- ret_tib %>%
      purrr::set_names(paste(names(.),'basic',idx, sep='_'))
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_tib={RET}"))
      print(ret_tib)
      
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_tab={RET}"))
      print(ret_tab)
    }

    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Sesame SSET Prediction Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToPredictions = function(sset, del='_', platform=NULL, fresh=FALSE,
                             quality.mask = FALSE, sum.TypeI = FALSE,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPredictions'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (is.null(platform)) platform <- sset@platform
    
    if (verbose>=vt+1)
      cat(glue::glue("[{funcTag}]:{tabsStr} Using platform={platform} for Inference/Prediction calls.{RET}{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Sesame Channel Swap Methods:: Extra(sset)
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    igr_cnt <- 0
    irg_cnt <- 0
    igr_key <- 'SwapCntIGR'
    irg_key <- 'SwapCntIRG'
    if (!is.null(sesame::extra(sset)[['IRR']]) && 
        !is.null(sesame::extra(sset)[['IGG']])) {
      igr_cnt <- length(which(!sesame::extra(sset)[['IGG']]))
      irg_cnt <- length(which(!sesame::extra(sset)[['IRR']]))
    }
    ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!igr_key := !!igr_cnt))
    ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!irg_key := !!irg_cnt))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Sesame Inference Methods:: SSET
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (platform=='EPIC' || platform=='HM450') {
      gct_val <- NULL
      gct_key <- 'GCT'
      gct_val <- safeGCT(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!gct_key := !!gct_val))
      
      sex_val <- NULL
      sex_key <- 'Sex'
      sex_val <- safeSex(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!sex_key := !!sex_val))
      
      kar_val <- NULL
      kar_key <- 'SexKaryotype'
      kar_val <- safeSexKaryo(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!kar_key := !!kar_val))
      
      eth_val <- NULL
      eth_key <- 'Ethnicity'
      eth_val <- safeEthnicity(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!eth_key := !!eth_val))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      Sesame Predictions Methods:: Beta
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (platform=='EPIC' || platform=='HM450') {
      if (fresh || is.null(sesame::extra(sset)[['betas']]) ) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr} Generating beta calls...{RET}"))
        
        # sesame::extra(sset)[['betas']] <-
        #   sesame::getBetas(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
        sesame::extra(sset)[['betas']] <-
          getBetas2(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
      }
      
      if (platform=='EPIC') {
        skn_val <- NULL
        skn_key <- 'AgeSkinBlood'
        skn_val <- safeSkinAge(beta=sesame::extra(sset)[['betas']], 
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!skn_key := !!skn_val))
      }
      
      if (platform=='HM450') {
        phn_val <- NULL
        phn_key <- 'AgePheno'
        phn_val <- safePhenoAge(beta=sesame::extra(sset)[['betas']], 
                                verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!phn_key := !!phn_val))
      }
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Sesame SNP Methods:: VCF
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

safeVCF = function(sset, vcf, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeVCF'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; vcf={vcf}{RET}"))
  
  ret_val <- NA
  stime <- system.time({
    
    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressMessages(suppressWarnings(sesame::formatVCF(sset=sset, vcf=vcf) ))
      if (verbose>=vt) cat("\n")
      
      #
      # TBD:: Use bgzip and tabix to compress VCF
      #   - Need to install samtools, bgzip, tabix on Centos image first...
      #
      # cmd <- ""
      
      try_str
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Sesame Inference Methods:: SSET
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

safeGCT = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeGCT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({
    
    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::bisConversionControl(sset) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

safeEthnicity = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeEthnicity'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({
    
    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(inferEthnicity2(sset) )
      # suppressWarnings(sesame::inferEthnicity(sset) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

safeSexKaryo = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSexKaryo'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({
    
    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(inferSexKaryotypes2(sset) )
      # suppressWarnings(sesame::inferSexKaryotypes(sset) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

safeSex = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSex'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({
    
    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(inferSex2(sset) )
      # suppressWarnings(sesame::inferSex(sset) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame Predictions Methods:: Beta
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

safeSkinAge = function(beta, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSkinAge'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({
    
    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::predictAgeSkinBlood(betas=beta) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

safePhenoAge = function(beta, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safePhenoAge'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({
    
    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::predictAgePheno(betas=beta) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

expand_ids = function(tib,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'expand_ids'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib %>%
      dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>%
      tidyr::separate(Probe_ID, into=c("Seq_ID", "Seq_Str"), 
                      sep='_', remove=FALSE) %>% 
      tidyr::separate(Seq_Str, into=c("Blank","Strand_TB", "Strand_CO", "Design_Type", "Rep_Num"), 
                      sep="", convert=TRUE) %>% dplyr::select(-Blank)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Sesame Sample Sheet Methods:: Requeue
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

requeueFlagOpenSesame = function(tib, name, csv=NULL, idx=1,
                                 verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'requeueFlagOpenSesame'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({

    ret_tib <- tibble::tibble(
      Sentrix_Name  = name,
      Failed_QC     = tib$cg_Failed_QC_basic_0,
      Min_Pass_Perc = tib$cg_pass_perc_basic_0
    )

    if (!is.null(csv)) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing CSV={csv}.{RET}"))
      csv_dir <- base::dirname(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      readr::write_csv(ret_tib, csv)
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

requeueFlag = function(tib, name, csv=NULL, idx=1,
                       verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'requeueFlag'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    requeue_cnt <- 0
    requeue_cnt <- tib %>%
      dplyr::filter(Workflow_idx==idx & Probe_Type=='cg' & !is.na(Requeue)) %>%
      dplyr::filter(Requeue==TRUE) %>%
      base::nrow()
    
    min_pass_perc <- tib %>% 
      dplyr::filter(Workflow_idx==idx & Probe_Type=='cg' & !is.na(Requeue)) %>%
      dplyr::arrange(full_pass_perc) %>%
      head(n=1) %>% pull(full_pass_perc)
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} requeue_cnt={requeue_cnt}; ",
                     "min_pass_perc={min_pass_perc}; sel={RET}"))
      tib %>% 
        dplyr::filter(Workflow_idx==idx & Probe_Type=='cg' & !is.na(Requeue)) %>%
        print()
    }
    
    requeue_flag <- FALSE
    if (requeue_cnt>0) requeue_flag <- TRUE
    
    ret_tib <- tibble::tibble(Sentrix_Name = !!name,
                              Failed_QC = !!requeue_flag,
                              Min_Pass_Perc = !!min_pass_perc)
    
    if (!is.null(csv)) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing CSV={csv}.{RET}"))
      csv_dir <- base::dirname(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      readr::write_csv(ret_tib, csv)
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Sesame SSET Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sigsSumToSSheet = function(tib, metric=NULL,
                           by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsSumToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting (metric={metric})...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    by_sym   <- by %>% rlang::sym()
    des_sym  <- des %>% rlang::sym()
    type_sym <- type %>% rlang::sym()
    
    ret_tib <- tib %>% 
      tidyr::unite(!!type_sym, !!type_sym,!!des_sym, sep='_')
    
    if (!is.null(metric)) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Selecting only metric={metric}...{RET}"))
      ret_tib <- ret_tib %>% dplyr::select(!!type_sym, ends_with(!!metric))
    }
    
    ret_tib <- ret_tib %>%
      tidyr::gather(Metric, Value, -!!type_sym) %>%
      tidyr::unite(Key, !!type_sym,Metric, sep='_') %>%
      tidyr::spread(Key, Value)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame Tibs Summary Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetTibToSummary = function(tib, man=NULL,
                            pval=NULL, perc=NULL, percision=-1,
                            save=FALSE, csv=NULL,
                            by="Probe_ID", type="Probe_Type", des="Probe_Design",
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetTibToSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; by={by}, type={type}, des={des}...{RET}"))
  
  if (verbose>=vt) {
    if (!is.null(pval)) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval={pval}...{RET}"))
    if (!is.null(perc)) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} perc={perc}...{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    by_sym  <- by %>% rlang::sym()
    des_sym  <- des %>% rlang::sym()
    type_sym <- type %>% rlang::sym()
    
    if (!is.null(pval)) pval <- as.numeric(pval)
    if (!is.null(pval)) perc <- as.numeric(perc)
    
    if (!is.null(man)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Will use provided manifest...{RET}"))
      if (verbose>=vt+4) man %>% print()
      
      tib_cols <- names(tib)
      if (!(des %in% tib_cols)) {
        if (verbose>=vt+5)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} right_join3(by,type,des)...{RET}"))
        
        tib <- dplyr::select(man, !!by, !!type, !!des) %>% dplyr::right_join(tib, by=by)
      } else {
        if (verbose>=vt+5)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} right_join2(by,type)...{RET}"))
        
        tib <- dplyr::select(man, !!by, !!type) %>% dplyr::right_join(tib, by=by)
      }
      
    } else {
      if (verbose>=vt+1) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} No manifest provided...{RET}"))
    }
    tib_cols <- names(tib)
    if (!(type %in% tib_cols)) 
      tib <- tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(!!by_sym, 1,2))
    
    # Check for all fields::
    #
    tib_cols <- names(tib)
    if (!(by %in% tib_cols)) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr}{TAB} ERROR: Failed to find by={by}!!!{RET}{RET}"))
      return(ret_tib)
    }
    if (!(des %in% tib_cols)) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr}{TAB} ERROR: Failed to find des={des}!!!{RET}{RET}"))
      return(ret_tib)
    }
    if (!(type %in% tib_cols)) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr}{TAB} ERROR: Failed to find type={type}!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    ret_tib <- tib %>% 
      dplyr::select(-dplyr::all_of(!!by) ) %>%
      tidyr::gather(Metric, value, -dplyr::all_of(c(!!type, !!des)) ) %>%
      dplyr::group_by(!!type_sym, !!des_sym, Metric) %>%
      summarise_if(is.numeric, list(min=min, median=median, mean=mean, sd=sd, max=max), na.rm=TRUE) %>%
      dplyr::ungroup()
    
    if (!is.null(pval)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Applying pval cutoff={pval}...{RET}"))
      if (verbose>=vt+4) print(tib)
      
      cut2_tib <- tib %>% 
        dplyr::group_by(!!type_sym, !!des_sym) %>%
        dplyr::select(-dplyr::all_of(by)) %>%
        dplyr::summarise_all(list(pass_perc=cntPer_lte), min=pval)
      
      cut1_tib <- tib %>% 
        dplyr::select(-dplyr::all_of(des_sym)) %>%
        dplyr::group_by(!!type_sym) %>%
        dplyr::select(-dplyr::all_of(by)) %>%
        dplyr::summarise_all(list(full_pass_perc=cntPer_lte), min=pval)

      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} cut2_tib={RET}"))
        print(cut2_tib)
        cat(glue::glue("[{funcTag}]:{tabsStr} cut1_tib={RET}"))
        print(cut1_tib)
      }
      
      if (!is.null(perc)) {
        cut1_tib <- cut1_tib %>% 
        dplyr::mutate( Requeue=dplyr::case_when(
          full_pass_perc>=perc ~ "FALSE",
          full_pass_perc<perc ~ "TRUE",
          TRUE ~ "TRUE") )
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} cut1_tib={RET}"))
          print(cut1_tib)
        }
        
        # if (!is.null(perc))
        #   cut2_tib <- cut2_tib %>% 
        #   dplyr::mutate( Requeue=dplyr::case_when(
        #     pass_perc>=perc ~ "FALSE", 
        #     pass_perc<perc ~ "TRUE",
        #     TRUE ~ "TRUE") )
      }
      cut1_tib <- dplyr::ungroup(cut1_tib)
      cut2_tib <- dplyr::ungroup(cut2_tib)
      
      ret_tib <- dplyr::inner_join(ret_tib,cut1_tib, by=c(type))
      ret_tib <- dplyr::inner_join(ret_tib,cut2_tib, by=c(type, des))
    }
    
    if (percision >= 0)
      ret_tib <- dplyr::mutate_if(ret_tib, is.numeric, list(round), percision)
    
    if (!is.null(save) && save==TRUE && !is.null(csv)) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing (percision={percision}) CSV={csv}.{RET}"))
      csv_dir <- base::dirname(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      readr::write_csv(ret_tib, csv)
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Sesame SSET Transform to Tib Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToTib = function(sset, source, name=NULL, man=NULL, mask=NULL,
                     percision=-1, sort=FALSE, 
                     fresh=FALSE, save=FALSE, csv=NULL, del='_',
                     by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                     verbose=0,vt=3,tc=1,tt=NULL,
                     funcTag='ssetToTib') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    by_sym  <- by %>% rlang::sym()
    des_sym <- des %>% rlang::sym()
    
    man_tib <- NULL
    if (!is.null(man)) man_tib <- dplyr::select(man,!!by, !!type, !!des)
    
    if (source=='sigs') {
      ret_tib <- dplyr::bind_rows(
        sset@IG   %>% tibble::as_tibble(rownames=by,.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='IG'),
        sset@oobG %>% tibble::as_tibble(rownames=by,.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='OG'),
        
        sset@IR   %>% tibble::as_tibble(rownames=by,.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='IR'),
        sset@oobR %>% tibble::as_tibble(rownames=by,.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='OR'),
        
        sset@II   %>% tibble::as_tibble(rownames=by,.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='2'),
        
        if (sset@ctl %>% base::nrow() != 0) {
          sset@ctl  %>% tibble::as_tibble(rownames=by,.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='2') %>%
            dplyr::rename(M=G,U=R) %>% dplyr::select(Probe_ID,!!des_sym,M,U)
        }
      ) %>% dplyr::select(!!by, !!des, dplyr::everything())

      # Add sig to the names::
      #
      ret_tib <- ret_tib %>% dplyr::rename(sig_M=M, sig_U=U)
      
      # Old Code to record original manifest reference designs::
      #   man_tib <- man_tib %>% dplyr::rename(Manifest_Design=!!des)
      if (!is.null(man_tib)) man_tib <- man_tib %>% dplyr::select(!!by, !!type)

    } else {
      if (is.null(name)) {
        if (fresh || is.null(sesame::extra(sset)[[source]]) ) {
          sset <- mutateSset(sset=sset, method=source,
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        }
        ret_tib <- sesame::extra(sset)[[source]] %>% 
          tibble::enframe(name=by, value=source)
      } else {
        if (fresh || 
            is.null(sesame::extra(sset)[[source]]) || 
            is.null(sesame::extra(sset)[[source]][[name]]) ) {
          sset <- mutateSset(sset=sset, method=name,
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        }
        ret_tib <- sesame::extra(sset)[[source]][[name]] %>% 
          tibble::enframe(name=by, value=paste(source,name, sep=del))
      }
    }

    if (!is.null(mask))
      ret_tib <- ret_tib %>% dplyr::filter(! (!!by_sym %in% mask) )

    if (percision!=-1) 
      ret_tib <- dplyr::mutate_if(ret_tib, is.numeric, list(round), percision)
    if (!is.null(man)) ret_tib <- dplyr::right_join(man_tib, ret_tib, by=by)
    if (sort) ret_tib <- dplyr::arrange(ret_tib, !!by_sym)
    
    if (save && !is.null(csv)) {
      safe_write(x=ret_tib,type="csv",file=csv,funcTag=funcTag,
                 verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      # if (verbose>=vt) 
      #   cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing (percision={percision}) CSV={csv}.{RET}"))
      # csv_dir <- base::dirname(csv)
      # if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      # readr::write_csv(ret_tib, csv)
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame SSET Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mutateSset = function(sset, method, full=TRUE,
                      quality.mask = FALSE, nondetection.mask = FALSE, 
                      correct.switch = FALSE, mask.use.tcga = FALSE, 
                      pval.threshold = 1, force=TRUE,
                      pval.method = "pOOBAH", sum.TypeI = FALSE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; Mutate Sesame({method})...{RET}"))
  
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} sset={RET}"))
    print(sset)
  }
  
  ctl_cnt <- sset@ctl %>% base::nrow()
  if (ctl_cnt==0 && method=='detectionPnegEcdf') {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Unable to mutate ctl_cnt={ctl_cnt} ",
                     "for method={method}. Returning original sset...{RET}"))
    return(sset)
  }
  
  ret_cnt <- 0
  stime <- system.time({
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} sset(ctl={ctl_cnt})={RET}"))
      print(sset)
    }
    
    oobR_ids = NULL
    oobG_ids = NULL
    if (!is.null(sset@extra$IGG) && !is.null(sset@extra$IRR)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Will use inferred channels...{RET}"))
      
      pass_ids_g <- rownames( sset@oobG)[ sset@extra$IRR]
      if (verbose>=vt+5) {
        pass_cnt_g <- pass_ids_g %>% length()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pass_ids_g({pass_cnt_g})={RET}"))
        pass_ids_g %>% head() %>% print()
      }
      
      fail_ids_g <- rownames( sset@oobG)[!sset@extra$IRR]
      if (verbose>=vt+5) {
        fail_cnt_g <- fail_ids_g %>% length()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} fail_ids_g({fail_cnt_g})={RET}"))
        fail_ids_g %>% head() %>% print()
      }
      
      pass_ids_r <- rownames( sset@oobR)[ sset@extra$IGG]
      if (verbose>=vt+5) {
        pass_cnt_r <- pass_ids_r %>% length()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pass_ids_r({pass_cnt_r})={RET}"))
        pass_ids_r %>% head() %>% print()
      }
      
      fail_ids_r <- rownames( sset@oobR)[!sset@extra$IGG]
      if (verbose>=vt+5) {
        fail_cnt_r <- fail_ids_r %>% length()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} fail_ids_r({fail_cnt_r})={RET}"))
        fail_ids_r %>% head() %>% print()
      }
      
      if (full) {
        oobG_ids <- c(pass_ids_g,fail_ids_g)
        oobR_ids <- c(pass_ids_r,fail_ids_r)
      } else {
        oobG_ids <- c(pass_ids_g)
        oobR_ids <- c(pass_ids_r)
      }
      
      if (verbose>=vt+4) {
        inbR_cnt <- sset@IR %>% base::nrow()
        inbG_cnt <- sset@IG %>% base::nrow()
        inbT_cnt <- inbR_cnt+inbG_cnt
        
        oobR_cnt <- oobR_ids %>% length()
        oobG_cnt <- oobG_ids %>% length()
        oobT_cnt <- oobR_cnt+oobG_cnt
        
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} method={method}; ",
                       "G={inbG_cnt}/{oobG_cnt}, ",
                       "R={inbR_cnt}/{oobR_cnt}, ",
                       "T={inbT_cnt}/{oobT_cnt}.{RET}{RET}"))
      }
    } else {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Will NOT use inferred channels...{RET}"))
      
      stopifnot(is.null(oobR_ids))
      stopifnot(is.null(oobG_ids))
    }
    
    if (is.null(method)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Missing method!!!{RET}{RET}"))
    } else if (method=='open') {
      sset <- sset %>% 
        sesame::pOOBAH(force=force) %>% 
        # sesame::noob(oobRprobes=oobR_ids, oobGprobes=oobG_ids) %>% 
        sesame::noob(bgR=oobR_ids, bgG=oobG_ids) %>% 
        sesame::dyeBiasCorrTypeINorm()
    } else if (method=='dyeBiasCorrTypeINorm') {
      sset <- sset %>% sesame::dyeBiasCorrTypeINorm()
    } else if (method=='detectionPnegEcdf' || method=='PnegEcdf') {
      sset <- sset %>% sesame::detectionPnegEcdf(force=force)
    } else if (method=='pOOBAH') {
      sset <- sset %>% sesame::pOOBAH(force=force)
    } else if (method=='noob') {
      # sset <- noob2(sset=sset,
      #               oobRprobes=oobR_ids, oobGprobes=oobG_ids,
      #               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- sset %>% sesame::noob(bgR=oobR_ids, bgG=oobG_ids)
    } else if (method=='noobsb') {
      sset <- sset %>% sesame::noobsb()
    } else if (method=='inferTypeIChannel') {
      sset <- sset %>% sesame::inferTypeIChannel(switch_failed=correct.switch, verbose=FALSE)
    } else if (method=='betas') {
      # sesame::extra(sset)[[method]] <- sesame::getBetas(sset=sset, mask=quality.mask,sum.TypeI=sum.TypeI)
      sesame::extra(sset)[[method]] <- getBetas2(sset=sset, mask=quality.mask,sum.TypeI=sum.TypeI)
    } else if (method=='r' || method=='raw') {
      # sset <- sset
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported method={method}!!!{RET}{RET}"))
    }
    
    ret_cnt <- sset %>% slotNames() %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} sset(slots={ret_cnt})={RET}"))
      print(sset)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  sset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Sesame SSET Initialization Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

newSset = function(prefix, platform, manifest, 
                   load=FALSE, save=FALSE, rds=NULL,
                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'newSset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; platform={platform}, prefix={prefix}.{RET}"))
  
  ret_cnt <- 0
  stime <- system.time({
    
    if (load && !is.null(rds) && file.exists(rds)) {
      if (verbose>=vt+1) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading RDS={rds}.{RET}"))
      sset <- readr::read_rds(rds)
    } else {
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} manifest={RET}"))
        manifest %>% print()
      }
      
      sset <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest)
      if (save && !is.null(rds)) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
        if (verbose>=vt+4) sset %>% print()
        readr::write_rds(sset, rds, compress="gz")
      }
    }
    
    ret_cnt <- sset %>% slotNames() %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} sset(slots={ret_cnt})={RET}"))
      print(sset)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  sset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Extracted Sesame Functions:: Betas
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getBetas2 = function (sset, mask = FALSE, sum.TypeI = FALSE, extra=TRUE) 
{
  if (is(sset, "SigSetList")) {
    return(do.call(cbind, lapply(sset, getBetas2, mask = mask, 
                                 sum.TypeI = sum.TypeI)))
  }
  stopifnot(is(sset, "SigSet"))
  IGs <- sesame::IG(sset)
  IRs <- sesame::IR(sset)
  if (sum.TypeI) {
    IGs <- IGs + oobR2(sset)
    IRs <- IRs + oobG2(sset)
  }
  else if (extra && !is.null(sset@extra$IGG) && !is.null(sset@extra$IRR)) {
    IGs[!sset@extra$IGG, ] <- sset@oobR[!sset@extra$IGG, ]
    IRs[!sset@extra$IRR, ] <- sset@oobG[!sset@extra$IRR, ]
  }
  betas <- c(pmax(IGs[,"M"], 1)/pmax(IGs[,"M"] + IGs[, "U"], 2),
             pmax(IRs[,"M"], 1)/pmax(IRs[,"M"] + IRs[, "U"], 2),
             pmax(sesame::II(sset)[,"M"], 1)/
               pmax(sesame::II(sset)[,"M"] + sesame::II(sset)[,"U"], 2))
  if (mask) 
    betas[!is.na(match(names(betas), sset@extra$mask))] <- NA
  betas
}

sset_pOOBAH2 = function (sset, force = FALSE) 
{
  stopifnot(is(sset, "SigSet"))
  method <- "pOOBAH"
  if (!force && method %in% names(extra(sset)$pvals)) {
    # cat("Retuning original sset\n")
    return(sset)
  }
  funcG <- ecdf(sesame::oobG(sset))
  funcR <- ecdf(sesame::oobR(sset))
  # funcR <- ecdf(oobR(sset) %>% head())
  
  pIR <- 1 - apply(cbind(funcR(sesame::IR(sset)[,"M"]), 
                         funcR(sesame::IR(sset)[,"U"])), 1, max)
  pIG <- 1 - apply(cbind(funcG(sesame::IG(sset)[,"M"]), 
                         funcG(sesame::IG(sset)[,"U"])), 1, max)
  pII <- 1 - apply(cbind(funcG(sesame::II(sset)[,"M"]), 
                         funcR(sesame::II(sset)[,"U"])), 1, max)
  names(pIR) <- rownames(sesame::IR(sset))
  names(pIG) <- rownames(sesame::IG(sset))
  names(pII) <- rownames(sesame::II(sset))
  if (!("pvals" %in% names(sesame::extra(sset)))) {
    # cat("Creating fresh list\n")
    sesame::extra(sset)[["pvals"]] <- list()
  }
  # cat(glue::glue("Setting method={method}{RET}{RET}"))
  pIR %>% head() %>% print()
  sesame::IR(sset) %>% head() %>% print()
  # oobR(sset) %>% head() %>% print()
  print(funcR)
  
  sesame::extra(sset)[["pvals"]][[method]] <- c(pIR, pIG, pII)
  sset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Extracted Sesame Functions:: Sex
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

inferSexKaryotypes2 = function (sset) 
{
  stopifnot(is(sset, "SigSet"))
  sex.info <- getSexInfo2(sset)
  auto.median <- median(sex.info[paste0("chr", seq_len(22))], 
                        na.rm = TRUE)
  XdivAuto <- sex.info["medianX"]/auto.median
  YdivAuto <- sex.info["medianY"]/auto.median
  if (XdivAuto > 1.2) {
    if (sex.info["fracXlinked"] >= 0.5) 
      sexX <- "XaXi"
    else if (sex.info["fracXmeth"] > sex.info["fracXunmeth"]) 
      sexX <- "XiXi"
    else sexX <- "XaXa"
  }
  else {
    if (sex.info["fracXmeth"] > sex.info["fracXunmeth"]) 
      sexX <- "Xi"
    else sexX <- "Xa"
  }
  if ((sexX == "Xi" || sexX == "Xa") && XdivAuto >= 1 && sex.info["fracXlinked"] >= 
      0.5) 
    sexX <- "XaXi"
  if (YdivAuto > 0.3 || sex.info["medianY"] > 2000) 
    sexY <- "Y"
  else sexY <- ""
  karyotype <- paste0(sexX, sexY)
  karyotype
}

inferSex2 = function (sset) 
{
  stopifnot(is(sset, "SigSet"))
  sex.info <- getSexInfo2(sset)[seq_len(3)]
  as.character(predict(sesameDataGet("sex.inference"), sex.info))
}

getSexInfo2 = function (sset) 
{
  if (is(sset, "SigSetList")) 
    return(do.call(cbind, lapply(sset, getSexInfo)))
  stopifnot(is(sset, "SigSet"))
  cleanY <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrY.clean
  # cat("length(cleanY)=\n")
  # cleanY %>% length() %>% print()
  # cleanY %>% head() %>% print()
  
  xLinked <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrX.xlinked
  # cat("length(xLinked)=\n")
  # xLinked %>% length() %>% print()
  # xLinked %>% head() %>% print()
  # cat("\n\n\n")
  
  probe2chr <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$probe2chr.hg19
  # cat("length(probe2chr)=\n")
  # probe2chr %>% length() %>% print()
  # probe2chr %>% head() %>% print()
  # cat("\n\n\n")
  
  # xLinkedBeta <- sesame::getBetas(sset=sesame::subsetSignal(sset, xLinked))
  xLinkedBeta <- getBetas2(sset=sesame::subsetSignal(sset, xLinked), extra=FALSE)
  # cat("head(xLinkedBeta)=\n")
  # xLinkedBeta %>% length() %>% print()
  # xLinkedBeta %>% head() %>% print()
  
  intens <- sesame::totalIntensities(sset)
  probes <- intersect(names(intens), names(probe2chr))
  intens <- intens[probes]
  probe2chr <- probe2chr[probes]
  # print(probe2chr)
  # return( sesame::subsetSignal(sset, cleanY) )
  # return( median(sesame::totalIntensities(sesame::subsetSignal(sset, cleanY))) )
  # cat("\n\n\n")
  
  c(medianY = median(sesame::totalIntensities(sesame::subsetSignal(sset, cleanY)), na.rm=TRUE), 
    medianX = median(sesame::totalIntensities(sesame::subsetSignal(sset, xLinked)), na.rm=TRUE), fracXlinked = 
      (sum(xLinkedBeta > 0.3 & xLinkedBeta < 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta))) ), 
    fracXmeth = (sum(xLinkedBeta > 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    fracXunmeth = (sum(xLinkedBeta < 0.3, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    tapply(intens, probe2chr, median, na.rm=TRUE))
}

inferEthnicity2 = function (sset) 
{
  if (is(sset, "SigSetList")) 
    return(vapply(sset, inferEthnicity2, character(1)))
  stopifnot(is(sset, "SigSet"))
  stopifnot(sset@platform %in% c("EPIC", "HM450"))
  ethnicity.inference <- sesameDataGet("ethnicity.inference")
  ccsprobes <- ethnicity.inference$ccs.probes
  rsprobes <- ethnicity.inference$rs.probes
  ethnicity.model <- ethnicity.inference$model
  af <- c(getBetas2(sset, mask = FALSE)[rsprobes], 
          getAFTypeIbySumAlleles(subsetSignal(sset,ccsprobes)))
  as.character(predict(ethnicity.model, af))
}

noob2 = function (sset, oobRprobes = NULL, oobGprobes = NULL, offset = 15,
                  verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'noob2'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} sset={RET}"))
    print(sset)
    
    oobRprobes_cnt <- oobRprobes %>% length()
    cat(glue::glue("[{funcTag}]:{tabsStr} oobRprobes({oobRprobes_cnt})={RET}"))
    oobRprobes %>% unique() %>% head() %>% print()
    
    oobGprobes_cnt <- oobGprobes %>% length()
    cat(glue::glue("[{funcTag}]:{tabsStr} oobGprobes({oobGprobes_cnt})={RET}"))
    oobGprobes %>% unique() %>% head() %>% print()
  }
  
  stime <- system.time({
    
    if (all(sesame::oobG(sset) == 0) || all(oobR(sset) == 0)) {
      return(sset)
    }
    ibR <- c(sesame::IR(sset), sesame::II(sset)[, "U"])
    ibG <- c(sesame::IG(sset), sesame::II(sset)[, "M"])
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ibR(1)={RET}"))
      ibR %>% names() %>% unique() %>% head() %>% print()
      cat(glue::glue("[{funcTag}]:{tabsStr} ibG(1)={RET}"))
      ibG %>% names() %>% unique() %>% head() %>% print()
    }
    
    ibR[ibR == 0] <- 1
    ibG[ibG == 0] <- 1
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ibR(2)={RET}"))
      ibR %>% names() %>% unique() %>% head() %>% print()
      cat(glue::glue("[{funcTag}]:{tabsStr} ibG(2)={RET}"))
      ibG %>% names() %>% unique() %>% head() %>% print()
    }
    
    sesame::oobR(sset)[sesame::oobR(sset) == 0] <- 1
    sesame::oobG(sset)[sesame::oobG(sset) == 0] <- 1
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} oobR={RET}"))
      sesame::oobR(sset) %>% head() %>% print()
      cat(glue::glue("[{funcTag}]:{tabsStr} oobG={RET}"))
      sesame::oobG(sset) %>% head() %>% print()
    }
    
    if (is.null(oobRprobes)) {
      real_oobR <- sesame::oobR(sset)
    } else {
      real_oobR <- sesame::oobR(sset)[oobRprobes, ]
    }
    #
    # Red Verbose::
    #
    if (verbose>=vt+4) {
      real_oobR_cnt <- base::nrow(real_oobR)
      cat(glue::glue("real_oobR({real_oobR_cnt})={RET}"))
      real_oobR %>% head() %>% print()
    }
    
    if (is.null(oobGprobes)) {
      real_oobG <- sesame::oobG(sset)
    } else {
      real_oobG <- sesame::oobG(sset)[oobGprobes, ]
    }
    #
    # Grn Verbose::
    #
    if (verbose>=vt+4) {
      real_oobG_cnt <- base::nrow(real_oobG)
      cat(glue::glue("real_oobG({real_oobG_cnt})={RET}"))
      real_oobG %>% head() %>% print()
    }
    
    ibR.nl <- 
      backgroundCorrectionNoobCh1_2(ibR, real_oobR, sesame::ctl(sset)$R, offset = offset)
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ibR.nl={RET}"))
      ibR.nl %>% names() %>% print()
    }
    
    ibG.nl <- 
      backgroundCorrectionNoobCh1_2(ibG, real_oobG, sesame::ctl(sset)$G, offset = offset)
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ibG.nl={RET}"))
      ibG.nl %>% names() %>% print()
    }
    
    if (length(sesame::IG(sset)) > 0) {
      sesame::IG(sset) <- 
        matrix(ibG.nl$i[seq_along(sesame::IG(sset))], 
               nrow = nrow(sesame::IG(sset)), 
               dimnames = dimnames(sesame::IG(sset)))
    } else {
      sesame::IG(sset) <- 
        matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("M", "U")))
    }
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} IG={RET}"))
      sesame::IG(sset) %>% head() %>% print()
    }
    
    if (length(sesame::IR(sset)) > 0) {
      sesame::IR(sset) <- matrix(ibR.nl$i[seq_along(IR(sset))], nrow = nrow(sesame::IR(sset)), 
                                 dimnames = dimnames(sesame::IR(sset)))
    } else {
      sesame::IR(sset) <- matrix(ncol = 2, nrow = 0, dimnames = list(NULL, 
                                                                     c("M", "U")))
    }
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} IR={RET}"))
      sesame::IR(sset) %>% head() %>% print()
    }
    
    if (nrow(sesame::II(sset)) > 0) {
      if (verbose>=vt+4)
        cat(glue::glue("[{funcTag}]:{tabsStr} nrow(sesame::II(sset)) > 0 == TRUE{RET}"))
      
      sesame::II(sset) <- 
        as.matrix(data.frame(M = ibG.nl$i[(length(sesame::IG(sset)) + 1):length(ibG)], 
                             U = ibR.nl$i[(length(sesame::IR(sset)) + 1):length(ibR)], 
                             row.names = rownames(sesame::II(sset))))
    } else {
      if (verbose>=vt+4)
        cat(glue::glue("[{funcTag}]:{tabsStr} nrow(sesame::II(sset)) > 0 == FALSE{RET}"))
      
      sesame::II(sset) <- 
        matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("M", "U")))
    }
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} II={RET}"))
      sesame::II(sset) %>% head() %>% print()
    }
    
    sesame::ctl(sset)$G <- ibG.nl$c
    sesame::ctl(sset)$R <- ibR.nl$c
    sesame::oobR(sset)  <- ibR.nl$o
    sesame::oobG(sset)  <- ibG.nl$o
    
    ret_cnt <- sset %>% slotNames() %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} sset(slots={ret_cnt})={RET}"))
      print(sset)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  sset
}

backgroundCorrectionNoobCh1_2 = function(ib, oob, ctl, offset=15) {
  
  e <- MASS::huber(oob)
  mu <- e$mu
  sigma <- e$s
  alpha <- pmax(MASS::huber(ib)$mu-mu, 10)
  list(
    i=offset+normExpSignal2(mu, sigma, alpha, ib),
    c=offset+normExpSignal2(mu, sigma, alpha, ctl),
    o=offset+normExpSignal2(mu, sigma, alpha, oob))
}

normExpSignal2 <- function (mu, sigma, alpha, x)  {
  
  sigma2 <- sigma * sigma
  
  if (alpha <= 0)
    stop("alpha must be positive")
  if (sigma <= 0)
    stop("sigma must be positive")
  
  mu.sf <- x - mu - sigma2/alpha
  signal <- mu.sf + sigma2 * exp(
    dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
      pnorm(
        0, mean = mu.sf, sd = sigma,
        lower.tail = FALSE, log.p = TRUE))
  
  o <- !is.na(signal)
  if (any(signal[o] < 0)) {
    warning("Limit of numerical accuracy reached with very
low intensity or very high background:\nsetting adjusted intensities
to small value")
    signal[o] <- pmax(signal[o], 1e-06)
  }
  signal
}

# End of file

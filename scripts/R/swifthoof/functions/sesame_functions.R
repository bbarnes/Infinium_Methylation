
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
#                         Sesame Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToSummary = function(sset, man, idx, workflow, name, platform=NULL,
                         outDir=NULL, pre=NULL, ref=NULL, mask=NULL,
                         pvals=NULL, min_pvals=NULL, min_percs=NULL,
                         basic=NULL,
                         report_vec=c("Requeue","pass_perc","mean",
                                      "basis_r2_val", "basis_dB_val"),
                         
                         write_sset=FALSE, sset_rds=NULL, ret_sset=FALSE,
                         
                         write_beta=FALSE, beta_csv=NULL, ret_beta=FALSE,
                         write_bsum=FALSE, bsum_csv=NULL, ret_bsum=FALSE,
                         
                         write_pval=FALSE, pval_csv=NULL, ret_pval=FALSE,
                         write_psum=FALSE, psum_csv=NULL, ret_psum=FALSE,
                         
                         write_sigs=FALSE, sigs_csv=NULL, ret_sigs=FALSE,
                         write_ssum=FALSE, ssum_csv=NULL, ret_ssum=FALSE,
                         
                         write_call=FALSE, call_csv=NULL, ret_call=FALSE,
                         write_csum=FALSE, csum_csv=NULL, ret_csum=FALSE,
                         
                         write_snps=FALSE, snps_csv=NULL, ret_snps=FALSE,
                         write_auto=FALSE, plot_auto=FALSE, makek_pred=TRUE,
                         
                         percision_sigs=-1,percision_beta=-1,percision_pval=-1,
                         by="Probe_ID", type="Probe_Type", des="Probe_Design",
                         minDb=0.02, dpi=120, plotFormat="png", datIdx=5,
                         non_ref=FALSE, fresh=FALSE, del='_',
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='ssetToSummary') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}...{RET}"))
  
  ret_cnt <- 0
  ret_dat <- NULL
  
  # Validate Inputs::
  #
  pvals_cnt <- 0
  min_pvals_cnt <- 0
  min_percs_cnt <- 0
  if (!is.null(pvals)) pvals_cnt <- length(pvals)
  if (!is.null(min_pvals)) min_pvals_cnt <- length(min_pvals)
  if (!is.null(min_percs)) min_percs_cnt <- length(min_percs)
  
  if (pvals_cnt != min_pvals_cnt || min_pvals_cnt != min_percs_cnt) {
    stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR; pval vecs not equal; pvals_cnt={pvals_cnt}, ",
                    "min_pvals_cnt={min_pvals_cnt}, min_percs_cnt={min_percs_cnt}!!!{RET}"))
    return(ret_dat)
  }
  
  sigs_dat_tib <- NULL
  call_dat_tib <- NULL
  auto_ssh_tib <- NULL
  sums_dat_tib <- NULL
  sums_ssh_tib <- NULL
  pred_dat_tib <- NULL
  data_ssh_tib <- NULL
  beta_sum_tib <- NULL
  
  stime <- system.time({
    
    # Initialize Outputs::
    #
    if (!is.null(name) && !is.null(outDir)) {
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
      
      out_name <- paste(name,workflow, sep=del)
      if (is.null(sset_rds)) 
        sset_rds <- file.path(outDir, paste(out_name,'sset.rds', sep='.'))
      
      if (is.null(beta_csv))
        beta_csv <- file.path(outDir, paste(out_name,'beta.dat.csv.gz', sep='.'))
      if (is.null(bsum_csv))
        bsum_csv <- file.path(outDir, paste(out_name,'beta.sum.csv.gz', sep='.'))
      
      if (is.null(sigs_csv))
        sigs_csv <- file.path(outDir, paste(out_name,'sigs.dat.csv.gz', sep='.'))
      if (is.null(ssum_csv))
        ssum_csv <- file.path(outDir, paste(out_name,'sigs.sum.csv.gz', sep='.'))
      
      if (is.null(call_csv)) 
        call_csv <- file.path(outDir, paste(out_name,'call.dat.csv.gz', sep='.'))
      if (is.null(csum_csv)) 
        csum_csv <- file.path(outDir, paste(out_name,'call.sum.csv.gz', sep='.'))
      
      if (is.null(snps_csv))
        snps_csv <- file.path(outDir, paste(out_name,'snps.vcf', sep='.'))
      
      # if (write_sigs) sigs_csv <- clean_file(sigs_csv, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    # OLD CODE TO BE REMOVED::
    #  - This was left over from concerns about memory leaks early on...
    # Clear sset to ensure there's no code mistakes
    # sset_dat <- sset
    # sset <- NULL
    
    by_sym <- rlang::sym(by)
    
    # call_dat_tib <- man %>% dplyr::select(!!by)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Set/Update/Summarize:: Pvals
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_dat$pval <- NULL
    pval_cut_tib <- tibble::tibble(Method=pvals, min_pval=min_pvals, min_perc=min_percs)
    for (pval_key in pval_cut_tib$Method) {
      
      pval_out_str <- paste('pval',pval_key, sep='-')
      pval_csv <- file.path(outDir, paste(out_name,pval_out_str,'dat.csv.gz', sep='.'))
      psum_csv <- file.path(outDir, paste(out_name,pval_out_str,'sum.csv.gz', sep='.'))
      
      ret_pval_dat <- NULL
      pval_tib <- pval_cut_tib %>% dplyr::filter(Method==pval_key)
      min_pval <- pval_tib$min_pval[1]
      min_perc <- pval_tib$min_perc[1]
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval={pval_key}; min={min_pval}, perc={min_perc}.{RET}"))
      
      #
      # Need to repeat this part if mask is present and update names...
      #
      pval_dat_tib <- NULL
      pval_dat_tib <- ssetToTib(
        sset=sset, source='pvals', name=pval_key,
        percision=percision_pval, sort=FALSE,
        save=write_pval, csv=pval_csv,
        by=by, type=type, des=des,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      if (ret_pval) ret_pval_dat$pval_dat <- pval_dat_tib
      
      if (verbose>=vt+6) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}pval_sum_tab({pval_key},min={min_pval},per={min_perc})={RET}"))
        ret_cnt <- print_tib(pval_dat_tib,funcTag, verbose,vt+4,tc, n="pval_dat_tib")
        
        pval_sum_tab <- pval_dat_tib %>% 
          purrr::set_names(c("Probe_ID","Pval")) %>%
          dplyr::left_join(man, by="Probe_ID") %>%
          dplyr::group_by(Probe_Type,Probe_Design) %>%
          dplyr::summarise(Pass_Perc=cntPer_lte(Pval,min_pval),
                           .groups="drop")
        ret_cnt <- print_tib(pval_sum_tab,funcTag, verbose,vt+4,tc, n="pval_sum_tab")
        cat(glue::glue("{RET}{RET}{RET}"))
      }
      
      pval_sum_tib <- NULL
      pval_sum_tib <- ssetTibToSummary(
        tib=pval_dat_tib,man=man,
        pval=min_pval, perc=min_perc, percision=percision_pval,
        save=write_psum, csv=psum_csv,
        by=by, type=type, des=des,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      if (ret_psum) ret_pval_dat$pval_sum <- pval_sum_tib
      if (ret_pval || ret_psum) ret_dat$pval[[pval_key]] = ret_pval_dat
      
      if (is.null(call_dat_tib)) {
        call_dat_tib <- pval_dat_tib
      } else {
        call_dat_tib <- dplyr::left_join(call_dat_tib,pval_dat_tib, by=by)
      }
      sums_dat_tib <- sums_dat_tib %>% dplyr::bind_rows(pval_sum_tib)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}pval_dat_tib({pval_key},min={min_pval},per={min_perc})={RET}"))
        ret_cnt <- print_tib(pval_dat_tib,funcTag, verbose,vt+4,tc, n="pval_dat_tib")
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}pval_sum_tib({pval_key},min={min_pval},per={min_perc})={RET}"))
        ret_cnt <- print_tib(pval_sum_tib,funcTag, verbose,vt+4,tc, n="pval_sum_tib")
        cat(glue::glue("{RET}RET}"))
      }
      
      if (!is.null(mask)) {
        mask_pval_key <- paste('Mask',pval_key, sep='_')
        
        pval_dat_tib <- pval_dat_tib %>% dplyr::filter(! (!!by_sym %in% mask) )

        pval_sum_tib <- NULL
        pval_sum_tib <- ssetTibToSummary(
          tib=pval_dat_tib,man=man,
          pval=min_pval, perc=min_perc, percision=percision_pval,
          save=FALSE, csv=NULL,
          by=by, type=type, des=des,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        pval_sum_tib <- pval_sum_tib %>% dplyr::mutate(Metric=paste(Metric,'Mask',sep=del))
        sums_dat_tib <- sums_dat_tib %>% dplyr::bind_rows(pval_sum_tib)
        
        if (verbose>=vt+6) {
          cat(glue::glue("{RET} sums_dat_tib(mask_pval_key={mask_pval_key})={RET}"))
          ret_cnt <- print_tib(sums_dat_tib,funcTag, verbose,vt+4,tc, n="sums_dat_tib")
          cat(glue::glue("{RET}RET}"))
        }
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Summarize:: Betas
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_beta_dat <- NULL
    beta_dat_tib <- NULL
    beta_dat_tib <- ssetToTib(
      sset=sset, source='betas',
      percision=percision_beta, sort=TRUE,
      save=write_beta, csv=beta_csv,
      by=by, type=type, des=des,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (ret_beta) ret_beta_dat$beta_dat <- beta_dat_tib
    
    call_dat_tib <- dplyr::left_join(call_dat_tib,beta_dat_tib, by=by)
    sums_dat_tib <- sums_dat_tib %>% dplyr::bind_rows(beta_sum_tib)
    
    beta_sum_tib <- ssetTibToSummary(
      tib=beta_dat_tib,man=man,
      percision=percision_beta,
      save=write_bsum, csv=bsum_csv,
      by=by, type=type, des=des,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (ret_bsum) ret_beta_dat$beta_sum <- beta_sum_tib
    
    # Add open_sesame comparison if provided::
    r2_basic_val <- NULL
    dB_basic_val <- NULL
    base_dat_tib <- NULL
    if (!is.null(basic)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting Basic(Open) Comparison...{RET}"))
      
      r2_basic_val <- basic %>% 
        # tibble::enframe(name="Probe_ID", value="betas") %>%
        dplyr::inner_join(beta_dat_tib, 
                          by="Probe_ID", suffix=c("_ref","_can")) %>%
        tibble::column_to_rownames(var="Probe_ID") %>% 
        as.matrix() %>% cor() %>% as_tibble() %>% 
        head(n=1) %>% dplyr::pull(2)
      
      if (percision_pval != -1)
        r2_basic_val <- round(r2_basic_val, percision_pval)
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Basic(r2)={r2_basic_val}...{RET}"))
      
      dB_basic_val <- basic %>% 
        # tibble::enframe(name="Probe_ID", value="betas") %>%
        dplyr::inner_join(beta_dat_tib, 
                          by="Probe_ID", suffix=c("_ref","_can")) %>%
        dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
        dplyr::mutate(dB=base::abs(betas_ref-betas_can)) %>% 
        dplyr::summarise(pass_perc=cntPer_lte(dB,minDb)) %>%
        head(n=1) %>% dplyr::pull(1)
      
      if (percision_beta != -1)
        dB_basic_val <- round(dB_basic_val, percision_beta)
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Basic(dB)={dB_basic_val}...{RET}"))
      
      r2_basic_key <- paste('cg_basic_r2', sep=del)
      dB_basic_key <- paste('cg_basic_dB', sep=del)
      base_dat_tib <- tibble::tibble(
        !!r2_basic_key := r2_basic_val,
        !!dB_basic_key := dB_basic_val )

      if (verbose>=vt+6) {
        cat(glue::glue("{RET}"))
        ret_cnt <- print_tib(base_dat_tib,funcTag, verbose,vt+4,tc, n="base_dat_tib")
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ",
                       "base_dat_tib... done...{RET}{RET}"))
      }
    }
    if (ret_beta || ret_bsum) ret_dat$beta <- ret_beta_dat
    
    # Run Auto Sample Detection::
    #
    if (!is.null(ref) &&
        pvals_cnt != 0 && min_pvals_cnt != 0) {
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting Sample Auto Detection...{RET}"))
      
      auto_beta_key <- paste('betas', sep=del)
      auto_negs_key <- paste('pvals',pvals[pvals_cnt], sep=del)
      auto_min_pval <- min_pvals[min_pvals_cnt]
      
      if (verbose>=vt) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}{TAB} auto_beta_key={auto_beta_key}.{RET}"))
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}{TAB} auto_negs_key={auto_negs_key}.{RET}"))
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}{TAB} auto_min_pval={auto_min_pval}.{RET}"))
      }
      
      # Simplifying by_sym to only CG#'s::
      #   TBD:: Should sort by detection p-value...
      #
      auto_man_tib <- man %>%
        dplyr::mutate(!!by_sym := stringr::str_remove(!!by_sym, "_.*$")) %>%
        dplyr::distinct(!!by_sym, .keep_all=TRUE)
      ret_cnt <- print_tib(auto_man_tib,funcTag, verbose,vt+4,tc, n="auto_man_tib")
      ret_cnt <- print_tib(call_dat_tib,funcTag, verbose,vt+4,tc, n="call_dat_tib")
      
      auto_dat_tib <- call_dat_tib %>% 
        dplyr::mutate(!!by_sym := stringr::str_remove(!!by_sym, "_.*$")) %>%
        dplyr::distinct(!!by_sym, .keep_all=TRUE)
      ret_cnt <- print_tib(auto_dat_tib,funcTag, verbose,vt+4,tc, n="auto_dat_tib")
      
      auto_ssh_tib <- autoDetect_Wrapper(
        can=auto_dat_tib, ref=ref, man=auto_man_tib, # mask=mask,
        minPval=auto_min_pval, minDelta=minDb,
        dname='Design_Type', pname=type, ptype='cg', jval=by, 
        field=auto_beta_key, pval=auto_negs_key, suffix='beta', del=del,
        outDir=outDir, sname=out_name, 
        plotMatrix=plot_auto, writeMatrix=write_auto,
        dpi=dpi, format=plotFormat, datIdx=datIdx, non.ref=non_ref,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      #
      # Add dB against best auto detect::
      #
      # auto_sam_tib %>% dplyr::select(dplyr::all_of(c("Probe_ID", rdat$ssheet_tib$AutoSample_dB_1_Key_1)))
      #
      dB_ref_tib <- ref %>% dplyr::select(dplyr::all_of(c("Probe_ID", auto_ssh_tib$AutoSample_dB_1_Key))) %>%
        purrr::set_names(c("Probe_ID","ref_beta")) %>%
        dplyr::inner_join(auto_dat_tib, by="Probe_ID") %>%
        dplyr::mutate(dB_ref=ref_beta-betas) %>%
        dplyr::select(Probe_ID, dB_ref)
        # dplyr::mutate(dB_ref=abs(ref_beta-betas))
      ret_cnt <- print_tib(dB_ref_tib,funcTag, verbose,vt-4,tc, n="dB_ref_tib")

      call_dat_tib <- call_dat_tib %>% dplyr::left_join(dB_ref_tib, by="Probe_ID")
      ret_cnt <- print_tib(call_dat_tib,funcTag, verbose,vt-4,tc, n="call_dat_tib")
      
      ret_cnt <- print_tib(auto_ssh_tib,funcTag, verbose,vt+4,tc, n="auto_ssh_tib")
      
    } else {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Skipping Sample Auto Detection...{RET}"))
    }
    
    # Write Updated SSET
    #
    if (write_sset && !is.null(sset_rds) && !file.exists(sset_rds))
      readr::write_rds(sset, sset_rds, compress="gz")
    
    # Write Calls CSV
    #
    if (write_call && !is.null(call_csv))
      readr::write_csv(call_dat_tib, call_csv)
    
    call_sum_tib <- NULL
    call_sum_tib <-
      callToPassPerc(file=call_csv, key="pvals_pOOBAH",
                     name=NULL, idx=NULL, min=min_pvals[1], type='cg',
                     verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # Write VCF SNPs calls::
    #
    if (write_snps && !is.null(snps_csv))
      vcf_ret <- safeVCF(sset=sset, vcf=snps_csv,
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Set/Summarize:: Sigs
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_sigs_dat <- NULL
    sigs_dat_tib <- ssetToTib(
      sset=sset, man=man, 
      source='sigs', name=pval,
      percision=percision_sigs, sort=FALSE, 
      save=write_sigs, csv=sigs_csv, 
      by=by, type=type, des=des,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (ret_sigs) ret_sigs_dat$sigs_dat <- sigs_dat_tib
    
    sigs_sum_tib <- ssetTibToSummary(
      tib=sigs_dat_tib,
      percision=percision_sigs,
      save=write_ssum, csv=ssum_csv, 
      by=by, type=type, des=des,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (ret_ssum) ret_sigs_dat$sigs_sum <- sigs_sum_tib
    
    if (ret_sigs || ret_ssum) ret_dat$sigs <- ret_sigs_dat
    
    sums_dat_tib <- sums_dat_tib %>% dplyr::bind_rows(sigs_sum_tib)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                Update Calls, Signal and Summary Headers::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    fix_cols <- c(by,type,des)
    sigs_dat_tib <- addColNames(sigs_dat_tib, add=workflow, fix=fix_cols,
                                verbose=verbose,vt=vt+5,tc=tc+1,tt=tt)
    
    fix_cols <- c(by)
    call_dat_tib <- addColNames(call_dat_tib, add=workflow, fix=fix_cols,
                                verbose=verbose,vt=vt+5,tc=tc+1,tt=tt)
    
    sums_dat_tib <- sums_dat_tib %>%
      dplyr::mutate(Workflow_key=workflow, Workflow_idx=idx) %>% 
      dplyr::select(Workflow_key, Workflow_idx, dplyr::everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Predicted and Inferred Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Initialize Sample Sheet::
    #
    sums_ssh_tib <- tibble::tibble(Method_Key=workflow,Method_Idx=idx)
    
    if (makek_pred) pred_dat_tib <- 
      ssetToPredictions(
        sset=sset, platform=platform,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Make Sample Sheet::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    data_ssh_tib <- sums_dat_tib %>% 
      tidyr::unite(key, Probe_Type,Probe_Design,Metric, sep=del) %>% 
      tidyr::gather(metric, value, -key, -Workflow_key, -Workflow_idx) %>% 
      dplyr::filter(!is.na(value)) %>% 
      dplyr::filter(metric %in% report_vec) %>%
      tidyr::unite(key, key,metric, sep=del) %>% 
      dplyr::select(-Workflow_key, -Workflow_idx) %>% 
      tidyr::spread(key, value) %>%
      dplyr::select(dplyr::contains("_Requeue_"), 
                    dplyr::contains("_pass_perc_"),
                    dplyr::everything())
    
    ret_cnt <- print_tib(data_ssh_tib,funcTag, verbose,vt+4,tc, n="data_ssh_tib")
    # pvals_pOOBAH
    
    # NOTE:: Need to extract full p-values seperate from the rest of the stats
    #  to avoid duplicate entries...
    #
    # TBD:: FIX THIS ISSUE
    pval_ssh_tib <- sums_dat_tib %>%
      tidyr::unite(key, Probe_Type,Metric, sep=del) %>% 
      tidyr::gather(metric, value, -key, -Workflow_key, -Workflow_idx, -Probe_Design) %>% 
      dplyr::filter(!is.na(value)) %>% 
      dplyr::filter(metric %in% c('full_pass_perc')) %>%
      dplyr::mutate(metric='pass_perc') %>%
      tidyr::unite(key, key,metric, sep=del) %>% 
      dplyr::select(-Workflow_key, -Workflow_idx, -Probe_Design) %>%
      dplyr::distinct() %>%
      tidyr::spread(key, value) %>%
      dplyr::select(dplyr::contains('_pass_perc_'), 
                    dplyr::everything())
    
    if (verbose>=vt+4) {
      cat(glue::glue("{RET}"))
      ret_cnt <- print_tib(pval_ssh_tib,funcTag, verbose,vt+4,tc, n="pval_ssh_tib")
      cat(glue::glue("{RET}{RET}"))
    }
    
    sums_ssh_tib <- 
      dplyr::bind_cols(sums_ssh_tib,auto_ssh_tib,pred_dat_tib,call_sum_tib,
                       pval_ssh_tib,data_ssh_tib,base_dat_tib) %>%
      purrr::set_names(paste(names(.),idx, sep=del))
    
    #
    #  NOT NEEDED::
    #    tidyr::pivot_wider(id_cols=c(Workflow_key, Workflow_idx), 
    #                       names_from="key", values_from="value")
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Gather/Merge Results:: Sample Sheet
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (is.null(pre)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Building new data list...{RET}"))
      
      ret_dat$sigs_dat <- sigs_dat_tib
      ret_dat$call_dat <- call_dat_tib
      ret_dat$sums_dat <- sums_dat_tib
      ret_dat$sums_ssh <- sums_ssh_tib
      # ret_dat$pred_dat <- pred_dat_tib
    } else {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Join with previous data...{RET}"))
      
      ret_dat$sigs_dat <- dplyr::left_join(pre$sigs_dat, sigs_dat_tib, by=c(by,type,des))
      ret_dat$call_dat <- dplyr::left_join(pre$call_dat, call_dat_tib, by=by)
      ret_dat$sums_dat <- dplyr::bind_rows(pre$sums_dat, sums_dat_tib)
      ret_dat$sums_ssh <- dplyr::bind_cols(pre$sums_ssh, sums_ssh_tib)
      # ret_dat$pred_dat <- dplyr::bind_cols(pre$pred_dat, pred_dat_tib)
    }
    
    ret_cnt <- ret_dat %>% names %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_dat({ret_cnt})={RET}"))
      print(ret_dat)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_dat
}

mutateSSET_workflow = function(sset, workflow, pvals=NULL,
                               save=FALSE, rds=NULL, 
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSSET_workflow'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
  
  ret_cnt <- 0
  stime <- system.time({
    
    # Clear extra fields
    sesame::extra(sset)[['pvals']] <- NULL
    sesame::extra(sset)[['betas']] <- NULL
    
    #
    # TBD:: This could be parsed into individual letters and then processed in that order...
    #
    if (workflow=='r' || workflow=='raw') {
      # Do nothing...
    } else if (workflow=='i') {
      sset <- mutateSset(sset=sset, method = 'inferTypeIChannel',
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='n') {
      sset <- mutateSset(sset=sset, method = 'noob', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='d') {
      sset <- mutateSset(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='nd') {
      sset <- mutateSset(sset=sset, method = 'noob', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='dn') {
      sset <- mutateSset(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'noob', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='ind') {
      sset <- mutateSset(sset=sset, method = 'inferTypeIChannel', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'noob', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='ndi') {
      sset <- mutateSset(sset=sset, method = 'noob', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'inferTypeIChannel', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='din') {
      sset <- mutateSset(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'inferTypeIChannel', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'noob', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='dni') {
      sset <- mutateSset(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'noob', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSset(sset=sset, method = 'inferTypeIChannel', 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else {
      stop(glue::glue("[{funcTag}]: ERROR: Unsupported workflow={workflow}!{RET}{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Set:: Pvals
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    for (pval in pvals) {
      if (verbose>=vt+2) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Mutating pval={pval}...{RET}"))
      
      sset <- mutateSset(sset=sset, method=pval,
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Set:: Betas
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt+2) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Mutating betas...{RET}"))
    
    sset <- mutateSset(sset=sset, method="betas",
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Save Output::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (save && !is.null(rds)) {
      if (verbose>=vt+1) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
      readr::write_rds(sset, rds, compress="gz")
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
#                        Quick Sesame File Conversion::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

build_mask_list = function(verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'build_mask_list'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  mask_cpg_tib <- NULL
  stime <- system.time({
    
    hg19_man_grs  <- sesameData::sesameDataGet("EPIC.hg19.manifest")
    hg19_mask_tib <- hg19_man_grs %>% as.data.frame() %>% 
      tibble::rownames_to_column(var="Probe_ID") %>% 
      dplyr::filter(probeType=='cg') %>%
      dplyr::select(Probe_ID,MASK_general) %>% 
      tibble::as_tibble() %>% 
      dplyr::filter(MASK_general) %>% 
      dplyr::distinct(Probe_ID)
    
    hg38_man_grs  <- sesameData::sesameDataGet("EPIC.hg38.manifest")
    hg38_mask_tib <- hg38_man_grs %>% as.data.frame() %>% 
      tibble::rownames_to_column(var="Probe_ID") %>%
      dplyr::filter(probeType=='cg') %>%
      dplyr::select(Probe_ID,MASK_general) %>% 
      tibble::as_tibble() %>% 
      dplyr::filter(MASK_general) %>% 
      dplyr::distinct(Probe_ID)
    
    mask_cpg_csv <- file.path(par$datDir, 'manifest/mask/sesame-general-mask.cpg.csv.gz')
    mask_cpg_tib <- dplyr::bind_rows( hg19_mask_tib, hg38_mask_tib ) %>% 
      dplyr::distinct() %>% dplyr::arrange(Probe_ID)
    readr::write_csv(mask_cpg_tib, mask_cpg_csv)
    
    ret_cnt <- mask_cpg_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  mask_cpg_tib
}

# End of file

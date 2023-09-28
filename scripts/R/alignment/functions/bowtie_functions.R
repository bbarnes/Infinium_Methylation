
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Bowtie Alignment Methods::
#
#  NOTE:: This is pretty old and hasn't been used in a long time. It
#    should probably be move to the graveyard and re-written when its's
#    needed...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

COM <- ","
TAB <- "\t"
RET <- "\n"
BNG <- "|"
BRK <- paste0("# ",
              paste(rep("-----",6),collapse=" "),"|",
              paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # ret_cnt <- ret_tib %>% base::nrow()
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Bowtie Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadProbeAlignBowtieInfI = function(sam, reduced=FALSE, filtered=FALSE, flipSeq=FALSE,
                                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadProbeAlignBowtieInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} sam={sam}.{RET}"))
  
  col_vec <- c('QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL',
               'AS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YT')
  
  stime <- system.time({
    
    snp_raw_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file=sam, col_names=col_vec, comment='@') ))
    
    if (reduced)  snp_raw_tib <- snp_raw_tib %>% dplyr::select(QNAME:POS,SEQ,MD)
    if (filtered) snp_raw_tib <- snp_raw_tib %>% dplyr::filter(FLAG==0 | FLAG==16)
    if (flipSeq)  snp_raw_tib <- snp_raw_tib %>% dplyr::mutate(SEQ=case_when( FLAG==16 ~ revCmp(SEQ),TRUE ~ SEQ ))
    
    # Split by Infinium I Probe Design::
    sam_1A_tib <- snp_raw_tib %>% dplyr::filter(stringr::str_ends(QNAME,'_IA')) %>% dplyr::mutate(QNAME=stringr::str_remove(QNAME, '_IA$'))
    sam_1B_tib <- snp_raw_tib %>% dplyr::filter(stringr::str_ends(QNAME,'_IB')) %>% dplyr::mutate(QNAME=stringr::str_remove(QNAME, '_IB$'))
    sam_tib <- dplyr::inner_join(sam_1A_tib,sam_1B_tib, by=c("QNAME","FLAG","RNAME","POS"), suffix=c("_IA", "_IB")) # %>%
    
    MD_IA_cnt <- sam_tib %>% dplyr::filter(MD_IA %in% art_snp_md_vec) %>% base::nrow()
    MD_IB_cnt <- sam_tib %>% dplyr::filter(MD_IB %in% art_snp_md_vec) %>% base::nrow()
    if (verbose>=vt) cat(glue::glue("[{par$prgmTag}]: IA/IB = {MD_IA_cnt}, {MD_IB_cnt}{RET}"))
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  sam_tib
}

bowtieProbeAlign = function(exe, fas, gen, dir, lan=NULL, run=FALSE,
                            verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'bowtieProbeAlign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} dir={dir}.{RET}"))
  
  stime <- system.time({
    fas_name <- fas %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    gen_name <- gen %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    gen_file <- gen %>% stringr::str_remove('.gz$')
    
    # Create Shell and Alignment directories::
    sh_dir <- file.path(dir, 'shells')
    al_dir <- file.path(dir, 'align')
    
    if (!dir.exists(sh_dir)) dir.create(sh_dir, recursive=TRUE)
    if (!dir.exists(al_dir)) dir.create(al_dir, recursive=TRUE)
    
    aln_ssh <- file.path(sh_dir, paste0('run_bow-',fas_name,'-',gen_name,'.sh') )
    aln_sam <- file.path(al_dir, paste0(fas_name,'-',gen_name,'.bowtie.sam.gz') )
    aln_cmd <- paste(exe, '-f -x',gen_file, '-U',fas, '| gzip -c ->',aln_sam, sep=' ')
    
    cat(glue::glue("[{funcTag}]:{TAB} Launching bowtie alignments: {fas_name} vs. {gen_name}...{RET}"))
    readr::write_lines(x=aln_cmd, file=aln_ssh, append=FALSE)
    Sys.chmod(paths=aln_ssh, mode="0777")
    if (!is.null(lan)) aln_ssh <- paste(lan,aln_ssh, sep=' ')
    if (run) base::system(aln_ssh)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (run) return(aln_sam)
  aln_ssh
}

# End of file

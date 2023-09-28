
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Probe Order Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Re-Order EPIC Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

formatReorderEPIC = function(tib1, tib2, dir, name='', verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'formatReorderEPIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  
  ann_csv <- file.path(dir, paste(name,'annotation.csv.gz', sep='.') )
  ord_csv <- file.path(dir, paste(name,'order.csv.gz', sep='.') )
  sum_csv <- file.path(dir, paste(name,'summary.csv.gz', sep='.') )
  
  ann_tib <- NULL
  ord_tib <- NULL
  sum_tib <- NULL
  
  # TBD:: Add chr,pos,strands,scores and other fields
  add_fields <- c('Forward_Sequence_Man','DesSeqN_Brac',
                  'Chromosome_Des','Coordinate_Des','IlmnID',
                  'Methyl_Allele_FR_Strand','Methyl_Allele_CO_Strand','Methyl_Allele_TB_Strand',
                  'Methyl_Next_Base','Infinium_Design',
                  'Min_Final_Score','Methyl_Final_Score','UnMethyl_Final_Score')
  
  ann_tib <- dplyr::bind_rows(
    prbs2order(tib1, addCols=add_fields, retGrp=1, blank=TRUE, verbose=opt$verbose),
    prbs2order(tib2, addCols=add_fields, retGrp=2, blank=TRUE, verbose=opt$verbose)
  ) %>% dplyr::select(-Next_Base) %>%
    dplyr::rename(Forward_Sequence=Forward_Sequence_Man, Design_Sequence=DesSeqN_Brac,
                      Chromosome=Chromosome_Des,Coordinate=Coordinate_Des,
                      Next_Base=Methyl_Next_Base) %>%
    dplyr::mutate(Color_Channel=dplyr::case_when(
      Next_Base=='A' | Next_Base=='T' ~ 'R',
      Next_Base=='C' | Next_Base=='G' ~ 'G',
      TRUE ~ NA_character_
    )) %>%
    dplyr::arrange(Assay_Design_Id)
  
  if (verbose>=vt+4) ann_tib %>% head(n=3) %>% as.data.frame() %>% print()
  
  ord_tib <- ann_tib %>% select(1:6)
  sum_tib <- order2stats(ann_tib)
  
  # Write Order::
  readr::write_csv(ann_tib,ann_csv)
  readr::write_csv(ord_tib,ord_csv)
  readr::write_csv(sum_tib,sum_csv)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  # sum_tib
  # ord_tib
  ann_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Probe Order Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

order2stats = function(tib, verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'order2stats'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  prb_seqs_tib <- NULL
  prb_seqs_tib <- tib %>% dplyr::select(AlleleA_Probe_Sequence, AlleleB_Probe_Sequence) %>% 
    tidyr::gather(Probe_Src, Probe_Seq) %>% 
    dplyr::filter(!is.na(Probe_Seq)) %>% 
    dplyr::filter(stringr::str_length(Probe_Seq)>0) %>%
    dplyr::mutate(Probe_Seq=stringr::str_to_upper(Probe_Seq))
  if (verbose>=vt+4) print(prb_seqs_tib)
  
  prb_sum_tib <- tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Assay_Design_Id, 1,2)) %>% 
    dplyr::group_by(Probe_Type) %>% dplyr::summarise(Type_Count=n()) %>% tidyr::spread(Probe_Type, Type_Count) %>%
    purrr::set_names(paste(names(.),"Count", sep='_'))
  if (verbose>=vt+4) print(prb_sum_tib)
  
  tot_ids_cnt <- tib %>% base::nrow()
  unq_ids_cnt <- tib %>% dplyr::distinct(Assay_Design_Id) %>% base::nrow()
  
  tot_seq_cnt <- prb_seqs_tib %>% base::nrow()
  unq_seq_cnt <- prb_seqs_tib %>% dplyr::distinct(Probe_Seq) %>% base::nrow()
  
  inf1_cnt <- tib %>% dplyr::filter(!is.na(AlleleB_Probe_Id) & stringr::str_length(AlleleB_Probe_Id)!=0) %>% base::nrow()
  inf2_cnt <- tib %>% dplyr::filter( is.na(AlleleB_Probe_Id) | stringr::str_length(AlleleB_Probe_Id)==0) %>% base::nrow()
  inf_tot_cnt <- (2*inf1_cnt) + inf2_cnt
  
  unq_cgn_cnt <- NA
  if(length(grep('IlmnID', names(tib)))>0) unq_cgn_cnt <- tib %>% dplyr::distinct(IlmnID) %>% base::nrow()
  
  stats <- tibble::tibble(tot_ids_cnt := tot_ids_cnt,
                          unq_ids_cnt := unq_ids_cnt,
                          inf1_cnt    := inf1_cnt,
                          int2_cnt    := inf2_cnt,
                          inf_tot_cnt := inf_tot_cnt,
                          tot_seq_cnt := tot_seq_cnt,
                          unq_seq_cnt := unq_seq_cnt,
                          unq_cgn_cnt := unq_cgn_cnt) %>%
    dplyr::bind_cols(prb_sum_tib)
  
  if (verbose>=vt) print(stats)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} DONE...{RET}{RET}"))
  
  
  stats
}

prbs2order = function(tib, addCols=NULL, retGrp=0, blank=FALSE, 
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'prbs2order'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  # Filter::
  #  - Flag Infinium I where PRB1_U != PRB1_M
  #  - Flag Infinium II probes with same color degenerate bases (e.g. w={a,t}, etc.)
  #  - Flag user selected probes
  #
  
  # First build Assay_Design_Id and then add replicates to ids
  tib <- tib %>% dplyr::group_by(Seq_ID_Uniq) %>% dplyr::add_count() %>% dplyr::mutate(Probe_Rep=n) %>% dplyr::ungroup()
  
  # Add CGN Replicates
  # selPrbA_tib %>% dplyr::group_by(Seq_ID_Uniq) %>% dplyr::summarise(Rep=n())
  
  sel_fields <- c("Assay_Design_Id","AlleleA_Probe_Id","AlleleA_Probe_Sequence",
                  "AlleleB_Probe_Id","AlleleB_Probe_Sequence","Normalization_Bin",
                  "Valid_Design", "Valid_Design_Bool")
  
  if (!is.null(addCols)) sel_fields <- c(sel_fields, addCols)
  # print(sel_fields)
  
  # Infinium I probes::
  inf1_tib <- tib %>% dplyr::mutate(Assay_Design_Id=paste(Seq_ID_Uniq,'1','T',Probe_Rep, sep=''),
                                    AlleleA_Probe_Id=paste(Assay_Design_Id,'A',sep='_'),
                                    AlleleA_Probe_Sequence=PRB1_U, 
                                    AlleleB_Probe_Id=paste(Assay_Design_Id,'B',sep='_'),
                                    AlleleB_Probe_Sequence=PRB1_M,
                                    Normalization_Bin=case_when(is.na(AlleleB_Probe_Sequence) | AlleleB_Probe_Sequence=='' ~ 'C',
                                                                NXB_M=='A' | NXB_M=='T' ~ 'A',
                                                                NXB_M=='a' | NXB_M=='t' ~ 'A',
                                                                NXB_M=='C' | NXB_M=='G' ~ 'B',
                                                                NXB_M=='c' | NXB_M=='g' ~ 'B',
                                                                TRUE ~ NA_character_),
                                    Valid_Design=case_when(
                                      stringr::str_to_upper(PRB1_U)==stringr::str_to_upper(PRB1_M) ~ paste('Fail_MatchInfI',NXB_M, sep='_'),
                                      TRUE ~ 'Pass'),
                                    Valid_Design_Bool=case_when(
                                      stringr::str_to_upper(PRB1_U)==stringr::str_to_upper(PRB1_M) ~ FALSE,
                                      TRUE ~ TRUE)
  ) %>% 
    dplyr::select(sel_fields, everything())
    # dplyr::select(dplyr::all_of(sel_fields), everything() )
  
  # Infinium II probes::
  inf2_tib <- tib %>% dplyr::mutate(Assay_Design_Id=paste(Seq_ID_Uniq,'2','T',Probe_Rep, sep=''),
                                    AlleleA_Probe_Id=paste(Assay_Design_Id,'A',sep='_'),
                                    AlleleA_Probe_Sequence=PRB2_D, 
                                    AlleleB_Probe_Id=NA_character_,
                                    AlleleB_Probe_Sequence=NA_character_,
                                    Normalization_Bin=case_when(is.na(AlleleB_Probe_Sequence) | AlleleB_Probe_Sequence=='' ~ 'C',
                                                                NXB_M=='A' | NXB_M=='T' ~ 'A',
                                                                NXB_M=='a' | NXB_M=='t' ~ 'A',
                                                                NXB_M=='C' | NXB_M=='G' ~ 'B',
                                                                NXB_M=='c' | NXB_M=='g' ~ 'B',
                                                                TRUE ~ NA_character_),
                                    Valid_Design=case_when(
                                      stringr::str_to_upper(CPN_D)=='R' | stringr::str_to_upper(CPN_D)=='Y' |
                                        stringr::str_to_upper(CPN_D)=='K' | stringr::str_to_upper(CPN_D)=='M' |
                                        stringr::str_to_upper(CPN_D)=='B' | stringr::str_to_upper(CPN_D)=='D' |
                                        stringr::str_to_upper(CPN_D)=='H' | stringr::str_to_upper(CPN_D)=='V' ~ paste('Pass_Channel',CPN_D, sep='_'),
                                      
                                      # stringr::str_to_upper(CPN_D)=='S' | stringr::str_to_upper(CPN_D)=='W' ~ 'Fail_Channel',
                                      # stringr::str_to_upper(PRB1_U)==stringr::str_to_upper(PRB1_M) ~ 'Fail_MatchInfI',
                                      TRUE ~ paste('Fail_Channel',CPN_D, sep='_')),
                                    Valid_Design_Bool=case_when(
                                      stringr::str_to_upper(CPN_D)=='R' | stringr::str_to_upper(CPN_D)=='Y' |
                                        stringr::str_to_upper(CPN_D)=='K' | stringr::str_to_upper(CPN_D)=='M' |
                                        stringr::str_to_upper(CPN_D)=='B' | stringr::str_to_upper(CPN_D)=='D' |
                                        stringr::str_to_upper(CPN_D)=='H' | stringr::str_to_upper(CPN_D)=='V' ~ TRUE,
                                      TRUE ~ FALSE)
  ) %>% 
    dplyr::select(sel_fields, everything())
    # dplyr::select(dplyr::all_of(sel_fields), everything() )
  
  if (blank) {
    inf2_tib <- inf2_tib %>%
      dplyr::mutate(AlleleB_Probe_Id=case_when(is.na(AlleleB_Probe_Id) ~ '', TRUE ~ AlleleB_Probe_Id),
                    AlleleB_Probe_Sequence=case_when(is.na(AlleleB_Probe_Sequence) ~ '', TRUE ~ AlleleB_Probe_Sequence))
  }
  
  ret <- NULL
  if (retGrp==1) {
    ret <- inf1_tib
  } else if (retGrp==2) {
    ret <- inf2_tib
  } else {
    ret$inf1 <- inf1_tib
    ret$inf2 <- inf2_tib
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))

  ret
}

# End of file

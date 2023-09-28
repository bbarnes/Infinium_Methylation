
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Basic Controls Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Global Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Non-CpG Adhoc I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadGRCm38_SNP = function(file, outDir,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadGRCm38_SNP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    
    probeType <- 'rs'
    genomeBuild <- 'GRCm38'
    dat_col <- c('Seq_ID','Probe_Type_Str','Design_Org_Seq',
                 'Genome_Build','Chromosome','Coordinate','AlleleA','AlleleB',
                 'IntA','IntB','Ann_Str')
    
    src_fon <- paste(genomeBuild,'improbeSourceInput.tsv.gz', sep='.')
    imp_fon <- paste(genomeBuild,'improbeDesignInput.tsv.gz', sep='.')
    
    src_tsv <- file.path(outDir, src_fon)
    imp_tsv <- file.path(outDir, imp_fon)
    
    src_tib <- 
      suppressMessages(suppressWarnings( readr::read_tsv(file, col_names=dat_col) )) %>% 
      dplyr::mutate(
        Seq_ID=stringr::str_remove(Seq_ID,'^snp_') %>% 
          stringr::str_replace_all('_','-'),
        Seq_ID=dplyr::case_when(
          stringr::str_starts(Seq_ID,'chr') ~ paste(genomeBuild,Seq_ID,sep='-'), 
          TRUE ~ Seq_ID),
        Iupac_Allele=mapDIs(paste0(AlleleA,AlleleB)),
        Pre_Seq=stringr::str_remove(Design_Org_Seq, '\\[.*$'),
        Pos_Seq=stringr::str_remove(Design_Org_Seq, '^.*\\]'),
        Beg_Nuc=stringr::str_sub(Pre_Seq,1,1),
        End_Nuc=dplyr::case_when(
          Beg_Nuc=='A' || Beg_Nuc=='a' || Beg_Nuc=='T' || Beg_Nuc=='t' ~ 'G',
          Beg_Nuc=='C' || Beg_Nuc=='c' || Beg_Nuc=='G' || Beg_Nuc=='g' ~ 'T',
          TRUE ~ 'N'
        ),
        Pos_Nuc=stringr::str_sub(Pos_Seq,1,1),
        Pos_Seq=stringr::str_sub(Pos_Seq,2),
        Sequence=stringr::str_to_upper(paste0(Pre_Seq,'[CG]',Pos_Seq,End_Nuc)),
        IUPAC_Sequence=paste0(Pre_Seq,'[',Iupac_Allele,Pos_Nuc,']',Pos_Seq,End_Nuc),
        Chromosome=stringr::str_remove(Chromosome,'^chr') %>% as.character(),
        Genome_Build=dplyr::case_when(
          Genome_Build=='hg18' ~ "GRCh36",
          Genome_Build=='hg19' ~ "GRCh37",
          Genome_Build=='hg38' ~ "GRCh38",
          Genome_Build=='mm10' ~ "GRCm38",
          TRUE ~ Genome_Build
        ),
        Probe_Type=probeType,
        CpG_Island="FALSE"
      ) %>% 
      dplyr::select(Seq_ID,Sequence,Genome_Build,Chromosome,Coordinate,
                    CpG_Island,IUPAC_Sequence,Probe_Type) %>%
      dplyr::arrange(Seq_ID)
    
    ret_tib <- src_tib %>% dplyr::select(Seq_ID:CpG_Island)
    
    #
    # Write Design Input and Annotation File::
    #
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Source (TSV)={imp_fon}...{RET}"))
    readr::write_tsv(src_tib,src_tsv)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Source (TSV)={imp_tsv}...{RET}"))
    readr::write_tsv(ret_tib,imp_tsv)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_files <- NULL
  ret_files[['src']] <- src_tsv
  ret_files[['imp']] <- imp_tsv
  
  ret_files
}

loadGRCh37_SNP = function(file, outDir,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadGRCh37_SNP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    
    probeType <- 'rs'
    genomeBuild <- 'GRCh37'
    dat_col <- c('Seq_ID_Long','Seq_ID','Design_Org_Seq',
                 'Genome_Build','Chromosome','Coordinate')
    
    src_fon <- paste(genomeBuild,'improbeSourceInput.tsv.gz', sep='.')
    imp_fon <- paste(genomeBuild,'improbeDesignInput.tsv.gz', sep='.')
    
    src_tsv <- file.path(outDir, src_fon)
    imp_tsv <- file.path(outDir, imp_fon)
    
    src_tib <- 
      suppressMessages(suppressWarnings( readr::read_tsv(file, col_names=dat_col) )) %>%
      tidyr::separate(Seq_ID_Long, into=c('Seq_ID','Allele_Count','SNP_Score','DiNuc'), sep=',') %>% 
      dplyr::mutate(
        Seq_ID=stringr::str_remove(Seq_ID,'^snp_') %>% 
          stringr::str_replace_all('_','-'),
        Seq_ID=dplyr::case_when(
          stringr::str_starts(Seq_ID,'chr') ~ paste(genomeBuild,Seq_ID,sep='-'), 
          TRUE ~ Seq_ID),
        Ref_Nuc=stringr::str_remove(Design_Org_Seq,'^.*\\[') %>%
          stringr::str_remove('[A-Za-z]\\].*$'),
        Iupac_Allele=mapDIs(DiNuc),
        Iupac_Allele=dplyr::case_when(
          is.na(Iupac_Allele) & Ref_Nuc=='A' ~ 'R',
          is.na(Iupac_Allele) & Ref_Nuc=='G' ~ 'R',
          is.na(Iupac_Allele) & Ref_Nuc=='T' ~ 'Y',
          is.na(Iupac_Allele) & Ref_Nuc=='C' ~ 'Y',
          TRUE ~ Iupac_Allele
        ),
        Pre_Seq=stringr::str_remove(Design_Org_Seq, '\\[.*$'),
        Pos_Seq=stringr::str_remove(Design_Org_Seq, '^.*\\]'),
        Beg_Nuc=stringr::str_sub(Pre_Seq,1,1),
        End_Nuc=dplyr::case_when(
          Beg_Nuc=='A' || Beg_Nuc=='a' || Beg_Nuc=='T' || Beg_Nuc=='t' ~ 'G',
          Beg_Nuc=='C' || Beg_Nuc=='c' || Beg_Nuc=='G' || Beg_Nuc=='g' ~ 'T',
          TRUE ~ 'N'
        ),
        Pos_Nuc=stringr::str_sub(Pos_Seq,1,1),
        Pos_Seq=stringr::str_sub(Pos_Seq,2),
        Sequence=stringr::str_to_upper(paste0(Pre_Seq,'[CG]',Pos_Seq,End_Nuc)),
        IUPAC_Sequence=paste0(Pre_Seq,'[',Iupac_Allele,Pos_Nuc,']',Pos_Seq,End_Nuc),
        Chromosome=stringr::str_remove(Chromosome,'^chr') %>% as.character(),
        Genome_Build_Str=dplyr::case_when(
          Genome_Build==18 ~ "GRCh36",
          Genome_Build==36 ~ "GRCh36",
          Genome_Build==19 ~ "GRCh37",
          Genome_Build==37 ~ "GRCh37",
          Genome_Build==38 ~ "GRCh38",
          Genome_Build==10 ~ "GRCm38",
          TRUE ~ NA_character_
        ),
        Genome_Build=Genome_Build_Str,
        Probe_Type=probeType,
        CpG_Island="FALSE"
      ) %>% 
      dplyr::select(Seq_ID,Sequence,Genome_Build,Chromosome,Coordinate,
                    CpG_Island,IUPAC_Sequence,Probe_Type) %>%
      dplyr::arrange(Seq_ID)
    
    ret_tib <- src_tib %>% dplyr::select(Seq_ID:CpG_Island)
    
    #
    # Write Design Input and Annotation File::
    #
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Source (TSV)={imp_fon}...{RET}"))
    readr::write_tsv(src_tib,src_tsv)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Source (TSV)={imp_tsv}...{RET}"))
    readr::write_tsv(ret_tib,imp_tsv)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_files <- NULL
  ret_files[['src']] <- src_tsv
  ret_files[['imp']] <- imp_tsv
  
  ret_files
}


# End of file

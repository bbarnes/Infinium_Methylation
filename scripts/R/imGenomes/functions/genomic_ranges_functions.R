
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Genomic Range Methods:: Intersection
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersect_GRS = function(ref, can,
                         
                         ref_key=NULL,
                         ref_col=NULL,
                         ref_prefix=NULL,
                         ref_red=TRUE,
                         
                         can_key=NULL,
                         can_col=NULL,
                         can_prefix=NULL, 
                         
                         ret_off = FALSE,
                         
                         out_dir,
                         run_tag,
                         pre_tag = NULL,
                         
                         reload     = 0,
                         reload_min = 2,
                         ret_data   = FALSE,
                         parallel   = FALSE,
                         
                         vb=0, vt=3, tc=1, tt=NULL,
                         fun_tag='intersect_GRS') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ),
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )

  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}      ret_off = '{ret_off}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_tag = '{out_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, aux_csv, out_csv, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Intersect Ranges::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( vb >= vt+2 ) cat(glue::glue("{mssg} Finding Overlaps...{RET}"))
    map_tib <- 
      GenomicRanges::findOverlaps(can,ref, ignore.strand=TRUE) %>%
      as.data.frame() %>% tibble::as_tibble()
    map_key <- glue::glue("map-tib")
    map_cnt <- print_tib( map_tib, fun_tag = fun_tag, name = map_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Extract Candidate Matches::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( vb >= vt+2 )
      cat(glue::glue("{mssg} Extracting Candidate Data...{RET}"))
    
    can_tib <- as.data.frame(can) %>% tibble::as_tibble(rownames=can_key)
    
    if (!is.null(can_col) && length(can_col)!=0) can_tib <- can_tib %>%
      dplyr::select(dplyr::all_of(can_col))
    
    if (!is.null(can_prefix)) can_tib <- can_tib %>% 
      purrr::set_names(paste(can_prefix,names(.), sep="_"))
    
    can_key <- glue::glue("can-tib")
    can_cnt <- print_tib( can_tib, fun_tag = fun_tag, name = can_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Extract Reference Matches::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( vb >= vt+2 )
      cat(glue::glue("{mssg} Extracting Reference Data...{RET}"))
    
    ref_tib <- as.data.frame(ref) %>% tibble::as_tibble(rownames=ref_key)
    
    if (ref_red) ref_tib <- ref_tib %>% 
      dplyr::mutate(seqnames=as.character(seqnames),
                    strand=as.character(strand)) %>%
      dplyr::rename(chr=seqnames,
                    pos=start,
                    srd=strand)
    
    if (!is.null(ref_col) && length(ref_col)!=0) ref_tib <- ref_tib %>%
      dplyr::select(dplyr::all_of(ref_col))
    
    if (!is.null(ref_prefix)) ref_tib <- ref_tib %>% 
      purrr::set_names(paste(ref_prefix,names(.), sep="_"))
    
    if (FALSE) {
      ref_tib <- ref %>% as.data.frame() %>%
        tibble::as_tibble(rownames=ref_key) %>%
        dplyr::mutate(seqnames=as.character(seqnames),
                      strand=as.character(strand)) %>%
        dplyr::rename(chr=seqnames,
                      pos=start,
                      top_srd=strand)
      
      if (!is.null(ref_prefix)) ref_tib <- ref_tib %>% 
          purrr::set_names(paste(ref_prefix,names(.), sep="_"))
    }
    
    ref_str <- glue::glue("ref-tib")
    ref_cnt <- print_tib( ref_tib, fun_tag = fun_tag, name = ref_str, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Bind Reference/Candidate Matches::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( vb >= vt+2 )
      cat(glue::glue("{mssg} Binding Candidate/Reference Data...{RET}"))
    
    # print(can_tib[map_tib$queryHits,] )
    
    ret_tib <- dplyr::bind_cols(
      can_tib[map_tib$queryHits, ],
      ref_tib[map_tib$subjectHits,]
    )
    print(ret_tib)
    
    # NOTE:: Returning the unmatched candidates does NOT currently work!
    #
    # if ( ret_off ) ret_tib <- ret_tib %>% 
    #   dplyr::bind_rows( can_tib[!map_tib$queryHits, ] )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Generate Reference Coverage Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # if ( vb >= vt+2 )
    #   cat(glue::glue("{mssg} Merging back full Reference...{RET}"))
    # ref_cov_tib <- ret_tib %>% dplyr::right_join( ref_tib )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv",
    #                        done = TRUE, write_spec = TRUE, append = FALSE,
    #                        fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  # ref_cov_tib
  ret_tib
}

intersect_GRS_old = function(ref,can,
                             ref_key=NULL,ref_col=NULL,ref_prefix=NULL,ref_red=TRUE,
                             can_key=NULL,can_col=NULL,can_prefix=NULL, 
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersect_GRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    #
    # Perform Intersection::
    #
    map_tib <- 
      GenomicRanges::findOverlaps(can,ref, ignore.strand=TRUE) %>%
      as.data.frame() %>% tibble::as_tibble()
    map_key <- glue::glue("map_tib({funcTag})")
    map_cnt <- print_tib(map_tib,funcTag, verbose,vt+4,tc, n=map_key)
    
    #
    # Return Tibble from Genomic Range: Candidate
    #
    can_tib <- as.data.frame(can) %>% tibble::as_tibble(rownames=can_key)
    
    if (!is.null(can_col) && length(can_col)!=0) can_tib <- can_tib %>%
      dplyr::select(dplyr::all_of(can_col))
    
    if (!is.null(can_prefix)) can_tib <- can_tib %>% 
      purrr::set_names(paste(can_prefix,names(.), sep="_"))
    
    can_key <- glue::glue("can_tib({funcTag})")
    can_cnt <- print_tib(can_tib,funcTag, verbose,vt+4,tc, n=can_key)
    
    #
    # Return Tibble from Genomic Range: Reference
    #
    ref_tib <- as.data.frame(ref) %>% tibble::as_tibble(rownames=ref_key)
    
    if (ref_red) ref_tib <- ref_tib %>% 
      dplyr::mutate(seqnames=as.character(seqnames),
                    strand=as.character(strand)) %>%
      dplyr::rename(chr=seqnames,
                    pos=start,
                    srd=strand)
    
    if (!is.null(ref_col) && length(ref_col)!=0) ref_tib <- ref_tib %>%
      dplyr::select(dplyr::all_of(ref_col))
    
    if (!is.null(ref_prefix)) ref_tib <- ref_tib %>% 
      purrr::set_names(paste(ref_prefix,names(.), sep="_"))
    
    if (FALSE) {
      ref_tib <- ref %>% as.data.frame() %>%
        tibble::as_tibble(rownames=ref_key) %>%
        dplyr::mutate(seqnames=as.character(seqnames),
                      strand=as.character(strand)) %>%
        dplyr::rename(chr=seqnames,
                      pos=start,
                      top_srd=strand)
      
      if (!is.null(ref_prefix)) ref_tib <- ref_tib %>% 
          purrr::set_names(paste(ref_prefix,names(.), sep="_"))
    }
    
    ref_key <- glue::glue("ref_tib({funcTag})")
    ref_cnt <- print_tib(ref_tib,funcTag, verbose,vt+4,tc, n=ref_key)
    
    #
    # Bind Candidate and Reference Tibs::
    #
    ret_tib <- dplyr::bind_cols(
      can_tib[map_tib$queryHits, ],
      ref_tib[map_tib$subjectHits,]
    )
    # This will be accounted for in the summary step::
    #   ret_tib <- ret_tib %>% dplyr::right_join( ref_tib )
    
    ret_key <- glue::glue("ret_FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Intersection Method:: 
#                                 OLD CODE!!!
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersectGranges = function(man,ref,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersectGranges'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (! purrr::is_list(ref)) {
      map_tib <- 
        GenomicRanges::findOverlaps(man,ref, ignore.strand=TRUE) %>%
        as.data.frame() %>% tibble::as_tibble()
      
      mani_tib <- man %>% as.data.frame() %>%
        rownames_to_column(var='Seq_ID') %>% tibble::as_tibble() 
      mani_len <- mani_tib %>% base::nrow()
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} mani_tib({mani_len})={RET}"))
        print(mani_tib)
      }
      
      gene_tib <- ref %>% as.data.frame() %>%
        rownames_to_column(var='Name') %>% tibble::as_tibble() # %>% 
      # purrr::set_names(c('Gene','chrom','chromStart','chromEnd','chromLength','chromStrand'))
      gene_len <- gene_tib %>% base::nrow()
      
      #
      # TBD::
      # TBD:: Lame way to fix this; the code commented out above:
      # TBD::
      #
      colnames(gene_tib)[1] <- "Gene"
      colnames(gene_tib)[2] <- "chrom"
      colnames(gene_tib)[3] <- "chromStart"
      colnames(gene_tib)[4] <- "chromEnd"
      colnames(gene_tib)[5] <- "chromLength"
      colnames(gene_tib)[6] <- "chromStrand"
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} gene_tib({gene_len})={RET}"))
        print(gene_tib)
      }
      
      #
      # Last change:: # dplyr::select(Seq_ID),
      #
      cur_tib <- dplyr::bind_cols(
        mani_tib[map_tib$queryHits, ], # %>% dplyr::select(Seq_ID),
        gene_tib[map_tib$subjectHits,] ) # %>% dplyr::mutate(Feature=feat_key)
      cur_len <- cur_tib %>% base::nrow()
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} cur_tib({cur_len})={RET}"))
        print(cur_tib)
      }
      ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
    } else {
      for (feat_key in names(ref)) {
        if (verbose>=vt)
          cat(glue::glue("[{funcTag}]:{tabsStr} GRange Overlap; feature={feat_key}...{RET}"))
        
        map_tib <- 
          GenomicRanges::findOverlaps(man,ref[[feat_key]], ignore.strand=TRUE) %>%
          as.data.frame() %>% tibble::as_tibble()
        map_len <- base::nrow(map_tib)
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; man={RET}"))
          print(man)
          
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; ref={RET}"))
          print(ref[[feat_key]])
          
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; map_tib={map_len}{RET}"))
          print(map_tib)
          
          cat(glue::glue("[{funcTag}]:{tabsStr}{RET}{RET}"))
        }
        
        mani_tib <- man %>% as.data.frame() %>%
          rownames_to_column(var='Seq_ID') %>% tibble::as_tibble() 
        mani_len <- mani_tib %>% base::nrow()
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; mani_tib({mani_len})={RET}"))
          print(mani_tib)
        }
        
        gene_tib <- ref[[feat_key]] %>% as.data.frame() %>% 
          rownames_to_column(var='Name') %>% tibble::as_tibble() %>% 
          purrr::set_names(c('Gene','chrom','chromStart','chromEnd','chromLength','chromStrand'))
        gene_len <- gene_tib %>% base::nrow()
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; gene_tib({gene_len})={RET}"))
          print(gene_tib)
        }
        
        cur_tib <- dplyr::bind_cols(
          # tib[map_tib$queryHits, ] %>% dplyr::select(Seq_ID),
          mani_tib[map_tib$queryHits, ] %>% dplyr::select(Seq_ID),
          gene_tib[map_tib$subjectHits,] ) %>%
          dplyr::mutate(Feature=feat_key)
        cur_len <- cur_tib %>% base::nrow()
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; cur_tib({cur_len})={RET}"))
          print(cur_tib)
        }
        ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
      }
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}



# End of file

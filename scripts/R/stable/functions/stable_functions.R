
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               TripleCrown::
#                         CG# Database Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Tidyverse Core Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("tidyverse",  quietly = TRUE) ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Common Human Abbreviations for String Variables::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
TAB2 <- "\t\t"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Basic R-Squared Method::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

r2_matrix = function( tib,
                      # file,
                      id,
                      colsA,
                      colsB,
                      
                      out_dir,
                      run_tag,
                      
                      reload     = 0,
                      reload_min = 2,
                      reload_pre = NULL,
                      
                      ret_data   = FALSE,
                      parallel   = FALSE,
                      
                      vb=0, vt=3, tc=1, tt=NULL,
                      fun_tag='r2_matrix')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt
  p1  <- vb > vt + 1
  p1  <- vb > vt + 2
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    # cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
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
  
  # errs_mssg <- glue::glue("File file='{file}' does not exist")
  # if ( !file.exists( file) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    for ( colA in colsA ) {
      
      # for ( colB in colsB ) {
      # }
      
      matA = tib %>% 
        dplyr::select( dplyr::all_of( c(id,colA) ) ) %>%
        tibble::column_to_rownames( id ) %>%
        as.matrix()
      
      matA %>% head() %>% print()
      
      matB = tib %>% 
        dplyr::select( dplyr::all_of( c(id,colsB) ) ) %>%
        tibble::column_to_rownames( id ) %>%
        as.matrix()
      
      matB %>% head() %>% print()
      
      ret_tib <- ret_tib %>% dplyr::bind_rows(
        stats::cor( x = matA, y = matB, 
                    method = "pearson", 
                    use = "pairwise.complete.obs" ) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column( var = "Meteric" )
      )
    }
    
    # warn_mssg <- glue::glue("WARN_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # wflag <- FALSE
    # 
    # errs_mssg <- glue::glue("ERROR_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    # if ( eflag ) return( NULL )
    
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
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   Evonik Chromosome/Contig BSP Analysis::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NOTE:: Will need to make a function that is used by lapply()
#
# bsp_tib/can = bsp_source_all_alleleA_list[[ref_key]]
# ann_grs/ref = pos_reg_all_grs[[ref_key]]
# fai_tib/fai = fai_source_list[[ref_key]]
# run_tag = ref_key
#
esitmate_core_contigs = function( can_tib,
                                  ref_grs,
                                  fai_tib,
                                  
                                  max_can = 0,
                                  max_ref = 0,
                                  max_fai = 0,
                                  
                                  out_dir,
                                  run_tag,
                                  
                                  reload     = 0,
                                  reload_min = 2,
                                  reload_pre = NULL,
                                  
                                  ret_data   = FALSE,
                                  parallel   = FALSE,
                                  
                                  vb=0, vt=3, tc=1, tt=NULL,
                                  fun_tag='esitmate_core_contigs' )
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt
  p1  <- vb > vt + 1
  p1  <- vb > vt + 2
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      max_can = '{max_can}'.{RET}"))
    cat(glue::glue("{mssg}      max_ref = '{max_ref}'.{RET}"))
    cat(glue::glue("{mssg}      max_fai = '{max_fai}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
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
  
  # errs_mssg <- glue::glue("File file='{file}' does not exist")
  # if ( !file.exists( file) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( max_can > 0 ) can_tib <- can_tib %>% head( n=max_can )
  if ( max_ref > 0 ) ref_grs <- can_tib %>% head( n=max_ref )
  if ( max_fai > 0 ) fai_tib <- can_tib %>% head( n=max_fai )
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                         Genomic Ranges:: All BSP
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( p1 ) cat( glue::glue("{pmssg} Joining Alignment with Source_Types; ",
                              "Ref Key='{run_tag}'...{RET}") )
    
    # A tibble: 2,755,670 × 9:: run_tag="galGal5"
    can_mat_grs <- NULL
    can_mat_grs <- 
      GenomicRanges::GRanges(
        seqnames = Rle( can_tib$Chr ),
        # strand = Rle(can_tib$srd),
        
        Ilmn_ID = can_tib$Ilmn_ID,
        Tag     = can_tib$Tag,
        Mis_Cnt = can_tib$Mis_Cnt,
        Map_Idx = can_tib$Map_Idx,
        # NOTE:: Map_Idx2 is needed for uniqueness::
        #  - A tibble: 43 × 9 out of A tibble: 35,331,564 × 9 have Map_Idx2 != 1
        Map_Idx2 = can_tib$Map_Idx2,
        Map_Count = can_tib$Map_Count,
        Species = can_tib$Ref_Source,
        CpG_Pos = can_tib$CpG_Pos,
        
        IRanges(start = can_tib$CpG_Pos - 1, 
                end = can_tib$CpG_Pos + 1,
                names = paste( can_tib$Ilmn_ID,
                               can_tib$Chr,
                               can_tib$CpG_Pos,
                               can_tib$Map_Idx2,
                               sep="-") )
      )
    
    if ( p4 ) cat( glue::glue("{pmssg} can_mat_grs::{RET}") )
    if ( p4 ) can_mat_grs %>% print()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                 Intersection:: BSP vs. Target Regions::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    bsp_reg_int_tib <- NULL
    bsp_reg_int_tib <- 
      intersect_GRS( ref = ref_grs, 
                     can = can_mat_grs, 
                     ref_prefix = "ref", 
                     can_prefix = "can", 
                     out_dir = file.path( out_dir,"intersect/bsp-vs-reg/on"), 
                     run_tag = run_tag, 
                     reload = reload ,
                     reload_min = reload_min,
                     parallel = parallel,
                     ret_data = FALSE,
                     vb=vb,vt=vt+1,tc=tc,tt=tt )
    
    out_key <- glue::glue("bsp_reg_int")
    out_cnt <- print_tib( bsp_reg_int_tib, fun_tag=fun_tag, name=out_key,
                          vb=vb,vt=vt+4,tc=tc+1 )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #          Reverse look up any probes that did not get matched 
    #                           to target regions::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    bsp_source_mis_alleleA_tib <- NULL
    bsp_source_mis_alleleA_tib <- 
      can_tib %>% 
      dplyr::anti_join( bsp_reg_int_tib, 
                        by=c("Chr"="can_seqnames", "CpG_Pos"="can_CpG_Pos", 
                             "Ilmn_ID"="can_Ilmn_ID", "Map_Idx2"="can_Map_Idx2",
                             "Tag"="can_Tag", "Ref_Source"="can_Species" )
      )
    
    can_mis_grs <- NULL
    can_mis_grs <- 
      GenomicRanges::GRanges(
        seqnames = Rle( bsp_source_mis_alleleA_tib$Chr ),
        # strand = Rle(bsp_source_mis_alleleA_tib$srd),
        
        # Tag   Mis_Cnt Map_Idx Map_Idx2 Map_Count
        Ilmn_ID = bsp_source_mis_alleleA_tib$Ilmn_ID,
        Tag     = bsp_source_mis_alleleA_tib$Tag,
        Mis_Cnt = bsp_source_mis_alleleA_tib$Mis_Cnt,
        Map_Idx = bsp_source_mis_alleleA_tib$Map_Idx,
        # NOTE:: Map_Idx2 is needed for uniqueness::
        #  - A tibble: 43 × 9 out of A tibble: 35,331,564 × 9 have Map_Idx2 != 1
        Map_Idx2 = bsp_source_mis_alleleA_tib$Map_Idx2,
        Map_Count = bsp_source_mis_alleleA_tib$Map_Count,
        Species = bsp_source_mis_alleleA_tib$Ref_Source,
        CpG_Pos = bsp_source_mis_alleleA_tib$CpG_Pos,
        
        IRanges(start = bsp_source_mis_alleleA_tib$CpG_Pos - 1, 
                end = bsp_source_mis_alleleA_tib$CpG_Pos + 1,
                names = paste( bsp_source_mis_alleleA_tib$Ilmn_ID,
                               bsp_source_mis_alleleA_tib$Chr,
                               bsp_source_mis_alleleA_tib$CpG_Pos,
                               bsp_source_mis_alleleA_tib$Map_Idx2,
                               sep="-") )
      )
    
    if ( p4 ) cat( glue::glue("{pmssg} can_mis_grs::{RET}") )
    if ( p4 ) can_mis_grs %>% print()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                  Build fai_grs from fai_tab to look like::
    #                           pos_reg_all_grs[[run_tag]]
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    fai_grs <- NULL
    fai_grs <- 
      GenomicRanges::GRanges(
        seqnames = Rle( fai_tib$chr %>% stringr::str_remove("^chr") ),
        # strand = Rle(fai_tib$srd),
        
        Species = fai_tib$source,
        Source_Type = "Off",
        Group = "Off",
        
        IRanges(start = 1, 
                end = fai_tib$Chrom_Length,
                names=paste( fai_tib$source,
                             fai_tib$chr,
                             fai_tib$Chrom_Length,
                             sep="-") )
      )
    
    if ( p4 ) cat( glue::glue("{pmssg} fai_grs::{RET}") )
    if ( p4 ) fai_grs %>% print()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                 Intersection:: BSP vs. Target Regions::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    bsp_reg_off_tib <- NULL
    bsp_reg_off_tib <- 
      intersect_GRS( ref = fai_grs, 
                     can = can_mis_grs, 
                     ref_prefix = "ref", 
                     can_prefix = "can", 
                     out_dir = file.path( out_dir,"intersect/bsp-vs-reg/off"), 
                     run_tag = run_tag, 
                     reload = reload,
                     reload_min = reload_min,
                     parallel = parallel,
                     ret_data = FALSE,
                     vb=vb,vt=vt+1,tc=tc,tt=tt )
    
    out_key <- glue::glue("bsp_reg_off")
    out_cnt <- print_tib( bsp_reg_off_tib, fun_tag=fun_tag, name=out_key,
                          vb=vb,vt=vt+4,tc=tc+1 )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                        Merge Results Back Together::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    bsp_reg_all_tib <- NULL
    bsp_reg_all_tib <- dplyr::bind_rows( bsp_reg_int_tib, 
                                         bsp_reg_off_tib ) %>%
      dplyr::distinct( can_seqnames,can_CpG_Pos,can_Ilmn_ID,
                       can_Map_Idx2,can_Tag,can_Species, .keep_all = TRUE ) %>%
      dplyr::arrange( can_Ilmn_ID )
    
    out_key <- glue::glue("bsp_reg_all")
    out_cnt <- print_tib( bsp_reg_all_tib, fun_tag=fun_tag, name=out_key,
                          vb=vb,vt=vt+4,tc=tc+1 )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #               Validate ALl BSP Alignments are Accounted for::
    #
    #  NOTE:: Also check the number of chromosome not covered (not important...)
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    mis_ref_cnt <- can_tib %>% dplyr::anti_join( 
      bsp_reg_all_tib, 
      by=c("Chr"="can_seqnames", "CpG_Pos"="can_CpG_Pos", 
           "Ilmn_ID"="can_Ilmn_ID", "Map_Idx2"="can_Map_Idx2",
           "Tag"="can_Tag", "Ref_Source"="can_Species" ) ) %>% base::nrow()
    
    mis_can_cnt <- bsp_reg_all_tib %>% dplyr::anti_join( 
      can_tib, 
      by=c("can_seqnames"="Chr", "can_CpG_Pos"="CpG_Pos", 
           "can_Ilmn_ID"="Ilmn_ID", "can_Map_Idx2"="Map_Idx2",
           "can_Tag"="Tag", "can_Species"="Ref_Source" ) ) %>% base::nrow()
    
    if ( mis_can_cnt != 0 ) {
      bsp_reg_all_tib %>% dplyr::anti_join( 
        can_tib, 
        by=c("can_seqnames"="Chr", "can_CpG_Pos"="CpG_Pos", 
             "can_Ilmn_ID"="Ilmn_ID", "can_Map_Idx2"="Map_Idx2",
             "can_Tag"="Tag", "can_Species"="Ref_Source" ) ) %>%
        print (n=mis_can_cnt )
      
      can_tib %>% dplyr::anti_join( 
        bsp_reg_all_tib, 
        by=c("Chr"="can_seqnames", "CpG_Pos"="can_CpG_Pos", 
             "Ilmn_ID"="can_Ilmn_ID", "Map_Idx2"="can_Map_Idx2",
             "Tag"="can_Tag", "Ref_Source"="can_Species" ) ) %>%
        print( n=mis_ref_cnt )
      
      eflag <- TRUE
      errs_mssg <- glue::glue("Failed to recover all BSP Alignments: ",
                              "mis_can_cnt::{mis_can_cnt} != 0, ",
                              "mis_ref_cnt={mis_ref_cnt}")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    if ( p1 ) cat( glue::glue("{pmssg} Validated BSP::{RET}",
                              "{pmssg}{TAB} mis_ref_cnt='{mis_ref_cnt}'{RET}",
                              "{pmssg}{TAB} mis_can_cnt='{mis_can_cnt}'{RET2}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                             Quality Control::
    #                   Validate Basic Intersection Stats::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # TBD:: Fix Genomic Ranges so we don't have one-off issues...
    # Coordinate Validation
    #  Example::  bsp_reg_all_tibs[[1]] %>% dplyr::filter( can_CpG_Pos < ref_pos-1 | can_CpG_Pos > ref_end+1 )
    #
    off_one_tib <- bsp_reg_all_tib %>% 
      dplyr::filter( can_CpG_Pos < ref_pos-1 | can_CpG_Pos > ref_end+1 )
    off_one_cnt <- off_one_tib %>% base::nrow()
    
    if ( off_one_cnt != 0 ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("off_one_cnt::{off_one_cnt} != 0")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    # Species Validation::
    #  Example::  bsp_reg_all_tibs[[3]] %>% dplyr::filter( can_Species != ref_Species )
    #
    spe_mis_tib <- bsp_reg_all_tib %>% 
      dplyr::filter( can_Species != ref_Species )
    spe_mis_cnt <- spe_mis_tib %>% base::nrow()
    
    if ( spe_mis_cnt != 0 ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("spe_mis_cnt::{spe_mis_cnt} != 0")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    # TBD:: Chromosome Validation
    #  Example::  bsp_reg_all_tibs[[1]] %>% dplyr::mutate( can_seqnames=can_seqnames %>% as.character() ) %>% dplyr::filter( can_seqnames != ref_chr )
    #
    chr_mis_tib <- bsp_reg_all_tib %>% 
      dplyr::mutate( can_seqnames=can_seqnames %>% as.character() ) %>% 
      dplyr::filter( can_seqnames != ref_chr )
    chr_mis_cnt <- chr_mis_tib %>% base::nrow()
    
    if ( chr_mis_cnt != 0 ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("chr_mis_cnt::{chr_mis_cnt} != 0")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    # TBD:: Unique Probe_ID/Tag Pair::
    #  Example::  bsp_reg_all_tibs[[1]] %>% dplyr::distinct( can_Ilmn_ID,ref_Species,can_Tag ) %>% dplyr::add_count( can_Ilmn_ID, name="Dup_Cnt") %>% dplyr::filter( Dup_Cnt != 1 )
    #
    tag_mis_tib <- bsp_reg_all_tib %>% 
      dplyr::distinct( can_Ilmn_ID,ref_Species,can_Tag ) %>% 
      dplyr::add_count( can_Ilmn_ID, name="Dup_Cnt") %>% 
      dplyr::filter( Dup_Cnt != 1 )
    tag_mis_cnt <- tag_mis_tib %>% base::nrow()
    
    if ( tag_mis_cnt != 0 ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("tag_mis_cnt::{tag_mis_cnt} != 0")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    if ( p1 ) cat( glue::glue("{pmssg} Validated Intersection Stats::{RET}",
                              "{pmssg}{TAB} off_one_cnt='{off_one_cnt}'{RET}",
                              "{pmssg}{TAB} spe_mis_cnt='{spe_mis_cnt}'{RET}",
                              "{pmssg}{TAB} tag_mis_cnt='{tag_mis_cnt}'{RET2}") )
    
    
    #
    # TBD:: Add More Summary Stats::
    #  Examples Below::
    #    bsp_reg_all_tibs[[1]] %>% dplyr::group_by(can_Mis_Cnt) %>% dplyr::summarise( Count=n(), .groups = "drop" )
    #    contig_ret_tib[[ref_key]] %>% dplyr::group_by( can_Tag,can_Mis_Cnt ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
    #
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                       Add Back Chromosome Length::
    #             Calculate chrDb (Chromosome Length Delta Beta)
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Clean up names as well...
    #
    bsp_reg_sel_tib <- NULL
    bsp_reg_sel_tib <- 
      bsp_reg_all_tib %>% 
      dplyr::select( can_Ilmn_ID, 
                     can_Tag,can_Mis_Cnt,can_Map_Idx,can_Map_Idx2,can_Map_Count, 
                     ref_chr,can_CpG_Pos,ref_srd,
                     ref_Species,ref_Source_Type,ref_Group ) %>%
      dplyr::rename( Ilmn_ID=can_Ilmn_ID, Tag=can_Tag, Mis_Cnt=can_Mis_Cnt, 
                     Map_Idx=can_Map_Idx, Map_Idx2=can_Map_Idx2, 
                     Map_Count=can_Map_Count, Chromosome=ref_chr, 
                     CpG_Pos=can_CpG_Pos, Strand=ref_srd, 
                     Species=ref_Species, Source_Type=ref_Source_Type, 
                     Group=ref_Group ) %>%
      # dplyr::inner_join( fai_source_list[[ref_key]] %>%
      dplyr::inner_join( fai_tib %>%
                           dplyr::mutate( chr=chr %>% stringr::str_remove("^chr") ), 
                         by=c("Species"="source","Chromosome"="chr") ) %>%
      dplyr::arrange( Ilmn_ID,-Chrom_Length,Mis_Cnt ) %>%
      dplyr::group_by( Ilmn_ID ) %>%
      dplyr::mutate( Order_Rank = dplyr::row_number() ) %>%
      dplyr::ungroup()
    
    bsp_reg_bad_tib <- NULL
    bsp_reg_bad_tib <- bsp_reg_sel_tib %>% 
      dplyr::filter( Order_Rank == 1 ) %>%
      dplyr::filter( Mis_Cnt != 0 | Source_Type == "Off" )
    
    #
    # Remove Bad Top Probes and all Others
    #  - Calculate summary stats
    #  - Select only Top Probes
    #  - Join Top with Summary
    #
    sel_all_tib <- bsp_reg_sel_tib %>% 
      dplyr::anti_join( bsp_reg_bad_tib, by=c("Ilmn_ID") ) 
    top_all_tib <- sel_all_tib %>% dplyr::filter( Order_Rank == 1 )
    
    #
    # Sanity Check top_all_tib == sel_all_tib %>% dplyr::distinct( Ilmn_ID )
    #
    sel_unq_cnt <- sel_all_tib %>% dplyr::distinct( Ilmn_ID ) %>% base::nrow()
    top_unq_cnt <- top_all_tib %>% base::nrow()
    if ( sel_unq_cnt != top_unq_cnt ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("Failed to recover all BSP Alignments: ",
                              "sel_unq_cnt={sel_unq_cnt} != 0, ",
                              "top_unq_cnt={top_unq_cnt}")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    qc_tib <- top_all_tib %>% dplyr::add_count( Ilmn_ID, name="Dup_Cnt" ) %>% 
      dplyr::filter( Dup_Cnt != 1 )
    qc_cnt <- qc_tib %>% base::nrow()
    
    if ( p1 ) cat( glue::glue("{pmssg} Quality Control Check::{RET}",
                              "{pmssg}{TAB} qc_cnt='{qc_cnt}'{RET2}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                        Calculate Summary Stats::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    sum_all_tib <- NULL
    sum_all_tib <- sel_all_tib %>% 
      dplyr::group_by( Ilmn_ID,Tag,Species,Source_Type ) %>% 
      dplyr::summarise( Cnt=n(),
                        Sum = sum( Chrom_Length, na.rm = TRUE ),
                        Min = min( Mis_Cnt, na.rm = TRUE ),
                        # Avg = mean( Mis_Cnt, na.rm = TRUE ),
                        Max  = max( Mis_Cnt, na.rm = TRUE ),
                        .groups = "drop" ) %>% 
      tidyr::pivot_wider( id_cols = c(Ilmn_ID,Species,Tag), 
                          names_from = c(Source_Type), 
                          values_from = c(Cnt,Sum,Min,Max), 
                          values_fill = 0 )
    
    sum_all_tib <- sum_all_tib %>%
      dplyr::mutate( Sum_All_Chrom_Length = sum_all_tib %>% 
                       dplyr::select( dplyr::starts_with("Sum_") ) %>% 
                       as.matrix() %>% matrixStats::rowSums2() %>% as.vector() ) %>%
      dplyr::left_join( sel_all_tib %>% dplyr::select( Ilmn_ID,Chrom_Length ),
                        by=c("Ilmn_ID") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                   END::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- bsp_reg_sel_tib
    
    if ( FALSE ) {
      warn_mssg <- glue::glue("WARN_MESSAGE")
      if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
      if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
      wflag <- FALSE
      
      errs_mssg <- glue::glue("ERROR_MESSAGE")
      if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
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
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if ( ret_data ) {
    ret_dat <- NULL
    ret_dat$bsp_reg_sel_tib <- bsp_reg_sel_tib
    ret_dat$sel_all_tib <- sel_all_tib
    ret_dat$top_all_tib <- top_all_tib
    ret_dat$sum_all_tib <- sum_all_tib
    
    return( ret_dat )
  }
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               Evonik IO:: 
#                       Target Positions & Regions
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Evonik:: Load Chromosome Lengths
#
# fai_tab
# A tibble: 302,145 × 7
# source  chr   Chrom_Length Ord_Rank Len_Rank Len_Count Chr_Count
# <chr>   <chr>        <int>    <int>    <int>     <int>     <int>
# 1 galGal5 chr1     196202544        1        1         1     23475
# 2 galGal5 chr2     149560735     2852        2         1     23475
# 3 galGal5 chr3     111302122     4636        3         1     23475
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
load_fais_list = function( names_vec,
                           files_vec,
                           
                           vb=0, vt=6, tc=1, tt=NULL,
                           fun_tag='load_fais_list')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  cat(glue::glue("{mssg} NOT TESTED YET!!!{RET2}"))
  return(NULL)
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  dat_cols <- cols(
    chr = col_character(),
    len = col_integer(),
    off = col_double(),
    base = col_integer(),
    width = col_integer()
  )
  
  files <- list()
  for ( idx in c(1:length(names_vec) ) ) {
    name <- names_vec[idx]
    file <- files_vec[idx]
    files[[name]] <- file
  }
  
  ret_tib <- files %>% lapply( function(x) { 
    readr::read_tsv( x, 
                     col_names = names(dat_cols$cols), 
                     col_types = dat_cols ) %>%
      dplyr::mutate( Ord_Rank = dplyr::row_number() ) %>%
      dplyr::arrange( -len ) %>%
      dplyr::mutate( Len_Rank = dplyr::row_number() )
  } ) %>% 
    dplyr::bind_rows( .id = "source" ) %>% 
    dplyr::add_count( source,len, name="Len_Count" ) %>%
    dplyr::add_count( source, name="Chr_Count" ) %>% 
    dplyr::distinct( source,chr,len,Ord_Rank,Len_Rank,Len_Count,Chr_Count ) %>% 
    dplyr::rename( Chrom_Length = len )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

evonik_load_fais = function( top_path,
                             
                             vb=0, vt=6, tc=1, tt=NULL,
                             fun_tag='evonik_load_fais')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  dat_cols <- cols(
    chr = col_character(),
    len = col_integer(),
    off = col_double(),
    base = col_integer(),
    width = col_integer()
  )
  
  dat_files <- list()
  dat_files[["galGal5"]] <- ref_seqs$galGal5 %>% stringr::str_replace(".gz$", ".fai")
  dat_files[["criGriChoV1"]] <- ref_seqs$criGriChoV1 %>% stringr::str_replace(".gz$", ".fai")
  dat_files[["pvirPacbioDovetail"]] <- ref_seqs$pvirPacbioDovetail %>% stringr::str_replace(".gz$", ".fai")
  
  ret_tib <- dat_files %>% lapply( function(x) { 
    readr::read_tsv( x, 
                     col_names = names(dat_cols$cols), 
                     col_types = dat_cols ) %>%
      dplyr::mutate( Ord_Rank = dplyr::row_number() ) %>%
      dplyr::arrange( -len ) %>%
      dplyr::mutate( Len_Rank = dplyr::row_number() )
  } ) %>% 
    dplyr::bind_rows( .id = "source" ) %>% 
    dplyr::add_count( source,len, name="Len_Count" ) %>%
    dplyr::add_count( source, name="Chr_Count" ) %>% 
    dplyr::distinct( source,chr,len,Ord_Rank,Len_Rank,Len_Count,Chr_Count ) %>% 
    dplyr::rename( Chrom_Length = len )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Evonik:: Load Target CpG Sites
#
# pos_tab
# A tibble: 216,505 × 7
# source      source_type chr                beg     end srd   Group         
# <chr>       <chr>       <chr>            <int>   <int> <chr> <chr>         
# 1 criGriChoV1 Site        NW_003613580v1  537255  537256 -     choSMD15_NA_1 
# 2 criGriChoV1 Site        NW_003613580v1  741767  741768 +     choSMD15_NA_2 
# 3 criGriChoV1 Site        NW_003613580v1  785510  785511 +     choSMD15_NA_3 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
evonik_load_positions = function( top_path,
                                  
                                  vb=0, vt=6, tc=1, tt=NULL,
                                  fun_tag='evonik_load_positions')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  dat_cols <- cols(
    chr  = col_character(),
    beg  = col_integer(),
    end  = col_integer(),
    srd  = col_character(),
    grp1 = col_character(),
    grp2 = col_character(),
    grp3 = col_character(),
    grp4 = col_character()
  )
  
  dat_files <- list()
  dat_files[["galGal5"]] <- file.path( top_path, "Projects/Envonik/positions/gal/Epitrace_positions_with_category.formatted.group-sorted.txt" )
  dat_files[["criGriChoV1"]] <- file.path( top_path, "Projects/Envonik/positions/cho/SigMethDiff15cho_category.formatted.txt" )
  dat_files[["pvirPacbioDovetail"]] <- file.path( top_path, "Projects/Envonik/positions/crw/ProbesList_with_category.formatted.txt" )
  
  ret_tib <- dat_files %>% lapply( function(x) { 
    readr::read_tsv( x, 
                     col_names = names(dat_cols$cols), 
                     col_types = dat_cols )
  } ) %>% 
    dplyr::bind_rows( .id = "source") %>% 
    dplyr::arrange( source,chr,beg,end ) %>% 
    tidyr::pivot_longer( cols = c(grp1,grp2,grp3,grp4), 
                         names_to = "Group_Idx", 
                         values_to = "Category", 
                         values_drop_na = TRUE ) %>%
    dplyr::distinct() %>%
    dplyr::group_by( Category ) %>%
    dplyr::mutate( Category_Rep = dplyr::row_number() ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Category =  dplyr::case_when(
        source == "pvirPacbioDovetail" ~ paste(Category,Category_Rep, sep="_"),
        TRUE ~ Category
      ),
      Category_Pre = Category %>%  stringr::str_remove("_.*$") %>%
        stringr::str_remove("G$"),
      # Category_Pre = paste0(Category_Pre,"G"),
      Category_Idx = Category %>% stringr::str_remove("^.*_") %>% as.integer(),
      Category_Mid = Category %>% stringr::str_remove("^[^_]+_") %>%
        stringr::str_remove("[0-9]+$") %>%
        stringr::str_remove("_$") %>%
        stringr::str_replace_all("_","-"),
      Category_Mid = dplyr::case_when(
        stringr::str_length(Category_Mid) == 0 ~ "NA",
        TRUE ~ Category_Mid
      ),
      Group = paste( Category_Pre,Category_Mid,Category_Idx, sep="_")
    ) %>%
    dplyr::distinct( source,chr,beg,end,srd,Group ) %>%
    dplyr::arrange( source,chr,beg,end,srd,Group ) %>%
    dplyr::mutate( source_type = "Site" ) %>%
    dplyr::select( source,source_type,chr,beg,end,srd,Group )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Evonik:: Load Target Regions
#
# reg_tab
# A tibble: 52,502 × 7
# source      source_type chr     beg   end srd   Group     
# <chr>       <chr>       <chr> <int> <int> <chr> <chr>     
# 1 criGriChoV1 Region      chrM    732  2732 +     choG_NA_1 
# 2 criGriChoV1 Region      chrM   1897  3897 +     choG_NA_2 
# 3 criGriChoV1 Region      chrM   3309  5309 +     choG_NA_3 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
evonik_load_regions = function( top_path,
                                
                                vb=0, vt=6, tc=1, tt=NULL,
                                fun_tag='evonik_load_regions')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  dat_cols <- cols(
    chr  = col_character(),
    beg  = col_integer(),
    end  = col_integer(),
    srd  = col_character(),
    grpA = col_character(),
    grpB = col_character(),
    grpC = col_character()
  )
  
  dat_files <- list()
  dat_files[["galGal5_1"]] <- file.path( top_path, "Projects/Envonik/regions/gal/Genelist_all_only_promoter_with_category.txt" )
  dat_files[["galGal5_2"]] <- file.path( top_path, "Projects/Envonik/regions/gal/LMR_regions_with_category.formatted.txt" )
  dat_files[["criGriChoV1_1"]] <- file.path( top_path, "Projects/Envonik/regions/cho/genelist_all_only_promoter_2k_filtered_with_category_corrected.txt" )
  dat_files[["pvirPacbioDovetail_1"]] <- file.path( top_path, "Projects/Envonik/regions/crw/RegionsOfInterest_with_category.formatted.txt" )
  
  ret_tib <- dat_files %>% lapply( function(x) { 
    readr::read_tsv( x, 
                     col_names = names(dat_cols$cols), 
                     col_types = dat_cols )
  } ) %>% 
    dplyr::bind_rows( .id = "source" ) %>%
    dplyr::arrange( source,chr,beg,end,srd ) %>%
    tidyr::separate( source, into=c("source", "source_rep"), sep="_", convert = TRUE )  %>% 
    dplyr::add_count( source, name="Reg_Count" ) %>% 
    dplyr::mutate(
      Category = dplyr::case_when(
        is.na(grpC) ~ grpA,
        !is.na(grpC) ~ grpC,
        TRUE ~ NA_character_
      ),
      Category = stringr::str_remove( Category, "^[Cc]ategory_")
    ) %>% 
    dplyr::group_by( Category ) %>%
    dplyr::mutate( Category_Rep = dplyr::row_number() ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Category =  dplyr::case_when(
        # source == "pvirPacbioDovetail" ~ paste("crwG",Category_Rep, sep="_"),
        source == "pvirPacbioDovetail" ~ paste(Category,Category_Rep, sep="_"),
        TRUE ~ Category
      ),
      Category_Pre = Category %>%  stringr::str_remove("_.*$") %>%
        stringr::str_remove("G$"),
      Category_Pre = paste0(Category_Pre,"G"),
      Category_Idx = Category %>% stringr::str_remove("^.*_") %>% as.integer(),
      Category_Mid = Category %>% stringr::str_remove("^[^_]+_") %>%
        stringr::str_remove("[0-9]+$") %>%
        stringr::str_remove("_$") %>%
        stringr::str_replace_all("_","-"),
      Category_Mid = dplyr::case_when(
        stringr::str_length(Category_Mid) == 0 ~ "NA",
        TRUE ~ Category_Mid
      ),
      Group = paste( Category_Pre,Category_Mid,Category_Idx, sep="_")
    ) %>%
    dplyr::distinct( source,chr,beg,end,srd,Group ) %>%
    dplyr::arrange( source,chr,beg,end,srd,Group ) %>%
    dplyr::mutate( source_type = "Region" ) %>%
    dplyr::select( source,source_type,chr,beg,end,srd,Group )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Evonik:: Join all Region Data
#
# A tibble: 269,007 × 7
# source      source_type chr     beg   end srd   Group     
# <chr>       <chr>       <chr> <int> <int> <chr> <chr>     
# 1 criGriChoV1 Region      chrM    732  2732 +     choG_NA_1 
# 2 criGriChoV1 Region      chrM   1897  3897 +     choG_NA_2 
# 3 criGriChoV1 Region      chrM   3309  5309 +     choG_NA_3 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
evonik_load_targets = function( top_path,
                                
                                out_dir,
                                run_tag,
                                
                                reload     = 0,
                                reload_min = 2,
                                reload_pre = NULL,
                                
                                ret_data   = FALSE,
                                parallel   = FALSE,
                                
                                vb=0, vt=3, tc=1, tt=NULL,
                                fun_tag='evonik_load_targets')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt
  p1  <- vb > vt + 1
  p1  <- vb > vt + 2
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}     top_path = '{top_path}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
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
  
  errs_mssg <- glue::glue("Top Path='{top_path}' does not exist")
  if ( !dir.exists( top_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    pos_tib <- NULL
    pos_tib <- evonik_load_positions( top_path = top_path, 
                                      vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    reg_tib <- NULL
    reg_tib <- evonik_load_regions( top_path = top_path, 
                                    vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    ret_tib <- dplyr::bind_rows(
      dplyr::select( pos_tib,source,source_type,chr,beg,end,srd,Group ),
      dplyr::select( reg_tib,source,source_type,chr,beg,end,srd,Group )
    ) %>%
      dplyr::distinct( source,source_type,chr,beg,end,srd,Group ) %>%
      dplyr::arrange( source,chr,beg,end,srd,source_type,Group ) %>%
      dplyr::select( source,source_type,chr,beg,end,srd,Group )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Delta-Beta Functions:: Tibble
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

delta_beta_across_wrap = function( pair_list, call_listX, call_listY, pval_min, 
                                   prefix = "dB", del = "_",
                                   
                                   vb=0, vt=3, tc=1, tt=NULL,
                                   fun_tag='delta_beta_across_wrap') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  ret_cnt <- 0
  ret_tib <- NULL
  dB_tab  <- NULL
  
  true_n <- pair_list %>% length()
  
  # pair_list %>% head( n=2 ) %>% print()
  # call_listA %>% head( n=2 ) %>% print()
  # call_listB %>% head( n=2 ) %>% print()
  
  ftime <- system.time({
    dB_tab <- foreach ( idx_pair = pair_list, .combine = rbind ) %dopar% {
      # print(idx_pair)
      
      delta_beta_across_vec2( datX = call_listX,
                              datY = call_listY,
                              x = idx_pair[[2]],
                              y = idx_pair[[3]],
                              name = "dB",
                              pids="Probe_ID", beta="beta", pval="pval",
                              min = pval_min,
                              prefix = "dB", del = "_" )
    }
    
    ret_tib <- dB_tab %>% dplyr::group_by(Probe_ID) %>%
      dplyr::summarise( pass_comb_cnt = n(),
                        full_comb_cnt = true_n,
                        db_miss_cnt = full_comb_cnt - pass_comb_cnt,
                        dB_fail_cnt = sum(dB >  0.2, na.rm = TRUE),
                        dB_pass_cnt = sum(dB <= 0.2, na.rm = TRUE),
                        dB_miss_per = round( 100*db_miss_cnt / full_comb_cnt ),
                        dB_comb_per = round( 100*dB_pass_cnt / full_comb_cnt ),
                        dB_pass_per = round( 100*dB_pass_cnt / pass_comb_cnt ),
                        
                        
                        dB_min = min( dB, na.rm = TRUE),
                        dB_max = max( dB, na.rm = TRUE),
                        dB_avg = mean( dB, na.rm = TRUE),
                        dB_med = median( dB, na.rm = TRUE),
                        dB_sds = sd( dB, na.rm = TRUE),
                        
                        .groups = "drop" )
  })
  
  ret_tib %>% dplyr::group_by(pass_comb_cnt, dB_fail_cnt ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" ) %>%
    print(n=3)
  
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

delta_beta_inside_wrap = function( pair_list, call_list, pval_min, 
                                   prefix = "dB", del = "_",
                                   
                                   vb=0, vt=3, tc=1, tt=NULL,
                                   fun_tag='delta_beta_inside_wrap') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  ret_cnt <- 0
  ret_tib <- NULL
  dB_tab  <- NULL
  
  true_n <- pair_list %>% length()
  
  # pair_list %>% head( n=2 ) %>% print()
  # call_list %>% head( n=2 ) %>% print()
  
  ftime <- system.time({
    dB_tab <- foreach ( idx_pair = pair_list, .combine = rbind ) %dopar% {
      # print(idx_pair)
      
      delta_beta_inside_vec2( dat = call_list,
                              x = idx_pair[[2]],
                              y = idx_pair[[3]],
                              name = "dB",
                              pids="Probe_ID", beta="beta", pval="pval",
                              min = pval_min,
                              prefix = "dB", del = "_" )
    }
    
    ret_tib <- dB_tab %>% dplyr::group_by(Probe_ID) %>%
      dplyr::summarise( pass_comb_cnt = n(),
                        full_comb_cnt = true_n,
                        db_miss_cnt = full_comb_cnt - pass_comb_cnt,
                        dB_fail_cnt = sum(dB >  0.2, na.rm = TRUE),
                        dB_pass_cnt = sum(dB <= 0.2, na.rm = TRUE),
                        dB_miss_per = round( 100*db_miss_cnt / full_comb_cnt ),
                        dB_comb_per = round( 100*dB_pass_cnt / full_comb_cnt ),
                        dB_pass_per = round( 100*dB_pass_cnt / pass_comb_cnt ),
                        
                        dB_min = min( dB, na.rm = TRUE),
                        dB_max = max( dB, na.rm = TRUE),
                        dB_avg = mean( dB, na.rm = TRUE),
                        dB_med = median( dB, na.rm = TRUE),
                        dB_sds = sd( dB, na.rm = TRUE),
                        
                        .groups = "drop" )
  })
  
  ret_tib %>% dplyr::group_by(pass_comb_cnt, dB_fail_cnt ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" ) %>%
    print(n=3)
  
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

beta_inside_wrap = function( call_list, pval_min, 
                             
                             vb=0, vt=3, tc=1, tt=NULL,
                             fun_tag='beta_inside_wrap') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  ret_cnt  <- 0
  ret_tib  <- NULL
  beta_tab <- NULL
  
  # pair_list %>% head( n=2 ) %>% print()
  # call_list %>% head( n=2 ) %>% print()
  
  ftime <- system.time({
    
    beta_tab <- call_list %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate( beta = dplyr::case_when( 
        pval > pval_min ~ NA_real_, 
        TRUE ~ beta ) )
    
    ret_tib <- beta_tab %>% 
      dplyr::group_by(Probe_ID) %>% 
      dplyr::summarise( sample_cnt = n(), 
                        pval_fail_cnt = sum( pval >  pval_min, na.rm = TRUE), 
                        pval_pass_cnt = sum( pval <= pval_min, na.rm = TRUE),
                        pval_pass_per = round( 100*pval_pass_cnt / sample_cnt ),
                        
                        beta_min = min( beta, na.rm = TRUE),
                        beta_max = max( beta, na.rm = TRUE),
                        beta_avg = mean( beta, na.rm = TRUE),
                        beta_med = median( beta, na.rm = TRUE),
                        beta_sds = sd( beta, na.rm = TRUE),
                        
                        .groups = "drop" )
    
  })
  
  ret_tib %>% dplyr::group_by(sample_cnt, pval_fail_cnt ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" ) %>%
    print(n=3)
  
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

delta_beta_inside_vec2 = function( dat, x, y, 
                                   name=NULL,
                                   pids="Probe_ID", beta="beta", pval="pval",
                                   min = 0.05,
                                   prefix = "dB", del = "_",
                                   
                                   vb=0, vt=3, tc=1, tt=NULL,
                                   fun_tag='delta_beta_inside_vec2') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  col_str <- paste( prefix, x, y, sep=del )
  if ( !is.null(name) ) col_str <- name
  
  abs_str <- paste(col_str,"abs", sep=del)
  col_sym <- rlang::sym(col_str)
  abs_sym <- rlang::sym(abs_str)
  
  pids_sym <- rlang::sym(pids)
  beta_sym <- rlang::sym(beta)
  pval_sym <- rlang::sym(pval)
  
  ftime <- base::system.time({
    
    datX <- dat[[x]] %>% dplyr::filter( !!pval_sym <= min ) %>%
      dplyr::select( !!pids_sym, !!beta_sym )
    
    datY <- dat[[y]] %>% dplyr::filter( !!pval_sym <= min ) %>%
      dplyr::select( !!pids_sym, !!beta_sym )
    
    datB <- dplyr::inner_join( datX,datY, by=pids, suffix=c("_A","_B"))
    
    datB <- datB %>% dplyr::mutate(
      # !!col_sym := datB[[2]] - datB[[3]],
      # !!abs_sym := base::abs( !!col_sym )
      !!col_sym := base::abs( datB[[2]] - datB[[3]] )
    ) %>% 
      # dplyr::select( 1, 4,5 ) %>% dplyr::distinct()
      dplyr::select( 1, 4 ) %>% dplyr::distinct()
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  datB
}

delta_beta_across_vec2 = function( datX, datY, x, y, 
                                   name=NULL,
                                   pids="Probe_ID", beta="beta", pval="pval",
                                   min = 0.05,
                                   prefix = "dB", del = "_",
                                   
                                   vb=0, vt=3, tc=1, tt=NULL,
                                   fun_tag='delta_beta_across_vec2') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  col_str <- paste( prefix, x, y, sep=del )
  if ( !is.null(name) ) col_str <- name
  
  abs_str <- paste(col_str,"abs", sep=del)
  col_sym <- rlang::sym(col_str)
  abs_sym <- rlang::sym(abs_str)
  
  pids_sym <- rlang::sym(pids)
  beta_sym <- rlang::sym(beta)
  pval_sym <- rlang::sym(pval)
  
  ftime <- base::system.time({
    
    datX <- datX[[x]] %>% dplyr::filter( !!pval_sym <= min ) %>%
      dplyr::select( !!pids_sym, !!beta_sym )
    
    datY <- datY[[y]] %>% dplyr::filter( !!pval_sym <= min ) %>%
      dplyr::select( !!pids_sym, !!beta_sym )
    
    datB <- dplyr::inner_join( datX,datY, by=pids, suffix=c("_A","_B"))
    
    datB <- datB %>% dplyr::mutate(
      # !!col_sym := datB[[2]] - datB[[3]],
      # !!abs_sym := base::abs( !!col_sym )
      !!col_sym := base::abs( datB[[2]] - datB[[3]] )
    ) %>% 
      # dplyr::select( 1, 4,5 ) %>% dplyr::distinct()
      dplyr::select( 1, 4 ) %>% dplyr::distinct()
    
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  datB
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Delta-Beta Functions:: Matrix
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

call_list_to_matrix = function( calls, out_path = NULL, out_name = "alpha-testing" ) {
  
  ids_vec <- calls %>%
    dplyr::bind_rows() %>% 
    dplyr::distinct( Probe_ID ) %>%
    dplyr::pull( Probe_ID )
  
  temp_beta_tib <- tibble::tibble( Probe_ID = ids_vec )
  temp_pval_tib <- tibble::tibble( Probe_ID = ids_vec )
  
  mast_beta_tib <- tibble::tibble( Probe_ID = ids_vec )
  mast_pval_tib <- tibble::tibble( Probe_ID = ids_vec )
  
  for ( id in names(calls) ) {
    # check_unique_pid( calls[[id]], paste(id,"calls",sep="_" ) )
    
    cur_beta_tib <- NULL
    cur_pval_tib <- NULL
    
    cur_beta_tib <- calls[[id]] %>% dplyr::select( Probe_ID, beta ) %>% dplyr::distinct()
    cur_pval_tib <- calls[[id]] %>% dplyr::select( Probe_ID, pval ) %>% dplyr::distinct()
    success <- check_unique_pid( cur_beta_tib, paste(id,"Merge1",sep="_" ) )
    if ( !success ) return(NULL)
    
    #
    #  1. Check uniquenes
    #  2. Merge to template
    #  3. Merge template+sample to master
    #
    cur_beta_tib <- temp_beta_tib %>% dplyr::left_join( cur_beta_tib, by="Probe_ID" )
    cur_pval_tib <- temp_pval_tib %>% dplyr::left_join( cur_pval_tib, by="Probe_ID" )
    success <- check_unique_pid( cur_beta_tib, paste(id,"Merge2",sep="_" ) )
    if ( !success ) return(NULL)
    
    #
    # Now join to master::
    #
    mast_beta_tib <- mast_beta_tib %>% dplyr::left_join(cur_beta_tib, by="Probe_ID" )
    mast_pval_tib <- mast_pval_tib %>% dplyr::left_join(cur_pval_tib, by="Probe_ID" )
    success <- check_unique_pid( mast_beta_tib, paste(id,"Master",sep="_" ) )
    if ( !success ) return(NULL)
  }
  success <- check_unique_pid( mast_beta_tib, "Master-Done" )
  if ( !success ) return(NULL)
  
  mast_beta_tib <- mast_beta_tib %>% purrr::set_names( c( "Probe_ID", names(calls) ) )
  mast_pval_tib <- mast_pval_tib %>% purrr::set_names( c( "Probe_ID", names(calls) ) )
  
  # test_beta_rds <- file.path( out_path, paste(out_name,"beta-matrix.rds", sep="." ) )
  beta_mat <- mast_beta_tib %>% 
    tibble::column_to_rownames( var = "Probe_ID" ) %>% data.matrix()
  beta_mat %>% head() %>% print()
  # readr::write_rds( beta_mat, test_beta_rds )
  
  # test_pval_rds <- file.path( out_path, paste(out_name,"pval-matrix.rds", sep="." ) )
  pval_mat <- mast_pval_tib %>% 
    tibble::column_to_rownames( var = "Probe_ID" ) %>% data.matrix()
  pval_mat %>% head() %>% print()
  # readr::write_rds( pval_mat, test_pval_rds )
  
  ret_dat <- NULL
  ret_dat$beta <- beta_mat
  ret_dat$pval <- pval_mat
  
  ret_dat
}

add_delta_beta_mat = function( tib, dat, x, y, prefix = "dB", del = "_" ) {
  
  # col_str <- paste( prefix, x, y, sep=del )
  col_str <- prefix
  col_sym <- rlang::sym(col_str)
  
  tib <- tib %>% dplyr::mutate(
    # !!col_sym := dat[[x]] - dat[[y]]
    # !!col_sym := dat[ , x] - dat[ , y]
    !!col_sym := base::abs(dat[ , x] - dat[ , y])
  )
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Simple Off Hand Functions::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

check_unique_pid = function( tib, id = "Sample" ) {
  
  success <- "Success!"
  tot_cnt <- tib %>% base::nrow()
  unq_cnt <- tib %>% dplyr::distinct( Probe_ID ) %>% base::nrow()
  if ( tot_cnt != unq_cnt ) success <- "Failure!"
  cat(glue::glue("[{id}]: tot={tot_cnt}, unq={unq_cnt}: {success}{RET}"))
  if ( tot_cnt != unq_cnt ) return(FALSE)
  
  return(TRUE)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Mutate UCSC dbSNP Functions::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

wash_snvs = function( imp_tib,
                      snv_tib,
                      
                      out_dir,
                      run_tag,
                      
                      reload     = 0,
                      reload_min = 2,
                      reload_pre = NULL,
                      
                      ret_data   = FALSE,
                      parallel   = FALSE,
                      
                      vb=0, vt=3, tc=1, tt=NULL,
                      fun_tag='wash_snvs') {
  
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
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
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
  
  # errs_mssg <- glue::glue("File file='{file}' does not exist")
  # if ( !file.exists( file) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    off_tib <- imp_tib %>% 
      dplyr::distinct( Ilmn_ID, Loci_Num, Rsn, Var_Type, 
                       Probe_Type, Strand_FR, Strand_CO, 
                       Forward_Sequence, Offset, Distance,
                       Chromosome, Coordinate ) %>% 
      dplyr::mutate( Ilmn_ID = paste0(Ilmn_ID, Strand_FR) ) %>%
      dplyr::inner_join( snv_tib, by=c("Rsn", "Chromosome", "Var_Type") ) %>%
      dplyr::arrange( Ilmn_ID, End, Distance, Rsn ) %>% 
      dplyr::add_count( Ilmn_ID, name="Var_Cnt" ) %>%
      dplyr::group_by( Ilmn_ID ) %>%
      dplyr::mutate( Var_Idx = dplyr::row_number() ) %>%
      dplyr::ungroup() %>%
      # dplyr::filter( Var_Cnt == 2 ) %>% head() %>%
      dplyr::distinct( Ilmn_ID, Var_Cnt, Forward_Sequence, Distance, Offset, IUPAC, Ref ) %>%
      # dplyr::distinct( Loci_Num, Probe_Type, Strand_FR, Strand_CO, Var_Cnt, Forward_Sequence, Distance, Offset, IUPAC, Ref ) %>%
      dplyr::mutate( Forward_Sequence = Forward_Sequence %>% shear_brac() )
    
    seq_tib <- off_tib %>% 
      dplyr::distinct( Ilmn_ID, Var_Cnt, Forward_Sequence )
    # dplyr::distinct( Loci_Num, Probe_Type, Strand_FR, Strand_CO, Var_Cnt, Forward_Sequence )
    
    #
    # TBD:: Investigate offsets with lower MAF to see that they are not all -3 to 0 offset...
    #
    fwd_tib <- replace_nucs_at_cpp( 
      cgn_vec_r = seq_tib$Ilmn_ID, 
      cnt_vec_r = seq_tib$Var_Cnt, 
      seq_vec_r = seq_tib$Forward_Sequence, 
      pos_vec_r = off_tib$Offset, 
      alt_vec_r = off_tib$IUPAC, 
      ref_vec_r = off_tib$Ref, 
      vb=0, vt=10 ) %>% 
      tibble::as_tibble() %>%
      dplyr::mutate( Probe_Type   = stringr::str_sub(Probe_ID, 1,2),
                     Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"),
                     Strand_TB = Probe_ID %>% stringr::str_remove("^.*_") %>% stringr::str_sub( 1,1 ),
                     Strand_CO = Probe_ID %>% stringr::str_remove("^.*_") %>% stringr::str_sub( 2,2 ),
                     Strand_FR = Probe_ID %>% stringr::str_remove("^.*_") %>% stringr::str_sub( 3,3 ),
                     SNP_Sequence = SNP_Sequence, 
                     RAW_Sequence = RAW_Sequence )
    
    ref_tib <- NULL
    ref_tib <- improbe_seqs_cpp( 
      # ids_vec_r = fwd_tib$Loci_ID,
      ids_vec_r = fwd_tib$Probe_ID,
      fwd_vec_r = fwd_tib$RAW_Sequence, 
      din_vec_r = fwd_tib$Probe_Type,
      bsc_str_r = "numd", 
      frs_vec_r_ = fwd_tib$Strand_FR, 
      cos_vec_r_ = fwd_tib$Strand_CO,
      prb_len = 50, vb = 1 ) %>%
      dplyr::bind_rows() %>% dplyr::filter( Probe_ID != "" ) %>%
      dplyr::mutate( PRB_U=stringr::str_to_upper(PRB_U), 
                     PRB_M=stringr::str_to_upper(PRB_M), 
                     PRB_D=stringr::str_to_upper(PRB_D),
                     PRB_N=stringr::str_to_upper(PRB_N) )
    
    snp_tib <- NULL
    snp_tib <- improbe_seqs_cpp( 
      # ids_vec_r = fwd_tib$Loci_ID,
      ids_vec_r = fwd_tib$Probe_ID,
      fwd_vec_r = fwd_tib$SNP_Sequence, 
      din_vec_r = fwd_tib$Probe_Type,
      bsc_str_r = "numd", 
      frs_vec_r_ = fwd_tib$Strand_FR, 
      cos_vec_r_ = fwd_tib$Strand_CO,
      prb_len = 50, vb = 1 ) %>%
      dplyr::bind_rows() %>% dplyr::filter( Probe_ID != "" ) %>%
      dplyr::mutate( PRB_U=stringr::str_to_upper(PRB_U), 
                     PRB_M=stringr::str_to_upper(PRB_M), 
                     PRB_D=stringr::str_to_upper(PRB_D),
                     PRB_N=stringr::str_to_upper(PRB_N) )
    
    ret_tib <- dplyr::left_join(
      ref_tib, snp_tib,
      by=c("Probe_ID", "Strand_SR", "Strand_CO"),
      suffix=c("_ref", "_snp")
    ) %>%
      dplyr::mutate( 
        Infinum1_Affected = dplyr::case_when(
          is.na( PRB_U_snp ) | is.na( PRB_M_snp ) ~ -1,
          PRB_U_snp == PRB_U_ref | PRB_M_snp == PRB_M_ref ~ 0,
          PRB_U_snp != PRB_U_ref | PRB_M_snp != PRB_M_ref ~ 1,
          TRUE ~ NA_real_ ) %>% as.integer(),
        Infinum2_Affected = dplyr::case_when(
          is.na( PRB_D_snp ) ~ -1,
          PRB_D_snp == PRB_D_ref ~ 0,
          PRB_D_snp != PRB_D_ref ~ 1,
          TRUE ~ NA_real_ ) %>% as.integer()
      )
    
    # warn_mssg <- glue::glue("WARN_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # wflag <- FALSE
    # 
    # errs_mssg <- glue::glue("ERROR_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    # if ( eflag ) return( NULL )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
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
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Format UCSC dbSNP Functions::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

dbsnp_table_to_IUPAC = function( tib,
                                 maj = NULL,
                                 
                                 exp_max = 12,
                                 
                                 out_dir,
                                 run_tag,
                                 
                                 reload     = 0,
                                 reload_min = 2,
                                 reload_pre = NULL,
                                 
                                 ret_data   = FALSE,
                                 parallel   = FALSE,
                                 
                                 vb=0, vt=3, tc=1, tt=NULL,
                                 fun_tag='dbsnp_table_to_IUPAC') {
  
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
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      exp_max = '{exp_max}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
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
  
  # errs_mssg <- glue::glue("File file='{file}' does not exist")
  # if ( !file.exists( file) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    ret_tib <- tib %>%
      dplyr::select( Chromosome, Beg, End, Rsn, Var_Type, Ref, Major, Minor, Experiment_Idx ) %>% 
      dplyr::filter( Var_Type == "snv" ) %>% 
      dplyr::filter( Experiment_Idx <= exp_max ) %>%
      # dplyr::filter( Rsn == "rs9256983" | Rsn == "rs9256982" ) %>%
      dplyr::distinct() %>% 
      tidyr::pivot_longer( cols = c(Major,Minor) , names_to = "Var" ) %>% 
      dplyr::select( -Var, -Experiment_Idx )
    
    if ( !is.null(maj) ) ret_tib <- ret_tib %>% 
      dplyr::bind_rows( maj_imp_snv_tib )
    
    ret_tib <- ret_tib %>% 
      dplyr::arrange( Chromosome, Beg, End, Rsn ) %>%
      dplyr::distinct() %>% 
      tidyr::pivot_wider( names_from = value ) %>% 
      tidyr::unite( nucs,  A,C,T,G, na.rm = TRUE, sep = "" ) %>% 
      dplyr::mutate( IUPAC = mapDIs(nucs) ) %>% 
      dplyr::filter( Ref != IUPAC ) %>%
      dplyr::distinct()
    
    # warn_mssg <- glue::glue("WARN_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # wflag <- FALSE
    # 
    # errs_mssg <- glue::glue("ERROR_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    # if ( eflag ) return( NULL )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
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
  
  ret_tib
}

ucsc_dbsnp_to_table = function( tib,
                                
                                out_dir,
                                run_tag,
                                
                                reload     = 0,
                                reload_min = 2,
                                reload_pre = NULL,
                                
                                ret_data   = FALSE,
                                parallel   = FALSE,
                                
                                vb=0, vt=3, tc=1, tt=NULL,
                                fun_tag='ucsc_dbsnp_to_table')
{
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
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
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
  
  # errs_mssg <- glue::glue("File file='{file}' does not exist")
  # if ( !file.exists( file) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    alt_col_vec <- paste0( "XALT", c(1:3),"_XXALT" )
    maf_col_vec <- paste0( "IDX", c(1:16),"_MAF" )
    maj_col_vec <- paste0( "IDX", c(1:16),"_Major" )
    min_col_vec <- paste0( "IDX", c(1:16),"_Minor" )
    
    com_tib <- tib %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        alts = alts %>% stringr::str_remove(",$"),
        mafs = mafs %>% stringr::str_remove(",$"),
        maxs = maxs %>% stringr::str_remove(",$"),
        mins = mins %>% stringr::str_remove(",$"),
        var_note = var_note %>% stringr::str_remove(",$")
      )
    
    com_cols <- names(com_tib)
    
    com_sum <- com_tib %>% 
      dplyr::group_by( var_type ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    if ( vb >= vt+ 1 ) com_sum %>% print( n=base::nrow(com_sum) )
    
    ret_tib <- dplyr::bind_rows(
      
      com_tib %>%
        dplyr::distinct() %>% 
        dplyr::filter( var_type != "snv" ),
      com_tib %>%
        dplyr::distinct() %>% 
        dplyr::filter( var_type == "snv" ) %>%
        tidyr::separate( maxs, into=maj_col_vec, sep=',', convert = TRUE, remove=FALSE ) %>%
        tidyr::separate( mins, into=min_col_vec, sep=',', convert = TRUE, remove=FALSE ) %>%
        dplyr::mutate( dplyr::across( dplyr::everything(.) ), na_if(., "") ) %>%
        tidyr::unite( MAJ_Str, dplyr::all_of( maj_col_vec ), sep=",", na.rm = TRUE ) %>%
        tidyr::unite( MIN_Str, dplyr::all_of( min_col_vec ), sep=",", na.rm = TRUE ) %>%
        tidyr::separate( alts, into=alt_col_vec, sep=',', convert = TRUE, remove = FALSE ) %>%
        tidyr::pivot_longer( cols = dplyr::starts_with("XALT"), 
                             names_to = c("XALT", ".value"), 
                             names_sep = "_", 
                             values_drop_na = TRUE ) %>%
        dplyr::mutate(
          Match_Ref_Maj = stringr::str_detect(MAJ_Str, ref),
          Match_Ref_Min = stringr::str_detect(MIN_Str, ref),
          Match_Alt_Maj = stringr::str_detect(MAJ_Str, XXALT),
          Match_Alt_Min = stringr::str_detect(MIN_Str, XXALT),
          mafs = dplyr::case_when(
            !Match_Ref_Maj & !Match_Ref_Min ~ paste0(mafs,",-0.01"),
            TRUE ~ paste0(mafs,",") ),
          maxs = dplyr::case_when(
            !Match_Ref_Maj & !Match_Ref_Min ~ paste0(maxs,",",stringr::str_sub( maxs, 1,1 ) ),
            TRUE ~ paste0(maxs,",") ),
          mins = dplyr::case_when(
            !Match_Ref_Maj & !Match_Ref_Min ~ paste0(mins,",",ref ),
            TRUE ~ paste0(mins,",") ),
          
          mafs = dplyr::case_when(
            !Match_Alt_Maj & !Match_Alt_Min ~ paste0(mafs,",-0.02"),
            TRUE ~ paste0(mafs,",") ),
          maxs = dplyr::case_when(
            !Match_Alt_Maj & !Match_Alt_Min ~ paste0(maxs,",",stringr::str_sub( maxs, 1,1 ) ),
            TRUE ~ paste0(maxs,",") ),
          mins = dplyr::case_when(
            !Match_Alt_Maj & !Match_Alt_Min ~ paste0(mins,",",XXALT ),
            TRUE ~ paste0(mins,",") ),
          mafs_len = stringr::str_length( mafs )
        ) %>% dplyr::arrange( -mafs_len ) %>% 
        dplyr::distinct(rsn, .keep_all = TRUE ) %>%
        dplyr::select( dplyr::all_of( com_cols ) )
      # dplyr::filter( !Match_Ref_Maj & !Match_Ref_Min ) %>% as.data.frame()
      # dplyr::filter( !Match_Alt_Maj & !Match_Alt_Min ) %>% head() %>% as.data.frame()
    ) %>% 
      dplyr::distinct() %>% 
      tidyr::separate( mafs, into=maf_col_vec, sep=',', convert = TRUE ) %>%
      tidyr::separate( maxs, into=maj_col_vec, sep=',', convert = TRUE ) %>%
      tidyr::separate( mins, into=min_col_vec, sep=',', convert = TRUE ) %>%
      dplyr::mutate( dplyr::across( dplyr::everything(.) ), na_if(., -Inf) ) %>%
      dplyr::mutate( dplyr::across( dplyr::everything(.) ), na_if(., "") ) %>%
      dplyr::mutate( dplyr::across( dplyr::all_of( maf_col_vec ), ~ .x * 100 ) ) %>%
      tidyr::unite( mafs, dplyr::all_of( maf_col_vec ), sep=",", na.rm = TRUE, remove = FALSE ) %>%
      tidyr::unite( maxs, dplyr::all_of( maj_col_vec ), sep=",", na.rm = TRUE, remove = FALSE ) %>%
      tidyr::unite( mins, dplyr::all_of( min_col_vec ), sep=",", na.rm = TRUE, remove = FALSE ) %>%
      # dplyr::mutate( alts = alts %>% stringr::str_remove_all(",") ) %>%
      tidyr::pivot_longer( cols = dplyr::starts_with("IDX"), 
                           names_to = c("IDX", ".value"), 
                           names_sep = "_", 
                           values_drop_na = TRUE ) %>%
      dplyr::mutate(
        IDX = IDX %>% stringr::str_remove("^IDX") %>% as.integer(),
      )
    
    
    # warn_mssg <- glue::glue("WARN_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # wflag <- FALSE
    # 
    # errs_mssg <- glue::glue("ERROR_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    # if ( eflag ) return( NULL )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
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
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Tabix Intersection Function::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersect_bed_tabix = function( tib,
                                file,
                                cols = NULL,
                                
                                ref_vcf   = FALSE,
                                com_len   = 0,
                                rm_chr0   = TRUE,
                                strip_chr = FALSE,
                                snps_only = FALSE,
                                match_src = FALSE,
                                parse_caf = FALSE,
                                min_maf = 0.0,
                                
                                ids_key,
                                chr_key,
                                beg_key,
                                end_key,
                                
                                beg_buf = 0,
                                end_buf = 0,
                                beg_off = 1,
                                
                                out_dir,
                                run_tag,
                                
                                reload     = 0,
                                reload_min = 2,
                                reload_pre = NULL,
                                
                                ret_data   = FALSE,
                                parallel   = FALSE,
                                write_out  = FALSE,
                                
                                vb=0, vt=3, tc=1, tt=NULL,
                                fun_tag='intersect_bed_tabix')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  p0  <- vb > vt + 0
  p1  <- vb > vt + 1
  p2  <- vb > vt + 2
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL

  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  reg_bed <- file.path( out_dir, paste(out_tag, 'reg.bed.gz', sep='.') )
  tab_tsv <- file.path( out_dir, paste(out_tag, 'tab.bed.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}      ref_vcf = '{ref_vcf}'.{RET}"))
    cat(glue::glue("{mssg}      com_len = '{com_len}'.{RET}"))
    cat(glue::glue("{mssg}      rm_chr0 = '{rm_chr0}'.{RET}"))
    cat(glue::glue("{mssg}    strip_chr = '{strip_chr}'.{RET}"))
    cat(glue::glue("{mssg}    snps_only = '{snps_only}'.{RET}"))
    cat(glue::glue("{mssg}    match_src = '{match_src}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}      ids_key = '{ids_key}'.{RET}"))
    cat(glue::glue("{mssg}      chr_key = '{chr_key}'.{RET}"))
    cat(glue::glue("{mssg}      beg_key = '{beg_key}'.{RET}"))
    cat(glue::glue("{mssg}      end_key = '{end_key}'.{RET}"))
    cat(glue::glue("{mssg}      beg_buf = '{beg_buf}'.{RET}"))
    cat(glue::glue("{mssg}      end_buf = '{end_buf}'.{RET}"))
    cat(glue::glue("{mssg}      beg_off = '{beg_off}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid  = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload  = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min  = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}   reload_pre  = '{reload_pre}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data  = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel  = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     write_out = '{write_out}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_tag = '{out_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      reg_bed = '{reg_bed}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c( reg_bed, tab_tsv, out_csv, end_txt) )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("File file='{file}' does not exist")
  if ( !file.exists( file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ids_sym <- rlang::sym( ids_key )
  chr_sym <- rlang::sym( chr_key )
  beg_sym <- rlang::sym( beg_key )
  end_sym <- rlang::sym( end_key )
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # Sort tib by chr, beg, end::
    if ( rm_chr0 ) tib <- tib %>%
        dplyr::filter( (!!chr_sym) != "chr0" ) %>%
        dplyr::filter( (!!chr_sym) != "0" )

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Write BED FIle::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( strip_chr ) tib <- tib %>% 
        dplyr::mutate( !!chr_sym := !!chr_sym %>% stringr::str_remove("^chr") )

    reg_tib <- NULL
    reg_tib <- tib %>% 
      dplyr::select( dplyr::all_of( c(chr_key, beg_key, end_key, ids_key) ) ) %>%
      dplyr::mutate( REG_BEG = !!beg_sym - beg_buf - beg_off,
                     REG_END = !!end_sym + end_buf ) %>%
      dplyr::arrange(  !!chr_sym, REG_BEG, REG_END, !!ids_sym ) %>%
      dplyr::distinct( !!chr_sym, REG_BEG, REG_END, !!ids_sym ) %>%
      dplyr::select( !!chr_sym, REG_BEG, REG_END, !!ids_sym )
    
    readr::write_tsv( x = reg_tib, file = reg_bed, col_names = FALSE )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                Run Tabix::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cmd = paste0( "tabix -R ",reg_bed," ",file," | gzip -c - > ",tab_tsv )
    if ( ref_vcf ) 
      cmd = paste0( "tabix -h -R ",reg_bed," ",file," | gzip -c - > ",tab_tsv )
    if ( p2 ) cat(glue::glue("{mssg} Running Cmd='{cmd}'...{RET2}"))
    base::system( command = cmd, ignore.stdout = FALSE, ignore.stderr = TRUE )
    
    errs_mssg <- glue::glue("Intersection file='{tab_tsv}' does not exist")
    if ( !file.exists( tab_tsv) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Load Tabix Results::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    tab_tib <- NULL
    if ( !is.null(cols) ) {
      tab_tib <- readr::read_tsv( tab_tsv, 
                                  col_names = names(cols$cols), 
                                  col_types = cols )
    } else if ( ref_vcf) {
      tab_tib <- readr::read_tsv( tab_tsv,
                                  skip = com_len, show_col_types = FALSE  ) %>% 
        magrittr::set_names( c("Chromosome","Coordinate_SNP","SNP_ID",
                               "REF_SNP","ALT_SNP","QUAL","FILTER","INFO") )
      
      if ( snps_only ) tab_tib <- tab_tib %>%
          dplyr::filter( stringr::str_length(REF_SNP) == 1 ) %>%
          # This removes most large indels::
          dplyr::filter( ALT_SNP %>% stringr::str_length() <= 5 ) %>%
          dplyr::filter( REF_SNP != ALT_SNP ) %>%
          dplyr::filter( !SNP_ID %>% stringr::str_detect(",") ) %>% 
          dplyr::distinct( Chromosome,Coordinate_SNP,SNP_ID,REF_SNP,ALT_SNP, .keep_all = TRUE ) %>%
          dplyr::select( Chromosome, Coordinate_SNP, SNP_ID,REF_SNP,ALT_SNP, INFO ) %>%
          clean_tib() %>% 
          dplyr::mutate( ALT_LenT = ALT_SNP %>% stringr::str_length(), 
                         ALT_LenC = ALT_SNP %>% stringr::str_remove_all(",") %>% stringr::str_length(),
                         ALT_LenD = ALT_LenT - ALT_LenC ) %>% 
          dplyr::filter( !(ALT_LenT > 1 & ALT_LenD == 0) ) %>%
          dplyr::distinct( SNP_ID, .keep_all = TRUE ) %>%
          dplyr::select( Chromosome, Coordinate_SNP, SNP_ID, REF_SNP, ALT_SNP, INFO )
      
      if ( parse_caf ) {
        tab_tib <- tab_tib %>%
          dplyr::mutate( 
            INFO    = INFO %>% stringr::str_remove("^.*CAF=") %>% stringr::str_remove(";.*$"), 
            REF_AF  = INFO %>% stringr::str_remove(",.*$") %>% as.double(),
            ALT_AF1 = INFO %>% stringr::str_remove("^[^,]+,") %>% stringr::str_remove(",.*$") %>% as.double(),
            ALT_AF2 = INFO %>% stringr::str_remove("^[^,]+,") %>% stringr::str_remove(",.*$") %>% as.double(),
            ALT_AF3 = INFO %>% stringr::str_remove("^[^,]+,") %>% stringr::str_remove(",.*$") %>% as.double()
          ) %>%
          dplyr::mutate( 
            INFO    = INFO %>% stringr::str_remove("^.*CAF=") %>% stringr::str_remove(";.*$")
          ) %>%
          tidyr::separate( INFO, into = c("AF0","AF1","AF2","AF3"), sep = ",", 
                           convert = TRUE, fill = "right", remove = FALSE, extra = "drop" ) %>%
          tidyr::separate( ALT_SNP, into = c("ALT1","ALT2","ALT3"), sep = ",", 
                           convert = TRUE, fill = "right", remove = FALSE, extra = "drop" )
      
        tab_tib <- dplyr::bind_rows(
          tab_tib %>% dplyr::mutate( MAF = AF1, ALT = as.character(ALT1) ) %>% dplyr::filter( !is.na(ALT) ),
          tab_tib %>% dplyr::mutate( MAF = AF2, ALT = as.character(ALT2) ) %>% dplyr::filter( !is.na(ALT) ),
          tab_tib %>% dplyr::mutate( MAF = AF3, ALT = as.character(ALT3) ) %>% dplyr::filter( !is.na(ALT) )
        ) %>%
          dplyr::filter( !is.na(ALT) ) %>%
          dplyr::mutate( 
            MAF = dplyr::case_when(
              MAF == "." ~ "0",
              is.na(MAF) ~ "0",
              TRUE ~ MAF ) %>% as.double()
          ) %>% 
          dplyr::filter( MAF > min_maf ) %>%
          dplyr::arrange( Chromosome,Coordinate_SNP,SNP_ID, -MAF )
      }
      
    } else {
      tab_tib <- readr::read_tsv( tab_tsv )
    }

    # Default Return Type: Tabix Results
    ret_tib <- tab_tib
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Analyze Match Types::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( match_src ) {
      ret_tib <- dplyr::bind_rows(
        dplyr::inner_join( 
          tab_tib, tib %>% dplyr::mutate( Coordinate_SNP=Nxb_Pos_Up ),
          by=c("Chromosome","Coordinate_SNP"), suffix=c("_SNP","_CPG"),
          multiple = "all" ) %>% 
          dplyr::mutate( SNP_Match_Type="Nxb_Up"),
        
        dplyr::inner_join( 
          tab_tib, tib %>% dplyr::mutate( Coordinate_SNP=Cpg_Pos_Up ),
          by=c("Chromosome","Coordinate_SNP"), suffix=c("_SNP","_CPG"),
          multiple = "all"  ) %>% 
          dplyr::mutate( SNP_Match_Type="Cpg_Up"),
        
        dplyr::inner_join( 
          tab_tib, tib %>% dplyr::mutate( Coordinate_SNP=Cpg_Pos_Dn ),
          by=c("Chromosome","Coordinate_SNP"), suffix=c("_SNP","_CPG"),
          multiple = "all"  ) %>% 
          dplyr::mutate( SNP_Match_Type="Cpg_Dn"),
        
        dplyr::inner_join( 
          tab_tib, tib %>% dplyr::mutate( Coordinate_SNP=Nxb_Pos_Dn ),
          by=c("Chromosome","Coordinate_SNP"), suffix=c("_SNP","_CPG"),
          multiple = "all"  ) %>% 
          dplyr::mutate( SNP_Match_Type="Nxb_Dn")
      )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( write_out )
      out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv",
                             done = TRUE, write_spec = TRUE, append = FALSE,
                             fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Genomic Range Intersection Function::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersect_vars_GRS = function( ref,
                               can,
                               file = NULL,
                               
                               ref_key = NULL,
                               ref_col = NULL,
                               ref_prefix = NULL,
                               ref_red = TRUE,
                               
                               can_key = NULL,
                               can_col = NULL,
                               can_prefix = NULL, 
                               can_red = TRUE,
                               
                               out_dir,
                               run_tag,
                               pre_tag = NULL,
                               
                               reload     = 0,
                               reload_min = 2,
                               ret_data   = FALSE,
                               parallel   = FALSE,
                               
                               vb=0, vt=3, tc=1, tt=NULL,
                               fun_tag='intersect_vars_GRS' ) {
  
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
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  # is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
  #                               out_dir = out_dir,
  #                               vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  # 
  # if ( reload >= reload_min && is_valid )
  #   return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
  #                      vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
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
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c( out_csv, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( !is.null(file) ) {
    errs_mssg <- glue::glue("File file='{file}' does not exist")
    if ( !file.exists( file) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
  }
  
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
    
    can_tib <- as.data.frame(can) %>% tibble::as_tibble( rownames = can_key )
    # print( can_tib, n=20 )
    
    ref_tib <- as.data.frame(ref) %>% tibble::as_tibble(rownames=ref_key)
    # print( ref_tib, n=20 )
    
    tmp_int_tib <- intersect(  names(can_tib), names(ref_tib) )
    print(tmp_int_tib)
    
    can_tib <- can_tib %>% dplyr::select( - dplyr::all_of( tmp_int_tib ) )
    ref_tib <- ref_tib %>% dplyr::select( - dplyr::all_of( tmp_int_tib ) )
    
    ret_tib <- dplyr::bind_cols(
      can_tib[map_tib$queryHits, ],
      ref_tib[map_tib$subjectHits,]
    )
    # print(ret_tib)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if ( vb >= vt )
      cat(glue::glue("{mssg} Pre ret_cnt = '{ret_cnt}'.{RET}"))
    
    ret_cnt <- ret_tib %>% 
      dplyr::distinct( Rsn,Cgn, .keep_all = TRUE ) %>% 
      base::nrow()
    
    if ( vb >= vt )
      cat(glue::glue("{mssg} Post ret_cnt = '{ret_cnt}'.{RET}"))
    
    if (FALSE) {
      
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Extract Candidate Matches::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( vb >= vt+2 )
        cat(glue::glue("{mssg} Extracting Candidate Data...{RET}"))
      
      if (can_red) can_tib <- can_tib %>% 
          dplyr::mutate( seqnames=as.character(seqnames),
                         strand=as.character(strand)) %>%
          dplyr::rename( chr=seqnames,
                         pos=start,
                         srd=strand)
      
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
      
      ref_key <- glue::glue("ref-tib")
      ref_cnt <- print_tib( ref_tib, fun_tag = fun_tag, name = ref_key, 
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                     Bind Reference/Candidate Matches::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( vb >= vt+2 )
        cat(glue::glue("{mssg} Binding Candidate/Reference Data...{RET}"))
      
      ret_tib <- dplyr::bind_cols(
        can_tib[map_tib$queryHits, ],
        ref_tib[map_tib$subjectHits,]
      )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                  Generate Reference Coverage Summary::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # if ( vb >= vt+2 )
      #   cat(glue::glue("{mssg} Merging back full Reference...{RET}"))
      # ref_cov_tib <- ret_tib %>% dplyr::right_join( ref_tib )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv",
                           done = TRUE, write_spec = TRUE, append = FALSE,
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Cpp Improbe Subset IO Function::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_infinium_design = function( tib, 
                                scrU_key = "Score_U",
                                scrM_key = "Score_M", 
                                strand_CO = "Strand_CO",
                                cpg_count = "Cpg_Count",
                                vb=0, vt=3, tc=1, tt=NULL,
                                fun_tag='load_improbe_subset' ) {
  
  scrU_sym <- rlang::sym(scrU_key)
  scrM_sym <- rlang::sym(scrM_key)
  strand_CO_sym <- rlang::sym(strand_CO)
  cpg_count_sym <- rlang::sym(cpg_count)
  
  ret_tib <- NULL
  ret_tib <- tib %>%
    dplyr::mutate(
      Score_Class = base::round(10 * pmin(!!scrU_sym, !!scrM_sym), 0 ) %>% 
        as.integer(),
      Infinium_Design = dplyr::case_when(
        Score_Class < 3 & !!strand_CO_sym == "C" ~ 0,
        Score_Class < 2 & !!strand_CO_sym == "O" ~ 0,
        
        !!cpg_count_sym > 3 ~ 1,
        
        Score_Class < 4 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 1 ~ 1,
        Score_Class < 5 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 2 ~ 1,
        Score_Class < 6 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 3 ~ 1,
        
        Score_Class < 4 & !!strand_CO_sym == "O" & !!cpg_count_sym <= 1 ~ 1,
        Score_Class < 5 & !!strand_CO_sym == "O" & !!cpg_count_sym <= 2 ~ 1,
        Score_Class < 6 & !!strand_CO_sym == "O" & !!cpg_count_sym <= 3 ~ 1,
        
        Score_Class >= 3 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 0 ~ 2,
        Score_Class >= 4 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 1 ~ 2,
        Score_Class >= 5 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 2 ~ 2,
        Score_Class >= 6 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 3 ~ 2,
        
        Score_Class >= 2 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 0 ~ 2,
        Score_Class >= 4 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 1 ~ 2,
        Score_Class >= 5 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 2 ~ 2,
        Score_Class >= 6 & !!strand_CO_sym == "C" & !!cpg_count_sym <= 3 ~ 2,
        
        TRUE ~ -1
      ) %>% as.integer()
    )
  
  ret_tib
}

load_improbe_subset = function( file,
                                vb=0, vt=3, tc=1, tt=NULL,
                                fun_tag='load_improbe_subset') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}      file = '{file}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  errs_mssg <- glue::glue("File '{file}' does not exist")
  if ( !file.exists(file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  ftime <- base::system.time({
    
    ret_tib <- safe_read( file = file, vb=vb, tt=tt ) %>% 
      dplyr::rename(
        Loci_ID = Seq_ID,
        
        Design_Seq  = Methyl_Probe_Covered_Top_Sequence,
        Probe_Seq_U = UnMethyl_Probe_Sequence,
        Probe_Seq_M = Methyl_Probe_Sequence,
        
        Strand_FR = Methyl_Allele_FR_Strand,
        Strand_TB = Methyl_Allele_TB_Strand,
        Strand_CO = Methyl_Allele_CO_Strand,
        Next_Base = Methyl_Next_Base,
        
        Score_U   = UnMethyl_Final_Score,
        Score_M   = Methyl_Final_Score,
        Cpg_Count = Methyl_Underlying_CpG_Count
      ) %>%
      dplyr::mutate(
        Chromosome = 
          paste0("chr", stringr::str_remove( Chromosome, "^chr" ) ),
        Coordinate = as.integer(Coordinate),
        Strand_TB  = Strand_TB %>% stringr::str_sub(1,1),
        Probe_49M_U = Probe_Seq_U %>% stringr::str_sub(1,49),
        Probe_49M_M = Probe_Seq_M %>% stringr::str_sub(1,49),
        Probe_48M_U = Probe_Seq_U %>% stringr::str_sub(2,49),
        Probe_48M_M = Probe_Seq_M %>% stringr::str_sub(2,49),
        Probe_Type  = stringr::str_sub(Loci_ID, 1,2) ) %>%
      add_infinium_design( scrU_key = "Score_U",
                           scrM_key = "Score_M", strand_CO = "Strand_CO", cpg_count = "Cpg_Count", 
                           vb=vb, vt=vt+1, tc=tc+1, tt=tt )
    
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
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Original r_improbe Function::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

r_improbe_original = function(
    tib,
    
    ids_key = "Loci_ID",
    seq_key = "Top_Sequence",
    din_key = "Probe_Type",
    srd_str = "TB",
    
    src_org_path = "/Users/bretbarnes/Documents/tools/backup/Infinium_Methylation_Workhorse.04012020/scripts/R/probe_design/functions/improbe_functions.R",
    src_cur_path = "/Users/bretbarnes/Documents/tools/Workhorse-Unstained/scripts/R/trifecta/functions/probe_design_functions.R",
    
    out_dir,
    run_tag,
    pre_tag = NULL,
    
    reload     = 0,
    reload_min = 2,
    ret_data   = FALSE,
    parallel   = FALSE,
    
    vb=0, vt=3, tc=1, tt=NULL,
    fun_tag='r_improbe_original') {
  
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
    cat(glue::glue("{mssg}     func_tag = '{func_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      ids_key = '{ids_key}'.{RET}"))
    cat(glue::glue("{mssg}      seq_key = '{seq_key}'.{RET}"))
    cat(glue::glue("{mssg}      din_key = '{din_key}'.{RET}"))
    cat(glue::glue("{mssg}      srd_str = '{srd_str}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      pre_tag = '{pre_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_tag = '{out_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{mssg} src_org_path = '{src_org_path}'.{RET}"))
    cat(glue::glue("{mssg} src_cur_path = '{src_cur_path}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c( out_csv, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # errs_mssg <- glue::glue("File file='{file}' does not exist")
  # if ( !file.exists( file) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("File src_org_path='{src_org_path}' does not exist")
  if ( !file.exists( src_org_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("File src_cur_path='{src_cur_path}' does not exist")
  if ( !file.exists( src_cur_path) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    ids_sym <- rlang::sym( ids_key )
    seq_sym <- rlang::sym( seq_key )
    din_sym <- rlang::sym( din_key )
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    source(src_org_path)
    
    if ( vb >= 1 )
      cat(glue::glue("{pmssg} Building r_imp_des_TB_tib...{RET}"))
    
    # TBD:: Add an option to force unique::
    #
    ret_tib <- desSeq_to_prbs( 
      tib = tib, # %>% dplyr::distinct( !!ids_sym, !!seq_sym, !!din_sym ),
      idsKey = ids_key, 
      seqKey = seq_key,
      prbKey = din_key,
      strsSR = srd_str,
      addMatSeq = TRUE, 
      parallel  = parallel,
      verbose   = vb, tt=tt )
    
    # source(src_cur_path)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
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
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Binary IO Benchmarking::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

binary_IO_test_rcpp = function(top, cgn,
                               
                               all_int,
                               sep_val,
                               
                               cgn_max = 0,
                               
                               out_dir,
                               run_tag,
                               pre_tag = NULL,
                               
                               reload     = 0,
                               reload_min = 2,
                               ret_data   = FALSE,
                               parallel   = FALSE,
                               
                               vb=0, vt=3, tc=1, tt=NULL,
                               fun_tag='binary_IO_test_rcpp') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      all_int = '{all_int}'.{RET}"))
    cat(glue::glue("{mssg}      sep_val = '{sep_val}'.{RET}"))
    cat(glue::glue("{mssg}      cgn_max = '{cgn_max}'.{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      pre_tag = '{pre_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    all_val = TRUE
    if ( all_int == 0 ) all_val = FALSE
    
    out_val = ""
    if (  all_val && sep_val == 'b' ) out_val <- "stack.bin"
    if ( !all_val && sep_val == 'b' ) out_val <- "index.bin"
    if (  all_val && sep_val != 'b' ) out_val <- "stack.tsv"
    if ( !all_val && sep_val != 'b' ) out_val <- "index.tsv"
    cat( glue::glue("{pmssg}{TAB} out_val='{out_val}'{RET}") )
    
    if ( out_val %>% length() ==  0 ) {
      cat( glue::glue("{pmssg} Out Val is length zero = '{out_val}'. Next ...{RET}") )
      next
    }
    
    stack_file <- file.path(out_dir, paste(run_tag, out_val, sep='.') )
    index_file <- file.path(out_dir, paste(run_tag, out_val, "idx", sep='.') )
    
    both = TRUE
    
    if (both) {
      ret_tib <- init_cgDb( top_vec = top, # %>% head(),
                            cgn_vec = cgn, # %>% head(),
                            
                            path = out_dir,
                            name = run_tag,
                            
                            all = all_val,
                            sep = sep_val,
                            
                            cgn_max = cgn_max,
                            
                            trim = TRUE,
                            vb = VB_VAL, vt = VT_VAL )
    }
    stack_present <- stack_file %>% file.exists()
    
    if ( stack_present ) {
      cat( glue::glue("{RET2}{pmssg} Loading stack_file = {stack_file}.{RET2}") )
      
      if (both) {
        new_stack_TB <- load_cgDb_cpp( file = stack_file,
                                       vb = VB_VAL, vt = VT_VAL )
      }
    } else {
      cat( glue::glue("{pmssg} Failed to find stack_file = {stack_file}.{RET}") )
      next
    }
    
    top_length <- ret_tib %>% length()
    
    if ( vb >= 3 )
    {
      cat( glue::glue("{pmssg}      stack_file = {stack_file}.{RET}") )
      cat( glue::glue("{pmssg}      top_length = {top_length}.{RET}") )
      cat( glue::glue("{pmssg}         out_val = {out_val}.{RET}") )
      cat( glue::glue("{pmssg} head(ret_tib) ...{RET}") )
      ret_tib %>% head() %>% print()
      cat( glue::glue("{pmssg} tail(ret_tib) ...{RET}") )
      ret_tib %>% tail() %>% print()
      cat( glue::glue("{pmssg}{RET2}") )
      cat( glue::glue("# - ----- ----- ----- -----| Done |----- ----- ----- -----{RET2}") )
    }
    
    # warn_mssg <- glue::glue("WARN_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # wflag <- FALSE
    # 
    # errs_mssg <- glue::glue("ERROR_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    # if ( eflag ) return(NULL)
    
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
  
  ret_tib
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            File IO Functions::
#  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Local Run Time Defaults::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

stable_options = function( pars, args,
                           
                           vb=0, vt=4, tc=1, tt=NULL,
                           fun_tag='stable_options') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  eflag <- FALSE
  wflag <- FALSE
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb>=vt+2) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   args={args}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Parse Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts <- base::list()
  
  # opts$run_name     <- NULL
  opts$out_path     <- NULL
  
  opts$sample_max   <- 0
  opts$platform     <- "EPIC"
  
  opts$ref_path     <- NULL
  opts$ref_file     <- NULL
  opts$ref_build    <- NULL
  opts$ref_source   <- NULL
  opts$ref_species  <- NULL
  
  opts$canonical_csv <- file.path( pars$aux_path, "canonical_cgn_top_grp.csv.gz")
  opts$temp_off_csv  <- file.path( pars$prgm_aux_path, "template_offset.csv.gz")
  
  opts$Rscript      <- "Rscript"
  opts$single       <- FALSE
  opts$cluster      <- FALSE
  opts$parallel     <- FALSE
  opts$track_time   <- TRUE
  opts$clean        <- FALSE
  opts$reload       <- 0
  opts$verbose      <- 3
  
  if (opts$verbose > 0)
    cat(glue::glue("[{pars$prgm_tag}]: Starting; {pars$prgm_tag}.{RET2}"))
  
  if (pars$run_mode == 'RStudio') {
    
    pars$local_run_type <- NULL
    
    # Read Options and Test Cases from pre-defined function
    # loc_dat <- lighthoof_local_defaults( opts = opt,
    #                                      pars = par,
    #                                      vb = opts$verbose )
    # opt <- loc_dat$opt
    # par <- loc_dat$par
    
    pars$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/imSuite/scripts/R") %>%
      stringr::str_remove("/tools/imProbeQC/scripts/R") %>%
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R") %>%
      stringr::str_remove("/tools/git-repos/imSuite/scripts/R")
    opts$top_path <- pars$top_path
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Params Defaults::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( is.null(pars$version) ) pars$version <- '1'
    if ( is.null(pars$version_key) ) pars$version_key <- "v"
    opts$version  <- paste0(pars$version_key,pars$version)
    
    opts$out_path  <- file.path( pars$top_path, 'scratch' )
    opts$imp_path  <- file.path( pars$top_path, 'data/improbe' )
    opts$ann_path  <- file.path( pars$top_path, 'data/annotation' )
    opts$man_path  <- file.path( pars$top_path, 'data/manifests' )
    opts$idat_path <- file.path( pars$top_path, 'data/idats' )
    
    if (  is.null(pars$run_name) ) opts$run_name = "SNP"
    if ( !is.null(pars$run_name) ) opts$run_name = pars$run_name
    
    # opts$ref_source <- paste( "NCBI", sep = ',' )
    opts$ref_source <- paste( "UCSC", sep = ',' )
    
    opts$run_name <- paste(opts$run_name,opts$ref_source,opts$version, sep='-')
    if ( opts$ref_source == "UCSC" ) {
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$chr_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "hg38.fa.gz",
                              "hg19.fa.gz",
                              "hg18.fa.gz",
                              "mm10.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 sep = ',' )
      
      opts$ref_build <- paste( "hg38",
                               "hg19",
                               "hg18",
                               "mm10",
                               sep = ',' )
      
      
    } else if ( opts$ref_source == "NCBI" ) {
      
      # opts$ref_path <- paste(
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta' ),
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta' ),
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh36/Sequence/WholeGenomeFasta' ),
      #   file.path( pars$top_path, 'data/imGenomes/Mus_musculus/NCBI/GRCm10/Sequence/WholeGenomeFasta' ),
      #   sep = ',' )
      # 
      # opts$chr_path <- paste(
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes' ),
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/Chromosomes' ),
      #   file.path( pars$top_path, 'data/imGenomes/Homo_sapiens/NCBI/GRCh36/Sequence/Chromosomes' ),
      #   file.path( pars$top_path, 'data/imGenomes/Mus_musculus/NCBI/GRCm10/Sequence/Chromosomes' ),
      #   sep = ',' )
      
      opts$ref_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$chr_path <- paste(
        file.path( pars$top_path, 'data/imGenomes' ),
        sep = ',' )
      
      opts$ref_file <- paste( "GRCh38.genome.fa.gz",
                              "GRCh37.genome.fa.gz",
                              "GRCh36.genome.fa.gz",
                              "GRCm38.genome.fa.gz",
                              sep = ',' )
      
      opts$ref_species <- paste( "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Homo_sapiens",
                                 "Mus_musculus",
                                 sep = ',' )
      
      opts$ref_build <- paste( "GRCh38",
                               "GRCh37",
                               "GRCh36",
                               "GRCm10",
                               sep = ',' )
      
    } else {
      eflag <- TRUE
      errs_mssg <- glue::glue("Unsupported default ref_source = {ref_source}")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    # Define tru-cgn-map csv::
    opts$cgn_tru_csv <- file.path( 
      pars$top_path, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-top/bins-100/map.csv.gz" )
    
    # opts$single       = TRUE
    opts$single       = FALSE
    
    opts$cluster      = FALSE
    opts$parallel     = TRUE
    opts$track_time   = TRUE
    opts$clean        = FALSE
    opts$reload       = 1
    opts$verbose      = 3
    
    opts$workflow = "parse_genome,parse_chromosome,update_database"
    
    #
    # ewas_loci_variation::
    #
    if ( pars$prgm_tag == "ewas_loci_variation" ||
         pars$prgm_tag == "ewas_swifthoof" ||
         pars$prgm_tag == "ewas_loci_variation_fileBased" ||
         pars$prgm_tag == "stable_loci_variation_anaysis_scratch" ||
         pars$prgm_tag == "stable_registration_score_scratch" ) {
      
      opts$sample_max <- 0
      opts$platform <- "EPIC"
      opts$manifest <- file.path( opts$top_path, "Projects/EWAS/data/manifests/EWAS_PQC122021-NA-NA-GRCh37_sesame.beta.historic-beadPool.csv.gz" )
      opts$sample_sheet <- file.path( opts$top_path, "Projects/EWAS/data/sample_sheets/EWAS-alpha-cross-product-testing.07012022.sampleSheet.csv.gz" )
    }
    if ( pars$prgm_tag == "ewas_loci_variation"  ||
         pars$prgm_tag == "ewas_loci_variation_fileBased" ||
         pars$prgm_tag == "stable_loci_variation_anaysis_scratch" ) {
      opts$sdf_path  <- file.path( opts$top_path, "scratch/stable_ewas_swifthoof_scratch/sesame-UCSC-v3/prefix_to_sdf" )
      # opts$sdf_path  <- file.path( opts$top_path, "scratch/stable_loci_variation_R_scratch/CEPH/delta_beta/CEPH" )
    }
    
  } else if (pars$run_mode == 'Command_Line') {
    
    options_list <- list(
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Run Time Version Options:: 
      #                       Platform, Genome Build, etc
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--run_name"), type="character", default=opts$run_name, 
        help=glue::glue(
          "Build run name.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--workflow"), type="character", default=opts$workflow, 
        help=glue::glue(
          "Workflow(s) to be executed.",
          "{RET}{TAB2} e.g. parse_genome,parse_chromosome,update_database, etc.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--top_path"), type="character", default=opts$top_path,
        help=glue::glue(
          "Top directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--out_path"), type="character", default=opts$out_path,
        help=glue::glue(
          "Build output directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                 Options:: Loci Variation/EWAS-Swifthoof
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--sdf_path"), type="character", default=opts$sdf_path,
        help=glue::glue(
          "Pre-built SDF path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--idat_path"), type="character", default=opts$idat_path,
        help=glue::glue(
          "Idats path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--sample_sheet"), type="character", default=opts$sample_sheet,
        help=glue::glue(
          "Sample Sheet path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--manifest"), type="character", default=opts$manifest,
        help=glue::glue(
          "Manifest path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--platform"), type="character", default=opts$platform,
        help=glue::glue(
          "Infinium Methylation Platform path (i.e. EPIC, HM450, MM285, etc.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--sample_max"), type="integer", default=opts$sample_max, 
        help=glue::glue(
          "Maximum Number of Samples to Process.{RET}",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Program Specific Options:
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--ref_path"), type="character", default=opts$ref_path, 
        help=glue::glue(
          "Reference Genome directory path(s).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--chr_path"), type="character", default=opts$chr_path, 
        help=glue::glue(
          "Reference Genome directory path(s) for individual chromosomes.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_file"), type="character", default=opts$ref_file, 
        help=glue::glue(
          "Reference Genome fasta file name(s).",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_build"), type="character", default=opts$ref_build, 
        help=glue::glue(
          "Reference Genome build names(s).",
          "{RET}{TAB2} e.g. GRch38,GRCh37,GRCh36,GRCm37.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_source"), type="character", default=opts$ref_source, 
        help=glue::glue(
          "Reference Source(s).",
          "{RET}{TAB2} e.g. UCSC, NCBI.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--ref_species"), type="character", default=opts$ref_species, 
        help=glue::glue(
          "Reference Specie(s).",
          "{RET}{TAB2} e.g. Homo_sapiens, Mus_musculus, Rattus_norvegicus, ",
          "SARS-CoV-2.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--canonical_csv"), type="character", default=opts$canonical_csv, 
        help=glue::glue(
          "CSV file(s) containg canonical (already defined) CG#/Top Sequence ",
          "Template definitions.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--temp_off_csv"), type="character", default=opts$temp_off_csv, 
        help=glue::glue(
          "CSV file containg Top Sequence Template start and end offsets from ",
          "CG# genomic postion.",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Run Time Mode Options::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      optparse::make_option(
        c("--Rscript"), type="character", default=opts$Rscript,
        help=glue::glue(
          "Rscript executable path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      # Process Parallel/Cluster Parameters::
      optparse::make_option(
        c("--single"), action="store_true", default=opts$single, 
        help=glue::glue(
          "Boolean variable to run a single sample on a single-core.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--parallel"), action="store_true", default=opts$parallel, 
        help=glue::glue(
          "Boolean variable to run parallel on multi-core.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--cluster"), action="store_true", default=opts$cluster,
        help=glue::glue(
          "Boolean variable to run jobs on cluster.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      # Run=time Options::
      optparse::make_option(
        c("--track_time"), action="store_true", default=opts$track_time,
        help=glue::glue(
          "Boolean variable tack run times.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--clean"), action="store_true", default=opts$clean, 
        help=glue::glue(
          "Boolean variable to run a clean build.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="boolean"),
      
      optparse::make_option(
        c("--reload"), type="integer", default=opts$reload, 
        help=glue::glue(
          "Integer value to reload intermediate files (for testing).{RET}",
          "{TAB2} Zero indicates no-reloading, higher numbers more reloads.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer"),
      
      # Verbosity level::
      optparse::make_option(
        c("-v", "--verbose"), type="integer", default=opts$verbose, 
        help=glue::glue(
          "Verbosity level: 0-5 (5 is very verbose).",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="integer")
    )
    
    opt_parser = optparse::OptionParser( option_list = options_list )
    opts = optparse::parse_args( opt_parser )
    
    if ( is.null(opts$top_path) && !is.null(pars$top_path) )
      opts$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R")
    
  } else {
    stop( glue::glue("{RET}[{pars$prgm_tag}]: ERROR: Unrecognized run_mode = ",
                     "'{pars$run_mode}'!{RET2}") )
  }
  
  #
  # Parameter Validation:: ewas_loci_variation
  #
  if ( pars$prgm_tag == "ewas_loci_variation" ||
       pars$prgm_tag == "ewas_swifthoof"  ||
       pars$prgm_tag == "ewas_loci_variation_fileBased" ||
       pars$prgm_tag == "stable_loci_variation_anaysis_scratch" ) {
    
    if ( !file.exists(opts$manifest) ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("Manifest does not exist: '{opts$manifest}'")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
    if ( !file.exists(opts$sample_sheet) ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("Sample Sheet does not exist: '{opts$sample_sheet}'")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
    
  }
  if ( pars$prgm_tag == "ewas_loci_variation" ) {
    if ( !dir.exists(opts$sdf_path) ) {
      eflag <- TRUE
      errs_mssg <- glue::glue("SDF Directory does not exist: '{opts$sdf_path}'")
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return( NULL )
    }
  }
  
  if ( opts$clean ) opts$reload <- -1
  
  ret_cnt <- length(opts)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

# End of file

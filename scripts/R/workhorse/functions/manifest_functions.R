
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Manifest Methods::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Options, Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

# Additional Tidy Practices Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("stringr", quietly = TRUE) ))
suppressWarnings(suppressPackageStartupMessages( 
  base::require("glue",    quietly = TRUE) ) )

# Matrix Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("matrixStats", quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("scales",      quietly = TRUE) ))

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

#
# TBD:: All these functions need to re-formalized with new function structure!!
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Manifest Comparison Methods::
#
#  manifest_column_agreement()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

manifest_column_agreement = function(
  a, b, source_a, source_b,
  join_cols = c("Probe_ID","Prb_Des","Prb_Din"),
  
  vb=0, vt=3, tc=1, tt=NULL,
  fun_tag='manifest_column_agreement') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}    source_a={source_a}.{RET}"))
    cat(glue::glue("{mssg}    source_b={source_b}.{RET}"))
    cat(glue::glue("{mssg}   join_cols={join_cols}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    int_cols <- 
      intersect( names(a), names(b) )
    
    dif_cols <-
      setdiff( names(a), names(b) )
    join_cols <- intersect(join_cols, int_cols)
    
    if (vb >= vt) {
      cat(glue::glue("{mssg} SetDiff Columns::{RET}"))
      cat(glue::glue("{mssg}   dif_cols={dif_cols}.{RET}"))
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} Intersect Columns::{RET}"))
      cat(glue::glue("{mssg}   dif_cols={dif_cols}.{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    int_tibs <- dplyr::inner_join(
      a %>% dplyr::select(dplyr::all_of(int_cols)),
      b %>% dplyr::select(dplyr::all_of(int_cols)),
      by=join_cols,
      suffix=c("_a","_b")
    )
    int_key <- glue::glue("int-tib({fun_tag})")
    int_cnt <- print_tib(int_tibs, fun_tag = fun_tag, name = int_key,
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    ret_tib <- NULL
    for (col_key in int_cols) {
      if (col_key %in% join_cols) next
      
      if (vb >= vt+1)
        cat(glue::glue("{mssg} Comparing column = {col_key}...{RET}"))
      
      colA_key <- paste0(col_key,"_a")
      colB_key <- paste0(col_key,"_b")
      
      colA_sym <- rlang::sym(colA_key)
      colB_sym <- rlang::sym(colB_key)
      
      man_tib <- int_tibs %>% 
        dplyr::select(dplyr::all_of( c(colA_key,colB_key) ) ) %>%
        dplyr::summarise(
          Mat_Cnt=sum(!!colA_sym == !!colB_sym, na.rm = TRUE),
          Mis_Cnt=sum(!!colA_sym != !!colB_sym, na.rm = TRUE),
          Total=sum(!is.na(!!colA_sym) & !is.na(!!colB_sym) ),
          Mat_Per=round(100*Mat_Cnt/Total, 2),
          Mis_Per=round(100*Mis_Cnt/Total, 2)
        ) %>% 
        dplyr::mutate(Column = col_key,
                      Source_A = source_a,
                      Source_B = source_b)
      
      cur_key <- glue::glue("cur-key({fun_tag})")
      cur_cnt <- print_tib(man_tib, fun_tag = fun_tag, name = cur_key, 
                           vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      
      ret_tib <- ret_tib %>% dplyr::bind_rows(man_tib)
    }
    
    ret_key <- glue::glue("ret-FIN({fun_tag})")
    ret_cnt <- print_tib(ret_tib, fun_tag = fun_tag, name = ret_key,
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Sesame Manifest IO Functions::
#
#  sesame_repo_address_workflow()
#  load_sesame_repo_address()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sesame_repo_address_workflow = function(
  manifests,
  
  out_dir,
  run_tag,
  pre_tag = NULL,
  
  reload     = 0,
  reload_min = 2,
  ret_data    = FALSE,
  write_local = FALSE,
  overwrite   = FALSE,
  
  vb=0, vt=3, tc=1, tt=NULL,
  fun_tag='sesame_repo_address_workflow') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  out_dir <- file.path( out_dir,fun_tag )
  # beg_txt <- file.path( out_dir, paste(run_tag, fun_tag,'begin.txt', sep='.') )
  # sum_csv <- file.path( out_dir, paste(run_tag, fun_tag,'sum.csv.gz', sep='.') )
  # aux_csv <- file.path( out_dir, paste(run_tag, fun_tag,'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(run_tag, fun_tag,'csv.gz', sep='.') )
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  is_valid <- valid_time_stamp( c(pre_tag, out_csv, end_txt ), 
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload != 0 && reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag=glue::glue("Reloading({fun_tag})"), 
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}       pre_tag = '{pre_tag}'.{RET}"))
    cat(glue::glue("{mssg}        reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}    reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}      is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}   write_local = '{write_local}'.{RET}"))
    cat(glue::glue("{mssg}     overwrite = '{overwrite}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}       out_dir = '{out_dir}'.{RET}"))
    # cat(glue::glue("{mssg}       beg_txt = '{beg_txt}'.{RET}"))
    # cat(glue::glue("{mssg}       sum_csv = '{sum_csv}'.{RET}"))
    # cat(glue::glue("{mssg}       aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}       out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}       end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  e_time <- 0
  f_time <- 0
  f_time <- base::system.time({
    
    sesame_address_list <- get_file_list(files=manifests,
                                         alpha_numeric = TRUE, del = COM)
    found_all <- TRUE
    for (sesame_key in names(sesame_address_list)) {
      cur_fns <- 
        paste(run_tag, sesame_key, 'load_sesame_repo_address',"csv.gz", sep='.')
      cur_dir <- file.path(out_dir, sesame_key, 'load_sesame_repo_address')
      cur_csv <- file.path(cur_dir, cur_fns)
      if (vb >= vt + 3)
        cat(glue::glue("{mssg} Checking for intermediate file = ",
                       "{cur_csv}...{RET}"))
      
      if (!file.exists(cur_csv)) found_all <- FALSE
    }
    if (vb >= vt + 2)
      cat(glue::glue("{mssg} found_all intermediate files = {found_all}.{RET}"))
    
    sesame_address_dat <- NULL
    if (!found_all || reload==0) {
      
      for (sesame_key in names(sesame_address_list)) {
        cur_out_dir <- file.path(out_dir, sesame_key)
        cur_run_tag <- paste(run_tag, sesame_key, sep='.')
        min_reload <- reload_min - 1
        if (min_reload < 1) min_reload = 1
        
        sesame_address_dat[[sesame_key]] <- load_sesame_repo_address(
          sesame_address_list[[sesame_key]],
          add_decoy = TRUE,
          add_masks = TRUE,
          
          out_dir = cur_out_dir,
          run_tag = cur_run_tag,
          pre_tag = pre_tag,
          
          reload = reload,
          reload_min = min_reload,
          
          vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      }
    } else {
      for (sesame_key in names(sesame_address_list)) {
        cur_fns <- 
          paste(run_tag, sesame_key, 'load_sesame_repo_address',"csv.gz", sep='.')
        cur_dir <- file.path(out_dir, sesame_key, 'load_sesame_repo_address')
        cur_csv <- file.path(cur_dir, cur_fns)
        
        sesame_address_dat[[sesame_key]] <-
          safe_read(cur_csv, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Calculate Comparisons::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # TBD:: Use manifest_column_agreement() to compare manifests and write
    #   reuslts to summary file...
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Merge All Manifests::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- 
      dplyr::bind_rows(sesame_address_dat, .id = "Manifest_Source") %>%
      dplyr::distinct(Probe_ID, Prb_Des, .keep_all = TRUE)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (overwrite || ( write_local && !file.exists(out_csv) ) ) 
      out_cnt <- safe_write( ret_tib, out_csv, done = TRUE, write_spec = TRUE, 
                             vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib({fun_tag})")
    ret_cnt <- 
      print_tib( ret_tib, func=fun_tag, name = ret_key, vb,vt=vt+4,tc=tc+1 )
  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

load_sesame_repo_address = function(
  name,
  add_decoy = FALSE,
  add_masks = FALSE,
  normalize = TRUE,
  old_cols = c("Probe_ID", "seqnames", "start", "end", "strand", "designType", 
               "channel", "nextBase", "nextBaseRef", "probeType", "gene",
               "gene_HGNC"),
  new_cols = c("Probe_ID", "Chromosome", "Coordinate", "CoordinateG", 
               "Strand_FR", "Infinium_Design_Type", "Color", "Prb_Nxb", 
               "Prb_Nxb_Ref", "Prb_Din", "gene", "gene_HGNC"),
  
  out_dir,
  run_tag,
  pre_tag = NULL,
  
  reload     = 0,
  reload_min = 1,
  ret_data   = FALSE,
  
  vb=0, vt=3, tc=1, tt=NULL,
  fun_tag='load_sesame_repo_address') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  out_dir <- file.path( out_dir,fun_tag )
  beg_txt <- file.path( out_dir, paste(run_tag, fun_tag,'begin.txt', sep='.') )
  sum_csv <- file.path( out_dir, paste(run_tag, fun_tag,'sum.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(run_tag, fun_tag,'csv.gz', sep='.') )
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
                                vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  if ( reload != 0 && reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag=glue::glue("Reloading({fun_tag})"), 
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}        name = '{name}'.{RET}"))
    cat(glue::glue("{mssg}   add_decoy = '{add_decoy}'.{RET}"))
    cat(glue::glue("{mssg}   add_masks = '{add_masks}'.{RET}"))
    cat(glue::glue("{mssg}      reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}  reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}    is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}    ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}     out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}     beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}     sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}     out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}     end_txt = '{end_txt}'.{RET}"))
    if (vb >= vt+4) {
      cat(glue::glue("{mssg}{RET}"))
      cat(glue::glue("{mssg} Old Sesame Columns::{RET}"))
      cat(glue::glue("{mssg}    old_cols={old_cols}.{RET}"))
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} New Sesame Columns::{RET}"))
      cat(glue::glue("{mssg}    new_cols={new_cols}.{RET}"))
    }
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  e_time <- 0
  f_time <- 0
  f_time <- base::system.time({
    
    # Cache Sesame Manifest Data::
    data_key <- name %>% stringr::str_remove("\\..*$")
    sesameData::sesameDataCache(data_key)
    
    man_tib <- sesameData::sesameDataGet( name ) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var="Probe_ID") %>% 
      tibble::as_tibble()
    man_key <- glue::glue("man-tib({fun_tag})")
    man_cnt <- print_tib(man_tib, fun_tag = fun_tag, name = man_key, 
                         vb,vt=vt+4,tc=tc+1 )
    
    mask_cols <- NULL
    if (add_masks) mask_cols <- man_tib %>% 
      dplyr::select( dplyr::starts_with("MASK_") ) %>% names()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract/Format Probe A::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_key <- "A"
    old_prbA_cols <- man_tib %>% 
      dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
      dplyr::select(-dplyr::starts_with("wDecoy_")) %>% 
      names()
    new_prbA_cols <- old_prbA_cols %>% 
      stringr::str_remove(paste0("_",prb_key,"$")) %>% 
      stringr::str_to_title() %>%
      paste("Prb",., sep="_") %>%
      stringr::str_replace("Prb_Address", "Address") %>%
      stringr::str_replace("Prb_Probeseq", "Prb_Seq") %>%
      stringr::str_replace("Prb_Chrm", "Prb_Chr")
    
    old_decoyA_cols <- NULL
    new_decoyA_cols <- NULL
    if (add_decoy) {
      old_decoyA_cols <- man_tib %>%
        dplyr::select(dplyr::starts_with("wDecoy_")) %>%
        dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>%
        names()
      
      new_decoyA_cols <- old_decoyA_cols %>%
        stringr::str_remove(paste0("_",prb_key,"$")) %>% 
        stringr::str_remove("^wDecoy_") %>% 
        stringr::str_to_title() %>%
        stringr::str_replace("Chrm", "Chr") %>%
        paste("Prb_Decoy",., sep="_")
    }
    
    oldA_cols <- c(old_cols, old_prbA_cols, old_decoyA_cols, mask_cols)
    newA_cols <- c(new_cols, new_prbA_cols, new_decoyA_cols, mask_cols)
    
    oldA_cnt <- oldA_cols %>% length()
    newA_cnt <- newA_cols %>% length()
    if (vb >= vt+2)
      cat(glue::glue("{mssg} Old Sesame Columns:: Probe A({oldA_cnt}).{RET}"))
    if (vb >= vt+4)
      cat(glue::glue("{mssg}    oldA_cols={oldA_cols}.{RET}"))
    if (vb >= vt+2)
      cat(glue::glue("{mssg} New Sesame Columns:: Probe A({newA_cnt}).{RET}"))
    if (vb >= vt+4)
      cat(glue::glue("{mssg}    newA_cols={newA_cols}.{RET}"))
    
    addA_tib <- NULL
    addA_tib <- man_tib %>% 
      dplyr::select( dplyr::all_of(oldA_cols) ) %>%
      purrr::set_names( newA_cols ) %>% 
      dplyr::mutate( 
        Prb_Des=dplyr::case_when(
          Infinium_Design_Type=="I"  ~ "U",
          Infinium_Design_Type=="II" ~ "2",
          TRUE ~ NA_character_ )
      )
    
    addA_cnt <- addA_tib %>% base::nrow()
    if (vb  >= vt+2)
      cat(glue::glue("{mssg} Probe Set A: I/U Count={addA_cnt}.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract/Format Probe B::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_key <- "B"
    old_prbB_cols <- man_tib %>% 
      dplyr::filter( designType=="I" ) %>%
      dplyr::select( dplyr::ends_with(paste0("_",prb_key) ) ) %>% 
      dplyr::select( -dplyr::starts_with("wDecoy_" ) ) %>% 
      names()
    new_prbB_cols <- old_prbB_cols %>% 
      stringr::str_remove( paste0("_",prb_key,"$" ) ) %>% 
      stringr::str_to_title() %>%
      paste( "Prb",., sep="_" ) %>%
      stringr::str_replace( "Prb_Address", "Address" ) %>%
      stringr::str_replace( "Prb_Probeseq", "Prb_Seq" ) %>%
      stringr::str_replace( "Prb_Chrm", "Prb_Chr" )
    
    old_decoyB_cols <- NULL
    new_decoyB_cols <- NULL
    if (add_decoy) {
      old_decoyB_cols <- man_tib %>%
        dplyr::filter( designType=="I" ) %>%
        dplyr::select( dplyr::starts_with("wDecoy_" ) ) %>%
        dplyr::select( dplyr::ends_with(paste0("_",prb_key ) ) ) %>%
        names()
      
      new_decoyB_cols <- old_decoyB_cols %>%
        stringr::str_remove( paste0("_",prb_key,"$" ) ) %>% 
        stringr::str_remove( "^wDecoy_" ) %>% 
        stringr::str_to_title() %>%
        stringr::str_replace( "Chrm", "Chr" ) %>%
        paste( "Prb_Decoy",., sep="_" )
    }
    
    oldB_cols <- c(old_cols, old_prbB_cols, old_decoyB_cols, mask_cols)
    newB_cols <- c(new_cols, new_prbB_cols, new_decoyB_cols, mask_cols)
    
    oldB_cnt <- oldB_cols %>% length()
    newB_cnt <- newB_cols %>% length()
    if (vb  >= vt+2)
      cat(glue::glue( "{mssg} Old Sesame Columns:: Probe B({oldB_cnt}){RET}" ) )
    if (vb  >= vt+4)
      cat(glue::glue( "{mssg}    oldB_cols={oldB_cols}.{RET}" ) )
    if (vb  >= vt+2)
      cat(glue::glue( "{mssg} New Sesame Columns:: Probe B({newB_cnt}){RET}" ) )
    if (vb  >= vt+4)
      cat(glue::glue( "{mssg}    newB_cols={newB_cols}.{RET}" ) )
    
    addB_tib <- NULL
    addB_tib <- man_tib %>% 
      dplyr::filter( designType=="I" ) %>%
      dplyr::select( dplyr::all_of(oldB_cols ) ) %>%
      purrr::set_names( newB_cols ) %>% 
      dplyr::mutate( Prb_Des="M" )
    
    addB_cnt <- addB_tib %>% base::nrow()
    if (vb  >= vt+2)
      cat(glue::glue( "{mssg} Probe Set B: M Count={addB_cnt}.{RET2}" ) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Stack Probe Designs:: A/B
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- dplyr::bind_rows(addA_tib,addB_tib) %>%
      dplyr::mutate(
        Strand_FR=dplyr::case_when(
          Strand_FR=="+" ~ "F",
          Strand_FR=="-" ~ "R",
          TRUE ~ NA_character_),
        Prb_Nxb=stringr::str_remove_all(Prb_Nxb,"[^a-zA-Z]") %>% mapDIs(),
        Prb_Inf=dplyr::case_when(
          Infinium_Design_Type=="I"  ~ 1,
          Infinium_Design_Type=="II" ~ 2,
          TRUE ~ NA_real_
        ) %>% as.integer(),
        Color=dplyr::case_when(
          Color=="Both" ~ NA_character_,
          TRUE ~ Color)
      )  %>% dplyr::arrange(Probe_ID)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Summary and Sanity Checks::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # This better be zero::
    dup_cnt <- ret_tib %>% 
      dplyr::arrange(Address) %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (dup_cnt != 0) {
      stop(glue::glue("{mssg} ERROR: Duplicates Found = {dup_cnt}!{RET2}"))
      return(NULL)
    }
    if (vb  >= vt+1)
      cat(glue::glue("{mssg} No duplicates detected = {dup_cnt}.{RET}"))
    
    sum_tib <- ret_tib %>% 
      dplyr::group_by(Prb_Des) %>% 
      dplyr::summarise(Count=n(), .groups = "drop")
    
    if (vb  >= vt+2) {
      cat(glue::glue("{mssg} Probe Design Coverage={RET}"))
      sum_tib %>% print(n=base::nrow(sum_tib))
      cat(glue::glue("{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    sum_cnt <- safe_write( sum_tib, file = sum_csv, fun_tag = fun_tag,
                           vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    out_cnt <- safe_write( ret_tib, file = out_csv, fun_tag = fun_tag,
                           done = TRUE, write_spec = TRUE, 
                           vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib({fun_tag})")
    ret_cnt <- 
      print_tib( ret_tib, func=fun_tag, name = ret_key, vb,vt=vt+4,tc=tc+1 )
  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Genome Studio Manifest I/O::
#
#  genome_studio_to_mock_aqps()
#  load_genome_studio_address()
#  load_genome_studio_manifest()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

partition_Ewapic_Aqp = function(map_csv,
                                map_col,
                                out_csv,
                                
                                src_ord_dir,
                                src_mat_dir,
                                src_aqp_dir,
                                
                                out_aux_dir,
                                out_ord_dir,
                                out_mat_dir,
                                out_aqp_dir,
                                
                                out_dir,
                                run_tag,
                                pre_tag = NULL,
                                
                                reload     = 0,
                                reload_min = 2,
                                ret_data   = FALSE,
                                parallel   = FALSE,
                                
                                vb=0, vt=3, tc=1, tt=NULL,
                                fun_tag='partition_Ewapic_Aqp') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  out_dir <- file.path( out_dir, run_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )

  if ( is.null( out_csv ) )
    out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  ord_dir <- file.path(out_dir, "order")
  mat_dir <- file.path(out_dir, "match")
  aqp_dir <- file.path(out_dir, "aqpqc")
  
  safe_mkdir( out_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  safe_mkdir( ord_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  safe_mkdir( mat_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  safe_mkdir( aqp_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  
  ord_vec <- 
    sort_files_by_date( path = ord_dir, pattern = "\\.csv\\.gz$", ret = "vec")
  mat_vec <- 
    sort_files_by_date( path = mat_dir, pattern = "\\.match\\.gz$", ret = "vec")
  aqp_vec <- 
    sort_files_by_date( path = aqp_dir, pattern = "\\.txt\\.gz$", ret = "vec")

  is_valid <- valid_time_stamp( 
    c(pre_tag, beg_txt, ord_vec, mat_vec, aqp_vec, out_csv, end_txt ), 
    out_dir = out_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c( out_csv, end_txt ) )
  if ( length( ord_vec ) != 0 ) unlink( ord_vec )
  if ( length( mat_vec ) != 0 ) unlink( mat_vec )
  if ( length( aqp_vec ) != 0 ) unlink( aqp_vec )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  e_time <- 0
  f_time <- 0
  f_time <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Build Full Map Table::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    map_tib <- suppressMessages( suppressWarnings(
      readr::read_csv( map_csv ) ) ) %>% 
      purrr::set_names( map_col ) %>% 
      tibble::as_tibble() %>% 
      tidyr::pivot_longer(cols=c(AQP1_Num,AQP2_Num), 
                          names_to="AQP_Num", 
                          values_to = "AQP_Name") %>% 
      dplyr::filter(!is.na(AQP_Name)) %>% 
      dplyr::filter(stringr::str_starts(AQP_Name, "BS")) %>%
      dplyr::mutate(AQP_Num=AQP_Num %>%
                      stringr::str_remove("^AQP") %>% 
                      stringr::str_remove("_Num$") %>% 
                      as.integer(),
                    
                    Src_Order_Path=Order_Path %>% 
                      stringr::str_remove( "^.*\\\\" ) %>% 
                      paste0( src_ord_dir,'/',.,".gz" ),
                    Src_Match_Path=dplyr::case_when(
                      AQP_Num == 2 ~ paste0( 
                        src_mat_dir,"/AQP2-",Match_Num,"_probes.match.gz" ),
                      AQP_Num == 1 ~ paste0( 
                        src_mat_dir,"/",Match_Num,"_probes.match.gz" ),
                      TRUE ~ NA_character_ ),
                    Src_AQP_Path = paste0(
                      src_aqp_dir,"/",AQP_Name,".txt.gz" ),
                    
                    Out_Order_Path=Order_Path %>% 
                      stringr::str_remove( "^.*\\\\" ) %>% 
                      paste0( out_ord_dir,'/',.,".gz" ),
                    Out_Match_Path=dplyr::case_when(
                      AQP_Num == 2 ~ paste0( 
                        out_mat_dir,"/AQP2-",Match_Num,"_probes.match.gz" ),
                      AQP_Num == 1 ~ paste0( 
                        out_mat_dir,"/",Match_Num,"_probes.match.gz" ),
                      TRUE ~ NA_character_ ),
                    Out_AQP_Path = paste0( out_aqp_dir,"/",AQP_Name,".txt.gz" ),
                    
                    Order_File_Name = base::basename( Out_Order_Path ),
                    Match_File_Name = base::basename( Out_Match_Path ),
                    AQP_File_Name = base::basename( Out_AQP_Path ) ) %>%
      dplyr::select( Bead_Pool:Bucket_Name, AQP_Num, AQP_Name,
                     Order_File_Name ,Match_File_Name, AQP_File_Name,
                     Src_Order_Path, Src_Match_Path, Src_AQP_Path,
                     Out_Order_Path, Out_Match_Path, Out_AQP_Path ) %>%
      dplyr::arrange( AQP_Num, Bead_Pool )
    
    if (ret_data) ret_dat$map <- map_tib

    map_key <- glue::glue("built-map-tib")
    map_cnt <- print_tib( map_tib, fun_tag = fun_tag, name = map_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Write Updated Full Map Table::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    org_map_csv <- base::basename(map_csv) %>% 
      stringr::str_remove("\\.gz$") %>% 
      stringr::str_remove("\\.[ct]sv$") %>% 
      paste("formatted.csv.gz", sep='.')
    org_map_csv <- file.path( base::dirname(map_csv), org_map_csv )
    
    if ( vb >= vt+1 )
      cat(glue::glue("{mssg} Writing org_map_csv = '{org_map_csv}'...{RET}"))
    
    map_cnt <- safe_write( x = map_tib, file = org_map_csv, done = TRUE, 
                           write_spec = TRUE, append = FALSE, fun_tag = fun_tag, 
                           vb=vb,vt=vt+6,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Sanity Check:: 
    #                     Validate All Source Files Exist!!!
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    src_ord_val <-
      unlist( cbind( lapply(map_tib$Src_Order_Path, file.exists) ) ) %>% unique()
    src_mat_val <-
      unlist( cbind( lapply(map_tib$Src_Match_Path, file.exists) ) ) %>% unique()
    src_aqp_val <-
      unlist( cbind( lapply(map_tib$Src_AQP_Path, file.exists) ) ) %>% unique()
    
    if ( src_ord_val != TRUE || src_mat_val != TRUE || src_aqp_val != TRUE ) {
      fail_mssg <- glue::glue("File check failed on source input files")
      stop(glue::glue( "{error} src_ord_val = '{src_ord_val}'",
                       "{error} src_mat_val = '{src_mat_val}'",
                       "{error} src_aqp_val = '{src_aqp_val}'",
                       "{error} {fail_mssg}!{error} Exiting...{RET2}") )
      return(NULL)
    }
    if ( vb >= vt+1 )
      cat(glue::glue("{mssg} Passed Source File Sanity Check!{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Clean:: Output Directories
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # NOTE:: I'm pretty sure the unlink's at the begining of the function
    #  actually take care of this cleaning...
    # TBD:: Validate the code below is not needed...
    #
    list.files(out_aux_dir, pattern = "\\.csv\\.gz$", full.names = TRUE) %>% unlink()
    list.files(out_ord_dir, pattern = "\\.csv\\.gz$", full.names = TRUE) %>% unlink()
    list.files(out_mat_dir, pattern = "\\.match\\.gz$", full.names = TRUE) %>% unlink()
    list.files(out_aqp_dir, pattern = "\\.txt\\.gz$", full.names = TRUE) %>% unlink()
    
    if ( vb >= vt+1 )
      cat(glue::glue("{mssg} Cleaned Output Directories!{RET2}"))

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Select Target Map Table::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- NULL
    if ( run_tag == "Epic2" ) ret_tib <- map_tib %>% 
      dplyr::filter(stringr::str_detect(Bucket_Name, "EPIC") )
    if ( run_tag == "Ewas1" ) ret_tib <- map_tib %>% 
      dplyr::filter(stringr::str_detect(Bucket_Name, "EWAS") )
    
    if ( is.null(ret_tib) || base::nrow(ret_tib) == 0 ) {
      fail_mssg <- glue::glue("Failed to match any match tabble. Please ",
                              "check user run_name input")
      stop(glue::glue( "{error} {fail_mssg}!{error} Exiting...{RET2}") )
      return(NULL)
    }

    ret_key <- glue::glue("selected-ret-tib({run_tag})")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+5,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Copy:: Input to Output Directories
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    base::file.copy( 
      ret_tib$Src_Order_Path, ret_tib$Out_Order_Path, copy.date = TRUE )
    base::file.copy( 
      ret_tib$Src_Match_Path, ret_tib$Out_Match_Path, copy.date = TRUE )
    base::file.copy( 
      ret_tib$Src_AQP_Path, ret_tib$Out_AQP_Path, copy.date = TRUE )
    
    if ( vb >= vt+1 )
      cat(glue::glue("{mssg} Copied Selected Source Files to Output!{RET2}"))

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Sanity Check:: 
    #              Validate All Output Files Exist After Copying!!!
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_ord_val <-
      unlist( cbind( lapply(ret_tib$Out_Order_Path, file.exists) ) ) %>% unique()
    out_mat_val <-
      unlist( cbind( lapply(ret_tib$Out_Match_Path, file.exists) ) ) %>% unique()
    out_aqp_val <-
      unlist( cbind( lapply(ret_tib$Out_AQP_Path, file.exists) ) ) %>% unique()
    
    if ( out_ord_val != TRUE || out_mat_val != TRUE || out_aqp_val != TRUE ) {
      fail_mssg <- glue::glue("File check failed on new output files")
      stop(glue::glue( "{error} out_ord_val = '{out_ord_val}'",
                       "{error} out_mat_val = '{out_mat_val}'",
                       "{error} out_aqp_val = '{out_aqp_val}'",
                       "{error} {fail_mssg}!{error} Exiting...{RET2}") )
      return(NULL)
    }
    if ( vb >= vt+1 )
      cat(glue::glue("{mssg} Passed Output File Sanity Check!{RET2}"))

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Write Selected Map Table::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv,
                           done = TRUE, write_spec = TRUE, append = FALSE,
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    if (ret_data) ret_dat$ret <- ret_tib

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  if (ret_data) return(ret_dat)
  
  ret_tib
}


genome_studio_to_mock_aqps = function(man_dir,
                                      man_str,

                                      decoder_len = 45,
                                      
                                      out_dir,
                                      run_tag,
                                      pre_tag = NULL,
                                      
                                      reload     = 0,
                                      reload_min = 2,
                                      ret_data   = FALSE,
                                      parallel   = FALSE,
                                      
                                      vb=0, vt=3, tc=1, tt=NULL,
                                      fun_tag='genome_studio_to_mock_aqps') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      man_dir = '{man_dir}'.{RET}"))
    cat(glue::glue("{mssg}      man_str = '{man_str}'.{RET}"))
    cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  e_time <- 0
  f_time <- 0
  f_time <- base::system.time({
    
    ids_key  <- "IlmnID"
    addA_key <- "AddressA_ID"
    idsA_key <- "AlleleA_Probe_Id"
    prbA_key <- "AlleleA_ProbeSeq"
    addB_key <- "AddressB_ID"
    idsB_key <- "AlleleB_Probe_Id"
    prbB_key <- "AlleleB_ProbeSeq"
    norm_key <- "Color_Channel"
    bead_key <- "Beadpool_ID"
    
    ids_sym  <- rlang::sym(ids_key)
    addA_sym <- rlang::sym(addA_key)
    idsA_sym <- rlang::sym(idsA_key)
    prbA_sym <- rlang::sym(prbA_key)
    addB_sym <- rlang::sym(addB_key)
    idsB_sym <- rlang::sym(idsB_key)
    prbB_sym <- rlang::sym(prbB_key)
    norm_sym <- rlang::sym(norm_key)
    bead_sym <- rlang::sym(bead_key)
    
    # Define Order Columns/Header::
    ord_cols <- c( "Assay_Design_Id",
                   "AlleleA_Probe_Id", "AlleleA_Probe_Sequence",
                   "AlleleB_Probe_Id", "AlleleB_Probe_Sequence",
                   "Normalization_Bin" )
    
    # Define Match Header::
    mat_cols <- c( "probe_id", "sequence", "type_b", "address_name", "bo_seq" )
    
    # Define AQP Header::
    aqp_cols <- 
      c( "Address", "Decode_Status", "Decode_Error_Code", "Decode_Score", 
         "Func_Status", "Func_Error_Code", "QC_Action" )
    
    # Define PQC Header::
    pqc_cols <- 
      c( "Address", "Status", "Eval Code", "Average Rep", "Expected Rep" )
    
    man_fns <- as.list( file.path( man_dir, split_str_to_vec( man_str ) ) )
    names(man_fns) <- run_tag
    
    for (man_name in names( man_fns ) ) {
      if (vb >= vt) cat(glue::glue("{mssg} Starting {man_name}...{RET2}"))
      
      cur_dir <- file.path( out_dir, man_name )
      out_tag <- paste( man_name, fun_tag, sep='.' )
      out_csv <- file.path( cur_dir, paste(out_tag, sep='.') )
      beg_txt <- paste(out_csv, 'start.txt', sep='.')
      end_txt <- paste(out_csv, 'done.txt', sep='.')
      
      man_csv <- man_fns[[man_name]]
      cur_basename <- base::basename( man_csv ) %>% 
        stringr::str_remove("\\.gz$") %>%
        stringr::str_remove("\\.csv$")
      
      ord_dir <- file.path(cur_dir, "order")
      mat_dir <- file.path(cur_dir, "match")
      aqp_dir <- file.path(cur_dir, "aqpqc")
      
      safe_mkdir( out_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
      safe_mkdir( ord_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
      safe_mkdir( mat_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
      safe_mkdir( aqp_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
      
      ord_vec <- sort_files_by_date( 
        path = ord_dir, pattern = "\\.csv\\.gz$", ret = "vec")
      mat_vec <- sort_files_by_date( 
        path = mat_dir, pattern = "\\.tsv\\.gz$", ret = "vec")
      aqp_vec <- sort_files_by_date( 
        path = aqp_dir, pattern = "\\.tsv\\.gz$", ret = "vec")
      
      is_valid <- valid_time_stamp( 
        c( man_csv, pre_tag, beg_txt, ord_vec, mat_vec, aqp_vec, end_txt ), 
        out_dir = cur_dir, vb=vb,vt=vt+6,tc=tc+1,tt=tt )
      
      if ( reload >= reload_min && is_valid ) {
        if ( vb >= vt ) 
          cat(glue::glue("{mssg} Files exist ({man_name}). Skipping...{RET2}"))
        if ( vb >= vt+6 ) {
          cat(glue::glue("{mssg} Order Files ={RET}") )
          cat(glue::glue("{mssg}{TAB} '{ord_vec}'{RET}") )
          
          cat(glue::glue("{mssg} Match Files ={RET}") )
          cat(glue::glue("{mssg}{TAB} '{mat_vec}'{RET}") )
          
          cat(glue::glue("{mssg} AQP/PQC Files ={RET}") )
          cat(glue::glue("{mssg}{TAB} '{aqp_vec}'{RET}") )
          cat(glue::glue("{RET}") )
        }
        next
      }
      if (vb >= vt+3) {
        cat(glue::glue("{mssg} Function Parameters::{RET}"))
        cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
        cat(glue::glue("{mssg} File IO Parameters::{RET}"))
        cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
        cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
        cat(glue::glue("{mssg}      ord_dir = '{ord_dir}'.{RET}"))
        cat(glue::glue("{mssg}      mat_dir = '{mat_dir}'.{RET}"))
        cat(glue::glue("{mssg}      aqp_dir = '{aqp_dir}'.{RET}"))
        cat(glue::glue("{mssg}      man_csv = '{man_csv}'.{RET}"))
        cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
        cat(glue::glue("{RET}"))
      }
      
      unlink( c( end_txt ) )
      if ( length( ord_vec ) != 0 ) unlink( ord_vec )
      if ( length( mat_vec ) != 0 ) unlink( mat_vec )
      if ( length( aqp_vec ) != 0 ) unlink( aqp_vec )
      
      if ( !file.exists(beg_txt) )
        sys_ret <- base::system( glue::glue("touch {beg_txt}") )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                              Load and Format::
      #                       Original Genome Studio Manifest::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ann_csv <- file.path(
        man_dir, paste(cur_basename, "analytical.csv.gz", sep = '.') )
      
      man_tib <- NULL
      if ( file.exists( ann_csv ) ) {
        if ( vb >= vt ) cat(glue::glue( "{mssg} Will use anlaytical file = ",
                                        "'{ann_csv}'.{RET2}") )
        
        man_tib <- safe_read( file = ann_csv, use_spec = TRUE, 
                              spec_prefernce = "rds",
                              vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      } else {
        if ( vb >= vt )
          cat(glue::glue( "{mssg} Will use raw file = '{man_csv}'.{RET2}") )
        
        man_tib <- load_genome_studio_manifest( 
          file = man_csv, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      }
      
      if ( !bead_key %in% names(man_tib) ) {
        if ( vb >= vt+1 )
          cat(glue::glue( "{mssg} Adding artifical {bead_key} = 'BP1'.{RET2}"))
                          
        man_tib <- man_tib %>% dplyr::mutate( !!bead_sym := "BP1" )
      }
      if (ret_data) ret_dat$man <- man_tib
      
      man_key <- glue::glue("man-tib")
      man_cnt <- print_tib( man_tib, fun_tag = fun_tag, name = man_key,
                            vb=vb,vt=vt+4,tc=tc+1,tt=tt )

      bead_cnt <- man_tib %>% dplyr::distinct( !!bead_sym ) %>% base::nrow()
      if ( bead_cnt == 1 ) {
        if ( vb >= vt+1 )
          cat(glue::glue( "{mssg} Only {bead_cnt} {bead_key}. Switching ",
                          "parallel to FALSE.{RET2}"))
        parallel = FALSE
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Build Mock Order Files::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( vb >= vt+1 )
        cat(glue::glue( "{RET}{mssg}{TAB} Building Mock Order Files for each ",
                        "Beadpool...{RET}") )

      ord_tib <- NULL
      ord_tib <- man_tib %>% dplyr::select(
        dplyr::all_of( c(ids_key, prbA_key, prbB_key, norm_key, bead_key) ) ) %>%
        dplyr::mutate( !!idsA_sym := paste(!!ids_sym,"A", sep = '_' ),
                       !!idsB_sym := paste(!!ids_sym,"B", sep = '_' ),
                       !!norm_sym := dplyr::case_when(
                         is.na(!!norm_sym) ~ "C",
                         !!norm_sym == "Red" ~ "A",
                         !!norm_sym == "Grn" ~ "B",
                         TRUE ~ NA_character_ ) ) %>%
        dplyr::select( 
          dplyr::all_of( c(ids_key, idsA_key, prbA_key, idsB_key, prbB_key, 
                           norm_key, bead_key) ) ) %>%
        purrr::set_names( c( ord_cols, bead_key ) )
      
      if (ret_data) ret_dat$ord <- ord_tib

      ord_key <- glue::glue("ord-tib")
      ord_cnt <- print_tib( ord_tib, fun_tag = fun_tag, name = ord_key,
                            vb=vb,vt=vt+5,tc=tc+1,tt=tt )
      
      ord_tibs <- ord_tib %>% 
        base::split.data.frame( f = .[[ bead_key ]], drop = TRUE )
      
      if (ret_data) ret_dat$ords <- ord_tibs
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Write Mock Order Files::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( parallel ) {
        if ( vb >= vt+1 )
          cat(glue::glue( "{mssg}{TAB} Writing Mock Order Files for each ",
                          "Beadpool in parallel processing...{RET2}") )

        foreach (beadpool=names( ord_tibs ) ) %dopar% {
          ord_csv <- paste( cur_basename,beadpool,"csv.gz", sep = '.' )
          ord_csv <- file.path( ord_dir, ord_csv )
          
          cur_tib <- NULL
          cur_tib <- ord_tibs[[ beadpool ]] %>% 
            dplyr::select( -dplyr::all_of( c( bead_key ) ) )
          
          cur_key <- glue::glue("cur-ord-tib(BP = {beadpool})")
          cur_cnt <- print_tib( cur_tib, fun_tag = fun_tag, name = cur_key,
                                vb=vb,vt=vt+5,tc=tc+1,tt=tt )
          
          safe_write( x = cur_tib, file = ord_csv, done = TRUE, write_spec = TRUE,
                      vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        }
        
      } else {
        if ( vb >= vt+1 )
          cat(glue::glue( "{mssg}{TAB} Writing Mock Order Files for each ",
                          "Beadpool in linear processing...{RET2}") )

        for ( beadpool in names( ord_tibs ) ) {
          ord_csv <- paste(cur_basename,beadpool,"csv.gz", sep = '.')
          ord_csv <- file.path( ord_dir, ord_csv )
          
          cur_tib <- NULL
          cur_tib <- ord_tibs[[ beadpool ]] %>% 
            dplyr::select( -dplyr::all_of( c( bead_key ) ) )
          
          cur_key <- glue::glue("cur-ord-tib(BP = {beadpool})")
          cur_cnt <- print_tib( cur_tib, fun_tag = fun_tag, name = cur_key,
                                vb=vb,vt=vt+5,tc=tc+1,tt=tt )
          
          safe_write( x = cur_tib, file = ord_csv, done = TRUE, write_spec = TRUE,
                      vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        }
        
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Build Mock Match Files::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( vb >= vt+1 )
        cat(glue::glue( "{RET}{mssg}{TAB} Building Mock Match Files for each ",
                        "Beadpool...{RET}") )
      
      decoder_seq <- rep("N", decoder_len) %>% paste(collapse = "")
      
      mat_tib <- NULL
      mat_tib <- man_tib %>% 
        dplyr::select( IlmnID, 
                       AddressA_ID, AlleleA_ProbeSeq,
                       AddressB_ID, AlleleB_ProbeSeq, !!bead_sym ) %>% 
        dplyr::rename( Address_A = AddressA_ID, sequence_A = AlleleA_ProbeSeq,
                       Address_B = AddressB_ID, sequence_B = AlleleB_ProbeSeq
        ) %>%
        tidyr::pivot_longer( cols = c( Address_A, sequence_A,
                                       Address_B, sequence_B ),
                             names_to = c(".value", "AB"), 
                             names_sep = "_", 
                             values_drop_na = TRUE ) %>%
        dplyr::distinct() %>%
        dplyr::mutate( IlmnID = paste(IlmnID, AB, sep="_"),
                       type_b = "S",
                       bo_seq = paste0( decoder_seq, sequence ) ) %>%
        dplyr::rename( address_name = Address,
                       probe_id = IlmnID ) %>%
        dplyr::select( dplyr::all_of( c( mat_cols, bead_key ) ) )
      
      if (ret_data) ret_dat$mat <- mat_tib
      
      mat_key <- glue::glue("mat-tib")
      mat_cnt <- print_tib( mat_tib, fun_tag = fun_tag, name = mat_key,
                            vb=vb,vt=vt+5,tc=tc+1,tt=tt )
      
      mat_tibs <- NULL
      mat_tibs <- mat_tib %>% 
        base::split.data.frame( f = .[[ bead_key ]], drop = TRUE )
      
      if (ret_data) ret_dat$mats <- mat_tibs
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Write Mock Match Files::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( parallel ) {
        if ( vb >= vt+1 )
          cat(glue::glue( "{mssg}{TAB} Writing Mock Match Files for each ",
                          "Beadpool in parallel processing...{RET2}") )
        
        foreach (beadpool=names( mat_tibs ) ) %dopar% {
          mat_csv <- 
            paste(cur_basename, beadpool, "AQP-1.match.tsv.gz", sep = '.')
          mat_csv <- file.path( mat_dir, mat_csv )
          
          cur_tib <- NULL
          cur_tib <- mat_tibs[[ beadpool ]] %>% 
            dplyr::select( -dplyr::all_of( c( bead_key ) ) )
          
          cur_key <- glue::glue("cur-mat-tib(BP = {beadpool})")
          cur_cnt <- print_tib( cur_tib, fun_tag = fun_tag, name = cur_key,
                                vb=vb,vt=vt+5,tc=tc+1,tt=tt )
          
          safe_write( x = cur_tib, file = mat_csv, done = TRUE, 
                      write_spec = TRUE, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        }
        
      } else {
        
        if ( vb >= vt+1 )
          cat(glue::glue( "{mssg}{TAB} Writing Mock Match Files for each ",
                          "Beadpool in linear processing...{RET2}") )
        
        for ( beadpool in names( mat_tibs ) ) {
          mat_csv <- 
            paste(cur_basename, beadpool, "AQP-1.match.tsv.gz", sep = '.')
          mat_csv <- file.path( mat_dir, mat_csv )
          
          cur_tib <- NULL
          cur_tib <- mat_tibs[[ beadpool ]] %>% 
            dplyr::select( -dplyr::all_of( c( bead_key ) ) )
          
          cur_key <- glue::glue("cur-mat-tib(BP = {beadpool})")
          cur_cnt <- print_tib( cur_tib, fun_tag = fun_tag, name = cur_key,
                                vb=vb,vt=vt+5,tc=tc+1,tt=tt )
          
          safe_write( x = cur_tib, file = mat_csv, done = TRUE, 
                      write_spec = TRUE, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        }
        
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Build Mock PQC File::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( vb >= vt+1 )
        cat(glue::glue( "{RET}{mssg}{TAB} Building Single PQC File for ",
                        "Beadpool.{RET}"))

      pqc_tib <- NULL
      pqc_tib <- mat_tib %>% 
        dplyr::distinct( address_name, !!bead_sym ) %>% 
        dplyr::rename( Address = address_name ) %>% 
        dplyr::arrange( Address ) %>% 
        dplyr::mutate( Stat = 0, 
                       Eval = 0, 
                       Avg  = 300, 
                       Exp  = 0 ) %>% 
        dplyr::select( dplyr::all_of( 
          c( "Address", "Stat", "Eval", "Avg", "Exp", bead_key ) ) )
        purrr::set_names( c( pqc_cols, bead_key ) )
      
      pqc_key <- glue::glue("pqc-tib")
      pqc_cnt <- print_tib( pqc_tib, fun_tag = fun_tag, name = pqc_key,
                            vb=vb,vt=vt+5,tc=tc+1,tt=tt )
      
      aqp_tibs <- NULL
      aqp_tibs <- pqc_tib %>% 
        base::split.data.frame( f = .[[ bead_key ]], drop = TRUE )
      
      pqc_tib <- pqc_tib %>% 
        dplyr::select( -dplyr::all_of( c( bead_key ) ) )
      
      if (ret_data) ret_dat$pqc  <- pqc_tib
      if (ret_data) ret_dat$aqps <- aqp_tibs

      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Write Mock PQC File::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( vb >= vt+1 )
        cat(glue::glue( "{mssg}{TAB} Writing Single PQC File for ",
                        "Beadpool.{RET}"))

      seqs_cnt <- pqc_tib %>% dplyr::distinct( Address ) %>% base::nrow()
      bead_cnt <- pqc_tib %>% dplyr::distinct( Address ) %>% base::nrow()
      pass_cnt <- bead_cnt
      
      # Define PQC Header::
      pqc_header <- 
        glue::glue( "Software Version{TAB}2.5.12{RET}",
                    "Pool{TAB}Unknown123456789{RET}",
                    "Revision{TAB}A{RET}",
                    "Total Sequences{TAB}{seqs_cnt}{RET}",
                    "Total Bead Types{TAB}{bead_cnt}{RET}",
                    "Total Bead Types Passed{TAB}{pass_cnt}{RET}" )
      
      pqc_tsv <- paste(cur_basename, "final-product.pqc.tsv.gz", sep = '.')
      pqc_tsv <- file.path( aqp_dir, pqc_tsv )
      
      safe_write( x = pqc_header, file = pqc_tsv, type = "line", done = FALSE,
                  write_spec = FALSE, append = FALSE, 
                  vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      safe_write( x = pqc_tib, file = pqc_tsv, type = "tsv", done = TRUE,
                  write_spec = FALSE, append = TRUE,
                  vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Build Mock AQP Files::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( parallel ) {
        if ( vb >= vt+1 )
          cat(glue::glue( "{mssg}{TAB} Building/Writing Mock AQP Files for ",
                          "each Beadpool in parallel processing...{RET2}") )
        
        foreach (beadpool=names( mat_tibs ) ) %dopar% {
          
          if ( vb >= vt+1 )
            cat(glue::glue( "{mssg}{TAB} Building Mock AQP Files for ",
                            "Beadpool = '{beadpool}' in parallel ",
                            "processing...{RET}") )
          
          aqp_tsv <- paste( cur_basename,beadpool,"AQP-1.aqp.tsv.gz", sep = '.' )
          aqp_tsv <- file.path( aqp_dir, aqp_tsv )
          
          if ( vb >= vt+1 )
            cat(glue::glue( "{mssg}{TAB} Current aqp_tsv = {aqp_tsv}.{RET}") )
          
          aqp_tib <- NULL
          aqp_tib <- aqp_tibs[[ beadpool ]] %>% 
            dplyr::distinct( Address ) %>% 
            dplyr::arrange( Address ) %>%
            dplyr::mutate( Dec = 0, 
                           Err = 0, 
                           Scr = 300, 
                           Fun = "",
                           Frr = "",
                           QCA = 0 ) %>% 
            purrr::set_names( aqp_cols )
          
          aqp_key <- glue::glue("cur-aqp-tib (bp = {beadpool})")
          aqp_cnt <- print_tib( aqp_tib, fun_tag = fun_tag, name = aqp_key,
                                vb=vb,vt=vt+5,tc=tc+1,tt=tt )
          
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          #                          Write Mock AQP Files::
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          
          if ( vb >= vt+1 )
            cat(glue::glue( "{mssg}{TAB} Writing Mock AQP Files for ",
                            "Beadpool = '{beadpool}' in parallel ",
                            "processing...{RET2}") )
          
          seqs_cnt <- aqp_tib %>% dplyr::distinct( Address ) %>% base::nrow()
          bead_cnt <- aqp_tib %>% dplyr::distinct( Address ) %>% base::nrow()
          pass_cnt <- bead_cnt
          
          # Define AQP Header::
          aqp_header = 
            glue::glue( "Software Version{TAB}0.0.0{RET}",
                        "Slide{TAB}Multiple Used{RET}",
                        "Pool{TAB}{beadpool}.123456789-AQP{RET}",
                        "Total Sequences{TAB}{seqs_cnt}{RET}",
                        "Total Bead Types{TAB}{bead_cnt}{RET}",
                        "Total Bead Types Passed{TAB}{pass_cnt}{RET}" )
          
          if ( vb >= vt+2 )
            cat(glue::glue( "{mssg}{TAB} aqp_header={aqp_header}.{RET}") )
          
          safe_write( x = aqp_header, file = aqp_tsv, type = "line",
                      done = FALSE, write_spec = FALSE, append = FALSE,
                      vb=vb,vt=vt+4,tc=tc+1,tt=tt )
          
          safe_write( x = aqp_tib, file = aqp_tsv, type = "tsv",
                      done = TRUE, write_spec = FALSE, append = TRUE,
                      vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        }

      } else {
        
        if ( vb >= vt+1 )
          cat(glue::glue( "{mssg}{TAB} Building Mock AQP Files for Beadpool = ",
                          "'{beadpool}' in linear processing...{RET}") )
        
        for ( beadpool in names( aqp_tibs ) ) {
          aqp_tsv <- paste( cur_basename,beadpool,"AQP-1.aqp.tsv.gz", sep = '.' )
          aqp_tsv <- file.path( aqp_dir, aqp_tsv )
          
          if ( vb >= vt+1 )
            cat(glue::glue( "{mssg}{TAB} Current aqp_tsv = {aqp_tsv}.{RET}") )
          
          aqp_tib <- NULL
          aqp_tib <- aqp_tibs[[ beadpool ]] %>% 
            dplyr::distinct( Address ) %>% 
            dplyr::arrange( Address ) %>%
            dplyr::mutate( Dec = 0, 
                           Err = 0, 
                           Scr = 300, 
                           Fun = "",
                           Frr = "",
                           QCA = 0 ) %>% 
            purrr::set_names( aqp_cols )
          
          aqp_key <- glue::glue("cur-aqp-tib (bp = {beadpool})")
          aqp_cnt <- print_tib( aqp_tib, fun_tag = fun_tag, name = aqp_key,
                                vb=vb,vt=vt+5,tc=tc+1,tt=tt )
          
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          #                          Write Mock AQP Files::
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          
          if ( vb >= vt+1 )
            cat(glue::glue( "{mssg}{TAB} Writing AQP Files for Beadpool = ",
                            "'{beadpool}' in linear processing...{RET2}") )
          
          seqs_cnt <- aqp_tib %>% dplyr::distinct( Address ) %>% base::nrow()
          bead_cnt <- aqp_tib %>% dplyr::distinct( Address ) %>% base::nrow()
          pass_cnt <- bead_cnt
          
          # Define AQP Header::
          aqp_header = 
            glue::glue( "Software Version{TAB}0.0.0{RET}",
                        "Slide{TAB}Multiple Used{RET}",
                        "Pool{TAB}{beadpool}.123456789-AQP{RET}",
                        "Total Sequences{TAB}{seqs_cnt}{RET}",
                        "Total Bead Types{TAB}{bead_cnt}{RET}",
                        "Total Bead Types Passed{TAB}{pass_cnt}{RET}" )
          
          if ( vb >= vt+2 )
            cat(glue::glue( "{mssg}{TAB} aqp_header={aqp_header}.{RET}") )
          
          safe_write( x = aqp_header, file = aqp_tsv, type = "line",
                      done = FALSE, write_spec = FALSE, append = FALSE,
                      vb=vb,vt=vt+4,tc=tc+1,tt=tt )
          
          safe_write( x = aqp_tib, file = aqp_tsv, type = "tsv",
                      done = TRUE, write_spec = FALSE, append = TRUE,
                      vb=vb,vt=vt+4,tc=tc+1,tt=tt )
        }

      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                               Done:: Next...
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( !file.exists(end_txt) )
        sys_ret <- base::system( glue::glue("touch {end_txt}") )
      
      if ( vb >= vt+1 )
        cat(glue::glue( "{mssg}{TAB} Done. {man_name} .{RET2}") )
    }

  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  if (ret_data) return(ret_dat)
  
  ret_tib
}


load_genome_studio_address = function(
  file,
  load_clean     = TRUE,
  load_controls  = FALSE,
  write_clean    = TRUE,
  overwrite      = FALSE,
  add_annotation = FALSE,
  ret_data       = FALSE,
  old_cols = c("IlmnID", "Name", 
               "AddressA_ID", "AlleleA_ProbeSeq", 
               "AddressB_ID", "AlleleB_ProbeSeq", "Next_Base",
               "Color_Channel", "CHR", "MAPINFO", "Strand"),
  new_cols = c("Probe_ID", "Prb_Cgn", 
               "Address_A", "Prb_Seq_A", "Address_B", "Prb_Seq_B",
               "Prb_Nxb", "Color", "Chromosome", "Coordinate", "Strand_FR"),
  non_ann_cols = c("Infinium_Design_Type", "Forward_Sequence",
                   "Genome_Build", "SourceSeq"),
  
  vb=0, vt=3, tc=1, tt=NULL,
  fun_tag='load_genome_studio_address') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          file={file}.{RET}"))
    cat(glue::glue("{mssg}    load_clean={load_clean}.{RET}"))
    cat(glue::glue("{mssg}   write_clean={write_clean}.{RET}"))
    cat(glue::glue("{mssg}     overwrite={overwrite}.{RET}"))
    cat(glue::glue("{mssg}      ret_data={ret_data}.{RET}"))
    cat(glue::glue("{RET}"))
    if (vb >= vt+4) {
      old_cnt <- old_cols %>% length()
      new_cnt <- new_cols %>% length()
      cat(glue::glue("{mssg} Old Genome Studio Columns({old_cnt})::{RET}"))
      cat(glue::glue("{mssg}    old_cols={old_cols}.{RET}"))
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} New Genome Studio Columns({new_cnt})::{RET}"))
      cat(glue::glue("{mssg}    new_cols={new_cols}.{RET}"))
      cat(glue::glue("{RET}"))
    }
  }  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    man_tib <-
      load_genome_studio_manifest(file = file,
                                  load_clean    = load_clean,
                                  load_controls = load_controls,
                                  write_clean   = write_clean,
                                  overwrite     = overwrite,
                                  ret_data      = ret_data,
                                  vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    man_cnt <- man_tib %>% base::nrow()
    if (vb >= vt+1)
      cat(glue::glue("{mssg} Loaded Genome Studio manifest({man_cnt})!{RET2}"))
    
    if (is.null(man_tib) || !tibble::is_tibble(man_tib)) {
      stop(glue::glue("{RET}{mssg} ERROR: Failed to load Genome Studio ",
                      "manifest!{RET2}"))
      return(ret_tib)
    }
    man_key <- glue::glue("genome-studio-manifest-1({fun_tag})")
    man_cnt <- print_tib(man_tib, fun_tag = fun_tag, name = man_key,
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Swap Old/New Column Names::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    man_tib <- dplyr::rename_with(man_tib, ~ new_cols, dplyr::all_of(old_cols) )
    
    man_key <- glue::glue("genome-studio-manifest-renamed-2({fun_tag})")
    man_cnt <- print_tib(man_tib, fun_tag = fun_tag, name = man_key,
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    if (!add_annotation) {
      man_tib <- man_tib %>% 
        dplyr::select(dplyr::any_of( c(new_cols,non_ann_cols) ))
    }
    
    man_key <- glue::glue("genome-studio-manifest-annotation-3({fun_tag})")
    man_cnt <- print_tib(man_tib, fun_tag = fun_tag, name = man_key,
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract/Format Probe A::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_key <- "A"
    old_prbA_cols <- man_tib %>% 
      dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
      names()
    new_prbA_cols <- old_prbA_cols %>%
      stringr::str_remove(paste0("_",prb_key,"$"))
    
    oldA_cnt <- old_prbA_cols %>% length()
    newA_cnt <- new_prbA_cols %>% length()
    if (vb >= vt+2)
      cat(glue::glue("{mssg} Old Sesame Columns:: Probe A({oldA_cnt}){RET}"))
    if (vb >= vt+4)
      cat(glue::glue("{mssg}    old_prbA_cols={old_prbA_cols}.{RET2}"))
    if (vb >= vt+2)
      cat(glue::glue("{mssg} New Sesame Columns:: Probe A({newA_cnt}){RET}"))
    if (vb >= vt+4)
      cat(glue::glue("{mssg}    new_prbA_cols={new_prbA_cols}.{RET2}"))
    
    prb_key <- "B"
    addA_tib <- NULL
    addA_tib <- man_tib %>% 
      dplyr::select(-dplyr::ends_with(paste0("_",prb_key))) %>% 
      dplyr::rename_with( ~ new_prbA_cols, dplyr::all_of(old_prbA_cols) ) %>%
      dplyr::mutate(Prb_Des=dplyr::case_when(
        Infinium_Design_Type=="I"  ~ "U",
        Infinium_Design_Type=="II" ~ "2",
        TRUE ~ NA_character_)
      )
    
    addA_cnt <- addA_tib %>% base::nrow()
    if (vb >= vt+2)
      cat(glue::glue("{mssg} Probe Set A: I/U Count={addA_cnt}.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract/Format Probe B::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_key <- "B"
    old_prbB_cols <- man_tib %>% 
      dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
      names()
    new_prbB_cols <- old_prbB_cols %>%
      stringr::str_remove(paste0("_",prb_key,"$"))
    
    oldB_cnt <- old_prbB_cols %>% length()
    newB_cnt <- new_prbB_cols %>% length()
    if (vb >= vt+2)
      cat(glue::glue("{mssg} Old Sesame Columns:: Probe B({oldB_cnt}){RET}"))
    if (vb >= vt+4)
      cat(glue::glue("{mssg}    old_prbB_cols={old_prbB_cols}.{RET2}"))
    if (vb >= vt+2)
      cat(glue::glue("{mssg} New Sesame Columns:: Probe B({newB_cnt}){RET}"))
    if (vb >= vt+4)
      cat(glue::glue("{mssg}    new_prbB_cols={new_prbB_cols}.{RET2}"))
    
    prb_key <- "A"
    addB_tib <- NULL
    addB_tib <- man_tib %>% 
      dplyr::filter(Infinium_Design_Type=="I") %>%
      dplyr::select(-dplyr::ends_with(paste0("_",prb_key))) %>% 
      dplyr::rename_with( ~ new_prbB_cols, dplyr::all_of(old_prbB_cols) ) %>% 
      dplyr::mutate(Prb_Des="M")
    
    addB_cnt <- addB_tib %>% base::nrow()
    if (vb >= vt+2)
      cat(glue::glue("{mssg} Probe Set B: M Count={addB_cnt}.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Stack Probe Designs:: A/B
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- dplyr::bind_rows(addA_tib,addB_tib) %>%
      dplyr::mutate(
        Prb_Inf=dplyr::case_when(
          Infinium_Design_Type=="I"  ~ 1,
          Infinium_Design_Type=="II" ~ 2,
          TRUE ~ NA_real_
        ) %>% as.integer(),
        Prb_Din = Probe_ID %>%
          stringr::str_sub(1,2),
        Prb_Cgn = Prb_Cgn %>%
          stringr::str_remove("^ch.[0-9A-Z]+.") %>%
          stringr::str_remove("[A-Z]+$") %>%
          stringr::str_remove("^[^0-9]+") %>%
          stringr::str_remove("^0+") %>%
          as.integer(),
        Chromosome=Chromosome %>%
          stringr::str_remove("^chr") %>%
          paste0("chr",.)
      ) %>% dplyr::arrange(Probe_ID)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Summary and Sanity Checks::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # This better be zero::
    dup_cnt <- ret_tib %>% 
      dplyr::arrange(Address) %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (dup_cnt != 0) {
      stop(glue::glue("{mssg} ERROR: Duplicates Found = {dup_cnt}!{RET2}"))
      return(NULL)
    }
    if (vb >= vt+1)
      cat(glue::glue("{mssg} No duplicates detected = {dup_cnt}.{RET}"))
    
    if (vb >= vt+2) {
      cat(glue::glue("{mssg} Probe Design Coverage={RET}"))
      sum_tib <- ret_tib %>% 
        dplyr::group_by(Prb_Des) %>% 
        dplyr::summarise(Count=n(), .groups = "drop")
      sum_tib %>% print(n=base::nrow(sum_tib))
      cat(glue::glue("{RET}"))
    }
    
    ret_key <- glue::glue("ret-FIN({fun_tag})")
    ret_cnt <- print_tib(ret_tib, fun_tag = fun_tag, name = ret_key,
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

load_genome_studio_manifest = function(file,
                                       
                                       load_clean    = TRUE,
                                       load_controls = FALSE,
                                       cols_convert  = FALSE,
                                       write_clean   = TRUE,
                                       overwrite     = FALSE,
                                       ret_data      = FALSE,
                                       
                                       vb=0, vt=3, tc=1, tt=NULL,
                                       fun_tag='load_genome_studio_manifest') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+2) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          file={file}.{RET}"))
    cat(glue::glue("{mssg}    load_clean={load_clean}.{RET}"))
    cat(glue::glue("{mssg}   write_clean={write_clean}.{RET}"))
    cat(glue::glue("{mssg}     overwrite={overwrite}.{RET}"))
    cat(glue::glue("{mssg}      ret_data={ret_data}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  probes_tib   <- NULL
  controls_tib <- NULL
  
  stime <- base::system.time({
    
    # Set Column Repair Prefix::
    col_prefix <- "Val_"
    if (stringr::str_detect(file,"EPIC")) col_prefix <- "HM"
    
    gs_control_names <- 
      c("Address","Control_Group","GS_Color","Control_Type")
    
    analytical_suffix_csv <- ".analytical.csv.gz"
    analytical_suffix_rds <- ".analytical.rds"
    control_suffix_csv    <- ".controls.csv.gz"
    
    # if (!R.utils::isGzipped(file, ".gz"))
    #   file <- R.utils::gzip(file) %>% as.character()
    if (!stringr::str_ends(file,".gz")) {
      system(glue::glue("gzip {file}"))
      file <- paste0(file,".gz")
    }
    
    if (vb >= vt+3) {
      cat(glue::glue("{mssg} Setting suffix variables::{RET}"))
      cat(glue::glue("{mssg}       control_suffix_csv=",
                     "{control_suffix_csv}.{RET}"))
      cat(glue::glue("{mssg}    analytical_suffix_csv=",
                     "{analytical_suffix_csv}.{RET}"))
      cat(glue::glue("{mssg}    analytical_suffix_rds=",
                     "{analytical_suffix_rds}.{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    # If the input file is the analytical file load it!
    if (stringr::str_ends(file, analytical_suffix_csv)) {
      clean_prb_csv <- file
      clean_col_rds <- clean_basename %>%
        paste(analytical_suffix_csv)
      clean_prb_done <- paste(clean_prb_csv,'done.txt', sep='.')
      
      if (vb >= vt+2)
        cat(glue::glue("{mssg} Found clean manifest! Will load clean ",
                       "manifest={clean_prb_csv}...{RET}"))
      
      if (file.exists(clean_col_rds) &&
          file.mtime(clean_prb_csv) <= file.mtime(clean_col_rds)) {
        if (vb >= vt+2)
          cat(glue::glue("{mssg} Using column header={clean_col_rds}...{RET}"))
        
        probes_type_cols <- readr::read_rds(clean_col_rds)
        probes_tib  <- readr::read_csv(clean_prb_csv,
                                       # col_names = names(probes_type_cols$cols),
                                       col_types = probes_type_cols)
      } else {
        if (vb >= vt+2)
          cat(glue::glue("{mssg} Will build column header...{RET}"))
        
        probes_tib <- safe_read( clean_prb_csv, clean = TRUE, 
                                 vb=vb,vt=vt+1,tc=tc+1,tt=tt )
        
        # Get Column Types String::
        col_types_str <- sapply(probes_tib, typeof) %>% 
          tibble::as_tibble(rownames = "Key") %>% 
          dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
          dplyr::pull(value) %>% paste0(collapse = '')
        
        # Build cols() object for future loading::
        probes_type_cols <- spec(readr::read_csv(
          readr::format_csv(probes_tib), col_types = col_types_str) )
        
        if (write_clean) {
          if (vb >= vt+2) 
            cat(glue::glue("{mssg} Writing clean manifest cols RDS...{RET}"))
          
          safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                      vb=vb,vt=vt+1,tc=tc+1,tt=tt )
        }
      }
      
      if (load_controls) {
        controls_csv <- file %>% 
          stringr::str_replace(analytical_suffix_csv, control_suffix_csv)
        
        if (!file.exists(controls_csv)) {
          cat(glue::glue("{mssg} Warning: No controls file present. Try ",
                         "rebuilding genome studio manifest from scratch.{RET}"))
        } else {
          if (vb >= vt+2)
            cat(glue::glue("{mssg} Loading controls={controls_csv}...{RET}"))
          
          controls_tib <- 
            suppressMessages(suppressWarnings( readr::read_csv( controls_csv )))
        }
      }
      
    } else {
      
      # Define clean formatted files::
      clean_basename <- file %>% 
        stringr::str_remove(".gz$") %>%
        stringr::str_remove(".csv$")
      
      clean_prb_csv <- clean_basename %>%
        paste0(analytical_suffix_csv)
      clean_col_rds <- clean_basename %>%
        paste0(analytical_suffix_rds)
      clean_prb_done <- paste(clean_prb_csv,'done.txt', sep='.')
      controls_csv   <- clean_basename %>%
        paste0(control_suffix_csv)
      
      if (vb >= vt+3) {
        cat(glue::glue("{mssg} Setting output variables::{RET}"))
        cat(glue::glue("{mssg}    clean_basename={clean_basename}.{RET}"))
        cat(glue::glue("{mssg}      controls_csv={controls_csv}.{RET}"))
        cat(glue::glue("{mssg}     clean_col_rds={clean_col_rds}.{RET}"))
        cat(glue::glue("{mssg}     clean_prb_csv={clean_prb_csv}.{RET}"))
        cat(glue::glue("{mssg}    clean_prb_done={clean_prb_done}.{RET}"))
        cat(glue::glue("{RET}"))
      }
      
      if (load_clean && file.exists(clean_prb_csv) &&
          file.mtime(file) <= file.mtime(clean_prb_csv) &&
          file.mtime(clean_prb_csv) <= file.mtime(clean_prb_done) ) {
        
        if (vb >= vt+2)
          cat(glue::glue("{mssg} Loading clean manifest={clean_prb_csv}...{RET}"))
        
        if (file.exists(clean_col_rds) &&
            file.mtime(clean_prb_csv) <= file.mtime(clean_col_rds)) {
          if (vb >= vt+2)
            cat(glue::glue("{mssg} Using column header={clean_col_rds}...{RET}"))
          
          probes_type_cols <- readr::read_rds(clean_col_rds)
          
          probes_tib  <- readr::read_csv(clean_prb_csv,
                                         # col_names = names(probes_type_cols$cols),
                                         col_types = probes_type_cols)
        } else {
          if (vb >= vt+2) 
            cat(glue::glue("{mssg} Will build column header...{RET}"))
          
          probes_tib <- safe_read( clean_prb_csv, clean = TRUE, 
                                   vb=vb,vt=vt+1,tc=tc+1,tt=tt )
          
          # Get Column Types String::
          col_types_str <- sapply(probes_tib, typeof) %>% 
            tibble::as_tibble(rownames = "Key") %>% 
            dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
            dplyr::pull(value) %>% paste0(collapse = '')
          
          # Build cols() object for future loading::
          probes_type_cols <- spec(readr::read_csv(
            readr::format_csv(probes_tib), col_types = col_types_str) )
          
          if (write_clean) {
            if (vb >= vt+2) 
              cat(glue::glue("{mssg} Writing clean manifest cols RDS...{RET}"))
            
            safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                        vb=vb,vt=vt+1,tc=tc+1,tt=tt )
          }
        }
        
        # Load Controls if requested
        if (load_controls) {
          if (!file.exists(controls_csv))
            cat(glue::glue("{mssg} Warning: No controls file present. Try ",
                           "rebuilding genome studio manifest from scratch.{RET}"))
          else
            controls_tib <- 
              suppressMessages(suppressWarnings( readr::read_csv( controls_csv )))
        }
        
      } else {
        if (vb >= vt+2)
          cat(glue::glue("{mssg} Will load raw manifest={file}...{RET}"))
        
        if (vb >= vt+2)
          cat(glue::glue("{mssg} Will build column header...{RET}"))
        
        lines_vec <- readr::read_lines(file)
        if (ret_data) ret_dat$lines <- lines_vec
        
        # Get Analytical and Control Start & End Indexes::
        beg_idx <- head(
          which( stringr::str_starts(head(lines_vec, n=20), "IlmnID" ) ), n=1)
        con_idx <- tail(
          which( stringr::str_starts(lines_vec, "\\[Controls\\]" ) ), n=1)
        end_idx <- lines_vec %>% length()
        
        if (vb >= vt+2)
          cat(glue::glue("{mssg} Setting Indexes=({beg_idx}, {con_idx}, ",
                         "{end_idx}).{RET}"))
        
        # Get Analytical Column Names::
        probes_name_cols <- 
          lines_vec[beg_idx] %>% 
          stringr::str_split(pattern = ",", simplify = TRUE) %>% 
          as.vector() %>% 
          stringr::str_replace("^([0-9])", paste0(col_prefix, "\\$1") ) %>%
          stringr::str_remove_all("\\\\")
        
        probes_name_cnt <- probes_name_cols %>% length()
        if (vb >= vt+3) {
          name1_str <- probes_name_cols[1]
          name2_str <- probes_name_cols[probes_name_cnt]
          out_str <- glue::glue("({name1_str} ... {name2_str})")
          cat(glue::glue("{mssg} Extracted Probe Names Vec({probes_name_cnt})",
                         " = {out_str}.{RET2}"))
        }
        
        probes_name_cols <- stringi::stri_remove_empty( probes_name_cols )
        
        # Build Analytical Probes Tibble
        probes_tib <- 
          lines_vec[c( (beg_idx+1):(con_idx-1) ) ] %>%
          stringr::str_remove(",+$") %>%
          tibble::as_tibble( .name_repair = "unique" ) %>%
          tidyr::separate(value, 
                          into=c( probes_name_cols ), 
                          sep=',', extra = "drop", fill = "right" ) %>%
          dplyr::mutate_all( list(~na_if(.,"") ) ) %>%
          clean_tib()

        ret_key <- glue::glue("probes-tib({fun_tag})")
        ret_cnt <- print_tib(probes_tib, fun_tag = fun_tag, name = ret_key,
                             vb=vb,vt=vt+1,tc=tc+1,tt=tt ) 
        
        # Get Column Types String::
        col_types_str <- sapply(probes_tib, typeof) %>% 
          tibble::as_tibble(rownames = "Key") %>% 
          dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
          dplyr::pull(value) %>% paste0(collapse = '')
        
        # Build cols() object for future loading::
        probes_type_cols <- spec(readr::read_csv(
          readr::format_csv(probes_tib), col_types = col_types_str) )
        if (vb >= vt+4) {
          cat(glue::glue("{mssg} probes_type_cols={RET}"))
          probes_type_cols %>% print()
        }
        
        # Build Controls Tibble::
        controls_tib <- 
          lines_vec[c( (con_idx+1):(end_idx) ) ] %>%
          stringr::str_remove(",+$") %>%
          stringr::str_remove(",AVG$") %>%
          tibble::as_tibble() %>% 
          tidyr::separate(value, into=c(gs_control_names), sep=',') %>% 
          clean_tib()
        ret_key <- glue::glue("controls-tib({fun_tag})")
        ret_cnt <- print_tib(controls_tib, fun_tag = fun_tag, name = ret_key,
                             vb=vb,vt=vt+1,tc=tc+1,tt=tt )
        
        if (ret_data) ret_dat$probes   <- probes_tib
        if (ret_data) ret_dat$controls <- controls_tib
        
        if (ret_data) ret_dat$probes_key_str <- probes_name_cols
        if (ret_data) ret_dat$probes_col_str <- col_types_str
        if (ret_data) ret_dat$probes_cols    <- probes_type_cols
        
        # Write clean files if requested::
        if (write_clean) {
          if (vb >= vt+2) 
            cat(glue::glue("{mssg} Writing clean manifest files...{RET}"))
          
          out_tag <- glue::glue("write-controls-CSV({fun_tag})")
          safe_write( controls_tib, file = controls_csv, done = TRUE, 
                      fun_tag = out_tag, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
          
          out_tag <- glue::glue("write-probes-CSV({fun_tag})")
          safe_write( probes_tib, file = clean_prb_csv, done = TRUE,
                      fun_tag = out_tag, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
          
          out_tag <- glue::glue("write-probes-col-RDS({fun_tag})")
          safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                      fun_tag = out_tag, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
        }
      }
    }
    
    probes_tib <- probes_tib %>% dplyr::arrange(IlmnID)
    
    if ( cols_convert ) {
      if ( vb >= vt+1 )
        cat(glue::glue("{mssg} Converting to modern names...{RET}"))
      
      probes_tib <- probes_tib %>%
        dplyr::select( -Name, -AddressA_ID, -AddressB_ID ) %>%
        dplyr::rename(
          Loci_ID = IlmnID,
          Probe_Seq_A = AlleleA_ProbeSeq,
          Probe_Seq_B = AlleleB_ProbeSeq,
          Strand_FR = Strand,
          Chromosome = CHR,
          Coordinate = MAPINFO
        ) %>%
        dplyr::mutate( 
          Infinium_Design = dplyr::case_when(
            !is.na(Probe_Seq_A) & !is.na(Probe_Seq_B) ~ 1,
            !is.na(Probe_Seq_A) &  is.na(Probe_Seq_B) ~ 2,
            TRUE ~ NA_real_ ) %>% as.integer(),
          Probe_Type = Loci_ID %>% stringr::str_sub( 1,2 ),
          Probe_Seq_A = Probe_Seq_A %>% stringr::str_to_upper(),
          Probe_Seq_B = Probe_Seq_B %>% stringr::str_to_upper(),
          Align_Seq_A = mutate_seqs_cpp( seqs = Probe_Seq_A, m='X', uc = TRUE ),
          Align_Seq_B = mutate_seqs_cpp( seqs = Probe_Seq_B, m='X', uc = TRUE ),
          Probe_48M_U = dplyr::case_when(
            Infinium_Design == 1 ~ stringr::str_sub( Align_Seq_A, 2,49 ),
            Infinium_Design == 2 ~ stringr::str_sub( Align_Seq_A, 3,50 ),
            TRUE ~ NA_character_ )
        )
    }
    
    if (load_controls && is.null(controls_tib)) {
      ret_dat$probes   <- probes_tib
      ret_dat$controls <- controls_tib
      ret_data <- TRUE
    }
    
    ret_key <- glue::glue("ret-FIN({fun_tag})")
    ret_cnt <- print_tib(probes_tib, fun_tag = fun_tag, name = ret_key,
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if (ret_data) return(ret_dat)
  
  probes_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                             Noob Probe_ID Masking::
# noob_mask()
# noob_mask_manifest()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Should probably clean these notes and both funcions up...
#
# changing default mod=100000000 to mod=999999999
# changing default prefix="nb" to prefix="cg9999"
# NOTE: 5 9's is optimal because improbe can only handle 15 characters
#  as input 2 character (cg) + 8 digits + 2 character ([TB][CO]) = 12
# However, we can reduce characters to a single leter 1+8 = 9 
#  Using 5 would put us at 14 leaving one charcter to expand 8 digits to 9
#
# None of that really matters since we'll never use these probes in improbe
#  The issue is running into your growing cg# (8 to 9 to 10) space. 
#  We'll make 5 9's
# Wait why don't we just make cgX
#
noob_mask = function(x, seed=21L, mod=100000000, prefix="cgBK",
                     fun_tag='noob_mask') {
  
  # Clean input::
  x <- x %>% stringr:: str_remove_all("[^0-9]+") %>% as.integer()
  
  hash <- digest::digest2int(as.character(x), seed) %% mod
  
  # Waiting to see if this ever fails
  m_len <- stringr::str_length(format(mod, scientific = FALSE))
  h_len <- stringr::str_length(hash)
  
  if (h_len>=m_len) {
    stop(glue::glue("{RET}[{fun_tag}]: ERROR: ",
                    "h_len({h_len}) >= m_len({m_len})!!!{RET}",
                    "hash={hash}{RET}",
                    "mod={mod}{RET}{RET}"))
    return(NULL)
  }
  
  # Now remake cg style number::
  hash <- paste0(prefix,
                 stringr::str_pad(hash, width=m_len-1,side="left", pad="0"))
  
  hash
}

noob_mask_manifest = function(tib,
                              key = "Probe_ID", 
                              out = NULL, 
                              prefix = "cg",
                              
                              vb=0, vt=3, tc=1, tt=NULL,
                              fun_tag='noob_mask_manifest') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      key={key}.{RET}"))
    cat(glue::glue("{mssg}      out={out}.{RET}"))
    cat(glue::glue("{mssg}   prefix={prefix}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (is.null(out)) out <- key
    key_sym <- rlang::sym(key)
    out_sym <- rlang::sym(out)
    
    noob_vec <- tib %>% 
      dplyr::pull(key_sym) %>% 
      lapply( noob_mask) %>% unlist()
    
    ret_tib <- tib %>% dplyr::mutate(!!out_sym := noob_vec)
    
    ret_key <- glue::glue("ret-FIN({fun_tag})")
    ret_cnt <- print_tib(ret_tib, fun_tag = fun_tag, name = ret_key,
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file

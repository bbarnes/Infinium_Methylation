
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               TripleCrown::
#                         CG# Database Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Command Line Options Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )

# Tidyverse Core Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("stringr", quietly = TRUE) ))
suppressWarnings(suppressPackageStartupMessages( 
  base::require("readr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("glue",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("furrr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("purrr",    quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("magrittr", quietly = TRUE) ) )

# Parallel Processing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

# Plotting Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("ggplot2", quietly = TRUE) ) )

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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Probe Performance Analysis::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

replicate_analysis_r = function( mat,
                                 ssh,
                                 par,
                                 
                                 out_dir = NULL,
                                 
                                 sds_min = 0.1,
                                 mad_min = 0.1,
                                 
                                 screen_key = "rep",
                                 percision = 0,
                                 
                                 vb=0, vt=6, tc=1, tt=NULL,
                                 fun_tag='replicate_analysis_r')
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
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  # ssh_list <- NULL
  # ssh_list <- ssh %>% split(.$Concentration)
  
  stats_tib <- NULL
  stats_tib <- tibble::tibble(
    Probe_ID = mat %>% rownames(),
    Tot = mat %>% ncol(),
    Nan = mat %>% matrixStats::rowCounts( value = NA ),
    PdP = round( 100*(Tot-Nan)/Tot, 2 ),
    # Avg = mat %>% matrixStats::rowMeans2( na.rm = TRUE ),
    Sds = mat %>% matrixStats::rowSds( na.rm = TRUE ),
    # Med = mat %>% matrixStats::rowMedians( na.rm = TRUE ),
    Mad = mat %>% matrixStats::rowMads( na.rm = TRUE )
    # Var = mat %>% matrixStats::rowVars( na.rm = TRUE ),
  ) %>% 
    dplyr::filter( !Probe_ID %>% stringr::str_starts("^ctl_") ) %>%
    #
    # TBD:: Add !is.na()
    #
    dplyr::mutate( Sds_Pass = !is.na(Sds) & Sds <= sds_min,
                   Mad_Pass = !is.na(Mad) & Mad <= mad_min,
                   mask = !(Sds_Pass & Mad_Pass) )
  
  # print( stats_tib )
  
  # Build UHM Params::
  pars_tib <- NULL
  pars_tib <- par %>% dplyr::mutate( Screen = screen_key )
  pars_str <- pars_tib %>% as.vector() %>% paste( collapse = "." )
  
  # print( pars_tib )
  
  # Join UHM Params and Probe Results::
  mask_tib <- NULL
  mask_tib <- pars_tib %>% dplyr::bind_cols( stats_tib )
  
  # print( mask_tib )
  
  # Concalculate UHM Summary::
  mask_sum <- NULL
  mask_sum <- mask_tib %>% 
    dplyr::group_by( dplyr::all_of( pars_tib ) ) %>% 
    dplyr::summarise( Tot = n(),
                      Nan = sum( is.na(mask) ),
                      Mis = sum(  mask, na.rm=TRUE ),
                      Pas = sum( !mask, na.rm=TRUE ),
                      Per = round( 100*Pas/Tot, 3 ),
                      .groups = "drop" )
  
  # Concatenate UHM Summary::
  # all_uhm_sum <- all_uhm_sum %>% dplyr::bind_rows( mask_sum )
  ret_tib <- mask_sum
  
  # Write Probe UHM Summary::
  if ( !is.null(out_dir) ) {
    out_dir <- safe_mkdir( out_dir, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    mask_csv <- NULL
    mask_csv <- file.path( out_dir, paste( pars_str,"r2.mask.dat.csv.gz", sep='.') )
    
    if ( percision != 0 ) mask_tib <- mask_tib %>%
      dplyr::mutate( dplyr::across( where(is.double), round, digit=percision ) )
    
    readr::write_csv( x = mask_tib, file = mask_csv )
  }
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

titration_analysis_r = function( mat,
                                 ssh,
                                 par,
                                 
                                 out_dir = NULL,
                                 
                                 avgHU_min = 0.2,
                                 avgMH_min = 0.1,
                                 avgMU_min = 0.5,
                                 ord_mis_cnt = 3,
                                 
                                 screen_key = "uhm",
                                 percision = 0,
                                 
                                 vb=0, vt=6, tc=1, tt=NULL,
                                 fun_tag='titration_analysis_r')
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
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  ssh_list <- NULL
  ssh_list <- ssh %>% split(.$Concentration)
  
  stats_tibs <- NULL
  stats_tibs <- ssh_list %>% lapply( function(x) {
    beta_mat <- NULL
    beta_mat <- mat[ , x$Sentrix_Name ] %>% as.matrix()
    
    tibble::tibble(
      Probe_ID = beta_mat %>% rownames(),
      Tot = beta_mat %>% ncol(),
      Nan = beta_mat %>% matrixStats::rowCounts( value = NA ),
      PdP = round( 100*(Tot-Nan)/Tot, 2 ),
      Avg = beta_mat %>% matrixStats::rowMeans2( na.rm = TRUE ),
      # Sds = beta_mat %>% matrixStats::rowSds( na.rm = TRUE ),
      Med = beta_mat %>% matrixStats::rowMedians( na.rm = TRUE )
      # Mad = beta_mat %>% matrixStats::rowMads( na.rm = TRUE ),
      # Var = beta_mat %>% matrixStats::rowVars( na.rm = TRUE ),
    ) %>% dplyr::filter( !Probe_ID %>% stringr::str_starts("^ctl_") )
    
  })
  
  # Compare intensities in order of titration::
  ord_tib <- NULL
  ord_tib <- tibble::tibble(
    Probe_ID = stats_tibs[[1]]$Probe_ID,
    
    dAvg12 = stats_tibs[[2]]$Avg - stats_tibs[[1]]$Avg,
    dMed12 = stats_tibs[[2]]$Med - stats_tibs[[1]]$Med,
    
    dAvg23 = stats_tibs[[3]]$Avg - stats_tibs[[2]]$Avg,
    dMed23 = stats_tibs[[3]]$Med - stats_tibs[[2]]$Med,
    
    # dAvg34 = stats_tibs[[4]]$Avg - stats_tibs[[3]]$Avg,
    # dMed34 = stats_tibs[[4]]$Med - stats_tibs[[3]]$Med,
    # 
    # dAvg45 = stats_tibs[[5]]$Avg - stats_tibs[[4]]$Avg,
    # dMed45 = stats_tibs[[5]]$Med - stats_tibs[[4]]$Med,
    # 
    # dAvg56 = stats_tibs[[6]]$Avg - stats_tibs[[5]]$Avg,
    # dMed56 = stats_tibs[[6]]$Med - stats_tibs[[5]]$Med,
    # 
    # dAvg67 = stats_tibs[[7]]$Avg - stats_tibs[[6]]$Avg,
    # dMed67 = stats_tibs[[7]]$Med - stats_tibs[[6]]$Med,
  )
  
  # Use Matrix R Math to do some quick calculations::
  ord_mat <- NULL
  ord_mat <- ord_tib %>% 
    tibble::column_to_rownames( var = "Probe_ID" ) %>% as.matrix()
  ord_mat[ which( ord_mat < 0 ) ] <- NA_real_
  
  # Tally Results::
  cut_tib <- NULL
  cut_tib <- tibble::tibble(
    Probe_ID = ord_mat %>% rownames(),
    Tot = ord_mat %>% ncol(),
    Nan = ord_mat %>% matrixStats::rowCounts( value = NA ),
    Ord_Pas_Perc = round( 100*(Tot-Nan)/Tot, 2 ),
    Ord_Pass = Ord_Pas_Perc >= 100*(Tot-ord_mis_cnt)/Tot,
    
    dHU = pmax( (stats_tibs[["50"]]$Avg  - stats_tibs[["0"]]$Avg),
                (stats_tibs[["50"]]$Med  - stats_tibs[["0"]]$Med) ),
    dMH = pmax( (stats_tibs[["100"]]$Avg - stats_tibs[["50"]]$Avg),
                (stats_tibs[["100"]]$Med - stats_tibs[["50"]]$Med) ),
    dMU = pmax( (stats_tibs[["100"]]$Avg - stats_tibs[["0"]]$Avg),
                (stats_tibs[["100"]]$Med - stats_tibs[["0"]]$Med) ),
    
    mask = !(Ord_Pass & 
               dHU >= avgHU_min & 
               dMH >= avgMH_min & 
               dMU >= avgMU_min )
  )
  
  # Build UHM Params::
  pars_tib <- NULL
  pars_tib <- par %>% dplyr::mutate( Screen = screen_key )
  pars_str <- pars_tib %>% as.vector() %>% paste( collapse = "." )
  
  # Join UHM Params and Probe Results::
  mask_tib <- NULL
  mask_tib <- pars_tib %>% dplyr::bind_cols( cut_tib )
  
  # Concalculate UHM Summary::
  mask_sum <- NULL
  mask_sum <- mask_tib %>% 
    dplyr::group_by( dplyr::all_of( pars_tib ) ) %>% 
    dplyr::summarise( Tot = n(),
                      Nan = sum( is.na(mask) ),
                      Mis = sum(  mask, na.rm=TRUE ),
                      Pas = sum( !mask, na.rm=TRUE ),
                      Per = round( 100*Pas/Tot, 3 ),
                      .groups = "drop" )
  
  # Concatenate UHM Summary::
  # all_uhm_sum <- all_uhm_sum %>% dplyr::bind_rows( mask_sum )
  ret_tib <- mask_sum
  
  # Write Probe UHM Summary::
  if ( !is.null(out_dir) ) {
    mask_csv <- NULL
    mask_csv <- file.path( out_dir, paste( pars_str,"r2.mask.dat.csv.gz", sep='.') )
    
    if ( percision != 0 ) mask_tib <- mask_tib %>%
      dplyr::mutate( dplyr::across( where(is.double), round, digit=percision ) )
    
    readr::write_csv( x = mask_tib, file = mask_csv )
  }
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  Specialized Sample Sheet Preparation::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

build_epicv2_sample_sheet = function( sheet_str = NULL,
                                      sheet_vec = NULL,
                                      idat_path = NULL,
                                      
                                      conc = 500,
                                      sentrix_unq = TRUE,
                                      
                                      out_dir,
                                      run_tag,
                                      
                                      reload     = 0,
                                      reload_min = 2,
                                      reload_pre = NULL,
                                      
                                      ret_data   = FALSE,
                                      parallel   = FALSE,
                                      
                                      vb=0, vt=3, tc=1, tt=NULL,
                                      fun_tag='build_epicv2_sample_sheet')
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
  p4  <- vb > vt + 4
  
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
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}    sheet_str = '{sheet_str}'.{RET}"))
    cat(glue::glue("{mssg}    sheet_vec = '{sheet_vec}'.{RET}"))
    cat(glue::glue("{mssg}    idat_path = '{idat_path}'.{RET}"))
    cat(glue::glue("{mssg}  sentrix_unq = '{sentrix_unq}'.{RET}"))
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
  
  errs_mssg <- glue::glue("Outdir out='{out_dir}' does not exist")
  if ( !dir.exists( out_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              Loop Through Each Parameter:: Sample_Sheet
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Both Sheet String and Sheet Vector can't be null")
  if ( ( is.null(sheet_str) || stringr::str_length(sheet_str) == 0 ) &&
       ( is.null(sheet_vec) || length(sheet_vec) == 0 ) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return( NULL )
  
  if ( is.null(sheet_vec) ) {
    sheet_vec <- NULL
    sheet_vec <- stringr::str_split( 
      string = sheet_str, 
      pattern = ",", simplify = FALSE ) %>% 
      base::unlist() %>% as.vector()
  }
  sheet_cnt = sheet_vec %>% length()
  
  errs_mssg <- glue::glue("Sheet Vector is Emplty (size={sheet_cnt})")
  if ( length(sheet_vec) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return( NULL )
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}",
                           "{mssg} \tSheet Vec Size = {sheet_cnt}.{RET}") )
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  ret_tib <- NULL
  ret_tib <- file_list( 
    file = sheet_vec, 
    prefix = "SampleSheet-", 
    suffix=".csv.gz$" ) %>%
    lapply( readr::read_csv, skip = 7, show_col_types = FALSE ) %>%
    dplyr::bind_rows( .id = "File_Sheet_Key" ) %>% 
    dplyr::mutate(
      Sheet_Prep = dplyr::case_when(
        File_Sheet_Key == "LightningAuto"   ~ "Lightning",
        File_Sheet_Key == "LightningManual" ~ "Lightning",
        
        File_Sheet_Key == "LightningAuto_EPIC-only"   ~ "Lightning",
        File_Sheet_Key == "LightningManual_EPIC-only" ~ "Lightning",
        
        File_Sheet_Key == "LightningAuto_EPICv2-only"   ~ "Lightning",
        File_Sheet_Key == "LightningManual_EPICv2-only" ~ "Lightning",
        
        File_Sheet_Key == "EZ-DNA-Manual"             ~ "EZ",
        File_Sheet_Key == "EZ-DNA-Manual-EPIC-only"   ~ "EZ",
        File_Sheet_Key == "EZ-DNA-Manual_EPICv2-only" ~ "EZ",
        TRUE ~ NA_character_
      ),
      Sheet_Proc = dplyr::case_when(
        File_Sheet_Key == "LightningAuto"   ~ "Auto",
        File_Sheet_Key == "LightningManual" ~ "Manual",
        
        File_Sheet_Key == "LightningAuto_EPIC-only"   ~ "Auto",
        File_Sheet_Key == "LightningManual_EPIC-only" ~ "Manual",
        
        File_Sheet_Key == "LightningAuto_EPICv2-only"   ~ "Auto",
        File_Sheet_Key == "LightningManual_EPICv2-only" ~ "Manual",
        
        File_Sheet_Key == "EZ-DNA-Manual"             ~ "Manual",
        File_Sheet_Key == "EZ-DNA-Manual-EPIC-only"   ~ "Manual",
        File_Sheet_Key == "EZ-DNA-Manual_EPICv2-only" ~ "Manual",
        TRUE ~ NA_character_
      ),
      Sheet_Data = dplyr::case_when(
        File_Sheet_Key == "LightningAuto"   ~ "EPICv12",
        File_Sheet_Key == "LightningManual" ~ "EPICv12",
        
        File_Sheet_Key == "LightningAuto_EPIC-only"   ~ "EPICv1",
        File_Sheet_Key == "LightningManual_EPIC-only" ~ "EPICv1",
        
        File_Sheet_Key == "LightningAuto_EPICv2-only"   ~ "EPICv2",
        File_Sheet_Key == "LightningManual_EPICv2-only" ~ "EPICv2",
        
        File_Sheet_Key == "EZ-DNA-Manual"             ~ "EPICv12",
        File_Sheet_Key == "EZ-DNA-Manual-EPIC-only"   ~ "EPICv1",
        File_Sheet_Key == "EZ-DNA-Manual_EPICv2-only" ~ "EPICv2",
        TRUE ~ NA_character_
      ),
      Sheet_Group = paste( Sheet_Prep,Sheet_Proc,Sheet_Data, sep="_" ),
      
      Sample_Group = Sample_Group %>%
        stringr::str_replace_all(" ",""),
      Sample_Group = dplyr::case_when(
        Sample_Group == "MethylationControl" ~ "MeTritration",
        Sample_Group == "Negative-Mouse" ~ "NegativeMouse",
        Sample_Group == "CellLine" ~ "CellLine",
        Sample_Group == "Coriell" ~ "Coriell",
        Sample_Group == "Blank" ~ "Blank",
        TRUE ~ NA_character_
      ),
      Sentrix_Name = paste(Sentrix_ID,Sentrix_Position, sep="_"),
      Sample_Name = Sample_Name %>% 
        stringr::str_remove("_[0-9]+$") %>%
        stringr::str_remove_all(" ") %>%
        stringr::str_replace_all("-","_"),
      Concentration = dplyr::case_when( 
        Sample_Name %>% stringr::str_detect("_[0-9]+$") ~ Sample_Name %>% stringr::str_remove("^.*_"), 
        TRUE ~ as.character( conc ) ) %>% as.integer(),
      Sample_Base = Sample_Name %>%
        stringr::str_remove("_[0-9]+$"),
      Chip_Name =  dplyr::case_when(
        EPICv2_vs_EPIC == "EPIC" ~ "EPIC",
        EPICv2_vs_EPIC == "EPICv2" ~ "EPIC",
        TRUE ~ NA_character_
      ),
      Chip_Version = dplyr::case_when(
        EPICv2_vs_EPIC == "EPIC" ~ 1,
        EPICv2_vs_EPIC == "EPICv2" ~ 2,
        TRUE ~ NA_real_
      ) %>% as.integer(),
      Manifest_Key = paste0( Chip_Name,"_v",Chip_Version ),
      Sample_Base = Sample_Name %>% 
        stringr::str_remove_all(" ") %>%
        stringr::str_remove("_[0-9]+$"),
      Sample_Name = paste( Sample_Base,Concentration, sep="_"),
    ) %>%
    dplyr::arrange( Chip_Name,Chip_Version, Sample_Group,
                    Sample_Base,Concentration ) %>%
    dplyr::select( File_Sheet_Key,
                   Sheet_Group,Sheet_Prep,Sheet_Proc,Sheet_Data,
                   Sentrix_Name,
                   Chip_Name,Chip_Version,Manifest_Key,
                   Sample_Group,
                   Sample_Base,Concentration,Sample_Name,
                   dplyr::everything() )
  
  if ( sentrix_unq ) ret_tib <- ret_tib %>% 
    dplyr::distinct( Sheet_Prep,Sheet_Proc,Sentrix_Name, .keep_all = TRUE ) %>%
    dplyr::select( 
      Sheet_Prep,Sheet_Proc,Sentrix_Name,
      Chip_Name,Chip_Version,Manifest_Key,
      Sample_Group,Sample_Base,Concentration,Sample_Name,
      Sample_Well,Sample_Plate,Pool_ID )
  
  if ( !is.null(idat_path) ) {
    idat_tib <- sesame::searchIDATprefixes( 
      dir.name = file.path( idat_path ), 
      recursive = TRUE ) %>% 
      as.data.frame() %>% 
      tibble::as_tibble( rownames = "Sentrix_Name" ) %>% 
      purrr::set_names( c("Sentrix_Name","Prefix") )
    
    ret_tib <- ret_tib %>%
      dplyr::left_join( idat_tib, by=c("Sentrix_Name") )
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
  
  # Sample Sheet Summary::
  if ( p4 ) {
    sheet_sum <- NULL
    sheet_sum <- ret_tib %>% 
      dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
      dplyr::group_by( Chip_Name,Chip_Version, Sample_Group,
                       Sample_Base,Concentration ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    sheet_sum %>% print( n=base::nrow(sheet_sum) )
  }
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  Specialized Sample Sheet Preparation::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_epicv2_manifests = function( sesame_csv = NULL,
                                  genome_csv = NULL,
                                  name,
                                  ctls = NULL,
                                  min_cols = TRUE,
                                  write_readr = FALSE,

                                  out_dir,
                                  run_tag,
                                  
                                  reload     = 0,
                                  reload_min = 2,
                                  reload_pre = NULL,
                                  
                                  ret_data   = FALSE,
                                  parallel   = FALSE,
                                  
                                  vb=0, vt=3, tc=1, tt=NULL,
                                  fun_tag='load_epicv2_manifests')
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
  p4  <- vb > vt + 4
  
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
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}   genome_csv = '{genome_csv}'.{RET}"))
    cat(glue::glue("{mssg}   sesame_csv = '{sesame_csv}'.{RET}"))
    cat(glue::glue("{mssg}         name = '{name}'.{RET}"))
    cat(glue::glue("{mssg}         ctls = '{ctls}'.{RET}"))
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
  
  errs_mssg <- glue::glue("Outdir out='{out_dir}' does not exist")
  if ( !dir.exists( out_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              Loop Through Each Parameter:: Sample_Sheet
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Both Genome Studio and Sesame Manifest Strings can't be null")
  if ( ( is.null(genome_csv) || stringr::str_length(genome_csv) == 0 ) &&
       ( is.null(sesame_csv) || stringr::str_length(sesame_csv) == 0 ) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return( NULL )
  
  ses_manifest_vec <- NULL
  if ( !is.null(sesame_csv) && stringr::str_length(sesame_csv) != 0 ) {
    
    ret_tib <- readr::read_csv( file = sesame_csv, show_col_types = FALSE )
    
  } else {
    
    # A tibble: 865,918 × 23 (EPIC v1)
    # A tibble: 969,804 × 24 (EPIC v2)
    gst_tib <- NULL
    gst_tib <- load_genome_studio_manifest( file = genome_csv, 
                                            load_clean    = TRUE,
                                            load_controls = FALSE,
                                            cols_convert  = FALSE,
                                            write_clean   = FALSE,
                                            overwrite     = FALSE,
                                            ret_data      = FALSE,
                                            vb=vb,vt=vt+1,tc=tc,tt=tt )
    
    ret_tib <- gs_to_sesame( tib = gst_tib, 
                             min_cols = min_cols,
                             write_readr = write_readr,
                             
                             out_dir = out_dir, 
                             run_tag = name, 
                             reload  = reload, 
                             reload_min = reload_min,
                             parallel   = parallel,
                             vb=vb,vt=vt+1,tc=tc,tt=tt )
    
  }
  
  if ( !is.null(ctls) && base::nrow(ctls) != 0 )
    ret_tib <- ret_tib %>% dplyr::bind_rows(ctls)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Print Summary::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Local Run Time Defaults::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

swifthoof_sesame_options = function( pars, args,
                                     
                                     vb=0, vt=4, tc=1, tt=NULL,
                                     fun_tag='swifthoof_sesame_options')
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
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}{TAB} Local run_mode='{pars$run_mode}'{RET}"))
    cat(glue::glue("{mssg}{TAB} Local top_path='{pars$top_path}'{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  if ( p2 ) cat(glue::glue("{mssg}   args={args}.{RET}"))
  if ( p0 ) cat(glue::glue("{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Parse Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts <- base::list()
  
  opts$Rscript      <- "Rscript"
  opts$single       <- FALSE
  opts$cluster      <- FALSE
  opts$parallel     <- TRUE
  opts$track_time   <- TRUE
  opts$clean        <- FALSE
  opts$reload       <- 1

  if (pars$run_mode == 'RStudio') {
    
    pars$top_path <- pars$src_path %>% 
      stringr::str_remove("/tools/Workhorse-Unstained/scripts/R")
    opts$top_path <- pars$top_path
    
    if ( is.null(opts$verbose) ) opts$verbose <- pars$verbose
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Params Local Defaults::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Simple Version SetUp::
    if ( is.null(pars$version) ) pars$version <- '0'
    if ( is.null(pars$version_key) ) pars$version_key <- "v"
    opts$version  <- paste0(pars$version_key,pars$version)
    
    # Other common directory defaults::
    opts$out_path <- pars$out_path
    if ( is.null(opts$out_path) )  opts$out_path  <- file.path( pars$top_path, 'scratch' )
    
    opts$idat_path <- pars$idat_path
    if ( is.null(opts$idat_path) ) opts$idat_path <- file.path( pars$top_path, 'data/idats' )
    
    opts$run_name <- pars$run_name
    if ( is.null(opts$run_name) ) opts$run_name <- pars$prgm_tag
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Run_Name Based Options::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( pars$run_name == "GSIBIOINFO-597" ) {
      opts$idat_path <- file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"idats" )
      
      opts$genomestudio <- paste(
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"Manifests/MethylationEPIC_v-1-0_B4.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"Manifests/Interm-PQC-UsedCodesApplied-EPICv2-NoBP28_NoBP29.csv.gz" ),
        sep=COM
      )
      
      opts$sesame <- paste(
        file.path( pars$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/EPIC_v1.gs_to_sesame.csv.gz" ),
        file.path( pars$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/EPIC_v2.gs_to_sesame.csv.gz" ),
        sep=COM
      )
      
      opts$controls <- paste(
        file.path( pars$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" ),
        sep=COM
      )
      
      opts$manifest_name <- paste( "EPIC_v1","EPIC_v2", sep=COM )
      
      opts$sample_sheets <- paste(
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-EPIC.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-EPICv2-Interm-Alpha.csv.gz" ),
        sep=COM
      )
      
    } else if ( pars$run_name == "GSIBIOINFO-638" ) {
      
      opts$idat_path <- file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"idats" )
      
      opts$genomestudio <- paste(
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"Manifests/EPIC-with-ManifestGeneratorTool/EPICv1GRCh37-NA-NA-GRCh37.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"Manifests/EPICv2/EPICv2_SMG_2022_09_12_1125-NA-NA-GRCh37_1A.csv.gz" ),
        sep=COM
      )
      
      opts$sesame <- paste(
        file.path( pars$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full/EPIC_v1.gs_to_sesame.csv.gz" ),
        file.path( pars$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full/EPIC_v2.gs_to_sesame.csv.gz" ),
        sep=COM
      )
      
      opts$controls <- paste(
        file.path( pars$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" ),
        sep=COM
      )
      
      opts$manifest_name <- paste( "EPIC_v1","EPIC_v2", sep=COM )
      
      opts$sample_sheets <- paste(
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-LightningAuto.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-LightningAuto_EPICv2-only.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-LightningAuto_EPIC-only.csv.gz" ),
        
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-LightningManual.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-LightningManual_EPIC-only.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-LightningManual_EPICv2-only.csv.gz" ),
        
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-EZ-DNA-Manual.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-EZ-DNA-Manual-EPIC-only.csv.gz" ),
        file.path( pars$top_path, "Projects/EPIC_v2",pars$run_name,"SampleSheets/SampleSheet-EZ-DNA-Manual_EPICv2-only.csv.gz" ),
        sep=COM
      )
      
    } else if ( pars$run_name == "Embark_v3" ) {
      
      opts$idat_path <- file.path( pars$top_path, "Projects/Embark/Embark_Methylation_Internal_V3" )
      
      opts$sesame <- paste(
        file.path( pars$top_path, "Projects/Embark/Embark_Methylation_Internal_V3/Manifests/unmasked/embark_prbs.06092022.manifest.csv.gz" ),
        sep=COM
      )

      opts$controls <- paste(
        file.path( pars$top_path, "Projects/Embark/Embark_Methylation_Internal_V3/Manifests/unmasked/embark_ctls.06092022.manifest.rds" ),
        sep=COM
      )
      
      opts$manifest_name <- paste( "Embark_v3", sep=COM )
      
      opts$sample_sheets <- paste(
        file.path( pars$top_path, "Projects/Embark/Embark_Methylation_Internal_V3/Sample_Sheets/Formatted/Technical_Non_Technical_Replicates.sample_sheet.08.12.2022.updated-dog-names.MeTitration.csv.gz" ),
        sep=COM
      )
      # opts$sample_sheets <- paste(
      #   file.path( pars$top_path, "Projects/Embark/Embark_Methylation_Internal_V3/Experiment_MethylationTitration/MethylationTitration_PDR.csv.gz" ),
      #   file.path( pars$top_path, "Projects/Embark/Embark_Methylation_Internal_V3/Sample_Sheets/Formatted/Technical_Non_Technical_Replicates.sample_sheet.08.12.2022.updated-dog-names.csv.gz" ),
      #   sep=COM
      # )
      
    } else {
      stop( glue::glue("{errs} Unrecognized run_name = '{pars$run_name}'!{RET2}") )
      return(NULL)
    }
    opts$run_name <- paste(opts$run_name,opts$version, sep='-')
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Validate Options::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Ensure idats directory exists!!!
    errs_mssg <- glue::glue("Idats Directory does NOT exist: {opts$idat_path}")
    if ( !dir.exists( opts$idat_path ) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
    
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
        c("--out_path"), type="character", default=opts$out_path,
        help=glue::glue(
          "Build output directory path.",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--genomestudio"), type="character", default=opts$genomestudio, 
        help=glue::glue(
          "Genome Studio Manifests full path(s) ",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--sesame"), type="character", default=opts$sesame, 
        help=glue::glue(
          "Sesame Manifests full path(s) ",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--controls"), type="character", default=opts$controls, 
        help=glue::glue(
          "Sesame Controls Manifests full path(s) ",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--sesame"), type="character", default=opts$manifest_name, 
        help=glue::glue(
          "Manifest names(s) used as look up key in sample sheet",
          "{RET}{TAB2}[ Comma separated list ]",
          "{RET}{TAB2}[ default = %default ]" ),
        metavar="character"),
      
      optparse::make_option(
        c("--sample_sheets"), type="character", default=opts$sample_sheets, 
        help=glue::glue(
          "Directory with multiple sample_sheets to auto-detect or full path(s) ",
          "{RET}{TAB2} to sample_sheet(s).",
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
    
    # Read Options from pre-defined function::
    # option_list <- program_options(pars = par, opts = opt, vb = opts$verbose)
    
    opt_parser = optparse::OptionParser( option_list = options_list )
    opts = optparse::parse_args( opt_parser )
    
  } else {
    stop( glue::glue("{errs} Unrecognized run_mode = '{pars$run_mode}'!{RET2}") )
    return(NULL)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Done Parsing Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_cnt <- length(opts)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

# End of file

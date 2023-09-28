# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Signal Data Frame (sdf)::
#                      Base (Simple) Sesame Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Options, Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

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
S10  <- paste0("  ",
               paste(rep("     ",2),collapse=" "))
S15  <- paste0("  ",
               paste(rep("     ",3),collapse=" "))
S20  <- paste0("  ",
               paste(rep("     ",4),collapse=" "))
S25  <- paste0("  ",
               paste(rep("     ",5),collapse=" "))
S30  <- paste0("  ",
               paste(rep("     ",6),collapse=" "))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   GS <=> Sesame Manifest Conversions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

gs_to_sesame = function( tib,
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
                         fun_tag='gs_to_sesame')
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
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  if ( !file.exists(beg_txt) )
    sys_ret <- base::system( glue::glue("touch {beg_txt}") )
  
  ret_tib <- tib %>%
    dplyr::rename(
      Probe_ID = IlmnID,
      U = AddressA_ID,
      M = AddressB_ID
    ) %>%
    dplyr::mutate( 
      mask = FALSE,
      dplyr::across( col, ~tidyr::replace_na(.x, "2" ) )
    )
  
  if ( min_cols ) ret_tib <- ret_tib %>%
    dplyr::select( Probe_ID,U,M,col,mask, AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  if ( !min_cols ) ret_tib <- ret_tib %>%
    dplyr::select( Probe_ID,U,M,col,mask, AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                   Name,Next_Base,Color_Channel,col,Probe_Type,
                   Strand_FR,Strand_TB,Strand_CO,Infinium_Design,
                   Infinium_Design_Type,CHR,MAPINFO,Species,
                   Genome_Build,Source_Seq,Forward_Sequence,
                   Top_Sequence,Rep_Num )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Write Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  out_cnt <- ret_tib %>% base::nrow()
  
  if ( write_readr ) {
    readr::write_csv( x = ret_tib, file = out_csv )
  } else {
    out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  }
  
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
#                  Base Wrapper for Prefix to Raw SDFs::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prefixes_to_rds_raw = function( ssh_list,
                                man_tibs,
                                out_dir,
                                
                                suffix_str    = ".prefix_to_sdf.rds$",
                                prefix_key    = "Prefix",
                                chip_name_key = "Chip_Name",
                                manifest_key  = "Manifest_Key",
                                sentrix_key   = "Sentrix_Name",
                                
                                reload     = 0,
                                reload_min = 10,
                                reload_pre = NULL,
                                
                                ret_data   = FALSE,
                                parallel   = FALSE,
                                
                                vb=0, vt=6, tc=1, tt=NULL,
                                fun_tag='prefixes_to_rds_raw')
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
  
  prefix_cnt <- ssh_list %>% length()
  
  mc_valid <- validate_mclapply1()
  if ( !mc_valid ) {
    stop(glue::glue("{mssg} Failed validate_mclapply1()!{RET}"))
    return( ret_dat )
  }
  errs_mssg <- glue::glue("Prefix Count is zero({prefix_cnt})")
  if ( is.null(ssh_list) || prefix_cnt == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  sdf_dir  <- safe_mkdir( dir = file.path( out_dir,"prefix_to_sdf") )
  rds_list <- list.files( path = sdf_dir, pattern = ".prefix_to_sdf.rds$" )
  rds_cnt  <- length(rds_list)
  
  is_valid <- FALSE
  if ( reload >= reload_min && rds_cnt == prefix_cnt ) is_valid <- TRUE
    
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}   prefix_cnt = '{prefix_cnt}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      rds_cnt = '{rds_cnt}'.{RET}"))
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
  }
  
  if ( is_valid ) {
    if ( p2 ) cat(glue::glue("{mssg} Loading Raw SDFs={sdf_dir}...{RET}"))
    
    ftime <- base::system.time({
      rds_list <- NULL
      rds_list <- file_list( 
        path = sdf_dir,
        suffix = suffix_str, 
        pattern = suffix_str )
      
      ret_dat <- NULL
      ret_dat <- foreach::foreach( 
        sentrix_name = names(rds_list), 
        .inorder=TRUE, 
        .final = function(x) setNames(x, names(rds_list)) ) %dopar% {
          readr::read_rds( file = rds_list[[sentrix_name]] )
        }
    })
    tt$addTime( ftime, paste( "load_sdfs_time", sep="." ) )
    
  } else {
    #  TBD:: Replace prefix_to_sdf() with c++
    if ( opt$parallel ) {
      if ( p2 ) cat(glue::glue("{mssg} Building (parallel) Raw SDFs={sdf_dir}...{RET}"))
      
      ftime <- base::system.time({
        ret_dat <- NULL
        ret_dat <- foreach::foreach( sn = names(ssh_list), 
                                     .inorder=TRUE, 
                                     .final = function(x) setNames(x, names(ssh_list)) ) %dopar% {
          prefix_to_sdf( prefix     = ssh_list[[sn]][[prefix_key]],
                         platform   = ssh_list[[sn]][[chip_name_key]], 
                         manifest   = man_tibs[[ ssh_list[[sn]][[manifest_key]] ]],
                         out_dir    = out_dir,
                         run_tag    = ssh_list[[sn]][[sentrix_key]],
                         reload     = 0,
                         reload_min = 10, 
                         parallel   = FALSE,
                         vb=vb,vt=vt+3,tc=tc,tt=tt )
        }
      })
      tt$addTime( ftime, paste( "make_sdfs_time_par", sep="." ) )
      
    } else {
      if ( p2 ) cat(glue::glue("{mssg} Building (linear) Raw SDFs={sdf_dir}...{RET}"))
      
      ftime <- base::system.time({
        ret_dat <- NULL
        ret_dat <- ssh_list %>% # head(n=3) %>%
          lapply( function(x) {
            prefix_to_sdf( prefix     = x[[prefix_key]],
                           platform   = x[[chip_name_key]],
                           manifest   = man_tibs[[ x[[manifest_key]] ]],
                           out_dir    = out_dir,
                           run_tag    = x[[sentrix_key]],
                           reload     = 0,
                           reload_min = 10,
                           parallel   = FALSE,
                           vb=vb,vt=vt+3,tc=tc,tt=tt )
          } )
      })
      tt$addTime( ftime, paste( "make_sdfs_time_lin", sep="." ) )
    }
  }

  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_dat, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#              Another Implementation of Sesame-On-Demand::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#  
# > sesame::prepSesameList()
# code                  func                                description
# 1     0              resetMask                    Reset mask to all FALSE
# 2     Q            qualityMask                 Mask probes of poor design
# 3     G        prefixMaskButCG                    Mask all but cg- probes
# 4     H         prefixMaskButC             Mask all but cg- and ch-probes
# 5     C* inferInfiniumIChannel        Infer channel for Infinium-I probes  [ switch_failed = TRUE ]
# 5     c* inferInfiniumIChannel        Infer channel for Infinium-I probes  [ switch_failed = FALSE ]
# 6     D*             dyeBiasNL           Dye bias correction (non-linear)  [ mask = TRUE ]
# 6     d*             dyeBiasNL           Dye bias correction (non-linear)  [ mask = FALSE ]
# 7     E               dyeBiasL               Dye bias correction (linear)  
# 8     P*                pOOBAH        Detection p-value masking using oob  [ return.pval = FALSE, combine.neg = TRUE,  pval.threshold = min_pval ]
# 8     p*                pOOBAH        Detection p-value masking using oob  [ return.pval = FALSE, combine.neg = FALSE, pval.threshold = min_pval ]
# 9     I            detectionIB Mask detection by intermediate beta values
# 10    B                   noob           Background subtraction using oob  [ combine.neg = TRUE,  offset = 15 ]
# 10    b                   noob           Background subtraction using oob  [ combine.neg = FALSE, offset = 15 ]
# 11    S           inferSpecies                  Set species-specific mask
# 12    T            inferStrain           Set strain-specific mask (mouse)
# 13    M            matchDesign        Match Inf-I/II in beta distribution
#
# 14    N*    detectionPnegEcdf2        Detection p-value masking using neg  [ return.pval = FALSE, pval.threshold = min_pval, use_type = FALSE ]
# 15    V*              getBetas                         Return Beta Values  [ mask = TRUE,  sum.TypeI = FALSE, collapseToPfx = FALSE ]
# 15    v*              getBetas                         Return Beta Values  [ mask = FALSE, sum.TypeI = FALSE, collapseToPfx = FALSE ]
# 15    W*              getBetas                         Return Beta Values  [ mask = TRUE,  sum.TypeI = TRUE,  collapseToPfx = FALSE ]
# 15    w*              getBetas                         Return Beta Values  [ mask = FALSE, sum.TypeI = TRUE,  collapseToPfx = FALSE ]
#
# 16    O*                pOOBAH          Detection p-value return using oob [ return.pval = TRUE, combine.neg = TRUE,  pval.threshold = poob_min ]
# 16    o*                pOOBAH          Detection p-value return using oob [ return.pval = TRUE, combine.neg = FALSE, pval.threshold = poob_min ]
#

detectionPnegEcdf2 = function( sdf, return.pval = FALSE, pval.threshold=0.05, use_type = FALSE ) {
  
  if ( use_type ) {
    negctls <- sdf %>% sesame::controls() %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate( Type = Type %>% stringr::str_to_upper() ) %>%
      dplyr::filter( Type %>% stringr::str_detect("NEGATIVE") )
    # return( negctls )
  } else {
    negctls <- sdf %>% sesame::controls() %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate( Type = Probe_ID %>% stringr::str_to_upper() ) %>%
      dplyr::filter( Type %>% stringr::str_detect("NEGATIVE") )
  }
  
  funcG <- stats::ecdf( negctls$UG )
  funcR <- stats::ecdf( negctls$UR )
  
  ## p-value is the minimium detection p-value of the 2 alleles
  pvals <- stats::setNames(BiocGenerics::pmin(
    1-funcR(BiocGenerics::pmax(sdf$MR, sdf$UR, na.rm=TRUE)),
    1-funcG(BiocGenerics::pmax(sdf$MG, sdf$UG, na.rm=TRUE))), sdf$Probe_ID)
  
  if (return.pval) { return(pvals) }
  
  sesame::addMask(sdf, pvals > pval.threshold)
}

mutate_sdf_linear = function( sdf, 
                              steps,
                              mask_ids = NULL,
                              
                              min_pval = 0.05, 
                              off_set  = 15, 
                              ran_lapply = FALSE,
                              
                              out_dir,
                              run_tag,
                              
                              reload     = 0,
                              reload_min = 10,
                              reload_pre = NULL,
                              
                              ret_data   = FALSE,
                              parallel   = FALSE,
                              
                              vb=0, vt=3, tc=1, tt=NULL,
                              fun_tag='mutate_sdf_linear')
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
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}        steps = '{steps}'.{RET}"))
    cat(glue::glue("{mssg}     min_pval = '{min_pval}'.{RET}"))
    cat(glue::glue("{mssg}      off_set = '{off_set}'.{RET}"))
    cat(glue::glue("{mssg}   ran_lapply = '{ran_lapply}'.{RET}"))
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
  
  sdf_cnt <- sdf %>% base::nrow()
  errs_mssg <- glue::glue("SDF is empty: sdf_cnt='{sdf_cnt}'!{RET}")
  if ( is.null( sdf) || sdf_cnt == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  if ( p1 ) cat(glue::glue("{mssg} SDF Row Count='{sdf_cnt}'.{RET}"))
  
  errs_mssg <- glue::glue("Min Pval='{min_pval}' not defined")
  if ( is.null( min_pval) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  if ( p1 ) cat(glue::glue("{mssg} Set Min Pvalue='{min_pval}'.{RET}"))
  
  step_vec <- NULL
  step_vec <- stringr::str_split( string = steps, 
                                  pattern = "", 
                                  n = stringr::str_length(steps), 
                                  simplify = FALSE ) %>% 
    base::unlist() %>% as.vector()
  step_cnt = step_vec %>% length()
  
  errs_mssg <- glue::glue("Failed to parse workflow steps='{steps}'. Step Counts='{step_cnt}'{RET}")
  if ( is.null(step_vec) || step_cnt == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  if ( p1 ) cat(glue::glue("{mssg} Total Number of Workflow Steps='{step_cnt}'.{RET}"))
  if ( p1 ) cat( glue::glue("{mssg} Inputs are valid!{RET2}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    #
    # Validate sdf is legit::
    #  TBD:: Understand why this fails somestimes...
    #
    nan_dat <- NULL
    nan_dat <- sesame::sesameQC_calcStats( sdf = sdf, "dyeBias" )@stat$RGdistort
    
    warn_mssg <- glue::glue("SesameQC_caclcStats():: is.na()!")
    if ( is.null(nan_dat) || is.na( nan_dat ) ) wflag <- TRUE
    if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    if ( wflag ) nan_dat %>% print()
    if ( wflag ) return( NULL )
    
    for ( s in step_vec ) {
      if ( p2 ) cat(glue::glue("{mssg} work={steps}, step='{s}'{RET}"))
      
      if ( s == "A" ) sdf = sesame::addMask( sdf = sdf, probes = mask_ids )
      if ( s == "C" ) sdf = sesame::inferInfiniumIChannel( sdf = sdf, switch_failed = FALSE, verbose = FALSE )
      if ( s == "c" ) sdf = sesame::inferInfiniumIChannel( sdf = sdf, switch_failed = TRUE,  verbose = FALSE )
      # if ( s == "D" ) sdf = sesame::dyeBiasCorrTypeINorm(  sdf = sdf, mask = TRUE )
      # if ( s == "d" ) sdf = sesame::dyeBiasCorrTypeINorm(  sdf = sdf, mask = FALSE )
      if ( s == "D" ) sdf = sesame::dyeBiasNL( sdf = sdf, mask = TRUE )
      if ( s == "d" ) sdf = sesame::dyeBiasNL( sdf = sdf, mask = FALSE )
      if ( s == "P" ) sdf = sesame::pOOBAH(                sdf = sdf, return.pval = FALSE, combine.neg = TRUE,  pval.threshold = min_pval )
      if ( s == "P" ) sdf = sesame::pOOBAH(                sdf = sdf, return.pval = FALSE, combine.neg = FALSE, pval.threshold = min_pval )
      if ( s == "O" ) sdf = sesame::pOOBAH(                sdf = sdf, return.pval = TRUE,  combine.neg = TRUE,  pval.threshold = min_pval )
      if ( s == "o" ) sdf = sesame::pOOBAH(                sdf = sdf, return.pval = TRUE,  combine.neg = FALSE, pval.threshold = min_pval )
      if ( s == "B" ) sdf = sesame::noob(                  sdf = sdf, combine.neg = TRUE,  offset = off_set )
      if ( s == "b" ) sdf = sesame::noob(                  sdf = sdf, combine.neg = FALSE, offset = off_set )
      if ( s == "N" ) sdf = detectionPnegEcdf2(            sdf = sdf, return.pval = FALSE, pval.threshold = min_pval, use_type = FALSE )
      if ( s == "n" ) sdf = detectionPnegEcdf2(            sdf = sdf, return.pval = TRUE,  pval.threshold = min_pval, use_type = FALSE )
      # if ( s == "N" ) sdf = sesame::detectionPnegEcdf(    sdf = sdf, return.pval = FALSE, pval.threshold = min_pval, use_type = FALSE )
      if ( s == "V" ) sdf = sesame::getBetas( sdf = sdf, mask = TRUE,  sum.TypeI = FALSE, collapseToPfx = FALSE )
      if ( s == "v" ) sdf = sesame::getBetas( sdf = sdf, mask = FALSE, sum.TypeI = FALSE, collapseToPfx = FALSE )
      if ( s == "W" ) sdf = sesame::getBetas( sdf = sdf, mask = TRUE,  sum.TypeI = TRUE,  collapseToPfx = FALSE )
      if ( s == "w" ) sdf = sesame::getBetas( sdf = sdf, mask = FALSE, sum.TypeI = TRUE,  collapseToPfx = FALSE )
      
      if ( p2 ) cat(glue::glue("{mssg} Done.{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # out_cnt <- safe_write( x = sdf, file = out_csv, type = "csv", 
    #                        done = TRUE, write_spec = TRUE, append = FALSE, 
    #                        fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( sdf, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  sdf
}

mutate_sdf_simple = function( sdf, 
                              steps,
                              negs_min = 0.05, 
                              poob_min = 0.05, 
                              off_set  = 15,
                              
                              vb=0, vt=3, tc=1, tt=NULL,
                              fun_tag='mutate_sdf_simple')
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
  p3  <- vb > vt + 3
  
  step_vec <- NULL
  step_vec <- stringr::str_split( string = steps, 
                                  pattern = "", 
                                  n = stringr::str_length(steps), 
                                  simplify = FALSE ) %>% 
    base::unlist() %>% as.vector()
  step_cnt = step_vec %>% length()
  
  if ( p2 ) cat(glue::glue("{mssg} work={steps}, negs/poob='{negs_min}/{poob_min}'{RET}"))
  for ( s in step_vec ) {
    if ( p3 ) cat(glue::glue("{mssg}{TAB}step='{s}'{RET}"))
    
    if ( s == "A" ) sdf = sesame::addMask( sdf = sdf, probes = mask_ids )
    if ( s == "C" ) sdf = sesame::inferInfiniumIChannel( sdf = sdf, switch_failed = FALSE, verbose = FALSE )
    if ( s == "c" ) sdf = sesame::inferInfiniumIChannel( sdf = sdf, switch_failed = TRUE,  verbose = FALSE )
    # if ( s == "D" ) sdf = sesame::dyeBiasCorrTypeINorm(  sdf = sdf, mask = TRUE )
    # if ( s == "d" ) sdf = sesame::dyeBiasCorrTypeINorm(  sdf = sdf, mask = FALSE )
    if ( s == "D" ) sdf = sesame::dyeBiasNL( sdf = sdf, mask = TRUE )
    if ( s == "d" ) sdf = sesame::dyeBiasNL( sdf = sdf, mask = FALSE )
    if ( s == "P" ) sdf = sesame::pOOBAH(                sdf = sdf, return.pval = FALSE, combine.neg = TRUE,  pval.threshold = poob_min )
    if ( s == "P" ) sdf = sesame::pOOBAH(                sdf = sdf, return.pval = FALSE, combine.neg = FALSE, pval.threshold = poob_min )
    if ( s == "O" ) sdf = sesame::pOOBAH(                sdf = sdf, return.pval = TRUE,  combine.neg = TRUE,  pval.threshold = poob_min )
    if ( s == "o" ) sdf = sesame::pOOBAH(                sdf = sdf, return.pval = TRUE,  combine.neg = FALSE, pval.threshold = poob_min )
    if ( s == "B" ) sdf = sesame::noob(                  sdf = sdf, combine.neg = TRUE,  offset = off_set )
    if ( s == "b" ) sdf = sesame::noob(                  sdf = sdf, combine.neg = FALSE, offset = off_set )
    if ( s == "N" ) sdf = detectionPnegEcdf2(            sdf = sdf, return.pval = FALSE, pval.threshold = negs_min, use_type = FALSE )
    if ( s == "n" ) sdf = detectionPnegEcdf2(            sdf = sdf, return.pval = TRUE,  pval.threshold = negs_min, use_type = FALSE )
    if ( s == "V" ) sdf = sesame::getBetas( sdf = sdf, mask = TRUE,  sum.TypeI = FALSE, collapseToPfx = FALSE )
    if ( s == "v" ) sdf = sesame::getBetas( sdf = sdf, mask = FALSE, sum.TypeI = FALSE, collapseToPfx = FALSE )
    if ( s == "W" ) sdf = sesame::getBetas( sdf = sdf, mask = TRUE,  sum.TypeI = TRUE,  collapseToPfx = FALSE )
    if ( s == "w" ) sdf = sesame::getBetas( sdf = sdf, mask = FALSE, sum.TypeI = TRUE,  collapseToPfx = FALSE )
  }
  if ( p2 ) cat(glue::glue("{mssg} Done.{RET2}"))
  
  sdf
}

mutate_sdfs = function( sdfs,

                        work_str = "i",
                        negs_min = 1,
                        poob_min = 1,
                        
                        sdfs_max = 0,
                        
                        out_dir,
                        run_tag,
                        
                        reload     = 0,
                        reload_min = 2,
                        reload_pre = NULL,
                        
                        ret_data   = FALSE,
                        parallel   = FALSE,
                        
                        vb=0, vt=3, tc=1, tt=NULL,
                        fun_tag='mutate_sdfs')
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
  
  work_str <- paste0( "i", work_str %>% stringr::str_remove("^i") )
  work_tag <- paste( "work",work_str, sep='-' )
  negs_tag <- paste( "negs",negs_min, sep='-' )
  poob_tag <- paste( "poob",poob_min, sep='-' )
  
  # out_dir <- safe_mkdir( file.path( out_dir,fun_tag, work_str,negs_tag,poob_tag ), 
  #                        vb=vb,vt=vt+1,tc=tc,tt=tt )
  out_dir <- file.path( out_dir, fun_tag, run_tag, work_tag,negs_tag,poob_tag )
  out_tag <- paste( fun_tag,run_tag, work_tag,negs_tag,poob_tag, sep='.' )
  
  out_pre <- file.path( out_dir, paste(out_tag, sep='.') )
  out_rds <- file.path( paste(out_pre, 'rds', sep='.') )
  out_csv <- file.path( paste(out_pre, 'csv.gz', sep='.') )
  
  beg_txt <- file.path( paste(out_rds, 'start.txt', sep='.') )
  end_txt <- file.path( paste(out_rds, 'done.txt', sep='.') )
  
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_rds, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( readr::read_rds( file = out_rds) )
  # if ( reload >= reload_min && is_valid )
  #   return( safe_read( out_rds, fun_tag = fun_tag, head = "Reloading",
  #                      vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}     work_str = '{work_str}'.{RET}"))
    cat(glue::glue("{mssg}     negs_min = '{negs_min}'.{RET}"))
    cat(glue::glue("{mssg}     poob_min = '{poob_min}'.{RET}"))
    cat(glue::glue("{mssg}     sdfs_max = '{sdfs_max}'.{RET}"))
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
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET2}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      out_rds = '{out_rds}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(out_rds, end_txt) )
  
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
  
  errs_mssg <- glue::glue("SDFs list is empty")
  if ( is.null(sdfs) || length(sdfs)==0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    if ( is.numeric(sdfs_max) && sdfs_max > 0 ) sdfs <- sdfs %>% head( n=sdfs_max )
    
    if ( parallel ) {
      ret_dat <- foreach::foreach( sn = names(sdfs), .inorder=TRUE, .combine = cbind ) %dopar% {
        mutate_sdf_simple( sdf = sdfs[[sn]], 
                           steps = work_str,
                           negs_min = negs_min,
                           poob_min = poob_min,
                           vb=vb,vt=vt+3,tc=tc,tt=tt )
      }
    } else {
      ret_dat <- foreach::foreach( sn = names(sdfs), .inorder=TRUE, .combine = cbind ) %do% {
        mutate_sdf_simple( sdf = sdfs[[sn]], 
                           steps = work_str,
                           negs_min = negs_min,
                           poob_min = poob_min,
                           vb=vb,vt=vt+3,tc=tc,tt=tt )
      }
    }
    colnames( ret_dat ) <- names(sdfs)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_dat, file = out_rds, type = "rds",
                           done = TRUE, write_spec = FALSE, append = FALSE,
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_dat, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

train_test_sdfs = function( sdfs,
                        ssh,
                        
                        # rep_vec = NULL,
                        # uhm_vec = NULL,
                        # val_vec = NULL,
                        
                        grp_key = "Sample_Group",
                        rep_str = "TechnicalReplicates",
                        uhm_str = "MeTritration",
                        val_str = NULL,
                        # val_str = "Replicates",
                        
                        work_str = "i",
                        negs_min = 1,
                        poob_min = 1,
                        
                        out_dir,
                        run_tag,
                        
                        reload     = 0,
                        reload_min = 2,
                        reload_pre = NULL,
                        
                        ret_data   = FALSE,
                        parallel   = FALSE,
                        
                        vb=0, vt=3, tc=1, tt=NULL,
                        fun_tag='train_test_sdfs')
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
  
  negs_tag <- paste( "negs-min",negs_min, sep='-' )
  poob_tag <- paste( "poob-min",poob_min, sep='-' )
  
  out_dir <- safe_mkdir( file.path( out_dir,fun_tag, work_str,negs_tag,poob_tag ), 
                         vb=vb,vt=vt+1,tc=tc,tt=tt )
  out_tag <- paste( run_tag, fun_tag, work_str,negs_tag,poob_tag, sep='.' )
  
  out_pre <- file.path( out_dir, paste(out_tag, sep='.') )
  beg_txt <- file.path( paste(out_pre, 'beg.txt', sep='.') )
  end_txt <- file.path( paste(out_pre, 'end.txt', sep='.') )
  
  out_rds <- file.path( paste(out_pre, 'rds', sep='.') )
  out_csv <- file.path( paste(out_pre, 'csv.gz', sep='.') )
  
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_rds, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_rds, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}     work_str = '{work_str}'.{RET}"))
    cat(glue::glue("{mssg}     negs_min = '{negs_min}'.{RET}"))
    cat(glue::glue("{mssg}     poob_min = '{poob_min}'.{RET}"))
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
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET2}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      out_rds = '{out_rds}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(out_rds, end_txt) )
  
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
  
  errs_mssg <- glue::glue("SDFs list is empty")
  if ( is.null(sdfs) || length(sdfs)==0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("Sample Sheet is Empty")
  if ( is.null(ssh) || base::nrow(ssh)==0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  # negs_ssh_tib %>% tib_group_split( vec = "Sample_Group" )
  ssh_list <- ssh %>% tib_group_split( vec = grp_key )
  
  errs_mssg <- glue::glue("Failed to find Training Replicates Key='{rep_str}'")
  if ( (!is.null(rep_str) && stringr::str_length(rep_str) != 0) &&
       is.null(ssh_list[[rep_str]]) || base::nrow(ssh_list[[rep_str]]) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("Failed to find Training UHM Key='{rep_str}'")
  if ( (!is.null(uhm_str) && stringr::str_length(uhm_str) != 0) &&
       is.null(ssh_list[[uhm_str]]) || base::nrow(ssh_list[[uhm_str]]) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  # errs_mssg <- glue::glue("Training and Testing Vectors are All Emplty")
  # if ( (is.null(rep_vec) || length(rep_vec)==0) &&
  #      (is.null(uhm_vec) || length(uhm_vec)==0) &&
  #      (is.null(val_vec) || length(val_vec)==0) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
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



# End of file

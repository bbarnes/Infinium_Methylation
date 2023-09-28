
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         Signal Data Frame (sdf)::
#                         Latest Sesame Functions::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )

# Parallel Processing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

# Plotting Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("ggplot2", quietly = TRUE) ) )

# Matrix Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("matrixStats", quietly = TRUE) ) )

# Scales Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("scales",      quietly = TRUE) ))
# 

# [TBD]: REMOVE OLD STUFF:
#

# Command Line Options Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("optparse",   quietly = TRUE) ) )
#

# Tidyverse Core Packages::
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("tidyverse",  quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("stringr", quietly = TRUE) ))
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("readr",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("glue",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("furrr",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("purrr",    quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("magrittr", quietly = TRUE) ) )

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
#                     Pairs Plot (R-Squred/Detla-Beta)::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

beta_panel_gg = function( data, mapping, 
                          plot_type, # c( "dB","r2","dg" )
                          pval_min = 1.0,
                          beta_min = 0.3,
                          beta_max = 0.7,
                          dB_min   = 0.2,
                          
                          beta_data = NULL,
                          pval_data = NULL,
                          pval_cols = NULL,
                          full_plot = FALSE,
                          
                          grp_key = "Groups", 
                          grp_vec = NULL, 
                          
                          alpha = 0.3, alpha_lab = 0.3,
                          size = 0.5, wsize = 2, ticks = 3,
                          
                          # r2 defaults below, but they work for dB as well...
                          x_min = 0, x_max = 1,
                          y_min = 0, y_max = 1,
                          
                          out_dir,
                          run_tag,
                          
                          reload     = 0,
                          reload_min = 2,
                          reload_pre = NULL,
                          
                          ret_data   = FALSE,
                          parallel   = FALSE,
                          
                          vb=0, vt=3, tc=1, tt=NULL,
                          fun_tag='beta_panel_gg')
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
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1 )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1, tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
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
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Build Groups::
    #                        Build All dB Combinations::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    grp_sym <- rlang::sym( grp_key )
    
    plot_data <- data
    calc_data <- data
    
    if ( !is.null(beta_data) ) calc_data <- beta_data
    if ( !is.null(beta_data) && full_plot ) plot_data <- beta_data
    
    if ( !is.null(grp_vec) && grp_vec %>% length() != 0 ) {
      plot_data <- plot_data %>% 
        dplyr::mutate( !!grp_sym := paste( grp_vec, collapse = "-" ) )
      
      calc_data <- calc_data %>% 
        dplyr::mutate( !!grp_sym := paste( grp_vec, collapse = "-" ) )
    }
    
    if ( plot_type == "dB" ) {
      
      plot_data <- plot_data %>% dplyr::mutate( 
        dB = GGally::eval_data_col(plot_data, mapping$x) - GGally::eval_data_col(plot_data, mapping$y)
      )
      calc_data <- calc_data %>% dplyr::mutate( 
        dB = GGally::eval_data_col(calc_data, mapping$x) - GGally::eval_data_col(calc_data, mapping$y)
      )
      
    } else if ( plot_type == "r2" ) {
      
      calc_data <- calc_data %>% dplyr::mutate( 
        bt1 = GGally::eval_data_col(calc_data, mapping$x),
        bt2 = GGally::eval_data_col(calc_data, mapping$y)
      )
      
    } else if ( plot_type == "dg" ) {
      
    } else {
      errs_mssg <- glue::glue("Unsupported plot -type = '{plot_type}'")
      if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    print_tib( plot_data, name = "dB Plot Data", vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    print_tib( calc_data, name = "dB Calc Data", vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Calculate Stats to be Displayed in Text::
    #                                    &
    #                      Calculate Plotting Boundary Limits::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( plot_type == "dB" ) {
      
      calc_sum <- calc_data %>%
        dplyr::group_by( dplyr::across( dplyr::all_of( grp_vec ) ) ) %>%
        dplyr::summarise( Total_cnt = n(),
                          Total_str = number_as_commaK(Total_cnt),
                          
                          na_Fail_cnt = sum( is.na(dB) ),
                          dB_Fail_cnt = sum( dB >  dB_min, na.rm=TRUE ),
                          dB_Pass_cnt = sum( dB <= dB_min, na.rm=TRUE ),
                          
                          na_Fail_str = paste0( round( 100*na_Fail_cnt/Total_cnt, 0 ),"%" ),
                          dB_Fail_str = paste0( round( 100*dB_Fail_cnt/Total_cnt, 0 ),"%" ),
                          dB_Pass_str = paste0( round( 100*dB_Pass_cnt/Total_cnt, 0 ),"%" ),
                          
                          .groups = "drop" ) %>%
        dplyr::select( dplyr::all_of( grp_vec ), dplyr::ends_with("_str") )
      
      #  dplyr::select( dplyr::all_of(grp_vec), Total_cnt, dplyr::ends_with("_str") )
      
      x_min <- min( plot_data$dB, na.rm=TRUE )
      x_max <- max( plot_data$dB, na.rm=TRUE )
      
      y_idx  <- which.max( density( plot_data$dB, na.rm=TRUE )$y )
      y_max  <- density( plot_data$dB, na.rm=TRUE )$y[y_idx]
      
    } else if ( plot_type == "r2" ) {
      
      calc_sum <- calc_data %>%
        dplyr::group_by( dplyr::across( dplyr::all_of( grp_vec ) ) ) %>%
        dplyr::summarise( r2 = round( cor( 
          bt1, bt2, method='pearson', use='pairwise.complete.obs' ), 3 ),
          .groups = "drop" )
      
    } else if ( plot_type == "dg" ) {
      
      errs_mssg <- glue::glue("Pval Data cannot be null or length zero when ",
                              "plot_type='{plot_type}'")
      if ( is.null(pval_data) || base::nrow(pval_data) == 0 ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
      
      #
      # TBD:: Test this it shouldn't ever happen...
      #
      calc_size <- base::nrow(calc_data)
      pval_size <- base::nrow(pval_data)
      if ( calc_size != pval_size ) {
        errs_mssg <- glue::glue("We shouldn't be here; calc={calc_size}, ",
                                "pval={pval_size}")
        if ( TRUE ) eflag <- TRUE
        if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
        if ( eflag ) return(NULL)
        
        calc_data <- calc_data %>% dplyr::filter( Probe_ID %in% pval_data$Probe_ID )
        pval_data <- pval_data %>% dplyr::filter( Probe_ID %in% calc_data$Probe_ID )
      }
      
      calc_sum <- calc_data %>%
        dplyr::mutate( pv = GGally::eval_data_col(pval_data, mapping$x),
                       bt = GGally::eval_data_col(calc_data, mapping$x) ) %>%
        dplyr::group_by( dplyr::across( dplyr::all_of( grp_vec ) ) ) %>%
        dplyr::summarise( Total_cnt = n(),
                          Total_str = number_as_commaK(Total_cnt),
                          
                          na_Fail_cnt = sum( is.na(pv) ),
                          pv_Fail_cnt = sum( pv >  pval_min, na.rm=TRUE ),
                          pv_Pass_cnt = sum( pv <= pval_min, na.rm=TRUE ),
                          
                          bt_hypo_cnt = sum( bt < beta_min, na.rm=TRUE ),
                          bt_hypr_cnt = sum( bt > beta_max, na.rm=TRUE ),
                          
                          na_Fail_str = paste0( round( 100*na_Fail_cnt/Total_cnt, 0 ),"%" ),
                          pv_Fail_str = paste0( round( 100*pv_Fail_cnt/Total_cnt, 0 ),"%" ),
                          pv_Pass_str = paste0( round( 100*pv_Pass_cnt/Total_cnt, 0 ),"%" ),
                          
                          bt_hypo_str = paste0( round( 100*bt_hypo_cnt/Total_cnt, 0 ),"%" ),
                          bt_hypr_str = paste0( round( 100*bt_hypr_cnt/Total_cnt, 0 ),"%" ),
                          
                          .groups = "drop" ) %>%
        dplyr::select( dplyr::all_of( grp_vec ), dplyr::ends_with("_str") )
      
      x_min <- min( GGally::eval_data_col(plot_data, mapping$x), na.rm=TRUE)
      x_max <- max( GGally::eval_data_col(plot_data, mapping$x), na.rm=TRUE)
      
      y_min <- 0
      y_idx  <- which.max(density( GGally::eval_data_col(plot_data, mapping$x), na.rm=TRUE)$y)
      y_max  <- density( GGally::eval_data_col(plot_data, mapping$x), na.rm=TRUE)$y[y_idx]
      y_max  <- y_max + (y_max*0.3)
      
    } else {
      errs_mssg <- glue::glue("Unsupported plot -type = '{plot_type}'")
      if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    print_tib( calc_sum, name = "Summary Data", vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #               Set Plotting Boundary Limits and Text Box::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    calc_str <- NULL
    for ( ii in c(1:base::nrow(calc_sum))) calc_str <- paste(
      calc_str, calc_sum[ ii , ] %>% paste( collapse = ", " ), sep = "\n" )
    
    calc_tib <- tibble::tibble(
      xlabel = as.double(x_min),
      ylabel = as.double(y_max),
      clabel = calc_str )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Generate Plot::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( plot_type == "dB" ) {
      gg <- ggplot( data = plot_data ) +
        geom_density( aes( x = dB, color=!!grp_sym, fill=!!grp_sym, alpha = alpha ) ) +
        geom_vline( xintercept = -dB_min, color="red",  linetype="dotted" ) +
        geom_vline( xintercept =  dB_min, color="blue", linetype="dotted" )
      
    } else if ( plot_type == "r2" ) {
      gg <- ggplot( data = plot_data, mapping = mapping ) +
        geom_point( size = size ) # + geom_density2d( alpha = alpha )
      
    } else if ( plot_type == "dg" ) {
      gg <- ggplot( data = plot_data ) +
        geom_density( aes( x = GGally::eval_data_col(plot_data, mapping$x), 
                           color=!!grp_sym, fill=!!grp_sym, alpha = alpha ) )
      #  geom_density( aes( x = dB, color=!!grp_sym, fill=!!grp_sym, alpha = alpha ) )
      
    } else {
      errs_mssg <- glue::glue("Unsupported plot -type = '{plot_type}'")
      if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Add Scaling and Text Boxes::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    gg <- gg +
      scale_x_continuous( breaks = round(seq(x_min, x_max, length.out=ticks),1), limits=c(x_min, x_max) ) +
      scale_y_continuous( breaks = round(seq(y_min, y_max, length.out=ticks),1), limits=c(y_min, y_max) ) +
      ggplot2::geom_label( data = calc_tib,
                           mapping = ggplot2::aes(x = xlabel, y = ylabel, label = clabel),
                           hjust = "left", vjust = "top",
                           size = wsize, fontface = "bold",
                           alpha = alpha_lab, family='mono',
                           inherit.aes = FALSE ) # do not inherit anything from the ...
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Plot::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # TBD:: Write Plot and Auxillary File...
    #   - Possiblily Change the format of this function to template_2 to save
    #     regenerating plots...
    #
    
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
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  gg
}

plot_beta_gg = function( betas,
                         detps = NULL,
                         manifest = NULL,
                         
                         ids_key = "Probe_ID",
                         prb_key = "Probe_Type",
                         inf_key = "Infinium_Design",
                         
                         dB_min = 0.2,
                         pval_min = 1.0,
                         beta_min = 0.3,
                         beta_max = 0.7,

                         sub_per = 5,
                         grp_key = "Group",
                         
                         top_tag,
                         sub_tag,
                         par_tag,
                         
                         dpi_val = 320,
                         format  = "pdf",
                         alpha   = 0.2,
                         
                         out_dir,
                         run_tag,
                         
                         reload     = 0,
                         reload_min = 2,
                         reload_pre = NULL,
                         
                         ret_data   = FALSE,
                         parallel   = FALSE,
                         write_out  = FALSE,
                         
                         vb=0, vt=3, tc=1, tt=NULL,
                         fun_tag='plot_beta_gg')
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
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1, tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}      ids_key = '{ids_key}'.{RET}"))
    cat(glue::glue("{mssg}      prb_key = '{prb_key}'.{RET}"))
    cat(glue::glue("{mssg}      inf_key = '{inf_key}'.{RET}"))
    cat(glue::glue("{mssg}       dB_min = '{dB_min}'.{RET}"))
    cat(glue::glue("{mssg}     pval_min = '{pval_min}'.{RET}"))
    cat(glue::glue("{mssg}     beta_min = '{beta_min}'.{RET}"))
    cat(glue::glue("{mssg}     beta_max = '{beta_max}'.{RET}"))
    cat(glue::glue("{mssg}      sub_per = '{sub_per}'.{RET}"))
    cat(glue::glue("{mssg}      dpi_val = '{dpi_val}'.{RET}"))
    cat(glue::glue("{mssg}       format = '{format}'.{RET}"))
    cat(glue::glue("{mssg}        alpha = '{alpha}'.{RET}"))
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
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ids_sym <- rlang::sym( ids_key )
  prb_sym <- rlang::sym( prb_key )
  inf_sym <- rlang::sym( inf_key )
  
  # [TBD]: Combine prb_key,inf_key as a new key below...
  # grp_key <- c(prb_key,inf_key)
  # grp_key <- inf_key
  
  grp_vec <- NULL
  grp_vec <- c( prb_key,inf_key )

  generic_ncol <- 3
  
  plot_path <- safe_mkdir( file.path( out_dir, "plots" ) )
  plot_file <- file.path( plot_path, paste(run_tag,"beta-r2",format, sep=".") )
  # plot_file <- file.path( plot_path, paste(par_tag,"r-squared",format, sep=".") )
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Process Inputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( is.null(detps) ) detps = matrix( data = 0.0, 
                                          nrow = base::nrow(betas), 
                                          ncol = base::ncol(betas) )
    
    # [TBD]: Apply p-value masking...
    # betas %>% tibble::column_to_row() %>% as.data.frame() %>% as.matrix
    # detps %>% tibble::column_to_row() %>% as.data.frame() %>% as.matrix
    # pvals_beta_mat <- clean_beta_mat
    # pvals_beta_mat[ which( clean_detp_mat > pval_min ) ] <- NA_real_
    
    beta_mat <- NULL
    beta_mat <- betas %>% 
      tibble::column_to_rownames( var = "Probe_ID" ) %>%
      as.data.frame() %>%
      as.matrix()
    
    detp_mat <- NULL
    detp_mat <- detps %>% 
      tibble::column_to_rownames( var = "Probe_ID" ) %>%
      as.data.frame() %>%
      as.matrix()
    
    beta_mat[ which( detp_mat > pval_min ) ] <- NA_real_
    betas <- beta_mat %>% as.data.frame() %>%
      tibble::rownames_to_column( var = "Probe_ID" ) %>%
      tibble::as_tibble()
    
    beta_all_tib <- NULL
    pval_all_tib <- NULL
    if ( !is.null(manifest) ) {
      if ( p0 ) cat(glue::glue("{mssg} Using Manifest...{RET}"))
      
      beta_all_tib <- manifest %>% 
        dplyr::select( dplyr::all_of( c(ids_key,prb_key,inf_key) ) ) %>%
        dplyr::inner_join( betas, by=c(ids_key) )

      pval_all_tib <- manifest %>% 
        dplyr::select( dplyr::all_of( c(ids_key,prb_key,inf_key) ) ) %>%
        dplyr::inner_join( detps, by=c(ids_key) )

    } else {
      if ( p0 ) cat(glue::glue("{mssg} NOT Using Manifest...{RET}"))
      
      beta_all_tib <- betas %>% 
        dplyr::mutate( !!prb_sym := !!ids_key %>% stringr::str_sub(1,2),
                       !!inf_sym := "I" )
      pval_all_tib <- detps %>%
        dplyr::mutate( !!prb_sym := !!ids_key %>% stringr::str_sub(1,2),
                       !!inf_sym := "I" )
    }
    
    # Add Grouping...
    # beta_all_tib <- beta_all_tib %>%
    #   dplyr::group_by( dplyr::across( dplyr::all_of( grp_vec ) ) )
    # pval_all_tib <- pval_all_tib %>%
    #   dplyr::group_by( dplyr::across( dplyr::all_of( grp_vec ) ) )
    
    # Sub-setting Plot Data::
    #
    tot_idx_cnt <- beta_all_tib %>% base::nrow()
    sub_idx_vec <- which( c(1:tot_idx_cnt) %% base::floor( 100 / sub_per ) == 1 )
    sub_idx_len <- base::length(sub_idx_vec)
    
    sub_idx_ksz <- number_as_commaK( sub_idx_len )
    tot_idx_ksz <- number_as_commaK( tot_idx_cnt )
    
    sub_val_per <- base::round( 100*sub_idx_len / tot_idx_cnt, 2 )
    sub_per_str <- glue::glue( "DPI={dpi_val}, ",
                               "Plot displays downsampled percent = {sub_val_per}% ",
                               "({sub_idx_ksz}/{tot_idx_ksz})")
    
    # Below line is from other code, not sure why the divide by 4 is needed...
    # rank_ncol <- ( base::ncol(beta_all_tib) - generic_ncol ) / 4
    # rank_ncol <- base::ncol(beta_all_tib) - generic_ncol
    # rank_cols <- paste0("V", c(1:rank_ncol) )
    # rank_beg <- base::ncol(beta_all_tib) - generic_ncol + 1
    
    rank_beg <- generic_ncol + 1
    rank_end <- base::ncol(beta_all_tib)
    rank_col <- names(beta_all_tib)[ c(rank_beg:rank_end) ]
    
    beta_sub_tib <- NULL
    beta_sub_tib <- beta_all_tib %>% 
      dplyr::filter( dplyr::row_number() %in% sub_idx_vec )
    
    beta_sub_cnt <- beta_sub_tib %>% base::nrow()
    beta_all_cnt <- beta_all_tib %>% base::nrow()
    
    if ( p0 ) cat(glue::glue("{mssg} Building Plot Params: ",
                             "sub={beta_sub_cnt}, ",
                             "all={beta_all_cnt}, ",
                             "rank_beg={rank_beg}, ",
                             "rank_end={rank_end}, ",
                             "...{RET}"))
    pairs_gg <- NULL
    pairs_gg <- GGally::ggpairs( 
      data = beta_sub_tib,
      mapping = ggplot2::aes( color = grp_key,
                              fill  = grp_key,
                              alpha = alpha ),
      columns = rank_col,
      
      upper = list(
        combo = "box_no_facet",
        continuous = GGally::wrap( beta_panel_gg,
                                   plot_type = "dB",
                                   pval_min = pval_min,
                                   beta_min = beta_min,
                                   beta_max = beta_max,
                                   dB_min   = dB_min,
                                   
                                   beta_data = beta_all_tib,
                                   pval_data = pval_all_tib,
                                   pval_cols = rank_col,
                                   full_plot = FALSE,
                                   grp_key = grp_key,
                                   grp_vec = grp_vec,
                                   
                                   out_dir = plot_path,
                                   run_tag = "Upper", 
                                   reload  = reload,
                                   reload_min = reload_min, 
                                   ret_data = FALSE,
                                   parallel = parallel,
                                   
                                   vb=vb,vt=vt+1,tc=tc+1,tt=tt ) ),
      # binwidth = c(5, 0.5),
      # high = "red",
      
      diag = list(
        continuous = GGally::wrap( beta_panel_gg,
                                   plot_type = "dg",
                                   pval_min = pval_min,
                                   beta_min = beta_min,
                                   beta_max = beta_max,
                                   dB_min   = dB_min,
                                   
                                   beta_data = beta_all_tib,
                                   pval_data = pval_all_tib,
                                   pval_cols = rank_col,
                                   full_plot = FALSE,
                                   grp_key = grp_key,
                                   grp_vec = grp_vec,
                                   
                                   out_dir = plot_path,
                                   run_tag = "Diag", 
                                   reload  = reload,
                                   reload_min = reload_min, 
                                   ret_data = FALSE,
                                   parallel = parallel,
                                   
                                   vb=vb,vt=vt+1,tc=tc+1,tt=tt ) ),
      
      lower = list(
        combo = "box_no_facet",
        continuous = GGally::wrap( beta_panel_gg,
                                   plot_type = "r2",
                                   pval_min = pval_min,
                                   beta_min = beta_min,
                                   beta_max = beta_max,
                                   dB_min   = dB_min,
                                   
                                   beta_data = beta_all_tib,
                                   pval_data = pval_all_tib,
                                   pval_cols = rank_col,
                                   full_plot = FALSE,
                                   grp_key = grp_key,
                                   grp_vec = grp_vec,
                                   
                                   out_dir = plot_path,
                                   run_tag = "Lower", 
                                   reload  = reload,
                                   reload_min = reload_min, 
                                   ret_data = FALSE,
                                   parallel = parallel,
                                   
                                   vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
    )
    pairs_gg <- pairs_gg + 
      labs( title=top_tag, subtitle=sub_tag, caption=sub_per_str )
    
    #
    # Write plots::
    #
    if ( p0 ) cat(glue::glue("{mssg} Plotting[{format}] = '{plot_file}'...{RET}"))
    ggplot2::ggsave( filename = plot_file, plot = pairs_gg, device = format, dpi = dpi_val )

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( write_out )
      out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
                             done = TRUE, write_spec = TRUE, append = FALSE, 
                             fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if ( ret_data ) {
    
    ret_dat$beta_sub_tib <- beta_sub_tib
    ret_dat$beta_all_tib <- beta_all_tib
    ret_dat$pval_all_tib <- pval_all_tib
    
    ret_dat$rank_beg     <- rank_beg
    ret_dat$rank_end     <- rank_end
    ret_dat$rank_col     <- rank_col
    ret_dat$sub_per_str  <- sub_per_str

    return( ret_dat )
  }
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                 Analyze Cutlist Percent Improvement::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

consolidate_matricies = function( betas, 
                                  detps, 
                                  clean=NULL,
                                  sum_stats = TRUE,
                                  all_stats = FALSE,
                                  vb=0, vt=6, tc=1, tt=NULL,
                                  fun_tag='consolidate_matricies')
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
  
  beta_col_cnt <- betas %>% base::ncol()
  beta_row_cnt <- betas %>% base::nrow()
  
  detp_col_cnt <- detps %>% base::ncol()
  detp_row_cnt <- detps %>% base::nrow()
  
  errs_mssg <- glue::glue("Matricies are not the same size: ",
                          "{beta_row_cnt} x {beta_col_cnt} != ",
                          "{detp_row_cnt} x {detp_col_cnt}.")
  if ( beta_col_cnt != detp_col_cnt ||
       beta_row_cnt != detp_row_cnt ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  ret_tib <- tibble::tibble(
    Probe_ID   = betas %>% rownames(),
    Sample_Cnt = beta_col_cnt,
    Probe_Cnt  = beta_row_cnt,
    
    Beta_Med   = betas %>% matrixStats::rowMedians( na.rm = TRUE ),
    Beta_Mad   = betas %>% matrixStats::rowMads( na.rm = TRUE ),
    
    Beta_Avg   = betas %>% matrixStats::rowMeans2( na.rm = TRUE ),
    Beta_Sds   = betas %>% matrixStats::rowSds( na.rm = TRUE ),
    
    Pval_Med   = detps %>% matrixStats::rowMedians( na.rm = TRUE ),
    Pval_Avg   = detps %>% matrixStats::rowMeans2( na.rm = TRUE )
  )
  
  # Add Extra Not Necessary Stats::
  if ( all_stats ) ret_tib <- ret_tib %>% 
    dplyr::mutate(
      Beta_Min  = betas %>% matrixStats::rowMins( na.rm = TRUE ),
      Beta_Max  = betas %>% matrixStats::rowMaxs( na.rm = TRUE ),
      
      Pval_Min  = detps %>% matrixStats::rowMins( na.rm = TRUE ),
      Pval_Max  = detps %>% matrixStats::rowMaxs( na.rm = TRUE ),
      
      Pval_Mad  = detps %>% matrixStats::rowMads( na.rm = TRUE ),
      Pval_Sds  = detps %>% matrixStats::rowSds( na.rm = TRUE )
    )
  
  # Add Summary Stats::
  if ( sum_stats ) ret_tib <- ret_tib %>% 
    dplyr::mutate(
      Pval_Fail_Cnt = betas %>% matrixStats::rowCounts( value = NA_real_ ),
      Pval_Fail_Per = round( 100*Pval_Fail_Cnt / Sample_Cnt, 3 ),
      
      Pval_Miss_Cnt = detps %>% matrixStats::rowCounts( value = NA_real_ ),
      Pval_Miss_Per = round( 100*Pval_Miss_Cnt / Sample_Cnt, 3 )
    )
  
  # Add Original Beta Failure Rates::
  if ( !is.null(clean) && 
       base::ncol(clean) == beta_col_cnt &&
       base::nrow(clean) == beta_row_cnt ) ret_tib <- ret_tib %>% 
    dplyr::mutate(
      Beta_Miss_Cnt = clean %>% matrixStats::rowCounts( value = NA_real_ ),
      Beta_Miss_Per = round( 100*Beta_Miss_Cnt / Sample_Cnt, 3 ),
    )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

consolidate_matricies_wrapper = function( betas,
                                          detps,
                                          manifest,
                                          
                                          pval_min = 1.0,
                                          dB_min   = 0.2,
                                          plot_con = FALSE,
                                          mix_all  = FALSE,
                                          
                                          sam_rgq = NULL,
                                          top_tag,
                                          sub_tag,
                                          par_tag,
                                          
                                          dpi_val = 320,
                                          
                                          out_dir,
                                          run_tag,
                                          
                                          reload     = 0,
                                          reload_min = 2,
                                          reload_pre = NULL,
                                          
                                          ret_data   = FALSE,
                                          parallel   = FALSE,
                                          
                                          vb=0, vt=3, tc=1, tt=NULL,
                                          fun_tag='consolidate_matricies_wrapper')
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
  p6  <- vb > vt + 6
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1, tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}     pval_min = '{pval_min}'.{RET}"))
    cat(glue::glue("{mssg}       dB_min = '{dB_min}'.{RET}"))
    cat(glue::glue("{mssg}      mix_all = '{mix_all}'.{RET}"))
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
  
  errs_mssg <- glue::glue("Beta Matrix cannot have zero length")
  if ( base::nrow(betas) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( p2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Clean Matricies::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    nan_tib <- NULL
    nan_tib <- tibble::tibble(
      Probe_ID  = detps %>% rownames(),
      Pval_Miss_Cnt = detps %>% matrixStats::rowCounts( value = NA_real_ ),
      Sample_Cnt  = detps %>% base::ncol()
    ) %>% dplyr::filter( ( Pval_Miss_Cnt == Sample_Cnt ) )
    valid_cgn_vec <- setdiff( rownames(detps), nan_tib$Probe_ID )
    clean_beta_mat  <- betas[ valid_cgn_vec, ]
    clean_detp_mat  <- detps[ valid_cgn_vec, ]
    
    if ( p2 ) {
      cat(glue::glue("{mssg} clean_beta_mat.dim={RET}"))
      clean_beta_mat %>% dim() %>% print()
      cat(glue::glue("{mssg} clean_detp_mat.dim={RET}"))
      clean_detp_mat %>% dim() %>% print()
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Filter by Detection P-value::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    pvals_beta_mat <- clean_beta_mat
    pvals_beta_mat[ which( clean_detp_mat > pval_min ) ] <- NA_real_
    
    if ( p2 ) {
      cat(glue::glue("{mssg} pvals_beta_mat.dim={RET}"))
      pvals_beta_mat %>% dim() %>% print()
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #          Calculate Sample Percet Passing Detection P-Value::
    #                        And Re-order Matrices::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    mat_col_cnt <- pvals_beta_mat %>% base::ncol()
    mat_row_cnt <- pvals_beta_mat %>% base::nrow()
    mid_col_idx <- base::ceiling( mat_col_cnt/2 )
    
    sam_tib <- NULL
    sam_tib <- tibble::tibble(
      Sentrix_Name  = clean_detp_mat %>% colnames(),
      Sample_Cnt  = mat_col_cnt,
      Probe_Cnt   = mat_row_cnt,
      
      Beta_Fail_Cnt = clean_beta_mat %>% matrixStats::colCounts( value = NA_real_ ),
      Pval_Fail_Cnt = pvals_beta_mat %>% matrixStats::colCounts( value = NA_real_ ),
      Beta_Fail_Per = round( 100*Beta_Fail_Cnt / Probe_Cnt, 3 ),
      Pval_Fail_Per = round( 100*Pval_Fail_Cnt / Probe_Cnt, 3 )
    ) %>% 
      dplyr::arrange( Pval_Fail_Per )
    #  dplyr::arrange( -Pval_Fail_Per )
    
    if ( !is.null( sam_rgq ) ) {
      sam_tib <- sam_tib %>% 
        dplyr::inner_join( sam_rgq, by=c("Sentrix_Name") ) %>% 
        dplyr::arrange( Registration_Score ) %>%
        dplyr::select( -Pval_Fail_Per ) %>%
        dplyr::rename( Pval_Fail_Per = Registration_Score )
      if ( p1 ) cat(glue::glue("{mssg} Sample Swap Done!{RET2}"))
      if ( p2 ) print(sam_tib)
    }
    
    pvals_beta_mat <- pvals_beta_mat[ , sam_tib$Sentrix_Name ]
    clean_detp_mat <- clean_detp_mat[ , sam_tib$Sentrix_Name ]
    
    #
    # TBD:: All provide Samples against eachother::
    #
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Partition Samples::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    rank_cols <- c( "Top", "All", "Bot" )
    generic_cols <- c("Probe_ID","Probe_Type","Infinium","Probe_Group_Mask","BP_Group","Probe_Group_BP")
    generic_ncol <- length( generic_cols )
    
    sentrix_names <- list()
    sentrix_names[["con"]] <- list()
    sentrix_names[["mix"]] <- list()
    
    sentrix_lists <- list()
    sentrix_lists[["con"]] <- list()
    sentrix_lists[["mix"]] <- list()
    
    sentrix_lists[["con"]][["All"]] <- sam_tib$Sentrix_Name
    sentrix_lists[["mix"]][["All"]] <- sam_tib$Sentrix_Name
    if ( mat_col_cnt > 2 ) sentrix_lists[["mix"]][["All"]] <- sam_tib$Sentrix_Name[ c(1,mid_col_idx,mat_col_cnt) ]
    
    # print(sentrix_lists[["con"]][["All"]])
    # print(sentrix_lists[["mix"]][["All"]])
    
    if ( mat_col_cnt > 5 ) {
      sentrix_lists[["con"]] <- list()
      sentrix_lists[["con"]][["Top"]] <- sam_tib$Sentrix_Name[ 1:mid_col_idx ]
      sentrix_lists[["con"]][["All"]] <- sam_tib$Sentrix_Name
      sentrix_lists[["con"]][["Bot"]] <- sam_tib$Sentrix_Name[ (mid_col_idx+1):mat_col_cnt ]
      
      sentrix_lists[["mix"]] <- list()
      if ( mat_col_cnt == 6 ) {
        sentrix_lists[["mix"]][["Top"]] <- sam_tib$Sentrix_Name[ c(1:3) ]
        sentrix_lists[["mix"]][["All"]] <- sam_tib$Sentrix_Name
        sentrix_lists[["mix"]][["Bot"]] <- sam_tib$Sentrix_Name[ c(3:5) ]
      }
      if ( mat_col_cnt > 6 ) {
        qrt_idx <- base::ceiling(mid_col_idx/2)
        mix_vec <- c( 1,qrt_idx,mid_col_idx,
                      mid_col_idx+1,mid_col_idx+qrt_idx+1,mat_col_cnt ) %>% 
          unique() %>% sort()
        
        top_vec <- mix_vec[ c(1:3) ]
        all_vec <- c(1,mid_col_idx, mat_col_cnt)
        bot_vec <- mix_vec[ c(4:length(mix_vec)) ]
        
        sentrix_lists[["mix"]][["Top"]] <- sam_tib$Sentrix_Name[ top_vec ]
        sentrix_lists[["mix"]][["All"]] <- sam_tib$Sentrix_Name[ all_vec ]
        sentrix_lists[["mix"]][["Bot"]] <- sam_tib$Sentrix_Name[ bot_vec ]
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Consolidate Matrices::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( p0 ) cat(glue::glue("{mssg} Consolidate matricies...{RET2}"))
    
    cond_tab <- NULL
    cond_tab <- sentrix_lists[["con"]] %>% lapply(
      function(x,b,d,c=NULL,s=TRUE,a=FALSE, vb=0,vt=6,tc=1, tt=NULL ) {
        consolidate_matricies( betas=b[,x],detps=d[,x],clean=c[,x],
                               sum_stats=s,all_stats=a,
                               vb=vb,vt=vt,tc=tc ) },
      b = pvals_beta_mat, d = clean_detp_mat, c = clean_beta_mat, 
      s = TRUE, a = FALSE, vb=vb,vt=vt+2,tc=tc+1 ) %>% 
      dplyr::bind_rows( .id = "Partition" )
    
    #
    # ret_dat$pvals_beta_mat[ , c( ret_dat$sentrix_lists$mix$Top,ret_dat$sentrix_lists$mix$Bot )  ] %>% dim()
    #
    # tmp_lst <- c( sentrix_lists[["mix"]]$Top, sentrix_lists[["mix"]]$Bot ) %>% as.list()
    # names(tmp_lst) <- c( "T1", "T2", "T3",  "B1", "B2", "B3" )
    # 
    # cond_tab <- NULL
    # cond_tab <- tmp_lst %>% lapply(
    #   function(x,b,d,c=NULL,s=TRUE,a=FALSE, vb=0,vt=6,tc=1,tt=NULL ) {
    #     consolidate_matricies( betas=b[,x],detps=d[,x],clean=c[,x],
    #                            sum_stats=s,all_stats=a,
    #                            vb=vb+10,vt=vt,tc=tc ) },
    #   b = pvals_beta_mat, d = clean_detp_mat, c = clean_beta_mat, 
    #   s = TRUE, a = FALSE, vb=vb,vt=vt+2,tc=tc+1 ) %>% 
    #   dplyr::bind_rows( .id = "Partition" )
    
    #
    # Return To Normal Workflow::
    #
    
    ret_tib <- cond_tab
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Extract Real Matrices::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( p0 ) cat(glue::glue("{mssg} Extracting real(mix) matricies...{RET2}"))
    # cat("\n\n\n\n\n\n\n-----------------------------------------\n\n\n\n\n\n")
    
    mix_all_betas <- list()
    mix_all_pvals <- list()
    mix_opp_datas <- list()
    
    # mix_all <- TRUE
    # mix_all <- FALSE
    
    if (mix_all) {
      # print(sentrix_lists)
      for ( mix_key in names(sentrix_lists) ) {
        if ( mix_key != "mix" ) next
        if ( p1 ) cat(glue::glue("{mssg}{TAB} Extracting matrix[{mix_key}]...{RET}"))
        for ( rank_key in names(sentrix_lists[[mix_key]]) ) {
          if ( p2 ) cat(glue::glue("{mssg}{TAB2} Extracting matrix[{mix_key}][{rank_key}]...{RET}"))
          
          cur_rank_cols <- paste0( 
            "p", sam_tib %>% 
              dplyr::filter( Sentrix_Name %in% sentrix_lists[[mix_key]][[rank_key]] ) %>% 
              dplyr::pull( Pval_Fail_Per ) %>% as.integer()
          )
          
          b_tib <- pvals_beta_mat[ , sentrix_lists[[mix_key]][[rank_key]] ] %>% 
            as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% 
            tibble::as_tibble()
          print_tib( b_tib, name = "b_tib", vb=vb,vt=vt+4,tc=tc+1 )
          
          p_tib <- clean_detp_mat[ , sentrix_lists[[mix_key]][[rank_key]] ] %>% 
            as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% 
            tibble::as_tibble()
          print_tib( p_tib, name = "p_tib", vb=vb,vt=vt+4,tc=tc+1 )
          
          mix_opp_datas[[mix_key]][[rank_key]] <- dplyr::inner_join(
            b_tib, p_tib, 
            by=c("Probe_ID"),
            suffix=c("_beta","_pval")
          )
          
          # mix_all_betas[[mix_key]][[rank_key]] <- manifest %>%
          #   dplyr::right_join( b_tib, by=c("Probe_ID") ) %>% 
          #   purrr::set_names( c(generic_cols, cur_rank_cols ) )
          # print_tib( mix_all_betas[[mix_key]][[rank_key]], name = "mix_all_betas", vb=vb,vt=vt+4,tc=tc+1 )
          # 
          # mix_all_pvals[[mix_key]][[rank_key]] <- manifest %>%
          #   dplyr::right_join( p_tib, by=c("Probe_ID") ) %>% 
          #   purrr::set_names( c(generic_cols, cur_rank_cols ) )
          # print_tib( mix_all_pvals[[mix_key]][[rank_key]], name = "mix_all_pvals", vb=vb,vt=vt+4,tc=tc+1 )
          
          if ( p2 ) cat(glue::glue("{mssg}{TAB2} Done matrix[{mix_key}][{rank_key}].{RET2}"))
        }
        if ( p1 ) cat(glue::glue("{mssg}{TAB} Done matrix[{mix_key}].{RET2}"))
      }
    }
    if ( p0 ) cat(glue::glue("{mssg} Done. Extracting real(mix) matricies.{RET2}"))
    
    #
    # For Mix we simply need to build the split tables
    #
    # mix_opp_betas <- ewas_class_man_tib %>% dplyr::right_join( dplyr::inner_join( ret_dat$mix_opp_datas$mix$Top, ret_dat$mix_opp_datas$mix$Bot, by=c("Probe_ID") ), by=c("Probe_ID") ) %>% names()
    # dplyr::inner_join( ret_dat$mix_opp_datas$mix$Top, ret_dat$mix_opp_datas$mix$Bot, by=c("Probe_ID") ) %>% tidyr::pivot_longer( cols = c(-Probe_ID), names_to = c("Sentrix_ID","Sentrix_Pos","Data_Type"), names_sep = "_", values_to = "Metric" ) %>% tidyr::unite( Sentrix_Name, Sentrix_ID, Sentrix_Pos, sep="_" ) %>% tidyr::pivot_wider( id_cols = c(Probe_ID, Sentrix_Name ), names_from = c( Data_Type ), values_from = c( Metric ) ) %>% dplyr::left_join( ret_dat$sam_tib, by=c("Sentrix_Name") )
    
    #
    # This build a lot of it...
    #
    # dplyr::inner_join( mix_opp_datas$mix$Top, mix_opp_datas$mix$Bot, by=c("Probe_ID") ) %>% tidyr::pivot_longer( cols = c(-Probe_ID), names_to = c("Sentrix_ID","Sentrix_Pos","Data_Type"), names_sep = "_", values_to = "Metric" ) %>% tidyr::unite( Sentrix_Name, Sentrix_ID, Sentrix_Pos, sep="_" ) %>% tidyr::pivot_wider( id_cols = c(Probe_ID, Sentrix_Name ), names_from = c( Data_Type ), values_from = c( Metric ) ) %>% dplyr::left_join( sam_tib, by=c("Sentrix_Name") )
    #
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Plot Matrices:: Consolidated
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cond_tib <- NULL
    cond_tib <- cond_tab %>% 
      dplyr::select( 
        Probe_ID, Partition, Beta_Med, Beta_Mad, Pval_Med, Pval_Fail_Per ) %>% 
      tidyr::pivot_wider( 
        id_cols = c(Probe_ID), 
        names_from = c(Partition), 
        values_from = c(Beta_Med,Beta_Mad,Pval_Med,Pval_Fail_Per) )
    print_tib( cond_tib, name = "cond_tib", vb=vb,vt=vt+4,tc=tc+1 )
    
    #
    # Plotting:: ggpairs()
    #
    if ( p6 ) {
      cat(glue::glue("{mssg}{TAB2} ----- ----- ----- ----- ----- DONE ----- ----- ----- ----- -----{RET2}"))
      cat(glue::glue("{mssg}{TAB2} cond_tib={RET}"))
      print(cond_tib)
      cat(glue::glue("{mssg}{TAB2} ----- ----- ----- ----- ----- DONE ----- ----- ----- ----- -----{RET2}"))
    }
    
    cond_all_tib <- NULL
    if ( !is.null(manifest) ) {
      cond_all_tib <- manifest %>% 
        dplyr::right_join( cond_tib, by=c("Probe_ID") ) %>%
        dplyr::arrange( Pval_Med_All )
    } else {
      die( glue::glue("{mssg} FAILED! Manifest must be provided!{RET2}") )
      return(NULL)
      
      # cond_all_tib <- cond_tib %>%
      #   dplyr::mutate( Pval_Med = Beta_Med_All ) %>%
      #   dplyr::arrange( Pval_Med_All )
    }
    print_tib( cond_all_tib, name = "cond_all_tib", vb=vb,vt=vt+4,tc=tc+1 )
    
    pairs_gg <- NULL
    # plot_con <- FALSE
    # plot_con <- TRUE
    
    if ( plot_con ) {
      
      #
      # Sub-setting Plot Data::
      #
      sub_per_val <- 5
      tot_idx_cnt <- cond_all_tib %>% base::nrow()
      sub_idx_vec <- which( c(1:tot_idx_cnt) %% base::floor( 100 / sub_per_val ) == 1 )
      sub_idx_len <- base::length(sub_idx_vec)
      
      sub_idx_ksz <- number_as_commaK( sub_idx_len )
      tot_idx_ksz <- number_as_commaK( tot_idx_cnt )
      
      sub_val_per <- base::round( 100*sub_idx_len / tot_idx_cnt, 2 )
      sub_per_str <- glue::glue( "DPI={dpi_val}, ",
                                 "Plot displays downsampled percent = {sub_val_per}% ",
                                 "({sub_idx_ksz}/{tot_idx_ksz})")
      
      # [Skip]: Facet Grid Rows = Probe_Group_BP[1-7], Other
      cond_grp_vec <- c("Probe_Group_BP","Infinium","Probe_Group_Mask")
      cond_grp_vec <- c("Probe_Type","Infinium" )
      
      cond_sub_tib <- NULL
      if ( mat_col_cnt > 5 ) {
        cond_sub_tib <- cond_all_tib %>% 
          dplyr::filter( !( is.na(Beta_Med_All) & is.na( Beta_Med_Top) & is.na( Beta_Med_Bot) ) ) %>%
          dplyr::filter( dplyr::row_number() %in% sub_idx_vec )
      } else {
        cond_sub_tib <- cond_all_tib %>% 
          dplyr::filter( !( is.na(Beta_Med_All) ) ) %>%
          dplyr::filter( dplyr::row_number() %in% sub_idx_vec )
        
        rank_ncol <- ( base::ncol(cond_sub_tib) - generic_ncol ) / 4
        rank_cols <- paste0("V", c(1:rank_ncol) )
      }
      print_tib( cond_sub_tib, name = "cond_sub_tib", vb=vb,vt=vt+4,tc=tc+1 )
      
      beta_dat_cols <- cond_all_tib %>% dplyr::select( dplyr::starts_with("Beta_Med") ) %>% names()
      beta_var_cols <- cond_all_tib %>% dplyr::select( dplyr::starts_with("Beta_Mad") ) %>% names()
      pval_dat_cols <- cond_all_tib %>% dplyr::select( dplyr::starts_with("Pval_Med") ) %>% names()
      samp_per_cols <- cond_all_tib %>% dplyr::select( dplyr::starts_with("Pval_Fail_Per") ) %>% names()
      
      if ( p6 ) {
        cat(glue::glue("{RET2}{mssg} beta_dat_cols:{RET}"))
        print( beta_dat_cols )
        cat(glue::glue("{RET2}{mssg} samp_per_cols:{RET}"))
        print( samp_per_cols )
        cat(glue::glue("{RET2}{mssg} generic_cols:{RET}"))
        print( generic_cols )
        cat(glue::glue("{RET2}{mssg} beta_dat_cols:{RET}"))
        print( beta_dat_cols )
        cat(glue::glue("{RET2}{mssg} rank_cols:{RET}"))
        print( rank_cols )
        cat(glue::glue("{RET2}{mssg} cond_sub_tib:{RET}"))
        print_tib( cond_sub_tib, name = "cond_sub_tib", vb=vb,vt=vt+4,tc=tc+1 )
        
        cond_sub_tib %>% 
          dplyr::select( dplyr::any_of( c(generic_cols, beta_dat_cols) ) ) %>%
          print()
        # c(generic_cols, rank_cols ) %>% print()
      }
      
      # cond_all_tib %>% dplyr::select( dplyr::all_of( generic_cols, beta_dat_cols ) )
      
      if ( p6 ) {
        cat(glue::glue("{mssg}{TAB2} generic_cols::{RET}"))
        print( generic_cols )
        cat(glue::glue("{mssg}{TAB2} beta_dat_cols::{RET}"))
        print( beta_dat_cols )
        cat(glue::glue("{mssg}{TAB2} ----- ----- ----- ----- ----- DONE ----- ----- ----- ----- -----{RET2}"))
      }
      
      beta_sub_tib <- NULL
      beta_sub_tib <- cond_sub_tib %>% 
        dplyr::select( dplyr::any_of( c(generic_cols, beta_dat_cols) ) ) %>% 
        purrr::set_names( c(generic_cols, rank_cols ) )
      print_tib( beta_sub_tib, name = "beta_sub_tib", vb=vb,vt=vt+4,tc=tc+1 )
      
      beta_all_tib <- NULL
      beta_all_tib <- cond_all_tib %>% dplyr::select( dplyr::all_of( c(generic_cols, beta_dat_cols) ) ) %>% 
        purrr::set_names( c(generic_cols, rank_cols ) )
      pval_all_tib <- NULL
      pval_all_tib <- cond_all_tib %>% dplyr::select( dplyr::all_of( c(generic_cols, pval_dat_cols) ) ) %>% 
        purrr::set_names( c(generic_cols, rank_cols ) )
      
      # NOTE: Not sure what these are used for...
      #
      # vars_all_tib <- NULL
      # vars_all_tib <- cond_all_tib %>% dplyr::select( dplyr::all_of( c(generic_cols, beta_var_cols) ) ) %>% 
      #   purrr::set_names( c(generic_cols, rank_cols ) )
      # samp_all_tib <- NULL
      # samp_all_tib <- cond_all_tib %>% dplyr::select( dplyr::all_of( c(generic_cols, samp_per_cols) ) ) %>% 
      #   purrr::set_names( c(generic_cols, rank_cols ) )
      
      grp_key <- "Group"
      plot_dir <- file.path( out_dir, "plots" )
      
      if ( p6 ) cat(glue::glue("{mssg}{TAB2} Setting plot_dir='{plot_dir}'..{RET2}"))
      
      pairs_gg <- NULL
      pairs_gg <- GGally::ggpairs( 
        data = beta_sub_tib,
        mapping = ggplot2::aes( color = grp_key,
                                fill  = grp_key,
                                alpha = 0.2 ),
        columns = rank_cols,
        
        upper = list(
          combo = "box_no_facet",
          continuous = GGally::wrap( beta_panel_gg,
                                     plot_type = "dB",
                                     pval_min = pval_min,
                                     beta_min = 0.3,
                                     beta_max = 0.7,
                                     dB_min   = dB_min,
                                     
                                     beta_data = beta_all_tib,
                                     pval_data = pval_all_tib,
                                     pval_cols = rank_cols,
                                     full_plot = FALSE,
                                     grp_key = grp_key,
                                     grp_vec = cond_grp_vec,
                                     
                                     out_dir = plot_dir,
                                     run_tag = "Upper", 
                                     reload = reload,
                                     reload_min = reload_min, 
                                     ret_data = FALSE,
                                     parallel = parallel,
                                     
                                     vb=vb,vt=vt+1,tc=tc+1,tt=tt ) ),
        # binwidth = c(5, 0.5),
        # high = "red",
        
        diag = list(
          continuous = GGally::wrap( beta_panel_gg,
                                     plot_type = "dg",
                                     pval_min = pval_min,
                                     beta_min = 0.3,
                                     beta_max = 0.7,
                                     dB_min   = dB_min,
                                     
                                     beta_data = beta_all_tib,
                                     pval_data = pval_all_tib,
                                     pval_cols = rank_cols,
                                     full_plot = FALSE,
                                     grp_key = grp_key,
                                     grp_vec = cond_grp_vec,
                                     
                                     out_dir = plot_dir,
                                     run_tag = "Diag", 
                                     reload = reload,
                                     reload_min = reload_min, 
                                     ret_data = FALSE,
                                     parallel = parallel,
                                     
                                     vb=vb,vt=vt+1,tc=tc+1,tt=tt ) ),
        
        lower = list(
          combo = "box_no_facet",
          continuous = GGally::wrap( beta_panel_gg,
                                     plot_type = "r2",
                                     pval_min = pval_min,
                                     beta_min = 0.3,
                                     beta_max = 0.7,
                                     dB_min   = dB_min,
                                     
                                     beta_data = beta_all_tib,
                                     pval_data = pval_all_tib,
                                     pval_cols = rank_cols,
                                     full_plot = FALSE,
                                     grp_key = grp_key,
                                     grp_vec = cond_grp_vec,
                                     
                                     out_dir = plot_dir,
                                     run_tag = "Lower", 
                                     reload = reload,
                                     reload_min = reload_min, 
                                     ret_data = FALSE,
                                     parallel = parallel,
                                     
                                     vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
      )
      pairs_gg <- pairs_gg + 
        labs( title=top_tag, subtitle=sub_tag, caption=sub_per_str )
      
      #
      # Write plots::
      #
      
      # Small function to display plots only if it's interactive
      #  p_ <- GGally::print_if_interactive
      
      tmp_dir <- safe_mkdir(plot_dir)
      
      # plot_types <- c( "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" )
      # plot_types <- c( "pdf", "jpeg", "tiff", "png", "bmp", "svg", "wmf" )
      plot_types <- c( "pdf" )
      
      for ( plot_type in plot_types ) {
        plot_file <- file.path( plot_dir, paste(par_tag,"consolidation-split",plot_type, sep=".") )
        if ( p0 ) cat(glue::glue("{mssg} Plotting[{plot_type}] = 'plot_file'...{RET2}"))
        
        # try( ggplot2::ggsave( filename = plot_file, plot = pairs_gg, device = plot_type, dpi = dpi_val ) )
        ggplot2::ggsave( filename = plot_file, plot = pairs_gg, device = plot_type, dpi = dpi_val )
      }
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
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if ( ret_data ) {
    
    ret_dat$sam_tib <- sam_tib
    ret_dat$pvals_beta_mat <- pvals_beta_mat
    ret_dat$clean_detp_mat <- clean_detp_mat
    ret_dat$cond_tab       <- cond_tab
    ret_dat$cond_tib       <- cond_tib
    ret_dat$cond_all_tib   <- cond_all_tib
    ret_dat$sentrix_lists  <- sentrix_lists
    ret_dat$mix_all_betas  <- mix_all_betas
    ret_dat$mix_all_pvals  <- mix_all_pvals
    ret_dat$mix_opp_datas  <- mix_opp_datas
    
    return( ret_dat )
  }
  
  cond_tab
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Replicate Analysis::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This is being called incorrectly::
#
replicate_analysis = function( sdfs, 
                               steps,
                               mask_ids,
                               sample,
                               min_dB   = 0.20,
                               min_pval = 0.20,
                               sam_pval = 98,
                               
                               ran_lapply = FALSE,
                               out_dir,
                               run_tag,
                               
                               reload     = 0,
                               reload_min = 2,
                               reload_pre = NULL,
                               
                               ret_data   = FALSE,
                               parallel   = FALSE,
                               
                               vb=0, vt=3, tc=1, tt=NULL,
                               fun_tag='replicate_analysis')
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
  p8  <- vb > vt + 8
  
  ret_dat <- NULL
  ret_tib <- NULL
  beta_r2_tib <- NULL
  samp_stats_tib <- NULL
  
  ret_cnt <- 0
  ret_tib <- NULL
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sam_csv <- file.path( out_dir, paste(out_tag, 'sample_sheet.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1, tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}       sample = '{sample}'.{RET}"))
    cat(glue::glue("{mssg}       min_dB = '{min_dB}'.{RET}"))
    cat(glue::glue("{mssg}     min_pval = '{min_pval}'.{RET}"))
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
    cat(glue::glue("{mssg}      sam_csv = '{sam_csv}'.{RET}"))
    cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  step_vec <- NULL
  step_vec <- stringr::str_split( string = steps, 
                                  pattern = "", 
                                  n = stringr::str_length(steps), 
                                  simplify = FALSE ) %>% 
    base::unlist() %>% as.vector()
  step_cnt = step_vec %>% length()
  
  errs_mssg <- glue::glue("Failed to parse workflow steps='{steps}' not defined")
  if ( base::length( step_vec) == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  if ( p1 ) cat(glue::glue("{mssg} Total Number of Workflow Steps={step_cnt}.{RET}"))
  
  #
  # Validate lapply has been previously ran::
  #  TBD:: Understand why we need this stupid trick. It happended after the 
  #  R 4.2.1 update...
  #
  # if ( !ran_lapply ) ran_lapply = validate_lapply( sdfs = sdfs[[1]] )
  # errs_mssg <- glue::glue("Failed to run lapply (silly issue) ran_lapply={ran_lapply}!!!")
  # if ( !ran_lapply ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  # if ( p1 ) cat(glue::glue("{mssg} Ran Lapply Completed={ran_lapply}.{RET}"))
  
  #
  # Run replicate Analysis::
  #
  call_tib <- tibble::tibble()
  call_tib <- sdfs %>% 
    lapply( function(x) { is.na( sesame::sesameQC_calcStats( 
      sdf = sesame::inferInfiniumIChannel( sdf = x, switch_failed = FALSE), 
      funs = "dyeBias" )@stat$RGdistort ) }) %>% 
    dplyr::bind_rows()
  
  call_cnt = call_tib %>% base::ncol()
  errs_mssg <- glue::glue("Failed to build passing samples call_cnt='{call_cnt}'!!!")
  if ( is.null(call_tib) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(call_tib)
  print_tib( x = call_tib, fun_tag = fun_tag, name = "beta_sel", 
             vb=vb,vt=vt+6,tc=tc )
  if ( p1 ) cat(glue::glue("{mssg} Callable Samples Count={call_cnt}.{RET}"))
  
  beta_sel <- tibble::tibble()
  beta_sel <- call_tib[ which( call_tib== FALSE) ]
  
  beta_nrow = beta_sel %>% base::nrow()
  beta_ncol = beta_sel %>% base::ncol()
  errs_mssg <- glue::glue("Failed to build select beta values ",
                          "beta_nrow='{beta_nrow}', beta_ncol='{beta_ncol}'!!!")
  if ( is.null(beta_sel) || beta_nrow==0 || beta_ncol==0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(beta_sel)
  print_tib( x = beta_sel, fun_tag = fun_tag, name = "beta_sel", 
             vb=vb,vt=vt+6,tc=tc )
  if ( p1 ) cat(glue::glue("{mssg} Beta Call Samples r x c = {beta_nrow} x {beta_ncol}.{RET}"))
  
  beta_raw <- tibble::tibble()
  if ( parallel ) {
    beta_raw <- BiocGenerics::do.call( 
      BiocGenerics::cbind, 
      parallel::mclapply( sdfs[ beta_sel %>% colnames() ],
                          mutate_sdf_linear, 
                          steps = steps, 
                          mask_ids = mask_ids,
                          
                          min_pval = min_pval, 
                          off_set = 15,
                          ran_lapply = ran_lapply,
                          
                          out_dir = opt$out_path,
                          run_tag = opt$run_name,
                          
                          reload     = opt$reload,
                          reload_min = 3,
                          reload_pre = NULL,
                          
                          ret_data   = FALSE,
                          parallel   = opt$parallel,
                          
                          vb=vb,vt=vt+4,tc=tc ) )
  } else {
    beta_raw <- BiocGenerics::do.call( 
      BiocGenerics::cbind, 
      base::lapply( sdfs[ beta_sel %>% colnames() ],
                    mutate_sdf_linear, 
                    steps = steps, 
                    mask_ids = mask_ids,
                    
                    min_pval = min_pval, 
                    off_set = 15,
                    ran_lapply = ran_lapply,
                    
                    out_dir = opt$out_path,
                    run_tag = opt$run_name,
                    
                    reload     = opt$reload,
                    reload_min = 3,
                    reload_pre = NULL,
                    
                    ret_data   = FALSE,
                    parallel   = opt$parallel,
                    
                    vb=vb,vt=vt+4,tc=tc ) )
  }
  print_tib( x = beta_raw, fun_tag = fun_tag, name = "beta_raw", 
             vb=vb,vt=vt+6,tc=tc )
  beta_nrow <- beta_raw %>% base::nrow()
  beta_ncol <- beta_raw %>% base::ncol()
  if ( p1 ) cat(glue::glue("{mssg} Beta Raw Samples Count=r{beta_nrow} x c{beta_ncol}.{RET}"))
  
  #
  # Screen Samples by Percent Passing Detection P-value::
  #
  sentrix_pass_tib <- tibble::tibble( 
    Sentrix_Name = beta_raw %>% colnames(), 
    Tot_Cnt = beta_raw %>% base::nrow(), 
    Nan_Cnt = beta_raw %>% matrixStats::colCounts( value=NA ), 
    Pass_Perc = 100 - round( 100 * Nan_Cnt/Tot_Cnt, 2 )
  ) %>% dplyr::filter( Pass_Perc > sam_pval )
  print_tib( x = sentrix_pass_tib, fun_tag = fun_tag, name = "sentrix_pass_tib", 
             vb=vb,vt=vt+6,tc=tc )
  
  if ( base::nrow(sentrix_pass_tib) != 0 ) {
    
    beta_raw <- beta_raw[ , sentrix_pass_tib$Sentrix_Name ]
    # return( beta_raw )
    print_tib( x = beta_raw, fun_tag = fun_tag, name = "beta_raw_filtered", 
               vb=vb,vt=vt+6,tc=tc )
    beta_nrow <- beta_raw %>% base::nrow()
    beta_ncol <- beta_raw %>% base::ncol()
    if ( p1 ) cat(glue::glue("{mssg} Beta Sample Filtered Samples Count={beta_nrow} x {beta_ncol}.{RET}"))
    
    beta_r2_tib <- NULL
    beta_r2_tib <- beta_raw  %>% 
      cor( use = "pairwise.complete.obs", method = "pearson" ) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column( var = "Sentrix_Name_A" ) %>% 
      tidyr::pivot_longer( cols = c(-Sentrix_Name_A), 
                           names_to = "Sentrix_Name_B", 
                           values_to = "r2" ) %>%
      dplyr::filter( Sentrix_Name_A != Sentrix_Name_B ) %>% 
      dplyr::group_by( Sentrix_Name_A ) %>% 
      dplyr::summarise( dplyr::across( 
        where(is.numeric), 
        list( min = min, max = max, avg = mean, sds = sd, 
              med = median, mad = mad ) ), 
        .groups = "drop" ) %>%
      dplyr::rename( Sentrix_Name = Sentrix_Name_A )
    beta_r2_tib <- sentrix_pass_tib %>% 
      dplyr::inner_join( beta_r2_tib, by=c("Sentrix_Name") )
    
    print_tib( x = beta_r2_tib, fun_tag = fun_tag, name = "beta_r2_tib", 
               vb=vb,vt=vt+6,tc=tc )
    
    samp_stats_tib <- NULL
    samp_stats_tib <- 
      sample_performance_mat( beta_mat_r = beta_raw, 
                              pval_mat_r = beta_raw, 
                              min_dB = min_dB, 
                              vb=vb,vt=vt+10,tc=tc )
    samp_stats_tib <- dplyr::left_join( 
      beta_r2_tib, 
      samp_stats_tib %>% 
        dplyr::filter( Sentrix_Name_A != Sentrix_Name_B ) %>% 
        dplyr::group_by( Sentrix_Name_A,Min_Delta_Beta,Pre_Pval_Count) %>% 
        dplyr::summarise( dplyr::across( 
          where(is.numeric), 
          list( min = min, max = max, avg = mean, med = median ) ), 
          .groups = "drop" ) %>%
        dplyr::rename( Sentrix_Name = Sentrix_Name_A ) %>% 
        purrr::set_names( names(.) %>% stringr::str_replace("r2_","r2B_") ) %>%
        dplyr::select( Sentrix_Name,Min_Delta_Beta, 
                       Pval_Pass_Percent_avg,
                       Pval_Pass_Percent_med,
                       dB_Pass_Percent_med,r2B_med,r2B_med ),
      by=c("Sentrix_Name") ) %>%
      dplyr::mutate( 
        Min_dB := !!min_dB,
        Sample_Pval_Min := !!sam_pval,
        Sample_Name := !!sample
      ) %>%
      dplyr::select( Min_dB,Sample_Pval_Min,Sample_Name, dplyr::everything() )
    
    print_tib( x = samp_stats_tib, fun_tag = fun_tag, name = "samp_stats_tib", 
               vb=vb,vt=vt+6,tc=tc )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Write Sample Sheet::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = samp_stats_tib, file = sam_csv, type = "csv",
                           done = TRUE, write_spec = TRUE, append = FALSE,
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Generate Probe Variation Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    sum_tib <- NULL
    sum_tib <- loci_variation_mat(
      beta_mat_r = beta_raw,
      pval_mat_r = beta_raw,
      min_dB   = min_dB,
      min_pval = min_pval,
      vb=vb, vt=vt+1+10, tc=tc ) %>% 
      tibble::as_tibble()
    # dplyr::mutate( Sample_Name := !!sample ) %>%
    # dplyr::select( Sample_Name, dplyr::everything() )
    
    print_tib( x = sum_tib, fun_tag = fun_tag, name = "sum_tib", 
               vb=vb,vt=vt+6,tc=tc )
    
    ret_tib <- NULL
    ret_tib <- ret_tib %>% dplyr::bind_rows( sum_tib ) %>%
      dplyr::mutate( 
        Min_dB := !!min_dB,
        Sample_Pval_Min := !!sam_pval,
        Sample_Name := !!sample
      )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  }
  etime <- round( as.double(ftime[3]), 2 )
  
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if ( ret_data ) {
    ret_dat <- NULL
    ret_dat$beta_r2_tib <- beta_r2_tib
    ret_dat$samp_stats_tib <- samp_stats_tib
    ret_dat$ret_tib <- ret_tib
    
    return( ret_dat )
  }
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                 Analyze Cutlist Percent Improvement::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

compare_samples = function( tib,
                            prbs_key,
                            datA_key,
                            datB_key,
                            dB_min,
                            
                            vb=0, vt=6, tc=1, tt=NULL,
                            fun_tag='compare_samples')
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
  
  datA_sym = rlang::sym(datA_key)
  datB_sym = rlang::sym(datB_key)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             General Stats::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  tot_cnt <- tib %>% base::nrow()
  pas_cnt <- tib %>%
    dplyr::filter( !is.na(datA_key) & !is.na(datB_key) ) %>%
    base::nrow()
  pas_per <- round( 100 * pas_cnt / tot_cnt, 3 )
  
  if ( vb >= vt+3 ) cat(glue::glue("{mssg} Built general stats...{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Reproducibility Stats::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Delta Beta Stats::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  dB_tib <- tib %>% 
    dplyr::select( dplyr::all_of( c(prbs_key, datA_key, datB_key) ) ) %>%
    dplyr::mutate( dB = base::abs( !!datA_sym - !!datB_sym) )
  print_tib( x = dB_tib, name = "dB Tibble", vb=vb,vt=vt+6,tc=tc+1 )
  
  dB_na_cnt <- dB_tib %>% 
    dplyr::filter( is.na(dB) ) %>% base::nrow()
  dB_pass_cnt <- dB_tib %>%
    dplyr::filter( !is.na(dB) & dB <= dB_min ) %>%
    base::nrow()
  dB_pass_per <- round( 100 * dB_pass_cnt / tot_cnt )
  
  if ( vb >= vt+3 ) cat(glue::glue("{mssg} Built Delta Beta stats...{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           R-Squared Correlation::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  dat_mat <- tib %>%
    dplyr::select( dplyr::all_of( c(prbs_key, datA_key, datB_key) ) ) %>%
    tibble::column_to_rownames( var = prbs_key ) %>%
    as.data.frame() %>% as.matrix()
  if ( vb >= vt+6 ) dat_mat %>% head(n=3) %>% print()
  
  r2_val <- NULL
  r2_val <- dat_mat %>% 
    stats::cor( use = "pairwise.complete.obs", method = "pearson" ) %>% 
    as.data.frame() %>% head(n=1) %>% dplyr::pull(2)
  
  if ( vb >= vt+3 ) cat(glue::glue("{mssg} Built R-squared stats...{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Build Summary Sample Sheet::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_tib <- tibble::tibble(
    Total_cnt = tot_cnt,
    Valid_cnt = pas_cnt,
    Valid_per = pas_per,
    
    dB_NA_cnt = dB_na_cnt,
    dB_pass_cnt = dB_pass_cnt,
    dB_pass_per = dB_pass_per,
    
    r2 = r2_val
  )
  
  # warn_mssg <- glue::glue("WARN_MESSAGE")
  # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
  # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
  # wflag <- FALSE
  # 
  # errs_mssg <- glue::glue("ERROR_MESSAGE")
  # if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return(NULL)
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

analyze_cutlist = function( datA, datB, sel_vec,
                            
                            prbs_key = "Probe_ID",
                            beta_key = "Beta_Med", # "Beta_Avg"
                            vars_key = "Beta_Mad",
                            pval_key = "Pval_Med",
                            
                            samp_key,
                            samA_key = "EWAS",
                            samB_key = "EPIC",
                            
                            dB_min = 0.2,
                            var_min = 0.2,
                            
                            key_vec = NULL,
                            # val_vec = NULL, # c(),
                            # sam_vec = NULL,
                            # sam_min = NULL,
                            
                            out_dir,
                            run_tag = "",
                            
                            reload     = 0,
                            reload_min = 2,
                            reload_pre = NULL,
                            
                            ret_data   = FALSE,
                            parallel   = FALSE,
                            
                            vb=0, vt=3, tc=1, tt=NULL,
                            fun_tag='analyze_cutlist')
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
  
  if ( length(key_vec) != 0 &&
       length(key_vec) == length(val_vec) ) out_tag <- paste0( 
         out_tag, 
         ".F_", paste( key_vec, collapse = "." ),
         ".F_", paste( val_vec, collapse = "." ) )
  
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1, tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1, tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}     prbs_key = '{prbs_key}'.{RET}"))
    cat(glue::glue("{mssg}     beta_key = '{beta_key}'.{RET}"))
    cat(glue::glue("{mssg}     vars_key = '{vars_key}'.{RET}"))
    cat(glue::glue("{mssg}     pval_key = '{pval_key}'.{RET}"))
    cat(glue::glue("{mssg}     samp_key = '{samp_key}'.{RET}"))
    cat(glue::glue("{mssg}     samA_key = '{samA_key}'.{RET}"))
    cat(glue::glue("{mssg}     samB_key = '{samB_key}'.{RET}"))
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
  
  errs_mssg <- glue::glue("Data Set A is null or empty")
  if ( datA %>% base::nrow() == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("Data Set B is null or empty")
  if ( datB %>% base::nrow() == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  errs_mssg <- glue::glue("Selection Vector is null or empty")
  if ( sel_vec %>% length() == 0 ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  prbs_sym <- rlang::sym( prbs_key )
  beta_sym <- rlang::sym( beta_key )
  
  datA_key <- paste( beta_key,samA_key, sep="_" )
  datB_key <- paste( beta_key,samB_key, sep="_" )
  
  varA_key <- paste( vars_key,samA_key, sep="_" )
  varB_key <- paste( vars_key,samB_key, sep="_" )
  
  varA_sym <- rlang::sym( varA_key )
  varB_sym <- rlang::sym( varB_key )
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    print_tib( x = datA, name = "DataSetA", vb=vb,vt=vt+4,tc=tc+1 )
    print_tib( x = datB, name = "DataSetB", vb=vb,vt=vt+4,tc=tc+1 )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Join DataSets::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    join_tib1 <- dplyr::left_join( datA, datB, by=c(prbs_key),
                                   suffix=c(paste0("_",samA_key),
                                            paste0("_",samB_key) )
    )
    # join_tib1 <- dplyr::left_join(
    #   dplyr::select( datA, !!prbs_sym, !!beta_sym ),
    #   dplyr::select( datB, !!prbs_sym, !!beta_sym ),
    #   by=c(prbs_key),
    #   suffix=c(paste0("_",samA_key),paste0("_",samB_key) )
    # )
    
    ssh_tib1 <- compare_samples( tib = join_tib1, 
                                 prbs_key = prbs_key, 
                                 datA_key = datA_key, 
                                 datB_key = datB_key, 
                                 dB_min = dB_min,
                                 vb=vb,vt=vt+1,tc=tc+1 )
    # ssh_tib1 %>% print()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Filter DataSets:: User Vector
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    join_tib2 <- join_tib1 %>% 
      dplyr::filter( !!prbs_sym %in% sel_vec )
    
    ssh_tib2 <- compare_samples( tib = join_tib2, 
                                 prbs_key = prbs_key,
                                 datA_key = datA_key,
                                 datB_key = datB_key,
                                 dB_min = dB_min,
                                 vb=vb,vt=vt+1,tc=tc+1 )
    # ssh_tib2 %>% print()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Filter DataSets:: Variance
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    join_tib3 <- join_tib2 %>% 
      dplyr::filter( !!varA_sym <= var_min ) %>%
      dplyr::filter( !!varB_sym <= var_min )
    
    ssh_tib3 <- compare_samples( tib = join_tib2, 
                                 prbs_key = prbs_key,
                                 datA_key = datA_key,
                                 datB_key = datB_key,
                                 dB_min = dB_min,
                                 vb=vb,vt=vt+1,tc=tc+1 )
    # ssh_tib2 %>% print()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Build Summary Sample Sheet::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- tibble::tibble(
      Sample_Name = samp_key,
      Product_A = samA_key,
      Product_B = samB_key,
      Beta_Metric = !!beta_key,
      Variance_Metric = !!vars_key,
      Variance_Cutoff = var_min,
      
      Product_A_cnt = datA %>% base::nrow(),
      Product_B_cnt = datB %>% base::nrow(),
      
      # Product_Intersect_cnt = join_cnt1,
      
      dB_pass_cnt_delta1 = ssh_tib2$dB_pass_cnt - ssh_tib1$dB_pass_cnt,
      dB_pass_per_delta1 = ssh_tib2$dB_pass_per - ssh_tib1$dB_pass_per,
      r2_delta1 = ssh_tib2$r2 - ssh_tib1$r2,
      r2_delta_per1 = round( 100*r2_delta1 / ( ssh_tib1$r2 ), 3 ),
      
      dB_pass_cnt_delta2 = ssh_tib3$dB_pass_cnt - ssh_tib1$dB_pass_cnt,
      dB_pass_per_delta2 = ssh_tib3$dB_pass_per - ssh_tib1$dB_pass_per,
      r2_delta2 = ssh_tib3$r2 - ssh_tib1$r2,
      r2_delta_per2 = round( 100*r2_delta2 / ( ssh_tib1$r2 ), 3 ),
      
    ) %>%
      dplyr::bind_cols( ssh_tib1 %>% purrr::set_names(paste(names(.),"1", sep="_")) ) %>%
      dplyr::bind_cols( ssh_tib2 %>% purrr::set_names(paste(names(.),"2", sep="_")) ) %>%
      dplyr::bind_cols( ssh_tib3 %>% purrr::set_names(paste(names(.),"3", sep="_")) )
    
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
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                 Analyze Matrix Cluster Separation Score::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Pretty Sure this Function is not used and should be removed...
#

# Inputs: csv, out_dir, 
# Filter: key_vec, val_vec, sam_vec, sam_min,
# dB_col: dat_col, 
# 
analyze_css_matrix = function( csv,
                               
                               key_vec = NULL, # c(),
                               val_vec = NULL, # c(),
                               
                               sam_vec = NULL,
                               sam_min = NULL,
                               
                               out_dir,
                               run_tag = "",
                               
                               reload     = 0,
                               reload_min = 2,
                               reload_pre = NULL,
                               
                               ret_data   = FALSE,
                               parallel   = FALSE,
                               
                               vb=0, vt=3, tc=1, tt=NULL,
                               fun_tag='analyze_css_matrix')
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
  
  if ( length(key_vec) != 0 &&
       length(key_vec) == length(val_vec) ) out_tag <- paste0( 
         out_tag, 
         ".F_", paste( key_vec, collapse = "." ),
         ".F_", paste( val_vec, collapse = "." ) )
  
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1 )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_csv, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1, tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}          csv = '{csv}'.{RET}"))
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
  
  errs_mssg <- glue::glue("File file='{csv}' does not exist")
  if ( !file.exists( csv) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    ann_tib <- safe_read( csv, vb=vb,vt=vt+1,tc=tc, tt=tt )
    
    ret_tib <- ann_tib
    
    for ( f_idx in c(1:length(key_vec) ) ) {
      cur_key_sym <- rlang::sym( key_vec[f_idx] )
      cur_val_sym <- rlang::sym( val_vec[f_idx] )
      
      print_tib( ret_tib, fun_tag = fun_tag, 
                 name = paste( "pre",key_vec[f_idx],val_vec[f_idx], sep="." ),
                 vb=vb,vt=vt+3,tc=tc+1, tt=tt )
      
      ret_tib <- ret_tib %>% dplyr::filter( !!cur_key_sym == val_vec[f_idx] )
      
      print_tib( ret_tib, fun_tag = fun_tag, 
                 name = paste( "aft",key_vec[f_idx],val_vec[f_idx], sep="." ),
                 vb=vb,vt=vt+3,tc=tc+1, tt=tt )
    }
    # ret_tib <- ann_tib %>% head()
    
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
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Binary Gizipped IO::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

list_to_table = function( list, 
                          col = "Probe_ID", 
                          row = "Sentrix_Name", 
                          val = "pval",
                          precision = 0,
                          matrix = FALSE,
                          vb=0, vt=3, tc=1, tt=NULL ) {
  
  col_sym = rlang::sym( col )
  row_sym = rlang::sym( row )
  val_sym = rlang::sym( val )
  
  ret_tib <- NULL
  ret_tib <- list %>% 
    dplyr::bind_rows( .id = row ) %>%
    tidyr::pivot_wider( names_from = !!row_sym,
                        values_from = c( !!val_sym ) ) %>%
    tibble::column_to_rownames( var = col ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column( var = col ) %>%
    tibble::as_tibble()
  
  if ( precision > 0 ) ret_tib <- ret_tib %>% dplyr::mutate(
    dplyr::across( purrr::is_double, function(x, y) {
      x = as.integer( x * y) }, precision ) )
  
  if ( matrix ) ret_tib <- ret_tib %>% 
    tibble::column_to_rownames( var = col ) %>% as.matrix()
  
  if ( vb >= vt ) ret_tib %>% head( n=3 ) %>% print()
  
  ret_tib
}

sdf_to_pval_tib = function( sdf, detp, ids=NULL, pval = 1.0,
                            
                            vb=0, vt=3, tc=1, tt=NULL,
                            fun_tag='sdf_to_pval_tib')
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
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  
  if ( detp == "none_none" ) {
    det_vec = sdf %>% 
      sesame::getBetas( mask = FALSE,
                        sum.TypeI = FALSE )
    
  } else if ( detp == "none_dyen" ) {
    det_vec = sdf %>% 
      sesame::dyeBiasCorrTypeINorm() %>%
      sesame::getBetas( mask = FALSE,
                        sum.TypeI = FALSE )
  } else if ( detp == "none_noob" ) {
    det_vec = sdf %>% 
      sesame::noob() %>% 
      sesame::getBetas( mask = FALSE,
                        sum.TypeI = FALSE )
  } else if ( detp == "noob_dyen" ) {
    det_vec = sdf %>% 
      sesame::noob() %>% 
      sesame::dyeBiasCorrTypeINorm() %>%
      sesame::getBetas( mask = FALSE,
                        sum.TypeI = FALSE )
  } else if ( detp == "poob" ) {
    det_vec = sesame::pOOBAH( sdf = sdf, return.pval = TRUE,  combine.neg = FALSE, pval.threshold = pval )
  } else if ( detp == "boob" ) {
    det_vec = sesame::pOOBAH( sdf = sdf, return.pval = TRUE,  combine.neg = TRUE, pval.threshold = pval )
  } else if ( detp == "negs") {
    if ( is.null( attr( sdf, "controls" ) ) ) {
      errs_mssg <- glue::glue("Sample does not have controls!")
      if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    det_vec = sesame::detectionPnegEcdf( sdf = sdf, return.pval = TRUE,  pval.threshold = pval )
  } else {
    eflag <- TRUE
    errs_mssg <- glue::glue("Unsupported detection p-value method = '{detp}'")
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
  }
  ret_tib <- tibble::tibble(
    Probe_ID = names(det_vec),
    pval = det_vec
  )
  
  # Filter Probes::
  if ( !is.null(ids) ) ret_tib <- ret_tib %>% 
    dplyr::filter( Probe_ID %in% ids )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

write_precision_matrix = function( sdfs,
                                   out_dir,
                                   manifest,
                                   
                                   samp_tag,
                                   detp_types,
                                   workflowsA,
                                   workflowsB,
                                   single = FALSE,
                                   precision = 1000000,
                                   data_max = 0,
                                   
                                   vb=0, vt=3, tc=1, tt=NULL,
                                   fun_tag='write_precision_matrix')
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
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}     samp_tag = '{samp_tag}'.{RET}"))
    cat(glue::glue("{mssg}   detp_types = '{detp_types}'.{RET}"))
    cat(glue::glue("{mssg}   workflowsA = '{workflowsA}'.{RET}"))
    cat(glue::glue("{mssg}   workflowsB = '{workflowsB}'.{RET}"))
    cat(glue::glue("{mssg}       single = '{single}'.{RET}"))
    cat(glue::glue("{mssg}    precision = '{precision}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  out_dir <- file.path( out_dir, fun_tag )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    safe_mkdir( out_dir, vb=vb )
    if ( vb >= vt ) cat(glue::glue("{mssg} Built out_dir='{out_dir}'{RET}"))
    
    # blank <- foreach ( wrk2_type = wrk2_types, .combine = rbind ) %dopar% {
    # detp_dats <- foreach ( detp_type = detp_types, .combine = rbind ) %dopar% {
    
    beta_col_vec <- colnames( beta_tib %>% dplyr::select( -Probe_ID ) )
    beta_row_vec <- beta_tib$Probe_ID
    
    beta_row_csv <- file.path( cur_dir, paste0( out_str,".beta.row.csv.gz" ))
    beta_col_csv <- file.path( cur_dir, paste0( out_str,".beta.col.csv.gz" ))
    readr::write_lines( x = beta_col_vec, file = beta_col_csv )
    readr::write_lines( x = beta_row_vec, file = beta_row_csv )
    
    for ( workflowA in workflowsA ) {  # break }
      work_strA <- struct_to_str( workflowA )
      
      for ( detp_type in detp_types ) { # break }
        #
        # Generate P-values for Current Workflow::
        #
        detp_str <- struct_to_str( detp_type )
        detp_dir <- file.path( out_dir, work_strA, detp_str ) %>% safe_mkdir()
        
        pvals <- sdfs %>% parallel::mclapply( sdf_to_pval_tib, 
                                              detp = detp_type, 
                                              ids  = manifest$Probe_ID, 
                                              pval = 1.0, 
                                              vb=vb,vt=vt+1,tc=tc+1 )
        
        pval_tab <- list_to_table( pvals,
                                   col = "Probe_ID", 
                                   row = "Sentrix_Name", 
                                   val = "pval",
                                   precision = precision,
                                   vb=vb,vt=vt+1,tc=tc+1 )
        
        if ( p1 ) cat(glue::glue("{mssg} Generated Pvals.{RET}"))
        
        #
        # Write bgz table and auxilary files::
        #
        
        
        # TBD: pval_param_str
        # TBD: pval_dir
        pval_mat_csv <- file.path( cur_dir, paste0( out_str,".pval.mat.csv" ))
        readr::write_csv( x = pval_tib[,-1], file = pval_mat_csv, col_names = FALSE )
      }
      
      
      for ( workflowB in workflowsB ) {  # break }
        #
        # Output Strings::
        #
        work_strB <- struct_to_str( workflowB )
        
        out_str <- paste( detp_str, work_strA, work_strB, sep="." )
        if ( vb >= vt ) cat(glue::glue("{mssg} out_str='{out_str}'{RET}"))
        cur_dir <- file.path( out_dir, detp_str, work_strA, work_strB )
        safe_mkdir( cur_dir, vb=vb )
        if ( vb >= vt ) cat(glue::glue("{mssg} out_str='{out_str}'{RET}"))
        if ( vb >= vt ) cat(glue::glue("{mssg} cur_dir='{cur_dir}'{RET}"))
        
        #
        # Generate Calls for Current Workflow::
        #
        calls <- NULL
        calls <- sdfs %>% head() %>%
          parallel::mclapply( mutate_sdf_to_call,
                              ids   = manifest$Probe_ID,
                              detp  = detp_type,
                              work1 = workflowA,
                              work2 = workflowB,
                              pval = 1.00,
                              vb=vb, vt=vt+1 )
        
        calls_tib <- calls %>% dplyr::bind_rows( .id = "Sentrix_Name" )
        if ( vb >= vt ) cat(glue::glue("{mssg} Genereateed Calls.{RET}"))
        
        #
        # Generate Unfiltered Beta Calls::
        #
        beta_tib <- calls_tib %>% 
          dplyr::select( Probe_ID, Sentrix_Name, beta ) %>%
          tidyr::pivot_wider( names_from = Sentrix_Name,
                              values_from = c( beta ) ) %>% 
          tibble::column_to_rownames( var = "Probe_ID" ) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column( var = "Probe_ID" ) %>% 
          tibble::as_tibble() %>% 
          dplyr::mutate(
            dplyr::across( purrr::is_double, function(x, y) { 
              x = as.integer( x * y) }, precision ) )
        beta_tib %>% head( n=3 ) %>% print()
        if ( vb >= vt ) cat(glue::glue("{mssg} Genereateed Betas.{RET}"))
        
        #
        # Generate Params for Current Workflow::
        #
        param_tib <- tibble::tibble(
          sample_tag = samp_tag,
          nrows_cnt  = sample_sheet_tib %>% base::nrow(),
          ncols_cnt  = sample_sheet_tib %>% base::ncol(),
          detp_type  = detp_type,
          workflow1  = work_strA,
          workflow2  = work_strB,
          precision  = precision
        )
        param_tib %>% print()
        if ( vb >= vt ) cat(glue::glue("{mssg} Genereateed Params.{RET}"))
        
        #
        # Write All Current Outputs::
        #
        beta_mat_csv <- file.path( cur_dir, paste0( out_str,".beta.mat.csv" ))
        params_csv <- file.path( cur_dir, paste0( out_str,".params.csv.gz" ))
        
        readr::write_csv( x = beta_tib[,-1], file = beta_mat_csv, col_names = FALSE )
        readr::write_csv( x = param_tib, file = params_csv )
        
        cur_tib <- tibble::tibble(
          out_dir    = cur_dir,
          prefix     = out_str,
          beta_csv   = beta_mat_csv,
          pval_csv   = pval_mat_csv,
          rows_csv   = beta_row_csv,
          cols_csv   = beta_col_csv,
          pars_csv   = params_csv
        )
        
        ret_tib <- ret_tib %>% dplyr::bind_rows( cur_tib )
        
        if ( vb >= vt ) cat(glue::glue("{mssg} Done.{RET2}"))
      }
      if ( single ) break
    }
    if ( single ) break
    
    # success <- write_df_bgz_rcpp( df = ret_tib, 
    #                              spec_str = spec_str1, 
    #                              file = file, 
    #                              precision = precision, 
    #                              vb=vb,vt=vt+1,tc=tc+1 )
    # 
    # if ( !success ) {
    #   errs_mssg <- glue::glue("Failed duing write_df_bgz_rcpp()!")
    #   if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    #   if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    #   if ( eflag ) return(NULL)
    # }
    
    # warn_mssg <- glue::glue("WARN_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # wflag <- FALSE
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

write_df_bgz = function( tib,
                         file,
                         precision = 1000000,
                         
                         vb=0, vt=3, tc=1, tt=NULL,
                         fun_tag='write_df_bgz')
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
    cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}    precision = '{precision}'.{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ftime <- base::system.time({
    
    safe_mkdir( file )
    
    spec_str1 <- tib %>% 
      summarise_all(class) %>% 
      gather(col_name, col_type) %>% 
      dplyr::pull(col_type) %>% stringr::str_sub(1,1) %>% 
      paste0(collapse = '')
    
    ret_tib <- tib %>% dplyr::mutate( 
      dplyr::across( purrr::is_double, function(x, y) { 
        x = as.integer( x * y) }, precision ) )
    
    spec_str2 <- ret_tib %>% 
      summarise_all(class) %>% 
      gather(col_name, col_type) %>% 
      dplyr::pull(col_type) %>% stringr::str_sub(1,1) %>% 
      paste0(collapse = '')
    
    tib_spec <- 
      readr::read_csv(
        readr::format_csv( ret_tib, col_names = TRUE ), 
        col_types = spec_str1 ) %>%
      readr::spec()
    
    # names(tib_spec$cols)
    
    success <- write_df_bgz_rcpp( df = ret_tib, 
                                  spec_str = spec_str1, 
                                  file = file, 
                                  precision = precision, 
                                  vb=vb,vt=vt+1,tc=tc+1 )
    
    if ( !success ) {
      errs_mssg <- glue::glue("Failed duing write_df_bgz_rcpp()!")
      if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
      if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
      if ( eflag ) return(NULL)
    }
    
    # warn_mssg <- glue::glue("WARN_MESSAGE")
    # if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    # if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # wflag <- FALSE
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat$spec_str1 <- spec_str1
  ret_dat$spec_str2 <- spec_str2
  
  ret_dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Summarize Calls::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

qcdf_to_tib = function( qdf,
                        rnd = 3,
                        vb=0, vt=6, tc=1, tt=NULL,
                        fun_tag='qcdf_to_tib')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  ret_tib <- tibble::tibble( a = 0 )
  for ( qc_col in names(qdf) ) {
    # cat(glue::glue("Current name = '{qc_col}'{RET}"))
    qc_sym <- rlang::sym(qc_col)
    ret_tib <- ret_tib %>% tibble::add_column( !!qc_sym := qdf[[qc_col]] )
  }
  ret_tib <- ret_tib %>% dplyr::select( -a )
  ret_cnt <- ret_tib %>% base::nrow()
  
  if ( rnd != 0 ) ret_tib <- ret_tib %>% 
    dplyr::mutate( dplyr::across( purrr::is_numeric, round, rnd) )
  
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

sdf_to_summary = function( sdf,
                           types = "cg",
                           pvals = 0.05,
                           detps = "poob",
                           mask  = FALSE,
                           best  = FALSE,
                           rnd = 3,
                           
                           vb=0, vt=3, tc=1, tt=NULL,
                           fun_tag='sdf_to_summary')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  # pval_sym <- rlang::sym( paste("pval_",detp,type, sep="_" ) )
  # beta_sym <- rlang::sym( paste("beta_",detp,type, sep="_" ) )
  
  ret_tab <- NULL
  for ( detp in detps ) {
    for ( pval in pvals ) {
      pval_int <- pval * 100;
      
      cur_sdf <- NULL
      if ( detp == "poob" ) {
        det_vec = sesame::pOOBAH( sdf = sdf, return.pval = TRUE,  combine.neg = FALSE, pval.threshold = pval )
        cur_sdf = sesame::pOOBAH( sdf = sdf, return.pval = FALSE, combine.neg = FALSE, pval.threshold = pval )
      } else if ( detp == "boob" ) {
        det_vec = sesame::pOOBAH( sdf = sdf, return.pval = TRUE,  combine.neg = TRUE, pval.threshold = pval )
        cur_sdf = sesame::pOOBAH( sdf = sdf, return.pval = FALSE, combine.neg = TRUE, pval.threshold = pval )
      } else if ( detp == "negs") {
        
        if ( is.null( attr( sdf, "controls" ) ) ) next
        
        det_vec = sesame::detectionPnegEcdf( sdf = sdf, return.pval = TRUE,  pval.threshold = pval )
        cur_sdf = sesame::detectionPnegEcdf( sdf = sdf, return.pval = FALSE, pval.threshold = pval )
      } else {
        eflag <- TRUE
        errs_mssg <- glue::glue("Unsupported detection p-value method = '{detp}'")
        if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
        if ( eflag ) return( NULL )
      }
      
      for ( type in types ) {
        prb_tib <- tibble::as_tibble( det_vec ) %>% 
          dplyr::mutate( Probe_ID = names(det_vec), 
                         Probe_Type = stringr::str_sub(Probe_ID, 1,2) ) %>%
          dplyr::filter( Probe_Type == type )
        # print(prb_tib)
        
        num_probes <- prb_tib %>% base::nrow()
        # print(num_probes)
        if ( num_probes == 0 ) next
        num_nondt <- prb_tib %>% dplyr::filter( value > pval ) %>% base::nrow()
        # print(num_nondt)
        
        cur_tab <- NULL
        cur_tab <- tibble::tibble(
          probe_type = type,
          detp_type  = detp,
          pval_min   = pval_int,
          num_probes = num_probes,
          num_nondt  = num_nondt,
          per_nondt  = base::round( 100* num_nondt/num_probes, rnd )
        )
        # %>% purrr::set_names( paste(names(.),detp,type,pval_int, sep="_" ) ) %>%
        #   as.data.frame()
        
        ret_tab <- ret_tab %>% dplyr::bind_rows( cur_tab )
      }
    }
  }
  
  ret_tib <- NULL
  ret_tib <- ret_tab %>% tidyr::pivot_wider( 
    id_cols = c(detp_type, probe_type, pval_min), 
    names_from = c(detp_type, probe_type, pval_min), 
    names_sep = "_", 
    values_from = c(num_probes, num_nondt, per_nondt) )
  
  # qc_sym <- rlang::sym(qc_col)
  # ret_tib <- ret_tib %>% tibble::add_column( !!qc_sym := qdf[[qc_col]] )
  
  if ( rnd != 0 ) ret_tib <- ret_tib %>% 
    dplyr::mutate( dplyr::across( purrr::is_numeric, round, rnd) )
  
  #
  # SNP's
  #
  # annoS = sesameDataGetAnno("EPIC/EPIC.hg19.snp_overlap_b151.rds")
  # annoI = sesameDataGetAnno("EPIC/EPIC.hg19.typeI_overlap_b151.rds")
  # formatVCF(sdf, annoS=annoS, annoI=annoI)
  
  if (FALSE) {
    ret_tib <- NULL
    ret_tib <- tibble::tibble(
      Probe_ID = sdf$Probe_ID,
      Loci_ID  = sdf$Probe_ID %>% stringr::str_remove("_.*$"),
      !!pval_sym := det_vec,
      !!beta_sym := sesame::getBetas( sdf = sdf_dat, 
                                      mask = mask,
                                      sum.TypeI = FALSE )
    )
    
    if ( best ) {
      ret_tib <- ret_tib %>% 
        dplyr::arrange( !!pval_sym ) %>%
        dplyr::distinct( Loci_ID, .keep_all = TRUE )
    }
    # dplyr::arrange(pval_poob) %>% 
    #   dplyr::distinct(Loci_ID, .keep_all = TRUE ) %>% dplyr::arrange(Loci_ID)
  }
  ret_cnt <- ret_tib %>% base::nrow()
  
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                              Convert SDF::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mutate_sdf = function( sdf,
                       detp = NULL,
                       work = NULL,
                       pval = 1.00,
                       
                       vb=0, vt=3, tc=1, tt=NULL,
                       fun_tag='mutate_sdf')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  if ( !is.null(detp) ) {
    for ( d in detp ) {
      if ( d == "" ) {
        # Do nothing...
      } else if ( d == "poob" ) {
        sdf = sesame::pOOBAH( sdf = sdf, return.pval = FALSE, combine.neg = FALSE, pval.threshold = pval )
      } else if ( d == "boob" ) {
        sdf = sesame::pOOBAH( sdf = sdf, return.pval = FALSE, combine.neg = TRUE, pval.threshold = pval )
      } else if ( d == "negs") {
        sdf = sesame::detectionPnegEcdf( sdf = sdf, return.pval = FALSE, pval.threshold = pval )
      } else {
        eflag <- TRUE
        errs_mssg <- glue::glue("Unsupported sesame p-value detection method = '{d}'")
        if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
        if ( eflag ) return( NULL )
      }
    }
  }
  
  if ( !is.null(work) ) {
    for ( w in work ) {
      if ( vb >= vt+1 ) cat(glue::glue("{mssg}{TAB} Mutating with '{w}'...{RET}"))
      
      if ( w == "" ) {
        # Do nothing...
      } else if ( w == "dyel") {
        sdf = sesame::dyeBiasCorr( sdf = sdf )
      } else if ( w == "dyen") {
        sdf = sesame::dyeBiasCorrTypeINorm( sdf = sdf )
        
      } else if ( w == "swap") {
        sdf = sesame::inferInfiniumIChannel( sdf = sdf, switch_failed = FALSE, summary = FALSE )
        
      } else if ( w == "noob") {
        sdf = sesame::noob( sdf = sdf )
      } else if ( w == "neob") {
        sdf = sesame::noob( sdf = sdf )
      } else if ( w == "scrub") {
        sdf = sesame::noob() %>% sesame::scrub( sdf = sdf )
        
      } else {
        eflag <- TRUE
        errs_mssg <- glue::glue("Unsupported sesame method = '{w}'")
        if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
        if ( eflag ) return( NULL )
      }
    }
  }
  
  ret_cnt <- sdf %>% base::nrow()
  
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  sdf
}

sdf_to_pval = function( sdf,
                        detp,
                        pval = 1.00,
                        loci = FALSE,
                        best = FALSE,
                        suffix = NULL,
                        
                        vb=0, vt=3, tc=1, tt=NULL,
                        fun_tag='sdf_to_pval')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  pval_sym <- rlang::sym( "pval" )
  beta_sym <- rlang::sym( "beta" )
  
  if ( !is.null(suffix) ) {
    pval_sym <- rlang::sym( paste0("pval_",suffix ) )
    beta_sym <- rlang::sym( paste0("beta_",suffix ) )
  }
  
  det_vec <- NULL
  if ( detp == "poob" ) {
    det_vec = sesame::pOOBAH( sdf = sdf, return.pval = TRUE,  combine.neg = FALSE, pval.threshold = pval )
  } else if ( detp == "boob" ) {
    det_vec = sesame::pOOBAH( sdf = sdf, return.pval = TRUE,  combine.neg = TRUE, pval.threshold = pval )
  } else if ( detp == "negs") {
    det_vec = sesame::detectionPnegEcdf( sdf = sdf, return.pval = TRUE,  pval.threshold = pval )
  } else {
    eflag <- TRUE
    errs_mssg <- glue::glue("Unsupported detection p-value method = '{detp}'")
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
  }
  
  ret_tib <- NULL
  if ( loci ) {
    ret_tib <- tibble::tibble(
      Probe_ID = sdf$Probe_ID,
      Loci_ID  = sdf$Probe_ID %>% stringr::str_remove("_.*$"),
      !!pval_sym := det_vec )
  } else {
    ret_tib <- tibble::tibble(
      Probe_ID = sdf$Probe_ID,
      !!pval_sym := det_vec )
  }
  
  if ( best ) {
    ret_tib <- ret_tib %>% 
      dplyr::arrange( !!pval_sym ) %>%
      dplyr::distinct( Loci_ID, .keep_all = TRUE )
  }
  ret_cnt <- ret_tib %>% base::nrow()
  
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

sdf_to_beta = function( sdf,
                        pval = 1.00,
                        mask = FALSE,
                        loci = FALSE,
                        suffix = NULL,
                        
                        vb=0, vt=3, tc=1, tt=NULL,
                        fun_tag='sdf_to_beta')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  pval_sym <- rlang::sym( "pval" )
  beta_sym <- rlang::sym( "beta" )
  
  if ( !is.null(suffix) ) {
    pval_sym <- rlang::sym( paste0("pval_",suffix ) )
    beta_sym <- rlang::sym( paste0("beta_",suffix ) )
  }
  
  ret_tib <- NULL
  if ( loci ) {
    ret_tib <- tibble::tibble(
      Probe_ID = sdf$Probe_ID,
      Loci_ID  = sdf$Probe_ID %>% stringr::str_remove("_.*$"),
      !!beta_sym := sesame::getBetas( sdf = sdf, 
                                      mask = mask,
                                      sum.TypeI = FALSE )
    )
  } else {
    ret_tib <- tibble::tibble(
      Probe_ID = sdf$Probe_ID,
      !!beta_sym := sesame::getBetas( sdf = sdf, 
                                      mask = mask,
                                      sum.TypeI = FALSE )
    )
  }
  ret_cnt <- ret_tib %>% base::nrow()
  
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

sdf_to_call = function( sdf,
                        detp,
                        pval = 1.00,
                        mask = FALSE,
                        loci = FALSE,
                        best = FALSE,
                        suffix = FALSE,
                        
                        vb=0, vt=3, tc=1, tt=NULL,
                        fun_tag='sdf_to_call')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  pval_sym <- rlang::sym( "pval" )
  beta_sym <- rlang::sym( "beta" )
  
  if ( suffix ) {
    pval_sym <- rlang::sym( paste0("pval_",detp ) )
    beta_sym <- rlang::sym( paste0("beta_",detp ) )
  }
  
  det_vec <- NULL
  if ( detp == "poob" ) {
    det_vec = sesame::pOOBAH( sdf = sdf, return.pval = TRUE,  combine.neg = FALSE, pval.threshold = pval )
  } else if ( detp == "boob" ) {
    det_vec = sesame::pOOBAH( sdf = sdf, return.pval = TRUE,  combine.neg = TRUE, pval.threshold = pval )
  } else if ( detp == "negs") {
    det_vec = sesame::detectionPnegEcdf( sdf = sdf, return.pval = TRUE,  pval.threshold = pval )
  } else {
    eflag <- TRUE
    errs_mssg <- glue::glue("Unsupported detection p-value method = '{detp}'")
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
  }
  
  if ( loci ) {
    ret_tib <- NULL
    ret_tib <- tibble::tibble(
      Probe_ID = sdf$Probe_ID,
      Loci_ID  = sdf$Probe_ID %>% stringr::str_remove("_.*$"),
      !!pval_sym := det_vec,
      !!beta_sym := sesame::getBetas( sdf = sdf,
                                      mask = mask,
                                      sum.TypeI = FALSE )
    )
  } else {
    ret_tib <- NULL
    ret_tib <- tibble::tibble(
      Probe_ID = sdf$Probe_ID,
      !!pval_sym := det_vec,
      !!beta_sym := sesame::getBetas( sdf = sdf,
                                      mask = mask,
                                      sum.TypeI = FALSE )
    )
  }
  
  if ( best ) {
    ret_tib <- ret_tib %>%
      dplyr::arrange( !!pval_sym ) %>%
      dplyr::distinct( Loci_ID, .keep_all = TRUE )
  }
  ret_cnt <- ret_tib %>% base::nrow()
  
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

mutate_sdf_to_call = function( sdf,
                               
                               ids = NULL,
                               detp = NULL,
                               work1 = NULL,
                               work2 = NULL,
                               
                               pval = 1.00,
                               mask = FALSE,
                               loci = FALSE,
                               best = FALSE,
                               suffix = FALSE,
                               
                               vb=0, vt=3, tc=1, tt=NULL,
                               fun_tag='mutate_sdf_to_call')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  pval_sym <- rlang::sym( "pval" )
  beta_sym <- rlang::sym( "beta" )
  
  if ( suffix ) {
    pval_sym <- rlang::sym( paste0("pval_",detp ) )
    beta_sym <- rlang::sym( paste0("beta_",detp ) )
  }
  
  sdf <- sdf %>% mutate_sdf( detp = detp,
                             work = work1, 
                             pval = 1.00,
                             vb=vb,vt=vt+1,tc=tc+1 )
  
  pval_tib <- NULL
  pval_tib <- sdf %>% sdf_to_pval( detp = detp,
                                   pval = 1.00,
                                   loci = loci,
                                   best = best,
                                   vb=vb,vt=vt+1,tc=tc+1 )
  
  sdf <- sdf %>% mutate_sdf( detp = detp,
                             work = work2,
                             pval = pval, 
                             vb=vb,vt=vt+1,tc=tc+1 )
  
  beta_tib <- NULL
  beta_tib <- sdf %>% sdf_to_beta( pval = pval,
                                   mask = mask,
                                   loci = loci,
                                   suffix = NULL,
                                   vb=vb,vt=vt+1,tc=tc+1 )
  
  ret_tib <- dplyr::inner_join( pval_tib, beta_tib, by="Probe_ID" )
  # ret_cnt <- ret_tib %>% base::nrow()
  
  if ( !is.null(ids) ) ret_tib <- ret_tib %>% dplyr::filter( Probe_ID %in% ids )
  
  ret_key <- glue::glue("final-ret-tib")
  ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

clean_sdf = function( sdf,
                      vb=0, vt=6, tc=1, tt=NULL,
                      fun_tag='clean_sdf')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  ret_tib <- NULL
  ret_tib <- sdf %>% dplyr::filter( !is.na(UG) | !is.na(UR) ) %>%
    #
    # NOTE: Not totally sure if these check are need, but why not...
    #
    dplyr::filter( !(is.na(MG) & is.na(MR) & is.na(UG) & is.na(UR) ) ) %>%
    dplyr::filter( !(is.na(MG) & is.na(MR) & col != 2) ) %>%
    dplyr::filter( !(is.na(UG) & is.na(UR) & col == 2) )
  
  ret_cnt <- ret_tib %>% base::nrow()
  
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                             Load Raw SDF::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prefix_to_sdf2 = function( prefix,
                           
                           platform = NULL,
                           manifest = NULL,
                           controls = NULL,
                           
                           add_mask = FALSE,
                           
                           out_dir,
                           run_idat = TRUE,
                           # run_tag,
                           
                           run_qc     = TRUE,
                           run_ctl    = TRUE,
                           clean_sdf  = TRUE,
                           save_light = TRUE,
                           
                           round      = 3,
                           reload     = 0,
                           reload_min = 2,
                           reload_pre = NULL,
                           
                           ret_data   = FALSE,
                           parallel   = FALSE,
                           
                           vb=0, vt=3, tc=1, tt=NULL,
                           fun_tag='prefix_to_sdf2')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  ftime <- 0
  etime <- 0
  wflag <- FALSE
  eflag <- FALSE
  
  run_tag <- prefix %>% base::basename()
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_rds <- file.path( out_dir, paste(out_tag, 'rds', sep='.') )
  qdf_rds <- file.path( out_dir, paste(out_tag, 'qdf.rds', sep='.') )
  ctl_csv <- file.path( out_dir, paste(out_tag, 'ctl.csv.gz', sep='.') )
  pdf_csv <- file.path( out_dir, paste(out_tag, 'pdf.csv.gz', sep='.') )
  beg_txt <- paste(out_rds, 'start.txt', sep='.')
  end_txt <- paste(out_rds, 'done.txt', sep='.')
  
  sigs_csv <- file.path( out_dir, paste(out_tag, 'idat.sigs.csv.gz', sep='.') )
  info_csv <- file.path( out_dir, paste(out_tag, 'idat.info.csv.gz', sep='.') )
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1 )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, sigs_csv, info_csv, 
                                  qdf_rds, pdf_csv, out_rds, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
  if ( reload >= reload_min && is_valid ) {
    if ( vb >= vt+3 ) cat(glue::glue("{mssg} Reloading rds = '{out_rds}'...{RET}"))
    return( readr::read_rds( out_rds ) )
  }
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}       prefix = '{prefix}'.{RET}"))
    cat(glue::glue("{mssg}      run_qc  = '{run_qc}'.{RET}"))
    cat(glue::glue("{mssg}      run_ctl = '{run_ctl}'.{RET}"))
    cat(glue::glue("{mssg}   clean_sdf  = '{clean_sdf}'.{RET}"))
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
    cat(glue::glue("{mssg}      out_rds = '{out_rds}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c( sigs_csv, info_csv, qdf_rds, pdf_csv, ctl_csv, out_rds, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  prefix_dir <- base::dirname(prefix)
  
  errs_mssg <- glue::glue("Directory ='{prefix_dir}' does not exist")
  if ( !dir.exists( prefix_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Get Idat Stats and Subset Manifest::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # TBD:: Screen Control Addresses as well...
    #
    if ( run_idat ) {
      if ( vb >= vt+8 ) manifest %>% head() %>% print()
      
      raw_idat_pair_tib <- prefixToIdat( prefix = prefix, 
                                         load = TRUE, 
                                         save = save_light, 
                                         light = save_light,
                                         csv = sigs_csv,
                                         ssh = info_csv, 
                                         gzip = TRUE, validate = FALSE,
                                         verbose = vb, vt=vt+8,tc=tc+1, tt=tt )
      
      if ( vb >= vt+8 ) raw_idat_pair_tib %>% head() %>% print()
      raw_idat_pair_tib %>% head() %>% print()
      
      beg_man_cnt <- manifest %>% base::nrow()
      
      manifest <- manifest %>%
        dplyr::filter( U %in% raw_idat_pair_tib$Address ) %>% 
        dplyr::filter( M %in% raw_idat_pair_tib$Address | is.na(M) ) %>% 
        dplyr::filter( !(is.na(U) ) ) %>%
        dplyr::filter( !(Infinium_Design_Type == 1 & is.na(M) ) ) %>%
        dplyr::distinct( U,M, .keep_all = TRUE )
      
      end_man_cnt <- manifest %>% base::nrow()
      
      if ( vb >= vt+2 ) cat(glue::glue("{mssg} beg_man_cnt={beg_man_cnt}{RET}",
                                       "{mssg} end_man_cnt={end_man_cnt}{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Sesame Load Idats::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( vb >= vt+8 ) manifest %>% head() %>% print()
    if ( vb >= vt+8 ) controls %>% head() %>% print()
    
    ret_dat <- sesame::readIDATpair( prefix.path = prefix, 
                                     platform = platform, 
                                     manifest = manifest,
                                     controls = controls,
                                     verbose = vb >= vt )
    
    if ( add_mask ) {
      mask_vec <- manifest %>% 
        dplyr::filter( mask ) %>% dplyr::pull( Probe_ID ) %>% as.vector()
      mask_len <- mask_vec %>% length()
      
      if ( vb >= vt+2 ) 
        cat(glue::glue("{mssg} Adding mask len={mask_len}{RET}"))
      
      if ( vb >= vt+8 ) head( mask_vec ) %>% print()
      ret_dat  <- sesame::addMask( sdf = ret_dat, probes = mask_vec )
      if ( vb >= vt+8 ) print( ret_dat )
    }
    if ( vb >= vt+8 ) ret_dat %>% head() %>% print()
    
    if ( clean_sdf ) ret_dat <- ret_dat %>% clean_sdf( vb=vb,vt=vt,tc=tc )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Build QC Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( run_qc ) {
      
      qdf_dat <- sesame::sesameQC( ret_dat )
      pdf_tib <- qcdf_to_tib( qdf_dat, rnd = round )
      
      readr::write_rds( qdf_dat, file = qdf_rds, compress = "gz" )
      safe_write( pdf_tib, pdf_csv, vb=vb )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Write Control Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # TBD:: Write Controls::
    # sesame::controls( cur_sdfs$`206662930005_R01C01` )
    # sesame::controls( cur_sdfs$`206662930005_R01C01` ) %>% as.data.frame() %>% tibble::as_tibble( rownames = "Probe_ID" ) %>% dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace_all("\\.+", "_") )
    
    if ( run_ctl ) {
      ctl_tib <- NULL
      ctl_tib <- sesame::controls( ret_dat ) %>% 
        as.data.frame() %>% 
        tibble::as_tibble( rownames = "Probe_ID" ) %>% 
        dplyr::mutate( Probe_ID = Probe_ID %>% 
                         stringr::str_replace_all("\\.+", "_") %>% 
                         stringr::str_remove("_+$") ) %>% 
        dplyr::select( -col )
      
      safe_write( x = ctl_tib, file = ctl_csv, vb=vb )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # readr::write_rds( ret_dat, out_rds, compress = "gz" )
    
    out_cnt <- safe_write( x = ret_dat, file = out_rds, type = "rds",
                           done = TRUE, write_spec = TRUE, append = FALSE,
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_cnt <- ret_dat$Probe_ID %>% length()
    # ret_key <- glue::glue("final-ret-tib")
    # ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
    #                       vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( vb >= vt ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

prefix_to_sdf = function( prefix,
                          
                          platform = NULL,
                          manifest = NULL,
                          controls = NULL,
                          
                          out_dir,
                          run_tag,
                          
                          reload     = 0,
                          reload_min = 2,
                          reload_pre = NULL,
                          
                          ret_data   = FALSE,
                          parallel   = FALSE,
                          
                          vb=0, vt=3, tc=1, tt=NULL,
                          fun_tag='prefix_to_sdf' )
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
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_rds <- file.path( out_dir, paste(out_tag, 'rds', sep='.') )
  # qdf_rds <- file.path( out_dir, paste(out_tag, 'qdf.rds', sep='.') )
  # ctl_csv <- file.path( out_dir, paste(out_tag, 'ctl.csv.gz', sep='.') )
  # pdf_csv <- file.path( out_dir, paste(out_tag, 'pdf.csv.gz', sep='.') )
  # sigs_csv <- file.path( out_dir, paste(out_tag, 'idat.sigs.csv.gz', sep='.') )
  # info_csv <- file.path( out_dir, paste(out_tag, 'idat.info.csv.gz', sep='.') )
  beg_txt <- paste(out_rds, 'start.txt', sep='.')
  end_txt <- paste(out_rds, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1 )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_rds, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1, tt=tt )
  
  if ( reload >= reload_min && is_valid ) {
    if ( p3 ) cat(glue::glue("{mssg} Reloading rds = '{out_rds}'...{RET}"))
    return( readr::read_rds( out_rds ) )
  }
  # if ( reload >= reload_min && is_valid )
  #   return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
  #                      vb=vb,vt=vt+1,tc=tc+1, tt=tt ) )
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters (vb={vb},vt={vt},tc={tc})::{RET}"))
    cat(glue::glue("{mssg}       prefix = '{prefix}'.{RET}"))
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
    cat(glue::glue("{mssg}      out_rds = '{out_rds}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, aux_csv, out_rds, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  prefix_dir <- base::dirname(prefix)
  
  errs_mssg <- glue::glue("Directory ='{prefix_dir}' does not exist")
  if ( !dir.exists( prefix_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    ret_tib <- sesame::readIDATpair( prefix.path = prefix, 
                                     platform = platform, 
                                     manifest = manifest,
                                     controls = controls,
                                     verbose  = p2 )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_rds, type = "rds", 
                           done = TRUE, write_spec = TRUE, append = FALSE, 
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1, tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1, tt=tt )
  })
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file

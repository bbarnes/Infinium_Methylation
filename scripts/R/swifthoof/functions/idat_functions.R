
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Basic IDAT Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Load Core Packages::
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("illuminaio",quietly=TRUE) ))

# Load Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel",quietly=TRUE) ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Idat File I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prefixesToChipTib = function(prefixes) {
  basecodes <- names(prefixes)
  
  tib <- NULL
  for (basecode in basecodes) {
    sentrix_codes <- basecode %>% stringr::str_split(pattern='_', n=2, simplify=TRUE)
    tib <- tib %>% dplyr::bind_rows(
      tibble::tibble(barcode=sentrix_codes[1],
                     poscode=sentrix_codes[2],
                     basecode=basecode,
                     path=prefixes[[basecode]])
    )
  }
  
  tib
}

prefixToIdat = function(prefix, load=FALSE, save=FALSE, light=TRUE,
                        csv=NULL, ssh=NULL,
                        gzip=TRUE, validate=TRUE, 
                        verbose=0,vt=6,tc=1,tt=NULL) {
  funcTag <- 'prefixToIdat'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Loading prefix={prefix}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    if ( load && !is.null(csv) && file.exists(csv) ) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading CSV={csv}.{RET}"))
      ret_tib <- suppressMessages(suppressWarnings( readr::read_csv(csv) ))
    } else {
      grn_idat <- loadIdat(prefix, 'Grn', gzip=gzip, 
                           verbose=verbose,vt=vt+4,tc=tc+1,tt=tt)
      
      red_idat <- loadIdat(prefix, 'Red', gzip=gzip, 
                           verbose=verbose,vt=vt+4,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Extract Signal Data::
      grn_sig <- getIdatSignalTib(grn_idat, channel='Grn', 
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      red_sig <- getIdatSignalTib(red_idat, channel='Red', 
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      if (verbose>=vt+4)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining Grn/Red data...{RET}"))
      ret_tib <- dplyr::full_join(grn_sig, red_sig, by="Address")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
      if (verbose>=vt+4)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done.Joining Grn/Red data{RET}{RET}"))
      
      grn_ann <- dplyr::bind_cols(
        getIdatBarcodeTib(grn_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatFormatTib(grn_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatTimeStampTib(grn_idat, method='Decoding', 
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatTimeStampTib(grn_idat, method='Extract', 
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      )
      
      if ( validate ) {
        red_ann <- dplyr::bind_cols(
          getIdatBarcodeTib(red_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
          getIdatFormatTib(red_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
          getIdatTimeStampTib(red_idat, method='Decoding', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
          getIdatTimeStampTib(red_idat, method='Extract', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        )
        
        unique_cnt <- dplyr::bind_rows(grn_ann, red_ann) %>% dplyr::distinct() %>% base::nrow()
        stopifnot(unique_cnt==1)
        if (verbose>=vt+1) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Passed idat pair validation!{RET}"))
      }
      
      if (save && !is.null(csv)) {
        outDir <- base::dirname(csv)
        if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
        
        if (file.exists(csv)) unlink(csv)
        if (verbose>=vt+1) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Signal Table CSV={csv}.{RET}"))
        if (  light ) readr::write_csv( dplyr::select( ret_tib, Address), csv )
        if ( !light ) readr::write_csv( ret_tib, csv )
      }
      
      if ( save && !is.null(ssh) ) {
        outDir <- base::dirname(ssh)
        if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
        
        if (file.exists(ssh)) unlink(ssh)
        if (verbose>=vt+1) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Sample Sheet CSV={ssh}.{RET}"))
        readr::write_csv( grn_ann, ssh )
      }
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

loadIdat = function(prefix, col, gzip=TRUE, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadIdat'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading prefix={prefix}.{RET}"))
  
  ret_cnt <- 0
  ret_dat <- NULL
  stime <- system.time({
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Ensure Idats are Gzipped::
    idat_file <- paste0(prefix,"_",col,".idat")
    if (file.exists(idat_file)) {
      if (gzip) {
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Compressing file={idat_file}.{RET}"))
        
        system(paste0("gzip ",idat_file))
        idat_file <- paste0(prefix,"_",col,".idat.gz")
      }
    } else {
      idat_file <- paste0(prefix,"_",col,".idat.gz")
    }
    if (!file.exists(idat_file)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: idat_file={idat_file} does NOT exist!!!{RET}{RET}"))
      return(ret_dat)
    }
    stopifnot(file.exists(idat_file))
    
    if (verbose>=vt+4) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading IDAT={idat_file}.{RET}"))
    
    ret_dat <- suppressMessages(suppressWarnings( illuminaio::readIDAT(idat_file) ))
    ret_cnt <- ret_dat %>% names() %>% length()
    if (verbose>=vt+5) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} names(ret_dat({ret_cnt}))={RET}"))
      ret_dat %>% names() %>% print()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Idat Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getIdatSignalTib = function(idat, channel, del='_', verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'getIdatSignalTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; channel={channel}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    datTag   <- 'sig'
    meanName <- paste('Raw', channel,datTag, sep=del)
    sdName   <- paste('SD',  channel,datTag, sep=del)
    beadName <- paste('Bead',channel, sep=del)
    
    new_cnames <- c('Address',meanName,sdName,beadName)
    ret_tib <- idat$Quants %>% tibble::as_tibble(rownames="Address") %>%
      dplyr::mutate(Address=as.integer(Address)) %>%
      purrr::set_names(new_cnames)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+5) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

getIdatBarcodeTib = function(idat, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'getIdatBarcodeTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    sentrixName <- paste(idat$Barcode,idat$Unknowns$MostlyA, sep='_')
    rowcol_df <- idat$Unknowns$MostlyA %>% stringr::str_remove("^R") %>% stringr::str_split('C', simplify=TRUE, n=2) %>% as.numeric()
    ret_tib <- tibble::tibble(Sentrix_Name=sentrixName,
                              Sentrix_Barcode=idat$Barcode,
                              Sentrix_Poscode=idat$Unknowns$MostlyA,
                              Sentrix_Row=as.integer(rowcol_df[1]),
                              Sentrix_Col=as.integer(rowcol_df[2]) )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+5) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

idat_to_tib = function( idat,
                        
                        out_dir,
                        run_tag,
                        
                        reload     = 0,
                        reload_min = 2,
                        reload_pre = NULL,
                        
                        ret_data   = FALSE,
                        parallel   = FALSE,
                        
                        vb=0, vt=3, tc=1, tt=NULL,
                        fun_tag='idat_to_tib') {
  
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
    
    print(idat)
    
    chipType <- idat$ChipType
    if ( chipType == 'SLIDE.15028542.24X1X1' ) chipFormat  <- '24x1'
    else if ( chipType == 'BeadChip 8x5' )     chipFormat  <- '8x1'
    else if ( chipType == 'BeadChip 12x8' )    chipFormat  <- '12x1'
    else if ( chipType == 'BeadChip 24x1x4' )  chipFormat  <- '24x1'
    else if ( chipType == 'BeadChip 24x1x2' )  chipFormat  <- '24x1'
    else if ( chipType == 'Beadchip 24x1x2' )  chipFormat  <- '24x1'
    else if ( chipType == 'BeadChip 48x4' )    chipFormat  <- '48x1'
    else 
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Unrecognized ChipType='{chipType}'!{RET2}"))
    
    ret_tib <- tibble::tibble( 'Chip_Type'   = chipType,
                               'Chip_Format' = chipFormat )
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+5) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
    
    warn_mssg <- glue::glue("WARN_MESSAGE")
    if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
    if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    wflag <- FALSE
    
    errs_mssg <- glue::glue("ERROR_MESSAGE")
    if ( is.null(args) || length(args) == 0 ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
    
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

getIdatFormatTib = function( idat, 
                             verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'getIdatFormatTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    chipType <- idat$ChipType
    if ( chipType == 'SLIDE.15028542.24X1X1' ) chipFormat  <- '24x1'
    else if ( chipType == 'BeadChip 8x5' )     chipFormat  <- '8x1'
    else if ( chipType == 'BeadChip 12x8' )    chipFormat  <- '12x1'
    else if ( chipType == 'BeadChip 24x1x4' )  chipFormat  <- '24x1'
    else if ( chipType == 'BeadChip 24x1x2' )  chipFormat  <- '24x1'
    else if ( chipType == 'Beadchip 24x1x2' )  chipFormat  <- '24x1'
    else if ( chipType == 'BeadChip 48x4' )    chipFormat  <- '48x1'
    else 
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Unrecognized ChipType='{chipType}'!{RET2}"))
    
    ret_tib <- tibble::tibble( 'Chip_Type'   = chipType,
                               'Chip_Format' = chipFormat )
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+5) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

getIdatTimeStampTib = function( idat, 
                                method='Extract', 
                                sherlockID='sherlockID', 
                                order='latest', 
                                verbose=0,vt=6,tc=1,tt=NULL) {
  funcTag <- 'getIdatTimeStampTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    time_tib  <- idat$RunInfo %>% tibble::as_tibble()
    method_idxs <- grep(method, time_tib$BlockType)
    machine_idxs <- grep(sherlockID, time_tib$BlockPars)
    stopifnot(length(method_idxs)>0)
    if (verbose>=vt+10) print(method_idxs)
    
    stopifnot(order=='latest')
    if (order=='latest')  {
      metod_idx <- method_idxs %>% tail(n=1)
      machine_idx <- machine_idxs %>% tail(n=1)
    } else {
      metod_idx <- method_idxs %>% head(n=1)
      machine_idx <- machine_idxs %>% head(n=1)
    }
    mach_vec <- time_tib$BlockPars[machine_idx] %>% str_split('\\|', simplify=TRUE)
    mach_var <- mach_vec[2] %>% stringr::str_split('=', simplify=TRUE)
    
    # print(time_tib$RunTime[metod_idx])
    time_tib$RunTime[metod_idx] <- time_tib$RunTime[metod_idx] %>% 
      stringr::str_remove(' [APM]+$')
    date_str <- time_tib$RunTime[metod_idx] %>% 
      as.POSIXct(format="%m/%d/%Y %H:%M:%S")
    if (is.na(date_str)) date_str <- time_tib$RunTime[metod_idx] %>% 
      as.POSIXct(format="%d/%m/%Y %H:%M:%S")
    
    # Another weird failure case found in Evonik data:
    #  PASS = "06/15/2022 8:40:25 PM"
    #  FAIL = "20-Jul-22 1:31:37 PM"
    if (is.na(date_str)) {
      if (verbose>=vt+5)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} time_tib$RunTime[metod_idx]=({time_tib$RunTime[metod_idx]}){RET}"))
      
      # date_str <- time_tib$RunTime[metod_idx] %>% 
      #   as.POSIXct( format="%d-%m-%Y %H:%M:%S" )
      date_str <- time_tib$RunTime[metod_idx] %>% 
        as.POSIXct( format="%d-%B-%Y %H:%M:%S",
                    tryFormats = c("%Y-%m-%d %H:%M:%OS",
                                   "%Y/%m/%d %H:%M:%OS",
                                   "%Y-%m-%d %H:%M",
                                   "%Y/%m/%d %H:%M",
                                   "%Y-%m-%d",
                                   "%Y/%m/%d" ) )
      if (verbose>=vt+5)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} date_str=({date_str}){RET}"))
    }
    
    name_str <- paste('Iscan',method, sep='_')
    
    mach_key <- mach_var[1]
    mach_val <- mach_var[2]
    
    # print(date_str)
    ret_tib <- date_str %>% 
      stringr::str_split('[\\-\\/ \\:]') %>% 
      BiocGenerics::unlist() %>% 
      purrr::set_names('Year','Mon','Day','Hour','Min','Sec') %>% 
      tibble::enframe() %>% spread(name, value) %>% dplyr::mutate_all(.funs = (as.integer)) %>%
      tibble::add_column('Date' = date_str) %>%
      tibble::add_column(!!mach_key := !!mach_val) %>%
      dplyr::select(!!mach_key, Date, Year, Mon, Day, Hour, Min, Sec) %>%
      purrr::set_names(paste(name_str, names(.), sep='_'))
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+5) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}


# End of file


suppressWarnings(suppressPackageStartupMessages(require("R6")) )
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )

COM <- ","
TAB <- "\t"
RET <- "\n"

timeTracker <- R6Class("timeTracker", list(
  time  = NULL,
  begTime = NULL,
  runTime = 0,
  file_vec = NULL,
  
  vt = 4,
  verbose = 0,

  initialize = function(verbose=0,vt=8) {
    stime <- system.time({
      self$time     <- NULL
      self$begTime  <- Sys.time()
      self$runTime  <- 0
      self$file_vec <- NULL
      
      self$vt <- vt
      self$verbose <- verbose
    })
    self$addTime(stime, 'init')
  },
  
  
  addTime = function(time, key) {
    self$time <- self$time %>% dplyr::bind_rows( self$formatSysTime(time,key) )
    
    lastIdx <- self$time %>% base::nrow()
    self$runTime <- difftime(self$time$DateStamp[lastIdx], self$begTime)
    if (self$verbose>=self$vt) self$print()
    
    invisible(self)
  },
  
  addFile = function(path) {
    self$file_vec <- c(self$file_vec, path) %>% unique()

    invisible(self)
  },

  addSummary = function() {
    if (!is.null(self$time)) {
      self$time <- self$time %>% dplyr::bind_rows(
        self$time %>%
          base::replace(is.na(.), 0) %>%
          dplyr::summarise_if(is.numeric, list(sum)) %>% 
          dplyr::mutate(Method='Summary', DateStamp=Sys.time()) %>% 
          dplyr::select(Method, everything())
      )
    }
    lastIdx <- self$time %>% base::nrow()
    self$runTime <- difftime(self$time$DateStamp[lastIdx], 
                             self$time$DateStamp[1])
    # self$runTime <- difftime(self$time$DateStamp[lastIdx], self$begTime)
    # if (self$verbose>=self$vt) self$print()

    invisible(self)
  },
  
  formatSysTime = function(stime,key) {
    funcTag <- 'formatTime'
    
    ctime <- stime %>% tibble::enframe() %>% 
      dplyr::mutate(Method=key, name=stringr::str_replace_all(name,'\\.','_')) %>% 
      tidyr::spread('name', 'value') %>% dplyr::mutate(DateStamp=Sys.time())
    
    ctime
  },
  
  print = function(...) {
    self$addSummary()
    cat(glue::glue("Beg Time: {self$begTime}{RET}"))
    cat(glue::glue("Run Time: {self$runTime}{RET}"))
    if (!is.null(self$time)) {
      self$time %>% print()
      if (base::nrow(self$time)>12) self$time %>% tail %>% print()
    }
  }
))

# End of file


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Align All Probe Sequence:: BSMAP
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
  base::require("readr", quietly = TRUE) ))
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Static Structures & Variables::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

slim_bsp_cols <- readr::cols(
  Probe_ID = readr::col_character(),
  Aln_Seq  = readr::col_character(),
  
  Tag  = readr::col_character(),
  Chr  = readr::col_character(),
  Beg  = readr::col_integer(),
  Strand_BS  = readr::col_character(),
  
  Ins_Len  = readr::col_integer(),
  Ref_Seq  = readr::col_character(),
  Mis_Cnt  = readr::col_integer(),
  Mis_Str  = readr::col_character()
)

all_bsp_cols <- readr::cols(
  Probe_ID = readr::col_character(),
  Aln_Seq  = readr::col_character(),
  Qual     = readr::col_character(),
  
  Tag  = readr::col_character(),
  Chr  = readr::col_character(),
  Beg  = readr::col_integer(),
  Strand_BS  = readr::col_character(),
  
  Ins_Len  = readr::col_integer(),
  Ref_Seq  = readr::col_character(),
  Mis_Cnt  = readr::col_integer(),
  Mis_Str  = readr::col_character()
)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Align All Probe Sequence:: BSMAP
#                           Main Workflow Driver
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_mapping_workflow = function(ref_fas = NULL, 
                                ref_tib = NULL,
                                can_tib,
                                
                                map_tib = NULL,
                                map_dir = NULL,
                                map_tsv = NULL,
                                map_col = NULL,
                                
                                bsp_dir = NULL,
                                bsp_exe,
                                bsp_opt = "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R",
                                
                                # Look into these...
                                join_key  = NULL,
                                join_type = "inner",
                                sort      = TRUE, 
                                slim      = TRUE,
                                full      = FALSE, 
                                merge     = TRUE,
                                add_chr   = FALSE,
                                
                                out_dir,
                                run_tag,
                                pre_tag = NULL,
                                
                                reload     = 0,
                                reload_min = 2,
                                ret_data   = FALSE,
                                parallel   = FALSE,
                                use_rcpp   = TRUE,
                                run_join   = TRUE,
                                
                                vb=0, vt=3, tc=1, tt=NULL,
                                fun_tag='bsp_mapping_workflow')
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
  top_csv <- file.path( out_dir, paste(out_tag, 'top.csv.gz', sep='.'))
  aux_csv <- file.path( out_dir, paste(out_tag, 'aux.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ),
                                out_dir = out_dir,
                                vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+3) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      ref_fas = '{ref_fas}'.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      map_dir = '{map_dir}'.{RET}"))
    cat(glue::glue("{mssg}      map_tsv = '{map_tsv}'.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      top_csv = '{top_csv}'.{RET}"))
    cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} BSMAP Parameters::{RET}"))
    cat(glue::glue("{mssg}      bsp_dir = '{bsp_dir}'.{RET}"))
    cat(glue::glue("{mssg}      bsp_exe = '{bsp_exe}'.{RET}"))
    cat(glue::glue("{mssg}      bsp_opt = '{bsp_opt}'.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}         sort = '{sort}'.{RET}"))
    cat(glue::glue("{mssg}         slim = '{slim}'.{RET}"))
    cat(glue::glue("{mssg}         full = '{full}'.{RET}"))
    cat(glue::glue("{mssg}        merge = '{merge}'.{RET}"))
    cat(glue::glue("{mssg}      add_chr = '{add_chr}'.{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     use_rcpp = '{use_rcpp}'.{RET}"))
    cat(glue::glue("{mssg}     run_join = '{run_join}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, top_csv, aux_csv, out_csv, end_txt) )
  
  f_time <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    #
    # TBD:: Should Pre-Load 
    #   i.e. Add CGN Map Data:: file.path(map_dir,map_tsv)
    #
    # TBD:: Need a wrapper for:
    #   - run_bsmap()
    #   - load_bsmap()
    #
    #   - join_bsmap() or rename annotate_bsmap???
    #   - join_cgn_map() with pre-loaded table
    #     - test run times:: loading vs. sorting
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Pre-processing BSMAP Inputs:: 
    #                       Reference/Candidate Fasta Files
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    can_ids_key <- "Probe_ID"
    can_prb_key <- "Probe_Seq"
    
    can_ids_sym <- rlang::sym(can_ids_key)
    can_prb_sym <- rlang::sym(can_prb_key)
    
    ref_tib <- format_reference_bsmap_tib( fas = ref_fas,
                                           tib = ref_tib,
                                           
                                           dna_vec = "dna",
                                           bsc_vec  = "N", 
                                           frs_vec  = "F", 
                                           
                                           out_dir = out_dir,
                                           run_tag = run_tag,
                                           pre_tag = c(pre_tag, beg_txt),
                                           
                                           ret_data   = ret_data,
                                           reload     = reload,
                                           reload_min = reload_min-1,
                                           
                                           vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    if (ret_data) ret_dat$ref_tib <- ref_tib
    
    can_fas <- format_candidate_bsmap_fas( tib = can_tib,
                                           
                                           ids_key = can_ids_key,
                                           prb_key = can_prb_key,
                                           
                                           out_dir = out_dir,
                                           run_tag = run_tag,
                                           pre_tag = c(pre_tag, beg_txt),
                                           
                                           reload     = reload,
                                           reload_min = reload_min-1,
                                           ret_data   = ret_data,
                                           use_rcpp   = TRUE,
                                           
                                           vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    if (ret_data) ret_dat$can_fas <- can_fas
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #         Align Candidate Probes Against all Reference Fastas:: BSMAP
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # TBD:: Make this parallel on the cluster...
    #
    ref_cnt <- 0
    ref_cnt <- ref_tib %>% base::nrow()
    if (vb >=vt+1)
      cat(glue::glue("{mssg} Will align candidate fasta={can_fas} against ",
                     "{ref_cnt} reference fasta file(s)...{RET}"))
    
    if (FALSE && parallel) {
      if (vb >=vt+1)
        cat(glue::glue("{mssg} Will align in parallel mode...{RET}"))
      
      
      
    } else {
      if (vb >=vt+1)
        cat(glue::glue("{mssg} Will align in linear mode...{RET}"))
      
      for (ref_idx in c(1:ref_cnt)) {
        
        cur_tib <- bsmap_workflow(
          ref_fas = ref_tib$Path[ref_idx], can_fas = can_fas, can_tib = can_tib,
          out_dir = out_dir, out_sub = ref_tib$Out_Sub_Path[ref_idx], 
          map_tib = map_tib, map_dir = map_dir, map_tsv = map_tsv, map_col = map_col,
          slim = slim, sort = sort, add_chr = add_chr,
          bsp_dir = bsp_dir, bsp_exe = bsp_exe, bsp_opt = bsp_opt,
          run_tag = run_tag, pre_tag = c(pre_tag, beg_txt),
          ret_data = ret_data, reload = reload, reload_min = reload_min-1,
          parallel = parallel, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
        
        ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
      }
      
    }
    if (ret_data) ret_dat$bsp_dat <- ret_tib
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                Join Bsmap::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( run_join ) {
      
      if ( !is.null(map_tsv) && !is.null(map_dir) && 
           !file.exists(map_tsv) && dir.exists(map_dir) )
        map_tsv <- file.path( map_dir, map_tsv )
      
      ret_tib <- join_bsmap( bsp_tib = ret_tib,
                             can_tib = can_tib, 
                             
                             map_tib = map_tib,
                             map_tsv = map_tsv,
                             map_col = map_col,
                             
                             out_dir = out_dir,
                             run_tag = run_tag,
                             pre_tag = pre_tag,
                             
                             uc         = TRUE,
                             ret_data   = ret_data,
                             reload     = reload,
                             reload_min = reload_min-1,
                             parallel   = parallel,
                             use_rcpp   = use_rcpp,
                             
                             vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Alignment Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( run_join ) {
      # Multiple Hit Summary::
      bsp_hit_sum <- ret_tib %>% 
        dplyr::group_by( !!can_ids_sym ) %>%
        dplyr::summarise(Count=n(), .groups="drop") %>%
        dplyr::group_by(Count) %>% 
        dplyr::summarise(His_Count=n(), .groups="drop")
      hit_key <- glue::glue("bsp_hit_sum")
      hit_key <- glue::glue("bsp-hit-sum")
      hit_cnt <- print_tib( bsp_hit_sum, fun_tag = fun_tag, name = hit_key, 
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )
      
      # Top Ranked Offenders::
      top_hit_sum <- ret_tib %>% 
        dplyr::group_by( !!can_ids_sym ) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>%
        dplyr::filter(Count!=1) %>%
        dplyr::arrange(-Count)
      # print(top_hit_sum, n=base::nrow(top_hit_sum))
      
      # Write Summaries::
      sum_cnt <- safe_write(x = bsp_hit_sum, file = sum_csv, fun_tag = fun_tag, 
                            write_spec = FALSE, append = FALSE, done = FALSE, 
                            vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      top_cnt <- safe_write(x = top_hit_sum, file = top_csv, fun_tag = fun_tag,
                            write_spec = FALSE, append = FALSE, done = FALSE, 
                            vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      
      # NOTE: Auxiliary file is more or less already written by bsp_worflow()
      #   function. That's actually a better spot for it anyways...
      #
      # - Split Core/Auxiliary
      # top_tib <- ret_tib %>% dplyr::select( dplyr::any_of( fun_var$top_vec ) )
      # aux_tib <- ret_tib %>% dplyr::select(!names( top_tib ) )
      # 
      # if (ret_data) ret_dat$aux <- aux_tib
      # if (ret_data) ret_dat$top <- top_tib
      
      # Write Auxiliary Probe Order Output::
      # aux_cnt <- safe_write(x = aux_tib, file = aux_csv, fun_tag = fun_tag,
      #                       vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      
      # TBD:: Need to subset the full output to only the important columns...
      #
      # Write Top Probe Order Fields::
      # out_cnt <- safe_write(x = top_tib, file = out_csv, fun_tag = fun_tag,
      #                       vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                               Write Data::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # safe_write( sum_csv, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      # safe_write( aux_csv, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      # safe_write( out_csv, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      # out_cnt <- safe_write( x = ret_tib, file = out_csv, type = "csv", 
      #                        done = TRUE, write_spec = TRUE, append = FALSE, 
      #                        fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                               Print Summary::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_key <- glue::glue("final-ret-tib")
      ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    }
  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  if (ret_data) return(ret_dat)
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                             BSMAP Methods::
#
#  run_bsmap()
#  load_bsmap()
#  format_reference_bsmap_tib()
#  format_candidate_bsmap_fas()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Pre-processing BSMAP Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_reference_bsmap_tib = function(fas = NULL,
                                      tib = NULL,
                                      
                                      dna_vec = c("dna"),
                                      bsc_vec = c("N"),
                                      frs_vec = c("F"),
                                      
                                      out_dir,
                                      run_tag,
                                      pre_tag = NULL,
                                      
                                      ret_data   = FALSE,
                                      reload     = 0,
                                      reload_min = 2,
                                      
                                      vb=0, vt=3, tc=1, tt=NULL,
                                      fun_tag='format_reference_bsmap_tib')
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
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+3) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          fas = '{fas}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      bsc_vec = '{bsc_vec}'.{RET}"))
    cat(glue::glue("{mssg}      frs_vec = '{frs_vec}'.{RET}"))
    cat(glue::glue("{mssg}      dna_vec = '{dna_vec}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_key <- glue::glue("initial-ret-tib")
  ret_cnt <- print_tib( tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("File file='{fas}' does not exist")
  if ( !file.exists( fas) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)

  unlink( c(out_csv, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  e_time <- 0
  f_time <- 0
  f_time <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Format Genomes:: Direct User Input
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( !is.null(fas) && file.exists(fas) )
      ret_tib <- 
        tibble::tibble(
          Path=fas,
          Genome_Base_Name = base::basename(Path) %>% 
            stringr::str_remove(".fa.gz"),
          Genome_Version = "Unknown",
          Molecule_Type  = "Whole_Genome",
          Out_Sub_Path = file.path(
            # out_dir,
            Genome_Version,
            Molecule_Type,
            Genome_Base_Name) 
        )
    
    ret_key <- glue::glue("Post-reference-genome-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Filter and Format Genomes:: imGenomes
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( !is.null(tib) )
      ret_tib <- tib %>%
      dplyr::filter( Genome_Alphabet   %in% dna_vec ) %>%
      dplyr::filter( Genome_Strand_BSC %in% bsc_vec ) %>%
      dplyr::filter( Genome_Strand_FR  %in% frs_vec ) %>%
      dplyr::mutate(
        Out_Sub_Path=file.path(
          # out_dir,
          Genome_Version,
          Molecule_Type,
          Genome_Base_Name)
      ) %>% dplyr::bind_rows(ret_tib)
    
    ret_key <- glue::glue("Post-imGenomes-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Validate Each Fasta File Exists::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    for (cur_fas in ret_tib$Path)
      if (is.null(cur_fas) || !file.exists(cur_fas)) {
        stop(glue::glue("{RET}{mssg} ERROR: Failed to find fasta={cur_fas}! ",
                        "user input or reference genome tibble!{RET2}"))
        print(ret_tib)
        return(NULL)
      }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv,
                           done = TRUE, write_spec = FALSE, append = FALSE,
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
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
  
  ret_tib
}

format_candidate_bsmap_fas = function(tib,
                                      fas = NULL,
                                      
                                      ids_key,
                                      prb_key,
                                      
                                      out_dir,
                                      run_tag,
                                      pre_tag = NULL,
                                      
                                      ret_data   = FALSE,
                                      reload     = 0,
                                      reload_min = 2,
                                      use_rcpp   = TRUE,
                                      
                                      vb=0, vt=3, tc=1, tt=NULL,
                                      fun_tag='format_candidate_bsmap_fas')
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
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  out_fas <- file.path( out_dir, paste(out_tag, 'prb.fas.gz', sep='.') )
  if (!is.null(fas)) out_fas <- fas
  beg_txt <- paste(out_fas, 'start.txt', sep='.')
  end_txt <- paste(out_fas, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_fas, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid ) return( out_fas )
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+3) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     use_rcpp = '{use_rcpp}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_fas = '{out_fas}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_key <- glue::glue("initial-ret-tib")
  ret_cnt <- print_tib( tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  
  if ( !ids_key %in% names(tib) ) {
    fail_mssg <- glue::glue("Failed to find ids_key={ids_key} in tib")
    stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
    return(NULL)
  }
  if ( !prb_key %in% names(tib) ) {
    fail_mssg <- glue::glue("Failed to find prb_key={prb_key} in tib")
    stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
    return(NULL)
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
  
  f_time <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    ids_sym <- rlang::sym(ids_key)
    prb_sym <- rlang::sym(prb_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Set:: De-Methylated Align Sequence
    #                      Extract:: Fasta Line from tib
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (use_rcpp)
      fas_vec <- tib %>%
      # dplyr::select( dplyr::all_of( !!ids_sym, !!prb_sym ) ) %>%
      dplyr::mutate( !!prb_sym := deM_cpp( seqs = !!prb_sym, uc = TRUE ) ) %>%
      dplyr::mutate( fas_line = paste0(">",!!ids_sym,"\n", !!prb_sym ) ) %>%
      dplyr::pull( fas_line )
    else
      fas_vec <- tib %>%
      # dplyr::select( dplyr::all_of( !!ids_sym, !!prb_sym ) ) %>%
      dplyr::mutate( !!prb_sym := deMs( seqs = !!prb_sym, uc = TRUE ) ) %>%
      dplyr::mutate( fas_line = paste0(">",!!ids_sym,"\n", !!prb_sym ) ) %>%
      dplyr::pull( fas_line )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Write:: Fasta File
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = fas_vec, file = out_fas, type = "line", 
                           done = TRUE, write_spec = FALSE, append = FALSE, 
                           compress = NULL, permissions = NULL, 
                           fun_tag = fun_tag, vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib({fun_tag})")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt )
  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  out_fas
}

bsmap_workflow = function(ref_fas,
                          can_fas,
                          can_tib,
                          
                          map_tib = NULL,
                          map_dir = NULL,
                          map_tsv = NULL,
                          map_col = NULL,
                          
                          bsp_exe,
                          bsp_dir = NULL,
                          bsp_opt = "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R",
                          
                          sort = TRUE,
                          slim = TRUE,
                          add_chr = FALSE,
                          
                          out_dir,
                          out_sub,
                          run_tag,
                          pre_tag = NULL,
                          
                          ret_data   = FALSE,
                          reload     = 0,
                          reload_min = 2,
                          
                          parallel   = FALSE,
                          use_rcpp   = TRUE,
                          
                          vb=0, vt=3, tc=1, tt=NULL,
                          fun_tag='bsmap_workflow')
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
  
  out_dir <- file.path( out_dir, fun_tag, out_sub )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
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
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} bsmap Parameters::{RET}"))
    cat(glue::glue("{mssg}      ref_fas = '{ref_fas}'.{RET}"))
    cat(glue::glue("{mssg}      can_fas = '{can_fas}'.{RET}"))
    cat(glue::glue("{mssg}      bsp_exe = '{bsp_exe}'.{RET}"))
    cat(glue::glue("{mssg}      bsp_dir = '{bsp_dir}'.{RET}"))
    cat(glue::glue("{mssg}      bsp_opt = '{bsp_opt}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}         sort = '{sort}'.{RET}"))
    cat(glue::glue("{mssg}         slim = '{slim}'.{RET}"))
    cat(glue::glue("{mssg}      add_chr = '{add_chr}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}     use_rcpp = '{use_rcpp}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_key <- glue::glue("initial-can-tib")
  ret_cnt <- print_tib( can_tib, fun_tag = fun_tag, name = ret_key, 
                        vb=vb,vt=vt+6,tc=tc+1,tt=tt )
  
  
  if ( is.null(ref_fas) ) {
    fail_mssg <- glue::glue("Reference fasta variable ref_fas is NULL")
    stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
    return(NULL)
  }
  if ( !file.exists(ref_fas) ) {
    fail_mssg <- glue::glue("Failed to find ref_fas = {ref_fas}")
    stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
    return(NULL)
  }
  
  if ( is.null(can_fas) ) {
    fail_mssg <- glue::glue("Candidate fasta variable can_fas is NULL")
    stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
    return(NULL)
  }
  if ( !file.exists(can_fas) ) {
    fail_mssg <- glue::glue("Failed to find can_fas = {can_fas}")
    stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
    return(NULL)
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
  
  f_time <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                Run Bsmap::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Force Bsp Output fields to match bsp_var results::
    #
    out_bsp <- NULL
    out_bsp <- run_bsmap(ref_fas = ref_fas,
                         can_fas = can_fas,
                         out_dir = out_dir,
                         
                         slim    = slim, 
                         
                         bsp_dir = bsp_dir,
                         bsp_exe = bsp_exe,
                         bsp_opt = bsp_opt,
                         
                         run_tag = run_tag,
                         pre_tag = pre_tag,
                         
                         ret_data   = ret_data,
                         reload     = reload,
                         reload_min = reload_min-1,
                         
                         vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                Load Bsmap::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- NULL
    ret_tib <- load_bsmap( file = out_bsp,
                           sort = sort,
                           slim = slim,
                           add_chr = add_chr,
                           
                           out_dir = out_dir,
                           run_tag = run_tag,
                           pre_tag = pre_tag,
                           
                           ret_data   = ret_data,
                           reload     = reload,
                           reload_min = reload_min-1,
                           
                           vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = ret_tib, file = out_csv,
                           done = TRUE, write_spec = TRUE, append = FALSE,
                           fun_tag = fun_tag, vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

write_bsmap_fasta = function( tib,
                              
                              # Example: c("Loci_ID","Strand_FR","Strand_CO", 'Allele',...)
                              key_col,
                              key_sep = "_",
                              
                              # Example: c("Probe_Seq_A", "Probe_Seq_B") [ Vector ]
                              seq_col,
                              
                              # Example: c("A", "B")
                              seq_sep = NULL,
                              
                              unique = TRUE,
                              sort = FALSE,
                              
                              out_dir,
                              run_tag,
                              
                              reload     = 0,
                              reload_min = 2,
                              reload_pre = NULL,
                              
                              ret_data   = FALSE,
                              parallel   = FALSE,
                              
                              vb=0, vt=3, tc=1, tt=NULL,
                              fun_tag='write_bsmap_fasta') {
  
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
  out_fas <- file.path( out_dir, paste(out_tag, 'fas.gz', sep='.') )
  beg_txt <- paste(out_fas, 'start.txt', sep='.')
  end_txt <- paste(out_fas, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(reload_pre, beg_txt, out_fas, end_txt ), 
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid )
    return( safe_read( out_fas, fun_tag = fun_tag, head = "Reloading",
                       vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if ( vb >= vt ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
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
    cat(glue::glue("{mssg}      out_fas = '{out_fas}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(sum_csv, aux_csv, out_fas, end_txt) )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # errs_mssg <- glue::glue("Seq_Sep must be same length as Seq_Vec if provided!")
  # if ( length(seq_sep) != length(seq_col) ) eflag <- TRUE
  # if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  # if ( eflag ) return( NULL )
  # 
  # if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Validate Inputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    key_cnt <- key_col %>% length()
    int_cnt <- intersect( key_col, names(tib) ) %>% length()
    
    errs_mssg <- glue::glue("Failed to match all keys in tib {int_cnt} != {key_cnt}")
    if ( key_cnt != int_cnt ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
    
    seq_cnt <- seq_col %>% length()
    int_cnt <- intersect( seq_col, names(tib) ) %>% length()
    
    errs_mssg <- glue::glue("Failed to match all seqs in tib {int_cnt} != {seq_cnt}")
    if ( seq_cnt != int_cnt ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return( NULL )
    
    if ( vb >= vt+2 ) cat(glue::glue("{mssg}{TAB} Found all keys and seqs!{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Extract Inputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    let_col <- base::LETTERS[ c(1:length(seq_col) ) ]
    mul_col <- c( "Key", let_col )
    sig_col <- c( "Key", "Probe" )
    
    ret_tib <- tib %>%
      dplyr::select( dplyr::all_of( c(key_col, seq_col) ) ) %>%
      tidyr::unite( Key, dplyr::all_of(key_col), sep = key_sep )
    
    if ( length(seq_col) == 1 ) ret_tib <- ret_tib %>% purrr::set_names( sig_col )
    
    if ( length(seq_col) > 1 ) ret_tib <- ret_tib %>%
      purrr::set_names( mul_col ) %>%
      tidyr::pivot_longer( cols = dplyr::all_of( let_col ),
                           names_to = c("Allele"),
                           values_to = "Probe",
                           # names_sep = "_", 
                           values_drop_na = TRUE ) %>%
      tidyr::unite(Key, c("Key","Allele"), sep=key_sep )
    
    ret_tib <- ret_tib %>% dplyr::mutate(
      Probe = mutate_seqs_cpp( seqs = Probe, m="X", uc=TRUE, vb=vb, vt=vt+4 ) )
    
    if ( unique ) ret_tib <- ret_tib %>% dplyr::distinct()
    if ( sort ) ret_tib <- ret_tib %>% dplyr::arrange( Key )
    
    fas_vec <- ret_tib %>%
      dplyr::mutate( Fasta_Str = paste0(">",Key,"\n",Probe) ) %>%
      dplyr::pull( Fasta_Str )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- safe_write( x = fas_vec, file = out_fas, type = "line", 
                           done = TRUE, write_spec = FALSE, append = FALSE, 
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
  
  if ( ret_data == "file" ) return( out_fas )
  
  # ret_tib
}

run_bsmap = function( ref_fas,
                      can_fas,
                      
                      bsp_exe,
                      bsp_dir = NULL,
                      bsp_opt = "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R",
                      slim    = TRUE,
                      dry_run = FALSE,
                      
                      out_dir,
                      run_tag,
                      pre_tag = NULL,
                      add_tag = NULL,
                      
                      ret_data   = FALSE,
                      reload     = 0,
                      reload_min = 2,
                      
                      vb=0, vt=3, tc=1, tt=NULL,
                      fun_tag='run_bsmap')
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

  if (slim) add_tag <- ".slim"
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste(run_tag, fun_tag, sep='.')
  out_tag <- paste0(out_tag, add_tag)
  bsp_ssh <- file.path( out_dir, paste(out_tag, 'aln.bsp.sh', sep='.') )
  out_bsp <- file.path( out_dir, paste(out_tag, 'aln.bsp', sep='.') )
  out_tsv <- file.path( out_dir, paste(out_tag, 'aln.bsp.tsv.gz', sep='.') )
  beg_txt <- paste(out_tsv, 'start.txt', sep='.')
  end_txt <- paste(out_tsv, 'done.txt', sep='.')
  
  safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  is_valid <- valid_time_stamp( c(pre_tag, beg_txt, bsp_ssh,out_tsv, end_txt ),
                                out_dir = out_dir,
                                vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  if ( reload >= reload_min && is_valid ) return( out_tsv )
  # return( safe_read( out_tsv, fun_tag = fun_tag, head = "Reloading",
  #                    vb=vb,vt=vt+4,tc=tc+1,tt=tt ) )
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+3) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} bsmap Parameters::{RET}"))
    cat(glue::glue("{mssg}      ref_fas = '{ref_fas}'.{RET}"))
    cat(glue::glue("{mssg}      can_fas = '{can_fas}'.{RET}"))
    cat(glue::glue("{mssg}      bsp_exe = '{bsp_exe}'.{RET}"))
    cat(glue::glue("{mssg}      bsp_dir = '{bsp_dir}'.{RET}"))
    cat(glue::glue("{mssg}      bsp_opt = '{bsp_opt}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}        slim = '{slim}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      out_bsp = '{out_bsp}'.{RET}"))
    cat(glue::glue("{mssg}      out_tsv = '{out_tsv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  unlink( c(beg_txt, out_tsv, end_txt) )

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("File ref_fas='{ref_fas}' does not exist")
  if ( !file.exists( ref_fas) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)

  errs_mssg <- glue::glue("File can_fas='{can_fas}' does not exist")
  if ( !file.exists( can_fas) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  f_time <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    add_cmd <- ""
    if (slim) add_cmd <- "cut -f 1,2,4-11 | "
    
    if (!is.null(bsp_dir)) {
      bsp_exe_test <- file.path( bsp_dir, base::basename(bsp_exe) )
      if (file.exists(bsp_exe_test)) bsp_exe_test
    }
    
    # TBD:: Clean up this missing bsp_exe for docker...
    #
    if ( ( is.null(bsp_exe) || !file.exists(bsp_exe ) ) &&
         (!is.null(bsp_dir) &&  dir.exists(bsp_dir  ) ) )
      bsp_exe <- file.path(bsp_dir, bsp_exe)
    
    if ( is.null(bsp_exe) || !file.exists(bsp_exe) ) {
      cat(glue::glue("{mssg} Warning: Unable to locate bsp_exe={bsp_exe}. ",
                     "Will try docker version: ",
                     "2.90 '/repo/bsmap-2.90/bsmap'.{RET2}"))
      bsp_exe <- '/repo/bsmap-2.90/bsmap'
    }
    
    # if (!file.exists(bsp_exe)) {
    #   cat(glue::glue("{mssg} Warning: Unable to locate bsp_exe={bsp_exe}. ",
    #                  "Will try docker version: '/repo/BSMAPz/bsmapz'.{RET2}"))
    #   bsp_exe <- '/repo/BSMAPz/bsmapz'
    # }
    
    if (!file.exists(bsp_exe)) {
      fail_mssg <- glue::glue("bsp_exe='{bsp_exe}' does NOT exist")
      stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
      return(NULL)
    }
    
    bsp_cmd <- glue::glue(
      "{bsp_exe} -a {can_fas} -d {ref_fas} {bsp_opt} -o {out_bsp}{RET}",
      "cat {out_bsp} | {add_cmd} gzip -c -> {out_tsv}{RET}",
      "rm -f {out_bsp}{RET}" )
    
    if (vb >= vt)
      cat(glue::glue("{mssg} Writing BSMAP shell = {bsp_ssh}...{RET2}"))
    
    out_cnt <- safe_write( x = bsp_cmd, type = "line", file = bsp_ssh, 
                           done = TRUE, write_spec = FALSE, append = FALSE, 
                           permissions = "0777", fun_tag = fun_tag,
                           vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    if (vb >= vt)
      cat(glue::glue("{mssg} Running; CMD = '{bsp_cmd}'...{RET}"))
    
    sys_ret <- 0
    if ( !dry_run ) sys_ret <- base::system(bsp_ssh)
    
    if ( sys_ret != 0 ) {
      fail_mssg <- glue::glue("sys_ret({sys_ret}) != 0!")
      stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
      return(NULL)
    } else {
      if (vb >= vt)
        cat(glue::glue("{mssg} BSMAP Completed Succesfully!{RET2}"))
      
      sys_ret <- base::system( glue::glue("touch {end_txt}") )
    }
    
    ret_cnt <- sys_ret
  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  out_tsv
}

load_bsmap = function( file,
                       sort = TRUE,
                       slim = TRUE,
                       allele = NULL,
                       ale_key = "Allele",
                       write_out = FALSE,
                       
                       add_col = NULL,
                       add_tib = NULL,
                       add_key = NULL,
                       
                       add_cnt = FALSE,
                       add_chr = FALSE,
                       parse_id = FALSE,
                       
                       out_dir,
                       run_tag,
                       
                       reload     = 0,
                       reload_min = 2,
                       reload_pre = NULL,
                       
                       ret_data   = FALSE,
                       parallel   = FALSE,
                       
                       vb=0, vt=3, tc=1, tt=NULL,
                       fun_tag='load_bsmap')
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
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}         file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}         sort = '{sort}.'{RET}"))
    cat(glue::glue("{mssg}         slim = '{slim}'.{RET}"))
    cat(glue::glue("{mssg}       allele = '{allele}'.{RET}"))
    cat(glue::glue("{mssg}      ale_key = '{ale_key}'.{RET}"))
    cat(glue::glue("{mssg}    write_out = '{write_out}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}      add_col = '{add_col}'.{RET}"))
    cat(glue::glue("{mssg}      add_tib = '{add_tib}'.{RET}"))
    cat(glue::glue("{mssg}      add_key = '{add_key}'.{RET}"))
    cat(glue::glue("{mssg}      add_cnt = '{add_cnt}'.{RET}"))
    cat(glue::glue("{mssg}      add_chr = '{add_chr}'.{RET}"))
    cat(glue::glue("{mssg}     parse_id = '{parse_id}'.{RET}"))
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
  
  errs_mssg <- glue::glue("File file='{file}' does not exist")
  if ( !file.exists( file) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if (slim) {
    bsp_cols <- slim_bsp_cols
    chr_key  <- names(bsp_cols$cols)[4]
    beg_key  <- names(bsp_cols$cols)[5]
    
  } else {
    bsp_cols <- all_cols
    chr_key  <- names(bsp_cols$cols)[5]
    beg_key  <- names(bsp_cols$cols)[6]
  }
  bsp_col_len <- length(bsp_cols$cols)

  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    chr_sym <- rlang::sym(chr_key)
    beg_sym <- rlang::sym(beg_key)
    ale_sym <- rlang::sym(ale_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Load Raw BSP Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (vb >= vt)
      cat(glue::glue("{mssg} Loading BSP={file}...{RET}"))
    
    ret_tib <- safe_read( file = file, type = "tsv",
                          clean = TRUE, 
                          use_spec = TRUE, spec_col = bsp_cols, 
                          has_head = FALSE, write_spec = TRUE, 
                          vb=vb,vt=vt+3,tc=tc+1,tt=tt ) %>%
      dplyr::select( -dplyr::any_of( c("Bsp_Qual") ) )
    
    ret_key <- glue::glue("after-loading-before-chr-modification")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    ret_tib <- ret_tib %>% dplyr::mutate(
      !!chr_sym := stringr::str_remove(!!chr_sym,'^chr') ) %>%
      clean_tib( fun_tag = fun_tag, name = "cleaning-bsp-tib",
                 vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Chromosome ID Modification::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # if (add_chr) ret_tib <- ret_tib %>% 
    #   dplyr::mutate( !!chr_sym:=paste0('chr', !!chr_sym,'^chr' ) )
    if (add_chr) ret_tib <- ret_tib %>%
      dplyr::mutate( !!chr_sym:=paste0('chr', !!chr_sym ) )
    
    ret_key <- glue::glue("after-chr-fix-and-clean-tibble")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    if (is.null(ret_tib) || length(ret_tib)==0) {
      fail_mssg <- glue::glue("bsmap results are NULL")
      stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
      return(NULL)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Probe ID Parsing Modification::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if ( parse_id ) {
      if (vb >= vt)
        cat(glue::glue("{mssg} Parsing ID...{RET}"))
      
      p1_cnt <- ret_tib %>% dplyr::filter( 
        stringr::str_detect( Probe_ID, "_[TBN][CON][12]$" ) ) %>% 
        base::nrow()
      
      p2_cnt <- ret_tib %>% dplyr::filter( 
        stringr::str_detect( Probe_ID, "_[TBN][CON][12]_[A-Z]$" ) ) %>% 
        base::nrow()
      
      if ( p1_cnt > 0 && p2_cnt == 0 ) {
        if (vb >= vt)
          cat(glue::glue("{mssg} Parsing V1: p1={p1_cnt}, p2={p2_cnt}...{RET}"))
        
        ret_tib <- ret_tib %>%
          # SRD_Str = Probe_ID %>% stringr::str_remove("^[^_]+") %>%
          tidyr::separate( Probe_ID, into=c("Probe_ID","SRD_Str"), sep='_' ) %>%
          dplyr::mutate( Ilmn_ID = paste( Probe_ID,SRD_Str, sep='_') ) %>%
          tidyr::separate( 
            SRD_Str, into=c( "Strand_TB", "Strand_CO", "Infinium_Design" ),
            sep = c(1,2,3), convert = TRUE ) %>%
          dplyr::mutate( !!ale_sym := 'N' )
        
      } else if ( p1_cnt == 0 && p2_cnt > 0 ) {
        if (vb >= vt)
          cat(glue::glue("{mssg} Parsing V2: p1={p1_cnt}, p2={p2_cnt}...{RET}"))
        
        ret_tib <- ret_tib %>%
          tidyr::separate( Probe_ID,
                           into = c( "Probe_ID","SRD_Str", ale_key ),
                           sep='_', remove = FALSE ) %>%
          dplyr::mutate( Ilmn_ID = paste( Probe_ID,SRD_Str, sep='_') ) %>%
          tidyr::separate(
            SRD_Str, into=c( "Strand_TB", "Strand_CO", "Infinium_Design" ),
            sep = c(1,2,3), convert = TRUE )
        
      } else {
        fail_mssg <- glue::glue("Failed to parse Probe_ID")
        stop(glue::glue("{errs} {fail_mssg}!{errs} Exiting...{RET2}"))
        return(NULL)
      }
      ret_key <- glue::glue("after-Parse-ID")
      ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                            vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      ret_tib <- ret_tib %>%
        dplyr::mutate( 
          Probe_Type = Probe_ID %>% stringr::str_sub( 1,2 ),
          Probe_Type = dplyr::case_when(
            Probe_Type %>% stringr::str_starts("[0-9]") ~ "cg",
            TRUE ~ Probe_Type
          )
        ) %>%
        dplyr::select( Ilmn_ID, Probe_ID, Probe_Type, Infinium_Design, 
                       !!ale_sym, Strand_TB, Strand_CO, dplyr::everything() )
      
      ret_key <- glue::glue("after-Adding-Probe_Type")
      ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                            vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                         CpG Offset Modification::
      #                                    &
      #                      Alignment Count Modification::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if ( !is.null(cpg_offset_tib) ) {
        if (vb >= vt)
          cat(glue::glue("{mssg} Adding CpG Postion Offset...{RET}"))
        
        off_tib <- cpg_offset_tib %>% 
          dplyr::rename( Strand_FR_Tbl = Strand_FR, Strand_CO_Tbl = Strand_CO )
        
        ret_tib <- ret_tib %>%
          dplyr::inner_join( off_tib,
                             by=c("Probe_Type", "Infinium_Design", "Strand_BS"),
                             suffix=c("_ord","_bsp") ) %>%
          dplyr::mutate( CpG_Pos = Beg + CG_Offset ) %>%
          dplyr::select( Ilmn_ID, Probe_ID, Probe_Type, Infinium_Design, 
                         !!ale_sym, Strand_TB, Strand_CO, 
                         Aln_Seq, Tag, Chr, Beg, CpG_Pos, dplyr::everything() )
        
        ret_key <- glue::glue("after-CpG-Position-Offset")
        ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                              vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      }
      
      if ( add_cnt )
        ret_tib <- ret_tib %>% 
        dplyr::group_by( Ilmn_ID,Aln_Seq ) %>%
        dplyr::mutate( Map_Count = dplyr::row_number() ) %>%
        dplyr::ungroup()
      
    } else {
      if ( add_cnt )
        ret_tib <- ret_tib %>% 
          dplyr::group_by( Probe_ID,Aln_Seq ) %>%
          dplyr::mutate( Map_Count = dplyr::row_number() ) %>%
          dplyr::ungroup()
    }
    
    ret_key <- glue::glue("after Adding CpG Postion Offset")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Additional Modifications::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Add Allele if provided::
    if ( !is.null(allele) ) ret_tib <- dplyr::mutate(ret_tib, !!ale_sym := allele)
    
    # Add additional columns if provided::
    if ( !is.null(add_col) ) ret_tib <- dplyr::bind_cols( ret_tib, add_col )
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(!!chr_sym, !!beg_sym)
    
    ret_key <- glue::glue("after Adding Additional Modifications ",
                          "[ Allele, Additiona-Cols, Genomic-Soriing ]")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key,
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    out_cnt <- 0
    if ( write_out ) out_cnt <- safe_write( 
      x = ret_tib, file = out_csv, type = "csv", 
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

match_tri_tag_improbe = function( tib,
                                  
                                  # Example: file.path( opt$top_path, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/prbs-B46-split3" )
                                  tag_dir,
                                  
                                  # Example:: ".b46.probe_cgn-table.raw.sorted.tsv.gz$"
                                  tag_str,
                                  
                                  # Example: file.path( opt$top_path, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/min/GRCh37.chr-pos-srd.slim.cgn-sorted.txt.gz" )
                                  imp_tsv = NULL,
                                  
                                  # Example: "ProbeSeq"
                                  seq_key,
                                  
                                  # Example: "Infinium_Design"
                                  des_key,
                                  
                                  out_dir,
                                  run_tag,
                                  
                                  reload     = 0,
                                  reload_min = 2,
                                  reload_pre = NULL,
                                  
                                  ret_data   = FALSE,
                                  parallel   = FALSE,
                                  
                                  vb=0, vt=3, tc=1, tt=NULL,
                                  fun_tag='match_tri_tag_improbe') {
  
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
    cat(glue::glue("{mssg}      tag_dir = '{tag_dir}'.{RET}"))
    cat(glue::glue("{mssg}      tag_str = '{tag_str}'.{RET}"))
    cat(glue::glue("{mssg}      imp_tsv = '{imp_tsv}'.{RET}"))
    cat(glue::glue("{mssg}      seq_key = '{seq_key}'.{RET}"))
    cat(glue::glue("{mssg}      des_key = '{des_key}'.{RET}"))
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
  
  tri_cols <- readr::cols(
    Probe_46M  = readr::col_character(),
    Last_Base  = readr::col_character(),
    Cgn_Num    = readr::col_integer(),
    Platform_Bits   = readr::col_character(),
    Degen_UM_Idx    = readr::col_integer(),
    Seq_Tag_Srd_Idx = readr::col_integer()
  )
  
  imp_cols <- readr::cols(
    Chromosome = readr::col_character(),
    Coordinate = readr::col_integer(),
    Cgn_Num    = readr::col_integer(),
    Top_Strand = readr::col_character()
  )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Validate Inputs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  errs_mssg <- glue::glue("Tag Directory='{tag_dir}' does not exist")
  if ( !dir.exists( tag_dir) ) eflag <- TRUE
  if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
  if ( eflag ) return(NULL)
  
  if ( vb >= vt+2 ) cat(glue::glue("{mssg} Inputs are valid!{RET}"))
  
  ftime <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    seq_sym = rlang::sym( seq_key )
    des_sym = rlang::sym( des_key )
    
    ret_dat <- tib %>%
      dplyr::mutate(
        ProbeAln = mutate_seqs_cpp( seqs = !!seq_sym, m="X", uc=TRUE, vb=vb, vt=8 ),
        
        # Prefix Sequence (default length = 3)
        Probe_49P = dplyr::case_when(
          !!des_sym == 1 ~ stringr::str_sub( ProbeAln, 1,3 ),
          !!des_sym == 2 ~ stringr::str_sub( ProbeAln, 2,4 ),
          TRUE ~ NA_character_ ),
        
        # Suffix Sequence (default length = 49 - 3)
        Probe_49S = dplyr::case_when(
          !!des_sym == 1 ~ stringr::str_sub( ProbeAln, 4,49 ),
          !!des_sym == 2 ~ stringr::str_sub( ProbeAln, 5,50 ),
          TRUE ~ NA_character_ )
      ) %>% 
      dplyr::select( Address, Probe_49P, Probe_49S ) %>%
      split( .$Probe_49P )
    
    tri_files <- file_list( path = tag_dir,
                            pattern = tag_str,
                            suffix  = tag_str,
                            vb=vb,vt=vt+3,tc=tc+1,tt=tt )
    
    tri_seqs <- base::names( ret_dat )
    for (tri in tri_seqs ) {
      
      tri_tsv <- tri_files[[tri]]
      
      if ( !file.exists(tri_tsv) ) {
        warn_mssg <- glue::glue("Failed to find({tri}): tri_tsv='{tri_tsv}'. Skipping...")
        if ( is.null(args) || length(args) == 0 ) wflag <- TRUE
        if ( wflag ) cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
        wflag <- FALSE
        next
      }
      
      tri_tib <- readr::read_tsv( file = tri_tsv, 
                                  col_names = names(tri_cols$cols), 
                                  col_types = tri_cols )
      
      ret_dat[[tri]] <- ret_dat[[tri]] %>%
        dplyr::left_join( tri_tib, by=c("Probe_49S"="Probe_46M") ) %>% 
        dplyr::left_join( seq_tag_srd_tib, by=c("Seq_Tag_Srd_Idx"="Tag_Idx") ) %>%
        # dplyr::add_count( Address, name="Address_Cgn_Tag_Count") %>%
        dplyr::select( -Probe_49P, -Probe_49S, -Seq_Tag_Srd_Idx )
      
      if ( vb >= vt+2 ) cat( glue::glue("{mssg} Done! tri={tri}...{RET2}") )
    }
    ret_tib  <- ret_dat %>% dplyr::bind_rows()
    
    if ( !is.null(imp_tsv) && file.exists(imp_tsv) ) {
      
      if ( vb >= vt+2 ) 
        cat(glue::glue("{mssg} Loading improbe [pos <=> cgn]: '{imp_tsv}'!{RET}"))
      
      # imp_tib <- safe_read( file = imp_tsv, spec_col = imp_cols, 
      #                       spec_prefernce = "col", 
      #                       vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      cgn_pos_tib <- readr::read_tsv( file = imp_tsv, 
                                      col_names = names(imp_cols$cols), 
                                      col_types = imp_cols )
      
      ret_tib <- ret_tib %>% dplyr::left_join( cgn_pos_tib, by="Cgn_Num")
    }
    
    ret_tib <- ret_tib %>%
      dplyr::add_count( Address, name="Address_Cgn_Tag_Count")
    
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

# TBD:: Re-qrite clean::
#
load_cgn_map = function(file,
                        dir,
                        cols = NULL,
                        names = NULL,
                        
                        vb=0, vt=3, tc=1, tt=NULL,
                        fun_tag='load_cgn_map')
{
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  warn <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  errs <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb >= vt+3) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}    file = '{file}'.{RET}"))
    cat(glue::glue("{mssg}   names = '{names}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  if (is.null(cols))
    cols <- readr::cols(
      Chr = readr::col_character(),
      Pos = readr::col_integer(),
      Cgn = readr::col_integer(),
      Top = readr::col_character()
    )
  
  # Update Names::
  if (!is.null(names) && length(names) == length(cols$cols))
    names(cols$cols) <- names
  
  if ( vt >= vt+4 ) {
    cat(glue::glue("{mssg} Using cols::{RET}"))
    print(cols)
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  e_time <- 0
  f_time <- 0
  f_time <- base::system.time({
    
    errs_mssg <- glue::glue("File file='{file}' does not exist")
    if ( !file.exists( file) ) eflag <- TRUE
    if ( eflag ) stop(glue::glue("{errs} {errs_mssg}!{errs} Exiting...{RET2}"))
    if ( eflag ) return(NULL)
    
    if ( vb >= vt+1 ) 
      cat(glue::glue("{mssg} Loading map_tsv = '{file}'...{RET2}"))
    
    ret_tib <- safe_read( file = file,
                          use_spec = TRUE, 
                          spec_col = cols,
                          has_head = FALSE,
                          write_spec = TRUE,
                          spec_prefernce = "col",
                          # over_write = TRUE, 
                          clean = TRUE,
                          fun_tag = fun_tag,
                          
                          vb=vb,vt=vt+4,tc=tc+1,tt=tt )
    
    if (vb >= vt+1) 
      cat(glue::glue("{mssg} Loading map_tsv = '{file}'...{RET2}"))
    
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
  
  ret_tib
}

# End of file

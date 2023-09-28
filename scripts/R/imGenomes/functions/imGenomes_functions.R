
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           imGenomes Functions::
#            Similar to iGenomes, but for infinium methylation 
#                          Hence the lower-case im....
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


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Locate and Load all 
#                             Normal, SNP and 
#                   Pre-Bisulfite Converted Genomes (bsc)::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_imGenomes_table = function(dir, 
                                genome_build,
                                ret_list = FALSE,
                                load_chroms = FALSE,
                                vb=0, vt=3, tc=1, tt=NULL,
                                fun_tag='load_imGenomes_table') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}            dir={dir}.{RET}"))
    cat(glue::glue("{mssg}       ret_list={ret_list}.{RET}"))
    cat(glue::glue("{mssg}    load_chroms={load_chroms}.{RET}"))
    cat(glue::glue("{mssg}   genome_build={genome_build}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  e_time <- 0
  f_time <- 0
  f_time <- base::system.time({
    
    gen_dir <- file.path(dir, genome_build,"Sequence/WholeGenomeFasta")
    gen_pattern <- paste0(genome_build,".genome.*.fa.gz$")
    
    if (vb >= vt) {
      cat(glue::glue("{mssg} Using Genome Build={genome_build}...{RET}"))
      cat(glue::glue("{mssg} Using Genome Build Directory={gen_dir}...{RET}"))
      cat(glue::glue("{mssg} Searching with Pattern={gen_pattern}...{RET}"))
    }
    fas_list <- list.files(gen_dir, pattern=gen_pattern, full.names=TRUE)
    fas_count <- fas_list %>% length()
    
    if (vb >= vt) {
      cat(glue::glue("{mssg} Found {fas_count} Fasta File(s)={RET}"))
      print(fas_list)
    }
    
    ret_tib <- fas_list %>% 
      tibble::as_tibble() %>% 
      purrr::set_names(c("Path")) %>%
      dplyr::mutate(Genome_Base_Name=base::basename(Path) %>% 
                      stringr::str_remove(".fa.gz$")) %>%
      dplyr::mutate(Unq_ID=stringr::str_remove(Genome_Base_Name, 
                                               paste(genome_build,"genome.",sep=".")), 
                    Unq_ID=stringr::str_replace(Unq_ID,"^F",
                                                paste(genome_build,"NCBI.dna.F", sep='.') ), 
                    Unq_ID=stringr::str_replace(Unq_ID,"^R",
                                                paste(genome_build,"NCBI.dna.R", sep='.') ),
                    Unq_ID=stringr::str_replace(
                      Unq_ID, paste(genome_build,"genome$", sep='.'), 
                      paste(genome_build,"NCBI.dna.FCN", sep='.') ),
                    Unq_ID=stringr::str_replace(Unq_ID,"dbSNP-151.iupac$",
                                                "dbSNP-151.iupac.FCN"), 
                    Unq_ID=stringr::str_replace(Unq_ID,"dbSNP-151.iupac",
                                                paste0(genome_build,".dbSNP-151.snp")) ) # %>%
    # dplyr::select(-Genome_Base_Name)
    
    ret_key <- glue::glue("mid-formatting-imGenomes")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+1,tc=tc+1,tt=tt )
    
    
    ret_tib <- ret_tib %>% 
      tidyr::separate(Unq_ID,
                      into=c("Genome_Version","Genome_Source",
                             "Genome_Alphabet","Genome_Key"), sep="\\.") %>%
      tidyr::separate(Genome_Key, 
                      into=c("Genome_Strand_FR","Genome_Strand_CO",
                             "Genome_Strand_BSC"), sep=c(1,2), 
                      remove=FALSE) %>%
      dplyr::mutate(Genome_Key=paste(Genome_Key,Genome_Alphabet, sep="_")) %>%
      dplyr::mutate(
        Genome_Alphabet_Int=dplyr::case_when(
          Genome_Alphabet=="dna" ~ 0,
          Genome_Alphabet=="snp" ~ 1,
          TRUE ~ 3,
        ) %>% as.integer(),
        Genome_Strand_BSC_Int=dplyr::case_when(
          Genome_Strand_BSC=="N" ~ 0,
          Genome_Strand_BSC=="U" ~ 1,
          Genome_Strand_BSC=="M" ~ 2,
          Genome_Strand_BSC=="D" ~ 3,
          TRUE ~ 4,
        ) %>% as.integer(),
        Genome_Strand_CO_Int=dplyr::case_when(
          Genome_Strand_CO=="C" ~ 0,
          Genome_Strand_CO=="O" ~ 1,
          TRUE ~ 2,
        ) %>% as.integer(),
        Genome_Strand_FR_Int=dplyr::case_when(
          Genome_Strand_FR=="F" ~ 0,
          Genome_Strand_FR=="R" ~ 1,
          TRUE ~ 2
        ) %>% as.integer()
      ) %>% dplyr::arrange(Genome_Alphabet_Int,
                           Genome_Strand_BSC_Int,
                           Genome_Strand_CO_Int,
                           Genome_Strand_FR_Int) %>%
      dplyr::mutate(Molecule_Type="Whole_Genome")
    
    chr_tib <- NULL
    chr_dir <- NULL
    
    if (load_chroms) {
      
      chr_dir <- file.path(dir, genome_build,"Sequence/Chromosomes")
      if (!is.null(chr_dir) && dir.exists(chr_dir)) {
        cat(glue::glue("{mssg} Will load individual chromosomes={chr_dir}...{RET}"))
        chr_pattern <- paste0(".fa.gz$")
        
        fas_list <- list.files(chr_dir, pattern=chr_pattern, full.names=TRUE)
        fas_count <- fas_list %>% length()
        
        chr_tib <- fas_list %>% 
          tibble::as_tibble() %>% 
          purrr::set_names(c("Path")) %>%
          dplyr::mutate(Genome_Base_Name=base::basename(Path) %>% 
                          stringr::str_remove(".fa.gz$")) %>% 
          dplyr::mutate(Chromosome=Genome_Base_Name,
                        Genome_Version=opt$genome_build,
                        Genome_Source="NCBI",
                        Genome_Alphabet="dna",Genome_Key="FCN_dna",
                        Genome_Strand_FR="F", Genome_Strand_CO="CO",
                        Genome_Strand_BSC="N",Molecule_Type="Chrom")
        
        chr_tib <- chr_tib %>%
          dplyr::mutate(
            
            Chr_Str_Len=Chromosome %>% stringr::str_remove("^chr") %>% stringr::str_remove("^[0-9XYM]+") %>% stringr::str_length(),
            Is_Full_Chromosome=dplyr::case_when(
              Chr_Str_Len==0 ~ TRUE,
              TRUE ~ FALSE),
            Molecule_Type=dplyr::case_when(
              Chr_Str_Len>0 ~ "Partial_Contigs",
              TRUE ~ "Whole_Chrom")
          )
        
        cat(glue::glue("{mssg} Found {fas_count} Fasta File(s)={RET}"))
        print(fas_list)
      }
      
    }
    ret_tib <- dplyr::bind_rows(ret_tib, chr_tib)
    
    if (ret_list) ret_dat <- ret_tib %>% split(f=ret_tib$Genome_Key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Print Summary::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_key <- glue::glue("final-ret-tib")
    ret_cnt <- print_tib( ret_tib, fun_tag = fun_tag, name = ret_key, 
                          vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  })
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Split Pos-Cgn Maps:: by Chromosome
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

split_pos_map_by_chrom = function(map_dir,
                                  map_tsv,

                                  gen_bld,
                                  imG_tib = NULL,
                                  chr_names = NULL,
                                  
                                  rm_chr = TRUE,
                                  list_id = NULL,
                                  
                                  out_dir,
                                  run_tag,
                                  pre_tag = NULL,
                                  
                                  reload     = 0,
                                  reload_min = 1,
                                  ret_data   = FALSE,
                                  parallel   = FALSE,
                                  
                                  vb=0, vt=3, tc=1, tt=NULL,
                                  fun_tag='split_pos_map_by_chrom') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("{RET}[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("{RET}[{fun_tag}]: ERROR:{tabs}")
  
  out_dir <- file.path( out_dir, fun_tag )
  out_tag <- paste( run_tag, fun_tag, sep='.' )
  sum_csv <- file.path( out_dir, paste(out_tag, 'sum.csv.gz', sep='.') )
  out_csv <- file.path( out_dir, paste(out_tag, 'csv.gz', sep='.') )
  beg_txt <- paste(out_csv, 'start.txt', sep='.')
  end_txt <- paste(out_csv, 'done.txt', sep='.')
  
  # safe_mkdir( out_dir, vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  # is_valid <- valid_time_stamp( c(pre_tag, beg_txt, out_csv, end_txt ), 
  #                               out_dir = out_dir,
  #                               vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  # 
  # if ( reload >= reload_min && is_valid )
  #   return( safe_read( out_csv, fun_tag = fun_tag, head = "Reloading",
  #                      vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  
  if (vb >= vt) cat(glue::glue("{mssg} Starting...{RET2}"))
  if (vb >= vt+3) {
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    # cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}      map_dir = '{map_dir}'.{RET}"))
    cat(glue::glue("{mssg}      map_tsv = '{map_tsv}'.{RET}"))
    cat(glue::glue("{mssg}      gen_bld = '{gen_bld}'.{RET}"))
    cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    cat(glue::glue("{RET}"))
  }
  # unlink( c(sum_csv, aux_csv, out_csv, end_txt) )

  map_tsv <- file.path(map_dir, map_tsv)
  
  # First check if files already exist and meet valid_time_stamps()
  # TBD:: Make these parameters::
  #
  if ( is.null( chr_names ) )
    chr_names <- imG_tib %>% dplyr::filter( 
      Molecule_Type == "Whole_Chrom" & 
        Is_Full_Chromosome == TRUE & 
        Genome_Strand_BSC == "N" & 
        Genome_Alphabet == "dna" &
        Chromosome != "chrM" ) %>% 
    dplyr::pull( Chromosome )
  
  # print(chr_names)
  
  chr_rds <- 
    paste(gen_bld, "chr-pos-srd.slim.pos-sorted", chr_names, "rds", sep = '.')
  chr_rds <- file.path( map_dir, "chromosomes", chr_rds)
  
  chr_cnt <- 0
  chr_cnt <- length(chr_rds)
  if ( vb >= vt+3 ) {
    cat(glue::glue("{mssg} Found {chr_cnt} chromosome split rds files!{RET}"))
    chr_rds %>% head() %>% print()
  }
  
  is_valid <- valid_time_stamp( c( map_tsv, chr_rds ),
                                out_dir = out_dir,
                                vb=vb,vt=vt+4,tc=tc+1,tt=tt )
  
  if ( is_valid ) {
    if ( vb >= vt+1 )
      cat(glue::glue("{mssg} {chr_cnt} chromosomes are valid!{RET}"))

    if (rm_chr) chr_names <- chr_names %>% stringr::str_remove("^chr")
    if ( !is.null( list_id ) && list_id == "chr" ) {
      chr_rds <- as.list( chr_rds )
      names(chr_rds) <- chr_names
    }
    
    return( chr_rds )
  }

  if ( vb >= vt+1 )
    cat(glue::glue("{mssg} Chromosomes are not valid!{RET}",
                   "{mssg} Will build new splits, but code is under",
                   "construction...{RET2}"))
  
  return( NULL )
  
  map_cols <-
    cols(
      Chromosome  = col_character(),
      Coordinate  = col_integer(),
      Bsp_Map_Cgn = col_integer(),
      Bsp_Top_Srd = col_character()
    )
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  e_time <- 0
  f_time <- 0
  f_time <- base::system.time({
    
    if ( !file.exists(beg_txt) )
      sys_ret <- base::system( glue::glue("touch {beg_txt}") )
    
    if (FALSE) {
      
      map_tib <- load_cgn_map( 
        tsv  = map_tsv,
        dir  = map_dir, 
        cols = map_cols, 
        vb=vb,vt=vt+4,tc=tc+1,tt=tt )
      
      prefix <- map_tsv %>% 
        stringr::str_remove(".gz$") %>% 
        stringr::str_remove(".txt")
      
      map_names <- c("Chr", "Pos", "Cgn", "Top")
      map_chrs <- map_tib %>% 
        purrr::set_names(map_names) %>%
        split( f = .$Chr )
      
      for (chr in names(map_chrs)) {
        chr_rds <- paste0(prefix,".chr",chr,".rds")
        chr_rds <- file.path( map_dir, "chromosomes", chr_rds )
        
        cat(glue::glue("{mssg} On chr{chr}, writing RDS = '{chr_rds}'.{RET}"))
        readr::write_rds( x = map_chrs[[chr]], chr_rds )
        cat(glue::glue("{mssg} Done...{RET2}"))
      }
    }
    
    # if (TRUE) {
    # } else if (FALSE) {
    #   warn_mssg <- glue::glue("WARN_MESSAGE")
    #   cat(glue::glue("{warn} {warn_mssg}!{RET2}"))
    # } else {
    #   fail_mssg <- glue::glue("FAIL_MESSAGE")
    #   stop(glue::glue("{error} {fail_mssg}!{error} Exiting...{RET2}"))
    #   return(NULL)
    # }
    
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
  e_time <- as.double(f_time[3]) %>% round(2)
  if (!is.null(tt)) tt$addTime(f_time,fun_tag)
  if (vb >= vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file

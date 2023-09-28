
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Workhorse Functions:: 
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
#                        Local Run Time Defaults::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

workhorse_local_defaults = function(opts, pars,
                                    vb=0, vt=6, tc=1, tt=NULL,
                                    fun_tag='workhorse_local_defaults') {
  
  tabs  <- paste0(rep(TAB, tc), collapse='')
  mssg  <- glue::glue("[{fun_tag}]:{tabs}")
  warn  <- glue::glue("[{fun_tag}]: Warning:{tabs}")
  error <- glue::glue("[{fun_tag}]: ERROR:{tabs}")
  
  if (vb>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (vb>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   fun_tag={fun_tag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  e_time  <- 0
  f_time  <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  if (vb >= vt) {
    cat(glue::glue("[{fun_tag}]: Local    run_mode={pars$run_mode}.{RET}"))
    cat(glue::glue("[{fun_tag}]: Local     top_dir={pars$top_dir}.{RET}"))
  }
  
  opts$Rscript <- "Rscript"
  
  # BSMAP Parameters::
  # opts$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  # opts$bsmap_opt <- "\"-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R\""
  opts$bsmap_opt <- "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
  opts$bsmap_dir <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz"
  opts$bsmap_exe <- "bsmapz"
  
  opts$out_dir <- file.path(pars$top_dir, 'scratch')
  opts$impDir  <- file.path(pars$top_dir, 'data/improbe')
  opts$ann_dir <- file.path(pars$top_dir, 'data/annotation')
  opts$manDir  <- file.path(pars$top_dir, 'data/manifests')
  opts$gen_dir <- file.path(pars$top_dir, 'data/imGenomes/Homo_sapiens/NCBI')
  opts$idatDir <- file.path(pars$top_dir, 'data/idats')
  
  opts$noob <- paste(
    file.path(pars$top_dir, "data/CustomContent/transfer/updated_manifest.csv.gz"),
    sep=',')
  
  opts$clean  <- FALSE
  opts$reload <- 0
  opts$reload <- 1
  
  opts$single   <- FALSE
  opts$cluster  <- FALSE
  opts$parallel <- TRUE
  opts$align_chroms <- FALSE
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Pre-defined local options run_types::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  pars$version_key <- "v"
  opts$version    <- '3'
  
  pars$cus_dir <- 
    file.path(pars$top_dir,"data/workhorse/custom_chips", pars$local_run_type)
  opts$ord_dir <- file.path(pars$cus_dir, "order")
  opts$mat_dir <- file.path(pars$cus_dir, "match")
  opts$aqp_dir <- file.path(pars$cus_dir, "aqpqc")
  opts$aux_dir <- file.path(pars$cus_dir, "auxdb")

  pars$is_pqc <- TRUE
  
  if (FALSE) {
    
  } else if (pars$local_run_type == 'Genknowme') {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                Genknowme::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    opts$genome_build <- 'GRCh37'
    opts$platform <- 'EPIC'
    opts$species  <- 'Homo_sapiens'
    opts$version  <- paste0(pars$version_key,opts$version)
    opts$version  <- 'v1'
    
    opts$ord_csv <- paste(
      "AQP1_NAremoved_GenKnowme_CpG_SNP_order.07082020.csv.gz",
      "AQP2_AP_Genknowme_AQP2_replicate_design_file2.csv.gz",
      sep = ',' )
    
    opts$mat_tsv <- paste(
      "20468029_AQP1_probes.match.gz",
      "20468029_AQP2_probes.match.gz",
      sep = ',' )
    
    if ( pars$is_pqc )
      opts$aqp_tsv <- paste(
        "20042793X371678_A_ProductQC_AP.txt.gz",
        sep = ',' )
    else
      opts$aqp_tsv <- paste(
        "BS0032678_AQP1-AQP.txt.gz",
        "20468029_AQP2_probes.match.gz",
        sep = ',' )
    
  } else if (pars$local_run_type == 'McMaster') {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                McMaster::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    opts$genome_build <- 'GRCh37'
    opts$platform <- 'EPIC'
    opts$species  <- 'Homo_sapiens'
    opts$version  <- paste0(pars$version_key,opts$version)
    opts$version  <- 'v2'

    opts$ord_csv <- paste(
      'McMaster_CpG_DesignFile_v4.csv.gz',
      sep=',' )
    
    opts$mat_tsv <- paste(
      '20532820_probes1.match.gz',
      '20532820_probes2.match.gz',
      sep=',' )
    
    if ( pars$is_pqc )
      opts$aqp_tsv <- paste(
        '20051339_A_ProductQC.txt.gz',
        sep=',')
    else
      opts$aqp_tsv <- paste(
        'BS0033057-AQP1.txt.gz',
        'BS0033090-AQP2.txt.gz',
        sep=',' )

  } else if (pars$local_run_type == 'Chicago') {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                 Chicago::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    opts$genome_build <- 'GRCh38'
    opts$genome_build <- 'GRCh37'
    opts$platform <- 'EPIC'
    opts$species  <- 'Homo_sapiens'
    opts$version  <- paste0(pars$version_key,opts$version)
    opts$version  <- 'v2'
    
    opts$ord_csv <- paste(
      'UofChicago-A_A_Array-CpG-order-FINAL.csv.gz',
      sep=',' )
    
    opts$mat_tsv <- paste(
      '20504790_probes.match.tsv.gz',
      sep=',' )
    
    opts$aqp_tsv <- paste(
      '329922X374054_A_ProductQC.txt.gz',
      sep=',' )

  } else if (pars$local_run_type == 'Epic2' ||
             pars$local_run_type == 'Ewas1') {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Epic2/Ewas1::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

    map_col <- c("Bead_Pool", "MN", "Match_Num", "Bead_Pool_Name", 
                 "Bucket_Name", "Order_Path", "AQP1_Num", "AQP2_Num")
    
    map_csv <- 
      file.path(opts$aux_dir, "Epic2-Ewas_aqp-file-mapping.v.4.csv.gz")
    
    map_tib <- readr::read_csv(map_csv) %>% 
      purrr::set_names(map_col) %>% 
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
                    Order_Path=Order_Path %>% 
                      stringr::str_remove("^.*\\\\") %>% 
                      paste0(opts$ord_dir,'/',.,".gz"),
                    Match_Path=dplyr::case_when(
                      AQP_Num == 2 ~ paste0(opts$mat_dir,"/AQP2-",Match_Num,"_probes.match.gz"),
                      AQP_Num == 1 ~ paste0(opts$mat_dir,"/",Match_Num,"_probes.match.gz"),
                      TRUE ~ NA_character_
                    ),
                    AQP_Path = paste0(opts$aqp_dir,"/",AQP_Name,".txt.gz"),
                    Order_File_Name = base::basename(Order_Path),
                    Match_File_Name = base::basename(Match_Path),
                    AQP_File_Name = base::basename(AQP_Path) ) %>%
      dplyr::select(Bead_Pool:Bucket_Name,AQP_Num,AQP_Name,
                    Order_File_Name,Match_File_Name,AQP_File_Name,
                    Order_Path,Match_Path,AQP_Path) %>%
      dplyr::arrange(AQP_Num, Bead_Pool)
    
    lapply(map_tib$Order_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    lapply(map_tib$Match_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    lapply(map_tib$AQP_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    
    epic_map_tib <- map_tib %>% dplyr::filter(stringr::str_detect(Bucket_Name, "EPIC") )
    ewas_map_tib <- map_tib %>% dplyr::filter(stringr::str_detect(Bucket_Name, "EWAS") )
    
    if (pars$local_run_type == 'Epic2') {
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                                  Epic2::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # EPIC v2:: 
      opts$genome_build <- 'GRCh37'
      opts$platform <- 'EPIC_v2'
      opts$species  <- 'Homo_sapiens'
      opts$version  <- paste0(pars$version_key,opts$version)
      
      opts$ord_dir <- ord_dir
      opts$mat_dir <- mat_dir
      opts$aqp_dir <- aqp_dir
      
      opts$ord_csv <- paste(epic_map_tib$Order_File_Name, collapse = ",")
      opts$mat_tsv <- paste(epic_map_tib$Match_File_Name, collapse = ",")
      opts$aqp_tsv <- paste(epic_map_tib$AQP_File_Name, collapse = ",")
      
    } else if (pars$local_run_type == 'Ewas1') {
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                                  Ewas1::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      opts$genome_build <- 'GRCh37'
      opts$platform <- 'EWAS'
      opts$species  <- 'Homo_sapiens'
      opts$version  <- paste0(pars$version_key,opts$version)
      
      opts$ord_dir <- ord_dir
      opts$mat_dir <- mat_dir
      opts$aqp_dir <- aqp_dir
      
      opts$ord_csv <- paste(ewas_map_tib$Order_File_Name, collapse = ",")
      opts$mat_tsv <- paste(ewas_map_tib$Match_File_Name, collapse = ",")
      opts$aqp_tsv <- paste(ewas_map_tib$AQP_File_Name, collapse = ",")
    }
    
  } else if (pars$local_run_type == 'MM285k') {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                 MM285k::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # This needs to be re-formatted to new dir/csv pair format!!!
    #
    
    # opts$genome_build <- 'GRCm38'
    # opts$genome_build <- 'GRCm10'
    # opts$platform <- 'MM285k'
    # opts$species  <- "Mus_musculus"
    # opts$version  <- 'M0'
    # 
    # opts$gen_dir <- file.path(pars$top_dir, 'data/iGenomes/Mus_musculus/NCBI')
    # pars$aqpDir <- file.path(pars$top_dir, 'data/CustomContent/LifeEpigentics/AQP')
    # opts$ord_csv <- paste(
    #   file.path(pars$aqpDir, 'orders/Mus_musculus.order_BP1.csv.gz'),
    #   file.path(pars$aqpDir, 'orders/Mus_musculus.order_BP2.csv.gz'),
    #   file.path(pars$aqpDir, 'orders/mm10_LEGX_nonCpG_probes.Jan16-2020.order.csv.gz'),
    #   file.path(pars$aqpDir, 'orders/LEGX_SpikeIn_Reorder-All-06052020.order.withHeader.csv.gz'),
    #   sep=',')
    # 
    # opts$mat_csv <- paste(
    #   file.path(pars$aqpDir, 'BP1/20420178_AQP1_LifeEpigen_BP1.txt.gz'),
    #   file.path(pars$aqpDir, 'BP2/20420260_AQP1_LifeEpigen_BP2.txt.gz'),
    #   file.path(pars$aqpDir, 'BP3/20420260_AQP2_LifeEpigen_BP2.txt.gz'),
    #   file.path(pars$aqpDir, 'BP4/20455357_AQP1_LifeEpigen_BP4.txt.gz'),
    #   sep=',')
    # 
    # if (opts$version == 'P0') {
    #   opts$aqp_csv <- paste(
    #     file.path(pars$aqpDir, 'PQC/20042400_A_ProductQC.txt.gz'),
    #     sep=',')
    #   opts$aqpn <- paste(1, sep=",")
    #   
    # } else {
    #   
    #   opts$aqp_csv <- paste(
    #     file.path(pars$aqpDir, 'AQP_Copy/BS0032527-AQP.txt.gz'),
    #     file.path(pars$aqpDir, 'AQP_Copy/BS0032533-AQP.txt.gz'),
    #     file.path(pars$aqpDir, 'AQP_Copy/BS0032545-AQP.txt.gz'),
    #     file.path(pars$aqpDir, 'AQP_Copy/BS0032636-AQP.txt.gz'),
    #     sep=',')
    #   opts$aqpn <- paste(1,1,2,1, sep=",")
    # }
    # opts$bpns <- paste(1,2,2,3, sep=",")
    
  } else if (pars$local_run_type == 'NZT') {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                   NZT::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # This needs to be re-formatted to new dir/csv pair format!!!
    #
    
    # opts$genome_build <- 'GRCh36'
    # opts$genome_build <- 'GRCh38'
    # opts$genome_build <- 'GRCh37'
    # opts$platform <- 'NZT'
    # opts$species  <- 'Homo_sapiens'
    # opts$version  <- 'N0'
    # opts$version  <- 'C4'
    # 
    # pars$combined <- FALSE
    # pars$combined <- TRUE
    # 
    # pars$aqpDir <- file.path(pars$top_dir, 'data/CustomContent/NZT/decode/combined')
    # opts$ord_csv <- paste(
    #   file.path(pars$aqpDir, 'order-12-combined.csv.gz'),
    #   sep=',')
    # 
    # opts$mat_csv <- paste(
    #   file.path(pars$aqpDir, 'match-12-combined.match.tsv.gz'),
    #   sep=',')
    # 
    # opts$aqp_csv <- paste(
    #   file.path(pars$aqpDir, 'aqp-12-combined.txt.gz'),
    #   sep=',')
    # 
    # opts$bpns <- paste(1, sep=",")
    # opts$aqpn <- paste(1, sep=",")
    # 
    # opts$idat <- NULL
    
  } else {
    local_run_type <- 
    stop(glue::glue("{RET}[{fun_tag}]: Unsupported pre-options local ",
                    "type: local_run_type = {pars$local_run_type}!{RET2}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Load improbe Pre-built Database Files::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts$tag_seq_dir <- 
    file.path(opts$imp_dir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-P49-split")
  
  opts$canonical_cgn_dir <- file.path(pars$dat_dir, "manifest/cgnDB")
  opts$canonical_cgn_csv <- "canonical.cgn-top-grp.csv.gz"
  
  opts$bsp_map_dir <- 
    file.path(pars$top_dir, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/min")
  opts$bsp_map_tsv <- paste(opts$genome_build,"chr-pos-srd.slim.pos-sorted.txt.gz", sep='.')
  
  opts$tag_map_dir <- 
    file.path(pars$top_dir, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/min")
  opts$tag_map_tsv <- paste(opts$genome_build,"chr-pos-srd.slim.cgn-sorted.txt.gz", sep='.')
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Load Pre-Built Confirmation Manifests::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Sesame Pre-defined manifests
  opts$sesame_manifest_dat <- "EPIC.hg19.manifest,HM450.hg19.manifest"
  
  # Genome Studio Pre-defined manifests
  genome_manifest_dir <- 
    file.path(pars$top_dir, "data/manifests/methylation/GenomeStudio")
  
  opts$genome_manifest_csv <- paste(
    file.path(genome_manifest_dir, "MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz"),
    file.path(genome_manifest_dir, "HumanMethylation450_15017482_v.1.2.csv.gz"),
    sep = ","
  )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                              Set Run Name::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  if ( is.null(opts$run_name) )
    opts$run_name <- paste( pars$local_run_type, 
                            opts$platform, 
                            opts$version, 
                            opts$genome_build,
                            sep='-' )
  
  # This should be done in program_init()
  #
  # opts$out_dir <- file.path( opts$out_dir, fun_tag, opts$run_name )
  # safe_mkdir( opts$out_dir )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Return Both Options and Parameters::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_dat$opt <- opts
  ret_dat$par <- pars
  
  ret_cnt <- length(opts) + length(pars)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                                  Done::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (vb>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={e_time}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

# End of file

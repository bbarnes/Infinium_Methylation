
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Basic Controls Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )

# Parallel Processing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )

# Illunia IO IDAT Processing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("illuminaio", quietly = TRUE) ) )

# [TBD]: REMOVE OLD STUFF:
# suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
# suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
# suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Controls Formatting:: HSA/New
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_controls_Unk = function(tib,
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_controls_Unk'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib %>% 
      # dplyr::filter(Probe_Type != 'cg') %>% 
      # dplyr::filter(Probe_Type != 'ch') %>% 
      # dplyr::filter(Probe_Type != 'rs') %>% 
      # dplyr::filter(Probe_Type != 'mu') %>% 
      # dplyr::filter(Probe_Type != 'rp') %>% 
      # dplyr::filter(Probe_Type == 'BS' | Probe_Type == 'NO' | Probe_Type == 'NE') %>% 
      dplyr::distinct(M,U, .keep_all=TRUE,Top_Sequence=NA) %>% 
      # dplyr::mutate(SR_Str=FR,CO_Str=CO,Rep_Num=NA,Probe_Source='MUS',NXB_D=NA) %>% 
      dplyr::mutate(Probe_Source='MUS',NXB_D=NA) %>% 
      dplyr::arrange(Probe_Type,Infinium_Design) %>%
      dplyr::mutate(
        Control_Group=dplyr::case_when(
          Probe_Type=='BS' & Infinium_Design==1 ~ "BISULFITE CONVERSION I",
          Probe_Type=='BS' & Infinium_Design==2 ~ "BISULFITE CONVERSION II",
          Probe_Type=='NO' ~ 'NON-POLYMORPHIC',
          Probe_Type=='ne' ~ 'NEGATIVE',
          TRUE ~ NA_character_
        ),
        Control_Group_Str=stringr::str_replace_all(Control_Group,' ','_') %>% 
          stringr::str_replace_all('-','_'),
        DiNuc=dplyr::case_when(
          Probe_Type=='BS' ~ stringr::str_replace(Seq_ID, '^.*-([ACTG][ACTG])-.*$', '\\$1') %>% 
            stringr::str_remove_all('\\\\'),
          Probe_Type=='NO' ~ stringr::str_replace(Seq_ID, '^.*_([ACTG][ACTG])_.*$', '\\$1') %>% 
            stringr::str_remove_all('\\\\'),
          TRUE ~ NA_character_
        ),
        N1=stringr::str_sub(DiNuc, 1,1),
        N2=stringr::str_sub(DiNuc, 2,2),
        Last_BaseA=stringr::str_sub(AlleleA_Probe_Sequence,-1),
        Last_BaseB=stringr::str_sub(AlleleB_Probe_Sequence,-1)
      ) %>% 
      dplyr::group_by(Probe_Type,Infinium_Design) %>%
      dplyr::mutate(
        Row_Idx=dplyr::row_number() + 100,
        Row_Str=Row_Idx,
        Col_Idx=Row_Idx %% (color_len-100)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Control_Name=dplyr::case_when(
          # Probe_Type=='BS' & Infinium_Design==1 & Last_BaseA=='G' ~ paste0(Control_Group_Str,'_U',Row_Str),
          # Probe_Type=='BS' & Infinium_Design==1 & Last_BaseA=='A' ~ paste0(Control_Group_Str,'_C',Row_Str),
          Probe_Type=='BS' & Infinium_Design==1 ~ paste0(Control_Group_Str,'_',Row_Str),
          Probe_Type=='BS' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
          Probe_Type=='NO' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
          Probe_Type=='ne' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
          TRUE ~ NA_character_ ),
        Seq_ID_Org=Seq_ID,
        Seq_ID=Control_Name
      ) %>%
      dplyr::select(Seq_ID,Strand_TB,Strand_TB,Infinium_Design,Probe_Type,U,M, 
                    dplyr::everything())
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

format_controls_New = function(tib,
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_controls_New'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib %>% 
      dplyr::filter(Probe_Type != 'cg') %>% 
      dplyr::filter(Probe_Type != 'ch') %>% 
      dplyr::filter(Probe_Type != 'rs') %>% 
      dplyr::filter(Probe_Type != 'mu') %>% 
      dplyr::filter(Probe_Type != 'rp') %>% 
      dplyr::filter(Probe_Type == 'BS' | Probe_Type == 'NO' | Probe_Type == 'NE') %>% 
      dplyr::distinct(M,U, .keep_all=TRUE,Top_Sequence=NA) %>% 
      dplyr::mutate(SR_Str=FR,CO_Str=CO,Rep_Num=NA,Probe_Source='MUS',NXB_D=NA) %>% 
      dplyr::arrange(Probe_Type,Infinium_Design) %>%
      dplyr::mutate(
        Control_Group=dplyr::case_when(
          Probe_Type=='BS' & Infinium_Design==1 ~ "BISULFITE CONVERSION I",
          Probe_Type=='BS' & Infinium_Design==2 ~ "BISULFITE CONVERSION II",
          Probe_Type=='NO' ~ 'NON-POLYMORPHIC',
          Probe_Type=='ne' ~ 'NEGATIVE',
          TRUE ~ NA_character_
        ),
        Control_Group_Str=stringr::str_replace_all(Control_Group,' ','_') %>% 
          stringr::str_replace_all('-','_'),
        DiNuc=dplyr::case_when(
          Probe_Type=='BS' ~ stringr::str_replace(Seq_ID, '^.*-([ACTG][ACTG])-.*$', '\\$1') %>% 
            stringr::str_remove_all('\\\\'),
          Probe_Type=='NO' ~ stringr::str_replace(Seq_ID, '^.*_([ACTG][ACTG])_.*$', '\\$1') %>% 
            stringr::str_remove_all('\\\\'),
          TRUE ~ NA_character_
        ),
        N1=stringr::str_sub(DiNuc, 1,1),
        N2=stringr::str_sub(DiNuc, 2,2),
        Last_BaseA=stringr::str_sub(AlleleA_Probe_Sequence,-1),
        Last_BaseB=stringr::str_sub(AlleleB_Probe_Sequence,-1)
      ) %>% 
      dplyr::group_by(Probe_Type,Infinium_Design) %>%
      dplyr::mutate(
        Row_Idx=dplyr::row_number() + 100,
        Row_Str=Row_Idx,
        Col_Idx=Row_Idx %% (color_len-100)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Control_Name=dplyr::case_when(
          # Probe_Type=='BS' & Infinium_Design==1 & Last_BaseA=='G' ~ paste0(Control_Group_Str,'_U',Row_Str),
          # Probe_Type=='BS' & Infinium_Design==1 & Last_BaseA=='A' ~ paste0(Control_Group_Str,'_C',Row_Str),
          Probe_Type=='BS' & Infinium_Design==1 ~ paste0(Control_Group_Str,'_',Row_Str),
          Probe_Type=='BS' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
          Probe_Type=='NO' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
          Probe_Type=='ne' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
          TRUE ~ NA_character_ ),
        Seq_ID_Org=Seq_ID,
        Seq_ID=Control_Name
      ) %>%
      dplyr::select(Seq_ID,SR_Str,CO_Str,Infinium_Design,Rep_Num,Probe_Type,U,M, 
                    dplyr::everything())
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

format_controls_HSA = function(file1, file2,
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_controls_HSA'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    std_ctl_seq_tib <- dplyr::inner_join(
      suppressMessages(suppressWarnings(readr::read_tsv(file1) )) %>% 
        dplyr::mutate(Address=stringr::str_remove(address_name, '^1') %>% as.integer()),
      suppressMessages(suppressWarnings(
        readr::read_csv(file2, col_names=c("Address","Probe_Type","COLOR_CHANNEL","Probe_ID")) )),
      by="Address") %>%
      dplyr::select(Address,Probe_Type,COLOR_CHANNEL,Probe_ID,probe_id,sequence, everything()) %>%
      dplyr::rename(Design_ID=probe_id) %>% dplyr::distinct(Address, .keep_all=TRUE) %>%
      dplyr::select(-type_b,-bo_seq,-address_name)
    
    std_ctl_tib <- std_ctl_seq_tib %>% 
      dplyr::mutate(PIDX=Probe_ID %>%
                      stringr::str_replace('^.*[^0-9]([0-9]+)$', '\\$1') %>% 
                      stringr::str_remove_all('\\\\')) %>%
      dplyr::mutate(Probe_ID=Probe_ID %>%
                      stringr::str_replace_all(' ', '_') %>%
                      stringr::str_replace_all('-', '_') %>% 
                      stringr::str_replace_all('\\(', '') %>% 
                      stringr::str_replace_all('\\)', '')) %>%
      dplyr::rename(U=Address) %>%
      dplyr::mutate(M=NA_real_,M=as.double(M),U=as.double(U),
                    DESIGN='II',Infinium_Design=2,
                    col=NA_character_,Probe_Source='HSA',
                    Next_Base=NA_character_,
                    Probe_ID=paste('ctl',Probe_ID, sep='_')) %>%
      dplyr::select(Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,Probe_Type,Probe_Source,Next_Base,
                    dplyr::everything()) %>%
      dplyr::arrange(Probe_Type,Probe_ID) %>%
      dplyr::distinct(M,U, .keep_all=TRUE) %>% 
      dplyr::mutate(
        Design_Base_ID=stringr::str_remove(Design_ID, '_[AB]$'),
        Design_Base_AB=stringr::str_replace(Design_ID, '^.*_([AB])$','\\$1') %>%
          stringr::str_remove_all('\\\\'),
        Control_Group=Probe_Type,
        Control_Group_Str=stringr::str_replace_all(Control_Group,' ','_') %>% 
          stringr::str_replace_all('-','_'),
        Probe_Type=dplyr::case_when(
          Control_Group=="BISULFITE CONVERSION I"  ~ 'BS',
          Control_Group=="BISULFITE CONVERSION II" ~ 'BS',
          Control_Group=="EXTENSION"       ~ 'EX',
          Control_Group=="HYBRIDIZATION"   ~ 'HB',
          Control_Group=="NEGATIVE"        ~ 'NG',
          Control_Group=="NON-POLYMORPHIC" ~ 'NP',
          
          Control_Group=="NORM_A"          ~ 'MA',
          Control_Group=="NORM_C"          ~ 'MC',
          Control_Group=="NORM_G"          ~ 'MG',
          Control_Group=="NORM_T"          ~ 'MT',
          
          Control_Group=="RESTORATION"     ~ 'RE',
          Control_Group=="SPECIFICITY I"   ~ 'S1',
          Control_Group=="SPECIFICITY II"  ~ 'S2',
          
          Control_Group=="TARGET REMOVAL"  ~ 'TR',
          TRUE ~ NA_character_
        ),
        # AlleleB_Probe_Sequence=NA,
        DiNuc=NA_character_,
        N1=NA_character_,
        N2=NA_character_,
        Seq_ID=paste0(Control_Group_Str,'_',PIDX) %>%
          stringr::str_replace_all(' ', '_') %>%
          stringr::str_replace_all('-', '_') %>% 
          stringr::str_replace_all('\\(', '') %>% 
          stringr::str_replace_all('\\)', '')
      )
    
    std_bsU_tib <- std_ctl_tib %>% 
      dplyr::filter(stringr::str_detect(Probe_ID,'ctl_BS_Conversion_I_U')) %>%
      dplyr::rename(Address=U,AlleleA_Probe_Sequence=sequence) %>% 
      dplyr::mutate(Probe_ID=stringr::str_remove(Probe_ID,'_[A-Z][0-9]+$')) %>%
      dplyr::select(-M)
    
    std_bsC_tib <- std_ctl_tib %>% 
      dplyr::filter(stringr::str_detect(Probe_ID,'ctl_BS_Conversion_I_C')) %>%
      dplyr::rename(Address=U,AlleleB_Probe_Sequence=sequence) %>% 
      dplyr::mutate(Probe_ID=stringr::str_remove(Probe_ID,'_[A-Z][0-9]+$')) %>%
      dplyr::select(-M)
    
    # Unique Addresses to Remove from Standard Controls
    std_bs1_add_vec <- dplyr::bind_rows(std_bsU_tib,std_bsC_tib) %>% 
      dplyr::distinct(Address) %>% dplyr::pull(Address)
    
    std_bs1_hum_tib <- 
      dplyr::inner_join(std_bsU_tib,std_bsC_tib, 
                        by=c("Probe_ID","Design_Base_ID","Probe_Type","Probe_Source","PIDX","Control_Group","Control_Group_Str",
                             "DESIGN","col","Next_Base","DiNuc","N1","N2","Infinium_Design","Seq_ID"), suffix=c("_A","_B") ) %>% 
      dplyr::rename(U=Address_A,M=Address_B) %>%
      dplyr::mutate(
        DESIGN='I',
        Infinium_Design=1,
        COLOR_CHANNEL='Both',
        col=COLOR_CHANNEL_A,
        Next_Base='N',
        Design_Base_AB="AB",
        Design_ID=Design_Base_ID
      )
    
    #
    # Build remainder of Standard EPIC Controls::
    #
    std_inf1_hum_tib <- std_ctl_tib %>% 
      dplyr::filter(!U %in% std_bs1_add_vec) %>% 
      dplyr::rename(AlleleA_Probe_Sequence=sequence) %>%
      dplyr::mutate(col=COLOR_CHANNEL)
    
    #
    # TBD:: This does not contain all controls and loses pairs...
    #
    ret_tib <- 
      dplyr::bind_rows(std_inf1_hum_tib,std_bs1_hum_tib) %>%
      dplyr::distinct(M,U, .keep_all=TRUE) %>% 
      dplyr::filter(!is.na(Probe_ID)) %>%
      dplyr::select(-Design_Base_AB_A,-Design_Base_AB_B,
                    -Design_ID_A,-Design_ID_B) %>%
      dplyr::group_by(Probe_ID) %>%
      dplyr::mutate(PID_CNT=dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        # Without unique Replicate IDs:: i.e. PID_CNT
        # Probe_ID=paste0(stringr::str_replace_all(Probe_ID,'_','-'),'-HSA_NN',Infinium_Design,'1'),
        
        # With unique Replicate IDs:: i.e. PID_CNT
        Probe_ID=paste0(stringr::str_replace_all(Probe_ID,'_','-'),'-HSA_NN',Infinium_Design,PID_CNT),
        Assay_Class='Control',
        AQP=1,
        Probe_Class=Probe_Type
      )
    
    if (verbose>=vt+4) {
      ret_1_cnt <- 0
      ret_2_cnt <- 0
      
      ret_1_tib <- ret_tib %>% dplyr::filter(Infinium_Design==1)
      ret_2_tib <- ret_tib %>% dplyr::filter(Infinium_Design==2)
      
      ret_1_cnt <- ret_1_tib %>% base::nrow()
      ret_2_cnt <- ret_2_tib %>% base::nrow()
      
      cat(glue::glue("[{funcTag}]:{tabsStr} Infinium 1 Controls({ret_1_cnt})={RET}"))
      ret_1_tib %>% print(base::nrow(ret_1_tib))
      cat(glue::glue("[{funcTag}]:{tabsStr} Infinium 2 Controls({ret_2_cnt})={RET}"))
      ret_2_tib %>% print(base::nrow(ret_2_tib))
      cat(glue::glue("[{funcTag}]:{RET}"))
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Controls I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addControlTypeSlim = function(tib,
                              ret_str = "Control_Class",
                              vb=0, vt=3, tc=1, tt=NULL,
                              fun_tag='addControlTypeSlim')
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
  
  if ( p0 ) cat(glue::glue("{mssg} Starting...{RET2}"))
  if ( p1 ) {
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      fun_tag = '{fun_tag}'.{RET}"))
    cat(glue::glue("{mssg}      ret_str = '{ret_str}'.{RET}"))
    cat(glue::glue("{mssg}{RET}"))
    # cat(glue::glue("{mssg} Function Parameters::{RET}"))
    # cat(glue::glue("{mssg}     is_valid = '{is_valid}'.{RET}"))
    # cat(glue::glue("{mssg}       reload = '{reload}'.{RET}"))
    # cat(glue::glue("{mssg}   reload_min = '{reload_min}'.{RET}"))
    # cat(glue::glue("{mssg}   reload_pre = '{reload_pre}'.{RET}"))
    # cat(glue::glue("{mssg}     ret_data = '{ret_data}'.{RET}"))
    # cat(glue::glue("{mssg}     parallel = '{parallel}'.{RET}"))
    # cat(glue::glue("{mssg}{RET}"))
    # cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    # cat(glue::glue("{mssg}      run_tag = '{run_tag}'.{RET}"))
    # cat(glue::glue("{mssg}      out_tag = '{out_tag}'.{RET}"))
    # cat(glue::glue("{mssg}      out_dir = '{out_dir}'.{RET}"))
    # cat(glue::glue("{mssg}      beg_txt = '{beg_txt}'.{RET}"))
    # cat(glue::glue("{mssg}      sum_csv = '{sum_csv}'.{RET}"))
    # cat(glue::glue("{mssg}      aux_csv = '{aux_csv}'.{RET}"))
    # cat(glue::glue("{mssg}      out_csv = '{out_csv}'.{RET}"))
    # cat(glue::glue("{mssg}      end_txt = '{end_txt}'.{RET}"))
    # cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  ret_sym <- rlang::sym( ret_str )

  ret_tib <- tib %>% dplyr::mutate(
    # Control_Type = case_when(
    !!ret_sym := case_when(
      stringr::str_starts(Probe_ID, 'neg')  ~ 'NG',
      # 
      # stringr::str_starts(Probe_ID, 'neg')  ~ 'NEGATIVE',
      
      stringr::str_starts(Probe_ID, 'Norm') ~ 'NR',
      # stringr::str_starts(Probe_ID, 'Norm') ~ 'Normalization',
      # 
      
      stringr::str_starts(Probe_ID, 'BS_Conversion_II') ~ 'B2',
      stringr::str_starts(Probe_ID, 'BS_Conversion_I_') ~ 'B1',
      # 
      # stringr::str_starts(Probe_ID, 'BS_Conversion_II') ~ 'BsConversion_II',
      # stringr::str_starts(Probe_ID, 'BS_Conversion_I_') ~ 'BsConversion_I',
      
      stringr::str_starts(Probe_ID, 'BN') ~ 'BsConversion_NoBiasHIII_ASPE',
      # 
      # stringr::str_starts(Probe_ID, 'BS_conversion_NoBiasHIII_ASPE') ~ 'BsConversion_NoBiasHIII_ASPE',
      
      stringr::str_starts(Probe_ID, 'Non_Specific_II_') ~ 'N2',
      stringr::str_starts(Probe_ID, 'Non_Specific_I_')  ~ 'N1',
      # 
      # stringr::str_starts(Probe_ID, 'Non_Specific_II_') ~ 'NonSpecific_II',
      # stringr::str_starts(Probe_ID, 'Non_Specific_I_')  ~ 'NonSpecific_I',
      
      stringr::str_starts(Probe_ID, 'A_Hairpin1')   ~ 'HA',
      stringr::str_starts(Probe_ID, 'C_Hairpin1')   ~ 'HC',
      stringr::str_starts(Probe_ID, 'G_Hairpin1')   ~ 'HG',
      stringr::str_starts(Probe_ID, 'T_Hairpin1')   ~ 'HT',
      # 
      # stringr::str_starts(Probe_ID, 'A_Hairpin1')   ~ 'Hairpin_A',
      # stringr::str_starts(Probe_ID, 'C_Hairpin1')   ~ 'Hairpin_C',
      # stringr::str_starts(Probe_ID, 'G_Hairpin1')   ~ 'Hairpin_G',
      # stringr::str_starts(Probe_ID, 'T_Hairpin1')   ~ 'Hairpin_T',
      # 
      
      stringr::str_starts(Probe_ID, 'A_Hairpin2')   ~ 'GA',
      stringr::str_starts(Probe_ID, 'C_Hairpin2')   ~ 'GC',
      stringr::str_starts(Probe_ID, 'G_Hairpin2')   ~ 'GG',
      stringr::str_starts(Probe_ID, 'T_Hairpin2')   ~ 'GT',
      # 
      # stringr::str_starts(Probe_ID, 'A_Hairpin2')   ~ 'Hairpin_A',
      # stringr::str_starts(Probe_ID, 'C_Hairpin2')   ~ 'Hairpin_C',
      # stringr::str_starts(Probe_ID, 'G_Hairpin2')   ~ 'Hairpin_G',
      # stringr::str_starts(Probe_ID, 'T_Hairpin2')   ~ 'Hairpin_T',
      
      stringr::str_starts(Probe_ID, 'GT_mismatch')  ~ 'MT',
      # 
      # stringr::str_starts(Probe_ID, 'GT_mismatch')  ~ 'Mismatch_GT',
      
      stringr::str_starts(Probe_ID, '2_HIGH_MM')      ~ 'MH',
      stringr::str_starts(Probe_ID, '3_HIGH_MM')      ~ 'MH',
      stringr::str_starts(Probe_ID, '18_MEDIUM_MM')   ~ 'MM',
      stringr::str_starts(Probe_ID, '39_MEDIUM_MM')   ~ 'MM',
      stringr::str_starts(Probe_ID, '60_LOW_MM')      ~ 'ML',
      stringr::str_starts(Probe_ID, '74_LOW_MM')      ~ 'ML',
      # 
      # stringr::str_starts(Probe_ID, '2_HIGH_MM')      ~ 'HIGH_MM_2',
      # stringr::str_starts(Probe_ID, '3_HIGH_MM')      ~ 'HIGH_MM_3',
      # stringr::str_starts(Probe_ID, '18_MEDIUM_MM')   ~ 'MEDIUM_MM_18',
      # stringr::str_starts(Probe_ID, '39_MEDIUM_MM')   ~ 'MEDIUM_MM_39',
      # stringr::str_starts(Probe_ID, '60_LOW_MM')      ~ 'LOW_MM_68',
      # stringr::str_starts(Probe_ID, '74_LOW_MM')      ~ 'LOW_MM_74',
      
      stringr::str_starts(Probe_ID, '74_YEAST_3MM')   ~ 'YS',
      stringr::str_starts(Probe_ID, '90_YEAST_3MM')   ~ 'YS',
      # 
      # stringr::str_starts(Probe_ID, '74_YEAST_3MM')   ~ 'YEAST_3MM_74',
      # stringr::str_starts(Probe_ID, '90_YEAST_3MM')   ~ 'YEAST_3MM_90',
      
      stringr::str_starts(Probe_ID, 'nonPolyA_ATG2')  ~ 'NA',
      stringr::str_starts(Probe_ID, 'nonPolyA_HK1')   ~ 'NA',
      stringr::str_starts(Probe_ID, 'nonPolyA_PPIH')  ~ 'NA',
      # 
      # stringr::str_starts(Probe_ID, 'nonPolyA_ATG2')  ~ 'NonPoly_A',
      # stringr::str_starts(Probe_ID, 'nonPolyA_HK1')   ~ 'NonPoly_A',
      # stringr::str_starts(Probe_ID, 'nonPolyA_PPIH')  ~ 'NonPoly_A',
      
      stringr::str_starts(Probe_ID, 'nonPolyC_HK2')   ~ 'NC',
      stringr::str_starts(Probe_ID, 'nonPolyC_PPID')  ~ 'NC',
      stringr::str_starts(Probe_ID, 'nonPolyC_PPIE')  ~ 'NC',
      #
      # NOT SURE WHY THIS IS NE and NC::
      # stringr::str_starts(Probe_ID, 'nonPolyC_PPIE')  ~ 'NonPoly_CE',
      # 
      # stringr::str_starts(Probe_ID, 'nonPolyC_HK2')   ~ 'NonPoly_C',
      # stringr::str_starts(Probe_ID, 'nonPolyC_PPID')  ~ 'NonPoly_C',
      # stringr::str_starts(Probe_ID, 'nonPolyC_PPIE')  ~ 'NonPoly_CE',
      
      stringr::str_starts(Probe_ID, 'nonPolyG_HK2')   ~ 'NG',
      stringr::str_starts(Probe_ID, 'nonPolyG_HK3')   ~ 'NG',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIE')  ~ 'NG',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIH')  ~ 'NG',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIG')  ~ 'NG',
      # 
      # stringr::str_starts(Probe_ID, 'nonPolyG_HK2')   ~ 'NonPoly_G',
      # stringr::str_starts(Probe_ID, 'nonPolyG_HK3')   ~ 'NonPoly_G',
      # stringr::str_starts(Probe_ID, 'nonPolyG_PPIE')  ~ 'NonPoly_G',
      # stringr::str_starts(Probe_ID, 'nonPolyG_PPIH')  ~ 'NonPoly_G',
      # stringr::str_starts(Probe_ID, 'nonPolyG_PPIG')  ~ 'NonPoly_G',
      
      stringr::str_starts(Probe_ID, 'nonPolyT_ALDOB') ~ 'NT',
      stringr::str_starts(Probe_ID, 'nonPolyT_HK2')   ~ 'NT',
      stringr::str_starts(Probe_ID, 'nonPolyT_PPIH')  ~ 'NT',
      # 
      # stringr::str_starts(Probe_ID, 'nonPolyT_ALDOB') ~ 'NonPoly_T',
      # stringr::str_starts(Probe_ID, 'nonPolyT_HK2')   ~ 'NonPoly_T',
      # stringr::str_starts(Probe_ID, 'nonPolyT_PPIH')  ~ 'NonPoly_T',
      
      stringr::str_starts(Probe_ID, 'CA8')  ~ 'CA',
      stringr::str_starts(Probe_ID, 'CA7')  ~ 'CA',
      stringr::str_starts(Probe_ID, 'CA6')  ~ 'CA',
      stringr::str_starts(Probe_ID, 'CA5')  ~ 'CA',
      stringr::str_starts(Probe_ID, 'CA1')  ~ 'CA',
      # 
      # stringr::str_starts(Probe_ID, 'CA8')  ~ 'CA8',
      # stringr::str_starts(Probe_ID, 'CA7')  ~ 'CA7',
      # stringr::str_starts(Probe_ID, 'CA6')  ~ 'CA6',
      # stringr::str_starts(Probe_ID, 'CA5')  ~ 'CA5',
      # stringr::str_starts(Probe_ID, 'CA1')  ~ 'CA1',
      
      stringr::str_starts(Probe_ID, 'CC8')  ~ 'CC',
      stringr::str_starts(Probe_ID, 'CC7')  ~ 'CC',
      stringr::str_starts(Probe_ID, 'CC5')  ~ 'CC',
      stringr::str_starts(Probe_ID, 'CC4')  ~ 'CC',
      stringr::str_starts(Probe_ID, 'CC3')  ~ 'CC',
      stringr::str_starts(Probe_ID, 'CC1')  ~ 'CC',
      # 
      # stringr::str_starts(Probe_ID, 'CC8')  ~ 'CC8',
      # stringr::str_starts(Probe_ID, 'CC7')  ~ 'CC7',
      # stringr::str_starts(Probe_ID, 'CC5')  ~ 'CC5',
      # stringr::str_starts(Probe_ID, 'CC4')  ~ 'CC4',
      # stringr::str_starts(Probe_ID, 'CC3')  ~ 'CC3',
      # stringr::str_starts(Probe_ID, 'CC1')  ~ 'CC1',
      
      stringr::str_starts(Probe_ID, 'CG12') ~ 'CG',
      stringr::str_starts(Probe_ID, 'CG10') ~ 'CG',
      stringr::str_starts(Probe_ID, 'CG8')  ~ 'CG',
      stringr::str_starts(Probe_ID, 'CG2')  ~ 'CG',
      stringr::str_starts(Probe_ID, 'CG1')  ~ 'CG',
      # 
      # stringr::str_starts(Probe_ID, 'CG12') ~ 'CG12',
      # stringr::str_starts(Probe_ID, 'CG10') ~ 'CG10',
      # stringr::str_starts(Probe_ID, 'CG8')  ~ 'CG8',
      # stringr::str_starts(Probe_ID, 'CG2')  ~ 'CG2',
      # stringr::str_starts(Probe_ID, 'CG1')  ~ 'CG1',
      
      stringr::str_starts(Probe_ID, 'CT12') ~ 'CT',
      stringr::str_starts(Probe_ID, 'CT10') ~ 'CT',
      stringr::str_starts(Probe_ID, 'CT9')  ~ 'CT',
      stringr::str_starts(Probe_ID, 'CT5')  ~ 'CT',
      stringr::str_starts(Probe_ID, 'CT2')  ~ 'CT',
      # 
      # stringr::str_starts(Probe_ID, 'CT12') ~ 'CT12',
      # stringr::str_starts(Probe_ID, 'CT10') ~ 'CT10',
      # stringr::str_starts(Probe_ID, 'CT9')  ~ 'CT9',
      # stringr::str_starts(Probe_ID, 'CT5')  ~ 'CT5',
      # stringr::str_starts(Probe_ID, 'CT2')  ~ 'CT2',
      
      stringr::str_starts(Probe_ID, 'TM183T')  ~ 'TM',
      stringr::str_starts(Probe_ID, 'TM182T')  ~ 'TM',
      stringr::str_starts(Probe_ID, 'TM169T')  ~ 'TM',
      stringr::str_starts(Probe_ID, 'TM167T')  ~ 'TM',
      stringr::str_starts(Probe_ID, 'TM150T')  ~ 'TM',
      stringr::str_starts(Probe_ID, 'TM148T')  ~ 'TM',
      # 
      # stringr::str_starts(Probe_ID, 'TM183T')  ~ 'TM183T',
      # stringr::str_starts(Probe_ID, 'TM182T')  ~ 'TM182T',
      # stringr::str_starts(Probe_ID, 'TM169T')  ~ 'TM169T',
      # stringr::str_starts(Probe_ID, 'TM167T')  ~ 'TM167T',
      # stringr::str_starts(Probe_ID, 'TM150T')  ~ 'TM150T',
      # stringr::str_starts(Probe_ID, 'TM148T')  ~ 'TM148T',
      
      # stringr::str_starts(Probe_ID, 'rs') ~ 'SNP',
      # 
      stringr::str_starts(Probe_ID, 'rs') ~ 'RS',
      
      TRUE ~ NA_character_)
  )
  
  if ( p1 ) {
    mat_cnt <- 0
    mat_cnt <- ret_tib %>% 
      dplyr::group_by(Control_Type) %>% 
      dplyr::summarise(Probe_Count=n()) %>% 
      base::nrow()
    
    unk_cnt <- 0
    unk_cnt <- ret_tib %>% 
      dplyr::filter(is.na(Control_Type) ) %>% 
      dplyr::select(Probe_ID) %>% base::nrow()
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Match Count={mat_cnt}, Unknown Type={unk_cnt}.{RET}"))
    sum_tib <- ret_tib %>% 
      dplyr::group_by(Control_Type) %>% 
      dplyr::summarise( Probe_Count=n(), .groups = "drop" )
    sum_tib %>% print( n=base::nrow( n=base::nrow(sum_tib) ) )
  }
  etime <- round( as.double(ftime[3]), 2 )
  if ( !is.null(tt) ) tt$addTime( ftime,fun_tag )
  if ( p0 ) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))

  ret_tib
}

addControlType = function(tib,
                          verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'addControlType'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- tib %>% dplyr::mutate(
    Control_Type=case_when(
      stringr::str_starts(Probe_ID, 'neg')  ~ 'NEGATIVE',
      stringr::str_starts(Probe_ID, 'Norm') ~ 'Normalization',
      
      stringr::str_starts(Probe_ID, 'BS_Conversion_II') ~ 'BsConversion_II',
      stringr::str_starts(Probe_ID, 'BS_Conversion_I_') ~ 'BsConversion_I',
      
      stringr::str_starts(Probe_ID, 'BS_conversion_NoBiasHIII_ASPE') ~ 'BsConversion_NoBiasHIII_ASPE',
      
      stringr::str_starts(Probe_ID, 'Non_Specific_II_') ~ 'NonSpecific_II',
      stringr::str_starts(Probe_ID, 'Non_Specific_I_')  ~ 'NonSpecific_I',
      
      stringr::str_starts(Probe_ID, 'A_Hairpin1')   ~ 'Hairpin_A_1',
      stringr::str_starts(Probe_ID, 'C_Hairpin1')   ~ 'Hairpin_C_1',
      stringr::str_starts(Probe_ID, 'G_Hairpin1')   ~ 'Hairpin_G_1',
      stringr::str_starts(Probe_ID, 'T_Hairpin1')   ~ 'Hairpin_T_1',
      
      stringr::str_starts(Probe_ID, 'A_Hairpin2')   ~ 'Hairpin_A_2',
      stringr::str_starts(Probe_ID, 'C_Hairpin2')   ~ 'Hairpin_C_2',
      stringr::str_starts(Probe_ID, 'G_Hairpin2')   ~ 'Hairpin_G_2',
      stringr::str_starts(Probe_ID, 'T_Hairpin2')   ~ 'Hairpin_T_2',
      
      stringr::str_starts(Probe_ID, 'GT_mismatch')  ~ 'Mismatch_GT',
      
      stringr::str_starts(Probe_ID, '2_HIGH_MM')      ~ 'HIGH_MM_2',
      stringr::str_starts(Probe_ID, '3_HIGH_MM')      ~ 'HIGH_MM_3',
      stringr::str_starts(Probe_ID, '18_MEDIUM_MM')   ~ 'MEDIUM_MM_18',
      stringr::str_starts(Probe_ID, '39_MEDIUM_MM')   ~ 'MEDIUM_MM_39',
      stringr::str_starts(Probe_ID, '60_LOW_MM')      ~ 'LOW_MM_68',
      stringr::str_starts(Probe_ID, '74_LOW_MM')      ~ 'LOW_MM_74',
      
      stringr::str_starts(Probe_ID, '74_YEAST_3MM')   ~ 'YEAST_3MM_74',
      stringr::str_starts(Probe_ID, '90_YEAST_3MM')   ~ 'YEAST_3MM_90',
      
      stringr::str_starts(Probe_ID, 'nonPolyA_ATG2')  ~ 'NonPoly_A_ATG2',
      stringr::str_starts(Probe_ID, 'nonPolyA_HK1')   ~ 'NonPoly_A_HK1',
      stringr::str_starts(Probe_ID, 'nonPolyA_PPIH')  ~ 'NonPoly_A_PPIH',
      
      stringr::str_starts(Probe_ID, 'nonPolyC_HK2')   ~ 'NonPoly_C_HK2',
      stringr::str_starts(Probe_ID, 'nonPolyC_PPID')  ~ 'NonPoly_C_PPID',
      stringr::str_starts(Probe_ID, 'nonPolyC_PPIE')  ~ 'NonPoly_C_PPIE',
      
      stringr::str_starts(Probe_ID, 'nonPolyG_HK2')   ~ 'NonPoly_G_HK2',
      stringr::str_starts(Probe_ID, 'nonPolyG_HK3')   ~ 'NonPoly_G_HK3',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIE')  ~ 'NonPoly_G_PPIE',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIH')  ~ 'NonPoly_G_PPIH',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIG')  ~ 'NonPoly_G_PPIG',
      
      stringr::str_starts(Probe_ID, 'nonPolyT_ALDOB') ~ 'NonPoly_T_ALDOB',
      stringr::str_starts(Probe_ID, 'nonPolyT_HK2')   ~ 'NonPoly_T_HK2',
      stringr::str_starts(Probe_ID, 'nonPolyT_PPIH')  ~ 'NonPoly_T_PPIH',
      
      stringr::str_starts(Probe_ID, 'CA8')  ~ 'CA8',
      stringr::str_starts(Probe_ID, 'CA7')  ~ 'CA7',
      stringr::str_starts(Probe_ID, 'CA6')  ~ 'CA6',
      stringr::str_starts(Probe_ID, 'CA5')  ~ 'CA5',
      stringr::str_starts(Probe_ID, 'CA1')  ~ 'CA1',
      
      stringr::str_starts(Probe_ID, 'CC8')  ~ 'CC8',
      stringr::str_starts(Probe_ID, 'CC7')  ~ 'CC7',
      stringr::str_starts(Probe_ID, 'CC5')  ~ 'CC5',
      stringr::str_starts(Probe_ID, 'CC4')  ~ 'CC4',
      stringr::str_starts(Probe_ID, 'CC3')  ~ 'CC3',
      stringr::str_starts(Probe_ID, 'CC1')  ~ 'CC1',
      
      stringr::str_starts(Probe_ID, 'CG12') ~ 'CG12',
      stringr::str_starts(Probe_ID, 'CG10') ~ 'CG10',
      stringr::str_starts(Probe_ID, 'CG8')  ~ 'CG8',
      stringr::str_starts(Probe_ID, 'CG2')  ~ 'CG2',
      stringr::str_starts(Probe_ID, 'CG1')  ~ 'CG1',
      
      stringr::str_starts(Probe_ID, 'CT12') ~ 'CT12',
      stringr::str_starts(Probe_ID, 'CT10') ~ 'CT10',
      stringr::str_starts(Probe_ID, 'CT9')  ~ 'CT9',
      stringr::str_starts(Probe_ID, 'CT5')  ~ 'CT5',
      stringr::str_starts(Probe_ID, 'CT2')  ~ 'CT2',
      
      stringr::str_starts(Probe_ID, 'TM183T')  ~ 'TM183T',
      stringr::str_starts(Probe_ID, 'TM182T')  ~ 'TM182T',
      stringr::str_starts(Probe_ID, 'TM169T')  ~ 'TM169T',
      stringr::str_starts(Probe_ID, 'TM167T')  ~ 'TM167T',
      stringr::str_starts(Probe_ID, 'TM150T')  ~ 'TM150T',
      stringr::str_starts(Probe_ID, 'TM148T')  ~ 'TM148T',
      
      stringr::str_starts(Probe_ID, 'rs') ~ 'SNP',
      
      TRUE ~ NA_character_)
  )
  
  if (verbose>=vt) {
    mat_cnt <- tib %>% dplyr::group_by(Control_Type) %>% dplyr::summarise(Probe_Count=n()) %>% base::nrow()
    unk_cnt <- tib %>% dplyr::filter(is.na(Control_Type)) %>% dplyr::select(Probe_ID) %>% base::nrow()
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Match Count={mat_cnt}, Unknown Type={unk_cnt}.{RET}"))
    tib %>% dplyr::group_by(Control_Type) %>% dplyr::summarise(Probe_Count=n()) %>% print()
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}"))
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       BSP Map Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

combineBspBed = function(bsp, bed, src,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'combineBspBed'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    bsp_cols <- c('Probe_ID', 'BSP_PrbSeq', 'BSP_Qual', 'BSP_Tag', 'BSP_Chrom', 'BSP_Pos', 'BSP_SRD', 'BSP_GapCnt', 'BSP_RefSeq', 'BSP_TopMisCnt', 'BSP_AllMisStr')
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading bsp={bsp}...{RET}"))
    bsp_tib <- suppressMessages(suppressWarnings( readr::read_tsv(bsp, col_names=bsp_cols) )) %>% dplyr::select(-BSP_Qual)
    #  dplyr::mutate(Design_Method="BSC") %>% 
    bsp_cnt <- bsp_tib %>% base::nrow()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} bsp_cnt={bsp_cnt}...{RET}{RET}"))
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading bed={bed}...{RET}"))
    bed_tib <- suppressMessages(suppressWarnings( readr::read_tsv(bed) ))
    
    tib <- bsp_tib %>% dplyr::filter(BSP_Tag=='UM') %>% dplyr::inner_join(bed_tib, by=c("Probe_ID"="Full_ID") ) %>% 
      dplyr::mutate(Design_Method=src,
                    Probe_New_Seq_A=paste0(PRB1_U_SS,'G'),
                    Probe_New_Seq_B=paste0(PRB1_M_SS,'A'),
                    Probe_New_Seq_D=PRB2_D_UC)
    bed_cnt <- bed_tib %>% base::nrow()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} bed_cnt={bed_cnt}...{RET}{RET}"))
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  tib
}

bsp2prbsControls = function(tib, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'bsp2prbsControls'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt) tib %>% dplyr::group_by(Control_Type) %>% dplyr::summarise(Count=n()) %>% as.data.frame() %>% print()
  
  prbs <- tib  %>%
    dplyr::mutate(
      BSP_RefBscU=case_when(
        # Control_Type=='BsConversion_II' & BSP_SRD=='--' ~ BSP_RefSeq,
        # Control_Type=='BsConversion_II' & BSP_SRD=='-+' ~ revCmp(BSP_RefSeq),
        # Control_Type=='BsConversion_II' & BSP_SRD=='+-' ~ revCmp(BSP_RefSeq),
        # Control_Type=='BsConversion_II' & BSP_SRD=='++' ~ BSP_RefSeq,
        
        stringr::str_starts(Control_Type, 'NonPoly_G') & BSP_SRD=='--' ~ BSP_RefSeq,          # NOT CHECKED
        stringr::str_starts(Control_Type, 'NonPoly_G') & BSP_SRD=='-+' ~ bscUs(revCmp(BSP_RefSeq)),
        stringr::str_starts(Control_Type, 'NonPoly_G') & BSP_SRD=='+-' ~ revCmp(BSP_RefSeq),  # NOT CHECKED
        stringr::str_starts(Control_Type, 'NonPoly_G') & BSP_SRD=='++' ~ bscUs(BSP_RefSeq),
        
        # TBD:: NORM_A/G Seems to be working now::
        stringr::str_starts(Control_Type, 'NORM_G') & BSP_SRD=='--' ~ BSP_RefSeq,          # Can Probably Remove
        stringr::str_starts(Control_Type, 'NORM_A') & BSP_SRD=='-+' ~ bscUs(revCmp(BSP_RefSeq)),
        stringr::str_starts(Control_Type, 'NORM_G') & BSP_SRD=='-+' ~ bscUs(revCmp(BSP_RefSeq)),
        stringr::str_starts(Control_Type, 'NORM_G') & BSP_SRD=='+-' ~ revCmp(BSP_RefSeq),  # Can Probably Remove
        stringr::str_starts(Control_Type, 'NORM_A') & BSP_SRD=='++' ~ bscUs(BSP_RefSeq),
        stringr::str_starts(Control_Type, 'NORM_G') & BSP_SRD=='++' ~ bscUs(BSP_RefSeq),
        
        BSP_SRD=='--' ~ revCmp(bscUs(revCmp(BSP_RefSeq))), # 728
        BSP_SRD=='-+' ~ revCmp(BSP_RefSeq),
        BSP_SRD=='+-' ~ revCmp(bscUs(BSP_RefSeq)), # 993
        BSP_SRD=='++' ~ BSP_RefSeq,  # Working So Far
        TRUE ~ NA_character_
      ),
      
      BSP_PrePosB= 2,
      BSP_NxbPosB= 3+stringr::str_length(BSP_PrbSeq),
      BSP_NxbPosE= 3+stringr::str_length(BSP_PrbSeq)+1,
      BSP_RefPreU=BSP_RefBscU %>% stringr::str_to_upper() %>% stringr::str_sub(BSP_PrePosB,BSP_PrePosB),
      BSP_RefNxbU=BSP_RefBscU %>% stringr::str_to_upper() %>% stringr::str_sub(BSP_NxbPosB,BSP_NxbPosB),
      BSP_RefPrbU=BSP_RefBscU %>% stringr::str_to_upper() %>% stringr::str_sub(3) %>% stringr::str_sub(1,stringr::str_length(BSP_PrbSeq)),
      BSP_CmpPrbU=BSP_RefPrbU %>% stringr::str_to_upper() %>% stringr::str_sub(1,stringr::str_length(BSP_RefPrbU)-1)
    )
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}"))
  
  prbs
}

loadControlProbes = function(name, dir, opt,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadControlProbes'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  # name <- '01152015_DarkMatterControls.unique.probe.match'
  ctl_tsv <- file.path(dir, name)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading ctl_tsv={ctl_tsv}...{RET}"))
  ctl_tib <- suppressMessages(suppressWarnings( readr::read_tsv(ctl_tsv) )) %>% dplyr::distinct()
  
  ctl_cnt <- ctl_tib %>% base::nrow()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loaded Control Sequences(n={ctl_cnt}).{RET}"))
  
  seq_tib <- ctl_tib %>% dplyr::rename(Probe_ID=probe_id, Probe_Seq=sequence, Address=address_name) %>% 
    dplyr::mutate(Control_ID=stringr::str_remove(Probe_ID, '_[AB]$'),
                  AB_Str=stringr::str_replace(Probe_ID, '^.*_([AB])$', '\\$1') %>% stringr::str_remove('\\\\'),
                  Last_Base_Src=stringr::str_sub(Probe_Seq, -1),
                  Address=stringr::str_remove(Address, '^1') %>% stringr::str_remove('^0+') %>% as.numeric(),
                  Address_Src=Address,
                  Probe_Seq_S=Probe_Seq %>% stringr::str_to_upper() %>% stringr::str_sub(1,stringr::str_length(Probe_Seq)-1)
    ) %>% addControlType() %>% 
    # dplyr::select(Control_Type, Control_ID, Probe_ID, AB_Str, Last_Base_Src, Address, Address_Src, Probe_Seq, Probe_Seq_S) %>%
    dplyr::arrange(Control_Type, Control_ID)
  
  # Format from stack to paired probes
  #
  prb_tib <- dplyr::full_join(seq_tib %>% dplyr::filter(AB_Str=='A'),
                              seq_tib %>% dplyr::filter(AB_Str=='B'),
                              by=c("Control_ID", "Control_Type"), suffix=c('_A', '_B')) %>%
    dplyr::mutate(Design_Type=case_when(!is.na(AB_Str_B) ~ 'I',
                                        is.na(AB_Str_B) ~ 'II',
                                        TRUE ~ NA_character_),
                  Last_Base_A=stringr::str_sub(Probe_Seq_A, -1),
                  Last_Base_B=stringr::str_sub(Probe_Seq_B, -1) )
  
  # Sanity Checks::
  #  Check that all Infinium I probes are equal except the last base::
  #  Check the purposely swapped base distributions::
  #
  if (verbose>=vt) {
    prb_tib %>% dplyr::filter(Design_Type=='I') %>% dplyr::filter(Probe_Seq_S_A!=Probe_Seq_S_B) %>% print()
    prb_tib %>%
      dplyr::filter(Design_Type=='I') %>%
      dplyr::group_by(AB_Str_A, Last_Base_Src_A, AB_Str_B, Last_Base_Src_B) %>% dplyr::summarise(Count=n()) %>% print()
    prb_tib %>% dplyr::filter(Design_Type=='I') %>%
      dplyr::group_by(AB_Str_A, Last_Base_A, AB_Str_B, Last_Base_B) %>% dplyr::summarise(Count=n()) %>% print()
  }
  
  # Last Base Mapping::
  last_base_map <- prb_tib %>% dplyr::group_by(Control_Type, AB_Str_A, AB_Str_B, Last_Base_A, Last_Base_B) %>% 
    dplyr::summarise(Count=n()) %>% arrange(AB_Str_A, AB_Str_B)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Write Original Manifest Fasta::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  opt$write_original_fasta <- FALSE
  if (opt$write_original_fasta) {
    ord_fas <- file.path(opt$outDir, paste(name,'bscU.fa.gz', sep='.') )
    fas_tib <- prb_tib %>% 
      dplyr::select(Probe_ID_A, Address_A, Probe_Seq_A) %>% 
      dplyr::mutate(line=paste0('>',Address_A,'_',Probe_ID_A,'\n',Probe_Seq_A)) %>% dplyr::pull(line)
    if (opt$writeOutput) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing ord_fas={ord_fas}...{RET}"))
      readr::write_lines(x=fas_tib, path=ord_fas)
    }
    
    # Print all Data::
    ord_fas <- file.path(opt$outDir, paste(name,'fa.gz', sep='.') )
    fas_tib <- c(prb_tib %>% dplyr::select(Probe_ID_A, Address_A, Probe_Seq_A) %>% dplyr::filter(!is.na(Probe_ID_A)) %>%
                   dplyr::mutate(line=paste0('>',Address_A,'_',Probe_ID_A,'\n',Probe_Seq_A)) %>% dplyr::pull(line),
                 prb_tib %>% dplyr::select(Probe_ID_B, Address_B, Probe_Seq_B) %>% dplyr::filter(!is.na(Probe_ID_B)) %>%
                   dplyr::mutate(line=paste0('>',Address_B,'_',Probe_ID_B,'\n',Probe_Seq_B)) %>% dplyr::pull(line) )
    if (opt$writeOutput) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing ord_fas={ord_fas}...{RET}"))
      readr::write_lines(x=fas_tib, path=ord_fas)
    }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}"))
  
  prb_tib
}

formatBSP = function(bsp) {
  bsp_bed_cols <- c('BSP_Chrom','BSP_Beg','BSP_End','Probe_ID','BSP_Tag','BSP_Srd','BSP_PrbSeq','BSP_RefSeq','BSP_GapCnt','BSP_MisCnt','BSP_MisStr',
                    'GEN_Chrom','GEN_Beg','GEN_End','GEN_Tran','GEN_Srd','GEN_Class')
  
  tib <- readr::read_tsv(bsp, col_names=bsp_bed_cols) %>% 
    tidyr::separate(BSP_MisStr, into=c('H0','H1','H2','H3','H4','H5'), sep=':', convert=TRUE) %>%
    dplyr::filter(H0>0) %>%
    tidyr::separate(Probe_ID, into=c('CGN', 'FR', 'TB', 'CO', 'ScrU', 'ScrM', 'ImpSumCnt', 'ImpCgnCnt'), sep='_') %>% 
    tidyr::unite(Probe_ID, CGN, FR, TB, CO, remove=FALSE, sep='_') %>%
    tidyr::separate(GEN_Tran, into=c('Gene', 'Transcript'), sep=':') %>%
    dplyr::mutate(MinScore=as.double(pmin(ScrU,ScrM)))
  
  tib
}


# End of file

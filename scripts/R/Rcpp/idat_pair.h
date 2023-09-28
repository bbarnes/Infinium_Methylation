#ifndef __IDAT_PAIR_H__
#define __IDAT_PAIR_H__

#include <stdio.h>

#include "idat.h"

using namespace idat;

namespace idat_pair {

class Idat_Pair
{
private:
  
  // Workflow Fields::
  // std::vector<_S> workflows;
  // _S work_str;
  // _S step_str;
  
  // Class Options Map::
  _S output_path;
  std::map< _S,_d > opts_dbl;
  std::map< _S,_S > opts_str;
  
  // Top Matched Manifest::
  std::pair<_S,_d> man_top_pair;
  
  // Maps:: idx_to_pair <=> pair_to_idx::
  std::vector< std::pair<_I,_I> > idx_to_pair;
  std::map< std::pair<_I,_I>,_I > pair_to_idx;
  
  // Keep track of dins(Probe_Type) to idx::
  std::map< _S, std::vector<_I> > dins_to_idx;
  
  // Data:: idx_to_data
  std::map< std::pair<_S,_S>, std::vector<_d> > pair_idx_to_data;
  // std::vector< _d > idx_to_beta;
  // std::vector< _d > idx_to_pval;
  
  // Stats:: key_to_stat
  std::map< std::pair<_S,_S>, std::map< _S,_d > > dins_to_sum;
  
  // Info:: key_to_info (strings/boolean)
  std::map< _S, std::vector< _S > > key_to_info;
  std::map< _S, std::vector< _Y > > key_to_bool;

  // Sentrix:: key => info/data (string/double)
  std::map< _S,_S > sentrix_info;
  std::map< _S,_d > sentrix_data;
  
public: 
  
  Idat_Pair( const _S prefix_path, const _S& output_path,
             std::vector<_I>& pval_add_vec,
             std::vector<_I>& addU_man_vec, std::vector<_I>& addM_man_vec,
             std::vector<_S>& cgns_man_vec, std::vector<_S>& cols_man_vec,
             std::vector<_S>& name_man_vec, std::vector<_S>& anno_man_vec,
             std::vector<_S>& chrs_man_vec,
             const _d min_pval = 0.05,
             const _d min_beta = 0.30, const _d max_beta = 0.70,
             const _d min_perI = 0.05, const _d min_perO = 0.75,
             const _Y read_bgz  = false, 
             const _Y write_bgz = false,
             const _Y rm_pval_outliers = false,
             const _I vb = 0, const _I vt = 1, const _I tc = 0,
             const _S ft = "Idat_Pair" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    // const 
    // this->opts[]
    
    const _Y p0 = vb > vt + 0;
    _Y success = true;
    
    _I pval_size = pval_add_vec.size();
    _I cgns_size = cgns_man_vec.size();
    
    if ( cgns_size == 0 ) success = false;
    if ( cgns_size != addU_man_vec.size() ) success = false;
    if ( cgns_size != addM_man_vec.size() ) success = false;
    if ( cgns_size != name_man_vec.size() ) success = false;
    if ( cgns_size != cols_man_vec.size() ) success = false;
    if ( cgns_size != chrs_man_vec.size() ) success = false;
    if ( cgns_size != anno_man_vec.size() ) success = false;
    if ( !success ) std::cerr
      <<_fe<<"Manifest Input Vector Size's do NOT Match! Exiting..,\n "
      <<_fe<<"\t cgns = '"<<cgns_size<<"'\n"
      <<_fe<<"\t addU = '"<<addU_man_vec.size()<<"'\n"
      <<_fe<<"\t addM = '"<<addM_man_vec.size()<<"'\n"
      <<_fe<<"\t name = '"<<name_man_vec.size()<<"'\n"
      <<_fe<<"\t cols = '"<<cols_man_vec.size()<<"'\n"
      <<_fe<<"\t chrs = '"<<chrs_man_vec.size()<<"'\n"
      <<_fe<<"\t anno = '"<<anno_man_vec.size()<<"'\n"
      <<std::endl;
    
    if ( success ) {
      
      if ( p0 ) std::cerr
        <<_fo<<"Statring...\n"
        <<_fo<<"\t       (vb,vt,tc): ("<<vb<<","<<vt<<","<<tc<<")\n"
        <<_fo<<"\t      prefix_path: '"<<prefix_path<<"'\n"
        <<_fo<<"\t      output_path: '"<<output_path<<"'\n"
        <<_fo<<"\t        pval_size: '"<<pval_size<<"'\n"
        <<_fo<<"\t        cgns_size: '"<<cgns_size<<"'\n\n"
        <<_fo<<"\t        pval_size: '"<<pval_size<<"'\n"
        <<_fo<<"\t         min_pval: '"<<min_pval<<"'\n"
        <<_fo<<"\t         min_beta: '"<<min_beta<<"'\n"
        <<_fo<<"\t         max_beta: '"<<max_beta<<"'\n"
        <<_fo<<"\t         min_per0: '"<<min_perO<<"'\n"
        <<_fo<<"\t         min_perI: '"<<min_perI<<"'\n\n"
        <<_fo<<"\t         read_bgz: '"<<read_bgz<<"'\n"
        <<_fo<<"\t        write_bgz: '"<<write_bgz<<"'\n"
        <<_fo<<"\t rm_pval_outliers: '"<<rm_pval_outliers<<"'\n"
        <<std::endl;
      
      if ( !read_bgz ) {
        success = this->parse_idat_pair( prefix_path, output_path, 
                                         pval_add_vec, 
                                         addU_man_vec, addM_man_vec,
                                         cgns_man_vec, cols_man_vec,
                                         name_man_vec, anno_man_vec,
                                         chrs_man_vec, 
                                         min_pval, min_beta, max_beta,
                                         min_perO, min_perI, 
                                         read_bgz, write_bgz,
                                         rm_pval_outliers,
                                         vb,vt+1,tc );
      } else {
        /*
         * BGZ IO Support is Priority. Removing this for now...
         *   Example code is in the Rcpp/core-bk directory!
         */
        std::cerr <<_fe<<"BGZ Input is NOT currently supported: '" 
                  <<prefix_path<<"' Returning fail code...\n"<< std::endl;
        success = false;
      }
    }
    if ( p0 ) std::cerr<<_fo<<"Done. Return: '"<<success<<"'\n"<<std::endl;
  };
  ~Idat_Pair() {};
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                 Validate Idat Pair:: Scalar/Vector Variables::
   *                         
   *  validate_pair()
   *     - validate_idat_version()
   *     - validate_address_count()
   *     - validate_Sentrix_Idx_count()
   *     - validate_sentrix_pos_count()
   *     - validate_chip_format()
   *     - validate_address_order()
   *     
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  _Y validate_idat_version( idat::Idat& grn, idat::Idat& red,
                            const _I vb = 0, const _I vt = 4, const _I tc = 0,
                            const _S ft="validate_idat_version" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    if ( grn.get_Scalar_I(IDATVERSION) != red.get_Scalar_I(IDATVERSION) ) {
      std::cerr <<_fe<<IDATVERSION<<" do NOT Match!\n"
                <<_fe<<"\t Grn "<<IDATVERSION<<" = '"<<grn.get_Scalar_I(IDATVERSION)<<"'\n"
                <<_fe<<"\t Red "<<IDATVERSION<<" = '"<<red.get_Scalar_I(IDATVERSION)<<"'\n"
                << std::endl;
      return( false );
    }
    return( true );
  }; // validate_idat_version()
  
  _Y validate_address_count( idat::Idat& grn, idat::Idat& red,
                             const _I vb = 0, const _I vt = 4, const _I tc = 0,
                             const _S ft="validate_address_count" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    if ( grn.get_Scalar_I(ADDRESS_CNT) != red.get_Scalar_I(ADDRESS_CNT) ) {
      std::cerr <<_fe<<"Idat "<<ADDRESS_CNT<<" do NOT Match!\n"
                <<_fe<<"\t Grn "<<ADDRESS_CNT<<" = '"<<grn.get_Scalar_I(ADDRESS_CNT)<<"'\n"
                <<_fe<<"\t Red "<<ADDRESS_CNT<<" = '"<<red.get_Scalar_I(ADDRESS_CNT)<<"'\n"
                << std::endl;
      return( false );
    }
    return( true );
  }; // validate_address_count()
  
  _Y validate_Sentrix_Idx_count( idat::Idat& grn, idat::Idat& red,
                                 const _I vb = 0, const _I vt = 4, const _I tc = 0,
                                 const _S ft="validate_Sentrix_Idx_count" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    if ( grn.get_Scalar_S(SENTRIX_ID) != red.get_Scalar_S(SENTRIX_ID) ) {
      std::cerr <<_fe<<"Idat "<<SENTRIX_ID<<" do NOT Match!\n"
                <<_fe<<"\t Grn "<<SENTRIX_ID<<" = '"<<grn.get_Scalar_S(SENTRIX_ID)<<"'\n"
                <<_fe<<"\t Red "<<SENTRIX_ID<<" = '"<<red.get_Scalar_S(SENTRIX_ID)<<"'\n"
                << std::endl;
      return( false );
    }
    return( true );
  }; // validate_Sentrix_Idx_count()
  
  _Y validate_sentrix_pos_count( idat::Idat& grn, idat::Idat& red,
                                 const _I vb = 0, const _I vt = 4, const _I tc = 0,
                                 const _S ft="validate_sentrix_pos_count" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    if ( grn.get_Scalar_S(SENTRIX_POS) != red.get_Scalar_S(SENTRIX_POS) ) {
      std::cerr <<_fe<<"Idat "<<SENTRIX_POS<<" do NOT Match!\n"
                <<_fe<<"\t Grn "<<SENTRIX_POS<<" = '"<<grn.get_Scalar_S(SENTRIX_POS)<<"'\n"
                <<_fe<<"\t Red "<<SENTRIX_POS<<" = '"<<red.get_Scalar_S(SENTRIX_POS)<<"'\n"
                << std::endl;
      return( false );
    }
    return( true );
  }; // validate_sentrix_pos_count()
  
  _Y validate_chip_format( idat::Idat& grn, idat::Idat& red,
                           const _I vb = 0, const _I vt = 4, const _I tc = 0,
                           const _S ft="validate_chip_format" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    if ( grn.get_Scalar_S(CHIP_FORMAT) != red.get_Scalar_S(CHIP_FORMAT) ) {
      std::cerr <<_fe<<"Idat "<<CHIP_FORMAT<<" do NOT Match!\n"
                <<_fe<<"\t Grn "<<CHIP_FORMAT<<" = '"<<grn.get_Scalar_S(CHIP_FORMAT)<<"'\n"
                <<_fe<<"\t Red "<<CHIP_FORMAT<<" = '"<<red.get_Scalar_S(CHIP_FORMAT)<<"'\n"
                << std::endl;
      return( false );
    }
    return( true );
  }; // validate_chip_format()
  
  _Y validate_address_order( idat::Idat& grn, 
                             idat::Idat& red,
                             const _S& add_key = ADD_STR,
                             const _I vb = 0, const _I vt = 4, const _I tc = 0,
                             const _S ft="validate_address_order" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    _Y success = true;
    if ( grn.get_Vector_I(add_key).size() == 0 ) {
      std::cerr <<_fe<<"Idat Grn Address Vector is Empty!\n"<<std::endl;
      success = false;
    }
    if ( red.get_Vector_I(add_key).size() == 0 ) {
      std::cerr <<_fe<<"Idat Red Address Vector is Empty!\n"<<std::endl;
      success = false;
    }
    if (!success ) return( success );
    
    _I idx = identical_vecs( grn.get_Vector_I(add_key), red.get_Vector_I(add_key) );
    if ( idx != 0 ) {
      std::cerr <<_fe<<"Idat Addresses(idx="<<idx<<") do NOT Match!\n"
                <<_fe<<"\t Grn["<<idx-1<<"] = '"<<grn.get_Vector_I(add_key)[idx-1]<<"'\n"
                <<_fe<<"\t Red["<<idx-1<<"] = '"<<red.get_Vector_I(add_key)[idx-1]<<"'\n"
                <<std::endl;
      return( false );
    }
    return( true );
  }; // validate_address_order()
  
  _Y validate_pair( idat::Idat& grn, idat::Idat& red,
                    const _I vb = 0,const _I vt = 4,const _I tc = 0,
                    const _S ft="validate_pair" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    
    if ( p0 ) std::cerr<<_fo<<"Comparing Idats:: Scalar Variables"<<std::endl;
    
    // Compare:: Scalars [ IdatVersion ]
    if ( !this->validate_idat_version( grn,red, vb,vt+1,tc ) ) return( false );
    
    // Compare:: Scalars [ Address_Cnt ]
    if ( !this->validate_address_count( grn,red, vb,vt+1,tc ) ) return( false );
    
    // Compare:: Scalars [ Sentrix_Idx ]
    if ( !this->validate_Sentrix_Idx_count( grn,red, vb,vt+1,tc ) ) return( false );
    
    // Compare:: Scalars [ Sentrix_Pos ]
    if ( !this->validate_sentrix_pos_count( grn,red, vb,vt+1,tc ) ) return( false );
    
    // Compare:: Scalars [ Chip_Format ]
    if ( !this->validate_chip_format( grn,red, vb,vt+1,tc ) ) return( false );
    
    if ( p0 ) std::cerr<<_fo<<"Comparing Idats:: Vector Variables"<<std::endl;
    
    // Compare:: Vectors [ Validate All Tango Addresses are in Order ]
    if ( !this->validate_address_order( grn,red, "Add", vb,vt+1,tc ) ) return( false );
    
    if ( p0 ) std::cerr<<_fo<<"Idat Pair Validation Successful. Done..\n"<<std::endl;
    
    return( true );
  }; // validate_pair()
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                        Set Best Matching Manifest::
   *
   *  set_manifest()
   *     - set_best_manfest()
   *     - set_pval()
   *     - set_beta()
   *     - set_manfest_data()
   *
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  _Y set_best_manfest( std::vector<_I>& addU_man_vec, std::vector<_I>& addM_man_vec,
                       std::vector<_I>& idat_add_vec, std::vector<_S>& name_man_vec,
                       
                       const _I vb = 0, const _I vt = 1, const _I tc = 0,
                       const _S ft="set_best_manfest" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    const _Y p3 = vb > vt + 3;
    
    _Y success = true;
    
    this->man_top_pair.first  = "";
    this->man_top_pair.second = 0.0;
    
    _I idat_add_cnt = idat_add_vec.size();
    _I cgns_man_cnt = addU_man_vec.size();
    
    std::map<_S,_I> man_cnt_map;
    std::map<_S,_I> man_tot_map;
    std::map<_I,std::vector<_S> > man_add_map;
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     * 
     *                  Initialize Manifest Address/Count Maps::
     *                 
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    for ( std::size_t ii=0; ii<cgns_man_cnt; ii++ ) {
      man_cnt_map[ name_man_vec[ii] ] = 0;
      
      if ( addU_man_vec[ii] != 0 ) {
        man_tot_map[ name_man_vec[ii] ]++;
        man_add_map[ addU_man_vec[ii] ].push_back( name_man_vec[ii] );
      }
      if ( addM_man_vec[ii] != 0 ) {
        man_tot_map[ name_man_vec[ii] ]++;
        man_add_map[ addM_man_vec[ii] ].push_back( name_man_vec[ii] );
      }
    }
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     * 
     *                 Count Manifest Support from Idat Addresses::
     *                 
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    _S max_key = "";
    _d max_per = 0.0;
    std::vector<_I> mis_vec;
    _I mat_cnt = 0; _I mis_cnt = 0; _I max_cnt = 0;
    for ( std::size_t ii=0; ii<idat_add_cnt; ii++ ) {
      _I cur_add = idat_add_vec[ii];
      
      if ( man_add_map.find(cur_add) != man_add_map.end() ) {
        for ( _I jj=0; jj<man_add_map[cur_add].size(); jj++ ) {
          man_cnt_map[ man_add_map[cur_add][jj] ]++;
          if ( man_cnt_map[ man_add_map[cur_add][jj] ] > max_cnt ) {
            max_cnt = man_cnt_map[ man_add_map[cur_add][jj] ];
            max_key = man_add_map[ cur_add ][jj];
          }
        }
        mat_cnt++;
      } else {
        // Keep track of Tango Addresses that are not in the manifest::
        mis_vec.push_back( cur_add );
        mis_cnt++;

        if ( p3 && mis_cnt < 5 ) std::cerr
        <<_fo<<"Background MIS ADDRESS["<<mis_cnt<<"]:: '"<<cur_add<<"'"
        <<std::endl;
      }
    }
    
    _S _xo = _fo;
    //
    // NOT Sure why we need this check:  mis_cnt != 0
    //
    // if ( mis_cnt != 0 && man_tot_map[max_key] != 0 ) {
    if ( man_tot_map[max_key] != 0 ) {
      max_per = (_d) 100 * max_cnt / (_d) man_tot_map[max_key];
      this->man_top_pair.first  = max_key;
      this->man_top_pair.second = max_per;
    } else {
      _xo = _fe;
      success = false;
    }
    
    if ( p0 ) std::cerr
      <<_xo<<"Done. Manifest Top Match Results:: success='"<<success<<"'\n"
      <<_xo<<"\t  idat_add_cnt: '"<<idat_add_cnt<<"'\n"
      <<_xo<<"\t  cgns_man_cnt: '"<<cgns_man_cnt<<"'\n"
      <<_xo<<"\t       mat_cnt: '"<<mat_cnt<<"'\n"
      <<_xo<<"\t       mis_cnt: '"<<mis_cnt<<"'\n"
      <<_xo<<"\t       max_cnt: '"<<max_cnt<<"'\n"
      <<_xo<<"\t       max_per: '"<<set_precision(max_per,2)<<"%'\n"
      <<_xo<<"\t       max_key: '"<<max_key<<"'\n"
      <<std::endl;
    
    return( success );
  }; // set_best_manfest()
  
  _Y set_pval( _d ug, _d mg, _d ur, _d mr, _S col, _d& pval )
  {
    pval = -1.0;
    _Y success = true;
    if ( col[0] == '2' ) {
      pval = std::min( ug,ur );
    } else if ( col[0] == 'G' || col[0] == 'g' ) {
      pval = std::min( ug,mg );
    } else if ( col[0] == 'R' || col[0] == 'r' ) {
      pval = std::min( ur,mr );
    } else {
      success = false;
    }
    return( success );
  }; // set_pval()
  
  _Y set_beta( _d ug, _d mg, _d ur, _d mr, _S col, _d& beta )
  {
    beta = -1.0;
    _Y success = true;
    if ( col[0] == '2' ) {
      if ( ug+ur != 0 ) beta = ug / ( ug+ur );
    } else if ( col[0] == 'G' || col[0] == 'g' ) {
      if ( ug+mg != 0 ) beta = mg / ( ug+mg );
    } else if ( col[0] == 'R' || col[0] == 'r' ) {
      if ( ur+mr != 0 ) beta = mr / ( ur+mr );
    } else {
      success = false;
    }
    return( success );
  }; // set_beta()
  
  // _Y set_mask( _d ug, _d mg, _d ur, _d mr, _S col, _d& beta )
  // {
  //   beta = -1.0;
  //   _Y success = true;
  //   if ( col[0] == '2' ) {
  //     if ( ug+ur != 0 ) beta = ug / ( ug+ur );
  //   } else if ( col[0] == 'G' ) {
  //     if ( ug+mg != 0 ) beta = mg / ( ug+mg );
  //   } else if ( col[0] == 'R' ) {
  //     if ( ur+mr != 0 ) beta = mr / ( ur+mr );
  //   } else {
  //     success = false;
  //   }
  //   return( success );
  // }; // set_mask()
  
  _Y set_manfest_data( idat::Idat& grn, idat::Idat& red,
                       std::vector<_I>& pval_add_vec,
                       std::vector<_I>& addU_man_vec, std::vector<_I>& addM_man_vec,
                       std::vector<_S>& cgns_man_vec, std::vector<_S>& cols_man_vec,
                       std::vector<_S>& name_man_vec, std::vector<_S>& anno_man_vec,
                       std::vector<_S>& chrs_man_vec,
                       const _d min_pval = 0.05,

                       const _I vb = 0, const _I vt = 1, const _I tc = 0,
                       const _S ft="set_manfest_data" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    const _Y p1 = vb > vt + 1;
    const _Y p2 = vb > vt + 2;
    const _Y p10 = vb > vt + 10;
    
    _Y success = true;
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     * 
     *               Build Manifest and Controls from Match Manifest::
     *
     * TBD::
     *   - Update negs addresses by color
     *   - Duplicate Infinium I G/R order
     *   
     *   - Set sigs/detP/beta/pval [0]
     *   - Set sigs/detP/beta/pval [1]
     *
     *
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    std::vector<_I> idat_add_vec = grn.get_Vector_I(ADD_STR);
    _I idat_add_cnt = idat_add_vec.size();
    _I cgns_man_cnt = cgns_man_vec.size();
    
    if ( p0 ) std::cerr
      <<_fo<<BRK
      <<_fo<<"Setting Manifest Data::\n"
      <<_fo<<"\t  (vb,vt,tc): (" <<vb<<","<<vt<<","<<tc<<")\n"
      <<_fo<<"\t  idat_add_cnt: '"<<idat_add_cnt<<"'\n"
      <<_fo<<"\t  cgns_man_cnt: '"<<cgns_man_cnt<<"'\n"
      <<std::endl;
    
    //
    // Clear Previous Data::
    //
    this->idx_to_pair.clear();
    this->pair_to_idx.clear();
    
    this->key_to_info.clear();
    this->key_to_bool.clear();
    
    this->pair_idx_to_data.clear();
    
    std::vector<_I> pval_grn_vec = pval_add_vec;
    std::vector<_I> pval_red_vec = pval_add_vec;

    _I pval_grn_cnt0 = pval_grn_vec.size();
    _I pval_red_cnt0 = pval_red_vec.size();
    
    _I man_bad_cnt = 0;
    _I man_mis_cnt = 0;
    for ( std::size_t idx=0; idx<cgns_man_cnt; idx++ ) {
      _Y valid_prb = true;
      
      // Skip Non-Target Genome::
      if ( name_man_vec[idx].compare(this->man_top_pair.first) != 0 &&
           name_man_vec[idx].compare("ctl") != 0 ) continue;
      
      // This should never happen:: Should throw an error!!!
      if ( addU_man_vec[idx] == 0 ) man_bad_cnt++;
      if ( addU_man_vec[idx] == 0 ) continue;

      // Original tango_pair
      std::pair<_I,_I> tango_pair( addU_man_vec[idx],addM_man_vec[idx] );
      
      // Symetric Map:: idx_to_pair <=> pair_to_idx
      // this->pair_col_to_mask[tango_pair][col] = mask;
      // this->pair_to_idx[tango_pair] = cgns_man_vec[idx];
      this->pair_to_idx[tango_pair] = idx;
      this->idx_to_pair.push_back( tango_pair );

      /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
       * 
       *                             Assign Sig Data::
       *                 
       * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
      
      std::pair<_S,_S> keyUG;
      std::pair<_S,_S> keyMG;
      std::pair<_d,_d> dataG;
      if ( valid_prb ) valid_prb = grn.get_Select_Data_Pair(
        SIG_STR, tango_pair,keyUG,keyMG,dataG, vb,vt+10,tc );
      
      std::pair<_S,_S> keyUR;
      std::pair<_S,_S> keyMR;
      std::pair<_d,_d> dataR;
      if ( valid_prb ) valid_prb = red.get_Select_Data_Pair(
        SIG_STR, tango_pair,keyUR,keyMR,dataR, vb,vt+10,tc );
      
      // Ensure Tango Pair Exists in Idats::
      if ( !valid_prb ) man_mis_cnt++;
      if ( !valid_prb ) continue;
      // Not sure the correct response: either skip (above), or fail (below)
      //   Also, don't like this whole pair<_S,_S> set up. Deal with it for now...
      // if ( !valid_prb ) {
      //   std::cerr <<_fe<<"Failed to find pair(idx="<<idx<<") in get_Select_Data_Pair()!\n"
      //             <<_fe<<"\t key[ UG,MG,UR,MR ] = ["<<keyUG.second<<", "<<keyMG.second<<", "<<keyUR.second<<", "<<keyMR.second<<"]\n"
      //             <<_fe<<"\t   tango_pair.first = '"<<tango_pair.first<<"'\n"
      //             <<_fe<<"\t  tango_pair.second = '"<<tango_pair.second<<"'\n"
      //             <<std::endl;
      //   return( false );
      // }

      this->pair_idx_to_data[ keyUG ].push_back( dataG.first );
      this->pair_idx_to_data[ keyMG ].push_back( dataG.second );
      
      this->pair_idx_to_data[ keyUR ].push_back( dataR.first );
      this->pair_idx_to_data[ keyMR ].push_back( dataR.second );
      
      // Set Relative CGN Index (probably slicker way to do this...)
      _I cgn_idx = this->pair_idx_to_data[ keyMR ].size() - 1;

      // Set Beta Value::
      std::pair<_S,_S> keyBeta( BETA_STR, "0" );
      _d beta = dNEG_ONE;
      if ( success ) success = set_beta(
        dataG.first,dataG.second, dataR.first,dataR.second,
        cols_man_vec[idx], beta );
      
      if ( !success ) {
        std::cerr
        <<_fe<<"Failed set_beta idx=["<<idx<<"]::\n"
        <<std::endl;
        success = false;
        break;
      }
      
      // this->idx_to_beta.push_back(beta);
      this->pair_idx_to_data[keyBeta].push_back(beta);
      
      /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
       * 
       *                             Assign Det Data::
       *                 
       * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
      
      // _S DETN_STR = DET_STR+"0";
      if ( valid_prb ) valid_prb = grn.get_Select_Data_Pair(
        DET_STR, tango_pair,keyUG,keyMG,dataG, vb,vt+10,tc );
      
      if ( valid_prb ) valid_prb = red.get_Select_Data_Pair(
        DET_STR, tango_pair,keyUR,keyMR,dataR, vb,vt+10,tc );
      
      // Ensure Tango Pair Exists in Idats::
      if ( !valid_prb ) man_mis_cnt++;
      if ( !valid_prb ) continue;
      
      keyUG.first = keyUG.first+"_0";
      keyMG.first = keyMG.first+"_0";
      this->pair_idx_to_data[ keyUG ].push_back( dataG.first );
      this->pair_idx_to_data[ keyMG ].push_back( dataG.second );
      
      keyUR.first = keyUR.first+"_0";
      keyMR.first = keyMR.first+"_0";
      this->pair_idx_to_data[ keyUR ].push_back( dataR.first );
      this->pair_idx_to_data[ keyMR ].push_back( dataR.second );
      
      // Set Pval Value::
      // std::pair<_S,_S> keyPval( PVAL_STR,cols_man_vec[idx] );
      std::pair<_S,_S> keyPval( PVAL_STR,"0" );
      _d pval = dNEG_ONE;
      
      if ( success ) success = set_pval(
        dataG.first,dataG.second, dataR.first,dataR.second,
        cols_man_vec[idx], pval );
      
      if ( !success ) {
        std::cerr
        <<_fe<<"Failed set_pval idx=["<<idx<<"]::\n"
        <<std::endl;
        success = false;
        break;
      }
      
      // this->idx_to_pval.push_back(pval);
      this->pair_idx_to_data[keyPval].push_back(pval);
      
      /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
       * 
       *               Add Re-calibration Background Intensities::
       *               
       *  TBD:: For now using 2 x min_pval; should parameterize instead...
       *
       *  NOTES:: Logic Seems to be working...
       *
       *   - Using 2.0 * min_pval::
       *   [set_manfest_data]: 		Background Recalibration Stats::
       *   [set_manfest_data]: 			 Address Count.0[G/R] = [880,880]
       *   [set_manfest_data]: 			 Address Count.1[G/R] = [15346,13419]
       *   
       *   - Using 3.0 * min_pval::
       *   [set_manfest_data]: 		Background Recalibration Stats::
       *   [set_manfest_data]: 			 Address Count.0[G/R] = [880,880]
       *   [set_manfest_data]: 			 Address Count.1[G/R] = [13258,12207]
       *
       *   - Using 4.0 * min_pval
       *   [set_manfest_data]: 		Background Recalibration Stats::
       *   [set_manfest_data]: 			 Address Count.0[G/R] = [880,880]
       *   [set_manfest_data]: 			 Address Count.1[G/R] = [12535,10891]
       *
       * TBD:: Count Outliers Pre/Post...
       *
       * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
      
      // _d min_pval2 = 2.0 * min_pval;
      // _d min_pval2 = 3.0 * min_pval;
      // _d min_pval2 = 4.0 * min_pval;
      // _d min_pval2 = 0.90;
      // _d min_pval2 = 1.00;
      _d min_pval2 = 0.95;
      if ( tango_pair.first != 0 && this->pair_idx_to_data[ keyUG ][cgn_idx] > min_pval2 ) 
        pval_grn_vec.push_back(tango_pair.first);
      if ( tango_pair.first != 0 && this->pair_idx_to_data[ keyUR ][cgn_idx] > min_pval2 ) 
        pval_red_vec.push_back(tango_pair.first);
      
      if ( tango_pair.second != 0 && this->pair_idx_to_data[ keyMG ][cgn_idx] > min_pval2 ) 
        pval_grn_vec.push_back(tango_pair.second);
      if ( tango_pair.second != 0 && this->pair_idx_to_data[ keyMR ][cgn_idx] > min_pval2 ) 
        pval_red_vec.push_back(tango_pair.second);
      
      /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
       * 
       *                          Assign Manifest Info::
       *                 
       * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
      
      // Idx to Meta Data::
      _S prb_din = get_probe_type( cgns_man_vec[idx] );
      this->key_to_info[PROBE_ID].push_back( cgns_man_vec[idx] );
      this->key_to_info[MANIFEST].push_back( name_man_vec[idx] );
      this->key_to_info[COL_STR].push_back( cols_man_vec[idx] );
      this->key_to_info[COL_REF].push_back( cols_man_vec[idx] );
      this->key_to_info[PROBE_DIN].push_back( prb_din );
      this->key_to_info[CHROMOSOME].push_back( chrs_man_vec[idx] );
      this->key_to_info[ANNOTATION].push_back( anno_man_vec[idx] );

      // Keep track of dins(Probe_Type) to idx::
      this->dins_to_idx[prb_din].push_back(cgn_idx);
      
      // Add Default Mask
      _Y is_masked = false;
      if ( this->pair_idx_to_data[keyPval][cgn_idx] > min_pval ) is_masked = true;
      this->key_to_bool[MASK_STR+"_0"].push_back( is_masked );
      
      // Record Summary Stats for Auto Sample Sheet::
      std::pair<_S,_S> keyStat( prb_din,"0" );
      this->dins_to_sum[keyStat]["Total_Cnt"]++;
      if ( !is_masked ) this->dins_to_sum[keyStat]["Pass_Cnt"]++;

      // if ( cgn_idx < 20 ) {
      //  std::cerr <<_fe<<"pair_idx_to_data["
      //            <<keyPval.first<<", "
      //            <<keyPval.second<<"]["<<cgn_idx<<",idx="<<idx<<"] = "
      //            <<this->pair_idx_to_data[keyPval][cgn_idx]<<" > "
      //            <<min_pval<<" = '"<<is_masked<<"'"
      //            <<std::endl;
      // }
      // if ( idx > 10 ) break;
    }
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     * 
     *            Recalibrate DetP[1] with Updated Neg Pval Address::
     *
     *  TBD:: Remove Outliers...
     * 
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    _I pval_grn_cnt1 = pval_grn_vec.size();
    _I pval_red_cnt1 = pval_red_vec.size();

    if ( p0 ) std::cerr
      <<_fo<<"Background Recalibration Stats::\n"
      <<_fo<<"\t Address Count.0[G/R] = ["<<pval_grn_cnt0<<","<<pval_red_cnt0<<"]\n"
      <<_fo<<"\t Address Count.1[G/R] = ["<<pval_grn_cnt1<<","<<pval_red_cnt1<<"]\n"
      <<std::endl;
    
    //
    // TBD: Investigate reset_pvals() function... LEFT OFF HERE
    //
    grn.reset_pvals( pval_grn_vec, false, vb,vt+2,tc );
    red.reset_pvals( pval_red_vec, false, vb,vt+2,tc );
    //
    // Previously the above calll was uncommented...
    //
    // grn.reset_pvals( pval_grn_vec, true, vb,vt+2,tc );
    // red.reset_pvals( pval_red_vec, true, vb,vt+2,tc );
    
    // TBD:: [Done]: in Idat() write function: reset_pval( adds_vec )
    // TBD:: [Done]: Run back through new idats and update ( mask_1,pval_1,dets_1 )
    //
    
    _I pval_grn_cnt2 = pval_grn_vec.size();
    _I pval_red_cnt2 = pval_red_vec.size();
    
    if ( p0 ) std::cerr
    <<_fo<<"P-value Recalibration Stats::\n"
    <<_fo<<"\t Address Count.0[G/R] = ["<<pval_grn_cnt0<<","<<pval_red_cnt0<<"]\n"
    <<_fo<<"\t Address Count.1[G/R] = ["<<pval_grn_cnt1<<","<<pval_red_cnt1<<"]\n"
    <<_fo<<"\t Address Count.1[G/R] = ["<<pval_grn_cnt2<<","<<pval_red_cnt2<<"]\n"
    <<std::endl;
    
    
    
    
    
    
    
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     * 
     *                      This is where the fix is needed. 
     *
     *  TBD:: Process all data in linear mode by index and nothing else...
     *    - Make sure we get G/R data in order...
     *    
     *    - Move prefix loop into Rcpp driver function()
     *      - At least run a full chip at a time
     *    - ecdf()
     *      - read negative signal into ordered stack (or something)
     *      - sored = true
     *      
     *    - Make the code below an option
     *    - Allow different methods
     *      - OOB
     *      - Any probe with pvalue < threshold (currently implemented )
     *    - Add number of background probes to stats
     *    
     * 
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    _I beta_mat_cnt = 0;
    _I beta_mis_cnt = 0;
    _I mask_mis_cnt = 0;
    _I mask_mat_cnt = 0;
    
    _I man_cgn_size = this->idx_to_pair.size();
    
    if ( p1 ) std::cerr
    <<_fo<<"P-value Recalibration Stats:: Analysis Loop\n"
    <<_fo<<"\t man_cgn_size = ["<<man_cgn_size<<"]\n"
    <<std::endl;
    
    if ( true ) {
    for ( std::size_t idx=0; idx<man_cgn_size; idx++ ) {
      _Y valid_prb = true;
      
      if ( p2 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop idx=["<<idx<<"/"<<man_cgn_size<<"]"
      <<std::endl;

      // Original tango_pair
      std::pair<_I,_I> tango_pair( this->idx_to_pair[idx] );
      // std::pair<_I,_I> tango_pair;
      // tango_pair = this->idx_to_pair[idx];
      
      if ( this->pair_to_idx.find(tango_pair) == this->pair_to_idx.end() ) {
        std::cerr
        <<_fe<<"Failed to find tango_pair["<<idx<<"]::\n"
        <<_fe<<"\t  tango_pair.first = ["<<tango_pair.first<<"]\n"
        <<_fe<<"\t tango_pair.second = ["<<tango_pair.second<<"]\n"
        <<std::endl;
        success = false;
        break;
      }

      /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
       * 
       *                             Assign Sig Data::
       *
       *  TBD:: Can clean up and functionalize code below...
       *
       * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
      
      _S prb_din = get_probe_type( cgns_man_vec[idx] );
      
      if ( p10 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop prb_din=["<<prb_din<<"]"
      <<std::endl;
      
      std::pair<_S,_S> keyUG;
      std::pair<_S,_S> keyMG;
      std::pair<_d,_d> dataG;
      if ( valid_prb ) valid_prb = grn.get_Select_Data_Pair(
        SIG_STR, tango_pair,keyUG,keyMG,dataG, vb,vt+10,tc );
      
      if ( p10 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop valid_prbG=["<<valid_prb<<"]"
      <<std::endl;

      std::pair<_S,_S> keyUR;
      std::pair<_S,_S> keyMR;
      std::pair<_d,_d> dataR;
      if ( valid_prb ) valid_prb = red.get_Select_Data_Pair(
        SIG_STR, tango_pair,keyUR,keyMR,dataR, vb,vt+10,tc );
      
      if ( p10 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop valid_prbR=["<<valid_prb<<"]"
      <<std::endl;
      
      // Ensure Tango Pair Exists in Idats::
      if ( !valid_prb ) man_mis_cnt++;
      if ( !valid_prb ) continue;
      
      if ( p10 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop man_mis_cnt=["<<man_mis_cnt<<"]"
      <<std::endl;

      // this->pair_idx_to_data[ keyUG ].push_back( dataG.first );
      // this->pair_idx_to_data[ keyMG ].push_back( dataG.second );
      // 
      // this->pair_idx_to_data[ keyUR ].push_back( dataR.first );
      // this->pair_idx_to_data[ keyMR ].push_back( dataR.second );
      
      // Set Relative CGN Index (probably slicker way to do this...)
      _I cgn_idx = this->pair_idx_to_data[ keyMR ].size() - 1;
      cgn_idx = idx;
      
      if ( p10 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop cgn_idx=["<<cgn_idx<<"]"
      <<std::endl;
      
      // Set Beta Value::
      std::pair<_S,_S> keyBeta( BETA_STR, "1" );
      _d beta = dNEG_ONE;
      if ( success ) success = set_beta(
        dataG.first,dataG.second, dataR.first,dataR.second,
        cols_man_vec[idx], beta );
      
      if ( !success ) {
        std::cerr
        <<_fe<<"Failed set_beta idx=["<<idx<<"]::\n"
        <<std::endl;
        success = false;
        break;
      }
      
      this->pair_idx_to_data[keyBeta].push_back(beta);
      
      if ( p10 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop beta=["<<beta<<", BETA_STR="<<BETA_STR<<"]"
      <<std::endl;

      // Line below was there previously, but I think that's an error::
      // std::pair<_S,_S> preBeta( BETA_STR, "1" );
      std::pair<_S,_S> preBeta( BETA_STR, "0" );
      
      if ( p10 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop beta=["<<beta<<", BETA_STR="<<BETA_STR<<"]\n"
      <<_fo<<"\t\t pair_idx_to_data[preBeta][cgn_idx] = ["<<this->pair_idx_to_data[preBeta][cgn_idx]<<"]"
      <<std::endl;
      
      continue;
      
      if ( p0 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop beta=["<<beta<<", BETA_STR="<<BETA_STR<<"]\n"
      <<_fo<<"\t\t pair_idx_to_data[keyBeta][cgn_idx] = ["<<this->pair_idx_to_data[keyBeta][cgn_idx]<<"]"
      <<std::endl;

      if ( this->pair_idx_to_data[preBeta][cgn_idx] - this->pair_idx_to_data[keyBeta][cgn_idx] != 0 )
        beta_mis_cnt++;
      
      if ( p0 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop beta_mis_cnt=["<<beta_mis_cnt<<"]"
      <<std::endl;
      
      if ( this->pair_idx_to_data[preBeta][cgn_idx] - this->pair_idx_to_data[keyBeta][cgn_idx] == 0 )
        beta_mat_cnt++;
      
      if ( p0 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop beta_mat_cnt=["<<beta_mat_cnt<<"]"
      <<std::endl;

      /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
       * 
       *                             Assign Det Data::
       *                 
       * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
      
      if ( valid_prb ) valid_prb = grn.get_Select_Data_Pair(
        DET_STR, tango_pair,keyUG,keyMG,dataG, vb,vt+10,tc );
      
      if ( p0 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop valid_prb1=["<<valid_prb<<"]"
      <<std::endl;

      if ( valid_prb ) valid_prb = red.get_Select_Data_Pair(
        DET_STR, tango_pair,keyUR,keyMR,dataR, vb,vt+10,tc );
      
      if ( p0 ) std::cerr
      <<_fo<<"\tP-value Recalibration Stats:: Analysis Loop valid_prb2=["<<valid_prb<<"]"
      <<std::endl;
      
      if ( p0 ) std::cerr <<std::endl;

      // Ensure Tango Pair Exists in Idats::
      if ( !valid_prb ) man_mis_cnt++;
      if ( !valid_prb ) continue;
      
      keyUG.first = keyUG.first+"_1";
      keyMG.first = keyMG.first+"_1";
      this->pair_idx_to_data[ keyUG ].push_back( dataG.first );
      this->pair_idx_to_data[ keyMG ].push_back( dataG.second );
      
      keyUR.first = keyUR.first+"_1";
      keyMR.first = keyMR.first+"_1";
      this->pair_idx_to_data[ keyUR ].push_back( dataR.first );
      this->pair_idx_to_data[ keyMR ].push_back( dataR.second );
      
      // Set Pval Value::
      std::pair<_S,_S> keyPval( PVAL_STR,"1" );
      _d pval = dNEG_ONE;
      if ( success) success = set_pval(
        dataG.first,dataG.second, dataR.first,dataR.second,
        cols_man_vec[idx], pval );
      
      if ( !success ) {
        std::cerr
        <<_fe<<"Failed set_pval idx=["<<idx<<"]::\n"
        <<std::endl;
        success = false;
        break;
      }
      
      this->pair_idx_to_data[keyPval].push_back(pval);
      
      // Add Default Mask 
      _Y is_masked = false;
      if ( this->pair_idx_to_data[keyPval][cgn_idx] > min_pval ) is_masked = true;
      this->key_to_bool[MASK_STR+"_1"].push_back( is_masked );
      
      // Record Summary Stats for Auto Sample Sheet::
      std::pair<_S,_S> keyStat( prb_din,"1" );
      this->dins_to_sum[keyStat]["Total_Cnt"]++;
      if ( !is_masked ) this->dins_to_sum[keyStat]["Pass_Cnt"]++;

      if ( this->key_to_bool[MASK_STR+"_0"][cgn_idx] != this->key_to_bool[MASK_STR+"_1"][cgn_idx] )
        mask_mis_cnt++;
      if ( this->key_to_bool[MASK_STR+"_0"][cgn_idx] == this->key_to_bool[MASK_STR+"_1"][cgn_idx] )
        mask_mat_cnt++;
      
      if ( is_masked && cgn_idx < 30 ) {
        std::cerr
        <<_fo<<"IS_MASKED_1::\n"
        <<_fo<<"\t          Probe_ID = '"<<this->key_to_info[PROBE_ID][cgn_idx]<<"'\n"
        <<_fo<<"\t          MANIFEST = '"<<this->key_to_info[MANIFEST][cgn_idx]<<"'\n"
        <<_fo<<"\t           COL_STR = '"<<this->key_to_info[COL_STR][cgn_idx]<<"'\n"
        <<_fo<<"\t           COL_REF = '"<<this->key_to_info[COL_REF][cgn_idx]<<"'\n"
        <<_fo<<"\t         is_masked = '"<<is_masked<<"'\n"
        <<_fo<<"\t     keyPval.first = '"<<keyPval.first<<"'\n"
        <<_fo<<"\t    keyPval.second = '"<<keyPval.second<<"'\n"
        <<_fo<<"\t           cgn_idx = '"<<cgn_idx<<"'\n"
        <<_fo<<"\t          min_pval = '"<<min_pval<<"'\n\n"
        <<_fo<<"\t      beta_mat_cnt = '"<<beta_mat_cnt<<"'\n"
        <<_fo<<"\t      beta_mis_cnt = '"<<beta_mis_cnt<<"'\n"
        <<_fo<<"\t      mask_mat_cnt = '"<<mask_mat_cnt<<"'\n"
        <<_fo<<"\t      mask_mis_cnt = '"<<mask_mis_cnt<<"'\n"
        <<_fo<<"\t  pair_idx_to_data = '"<<this->pair_idx_to_data[keyPval][cgn_idx]<<"'\n"
        <<std::endl;
      }
    }
    }
    
    //
    // Add Passing Percent::
    //
    for ( auto& keyStat : this->dins_to_sum ) {
      this->dins_to_sum[keyStat.first]["Pass_Per"] = set_precision( keyStat.second["Pass_Cnt"] / keyStat.second["Total_Cnt"], 6 );
    }
    
    if ( p1 ) {
      // std::pair<_I,_I> tango_pair( this->idx_to_pair[0] );
      std::pair<_S,_S> keyBeta( BETA_STR, "0" );
      // std::pair<_I,_I> keyBeta( this->idx_to_pair[0] );
      std::cerr <<_fo<<"KEY_BETA[0]=["<<this->pair_idx_to_data[keyBeta][0]<<"], "
                <<" first="<<keyBeta.first<<", "
                <<"second="<<keyBeta.second
                <<std::endl;
      std::pair<_S,_S> preBeta( BETA_STR, "1" );
      // std::pair<_I,_I> preBeta( this->idx_to_pair[1] );
      std::cerr <<_fo<<"PRE_BETA[1]=["<<this->pair_idx_to_data[preBeta][0]<<"], "
                <<" first="<<preBeta.first<<", "
                <<"second="<<preBeta.second<<"\n"
                <<std::endl;
    }
    
    // TBD:: [to_sample_sheet]: Calc Summary Stats (PdP_0, PdP_1, etc.)
    // TBD:: [to_sample_sheet]: Add Run_Time params (Machine_Name, etc.)
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                                    Done::
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    if ( p0 ) std::cerr
      <<_fo<<"Done. Setting Manifest Data:: success='"<<success<<"'\n"
      <<_fo<<"\t beta_mat_cnt = '"<<beta_mat_cnt<<"'\n"
      <<_fo<<"\t beta_mis_cnt = '"<<beta_mis_cnt<<"'\n"
      <<_fo<<"\t mask_mat_cnt = '"<<mask_mat_cnt<<"'\n"
      <<_fo<<"\t mask_mis_cnt = '"<<mask_mis_cnt<<"'\n"
      <<_fo<<"\t  man_bad_cnt = '"<<man_bad_cnt<<"'\n"
      <<_fo<<"\t  man_mis_cnt = '"<<man_mis_cnt<<"'\n"
      <<std::endl;
    
    return( success );
  }; // set_manfest_data()
  
  _Y set_manifest( idat::Idat& grn, idat::Idat& red, const _S& output_path,
                   std::vector<_I>& pval_add_vec,
                   std::vector<_I>& addU_man_vec, std::vector<_I>& addM_man_vec,
                   std::vector<_S>& cgns_man_vec, std::vector<_S>& cols_man_vec,
                   std::vector<_S>& name_man_vec, std::vector<_S>& anno_man_vec,
                   std::vector<_S>& chrs_man_vec,
                   const _d min_pval = 0.05,
                   const _d min_beta = 0.30, const _d max_beta = 0.70,
                   const _d min_perI = 0.05, const _d min_perO = 0.75,
                   const _Y read_bgz  = false, 
                   const _Y write_bgz = false,
                   const _I vb = 0, const _I vt = 1, const _I tc = 0,
                   const _S ft="set_manifest" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    
    _Y success = true;
    
    std::vector<_I> idat_add_vec = grn.get_Vector_I(ADD_STR);
    _I idat_add_cnt = idat_add_vec.size();
    _I cgns_man_cnt = cgns_man_vec.size();
    
    if ( p0 ) std::cerr
      <<_fo<<BRK
      <<_fo<<"Setting Manifest::\n"
      <<_fo<<"\t  (vb,vt,tc): (" <<vb<<","<<vt<<","<<tc<<")\n"
      <<_fo<<"\t  idat_add_cnt: '"<<idat_add_cnt<<"'\n"
      <<_fo<<"\t  cgns_man_cnt: '"<<cgns_man_cnt<<"'\n"
      <<std::endl;
    
    if ( success ) 
      success = set_best_manfest( addU_man_vec, addM_man_vec, 
                                  idat_add_vec, name_man_vec, 
                                  vb,vt+1,tc+1 );
    
    // After Best Manifest has been detected we can now expand the Infinium I
    //  probes.
    //  - Scan through data to calculate total expansion. Skip this just double
    //  - Re-allocate vectors
    //  - Assign expanded set
    // Then everything else runs the same...
    //
    
    // TBD:: [Done]: Reverse insert everything into the expanded vectors...
    // TBD:: Functionalize the code below::
    //
    std::size_t inc = cols_man_vec.size();
    for ( auto& col : cols_man_vec ) if ( col[0]=='R' || col[0]=='G' ) inc++;

    std::vector<_I> addU_man2_vec(inc,0);
    std::vector<_I> addM_man2_vec(inc,0);
    std::vector<_S> cgns_man2_vec(inc,"");
    std::vector<_S> cols_man2_vec(inc,"2");
    std::vector<_S> name_man2_vec(inc,"");
    std::vector<_S> anno_man2_vec(inc,"");
    std::vector<_S> chrs_man2_vec(inc,"");
    
    _S colG = "G";
    _S colR = "R";
    _I colS = cols_man_vec.size() - 1;
    std::size_t jj = inc - 1;
    
    // TBD:: Remove Extra::
    // std::cerr<<_fo<<BRK
    //          <<_fo<<"  inc = ["<<inc<<"]='"<<cols_man2_vec[inc-1]<<"'\n"
    //          <<_fo<<"ii,jj = ["<<colS<<","<<jj<<"]='"<<cgns_man_vec[0]<<"'\n"
    //          <<_fo<<"ii,jj = ["<<colS<<","<<jj<<"]='"<<cgns_man_vec[colS]<<"'\n"
    //          <<_fo<<BRK<<std::endl;

    _Y expand_inf1 = false;
    if ( expand_inf1 ) {
      for ( std::size_t ii=colS; ii>0; ii-- ) {
        if ( jj <= 0 ) {
          std::cerr<<_fe<<BRK
                   <<_fe<<"JJ <= 0\n"
                   <<_fe<<"ii,jj = ["<<colS<<","<<jj<<"]='"<<cgns_man_vec[0]<<"'\n"
                   <<_fe<<"ii,jj = ["<<colS<<","<<jj<<"]='"<<cgns_man_vec[colS]<<"'\n"
                   <<_fe<<BRK<<std::endl;
          success = false;
          break;
        }
        
        _B col = cols_man_vec[ii][0];
        
        addU_man2_vec[jj] = addU_man_vec[ii];
        addM_man2_vec[jj] = addM_man_vec[ii];
        
        cgns_man2_vec[jj] = cgns_man_vec[ii];
        cols_man2_vec[jj] = cols_man_vec[ii];
        name_man2_vec[jj] = name_man_vec[ii];
        anno_man2_vec[jj] = anno_man_vec[ii];
        chrs_man2_vec[jj] = chrs_man_vec[ii];
        
        if ( col=='2' ) {
        } else if ( col=='G' || col=='R' ) {
          colG = "G";
          colR = "R";
          if ( col=='G' ) colR = "r";
          if ( col=='R' ) colG = "g";
          
          cgns_man2_vec[jj] = cgns_man_vec[ii]+"_"+colR;
          cols_man2_vec[jj] = colR;
          
          jj--;
          addU_man2_vec[jj] = addU_man_vec[ii];
          addM_man2_vec[jj] = addM_man_vec[ii];
          
          cgns_man2_vec[jj] = cgns_man_vec[ii]+"_"+colG;
          cols_man2_vec[jj] = colG;
          name_man2_vec[jj] = name_man_vec[ii];
          anno_man2_vec[jj] = anno_man_vec[ii];
          chrs_man2_vec[jj] = chrs_man_vec[ii];
          
        } else {
          success = false;
          std::cerr<<_fe<<BRK
                   <<_fe<<"Unrecognized   cgn["<<ii<<","<<jj<<"]='"<<cgns_man_vec[ii]<<"'\n"
                   <<_fe<<"Unrecognized  addU["<<ii<<","<<jj<<"]='"<<addU_man_vec[ii]<<"'\n"
                   <<_fe<<"Unrecognized  addM["<<ii<<","<<jj<<"]='"<<addM_man_vec[ii]<<"'\n"
                   <<_fe<<BRK<<std::endl;
          std::cerr<<_fe<<BRK
                   <<_fe<<"Unrecognized color["<<ii<<","<<jj<<"]='"<<col<<"'\n"
                   <<_fe<<BRK<<std::endl;
          break;
        }
        jj--;
      }
      if ( success )
        success = set_manfest_data( grn,red,
                                    pval_add_vec,
                                    addU_man2_vec, addM_man2_vec,
                                    cgns_man2_vec, cols_man2_vec,
                                    name_man2_vec, anno_man2_vec, chrs_man2_vec,
                                    min_pval,
                                    vb,vt+1,tc+1 );
    } else {
      if ( success )
        success = set_manfest_data( grn,red,
                                    pval_add_vec,
                                    addU_man_vec, addM_man_vec,
                                    cgns_man_vec, cols_man_vec,
                                    name_man_vec, anno_man_vec, chrs_man_vec,
                                    min_pval,
                                    vb,vt+1,tc+1 );
    }
    

    // if ( success )
    //   success = set_manfest_data( grn,red,
    //                               addU_man_vec, addM_man_vec,
    //                               cgns_man_vec, cols_man_vec,
    //                               name_man_vec, anno_man_vec, chrs_man_vec,
    //                               vb,vt+1,tc+1 );

    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                                   Done::
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    if ( p0 ) std::cerr<<_fo<<"Done.\n"<<_fo<<BRK<<std::endl;
    
    return( true );
  }; // set_manifest()
  
  //
  // This is NOT Functional::
  //
  _Y set_summary( idat::Idat& grn, idat::Idat& red, const _S& output_path,
                   std::vector<_I>& pval_add_vec,
                   std::vector<_I>& addU_man_vec, std::vector<_I>& addM_man_vec,
                   std::vector<_S>& cgns_man_vec, std::vector<_S>& cols_man_vec,
                   std::vector<_S>& name_man_vec, std::vector<_S>& anno_man_vec,
                   std::vector<_S>& chrs_man_vec,
                   const _d min_pval = 0.05,
                   const _d min_beta = 0.30, const _d max_beta = 0.70,
                   const _d min_perI = 0.05, const _d min_perO = 0.75,
                   const _Y read_bgz  = false, 
                   const _Y write_bgz = false,
                   const _I vb = 0, const _I vt = 1, const _I tc = 0,
                   const _S ft="set_summary" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    
    _Y success = true;
    
    // this->key_to_bool[MASK_STR+"_0"].push_back( is_masked );
    _I data_man_cnt = this->key_to_bool[MASK_STR+"_0"].size();

    if ( p0 ) std::cerr
      <<_fo<<BRK
      <<_fo<<"Setting Manifest::\n"
      <<_fo<<"\t  (vb,vt,tc): (" <<vb<<","<<vt<<","<<tc<<")\n"
      <<_fo<<"\t  data_man_cnt: '"<<data_man_cnt<<"'\n"
      <<std::endl;
    
    // TBD:: End idx should be set...
    // _I grp_max = 1;
    
    // std::map< _S, std::vector<_I> > dins_to_idx;
    // std::map< std::pair<_S,_I>,_d > dins_to_sum;
    // for ( auto& din_dat : this->dins_to_idx ) {
    //   _S din_key = din_dat.first;
    //   for ( auto& idx : this->dins_to_idx[din_dat.first] ) {
    //     for ( auto& mask this->key_to_bool )
    //   }
    // }
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                                   Done::
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    if ( p0 ) std::cerr<<_fo<<"Done.\n"<<_fo<<BRK<<std::endl;
    
    return( success );
  }; // set_summary()

  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                      Load Idat Pair:: Version 3.0
   * 
   * parse_idat_pair()
   *   - validate_pair()
   *   - set_manifest()
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  _Y parse_idat_pair( 
      const _S prefix_path, const _S& output_path,
      std::vector<_I>& pval_add_vec,
      std::vector<_I>& addU_man_vec, std::vector<_I>& addM_man_vec,
      std::vector<_S>& cgns_man_vec, std::vector<_S>& cols_man_vec,
      std::vector<_S>& name_man_vec, std::vector<_S>& anno_man_vec,
      std::vector<_S>& chrs_man_vec,
      const _d min_pval = 0.05, 
      const _d min_beta = 0.30, const _d max_beta = 0.70,
      const _d min_perO = 0.75, const _d min_perI = 0.05,
      const _Y read_bgz  = false, 
      const _Y write_bgz = false,
      const _Y rm_pval_outliers = false,
      const _I vb = 0, const _I vt = 1, const _I tc = 0,
      const _S ft="parse_idat_pair" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    const _Y p4 = vb > vt + 4;
    
    _Y success = true;
    
    _S grn_path = prefix_path+"_Grn.idat.gz";
    _S red_path = prefix_path+"_Red.idat.gz";
    _I pval_cnt = pval_add_vec.size();
    _I cgns_cnt = cgns_man_vec.size();
    _I chrs_cnt = chrs_man_vec.size();
    
    if ( p0 ) std::cerr
      <<_fo<<BRK
      <<_fo<< "Statring...\n"
      <<_fo<< "\t  (vb,vt,tc): ("<<vb<<","<<vt<<","<<tc<<")\n"
      <<_fo<< "\t   prefix_path = '"<<prefix_path<<"'\n"
      <<_fo<< "\t   output_path = '"<<output_path<<"'\n\n"
      <<_fo<< "\t      pval_cnt = '"<<pval_cnt<<"'\n"
      <<_fo<< "\t      cgns_cnt = '"<<cgns_cnt<<"'\n"
      <<_fo<< "\t      chrs_cnt = '"<<chrs_cnt<<"'\n\n"
      <<_fo<< "\t      read_bgz = '"<<read_bgz<<"'\n"
      <<_fo<< "\t     write_bgz = '"<<write_bgz<<"'\n"
      <<std::endl;
    
    if ( output_path.size() > 0 ) 
      if ( !build_path( output_path ) ) return( false );
      if ( p4 ) std::cerr
        <<_fo<<"Built Out Path = '"<<output_path<<"'"<<std::endl;
      
      // Set local variables and then assign to this if successful...
      if ( success ) success = legal_file( grn_path, vb,vt+8,tc );
      if ( success ) success = legal_file( red_path, vb,vt+8,tc );
      
      // Load Idats::
      idat::Idat grn_idat( grn_path, 'G', pval_add_vec,rm_pval_outliers, vb,vt+5,tc );
      idat::Idat red_idat( red_path, 'R', pval_add_vec,rm_pval_outliers, vb,vt+5,tc );
      
      // Validate Idat Pair:: Scalar/Vector Variables::
      if ( success ) 
        success = this->validate_pair( grn_idat,red_idat, vb,vt+4,tc );
      
      // Set/Build Manifest::
      if ( success ) 
        success = this->set_manifest( grn_idat, red_idat, output_path,
                                      pval_add_vec,
                                      addU_man_vec, addM_man_vec,
                                      cgns_man_vec, cols_man_vec,
                                      name_man_vec, anno_man_vec,
                                      chrs_man_vec, 
                                      min_pval, min_beta, max_beta,
                                      min_perO, min_perI, 
                                      read_bgz, write_bgz,
                                      vb,vt+1,tc );
      //
      // Set Summary Stats::
      //
      
      //
      // TBD:: Write function to set sentrix_fields
      //
      if ( success ) {
        this->output_path = output_path;
        this->sentrix_info[ SENTRIX_NAME ] = 
          grn_idat.get_Scalar_S( SENTRIX_ID ) + "_" +
          grn_idat.get_Scalar_S( SENTRIX_POS );
        
        this->sentrix_info[ SENTRIX_ID ]  = grn_idat.get_Scalar_S( SENTRIX_ID );
        this->sentrix_info[ SENTRIX_POS ] = grn_idat.get_Scalar_S( SENTRIX_POS );
        
        this->sentrix_data[ MIN_PVAL_STR] = min_pval;
        this->sentrix_data[ MIN_BETA_STR] = min_beta;
        this->sentrix_data[ MAX_BETA_STR] = max_beta;
        this->sentrix_data[ MIN_PERI_STR] = min_perI;
        this->sentrix_data[ MIN_PERO_STR] = min_perO;
      }
      
      if ( success ) success = this->to_sample_sheet( _COM, vb,vt+1,tc );
      
      /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
       *                                   Done::
       * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
      
      if ( p0 ) 
        std::cerr<<_fo<<"Done, Success='"<<success<<"'\n"<<_fo<<BRK<<std::endl;
      
      return( success );
  }; // parse_idat_pair()
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                              Access Methods::
   * 
   * TBD:: *** OUT OF DATE ***
   * 
   * Scalar Methods::
   *   get_Scalar_I()
   *   get_Scalar_S()
   *   
   * Vector Methods::
   *   get_Info_Vector()
   *   get_Bool_Vector()
   *   get_Pair_Vector()
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  // Get Scalar Fields::
  // _I get_Field_I( _S field ) { return( this->sI_dat[ field ] ); };
  // _S get_Field_S( _S field ) { return( this->sS_dat[ field ] ); };
  // _d get_Field_d( _S field ) { return( this->sd_dat[ field ] ); };
  
  // Get Vector Fields::
  //  NOTE:: These are not used::
  // std::vector<_d> get_Beta_Vector() { return( this->idx_to_beta ); };
  // std::vector<_d> get_Pval_Vector() { return( this->idx_to_pval ); };
  
  std::vector<_S> get_Info_Vector( _S f1 ) { return( this->key_to_info[ f1 ] ); };
  std::vector<_Y> get_Bool_Vector( _S f1 ) { return( this->key_to_bool[ f1 ] ); };
  std::vector<_d> get_Pair_Vector( _S f1, _S f2 ) { 
    std::pair<_S,_S> key_pair( f1,f2 );
    return( this->pair_idx_to_data[ key_pair] );
  };
  
  // Run Info:: Matrix
  // std::vector< std::vector<_S> > get_runinfo_mat() const { return( this->runinfo_mat ); };
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                           Manipulation Methods::
   * 
   *   set_pvals()
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                               Output Methods::
   * 
   *   to_string()
   *   TBD:: to_sdf()
   *   TBD:: to_bgz()
   *   TBD:: to_sample_sheet()
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  _Y to_sample_sheet( const char sep = _COM,
                      
                      const _I vb = 0, const _I vt = 1, const _I tc = 0,
                      const _S ft="to_sample_sheet" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    const _Y p4 = vb > vt + 4;
    
    _Y success = true;
    
    fast_oss ssH; fast_oss ssD;
    fast_oss sentrixH;
    fast_oss sentrixD;
    
    fast_oss oss;
    for ( auto& d : this->sentrix_info ) {
      oss<<d.first<<sep;
      oss<<d.second<<sep;
    }
    if ( p4 ) std::cerr<<_fo<<"OSS::\n'"<<oss.str()<<"'\n"<<std::endl;
    
    // Sanity Check:: Sentrix Name defined
    if ( this->sentrix_info.find( SENTRIX_NAME ) == this->sentrix_info.end() ) {
      std::cerr<<_fe<<"Failed to Info Field: '"<<SENTRIX_NAME<<"'\n"<<std::endl;
      return( false );
    }

    // Open File Handel::
    //  TBD:: Write gzipped file...
    _S out_file = this->output_path+"/"+this->sentrix_info[ SENTRIX_NAME ]+".sample_sheet.csv";
    std::ofstream OUT_FH ( out_file, std::ios::out );
    
    if ( !OUT_FH.is_open() ) {
      std::cerr<<_fe<<"Failed to open idx file: '"<<out_file<<"'\n"<<std::endl;
      return( false );
    }
    if ( p0 ) std::cerr
      <<_fo<<"Writing Auto Sample Sheet: '"<<out_file<<"'"<<std::endl;
    
    // Build Sentrix Strings:: Info
    for ( auto& d : this->sentrix_info ) {
      sentrixH<<d.first<<sep;
      sentrixD<<d.second<<sep;
    }
    // Write Manifest Info::
    sentrixH<<"Manifest_Name"<<sep;
    sentrixD<<this->man_top_pair.first<<sep;
    sentrixH<<"Manifest_Score"<<sep;
    sentrixD<<this->man_top_pair.second<<sep;
    
    // Build Sentrix Strings:: Data
    for ( auto& d : this->sentrix_data ) {
      sentrixH<<d.first<<sep;
      sentrixD<<d.second<<sep;
    }
    
    // Write Summary Stats::
    for ( auto& keyStat : this->dins_to_sum ) {
      for ( auto& curStat : keyStat.second ) {
        sentrixH<<curStat.first+"_"+keyStat.first.first+"_"+keyStat.first.second<<sep;
        sentrixD<<curStat.second<<sep;
      }
    }
    
    // Write Sentrix Strings to StringStream Buffer::
    ssH <<sentrixH.str(); 
    ssD <<sentrixD.str();
    
    // Write StringStream Buffer to File::
    //  TBD:: Write gzipped file...
    ssH <<"Last_Field"<<std::endl;
    ssD <<"Last_Field"<<std::endl;
    
    OUT_FH<<ssH.str();
    OUT_FH<<ssD.str();
    OUT_FH.close();
    
    if ( p0 ) std::cerr
      <<_fo<<"Done. Writing Auto Sample Sheets!\n"<<std::endl;

    return( success );
  }; // to_sample_sheet()
  
  _S to_string( const _S format = "tab",
                const _S prefix = "[Idat_Pair]:",
                const char sep = ',',
                
                const _I vb = 0, const _I vt = 1, const _I tc = 0,
                const _S ft="to_string" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    // const _Y p0 = vb > vt + 0;
    
    stringstream ss;
    if ( format.compare("tab") == 0 )
    {
      
    }
    
    return( ss.str() );
  }; // to_string()

}; // Idat_Pair()

}; // namespace idat_pair;

#endif /* idat_pair.h */

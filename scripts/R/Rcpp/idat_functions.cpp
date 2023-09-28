
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppCommon.h>

// My headers::
#include "idat.h"

using namespace Rcpp;
using namespace idat;

//
// TBD:: This is pretty much junk and should be removed/replaced::
//
template <typename T>
std::vector<T> get_dB_stats5_vec( const std::vector<T>& X, 
                                  const std::vector<T>& Y, _Y use_abs = false )
{
  if ( !use_abs ) return( stats5_vec( X-Y, false ) );
  std::vector<T> D = X-Y;
  for ( size_t ii=0; ii<D.size(); ii++ ) D[ii] = std::abs(D[ii]);
  return( stats5_vec( D, false ) );
};

// [[Rcpp::export]]
_d pearsoncoeff_rcpp( const std::vector<_d>& X, const std::vector<_d>& Y )
{
  return( pearsoncoeff_vec(X,Y) );
};

// [[Rcpp::export]]
_I pass_count_rcpp( const std::vector<_d>& X, const std::vector<_d>& Y,
                    const _d cut = 0.2, const _i cmp = 0, 
                    const _Y is_abs = false )
{
  return( pass_count( X,Y, cut,cmp,is_abs ) );
};

/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *                    Rcpp to Cpp Local Conversion Function::
 * 
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

void rcpp_to_cpp_vec_dbl( Rcpp::NumericVector& vec_r,
                          std::vector<_d>& vec, _Y rm_na = false ) {
  _I s = vec_r.size();
  
  vec.clear();
  vec.resize( s );
  _I last_idx = 0;
  for ( _I ii = 0; ii < s; ii++ ) {
    if ( rm_na && Rcpp::NumericVector::is_na( vec_r(ii) ) ) continue;
    vec[ii] = vec_r(ii);
    last_idx++;
  }
  if ( rm_na ) vec.resize( last_idx );
}

void rcpp_to_cpp_vec_int( Rcpp::IntegerVector& vec_r,
                          std::vector<int>& vec ) {
  _I s = vec_r.size();
  
  vec.clear();
  vec.resize( s );
  for ( _I ii = 0; ii < s; ii++ ) vec[ii] = vec_r(ii);
}

void rcpp_to_cpp_vec_str( Rcpp::CharacterVector& vec_r,
                          std::vector<std::string>& vec ) {
  _I s = vec_r.size();
  
  vec.clear();
  vec.resize( s );
  for ( _I ii = 0; ii < s; ii++ ) vec[ii] = vec_r(ii);
}

void rcpp_to_cpp_vec_bool( Rcpp::LogicalVector& vec_r,
                           std::vector<bool>& vec ) {
  _I s = vec_r.size();
  
  vec.clear();
  vec.resize( s );
  for ( _I ii = 0; ii < s; ii++ ) vec[ii] = vec_r(ii);
}

/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *                        Loci Variation Function::
 * 
 * loci_variation_mat()
 *   - Description:: All-Against-All same probe across samples, r2/dB and ppp 
 *      (percent passing pval) 
 *      
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

// [[Rcpp::export]]
Rcpp::DataFrame loci_cross_variation_mat( Rcpp::NumericMatrix beta1_mat_r,
                                          Rcpp::NumericMatrix beta2_mat_r,
                                          const _d min_dB   = 0.2,
                                          const _d min_pval = 0.05,
                                          
                                          const unsigned int vb = 0,
                                          const unsigned int vt = 1,
                                          const unsigned int tc = 0,
                                          const std::string ft="loci_cross_variation_mat" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0 = vb > vt + 0;
  
  Rcpp::DataFrame df;
  
  _I mat1_nrow = beta1_mat_r.nrow();
  _I mat1_ncol = beta1_mat_r.ncol();

  _I mat2_nrow = beta2_mat_r.nrow();
  _I mat2_ncol = beta2_mat_r.ncol();
  
  _I min_nrow = std::min( mat1_nrow,mat2_nrow );
  
  CharacterVector rownames1_vec_r = rownames( beta1_mat_r );
  CharacterVector colnames1_vec_r = colnames( beta1_mat_r );

  CharacterVector rownames2_vec_r = rownames( beta2_mat_r );
  CharacterVector colnames2_vec_r = colnames( beta2_mat_r );
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<<"Statring...\n"
    <<_fo<<"\t mat1_nrow: '"<< mat1_nrow<<"'\n"
    <<_fo<<"\t mat1_ncol: '"<< mat1_ncol<<"'\n"
    <<_fo<<"\t mat2_nrow: '"<< mat2_nrow<<"'\n"
    <<_fo<<"\t mat2_ncol: '"<< mat2_ncol<<"'\n\n"
    <<_fo<<"\t  min_nrow: '"<< min_nrow<<"'\n"
    <<std::endl;
  
  /* ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                             Loop Over Rows!
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  for ( size_t row_idx=0; row_idx<min_nrow; row_idx++ ) {
    
    
    
  }
  
  
  /* ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                                  Done!
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  return( df );
};

// [[Rcpp::export]]
Rcpp::DataFrame loci_variation_mat( Rcpp::NumericMatrix beta_mat_r,
                                    Rcpp::NumericMatrix pval_mat_r,
                                    const _d min_dB   = 0.2,
                                    const _d min_pval = 0.05,
                                    
                                    const unsigned int vb = 0,
                                    const unsigned int vt = 1,
                                    const unsigned int tc = 0,
                                    const std::string ft="loci_variation_mat" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0 = vb > vt + 0;
  
  Rcpp::DataFrame df;
  
  _I mat_nrow = beta_mat_r.nrow();
  _I mat_ncol = beta_mat_r.ncol();
  
  CharacterVector rownames_vec_r = rownames( beta_mat_r );
  CharacterVector colnames_vec_r = colnames( beta_mat_r );
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<<"Statring...\n"
    <<_fo<<"\t mat_nrow: '" << mat_nrow<<"'\n"
    <<_fo<<"\t mat_ncol: '" << mat_ncol<<"'\n"
    <<std::endl;
  
  std::vector<_S> cg_vec;
  
  std::vector<_I> prb_tot_cnt_vec;
  std::vector<_I> prb_pas_cnt_vec;
  std::vector<_d> prb_pas_per_vec;
  
  std::vector<_I> cmb_tot_cnt_vec;
  std::vector<_I> cmb_pas_cnt_vec;
  std::vector<_d> cmb_pas_per_vec;
  
  std::vector<_I> r2N_tot_cnt_vec;
  
  std::vector<_I> dBc_pas_cnt_vec;
  std::vector<_d> dBp_pas_per_vec;
  
  std::vector<_d> r2_vec;
  
  std::vector< std::vector<_d> > dB_stats_vec(5);
  
  for ( _I row_idx=0; row_idx<mat_nrow; row_idx++ ) {
    _S cg( (_S) rownames_vec_r[row_idx] );
    
    _I prb_tot_cnt = 0; _I cmb_tot_cnt = 0;
    _I prb_pas_cnt = 0; _I cmb_pas_cnt = 0;
    _d prb_pas_per = 0; _d cmb_pas_per = 0;
    
    _I dBc_pas_cnt = 0; 
    _d dBp_pas_per = 0;
    _d r2 = 0;
    _I r2N_tot_cnt = 0;
    
    std::vector<_d> datA_vec;
    std::vector<_d> datB_vec;
    
    for ( _I colA_idx=0; colA_idx<mat_ncol; colA_idx++ ) {
      prb_tot_cnt++;
      if ( !Rcpp::NumericVector::is_na(  beta_mat_r(row_idx,colA_idx) )  ) 
        prb_pas_cnt++;
      
      for ( _I colB_idx=colA_idx+1; colB_idx<mat_ncol; colB_idx++ ) {
        cmb_tot_cnt++;
        // Skip NA/NaN Values::
        if ( Rcpp::NumericVector::is_na(  beta_mat_r(row_idx,colA_idx) ) || 
             Rcpp::NumericVector::is_na(  beta_mat_r(row_idx,colB_idx) )  ) continue;
        cmb_pas_cnt++;
        
        _d valA = (_d) beta_mat_r(row_idx,colA_idx);
        _d valB = (_d) beta_mat_r(row_idx,colB_idx);
        
        datA_vec.push_back( valA );
        datB_vec.push_back( valB );
      }
    } // Done Loading Data::
    r2N_tot_cnt = cmb_tot_cnt;
    
    if ( prb_tot_cnt != 0 ) prb_pas_per = (_d) prb_pas_cnt / (_d) prb_tot_cnt;
    
    if ( cmb_tot_cnt != 0 ) cmb_pas_per = (_d) cmb_pas_cnt / (_d) cmb_tot_cnt;
    if ( prb_tot_cnt != 0 && cmb_pas_cnt != 0 ) {
      
      if ( cmb_pas_cnt != 0 ) dBc_pas_cnt = pass_count( datA_vec,datB_vec,min_dB,-1,true );
      if ( cmb_pas_cnt != 0 ) dBp_pas_per = (_d) dBc_pas_cnt / (_d) cmb_pas_cnt;
      
      std::vector<_d> dB_stats_cur = get_dB_stats5_vec( datA_vec, datB_vec, true );
      dB_stats_vec.resize( dB_stats_cur.size() );
      for (size_t ii=0; ii<dB_stats_cur.size(); ii++ )
        dB_stats_vec[ii].push_back( dB_stats_cur[ii] );
      
      // Adding Corners::
      _d r0 = pearsoncoeff_vec( datA_vec, datB_vec );
      datA_vec.push_back(0.0);
      datB_vec.push_back(0.0);
      datA_vec.push_back(1.0);
      datB_vec.push_back(1.0);
      _d r1 = pearsoncoeff_vec( datA_vec, datB_vec );
      
      // Add additional corner cases porportional to sample size::
      for ( size_t jj=0; jj<int(prb_tot_cnt/4) - 1; jj++ ) {
        datA_vec.push_back(0.0);
        datB_vec.push_back(0.0);
        datA_vec.push_back(1.0);
        datB_vec.push_back(1.0);
      }
      r2 = pearsoncoeff_vec( datA_vec, datB_vec );
      
      if ( false && row_idx == 3 ) {
        std::cerr<<"R0='"<<r0<<"', "
                 <<"R1='"<<r1<<"', "
                 <<"R2='"<<r2<<"', "
                 <<"cmb_size='"<<datA_vec.size()<<"'"<<std::endl;
        for ( size_t jj=0; jj<datA_vec.size(); jj++ ) {
          std::cerr<<"\t row_idx["<<jj<<"]='"<<datA_vec[jj]<<"', '"<<datB_vec[jj]<<"'"<<std::endl;
        }
      }
      
      r2N_tot_cnt = datA_vec.size();
    } else {
      for (size_t ii=0; ii<5; ii++ )
        dB_stats_vec[ii].push_back( -1 );
    }
    
    cg_vec.push_back( cg );
    
    prb_tot_cnt_vec.push_back( prb_tot_cnt );
    prb_pas_cnt_vec.push_back( prb_pas_cnt );
    prb_pas_per_vec.push_back( prb_pas_per );
    
    cmb_tot_cnt_vec.push_back( cmb_tot_cnt );
    cmb_pas_cnt_vec.push_back( cmb_pas_cnt );
    cmb_pas_per_vec.push_back( cmb_pas_per );
    
    r2N_tot_cnt_vec.push_back( r2N_tot_cnt );

    dBc_pas_cnt_vec.push_back( dBc_pas_cnt );
    dBp_pas_per_vec.push_back( dBp_pas_per );
    
    r2_vec.push_back( r2 );
    
  } // End row ( probe evaluation )
  
  df.push_back( cg_vec,  "Probe_ID" );
  
  if ( true ) {
    df.push_back( prb_tot_cnt_vec, "rN" );
    df.push_back( prb_pas_cnt_vec, "pN" );
    // df.push_back( cmb_tot_cnt_vec, "cN" );
    df.push_back( r2N_tot_cnt_vec, "RN" );
    
    df.push_back( prb_pas_per_vec, "pD" );
    df.push_back( r2_vec,  "r2" );
    df.push_back( dBp_pas_per_vec, "dB" );
    df.push_back( dB_stats_vec[4], "vB" );
    df.push_back( dB_stats_vec[3], "aB" );
    df.push_back( dB_stats_vec[1], "mB" );
  } else {
    df.push_back( prb_tot_cnt_vec, "Probe_Total_Count" );
    df.push_back( prb_pas_cnt_vec, "Probe_Pass_Count" );
    df.push_back( prb_pas_per_vec, "Probe_Pass_Perc" );
    
    df.push_back( cmb_tot_cnt_vec, "Combn_Total_Count" );
    df.push_back( cmb_pas_cnt_vec, "Combn_Pass_Count" );
    df.push_back( cmb_pas_per_vec, "Combn_Pass_Perc" );
    
    df.push_back( dBc_pas_cnt_vec, "dB_Pass_Count" );
    df.push_back( dBp_pas_per_vec, "dB_Pass_Percent" );
    df.push_back( r2_vec,  "r2" );
    
    std::vector<_S> db_stat_Names;
    db_stat_Names.push_back("dB_Q1");
    db_stat_Names.push_back("dB_Q2");
    db_stat_Names.push_back("dB_Q3");
    db_stat_Names.push_back("dB_Avg");
    db_stat_Names.push_back("dB_Sds");
    
    for ( size_t ii=0; ii<dB_stats_vec.size(); ii++ )
      df.push_back( dB_stats_vec[ii], db_stat_Names[ii] );
  }
  
  /* ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                                  Done!
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

  return( df );
};

/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *                        Sample Performance Function::
 * 
 * sample_performance_mat()
 *   - Description:: Compare All-Against-All samples r2/dB and ppp (percent 
 *      passing pval) 
 *      
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

// [[Rcpp::export]]
Rcpp::DataFrame sample_performance_mat( Rcpp::NumericMatrix beta_mat_r,
                                        Rcpp::NumericMatrix pval_mat_r,
                                        const _d min_dB = 0.2,
                                        
                                        const unsigned int vb = 0,
                                        const unsigned int vt = 1,
                                        const unsigned int tc = 0,
                                        const std::string ft="sample_performance_mat" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0 = vb > vt + 0;
  
  Rcpp::DataFrame df;
  
  _I mat_nrow = beta_mat_r.nrow();
  _I mat_ncol = beta_mat_r.ncol();
  
  CharacterVector rownames_vec_r = rownames( beta_mat_r );
  CharacterVector colnames_vec_r = colnames( beta_mat_r );
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<<"Statring...\n"
    <<_fo<<"\t mat_nrow: '" << mat_nrow<<"'\n"
    <<_fo<<"\t mat_ncol: '" << mat_ncol<<"'\n"
    <<std::endl;
  
  /*
   * Delta Beta Calculation:: 
   *     std::map< _S, std::map< _S,_d > > dB_vec(mat_ncol);
   *
   */
  std::vector<_S> sentrixA_vec;
  std::vector<_S> sentrixB_vec;
  std::vector<_I> beta_len_vec;
  std::vector<_I> pval_len_vec;
  std::vector<_d> pval_pas_vec;
  std::vector<_d> min_dB_vec;
  
  std::vector<_I> dB_count_vec;
  std::vector<_d> dB_pperc_vec;
  std::vector<_d> r2_vec;
  
  for ( _I colA_idx=0; colA_idx<mat_ncol; colA_idx++ ) {
    for ( _I colB_idx=colA_idx; colB_idx<mat_ncol; colB_idx++ ) {
      _S nameA_str( (_S) colnames_vec_r[colA_idx] );
      _S nameB_str( (_S) colnames_vec_r[colB_idx] );
      
      Rcpp::NumericVector colA_vec_r = beta_mat_r( _ , colA_idx );
      Rcpp::NumericVector colB_vec_r = beta_mat_r( _ , colB_idx );
      
      size_t colR_len = std::min( colA_vec_r.size(),colB_vec_r.size() );
      
      std::vector<_d> colA_vec;
      rcpp_to_cpp_vec_dbl( colA_vec_r,colA_vec );
      
      std::vector<_d> colB_vec;
      rcpp_to_cpp_vec_dbl( colB_vec_r,colB_vec );
      
      _I dB_pass_cnt = pass_count( colA_vec,colB_vec,min_dB, -1,true );
      _d dB_pass_per = (_d) dB_pass_cnt / (_d) mat_nrow;
      
      /*
       * Calculate r2 (pearson r-squared correlation)::
       *   - NOTE:: Need to remove NA values...
       */
      colA_vec.clear();
      colB_vec.clear();
      colA_vec.resize( colR_len );
      colB_vec.resize( colR_len );
      
      size_t last_idx = 0;
      for ( size_t ii=0; ii <colR_len; ii++ ) {
        if ( Rcpp::NumericVector::is_na( colA_vec_r(ii) ) ||
             Rcpp::NumericVector::is_na( colB_vec_r(ii) ) ) continue;
        colA_vec[ii] = (_d) colA_vec_r(ii);
        colB_vec[ii] = (_d) colB_vec_r(ii);
        last_idx++;
      }
      colA_vec.resize( last_idx );
      colB_vec.resize( last_idx );
      // TBD:: This doesn't match the R version. Investigate later...
      //   Probably something to do with the pval failure removals. We should
      //   just use the pval matrix...
      // Calculate r2::
      _d r2 = pearsoncoeff_vec( colA_vec,colB_vec );
      
      // Calculate p-value passing percent::
      _d ppp_val = (_d) last_idx / (_d) colR_len;
      
      // Push Back All Data::
      sentrixA_vec.push_back( nameA_str );
      sentrixB_vec.push_back( nameB_str );
      min_dB_vec.push_back( min_dB );
      beta_len_vec.push_back( colR_len );
      pval_len_vec.push_back( last_idx );
      pval_pas_vec.push_back( ppp_val );
      dB_count_vec.push_back( dB_pass_cnt );
      dB_pperc_vec.push_back( dB_pass_per );
      r2_vec.push_back( r2 );
      
      // Push Back recirical non-identical comparisons::
      if ( nameA_str.compare(nameB_str) != 0 ) {
        sentrixA_vec.push_back( nameB_str );
        sentrixB_vec.push_back( nameA_str );
        min_dB_vec.push_back( min_dB );
        beta_len_vec.push_back( colR_len );
        pval_len_vec.push_back( last_idx );
        pval_pas_vec.push_back( ppp_val );
        dB_count_vec.push_back( dB_pass_cnt );
        dB_pperc_vec.push_back( dB_pass_per );
        r2_vec.push_back( r2 );
      }
    }
  }
  
  // Update DataFrame::
  df.push_back( sentrixA_vec, "Sentrix_Name_A" );
  df.push_back( sentrixB_vec, "Sentrix_Name_B" );
  df.push_back( min_dB_vec,   "Min_Delta_Beta" );
  
  df.push_back( beta_len_vec, "Pre_Pval_Count" );
  df.push_back( pval_len_vec, "Pos_Pval_Count" );
  
  df.push_back( dB_count_vec, "dB_Pass_Count" );
  df.push_back( pval_pas_vec, "Pval_Pass_Percent" );
  df.push_back( dB_pperc_vec, "dB_Pass_Percent" );
  df.push_back( r2_vec, "r2" );
  
  // for ( size_t ii=0; ii<dB_stats_vec.size(); ii++ )
  //   df.push_back( dB_stats_vec[ii], db_stat_Names[ii] );
  
  return( df );
};

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * External Rcpp Functions::
 * 
 *                              Load Pair of Idats::
 *                              
 *                                   Goals::
 *                                   
 *  1. Idat::
 *     - detectionPnegEcdf()
 *        Calculate Negative Control Back Grounds
 *        
 * Internal::
 *     read_pair_idats_rcpp()
 *     
 *     set_detection()
 *     set_manifest()
 *      
 *      - build manifest
 *      - write "Auto Sample Sheet"
 *      - write additional outputs
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// [[Rcpp::export]]
Rcpp::DataFrame read_pair_idats_rcpp( const std::string& prefix, 
                                      std::vector< unsigned int >& pval_add_vec,
                                      
                                      const std::string outdir = "",
                                      const double min_pval = 0.05,
                                      
                                      const bool read_bgz = false,
                                      const bool write_bgz = false,
                                      
                                      const unsigned int vb = 0,
                                      const unsigned int vt = 1,
                                      const unsigned int tc = 0,
                                      const std::string ft="read_pair_idats_rcpp" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0 = vb > vt + 0;
  // const bool p1 = vb > vt + 1;
  // const bool p2 = vb > vt + 2;
  const bool p8 = vb > vt + 8;
  
  Rcpp::DataFrame df;
  bool success = true;
  
  std::string grn_idat_path( prefix+"_Grn.idat.gz" );
  std::string red_idat_path( prefix+"_Red.idat.gz" );
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<< "Statring...\n"
    <<_fo<< "\t    Prefix: '" << prefix << "'\n"
    <<_fo<< "\t Input Grn: '" << grn_idat_path << "'\n"
    <<_fo<< "\t Input Red: '" << red_idat_path << "'\n"
    <<std::endl;
  
  if ( true ) {
    idat::Idat grn_idat( grn_idat_path, 'G', 
                         pval_add_vec,
                         vb,vt+1,tc+1 );
    if ( p0 ) Rcpp::Rcerr
      <<_fo<<"Grn Idat:: to_string():\n"
      <<_fo<<"\n"<<grn_idat.to_string()<<"\n"
      <<std::endl;
    
  }
  if ( false ) {
    idat::Idat red_idat( red_idat_path, 'R', 
                         vb,vt+1,tc+1 );
    if ( p8 ) Rcpp::Rcerr
      <<_fo<<"Red Idat:: to_string():\n"
      <<_fo<<"\n"<<red_idat.to_string()<<"\n"
      <<std::endl;
  }
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   *                                    DONE::
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  if ( p0 ) Rcpp::Rcerr 
    <<_fo<< "Done. Return: '" << success << "'\n" << std::endl;
  
  return( df );
}

// [[Rcpp::export]]
void test_math() {
  _S _fo="test_math: ";
  std::cerr<<"\n"<<_fo<<BRK<<std::endl;
  
  std::vector<_i> i1_vec({ 0,1,2,3,4,5,6,7,8,9 });
  std::vector<_i> i2_vec({ 2,2,2,2,2,2,2,2,2,2 });
  
  std::vector<_i> ti_dciv = i1_vec - 10;
  std::cerr<<_fo<<" delta(ti,10) = "<<std::endl;
  for ( auto& x : ti_dciv ) std::cerr<<_fo<<"\t x = '"<<x<<"'"<<std::endl;
  
  std::vector<_i> t1_dt2v = i1_vec - i2_vec;
  std::cerr<<_fo<<" delta(t1,t2) = "<<std::endl;
  for ( auto& x : t1_dt2v ) std::cerr<<_fo<<"\t x = '"<<x<<"'"<<std::endl;
  
  std::vector<_i> t1_pt2v = i1_vec * i2_vec;
  std::cerr<<_fo<<" prod(t1,t2) = "<<std::endl;
  for ( auto& x : t1_pt2v ) std::cerr<<_fo<<"\t x = '"<<x<<"'"<<std::endl;
  
  std::cerr<<"\n"<<_fo<<BRK<<std::endl;
  
  std::vector<_d> t1_vec({ 0,1,2,3,4,5,6,7,8,9 });
  std::vector<_d> t2_vec = t1_vec * -1.5;
  for ( size_t ii=0; ii<t2_vec.size(); ii+=2 ) t2_vec[ii]+=0.3;
  
  _i t_abs = simple_abs( t1_vec[0],t2_vec[0] );
  std::cerr<<_fo<<" simple_abs() = '"<<t_abs<<"'"<<std::endl;
  _i t_sum = sum_vec( t1_vec );
  std::cerr<<_fo<<" sum_vec() = '"<<t_sum<<"'"<<std::endl;
  _d t_avg = mean_vec( t1_vec );
  std::cerr<<_fo<<" mean_vec() = '"<<t_avg<<"'"<<std::endl;
  _d t_sqs = sqsum_vec( t1_vec );
  std::cerr<<_fo<<" sqsum_vec() = '"<<t_sqs<<"'"<<std::endl;
  _d t_std = stdev_vec( t1_vec );
  std::cerr<<_fo<<" stdev() = '"<<t_std<<"'"<<std::endl;
  std::cerr<<std::endl;
  
  std::vector<_d> t_q3v = quantile3_vec( t1_vec, false );
  std::cerr<<_fo<<" quantile3_vec() = "<<std::endl;
  for ( auto& x : t_q3v ) std::cerr<<_fo<<"\t x = '"<<x<<"'"<<std::endl;
  
  std::vector<_d> t_s5v = stats5_vec( t1_vec, false );
  std::cerr<<_fo<<" stats5_vec() = "<<std::endl;
  for ( auto& x : t_s5v ) std::cerr<<_fo<<"\t x = '"<<x<<"'"<<std::endl;
  std::cerr<<std::endl;
  
  _d r2 = pearsoncoeff_vec( t1_vec,t2_vec );
  std::cerr<<_fo<<" r2 = '"<<set_precision(r2,3)<<"'"<<std::endl;
  
  _I dB = pass_count( t1_vec,t2_vec, 4,-1,true );
  std::cerr<<_fo<<" dB = '"<<dB<<"'"<<std::endl;
  
  std::cerr<<_fo<<"Done."<<std::endl;
  std::cerr<<_fo<<BRK<<std::endl;
  
}


/*** R
# timesTwo(42)

run_math <- TRUE
run_math <- FALSE
if ( run_math ) test_math();

run_test <- TRUE
run_test <- FALSE

if ( run_test ) {
  suppressWarnings(suppressPackageStartupMessages( 
    base::require("tidyverse",  quietly = TRUE) ) )
  
  loc_opt <- NULL
  loc_opt$top_path <- "/Users/bbarnes/Documents"
  neg_ctl_rds  <- file.path( loc_opt$top_path, "data/manifests/methylation/bgz/all_negative_ctls.rds" )
  neg_ctl_tib  <- readr::read_rds( neg_ctl_rds )
  
  #
  # All out of date after new computer::
  #
  # loc_opt$prefix <- file.path( loc_opt$top_path, "Projects/Embark/Embark_Methylation_Internal_V3/Experiment_ConcentrationTitration/206712840006/206712840006_R01C01" )
  # loc_opt$prefix <- file.path( loc_opt$top_path, "data/idats/idats_NA12878/200348350023/200348350023_R02C01" )
  # loc_opt$prefix <- file.path( loc_opt$top_path, "data/idats/idats_Evonik/206662930005/206662930005_R06C02" )
  # loc_opt$prefix <- file.path( loc_opt$top_path, "data/idats/idats_EX_20210127_Meth_Ref_Material_amd_1um_Comp/204981580007/204981580007_R01C01" )

  loc_opt$prefix <- file.path( loc_opt$top_path, "data/idats/EPIC/202296710014/202296710014_R01C01" )
  
  loc_opt$output_dir <- "/Users/bbarnes/Documents/tmp/as_idat_bgz"
  if ( !dir.exists( loc_opt$output_dir ) ) base::dir.create( loc_opt$output_dir, recursive = TRUE )
  
  vt <- 1
  tc <- 0
  
  vb <- 4 # One level deep (p2)
  vb <- 3 # Standard
  vb <- 2 # Light
  
  success <- FALSE

  loc_opt$min_dB   <- 0.20
  loc_opt$min_pval <- 0.05
  
  loc_opt$run_idat <- FALSE
  loc_opt$run_idat <- TRUE
  
  loc_opt$write_bgz <- FALSE
  loc_opt$write_bgz <- TRUE
  
  idat_df <- NULL
  if ( loc_opt$run_idat )
    idat_df <- read_pair_idats_rcpp( prefix = loc_opt$prefix,
                                     pval_add_vec = neg_ctl_tib$Address %>% as.vector(), 
                                     outdir    = loc_opt$output_dir, 
                                     min_pval  = loc_opt$min_pval,
                                     read_bgz  = FALSE, 
                                     write_bgz = loc_opt$write_bgz, 
                                     vb=vb,vt=vt,tc=tc )
  if ( !is.null(idat_df) ) idat_df %>% print()
  
  # cat( glue::glue("\nDONE(test_read_idats): success='{success}'\n\n") )
}

*/

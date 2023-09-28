
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppCommon.h>

// My headers::
#include "bisulfite.h"

using namespace Rcpp;
using namespace bisulfite;

/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *                    Rcpp to Cpp Local Conversion Function::
 * 
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

void rcpp_to_cpp_vec_int( Rcpp::IntegerVector& vec_r,
                          std::vector<int>& vec ) {
  _I s = vec_r.size();
  
  vec.clear();
  vec.resize( s );
  for ( _I ii = 0; ii < s; ii++ ) vec[ii] = vec_r(ii);
}

void rcpp_to_cpp_vec_uint( Rcpp::IntegerVector& vec_r,
                           std::vector< _I>& vec ) {
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
 *                      Template Rcpp to c++ Function::
 * 
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

// [[Rcpp::export]]
Rcpp::DataFrame templace_rcpp_func( Rcpp::StringVector  aln_ids_vec_r,
                                    Nullable<Rcpp::StringVector>  aln_seq_vec_r_ = R_NilValue,
                                    Nullable<Rcpp::StringVector>  aln_tag_vec_r_ = R_NilValue,
                                    Nullable<Rcpp::StringVector>  aln_chr_vec_r_ = R_NilValue,
                                    Nullable<Rcpp::IntegerVector> aln_pos_vec_r_ = R_NilValue,
                                    
                                    const bool uc = false,
                                    const unsigned int vb = 0,
                                    const unsigned int vt = 1,
                                    const unsigned int tc = 0,
                                    const std::string ft="templace_rcpp_func" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p1 = vb > vt + 0;
  
  bool success = true;
  
  _I aln_ids_vec_size = aln_ids_vec_r.size();
  
  if ( p1 ) Rcpp::Rcerr
  <<_fo<< "Starting...\n"
  <<_fo<< "\t aln_ids_vec_size = "<<aln_ids_vec_size<<"\n"
  <<std::endl;
  
  Rcpp::DataFrame df = DataFrame::create();
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                 Allocate c++ vectors from Rcpp vectors::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  std::vector<std::string> aln_ids_vec( aln_ids_vec_size );
  std::vector<std::string> aln_seq_vec( aln_ids_vec_size );
  std::vector<std::string> aln_tag_vec( aln_ids_vec_size );
  
  // Copy Over Vectors::
  rcpp_to_cpp_vec_str( aln_ids_vec_r, aln_ids_vec );
  
  if ( aln_seq_vec_r_.isNotNull() ) {
    if ( p1 ) Rcpp::Rcerr <<_fo<< "aln_seq_vec_r is NOT NULL\n" <<std::endl;
    Rcpp::StringVector aln_seq_vec_r(aln_seq_vec_r_);
    rcpp_to_cpp_vec_str( aln_seq_vec_r, aln_seq_vec );
  } else {
    if ( p1 ) Rcpp::Rcerr <<_fo<< "aln_seq_vec_r is NULL\n" <<std::endl;
  }
  
  if ( p1 ) Rcpp::Rcerr
  <<_fo<< "Status...\n"
  <<_fo<< "\t aln_ids_vec_size = "<<aln_ids_vec.size()<<"\n"
  <<std::endl;
  
  if ( !success ) {
    Rcpp::Rcerr
    <<_fe<< "Failed during processing!\n"
    <<std::endl;
    return( df );
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *
   *                 Convert Final Output to Return Rcpp DataFrame::
   * 
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  // Direct Version::
  df = Rcpp::DataFrame::create(
    Rcpp::Named("Probe_ID") = aln_ids_vec,
    Rcpp::Named("Aln_Seq") = aln_seq_vec,
    Rcpp::Named("Aln_Tag") = aln_tag_vec
  );
  
  // Loop Version::
  // for ( _I ii = 0; ii < name_vec.size(); ii++ ) {
  //   df.push_back( int_vecs[ii], name_vec[ii] )
  //   
  // }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                                   Done::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  if ( p1 ) Rcpp::Rcerr
  <<_fo<< "Done!\n"
  <<std::endl;
  
  return( df );
}

/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *                          Load 122mer Function::
 * 
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

// [[Rcpp::export]]
Rcpp::DataFrame seqs_to_2bit( std::vector< _S >& seq_vec,
                              Nullable<Rcpp::IntegerVector> cgn_vec_r_ = R_NilValue,

                              const bool uc = false,
                              const unsigned int vb = 0,
                              const unsigned int vt = 1,
                              const unsigned int tc = 0,
                              const std::string ft="seqs_to_2bit" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p1 = vb > vt + 0;
  
  bool success = true;
  
  _I seq_vec_size = seq_vec.size();
  
  if ( p1 ) Rcpp::Rcerr
  <<_fo<< "Starting...\n"
  <<_fo<< "\t seq_vec_size = "<<seq_vec_size<<"\n"
  <<std::endl;
  
  Rcpp::DataFrame df = DataFrame::create();
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                 Allocate c++ vectors from Rcpp vectors::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  // std::vector<std::string> aln_ids_vec( aln_ids_vec_size );
  // std::vector<std::string> aln_seq_vec( aln_ids_vec_size );
  // std::vector<std::string> aln_tag_vec( aln_ids_vec_size );
  
  // Copy Over Vectors::
  // std::vector< _i > cgn_vec;
  // rcpp_to_cpp_vec_int( cgn_vec_r_, cgn_vec );
  // _I cgn_vec_size = cgn_vec.size();
  // 
  // if ( p1 ) Rcpp::Rcerr
  // <<_fo<< "Status...\n"
  // <<_fo<< "\t cgn_vec_size = "<<cgn_vec_size<<"\n"
  // <<std::endl;
  
  if ( !success ) {
    Rcpp::Rcerr
    <<_fe<< "Failed during processing!\n"
    <<std::endl;
    return( df );
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *
   *                         Loop through sequences::
   * 
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  for ( size_t ii=0; ii<seq_vec_size; ii++ ) {
    _S seq_str = seq_vec[ii];
    uint8_t* seq_bit = dna_to_bit( seq_str.c_str(), seq_str.size() );
    // _I seq_uint = (_I) seq_bit;
    
    Rcpp::Rcerr
    <<_fo<< "Seq='"<<seq_str<<"', uint='"<<seq_bit<<"'\n"
    <<std::endl;
    
    
    if ( ii > 3 ) break;
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *
   *                 Convert Final Output to Return Rcpp DataFrame::
   * 
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  // Direct Version::
  // df = Rcpp::DataFrame::create(
  //   Rcpp::Named("Probe_ID") = aln_ids_vec,
  //   Rcpp::Named("Aln_Seq") = aln_seq_vec,
  //   Rcpp::Named("Aln_Tag") = aln_tag_vec
  // );
  
  // Loop Version::
  // for ( _I ii = 0; ii < name_vec.size(); ii++ ) {
  //   df.push_back( int_vecs[ii], name_vec[ii] )
  //   
  // }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                                   Done::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  if ( p1 ) Rcpp::Rcerr
  <<_fo<< "Done!\n"
  <<std::endl;
  
  return( df );
}

/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *                 BSMAP Quick Parsing for Large Alignments::
 * 
 Nullable<Rcpp::StringVector>  aln_seq_vec_r_ = R_NilValue,
 Nullable<Rcpp::StringVector>  aln_tag_vec_r_ = R_NilValue,
 Nullable<Rcpp::StringVector>  aln_chr_vec_r_ = R_NilValue,
 Nullable<Rcpp::IntegerVector> aln_pos_vec_r_ = R_NilValue,
 Nullable<Rcpp::StringVector>  aln_srd_vec_r_ = R_NilValue,
 Nullable<Rcpp::IntegerVector> ins_cnt_vec_r_ = R_NilValue,
 Nullable<Rcpp::StringVector>  ref_seq_vec_r_ = R_NilValue,
 Nullable<Rcpp::IntegerVector> mis_cnt_vec_r_ = R_NilValue,
 Nullable<Rcpp::StringVector>  mis_str_vec_r_ = R_NilValue,
 * 
 * 
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

// [[Rcpp::export]]
Rcpp::DataFrame parse_bsmap( Rcpp::StringVector  aln_ids_vec_r,
                             Rcpp::StringVector  aln_tag_vec_r,
                             Rcpp::StringVector  aln_chr_vec_r,
                             Rcpp::IntegerVector aln_pos_vec_r,
                             Rcpp::StringVector  aln_srd_vec_r,
                             Rcpp::IntegerVector mis_cnt_vec_r,
                             
                             Rcpp::StringVector  fai_chr_vec_r,
                             Rcpp::IntegerVector fai_len_vec_r,
                             
                             const unsigned int max_cnt = 0,
                             const bool uc = false,
                             const unsigned int vb = 0,
                             const unsigned int vt = 1,
                             const unsigned int tc = 0,
                             const std::string ft="parse_bsmap" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0  = vb > vt + 0;
  const bool p1  = vb > vt + 1;
  const bool p10 = vb > vt + 10;
  
  bool success = true;
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<< "Starting...\n"
    <<std::endl;
  
  _I aln_ids_vec_size = aln_ids_vec_r.size();
  _I fai_chr_vec_size = fai_chr_vec_r.size();
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Input Parameters::\n"
    <<_fo<< "\t aln_ids_vec_size = "<<aln_ids_vec_size<<"\n"
    <<_fo<< "\t fai_chr_vec_size = "<<fai_chr_vec_size<<"\n"
    <<std::endl;
  
  Rcpp::DataFrame df = DataFrame::create();
  
  // std::vector<std::string> aln_ids_vec( aln_ids_vec_size );
  // rcpp_to_cpp_vec_str( aln_ids_vec_r, aln_ids_vec );
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                     Declare/Allocate Data Structures::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  // Option 1:: hit_vec[mud][ cgn_tbs_cos_inf ].push( struct )
  // std::vector<std::map<std::string, std::vector<std::vector<std::string>> > > hit_vec;
  
  // Option 2:: 
  //  str_map[ col_name ].vec[ridx] = col_value
  //  int_map[ col_name ].vec[ridx] = col_value
  //  top_map[ cgn_tbs_cos_inf_mud ] = ridx
  //  top_vec[ ridx ] = cgn_tbs_cos_inf_mud
  //
  // Init col_names_vec
  // Init each of the columns to the size of aln_ids_vec_size
  
  std::vector<std::string> col_str_names({
    "Probe_ID", "Loci_ID", 
    "Strand_TB", "Strand_CO", "Allele",
    "Bsp_Tag", "Bsp_Chr", "Bsp_Srd" });
  
  std::vector<std::string> col_int_names({
    "Infinium_Desigin", "Bsp_Beg_Pos","Bsp_CpG_Pos", "Bsp_Mis_Cnt",
    "Bsp_Chr_Len",
    "Bsp_Hit_Cnt", "Chr_Len_Sum", "Chr_Len_Min", "Chr_Len_Max",
    "Top_Bsp_Hit" });
  
  std::vector<std::string> col_dbl_names({
    "Chr_Len_Beta" });
  
  std::map<std::string, std::vector<std::string>> str_dat_map;
  std::map<std::string, std::vector<unsigned int>> int_dat_map;
  std::map<std::string, std::vector<double>> dbl_dat_map;
  
  for ( _I ii=0; ii<col_str_names.size(); ii++ ) {
    str_dat_map[col_str_names[ii]] = std::vector<std::string>(aln_ids_vec_size, "");
  }
  for ( _I ii=0; ii<col_int_names.size(); ii++ ) {
    int_dat_map[col_int_names[ii]] = std::vector<unsigned int>(aln_ids_vec_size, 0);
  }
  for ( _I ii=0; ii<col_dbl_names.size(); ii++ ) {
    dbl_dat_map[col_dbl_names[ii]] = std::vector<double>(aln_ids_vec_size, 0.0);
  }
  
  // TBD:: Initialze Summary map<Probe_ID, unsigned in>
  //   - int_sum_map = { "Bsp_Hit_Cnt", "Chr_Len_Sum", "Chr_Len_Min", "Chr_Len_Max", "Top_Bsp_Hit" }
  //   - top_sum_map = { "Top_Cgn_Idx", "Mis_Cnt_Min", "Chr_Len_Max" }
  std::map<std::string, std::map<std::string, unsigned int>> int_sum_map;
  std::map<std::string, std::map<std::string, unsigned int>> top_sum_map;
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Initialized Output Vectors::\n"
    <<_fo<< "\t str_vecs = "<<str_dat_map.size()<<"\n"
    <<_fo<< "\t int_vecs = "<<int_dat_map.size()<<"\n"
    <<std::endl;
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                     Build Chromosome Length Mappings::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  _I min_chr_len = 0;
  std::map<std::string, unsigned int> fai_map;
  for ( _I cidx = 0; cidx < fai_chr_vec_size; cidx++ ) {
    _S chr( fai_chr_vec_r(cidx) );
    _I len( fai_len_vec_r(cidx) );
    
    if ( len == 0 ) {
      Rcpp::Rcerr <<_fe<< "Current Chromosome["<<cidx<<"] ='"<<chr<<"' "
                  << " has length("<<len<<") eqal to zero!\n"
                  <<std::endl;
      return( df );
    }
    //
    // TESTING:: Dropping the last digit from Chromosome Length...
    //
    len = int( len/10 );
    
    if ( min_chr_len==0 || min_chr_len > len ) min_chr_len = len;
    
    fai_map[chr] = len;
  }
  _I fai_len = fai_map.size();
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "FAI Chromosome Stats::\n"
    <<_fo<< "\t     fai_len = "<<fai_len<<"\n"
    <<_fo<< "\t min_chr_len = "<<min_chr_len<<"\n"
    <<std::endl;
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *
   *  Workflow::
   *    - Parse each 'line'
   * 
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  _S cgn = "";
  _B tbs = '0';
  _B cos = '0';
  _I inf_int =  0;
  _B ale = '0';
  _I mud =  3;
  
  for ( _I ridx = 0; ridx < aln_ids_vec_size; ridx++ ) {
    mud = 3;
    
    /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
     *                Convert all the input data from Rcpp to c++::
     * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
    
    _S probe_id_str( aln_ids_vec_r(ridx) );
    _S tag_str( aln_tag_vec_r(ridx) );
    _S chr_str( aln_chr_vec_r(ridx) );
    _I beg_pos( aln_pos_vec_r(ridx) );
    _S srd_str( aln_srd_vec_r(ridx) );
    _I mis_cnt( mis_cnt_vec_r(ridx) );
    _S prb_type( probe_id_str.substr( 0, 2) );
    
    /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
     *                            Get Chromosome Length::
     * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
    
    if ( fai_map.find(chr_str) == fai_map.end() ) {
      Rcpp::Rcerr <<_fe<< "Parsing probe_id["<<ridx<<"] ='"<<probe_id_str<<"' "
                  << " failed to find chrom("<<chr_str<<") in fai_map!\n"
                  <<std::endl;
      return( df );
    }
    _I chr_len = fai_map[chr_str] - min_chr_len;
    
    /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
     *                              Parse Probe_ID::
     * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
    
    success = parse_probe_id( probe_id_str, cgn, tbs, cos, inf_int, ale );
    if ( !success || tbs=='0' || cos=='0' || inf_int==0 || ale=='0' ) {
      Rcpp::Rcerr <<_fe<< "Parsing probe_id["<<ridx<<"] ='"<<probe_id_str<<"' "
                  << "failed!\n"
                  <<std::endl;
      return( df );
    }
    if ( ale == 'A' && inf_int == 2 ) mud = 2;
    if ( ale == 'A' && inf_int == 1 ) mud = 0;
    if ( ale == 'B' && inf_int == 1 ) mud = 1;
    if ( mud > 2 ) {
      Rcpp::Rcerr <<_fe<< "Parsing probe_id["<<ridx<<"] ='"<<probe_id_str<<"' "
                  << " mud("<<mud<<") > 2!\n"
                  <<std::endl;
      return( df );
    }
    
    /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
     *                   Calculate CpG Postion from BSP Offset::
     * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
    
    int off_set = get_bsmap_offset( prb_type, inf_int, srd_str );
    if ( off_set < -2 || off_set > 50 ) {
      Rcpp::Rcerr <<_fe<< "Parsing probe_id["<<ridx<<"] ='"<<probe_id_str<<"' "
                  << " off_set("<<off_set<<") < -2 OR > 50!\n"
                  <<std::endl;
      return( df );
    }
    _I cpg_pos = beg_pos + off_set;
    
    if ( p10 ) Rcpp::Rcerr
      <<_fo<< "Status["<<ridx<<"].2::\n"
      <<_fo<< "\t probe_id_str    = '"<<probe_id_str<<"'\n"
      <<_fo<< "\t         cgn_str = '"<<cgn<<"'\n"
      <<_fo<< "\t       strand_TB = '"<<tbs<<"'\n"
      <<_fo<< "\t       strand_CO = '"<<cos<<"'\n"
      <<_fo<< "\t        infinium = '"<<inf_int<<"'\n"
      <<_fo<< "\t          allele = '"<<ale<<"'\n\n"
      <<_fo<< "\t         tag_str = '"<<tag_str<<"'\n"
      <<_fo<< "\t         chr_str = '"<<chr_str<<"'\n"
      <<_fo<< "\t         beg_pos = '"<<beg_pos<<"'\n"
      <<_fo<< "\t         srd_str = '"<<srd_str<<"'\n"
      <<_fo<< "\t         mis_cnt = '"<<mis_cnt<<"'\n\n"
      <<_fo<< "\t         off_set = '"<<off_set<<"'\n"
      <<_fo<< "\t         beg_pos = '"<<beg_pos<<"'\n"
      <<_fo<< "\t         cpg_pos = '"<<cpg_pos<<"'\n"
      <<_fo<< "\t         chr_str = '"<<chr_str<<"'\n"
      <<_fo<< "\t         chr_len = '"<<chr_len<<"'\n"
      <<std::endl;
    
    /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
     *                           Update Summary Maps::
     * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
    
    if ( int_sum_map.find(probe_id_str) == int_sum_map.end() ) {
      int_sum_map[probe_id_str]["Bsp_Hit_Cnt"] = 0;
      int_sum_map[probe_id_str]["Chr_Len_Sum"] = 0;
      int_sum_map[probe_id_str]["Chr_Len_Min"] = chr_len;
      int_sum_map[probe_id_str]["Chr_Len_Max"] = chr_len;
    }
    int_sum_map[probe_id_str]["Bsp_Hit_Cnt"]++;
    int_sum_map[probe_id_str]["Chr_Len_Sum"] += chr_len;
    if ( chr_len < int_sum_map[probe_id_str]["Chr_Len_Min"] )
      int_sum_map[probe_id_str]["Chr_Len_Min"] = chr_len;
    if ( chr_len > int_sum_map[probe_id_str]["Chr_Len_Max"] )
      int_sum_map[probe_id_str]["Chr_Len_Max"] = chr_len;
    
    /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
     *                           Update Top Hits Maps::
     * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
    
    //   - top_sum_map = { "Top_Cgn_Idx", "Mis_Cnt_Min", "Chr_Len_Max" }
    if ( top_sum_map.find(probe_id_str) == top_sum_map.end() ) {
      top_sum_map[probe_id_str]["Top_Cgn_Idx"] = ridx;
      top_sum_map[probe_id_str]["Mis_Cnt_Min"] = mis_cnt;
      top_sum_map[probe_id_str]["Chr_Len_Max"] = chr_len;
    }
    if ( mis_cnt <= top_sum_map[probe_id_str]["Mis_Cnt_Min"] &&
         chr_len >= top_sum_map[probe_id_str]["Chr_Len_Max"] ) {
      top_sum_map[probe_id_str]["Top_Cgn_Idx"] = ridx;
      top_sum_map[probe_id_str]["Mis_Cnt_Min"] = mis_cnt;
      top_sum_map[probe_id_str]["Chr_Len_Max"] = chr_len;
    }
    
    /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
     *                               Save Hit Data::
     * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
    
    // TBD:: Load these vectors now...
    //
    // str_dat_map =
    // std::vector<std::string> col_str_names({
    //   "Probe_ID", "Loci_ID", 
    //   "Strand_TB", "Strand_CO", "Allele",
    //   "Bsp_Tag", "Bsp_Chr", "Bsp_Srd" });
    std::string tbs_str( 1, static_cast<char>(tbs) );
    std::string cos_str( 1, static_cast<char>(cos) );
    std::string ale_str( 1, static_cast<char>(ale) );
    
    str_dat_map["Probe_ID"][ridx]  = probe_id_str;
    str_dat_map["Loci_ID"][ridx]   = cgn;
    str_dat_map["Strand_TB"][ridx] = tbs_str;
    str_dat_map["Strand_CO"][ridx] = cos_str;
    str_dat_map["Allele"][ridx]    = ale_str;
    str_dat_map["Bsp_Tag"][ridx]   = tag_str;
    str_dat_map["Bsp_Chr"][ridx]   = chr_str;
    str_dat_map["Bsp_Srd"][ridx]   = srd_str;
    
    // int_dat_map =
    // std::vector<std::string> col_int_names({
    //   "Infinium_Desigin", "Bsp_Beg_Pos","Bsp_CpG_Pos", "Bsp_Mis_Cnt",
    //   "Bsp_Chr_Len",
    //   "Bsp_Hit_Cnt", "Chr_Len_Sum", "Chr_Len_Min", "Chr_Len_Max", 
    //   "Top_Bsp_Hit" });
    int_dat_map["Infinium_Desigin"][ridx] = inf_int;
    int_dat_map["Bsp_Beg_Pos"][ridx]      = beg_pos;
    int_dat_map["Bsp_CpG_Pos"][ridx]      = cpg_pos;
    int_dat_map["Bsp_Mis_Cnt"][ridx]      = mis_cnt;
    int_dat_map["Bsp_Chr_Len"][ridx]      = chr_len;
    
    // NOTE:: Bsp_Hit_Cnt/Len_Sum must be calculated outside of this loop
    // TBD:: switch the order so Bsp_Hit_Cnt and Len_Sum are last...
    
    if ( max_cnt != 0 && ridx >= max_cnt ) break;
  }
  
  if ( !success ) {
    Rcpp::Rcerr
    <<_fe<< "Failed during processing!\n"
    <<std::endl;
    return( df );
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                        Calculate Summary Stats::
   *                         { Bsp_Hit_Cnt,Len_Sum }
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  for ( _I ridx=0; ridx<aln_ids_vec_size; ridx++ ) {
    _S prb_id = str_dat_map["Probe_ID"][ridx];
    
    int_dat_map["Bsp_Hit_Cnt"][ridx] = int_sum_map[prb_id]["Bsp_Hit_Cnt"];
    int_dat_map["Chr_Len_Sum"][ridx] = int_sum_map[prb_id]["Chr_Len_Sum"];
    int_dat_map["Chr_Len_Min"][ridx] = int_sum_map[prb_id]["Chr_Len_Min"];
    int_dat_map["Chr_Len_Max"][ridx] = int_sum_map[prb_id]["Chr_Len_Max"];
    
    if ( top_sum_map[prb_id]["Top_Cgn_Idx"] == ridx )
      int_dat_map["Top_Bsp_Hit"][ridx] = 1;
    
    if ( int_sum_map[prb_id]["Chr_Len_Sum"] != 0 )
      dbl_dat_map["Chr_Len_Beta"][ridx] = (double)int_dat_map["Bsp_Chr_Len"][ridx] / (double)int_sum_map[prb_id]["Chr_Len_Sum"];
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                 Convert Final Output to Return Rcpp DataFrame::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  // Direct Version::
  // df = Rcpp::DataFrame::create(
  //   Rcpp::Named("Probe_ID") = aln_ids_vec,
  //   Rcpp::Named("Aln_Seq") = aln_seq_vec,
  // );
  
  // Loop Version::
  for ( _I ii = 0; ii < col_str_names.size(); ii++ ) {
    df.push_back( str_dat_map[ col_str_names[ii] ], col_str_names[ii] );
  }
  for ( _I ii = 0; ii < col_int_names.size(); ii++ ) {
    df.push_back( int_dat_map[ col_int_names[ii] ], col_int_names[ii] );
  }
  for ( _I ii = 0; ii < col_dbl_names.size(); ii++ ) {
    df.push_back( dbl_dat_map[ col_dbl_names[ii] ], col_dbl_names[ii] );
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                                   Done::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<< "Done!\n"
    <<std::endl;
  
  return( df );
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * Exported Rcpp Functions::
 * 
 *   replace_nuc_at_cpp()
 *   replace_nucs_at_cpp()
 *   mutate_seqs_cpp()
 *   improbe_seqs_cpp()
 *   
 *   [Questionable]: load_tag_tsv()
 *   [Questionable]: load_ucsc_snps_bed()
 *   
 *   fwd2tops_cpp()
 *   expand_seqs_cpp()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// [[Rcpp::export]]
Rcpp::StringVector replace_nuc_at_cpp( Rcpp::StringVector  seq_vec_r,
                                       Rcpp::IntegerVector pos_vec_r,
                                       Rcpp::StringVector  rep_vec_r,
                                       Rcpp::StringVector  org_vec_r,
                                       const bool uc = false,
                                       unsigned int vb = 0,
                                       unsigned int vt = 1 ) {
  
  const std::string func_tag("[replace_nucs_at_cpp]: ");
  const std::string func_err("[replace_nucs_at_cpp]: ERROR: ");
  const std::string func_wrn("[replace_nucs_at_cpp]: Warning: ");
  
  if ( vb >= vt ) Rcpp::Rcerr << func_tag << "Starting..." << "\n";
  
  // Return Data Frame for R::
  Rcpp::StringVector ret_nil;
  
  // Sanity Check:: Quick input validation::
  unsigned int seq_vec_size = seq_vec_r.size();
  unsigned int pos_vec_size = pos_vec_r.size();
  unsigned int rep_vec_size = rep_vec_r.size();
  unsigned int org_vec_size = org_vec_r.size();
  
  if ( seq_vec_size == 0 || 
       seq_vec_size != pos_vec_size || 
       rep_vec_size != seq_vec_size )
  {
    Rcerr << func_err << "Input vector sizes do not match or are zero!" << "\n"
          << func_err << "\t seq_vec_size = " << seq_vec_size << "\n"
          << func_err << "\t pos_vec_size = " << pos_vec_size << "\n"
          << func_err << "\t rep_vec_size = " << rep_vec_size << "\n"
          << func_err << "\t org_vec_size = " << org_vec_size << "\n"
          << std::endl;
    return( ret_nil );
  }
  
  // Allow for zero length inputs::
  LogicalVector nan_vec = is_na( seq_vec_r );
  
  unsigned int idx = 0;
  std::string rep  = "";
  std::string org  = "";
  std::vector<std::string> seq_vec( seq_vec_size );
  for (unsigned int seq_idx = 0; seq_idx < seq_vec_size; seq_idx++) {
    if ( !nan_vec[seq_idx] )
    {
      seq_vec[seq_idx] = seq_vec_r(seq_idx);
      idx = pos_vec_r(seq_idx);
      rep  = rep_vec_r(seq_idx);
      org  = "";
      if ( org_vec_size == seq_vec_size ) org = org_vec_r(seq_idx);
      
      if ( idx > 0 ) replace_nuc_at( seq_vec[seq_idx], idx-1, rep, org, uc );
    }
  }
  
  // Convert back to Rcpp StringVector()
  Rcpp::StringVector seqs_ret( seq_vec.size() );
  seqs_ret = seq_vec;
  
  return seqs_ret;
}

// [[Rcpp::export]]
Rcpp::DataFrame replace_nucs_at_cpp( Rcpp::StringVector  cgn_vec_r,
                                     Rcpp::IntegerVector cnt_vec_r,
                                     Rcpp::StringVector  seq_vec_r,
                                     
                                     Rcpp::IntegerVector pos_vec_r,
                                     Rcpp::StringVector  alt_vec_r,
                                     Rcpp::StringVector  ref_vec_r,
                                     
                                     const bool uc = false,
                                     unsigned int vb = 0,
                                     unsigned int vt = 1 ) {
  
  const std::string func_tag("[replace_nucs_at_cpp]: ");
  const std::string func_err("[replace_nucs_at_cpp]: ERROR: ");
  const std::string func_wrn("[replace_nucs_at_cpp]: Warning: ");
  
  if ( vb >= vt ) Rcpp::Rcerr << func_tag << "Starting..." << "\n";
  
  // Return Data Frame for R::
  Rcpp::DataFrame ret_df;
  // Rcpp::StringVector ret_nil;
  
  // Sanity Check:: Quick input validation::
  unsigned int cgn_vec_size = cgn_vec_r.size();
  unsigned int cnt_vec_size = cnt_vec_r.size();
  unsigned int seq_vec_size = seq_vec_r.size();
  
  unsigned int pos_vec_size = pos_vec_r.size();
  unsigned int alt_vec_size = alt_vec_r.size();
  unsigned int ref_vec_size = ref_vec_r.size();
  
  if ( cgn_vec_size == 0 || 
       cgn_vec_size != cnt_vec_size ||
       cgn_vec_size != seq_vec_size ||
       pos_vec_size == 0 ||
       pos_vec_size != alt_vec_size )
  {
    Rcerr << func_err << "Input vector sizes do not match or are zero!" << "\n"
          << func_err << "\t cgn_vec_size = " << cgn_vec_size << "\n"
          << func_err << "\t cnt_vec_size = " << cnt_vec_size << "\n"
          << func_err << "\t seq_vec_size = " << seq_vec_size << "\n"
          << func_err << "\t pos_vec_size = " << pos_vec_size << "\n"
          << func_err << "\t alt_vec_size = " << alt_vec_size << "\n"
          << func_err << "\t ref_vec_size = " << ref_vec_size << "\n"
          << std::endl;
    return( ret_df );
  }
  
  // Allow for zero length inputs::
  LogicalVector nan_vec = is_na( seq_vec_r );
  
  // unsigned int var_cnt = 0;
  unsigned int vec_idx = 0;
  unsigned int var_pos = 0;
  std::string var_alt  = "";
  std::string var_ref  = "";
  
  std::vector<std::string> cgn_vec( cgn_vec_size );
  std::vector<unsigned int> cnt_vec( cnt_vec_size );
  std::vector<std::string> seq_vec( seq_vec_size );
  std::vector<std::string> raw_vec( seq_vec_size );
  
  for (unsigned int seq_idx = 0; seq_idx < seq_vec_size; seq_idx++) {
    if ( !nan_vec[seq_idx] )
    {
      // var_cnt = cnt_vec_r(seq_idx);
      cgn_vec[seq_idx] = cgn_vec_r(seq_idx);
      cnt_vec[seq_idx] = cnt_vec_r(seq_idx);
      seq_vec[seq_idx] = seq_vec_r(seq_idx);
      raw_vec[seq_idx] = seq_vec_r(seq_idx);
      
      if ( vb >= vt+1 )  Rcpp::Rcerr 
        << func_tag << "Current seq_idx = " << seq_idx << "!\n"
        << func_tag << "\t          seq_idx = " << seq_idx << "\n"
        << func_tag << "\t cgn_vec[seq_idx] = " << cgn_vec[seq_idx] << "\n"
        << func_tag << "\t cnt_vec[seq_idx] = " << cnt_vec[seq_idx] << "\n"
        << func_tag << "\t seq_vec[seq_idx] = " << seq_vec[seq_idx] << "\n"
        << std::endl;
      
      for (unsigned int var_idx = 0; var_idx < cnt_vec[seq_idx]; var_idx++)
      {
        
        if ( vec_idx >=  alt_vec_size )
        {
          Rcerr << func_err << "Vector idx >= alt_vec_size!\n"
                << func_err << "\t      vec_idx = " << vec_idx << "\n"
                << func_err << "\t pos_vec_size = " << pos_vec_size << "\n"
                << func_err << "\t alt_vec_size = " << alt_vec_size << "\n"
                << func_err << "\t ref_vec_size = " << ref_vec_size << "\n"
                << std::endl;
          return( ret_df );
        }
        
        var_pos = pos_vec_r(vec_idx);
        var_alt = alt_vec_r(vec_idx);
        var_ref = "";
        if ( alt_vec_size == ref_vec_size ) var_ref = ref_vec_r(vec_idx);
        
        if ( var_pos > 0 ) 
          replace_nuc_at( seq_vec[seq_idx], var_pos-1, var_alt, var_ref, uc );
        
        vec_idx++;
      }
      if ( cnt_vec[seq_idx] == 0 ) vec_idx++;
    }
  }
  
  // Convert back to Rcpp StringVector()
  // Rcpp::StringVector seqs_ret( seq_vec.size() );
  // seqs_ret = seq_vec;
  // return seqs_ret;
  
  ret_df = Rcpp::DataFrame::create(
    Rcpp::Named("Probe_ID") = cgn_vec,
    Rcpp::Named("SNP_Count") = cnt_vec,
    Rcpp::Named("SNP_Sequence") = seq_vec,
    Rcpp::Named("RAW_Sequence") = raw_vec
  );
  
  return( ret_df );
}

// [[Rcpp::export]]
Rcpp::StringVector mutate_seqs_cpp( Rcpp::StringVector seqs, 
                                    std::string m,
                                    const bool uc = false,
                                    unsigned int vb = 0,
                                    unsigned int vt = 1 ) {
  
  if ( !init_rvcp ) revcomp_init();
  
  const std::string func_tag("[mutate_seqs_cpp]: ");
  const std::string func_err("[mutate_seqs_cpp]: ERROR: ");
  const std::string func_wrn("[mutate_seqs_cpp]: Warning: ");
  
  if ( vb >= vt ) Rcpp::Rcerr << func_tag << "Starting..." << "\n";
  
  // Allow for zero length inputs::
  LogicalVector na_vec = is_na( seqs );
  
  std::vector<std::string> seq_vec( seqs.size() );
  for (int seq_idx = 0; seq_idx < seqs.size(); seq_idx++) {
    if ( !na_vec[seq_idx] )
    {
      seq_vec[seq_idx] = seqs(seq_idx);
      mutate_seq( seq_vec[seq_idx], m, uc );
    }
  }
  
  // Convert back to Rcpp StringVector()
  Rcpp::StringVector seqs_ret( seq_vec.size() );
  seqs_ret = seq_vec;
  
  return seqs_ret;
}

// [[Rcpp::export]]
Rcpp::DataFrame improbe_seqs_cpp( Rcpp::StringVector ids_vec_r,
                                  Rcpp::StringVector fwd_vec_r,
                                  Rcpp::StringVector din_vec_r,
                                  std::string bsc_str_r = "numd",
                                  std::string frs_str_r = "FR",
                                  std::string cos_str_r = "CO",
                                  Nullable<Rcpp::StringVector> frs_vec_r_ = R_NilValue,
                                  Nullable<Rcpp::StringVector> cos_vec_r_ = R_NilValue,
                                  
                                  const unsigned int prb_len = 50,
                                  const bool return_source = false,
                                  const bool uc = false,
                                  const unsigned int vb = 0, 
                                  const unsigned int vt = 1 )
{
  const std::string func_tag("[improbe_seqs_cpp]: ");
  const std::string func_err("[improbe_seqs_cpp]: ERROR: ");
  const std::string func_wrn("[improbe_seqs_cpp]: Warning: ");
  
  const bool p1 = vb >= vt + 1;
  // const bool p2 = vb >= vt + 2;
  
  if ( !init_rvcp ) revcomp_init();
  if ( !init_topbot ) topbot_init();
  
  if ( p1 ) Rcerr << func_tag << "Starting..." << "\n";
  
  // Return Data Frame for R::
  Rcpp::DataFrame ret_df;
  
  // Sanity Check:: Quick input validation::
  unsigned int ids_vec_size = ids_vec_r.size();
  unsigned int fwd_vec_size = fwd_vec_r.size();
  unsigned int din_vec_size = din_vec_r.size();
  
  if ( ids_vec_size == 0 || 
       ids_vec_size != fwd_vec_size ||
       fwd_vec_size != din_vec_size )
  {
    Rcpp::Rcerr 
    << func_err << "Input vector sizes do not match or are zero!" << "\n"
    << func_err << "\t ids_vec_size = " << ids_vec_size << "\n"
    << func_err << "\t fwd_vec_size = " << fwd_vec_size << "\n"
    << func_err << "\t din_vec_size = " << din_vec_size << "\n"
    << std::endl;
    return( ret_df );
  }
  
  if ( p1 ) Rcpp::Rcerr 
    << func_tag << "Input vector sizes match and not length zero!" << "\n"
    << func_tag << "\t ids_vec_size = " << ids_vec_size << "\n"
    << func_tag << "\t fwd_vec_size = " << fwd_vec_size << "\n"
    << func_tag << "\t din_vec_size = " << din_vec_size << "\n"
    << std::endl;
  
  // Allocate Conversion from R to c++ sequence vector::
  std::vector<std::string> ids_vec( ids_vec_size );
  std::vector<std::string> fwd_vec( fwd_vec_size );
  std::vector<std::string> din_vec( din_vec_size );
  
  Rcpp::StringVector frs_vec_r;
  Rcpp::StringVector cos_vec_r;
  
  std::vector< std::vector<std::string> > frs_vec;
  std::vector< std::vector<std::string> > cos_vec;
  std::vector<std::string> frd_vec;
  std::vector<std::string> cod_vec;
  
  unsigned int frs_vec_size = 0;
  unsigned int cos_vec_size = 0;
  unsigned int frd_vec_size = 0;
  unsigned int cod_vec_size = 0;
  unsigned int srd_cnt = 1;
  
  std::vector<char> FR_valid_chars = { 'F', 'R' };
  std::vector<char> CO_valid_chars = { 'C', 'O' };
  
  /*
   * Parse FR Strand Inputs::
   *   - Either Target Design State String is set or Vector is Provided
   */
  if ( frs_vec_r_.isNotNull() )
  {
    frs_vec_r = frs_vec_r_;
    frs_vec_size = frs_vec_r.size();
    if ( fwd_vec_size != frs_vec_size )
    {
      Rcpp::Rcerr << func_err << "FR Strand vector incorrect size!" << "\n"
                  << func_err << "\t fwd_vec_size = " << fwd_vec_size << "\n"
                  << func_err << "\t frs_vec_size = " << frs_vec_size << "\n"
                  << std::endl;
      return( ret_df );
    }
    frs_vec.resize( frs_vec_size );
    frs_vec_size = frs_vec.size();
    
  } else {
    
    if ( !chop_valid_str( frs_str_r, frd_vec, FR_valid_chars, vb, vt+3 ) )
      return( ret_df );
    frd_vec_size = frd_vec.size();
    
    srd_cnt *= frd_vec_size;
    if ( p1 ) Rcerr 
      << func_tag << "Status: srd_cnt = " << srd_cnt << "\n"
      << func_tag << "\t frd_vec_size = " << frd_vec_size << "\n"
      << std::endl;
  }
  
  /*
   * Parse CO Strand Inputs::
   *   - Either Target Design State String is set or Vector is Provided
   */
  if ( cos_vec_r_.isNotNull() )
  {
    cos_vec_r = cos_vec_r_;
    cos_vec_size = cos_vec_r.size();
    if ( fwd_vec_size != cos_vec_size )
    {
      Rcpp::Rcerr << func_err << "CO Strand vector incorrect size!" << "\n"
                  << func_err << "\t fwd_vec_size = " << fwd_vec_size << "\n"
                  << func_err << "\t cos_vec_size = " << cos_vec_size << "\n"
                  << std::endl;
      return( ret_df );
    }
    cos_vec.resize( cos_vec_size );
    cos_vec_size = cos_vec.size();
    
  } else {
    
    if ( !chop_valid_str( cos_str_r, cod_vec, CO_valid_chars, vb, vt+3 ) )
      return( ret_df );
    cod_vec_size = cod_vec.size();
    
    srd_cnt *= cod_vec_size;
    if ( p1 ) Rcerr 
      << func_tag << "Status: srd_cnt = " << srd_cnt << "\n"
      << func_tag << "\t cod_vec_size = " << cod_vec_size << "\n"
      << std::endl;
  }
  
  if ( frs_vec_size != 0 || cos_vec_size != 0 )
  {
    if ( frs_vec_size != 0 ) frd_vec.resize( frs_vec_size );
    if ( cos_vec_size != 0 ) cod_vec.resize( cos_vec_size );
    
    for ( unsigned int seq_idx = 0; seq_idx < fwd_vec_size; seq_idx++ )
    {
      if ( frs_vec_size != 0 ) {
        frd_vec[seq_idx] = frs_vec_r( seq_idx );
        ::str_to_upper( frd_vec[seq_idx] );
        if ( !chop_valid_str( frd_vec[seq_idx], frs_vec[seq_idx], FR_valid_chars, vb, vt+3 ) )
          return( ret_df );
      }
      if ( cos_vec_size != 0 ) {
        cod_vec[seq_idx] = cos_vec_r( seq_idx );
        ::str_to_upper( cod_vec[seq_idx] );
        if ( !chop_valid_str( cod_vec[seq_idx], cos_vec[seq_idx], CO_valid_chars, vb, vt+3 ) )
          return( ret_df );
      }
    }
    if ( frs_vec_size != 0 ) {
      std::vector<char> frc_unq_vec = vec_to_unique_chars(frd_vec);
      frd_vec.clear();
      srd_cnt *= frc_unq_vec.size();
      if ( p1 ) Rcerr 
        << func_tag << "Status: srd_cnt = " << srd_cnt << "\n"
        << func_tag << "\t frc_unq_vec = " << frc_unq_vec.size() << "\n"
        << std::endl;
    }
    if ( cos_vec_size != 0 ) {
      std::vector<char> coc_unq_vec = vec_to_unique_chars(cod_vec);
      cod_vec.clear();
      srd_cnt *= coc_unq_vec.size();
      if ( p1 ) Rcerr 
        << func_tag << "Status: srd_cnt = " << srd_cnt << "\n"
        << func_tag << "\t coc_unq_vec = " << coc_unq_vec.size() << "\n"
        << std::endl;
    }
  }
  
  /*
   * Build output header names::
   * 
   */
  const std::string bsc_str = bsc_str_r;
  const unsigned int bsc_len = bsc_str.size();
  std::vector<std::string> bsc_vec;
  for ( std::size_t ii = 0; ii < bsc_len; ii++ )
  {
    std::string cur_bsc = str_to_str_idx( bsc_str, false, true, ii );
    bsc_vec.push_back( cur_bsc );
  }
  unsigned int bsc_vec_size = bsc_vec.size();
  
  std::vector<std::string> out_names_vec;
  // User Provided Fields::
  out_names_vec.push_back("Probe_ID");   // 0
  out_names_vec.push_back("Probe_Type"); // 1
  out_names_vec.push_back("Strand_SR");  // 2
  out_names_vec.push_back("Strand_CO");  // 3
  out_names_vec.push_back("Strand_BC");  // 4
  out_names_vec.push_back("Template_Sequence_FWD"); // 5
  
  std::vector<std::string> bsc_names_vec;
  bsc_names_vec.push_back("Template_Sequence"); // 6:0
  bsc_names_vec.push_back("NXB"); // 7:1
  bsc_names_vec.push_back("CPN"); // 8:2
  bsc_names_vec.push_back("SEC"); // 9:3
  bsc_names_vec.push_back("BOD"); // 10:3
  bsc_names_vec.push_back("END"); // 11:5
  bsc_names_vec.push_back("PRB"); // 12:6
  
  unsigned int bsc_names_size = bsc_names_vec.size();
  for ( std::size_t ii = 0; ii < bsc_vec_size; ii++ )
    for ( std::size_t jj = 0; jj < bsc_names_size; jj++ )
      out_names_vec.push_back( bsc_names_vec[jj]+"_"+bsc_vec[ii] );
  unsigned int out_names_size = out_names_vec.size();
  
  // std::map< std::string, std::string > bsc_names_map;
  // for ( std::size_t ii = 0; ii < bsc_vec_size; ii++ )
  //   for ( std::size_t jj = 0; jj < bsc_names_size; jj++ )
  //     bsc_names_map.insert( pair<std::string, std::string>( (std::string) bsc_names_vec[jj], (std::string) bsc_vec[ii] ) );
  
  /*
   * Allocate Space::
   *   - bsc_vecs[ fields ][ seq_idx ] = Probe Designs
   *   - tri_maps< tri-seq >< n46-seq >< umd_str > = vector( seq_idx )
   */
  unsigned int long_des_cnt = fwd_vec_size * srd_cnt; // bsc_vec_size;
  std::vector< std::vector< std::string > > 
    bsc_vecs( out_names_size, std::vector<std::string>(long_des_cnt, "") );
  
  std::map< std::string, std::map< std::string, std::map< std::string, std::vector<unsigned int> > > > tri_maps;
  
  if ( p1 ) Rcpp::Rcerr 
    << func_tag << "Vector: out_names_size = " << out_names_size << std::endl
    << func_tag << "Vector:   long_des_cnt = " << long_des_cnt << "\n\n"
    << func_tag << "      Testing: bsc_vecs.size() = '" << bsc_vecs.size() << "'\n"
    << func_tag << "   Testing: bsc_vecs[0].size() = '" << bsc_vecs[0].size() << "'\n"
    << func_tag << "Testing: bsc_vecs[0][0].size() = '" << bsc_vecs[0][0].size() << "'\n"
    << func_tag << "       Testing: bsc_vecs[0][0] = '" << bsc_vecs[0][0] << "'\n"
    << std::endl;
  
  std::string ids_str;
  std::string fwd_seq;
  std::string frs_str;
  std::string cos_str;
  std::string din_str;
  
  unsigned int cur_frs_size;
  unsigned int cur_cos_size;
  std::vector<std::string> cur_frs_vec;
  std::vector<std::string> cur_cos_vec;
  
  bool success = false;
  
  unsigned int out_idx = 0;
  for ( unsigned int seq_idx = 0; seq_idx < fwd_vec_size; seq_idx++ )
  {
    if ( p1 ) Rcpp::Rcerr 
      << func_tag << "Status: seq_idx = " << seq_idx << std::endl;
    
    // R to c++ conversion:: FWD & REV
    ids_vec[seq_idx] = ids_vec_r( seq_idx );
    din_vec[seq_idx] = din_vec_r( seq_idx );
    fwd_vec[seq_idx] = fwd_vec_r( seq_idx );
    
    ids_str = ids_vec[seq_idx];
    din_str = din_vec[seq_idx];
    fwd_seq = fwd_vec[seq_idx];
    
    // Get FR Vector::
    if ( frs_vec_size != 0 ) cur_frs_vec = frs_vec[seq_idx];
    else cur_frs_vec = frd_vec;
    cur_frs_size = cur_frs_vec.size();
    
    // Get CO Vector::
    if ( cos_vec_size != 0 ) cur_cos_vec = cos_vec[seq_idx];
    else cur_cos_vec = cod_vec;
    cur_cos_size = cur_cos_vec.size();
    
    for ( std::size_t frs_idx = 0; frs_idx < cur_frs_size; frs_idx++ )
    {
      frs_str = cur_frs_vec[frs_idx];
      if ( p1 ) Rcpp::Rcerr 
        << func_tag << "Status:\t out_idx = " << out_idx << std::endl
        << func_tag << "Status:\t frs_str = " << frs_str << std::endl
        << std::endl;
      
      for ( std::size_t cos_idx = 0; cos_idx < cur_cos_size; cos_idx++ )
      {
        cos_str = cur_cos_vec[cos_idx];
        frs_str = cur_frs_vec[frs_idx];
        if ( p1 ) Rcpp::Rcerr 
          << func_tag << "Status:\t\t out_idx = " << out_idx << std::endl
          << func_tag << "Status:\t\t cos_str = " << cos_str << std::endl
          << func_tag << "Status:\t\t fwd_len = " << fwd_vec[seq_idx].size() << std::endl
          << std::endl;
        
        if ( fwd_vec[seq_idx].size() > 0 )
        {
          success = false;
          success = improbe_seq( bsc_vecs,
                                 ids_str, fwd_seq, out_idx, 
                                 out_names_vec, tri_maps,
                                 frs_str, cos_str, bsc_str, din_str,
                                 prb_len, vb+1, vt+1, 2 );
          
          if ( p1 ) Rcpp::Rcerr
            << func_tag << "Status:\t\t bsc_vecs = " << bsc_vecs[0][out_idx] << std::endl;
          
          if ( !success )
          {
            Rcpp::Rcerr << func_err << "improbe_seq() failure!!!" << "\n"
                        << func_err << "\t frs_str = " << frs_str << "\n"
                        << func_err << "\t cos_str = " << cos_str << "\n"
                        << func_err << "\t bsc_str = " << bsc_str << "\n"
                        << func_err << "\t din_str = " << din_str << "\n"
                        << func_err << "\t fwd_seq = " << fwd_seq << "\n"
                        << std::endl;
            return( ret_df );
          }
          
          out_idx++;
        }
        
      } // END:: for ( std::size_t cos_idx = 0; cos_idx < cos_vec_size; cos_idx++ )
    } // END:: for ( std::size_t frs_idx = 0; frs_idx < cur_frs_size; frs_idx++ )
  } // END:: for ( unsigned int seq_idx = 0; seq_idx < fwd_vec_size; seq_idx++ )
  
  // Trim any extra allocated space::
  if ( p1 ) Rcpp::Rcerr 
    << func_tag << "Done. Last Index = " << out_idx << "\n" << std::endl;
  for ( std::size_t col_idx = 0; col_idx < out_names_size; col_idx++ )
  {
    bsc_vecs[col_idx].resize(out_idx);
  }
  
  // ret_df = Rcpp::DataFrame::create();
  for ( std::size_t col_idx = 0; col_idx < out_names_size; col_idx++ )
  {
    std::string col_key = out_names_vec[col_idx];
    ret_df[col_key] = bsc_vecs[col_idx];
  }
  
  return( ret_df );
}

// TBD:: Add Column Types::
//   - Rcpp::List col_types,

// [[Rcpp::export]]
Rcpp::DataFrame load_tag_tsv( std::string file,
                              Rcpp::StringVector col_names_vec_r,
                              const unsigned int vb = 0,
                              const unsigned int vt = 1  )
{
  const std::string func_tag("[load_tag_tsv]: ");
  const std::string func_err("[load_tag_tsv]: ERROR: ");
  const std::string func_wrn("[load_tag_tsv]: Warning: ");
  
  const bool p1 = vb >= vt + 1;
  
  if ( p1 ) Rcerr << func_tag << "Loading: '" << file << "'\n";
  
  // Environment pkg = Environment::namespace_env("tibble");
  Environment pkg = Environment::namespace_env("readr");
  Function read_tsv_r("read_tsv");
  Function cols_r("cols");
  Function col_character_r("col_character");
  Function col_integer_r("col_integer");
  
  // Return Data Frame for R::
  // Rcpp::DataFrame ret_df = read_tsv_r( Named("file")=file );
  // col_names col_types
  
  // std::vector<std::string> col_names_vec;
  // col_names_vec.push_back("tag");
  // col_names_vec.push_back("end");
  // col_names_vec.push_back("cgn");
  // col_names_vec.push_back("man");
  // col_names_vec.push_back("mud_map");
  // col_names_vec.push_back("info");
  // unsigned int col_names_vec_size = col_names_vec.size();
  
  unsigned int col_names_vec_size = col_names_vec_r.size();
  std::vector<std::string> col_names_vec( col_names_vec_size );
  
  for ( std::size_t ii = 0; ii < col_names_vec_size; ii++ )
  {
    col_names_vec[ii] = col_names_vec_r( ii );
    Rcerr << func_tag << "col_names["<<ii<<"] = '"<< col_names_vec[ii] << "'\n"
          << std::endl;
  }
  
  // Rcpp::DataFrame ret_df;
  Rcpp::DataFrame ret_df = read_tsv_r( Named("file")=file,
                                       Named("col_names")=col_names_vec );
  
  Rcerr << func_tag << "ref_df[0]: '" << ret_df.nrows() << "'\n\n\n"
        << std::endl;
  
  return( ret_df );
}  

// [[Rcpp::export]]
Rcpp::DataFrame load_ucsc_snps_bed( std::string file,
                                    Rcpp::StringVector col_names_vec_r,
                                    const char sep,
                                    
                                    const bool uc = false,
                                    const unsigned int vb = 0,
                                    const unsigned int vt = 1  )
{
  const std::string func_tag("[load_ucsc_snps_bed]: ");
  const std::string func_err("[load_ucsc_snps_bed]: ERROR: ");
  const std::string func_wrn("[load_ucsc_snps_bed]: Warning: ");
  
  const bool p1 = vb >= vt + 1;
  
  if ( p1 ) Rcerr << func_tag << "Loading: '" << file << "'\n";
  
  // Environment pkg = Environment::namespace_env("tibble");
  Environment pkg = Environment::namespace_env("readr");
  Function read_tsv_r("read_tsv");
  
  // Return Data Frame for R::
  Rcpp::DataFrame df;
  
  unsigned int col_names_vec_size = col_names_vec_r.size();
  std::vector<std::string> col_names_vec( col_names_vec_size );
  
  for ( std::size_t ii = 0; ii < col_names_vec_size; ii++ )
  {
    col_names_vec[ii] = col_names_vec_r( ii );
    Rcerr << func_tag << "col_names["<<ii<<"] = '"<< col_names_vec[ii]
          << std::endl;
  }
  
  // Rcpp::DataFrame df;
  df = read_tsv_r( Named("file")=file,
                   Named("col_names")=col_names_vec );
  
  unsigned int df_nrows = df.nrows();
  unsigned int df_ncols = df.size();
  
  if ( p1 ) Rcpp::Rcerr 
    << func_tag << "Data Stats: \n"
    << func_tag << "\t nrows = '" << df_nrows << "'\n"
    << func_tag << "\t ncols = '" << df_ncols << "'\n"
    << std::endl;
  
  // std::vector<std::string>  chr_vec( df_ncols );
  // std::vector<unsigned int> beg_vec( df_ncols );
  // std::vector<unsigned int> end_vec( df_ncols );
  // std::vector<std::string>  cgn_vec( df_ncols );
  
  // std::vector<std::string>  alt_cnt_vec( df_ncols );
  std::vector<std::string>  ref_var_vec( df_ncols );
  std::vector<std::string>  alt_str_vec( df_ncols );
  
  // std::vector<std::string>  indel_shift_vec( df_ncols );
  // std::vector<std::string>  evidence_cnt_vec( df_ncols );
  
  std::vector<std::string>  MAF_str_vec( df_ncols );
  std::vector<std::string>  MAX_str_vec( df_ncols );
  std::vector<std::string>  MIN_str_vec( df_ncols );
  
  // std::vector<std::string>  var_type_vec( df_ncols );
  // std::vector<std::string>  comments_vec( df_ncols );
  
  // Rcpp::StringVector alt_cnt_r = df["alt_cnt"];
  Rcpp::StringVector ref_r  = df["ref"];
  Rcpp::StringVector alts_r = df["alts"];
  Rcpp::StringVector MAF_r  = df["mafs"];
  Rcpp::StringVector MAX_r  = df["maxs"];
  Rcpp::StringVector MIN_r  = df["mins"];
  
  std::vector< std::vector<std::string> > alt_vecs( df_nrows );
  std::vector< std::vector<std::string> > maf_vecs( df_nrows );
  std::vector< std::vector<std::string> > max_vecs( df_nrows );
  std::vector< std::vector<std::string> > min_vecs( df_nrows );
  
  // std::map<std::string, int> alt_map;
  // std::map<std::string, std::vector<int> > maf_map;
  // std::map<std::string, int> max_map;
  // std::map<std::string, int> min_map;
  
  df_ncols = 100;
  for ( std::size_t var_idx = 0; var_idx < df_nrows; var_idx++ )
  {
    // alt_map.clear();
    // maf_map.clear();
    // max_map.clear();
    // min_map.clear();
    
    std::vector<std::string> alt_vec;
    std::vector<std::string> maf_vec;
    std::vector<std::string> max_vec;
    std::vector<std::string> min_vec;
    
    ref_var_vec[var_idx] = ref_r( var_idx );
    alt_str_vec[var_idx] = alts_r( var_idx );
    MAF_str_vec[var_idx] = MAF_r( var_idx );
    MAX_str_vec[var_idx] = MAX_r( var_idx );
    MIN_str_vec[var_idx] = MIN_r( var_idx );
    
    // if ( p1 ) Rcpp::Rcerr
    //   << func_tag << "MAF["<<var_idx<<"]: '"<< MAF_str_vec[var_idx] << "'"
    //   << std::endl;
    
    parse_line( alt_str_vec[var_idx], ',', alt_vec );
    parse_line( MAF_str_vec[var_idx], ',', maf_vec );
    parse_line( MAX_str_vec[var_idx], ',', max_vec );
    parse_line( MIN_str_vec[var_idx], ',', min_vec );
    
    bool found_ref = false;
    bool found_alt = false;
    
    unsigned int maf_size = maf_vec.size();
    unsigned int alt_size = alt_vec.size();
    std::vector<bool> found_alts( alt_size, false );
    
    for ( std::size_t maf_idx = 0; maf_idx < maf_size; maf_idx++ )
    {
      
      if ( maf_vec[maf_idx].size() != 0 && maf_vec[maf_idx].compare("-inf") != 0 )
      {
        if ( ref_var_vec[var_idx].compare(max_vec[maf_idx]) == 0 ) found_ref = true;
        if ( ref_var_vec[var_idx].compare(min_vec[maf_idx]) == 0 ) found_ref = true;
        
        for ( std::size_t alt_idx = 0; alt_idx < alt_size; alt_idx++ )
        {
          if ( alt_vec[alt_idx].compare(max_vec[maf_idx]) == 0 ||
               alt_vec[alt_idx].compare(min_vec[maf_idx]) == 0 )
          {
            found_alt = true;
            found_alts[alt_idx] = true;
          }
        }
        maf_vecs[var_idx].push_back( maf_vec[maf_idx] );
        max_vecs[var_idx].push_back( max_vec[maf_idx] );
        min_vecs[var_idx].push_back( min_vec[maf_idx] );
        
        if ( p1 && !found_ref ) Rcpp::Rcerr
          << func_tag << "\tMAF_1["<<var_idx<<"]["<<maf_idx<<"]: '"
          << maf_vec[maf_idx] << "','"
          << max_vec[maf_idx] << "','"
          << min_vec[maf_idx] << "'," << maf_vec[maf_idx].size()
          << std::endl;
      }
    }
    
    bool found_any = false;
    if ( !found_ref && maf_size != 0 )
    {
      found_any = true;
      maf_vecs[var_idx].push_back( "-1" );
      max_vecs[var_idx].push_back( max_vec[0] );
      min_vecs[var_idx].push_back( ref_var_vec[var_idx] );
    }
    
    if ( !found_alt && maf_size != 0 )
    {
      found_any = true;
      for ( std::size_t alt_idx = 0; alt_idx < alt_size; alt_idx++ )
      {
        if ( !found_alts[alt_idx] ) 
        {
          maf_vecs[var_idx].push_back( "-1" );
          max_vecs[var_idx].push_back( max_vec[0] );
          min_vecs[var_idx].push_back( alt_vec[alt_idx] );
        }
      }
    }
    
    if ( found_any )
    {
      Rcpp::Rcerr << func_tag << "\n\n\n\nHEREREREERERERE\n\n\n\n" << std::endl;
      
      unsigned int maf_size = maf_vec.size();
      for ( std::size_t maf_idx = 0; maf_idx < maf_size; maf_idx++ )
      {
        
        if ( p1 ) Rcpp::Rcerr
          << func_tag << "\tMAF_2["<<var_idx<<"]["<<maf_idx<<"]: '"
          << maf_vecs[var_idx][maf_idx] << "','"
          << max_vecs[var_idx][maf_idx] << "','"
          << min_vecs[var_idx][maf_idx] << "'," 
          << maf_vecs[var_idx][maf_idx].size()
          << std::endl;
        
      } // END:: for ( std::size_t maf_idx = 0; maf_idx < maf_size; maf_idx++ )
    } // END:: if ( found_any )
    
  } // END:: for ( std::size_t var_idx = 0; var_idx < df_ncols; var_idx++ )
  
  return( df );
}

// [[Rcpp::export]]
Rcpp::DataFrame fwd2tops_cpp( Rcpp::StringVector  seq_vec_r,
                              Rcpp::IntegerVector pos_vec_r,
                              Nullable<Rcpp::StringVector> cgn_vec_r_ = R_NilValue,
                              
                              const std::string cgn_key = "Probe_ID",
                              const int var_len = 2,
                              const bool return_source = false,
                              const bool uc = false,
                              const unsigned int vb = 0,
                              const unsigned int vt = 1,
                              const unsigned int tc = 0,
                              const std::string ft="fwd2tops_cpp" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0  = vb > vt + 0;
  const bool p1  = vb > vt + 1;
  const bool p2  = vb > vt + 2;
  const bool p10 = vb > p1 + 9;
  
  bool success = true;
  
  if ( !init_rvcp ) revcomp_init();
  if ( !init_topbot ) topbot_init();
  
  if ( p0 ) Rcpp::Rcerr <<_fo<< "Starting..." << "\n";
  
  // Return Data Frame for R::
  Rcpp::DataFrame ret_df;
  
  // Sanity Check:: Quick input validation::
  unsigned int seq_vec_size = seq_vec_r.size();
  unsigned int pos_vec_size = pos_vec_r.size();
  
  if ( seq_vec_size == 0 || seq_vec_size != pos_vec_size )
  {
    Rcerr << _fe << "Input vector sizes do not match or are zero!" << "\n"
          << _fe << "\t seq_vec_size = " << seq_vec_size << "\n"
          << _fe << "\t pos_vec_size = " << pos_vec_size << "\n"
          << std::endl;
    return( ret_df );
  }
  
  if ( p0 ) {
    Rcpp::Rcerr <<_fo<< "Input vector sizes match and not length zero!" << "\n"
                <<_fo<< "\t seq_vec_size = " << seq_vec_size << "\n"
                <<_fo<< "\t pos_vec_size = " << pos_vec_size << "\n"
                << std::endl;
  }
  
  // Allocate Conversion from R to c++ sequence vector::
  std::vector<std::string> seq_vec( seq_vec_size );
  
  // NOTE:: Expand all forward sequences first::
  // NOTE:: Double the size off the bat for expanded seqs vec::
  std::vector<std::string>  fwd_vec( seq_vec_size * 2 );
  std::vector<std::string>  exp_vec( seq_vec_size * 2 );
  std::vector<unsigned int> pos_vec( pos_vec_size * 2 );
  std::vector<unsigned int> ord_vec( seq_vec_size * 2 );
  std::vector<std::string>  cgn_vec( seq_vec_size * 2 );
  
  
  bool add_cgn = false;
  std::vector<std::string> ids_vec;
  if ( cgn_vec_r_.isNotNull() ) {
    if ( p1 ) Rcpp::Rcerr <<_fo<< "cgn_vec_r_ is NOT NULL\n" <<std::endl;
    Rcpp::StringVector cgn_vec_r(cgn_vec_r_);
    rcpp_to_cpp_vec_str( cgn_vec_r, ids_vec );
    add_cgn = true;
    
    // Left off here. Need to manage cgn output to ret_df...
  } else {
    if ( p1 ) Rcpp::Rcerr <<_fo<< "cgn_vec_r_ is NULL\n" <<std::endl;
  }
  
  unsigned int exp_raw_cnts = 0;
  unsigned int exp_last_idx = 0;
  unsigned int exp_vec_size = exp_vec.size();
  std::string cgn_str = "";
  
  for ( unsigned int seq_idx = 0; seq_idx < seq_vec_size; seq_idx++ ) {
    // R to c++ conversion::
    seq_vec[seq_idx] = seq_vec_r( seq_idx );
    unsigned int cur_pos = pos_vec_r( seq_idx );
    
    if ( add_cgn ) cgn_str = ids_vec[seq_idx];
    
    if ( p1 ) {
      Rcpp::Rcerr <<_fo<< "Current (" << seq_idx << ")" << "\n"
                  <<_fo<< "\t seq_vec = " << seq_vec[seq_idx] << "\n"
                  <<_fo<< "\t cur_pos = " << cur_pos << "\n"
                  << std::endl;
    }
    
    std::vector<std::string> cur_vec = expand_iupac_seq( seq_vec[seq_idx], vb, vt+1 );
    unsigned int cur_vec_size = cur_vec.size();
    
    if ( p1 ) {
      Rcpp::Rcerr <<_fo<< "IUPAC Expanded: (seq_idx, cur_vec_size) = (" 
                  << seq_idx << ", " << cur_vec_size << ")" << "\n"
      
      <<_fo<< "\t seq_vec = " << seq_vec[seq_idx] << "\n"
      <<_fo<< "\t cur_pos = " << cur_pos << "\n"
      << std::endl;
    }
    
    for ( unsigned int exp_idx = 0; exp_idx < cur_vec_size; exp_idx++ )
    {
      // unsigned int new_idx = seq_idx + exp_idx;
      // unsigned int new_idx = exp_last_idx + exp_idx;
      unsigned int new_idx = exp_last_idx;
      
      if ( p10 ) {
        Rcpp::Rcerr <<_fo<< "\t new_idx, seq_idx, exp_idx ("
                    << new_idx << ", " << seq_idx << ", " << exp_idx << ")" << "\n"
                    <<_fo<< "\t\t seq_vec = " << seq_vec[ seq_idx ] << "\n"
                    <<_fo<< "\t\t cur_pos = " << cur_pos << "\n"
                    << std::endl;
      }
      
      // Re-size if needed, hopefully not since we did this upfront for 2x...
      if ( new_idx + 2 > exp_vec_size )
      {
        unsigned int new_vec_size = (new_idx + 1) * 2;
        if ( p2 ) {
          Rcpp::Rcerr <<_fo<< "Expanding Sequence/Position Vectors::" << "\n"
                      <<_fo<< "\t exp_vec_size = " << exp_vec_size << "\n"
                      <<_fo<< "\t new_vec_size = " << new_vec_size << "\n"
                      << std::endl;
        }
        fwd_vec.resize( new_vec_size );
        exp_vec.resize( new_vec_size );
        pos_vec.resize( new_vec_size );
        ord_vec.resize( new_vec_size );
        if ( add_cgn ) cgn_vec.resize( new_vec_size );
      }
      
      // Assign new expanded sequences and coordinates::
      fwd_vec[ new_idx ] = seq_vec[seq_idx];
      exp_vec[ new_idx ] = cur_vec[exp_idx];
      pos_vec[ new_idx ] = cur_pos;
      ord_vec[ new_idx ] = exp_idx;
      if ( add_cgn ) cgn_vec[ new_idx ] = cgn_str;
      exp_raw_cnts++;
      
      // if ( new_idx > exp_last_idx ) exp_last_idx = new_idx;
      exp_last_idx++;
      
      if ( p10 ) {
        Rcpp::Rcerr <<_fo<< "\t new_idx, seq_idx, exp_idx ("
                    << new_idx << ", " << seq_idx << ", " << exp_idx << ")" << "\n"
                    <<_fo<< "\t\t seq_vec = " << seq_vec[ seq_idx ] << "\n"
                    <<_fo<< "\t\t fwd_vec = " << fwd_vec[ new_idx ] << "\n"
                    <<_fo<< "\t\t exp_vec = " << exp_vec[ new_idx ] << "\n"
                    <<_fo<< "\t\t pos_vec = " << pos_vec[ new_idx ] << "\n"
                    <<_fo<< "\t\t ord_vec = " << ord_vec[ new_idx ] << "\n"
                    << std::endl;
      }
    }
  }
  
  // NOTE:: Might need to drop the + 1
  fwd_vec.resize( exp_last_idx );
  exp_vec.resize( exp_last_idx );
  pos_vec.resize( exp_last_idx );
  ord_vec.resize( exp_last_idx );
  if ( add_cgn ) cgn_vec.resize( exp_last_idx );
  
  unsigned int new_fwd_size = fwd_vec.size();
  unsigned int new_exp_size = exp_vec.size();
  unsigned int new_pos_size = pos_vec.size();
  unsigned int new_ord_size = ord_vec.size();
  
  // TBD:: Should check all three new sizes are actually the same...
  if ( p1 ) {
    Rcpp::Rcerr << "\n"
                <<_fo<< "Done Expanding Forward Sequences!" << "\n"
                <<_fo<< "\t seq_vec_size = " << seq_vec_size << "\n" 
                <<_fo<< "\t pos_vec_size = " << pos_vec_size << "\n"
                <<_fo<< "\t new_fwd_size = " << new_fwd_size << "\n"
                <<_fo<< "\t new_exp_size = " << new_exp_size << "\n"
                <<_fo<< "\t new_pos_size = " << new_pos_size << "\n"
                <<_fo<< "\t new_ord_size = " << new_ord_size << "\n"
                <<_fo<< "\t exp_last_idx = " << exp_last_idx << "\n"
                <<_fo<< "\t exp_raw_cnts = " << exp_raw_cnts << "\n"
                << std::endl;
  }
  
  /*
   * Now make top/bot assignments and write out data to file if provided...
   * 
   */
  
  // Allocate top/bot calling/sequence vector::
  std::vector<std::string> tbs_vec( new_exp_size );
  std::vector<std::string> top_vec( new_exp_size );
  
  for ( unsigned exp_idx = 0; exp_idx < new_exp_size; exp_idx++ ) {
    top_vec[exp_idx] = exp_vec[exp_idx];
    
    // if ( p10 || ord_vec[exp_idx] > 0 ) {
    if ( p10 ) {
      Rcpp::Rcerr <<_fo<< "Pre-Top/Bot::\n" 
                  <<_fo<< " Exp_Idx, Ord_Idx = (" << exp_idx << ", " << ord_vec[exp_idx] << ")\n"
                  <<_fo<< "\t fwd_vec = " << fwd_vec[exp_idx] << "\n"
                  <<_fo<< "\t exp_vec = " << exp_vec[exp_idx] << "\n"
                  <<_fo<< "\t top_vec = " << top_vec[exp_idx] << "\n"
                  <<_fo<< "\t pos_vec = " << pos_vec[exp_idx] << "\n"
                  <<_fo<< "\t tbs_vec = " << tbs_vec[exp_idx] << "\n"
                  << std::endl;
    }
    
    int tb = is_top_bot( top_vec[exp_idx], var_len, uc, vb, vt+4 );
    if ( tb < 0 ) {
      tbs_vec[exp_idx] = "B";
      mutate_seq( top_vec[exp_idx], "rc", uc );
    } else if ( tb > 0 ) {
      tbs_vec[exp_idx] = "T";
    } else {
      tbs_vec[exp_idx] = "N";
    }
    
    // if ( p10 || ord_vec[exp_idx] > 0 ) {
    if ( p10 ) {
      Rcpp::Rcerr <<_fo<< " Exp_Idx, Ord_Idx, tb = (" 
                  << exp_idx << ", " << ord_vec[exp_idx] << ", " << tb << ")\n"
      
      <<_fo<< "\t fwd_vec = " << fwd_vec[exp_idx] << "\n"
      <<_fo<< "\t exp_vec = " << exp_vec[exp_idx] << "\n"
      <<_fo<< "\t top_vec = " << top_vec[exp_idx] << "\n"
      <<_fo<< "\t pos_vec = " << pos_vec[exp_idx] << "\n"
      <<_fo<< "\t tbs_vec = " << tbs_vec[exp_idx] << "\n"
      << std::endl;
    }
  }
  
  if ( !success ) {
    Rcpp::Rcerr
    <<_fe<< "Failed during processing!\n"
    <<std::endl;
    return( ret_df );
  }
  
  if ( return_source ) {
    ret_df = Rcpp::DataFrame::create(
      Rcpp::Named("Fwd_Sequence") = fwd_vec,
      Rcpp::Named("Exp_Sequence") = exp_vec,
      Rcpp::Named("Top_Sequence") = top_vec,
      Rcpp::Named("Cpg_Position") = pos_vec,
      Rcpp::Named("Strand_TB")    = tbs_vec );
  } else {
    ret_df = Rcpp::DataFrame::create(
      Rcpp::Named("Top_Sequence") = top_vec,
      Rcpp::Named("Strand_TB")    = tbs_vec );
  }
  if ( add_cgn ) ret_df[cgn_key] = cgn_vec;
  
  return( ret_df );
}

// [[Rcpp::export]]
StringVector expand_seqs_cpp( Rcpp::StringVector seqs,
                              const int max = 0,
                              const int vb = 0, 
                              const int vt = 1 ) {
  
  const std::string func_tag("[expand_seqs_cpp]: ");
  const std::string func_err("[expand_seqs_cpp]: ERROR: ");
  const std::string func_wrn("[expand_seqs_cpp]: Warning: ");
  
  Rcpp::Rcerr << "This function is deprecated...\n" << std::endl;
  StringVector rt( 0 );
  return(rt);
  
  if ( !init_exp_iupac ) iupac_exp_init();
  
  unsigned int seqs_size = seqs.size();
  
  std::vector<std::string> seqs_vec(seqs_size);
  
  for (unsigned int ii = 0; ii < seqs_size; ii++) {
    seqs_vec[ii] = seqs(ii);
    
    unsigned int rep_num = 1;
    std::string cur_seq = seqs_vec[ii];
    // unsigned int cur_len = cur_seq.size();
    
    std::vector<std::string> exps_vec = expand_iupac_seq( cur_seq, vb, vt+1 );
    
    // exps_vec.push_back(cur_seq);
    // for (unsigned int jj = 0; jj < cur_len; jj++) {
    //   unsigned int exp_len = iupac_exp_len[ (unsigned char) cur_seq[jj] ];
    //   if ( exp_len > 1 ) {
    //     rep_num = rep_num * exp_len;
    //     std::vector<char> rep_nucs = iupac_exp_nuc[ (char) cur_seq[jj] ];
    //     expand_seq_vec( exps_vec, rep_nucs, jj );
    //   }
    // }
    Rcpp::Rcerr << "ii = " << ii << "; seq = " << cur_seq << "\n";
    Rcpp::Rcerr << "\t rep_num = " << rep_num << std::endl;
    
    for (unsigned int jj = 0; jj < exps_vec.size(); jj++) {
      Rcpp::Rcerr << "\t jj = " << jj << "; exps_vec = " << exps_vec[jj] << std::endl;
    }
    Rcpp::Rcerr << "\n" << std::endl;
  }
  
  // Convert back to R-base Sting Vector::
  StringVector ret_vec( seqs_size );
  ret_vec = seqs_vec;
  return ret_vec;
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                              File IO Functions::
 *
 * Exported Rcpp Functions::
 *   variableColumnList()
 *   variableColumnListAsTibble()
 *   read_file_fast_cpp()
 *   read_file_to_df_cpp()
 *   cgntop_file_to_bgz_rcpp()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// [[Rcpp::export]]
Rcpp::List variableColumnList( std::vector<std::vector<std::string>>& v, 
                               std::vector<std::string>& names) {
  List retval;
  for (int ii = 0; ii <names.size(); ii++) {
    retval.push_back( v[ii], names[ii] );
  }
  return retval;
}

// [[Rcpp::export]]
Rcpp::DataFrame variableColumnListAsTibble( std::vector<std::vector<std::string>>& v, 
                                            std::vector<std::string>& names) {
  Function asTibble("as_tibble");
  
  return asTibble(variableColumnList(v, names));
}

//[[Rcpp::export]]
std::string read_file_fast_cpp( std::string path,
                                const int vb = 0,
                                const int vt = 3 ) {
  
  const std::string func_tag("[read_file_fast_cpp]: ");
  const std::string func_err("[read_file_fast_cpp]: ERROR: ");
  const std::string func_wrn("[read_file_fast_cpp]: Warning: ");
  
  std::ifstream in(path.c_str());
  std::string contents;
  
  // Scan Complete File from Begining to End::
  in.seekg(0,std::ios::end);
  
  if ( vb >= vt )
    Rcpp::Rcerr << func_tag << "\tFile Size = " << in.tellg() << "\n\n";
  
  // Adjust the size of the resulting string with the file size(tellg function)
  contents.resize( in.tellg() );
  
  // Return to the beginning of the file::
  in.seekg(0,std::ios::beg);
  
  // Each value the pointer will copied into contents
  in.read(&contents[0],contents.size());
  
  // Close file handle::
  in.close();
  
  return( contents );
}

//[[Rcpp::export]]
Rcpp::DataFrame read_fast_files_cpp( Rcpp::StringVector paths_vec,
                                     const std::string header_line = "",
                                     const char sep = ',',
                                     const int size = 2,
                                     const int line = 0,
                                     const int vb = 0, 
                                     const int vt = 1 ) {
  
  const std::string func_tag("[read_fast_files_cpp]: ");
  const std::string func_err("[read_fast_files_cpp]: ERROR: ");
  const std::string func_wrn("[read_fast_files_cpp]: Warning: ");
  
  const bool p1 = vb >= vt + 1;
  const bool p2 = vb >= vt + 1;
  
  Rcpp::DataFrame errs_ret( 0 );
  Rcpp::DataFrame tib;
  
  Rcpp::Rcerr << func_tag << "Under Construction...\n\n";
  return( tib );
  
  std::vector<std::string> paths(paths_vec.size());
  int file_count = paths_vec.size();
  
  if ( p1 )
    Rcpp::Rcerr << func_tag << "Starting: file count = " << file_count << "\n\n";
  
  bool has_header = true;
  unsigned int data_idx = 1;
  std::vector<std::string> header;
  if ( header_line.size() != 0 )
  {
    if ( p1 ) std::cerr 
      << func_tag << "Will use header line: '" <<  header_line << "'...\n" 
      << std::endl;
    
    data_idx = 0;
    has_header = false;
    parse_line(header_line, sep, header);
  }
  else
  {
    if ( p1 ) std::cerr 
      << func_tag << "Expecting header to be provided in file...\n" 
      << std::endl;
  }
  
  // return( header );
  // int header_numrow = 0;
  int column_numrow = 0;
  int header_length = 0;
  int column_length = 0;
  
  std::vector<std::vector<std::string>> columns_vec;
  std::vector<std::string> headers_vec;
  
  for ( int file_idx = 0; file_idx < file_count; file_idx++ ) {
    paths[file_idx] = paths_vec(file_idx);
    bool gzipped = is_gzipped( paths[file_idx] );
    
    std::string contents;
    if ( gzipped ) {
      std::string gzip_cmd = "gzip -dc " + paths[file_idx];
      contents = exec_cmd( gzip_cmd, vb, vt+1 );
    } else {
      contents = read_file_fast( paths[file_idx], vb, vt+1 );
    }
    
    // Parse contents into Lines::
    std::vector<std::string> lines;
    parse_line( contents, '\n', lines );
    
    // Parse File Header into Vector::
    if ( has_header ) parse_line(lines[0], sep, header);
    
    // Parse Columns into Vectors
    std::vector< std::vector< std::string> > column( header.size() );
    for (int row = data_idx; row < lines.size(); row++) {
      parse_into_postion( lines[row], ',', column );
    }
    
    if ( p2 ) Rcpp::Rcerr 
      << func_tag << "* ***** ***** ***** *****|***** ***** ***** ** *\n"
      << func_tag << "Done. Reading File (" << file_idx << ") '"
      <<  paths[file_idx] << "'...\n"
      << func_tag << "\tContents Size = " << contents.size() << "\n"
      << func_tag << "\t Header Count = " << header.size() << "\n"
      << func_tag << "\t   Line Count = " << lines.size() << "\n"
      << func_tag << "\t    Header[0] = " << header[0] << "\n"
      << func_tag << "\t    Header[n] = " << header[header.size()-1] << "\n"
      << std::endl;
    
    const int header_size = header.size();
    const int column_size = column.size();
    
    if ( header_size != column_size ) {
      Rcerr << "\n"
            << func_err << "Header and column lengths are not equal!\n"
            << func_err << "\theader length = " << header_size << "\n"
            << func_err << "\tcolumn length = " << column_size << "\n"
            << func_err << "Exiting... " << "\n"
            << std::endl;
      return( errs_ret );
    }
    if ( header_size == 0 && column_size == 0 ) {
      Rcerr << "\n"
            << func_err << "Both Header and column lengths are zero!!\n"
            << func_err << "\theader length = " << header_size << "\n"
            << func_err << "\tcolumn length = " << column_size << "\n"
            << func_err << "Exiting... " << "\n"
            << std::endl;
      return( errs_ret );
    }
    
    // Initialize header and columns vectors::
    if ( file_idx == 0 ) {
      std::copy( header.begin(), header.end(), std::back_inserter(headers_vec) );
      columns_vec.resize( column_size, std::vector<std::string>(0) );
    } else {
      /*
       * TBD:: Extend Header by adding any new headers.
       *   - Also, add blank columns for missing data
       *   - And add blank columns for new data that is missing...
       */
    }
    
    // Extend Columns::
    for (int col = 0; col < column_size; col++) {
      columns_vec[col].insert( columns_vec[col].end(),
                               std::make_move_iterator( column[col].begin() ),
                               std::make_move_iterator( column[col].end() ) );
    }
    
    // Update sizes::
    column_numrow = columns_vec[0].size();
    column_length = columns_vec[0][0].size();
    // header_numrow = headers_vec.size();
    header_length = headers_vec[0].size();
    
    if ( p2 ) Rcpp::Rcerr 
      << func_tag << "Vector Join Stats\n"
      << func_tag << "\tColumn Num Rows0 = " << column_numrow << "\n"
      << func_tag << "\tColumn Length[0] = " << column_length << "\n"
      << func_tag << "\tHeader Length[0] = " << header_length << "\n\n"
      << func_tag << "* ***** ***** ***** *****|***** ***** ***** ** *\n"
      << std::endl;
  }
  
  // Format tibble output::
  tib = variableColumnListAsTibble( columns_vec, headers_vec );
  
  return( tib );
}

// [[Rcpp::export]]
Rcpp::DataFrame read_file_to_df_cpp( const std::string file,
                                     Nullable<Rcpp::StringVector> name_vec_r_ = R_NilValue,
                                     const char sep = ',',
                                     const char ret = '\n',
                                     
                                     const bool uc = false,
                                     const unsigned int vb = 0,
                                     const unsigned int vt = 1,
                                     const unsigned int tc = 0,
                                     const std::string ft="read_file_to_df_cpp" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0 = vb > vt + 0;
  const bool p1 = vb > vt + 1;
  const bool p8 = vb > vt + 8;
  
  bool success = true;
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<< "Starting...\n"
    <<_fo<< "\t file = '"<<file<<"'\n"
    <<std::endl;
  if ( p8 ) Rcpp::Rcerr
    <<_fo<< "Verbosity...\n"
    <<_fo<< "\t p0 = "<<p0<<"\n"
    <<_fo<< "\t p1 = "<<p1<<"\n"
    <<_fo<< "\t p8 = "<<p8<<"\n"
    <<std::endl;
  
  Rcpp::DataFrame df = DataFrame::create();
  std::vector<std::vector<std::string>> data_mat;
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   * 
   *                 Allocate c++ vectors from Rcpp vectors::
   *                 
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  std::vector<std::string> name_vec;
  if ( name_vec_r_.isNotNull() ) {
    if ( p1 ) Rcpp::Rcerr <<_fo<< "name_vec_r_ is NOT NULL\n" <<std::endl;
    Rcpp::StringVector name_vec_r(name_vec_r_);
    rcpp_to_cpp_vec_str( name_vec_r, name_vec );
    
    if ( p1 ) Rcpp::Rcerr
      <<_fo<< "Data Delimiter='"<<sep<<"'\n"
      <<_fo<< " Column Vector='"<<name_vec_r<<"'\n"
      <<std::endl;
    
  } else {
    if ( p1 ) Rcpp::Rcerr <<_fo<< "name_vec_r_ is NULL\n" <<std::endl;
  }
  
  /* ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   * 
   *                           Load file via zlib::
   *                           
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  success = read_file_zlib( file, name_vec, data_mat, sep, ret, vb,vt+1,tc+1 );
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Status...\n"
    <<_fo<< "\t data_mat_size = "<<data_mat.size()<<"\n"
    <<std::endl;
  
  if ( !success ) {
    Rcpp::Rcerr
    <<_fe<< "Failed during processing!\n"
    <<std::endl;
    return( df );
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *
   *                 Convert Final Output to Return Rcpp DataFrame::
   * 
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  //
  // TBD:: Investigate setting the variable type (int,double,string)
  //
  for ( _I ii = 0; ii < name_vec.size(); ii++ ) {
    df.push_back( data_mat[ii], name_vec[ii] );
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                                   Done::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Done!\n"
    <<std::endl;
  
  return( df );
}

// [[Rcpp::export]]
bool cgntop_file_to_bgz_rcpp( const std::string file,
                              const std::string bgz_path,
                              const unsigned int top_idx = 0,
                              const char sep = ',',
                              const char ret = '\n',
                              
                              const bool uc = false,
                              const unsigned int vb = 0,
                              const unsigned int vt = 1,
                              const unsigned int tc = 0,
                              const std::string ft="cgntop_file_to_bgz_rcpp" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0 = vb > vt + 0;
  bool success = true;
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<< "Starting...\n"
    <<_fo<< "\t    file = '"<<file<<"'\n"
    <<_fo<< "\t bgz_path = '"<<bgz_path<<"'\n"
    <<_fo<< "\t top_idx = '"<<top_idx<<"'\n"
    <<std::endl;
  
  /* ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                       Load/Write file via zlib::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  success = cgntop_file_to_bgz( file, bgz_path, top_idx,
                                sep, ret, vb,vt+1,tc+1 );
  
  
  if ( !success ) {
    Rcpp::Rcerr <<_fe<< "Failed during processing!\n" <<std::endl;
    return( false );
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                                   Done::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<< "Done!\n"
    <<std::endl;
  
  return( success );
}

// [[Rcpp::export]]
bool load_cgntop_bgz_rcpp( const std::string file,
                           
                           const bool uc = false,
                           const unsigned int vb = 0,
                           const unsigned int vt = 1,
                           const unsigned int tc = 0,
                           const std::string ft="load_cgntop_bgz_rcpp" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p0 = vb > vt + 0;
  const bool p1 = vb > vt + 1;
  bool success = true;
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<< "Starting...\n"
    <<_fo<< "\t    file = '"<<file<<"'\n"
    <<std::endl;
  
  /* ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                       Load/Write file via zlib::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  std::vector<std::string> top_vec;
  success = load_cgntop_bgz( file, top_vec, vb,vt+1,tc+1 );
  
  _I top_cnt = top_vec.size();
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Status...\n"
    <<_fo<< "\t top_cnt = "<<top_cnt<<"\n"
    <<std::endl;
  
  if ( !success ) {
    Rcpp::Rcerr <<_fe<< "Failed during processing!\n" <<std::endl;
    return( false );
  }
  
  /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
   *                                   Done::
   * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
  
  if ( p0 ) Rcpp::Rcerr
    <<_fo<< "Done!\n"
    <<std::endl;
  
  return( success );
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                              File IO Functions::
 *
 * Exported Rcpp Functions::
 *   [Questionable]: write_df_bgz_rcpp()
 *   [Questionable]: read_df_bgz_rcpp()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// Environment pkg_dplyr = Environment::namespace_env("tidyverse");
// Environment pkg_dplyr = Environment::namespace_env("dplyr");
// Function dplyr_mutate("mutate");
// Environment pkg_purrr = Environment::namespace_env("purrr");
// Function purrr_across("across");
// Function purrr_is_double("is_double");

// Environment pkg = Environment::namespace_env("readr");
// Function read_tsv_r("read_tsv");
// Function cols_r("cols");
// Function col_character_r("col_character");
// Function col_integer_r("col_integer");

// [[Rcpp::export]]
bool write_df_bgz_rcpp( Rcpp::DataFrame df,
                        std::string spec_str,
                        std::string file,
                        
                        const unsigned int precision = 1000000,
                        const bool write_tsv = false,
                        
                        const unsigned int vb = 0,
                        const unsigned int vt = 1,
                        const unsigned int tc = 0,
                        const std::string ft="write_df_bgz_rcpp" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p1 = vb > vt + 0;
  
  bool success = true;  
  const _I nrows = df.nrows();
  const _I ncols = df.size();
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Starting...\n"
    <<_fo<< "\t     nrows = "<<nrows<<"\n"
    <<_fo<< "\t     ncols = "<<ncols<<"\n"
    <<_fo<< "\t  spec_str = "<<spec_str<<"\n"
    <<_fo<< "\t precision = "<<precision<<"\n"
    <<_fo<< "\t      file = "<<file<<"\n"
    <<std::endl;
  
  std::vector<std::string> name_vec = df.names();
  std::vector<std::vector<std::string>> str_vecs;
  std::vector<std::vector<int>> int_vecs;
  std::vector<std::vector<bool>> bool_vecs;
  
  for ( _I ii = 0; ii < spec_str.size(); ii++ )
  {
    char char_spec = spec_str[ii];
    
    if ( p1 ) Rcpp::Rcerr
      <<_fo<< "Current Vector["<<ii<<"]: name='"<<name_vec[ii]<<"', char_spec='"
      <<char_spec<<"'\n";
    
    if ( char_spec == 'c' )
    {
      Rcpp::CharacterVector cur_vec_r = df[ii];
      std::vector<std::string> cur_vec(nrows);
      rcpp_to_cpp_vec_str( cur_vec_r, cur_vec );
      
      str_vecs.push_back( cur_vec );
      
    } else if ( char_spec == 'i' || char_spec == 'd' || char_spec == 'n' ) {
      Rcpp::IntegerVector cur_vec_r = df[ii];
      std::vector<int> cur_vec(nrows);
      rcpp_to_cpp_vec_int( cur_vec_r, cur_vec );
      
      int_vecs.push_back( cur_vec );
      
    } else if ( char_spec == 'l' ) {
      Rcpp::LogicalVector cur_vec_r = df[ii];
      std::vector<bool> cur_vec(nrows);
      rcpp_to_cpp_vec_bool( cur_vec_r, cur_vec );
      
      bool_vecs.push_back( cur_vec );
      
    } else {
      if ( p1 ) Rcpp::Rcerr
        <<_fe<< "Unsupported character spec["<<ii<<"] ='"<<char_spec<<"'\n"
        <<std::endl;
      return( false );
    }
  }
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Converted all data to c++!\n"
    <<_fo<< "\t    str size = "<<str_vecs.size()<<"\n"
    <<_fo<< "\t    int size = "<<int_vecs.size()<<"\n"
    <<_fo<< "\t   bool size = "<<bool_vecs.size()<<"\n"
    <<std::endl;
  
  _I sum_check_cnt = spec_str.size() - str_vecs.size() - int_vecs.size() - bool_vecs.size();
  
  if ( sum_check_cnt != 0 ) {
    if ( p1 ) Rcpp::Rcerr
      <<_fe<< "Specs and Data sizes do not match = "<<sum_check_cnt<<"!\n"
      <<_fe<< "\t  spec_str.size  = "<<spec_str.size()<<"\n"
      <<_fe<< "\t  str_vecs.size  = "<<str_vecs.size()<<"\n"
      <<_fe<< "\t  int_vecs.size  = "<<int_vecs.size()<<"\n"
      <<_fe<< "\t  bool_vecs.size = "<<bool_vecs.size()<<"\n"
      <<std::endl;
    return( false );
  }
  
  success = write_df_bgz( spec_str,
                          file,
                          name_vec,
                          str_vecs,
                          int_vecs,
                          bool_vecs,
                          precision,
                          write_tsv,
                          vb,vt+1,tc );
  
  if ( !success ) {
    if ( p1 ) Rcpp::Rcerr
      <<_fe<< "Failed during write_df_bgz()\n"
      <<std::endl;
    return( false );
  }
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Done!\n"
  // <<_fo<< "\t   spec_str = "<<spec_str<<"\n"
     <<std::endl;
  
  return( success );
}

// [[Rcpp::export]]
Rcpp::DataFrame read_df_bgz_rcpp( std::string file,
                                  
                                  unsigned int precision = 1000000,
                                  
                                  const unsigned int vb = 0,
                                  const unsigned int vt = 1,
                                  const unsigned int tc = 0,
                                  const std::string ft="read_df_bgz_rcpp" )
{
  const std::string tb(tc, '\t');
  const std::string _fo("["+ft+"]: "+tb);
  const std::string _fe("["+ft+"]: ERROR: ");
  const std::string _fw("["+ft+"]: "+tb+"Warning: ");
  
  const bool p1 = vb > vt + 0;
  
  bool success = true;
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Starting...\n"
    <<_fo<< "\t      file = "<<file<<"\n"
    <<std::endl;
  
  Rcpp::DataFrame df =DataFrame::create();
  
  std::vector<std::string> name_vec;
  std::vector<std::vector<std::string>> str_vecs;
  std::vector<std::vector<int>> int_vecs;
  std::vector<std::vector<bool>> bool_vecs;
  
  success = load_df_bgz( file,
                         name_vec,
                         str_vecs,
                         int_vecs,
                         bool_vecs,
                         
                         vb+100,vt+1,tc+1 );
  
  if ( !success ) {
    if ( p1 ) Rcpp::Rcerr
      <<_fe<< "Failed during load_df_bgz()\n"
      <<std::endl;
    return( false );
  }
  
  for ( _I ii = 0; ii < name_vec.size(); ii++ ) {
    if ( p1 ) Rcpp::Rcerr
      <<_fo<<"df.int["<<ii<<"]: '"<<name_vec[ii]<<"' = '"<<int_vecs[ii].size()<<"'\n"
      <<std::endl;
    
    df.push_back( int_vecs[ii], name_vec[ii] );
  }
  
  // df = Rcpp::DataFrame::create(
  //   Rcpp::Named("Fwd_Sequence") = fwd_vec,
  //   Rcpp::Named("Exp_Sequence") = exp_vec,
  //   Rcpp::Named("Top_Sequence") = top_vec,
  //   Rcpp::Named("Cpg_Position") = pos_vec,
  //   Rcpp::Named("Strand_TB")    = tbs_vec );
  
  if ( p1 ) Rcpp::Rcerr
    <<_fo<< "Done!\n"
    <<std::endl;
  
  return( df );
};

/*** R

#
# [TBD]: Fix the BGZ read/write functions...
#
run_test <- TRUE
run_test <- FALSE

if ( run_test ) {
  #
  # BGZ Spec String Codes:: https://docs.python.org/3/library/struct.html
  #  c = character/string
  #  d = double
  #  i = integer
  #  n = ssize_t (integer)
  #  l = logical
  #
  suppressWarnings(suppressPackageStartupMessages( 
    base::require("tidyverse",  quietly = TRUE) ) )
  
  precision <- 1000000
  precision <- 1
  
  bgz_specs <- "iiiiii"
  # bgz_specs <- "diiddd" # Working
  # bgz_specs <- "cdiiddd"
  
  org_tib <- NULL
  org_tib <- tibble::tibble(
    # Probe_ID = c( "cg1", "cg23", "cg345" ),
    beta  = as.integer( c( 10000, 950001, 6000006 ) ),
    rank  = as.integer( c( 1, 2, 3) ),
    negs  = as.integer( c( -10, -20, -30) ),
    rank1 = as.integer( c( 4, 5, 6) ),
    rank2 = as.integer( c( 7, 8, 9) ),
    rank3 = as.integer( c( 2, 2, -2) )
    # pass = c( TRUE, FALSE, TRUE )
  )
  
  bgz_specs <- "iiiiii"
  org_tib <- NULL
  org_tib <- tibble::tibble(
    # Probe_ID = c( "cg1", "cg23", "cg345" ),
    beta  = as.integer( c(  1,  2,  3 ) ),
    rank  = as.integer( c(  4,  5,  6) ),
    negs  = as.integer( c(  7,  8,  9) ),
    rank1 = as.integer( c( 10, 11, 12) ),
    rank2 = as.integer( c( 13, 14, 15) ),
    rank3 = as.integer( c( 16, 17, 18) )
    # pass = c( TRUE, FALSE, TRUE )
  )
  # print( org_tib )
  
  out_path   <- "/Users/bbarnes/Documents/tmp/bgz"
  if ( !dir.exists(out_path) )
    dir.create( path = out_path, recursive = TRUE )
  bgz_prefix <- file.path( out_path, paste( "test",precision,bgz_specs, sep="_" ) )
  
  vb <- 4
  vt <- 0
  tc <- 0
  
  write_tsv <- TRUE
  
  write_bgz <- FALSE
  write_bgz <- TRUE
  
  read_bgz <- FALSE
  read_bgz <- TRUE
  
  success <- FALSE
  if ( write_bgz )
    success <- write_df_bgz_rcpp( df = org_tib,
                                  spec_str = bgz_specs,
                                  file = bgz_prefix,
                                  precision = precision, 
                                  write_tsv = write_tsv,
                                  vb=vb,vt=vt,tc=tc )
  
  new_tib <- NULL
  if ( success && read_bgz )
    new_tib <- read_df_bgz_rcpp( file = bgz_prefix,
                                 precision = precision,
                                 vb=vb,vt=vt,tc=tc ) %>% tibble::as_tibble()
  
  cat("Original Tibble::\n")
  print( org_tib )
  
  if ( !is.null(new_tib) ) {
    cat("Reloaded Tibble::\n")
    print(new_tib)
  }
  
  cat("\n\nERROR: There is an issue with repeating indexes in either the ",
      "writing or loading of the BGZ data::\n\n\n")
  
  # print( tibble::as_tibble(new_df) )
}

*/

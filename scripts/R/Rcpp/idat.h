#ifndef __IDAT_H__
#define __IDAT_H__

#include "bisulfite.h"

using namespace bisulfite;

namespace idat {

static _B MAGIC_SHIFT = 7;

/*
 * Funny Support Functions::
 * 
 */
_S char_to_bin_str( _B c )
{
  // Convert a char to ASCII value
  _i val = int(c);
  
  // Convert ASCII value to binary string
  _S bin_str = "";
  while (val > 0)
  {
    (val % 2)? bin_str.push_back('1') :
    bin_str.push_back('0');
    val /= 2;
  }
  reverse( bin_str.begin(), bin_str.end() );
  
  return( bin_str );
}

_S str_to_bin_str( _S s )
{
  _S bin_str = "";
  _I n = s.length();
  for ( std::size_t ii = 0; ii <= n; ii++) bin_str += char_to_bin_str( s[ii] );
  return( bin_str );
}

void mssg_idat( _B n, _B m, 
                _I& idx,
                _S  tag = "", 
                _S  pre = "",
                
                const _I vb = 0,
                const _I vt = 1,
                const _I tc = 0,
                const _S ft="mssg_idat",
                
                const _S sep0 = ",",
                const _S sepS = ", ",
                const _S sepQ = "',",
                const _S sepQS = "', " )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb > vt + 0;
  const _Y p1 = vb > vt + 1;
  
  _I n_size_of = sizeof(n);
  _I m_size_of = sizeof(m);
  
  _S pre_bin_str = str_to_bin_str( pre );
  
  if ( p0 ) std::cerr
    <<_fo<<"Starting::\n"
    <<_fo<<"\t   idx = '"<<idx<<"'\n"
    <<_fo<<"\t   tag = '"<<tag<<"'\n"
    <<_fo<<"\t   pre = '"<<pre_bin_str<<"'\n"
    <<std::endl;
  
  _S n_bin_str = "";
  try { 
    n_bin_str = char_to_bin_str(n);
  } catch (const std::exception& e) {
    std::cerr<<_fe<<"Caught Bin ERROR: n '"<<e.what()<<"'\n"<<std::endl;
  }
  _S m_bin_str = "";
  try { 
    m_bin_str = char_to_bin_str(m);
  } catch (const std::exception& e) {
    std::cerr<<_fe<<"Caught Bin ERROR: n '"<<e.what()<<"'\n"<<std::endl;
  }
  
  std::intmax_t n_bin_int = 0;
  try {
    n_bin_int = strtoimax( n_bin_str.c_str(), nullptr, 2 );
  } catch (const std::exception& e) {
    std::cerr<<_fe<<"Caught Int ERROR: n '"<<e.what()<<"'\n"<<std::endl;
  }
  std::intmax_t m_bin_int = 0;
  try {
    m_bin_int = strtoimax( m_bin_str.c_str(), nullptr, 2 );
  } catch (const std::exception& e) {
    std::cerr<<_fe<<"Caught Int ERROR: n '"<<e.what()<<"'\n"<<std::endl;
  }
  
  if ( p1 ) std::cerr
    <<_fo<<"\t n::"
    <<"n_size_of = '"<<n_size_of<<sepQS
    <<"n_bin_str = '"<<n_bin_str<<sepQS
    <<"n_bin_int = '"<<n_bin_int<<sepQS
    <<std::endl
    <<_fo<<"\t m::"
    <<"m_size_of = '"<<m_size_of<<sepQS
    <<"m_bin_str = '"<<m_bin_str<<sepQS
    <<"m_bin_int = '"<<m_bin_int<<sepQS
    <<std::endl;
  
  idx++;
}

/*
 * Safe Read "Funny Illumina" String Length Template for Binary Gzipped Files::
 * * IDAT Seek & Read:: String (one byte at a time:: SLOW!!!!)
 * 
 * Source: From [1] 
 *  https://code.google.com/p/glu-genetics/source/browse/glu/lib/illumina.py#86:
 * 
 * String data are encoded as a sequence of one or more length bytes followed                                                                                         
 * by the specified number of data bytes.                                                                                                                             
 *
 * The lower 7 bits of each length byte encodes the bits that comprise the                                                                                            
 * length of the following byte string.  When the most significant bit it                                                                                             
 * set, then an additional length byte follows with 7 additional high bits to                                                                                         
 * be added to the current length.  The following string lengths are                                                                                                  
 * accommodated by increasing sequences of length bytes:                                                                                                              
 *
 * length  maximum                                                                                                                                                    
 * bytes   length                                                                                                                                                     
 * ------  --------                                                                                                                                                   
 * 1       127 B                                                                                                                                                    
 * 2        16 KB                                                                                                                                                   
 * 3         2 MB                                                                                                                                                   
 * 4       256 MB                                                                                                                                                   
 * 5        32 GB                                                                                                                                                   
 *
 * While this seems like a sensible progression, there is some uncertainty                                                                                            
 * about this interpretation, since the longest of string observed in the                                                                                             
 * wild has been of length 6,264 with two length bytes.  
 *
 *
 * Python::
 *   n = m = ord(read(1))
 *
 *   if not n:
 *     return ''
 *
 *   if m&0x80:
 *     shift = 7
 *     n     = m&0x7F
 *
 *     while m&0x80:
 *       m      = ord(read(1))
 *       n     += (m&0x7F)<<shift
 *       shift += 7
 *
 *   return read(n)
 *   
 *   Example::
 *   
 *   $RunInfo
 * RunTime                 BlockType  BlockPars
 * [1,] "12/3/2020 9:22:01 PM"  "Decoding" "CallsToUsed=2062647|CallsToUnused=23928|CallsToInvalid=200593"
 * [2,] "1/27/2021 12:56:00 PM" "Scan"     "sherlockID=N285|ScannerID=N285|Username=svc_meta|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.5.0.2"
 * [3,] "1/27/2021 12:56:00 PM" "Register" "Algorithm=StandardGeneric"
 * [4,] "1/27/2021 12:56:00 PM" "Extract"  "Algorithm=StandardWithBackground"
 *
 * [5,] "12/3/2020 9:27:18 PM"  "Decoding" "CallsToUsed=2281795|CallsToUnused=26773|CallsToInvalid=176484"
 * [6,] "1/27/2021 12:56:17 PM" "Scan"     "sherlockID=N285|ScannerID=N285|Username=svc_meta|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.5.0.2"
 * [7,] "1/27/2021 12:56:17 PM" "Register" "Algorithm=StandardGeneric"
 * [8,] "1/27/2021 12:56:17 PM" "Extract"  "Algorithm=StandardWithBackground"
 *
 * [9,] "12/3/2020 9:38:21 PM"  "Decoding" "CallsToUsed=2360595|CallsToUnused=27299|CallsToInvalid=164996"
 * [10,] "1/27/2021 12:56:32 PM" "Scan"     "sherlockID=N285|ScannerID=N285|Username=svc_meta|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.5.0.2"
 * [11,] "1/27/2021 12:56:32 PM" "Register" "Algorithm=StandardGeneric"
 * [12,] "1/27/2021 12:56:32 PM" "Extract"  "Algorithm=StandardWithBackground"
 *
 * [13,] "12/3/2020 9:48:45 PM"  "Decoding" "CallsToUsed=2310028|CallsToUnused=28396|CallsToInvalid=185835"
 * [14,] "1/27/2021 12:56:46 PM" "Scan"     "sherlockID=N285|ScannerID=N285|Username=svc_meta|FPGAVersion=4.0.20|SoftwareAppication=iScan Control Software|SoftwareVersion=3.5.0.2"
 * [15,] "1/27/2021 12:56:46 PM" "Register" "Algorithm=StandardGeneric"
 * [16,] "1/27/2021 12:56:46 PM" "Extract"  "Algorithm=StandardWithBackground"
 * 
 */

/*
 * TBD::
 *   safe_vload_bgz()
 *   safe_sload_bgz()
 *   OR GO BACK TO THE ORIGNAL METHOD...
 * 
 * 
 */
template <typename F>
int safe_strlen_bgz( F& fh,
                     const _B mag = MAGIC_SHIFT,
                     
                     const _I vb = 0,
                     const _I vt = 4,
                     const _I tc = 0,
                     const _S ft="safe_strlen_bgz" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb > vt + 0;
  
  if ( p0 ) std::cerr
    <<_fo<<"Starting::\n"
    <<_fo<<"\t   Magic = '"<<std::to_string( mag )<<"'\n"
    <<std::endl;
  
  _i rv = 0;
  _Y success = true;
  // buf.clear();
  // buf = "";
  
  /*
   * Read Funny Bits (described above) to determine the number of characters
   *  to read and then read characters...
   */
  _B n, m;
  _S bit_str = "";
  _I idx = 0;
  
  try {
    rv = gzread( fh, &n, sizeof( n ) );
  } catch (const std::exception& e) {
    std::cerr<<_fe<<"Caught ERROR: n '"<<e.what()<<"'\n"<<std::endl;
    success = false;
  }
  
  m = n;
  if ( n==0 ) return( 0 );
  
  mssg_idat( n, m, idx, "1", bit_str, vb,vt+1,tc,ft );
  if ( m & 0x80 ) {
    _B shift = mag;
    n = m & 0x7F;
    
    bit_str = bit_str + char_to_bin_str( n );
    mssg_idat( n, m, idx, "2", bit_str, vb,vt+1,tc,ft );
    
    while ( m & 0x80 )
    {
      try {
        rv = gzread( fh, &m, sizeof( _B ) );
      } catch (const std::exception& e) {
        std::cerr<<_fe<<"Caught ERROR: n '"<<e.what()<<"'\n"<<std::endl;
        success = false;
      }
      n  += (m & 0x7F) << shift;
      shift += mag;
      
      bit_str = bit_str + char_to_bin_str( n );
      mssg_idat( n, m, idx, "3", bit_str, vb,vt+1,tc,ft );
    }
  }
  mssg_idat( n, m, idx, "4", bit_str, vb,vt+1,tc,ft );
  
  // Final Read::
  int len;
  try {
    len = strtoimax( char_to_bin_str(n).c_str(), nullptr, 2 );
  } catch (const std::exception& e) {
    std::cerr<<_fe<<"Caught ERROR: n '"<<e.what()<<"'\n"<<std::endl;
  }
  
  if ( p0 ) std::cerr
    <<_fo<<"      rv = '"<<rv<<"'\n"
    <<_fo<<"     len = '"<<len<<"'\n"
    <<_fo<<" success = '"<<success<<"'\n"
    <<_fo<<"Done.\n"<<std::endl;
  
  return( len );
};

/*
 * Safe Seek Funny String Binary Gizzped Files:: Scalar
 * 
 */
template <typename F>
static _Y safe_load_str_bgz( F& fh,
                             _S& buf,
                             _I rep_cnt = 1,
                             const _S bin_str = "",
                             
                             const _Q off_set = 0,
                             
                             const _Y seek_cur = false,
                             const _Y ilmn_str = false,
                             
                             const _I vb = 0,
                             const _I vt = 4,
                             const _I tc = 0,
                             const _S ft="safe_load_str_bgz" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb > vt + 0;
  
  if ( p0 ) std::cerr
    <<_fo<<"Starting::\n"
    <<_fo<<"\t   rep_cnt = '"<< rep_cnt<<"'\n"
    <<_fo<<"\t   bin_str = '"<< bin_str<<"'\n"
    <<_fo<<"\t   off_set = '"<< off_set<<"'\n"
    <<_fo<<"\t  seek_cur = '"<< seek_cur<<"'\n"
    <<_fo<<"\t  ilmn_str = '"<< ilmn_str<<"'"
    <<std::endl;
  
  // 1. Seek to the position if off_set is non-zero
  if ( off_set > 0 )
    if ( !safe_seek_bgz( fh, off_set, seek_cur, vb,vt+1,tc ) ) 
      return( false );
    
    // 2. Read Funny "Magic" Illumina Binary Strings length and then string
    if ( ilmn_str ) rep_cnt = safe_strlen_bgz( fh, MAGIC_SHIFT, vb,vt+1,tc );
    if ( ilmn_str && rep_cnt > 1 ) buf.resize( rep_cnt );
    
    // 3. Finally Read and store the data
    if ( !safe_read_bgz( fh, buf[0], rep_cnt, bin_str, vb,vt+1,tc ) )
      return( false );
    
    if ( p0 ) std::cerr<<_fo<<"Done.\n"<<std::endl;
    
    return( true );
};

/*
 * Safe Seek Funny String Binary Gizzped Files:: Vector
 * 
 */
template <typename F>
static _Y safe_load_strv_bgz( F& fh,
                              std::vector<_S>& buf,
                              _I rep_cnt = 1,
                              const _S bin_str = "",
                              
                              const _Q off_set = 0,
                              
                              const _Y seek_cur = false,
                              
                              const _I vb = 0,
                              const _I vt = 4,
                              const _I tc = 0,
                              const _S ft="safe_load_strv_bgz" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb > vt + 0;
  
  if ( p0 ) std::cerr
    <<_fo<<"Starting::\n"
    <<_fo<<"\t   rep_cnt = '"<< rep_cnt<<"'\n"
    <<_fo<<"\t   bin_str = '"<< bin_str<<"'\n"
    <<_fo<<"\t   off_set = '"<< off_set<<"'\n"
    <<_fo<<"\t  seek_cur = '"<< seek_cur<<"'"
    <<std::endl;
  
  buf.clear();
  buf.resize(rep_cnt);
  
  // 1. Seek to the position if off_set is non-zero
  if ( off_set > 0 )
    if ( !safe_seek_bgz( fh, off_set, seek_cur, vb,vt+1,tc ) )
      return( false );
    
    for ( std::size_t  rep_idx=0; rep_idx<rep_cnt; rep_idx++ ) {
      
      // 2. Read Funny "Magic" Illumina Binary Strings length and then string
      _I rep_len = safe_strlen_bgz( fh, MAGIC_SHIFT, vb,vt+1,tc );
      buf[rep_idx].clear();
      buf[rep_idx].resize( rep_len );
      
      // 3. Finally Read and store the data
      if ( !safe_read_bgz( fh, buf[rep_idx][0], rep_len, bin_str, vb,vt+1,tc ) )
        return( false );
    }
    if ( p0 ) std::cerr<<_fo<<"Done.\n"<<std::endl;
    
    return( true );
};

/*
 * Safe Seek Funny String Binary Gizzped Files:: Matrix
 * 
 */
template <typename F>
static _Y safe_load_strm_bgz( F& fh,
                              std::vector<std::vector<_S>>& buf,
                              _I vec_cnt = 1,
                              _I rep_cnt = 1,
                              const _S bin_str = "",
                              
                              const _Q off_set = 0,
                              
                              const _Y seek_cur = false,
                              
                              const _I vb = 0,
                              const _I vt = 4,
                              const _I tc = 0,
                              const _S ft="safe_load_strm_bgz" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb > vt + 0;
  
  if ( p0 ) std::cerr
    <<_fo<<"Starting::\n"
    <<_fo<<"\t   rep_cnt = '"<< rep_cnt<<"'\n"
    <<_fo<<"\t   bin_str = '"<< bin_str<<"'\n"
    <<_fo<<"\t   off_set = '"<< off_set<<"'\n"
    <<_fo<<"\t  seek_cur = '"<< seek_cur<<"'"
    <<std::endl;
  
  buf.clear();
  buf.resize(vec_cnt);
  
  // 1. Seek to the position if off_set is non-zero
  if ( off_set > 0 )
    if ( !safe_seek_bgz( fh, off_set, seek_cur, vb, vt+1, tc+1 ) ) return( false );
    
    for ( std::size_t rep_idx=0; rep_idx<vec_cnt; rep_idx++ ) {
      buf[rep_idx].clear();
      buf[rep_idx].resize(rep_cnt);
      
      _S idx_str( bin_str );
      idx_str += "-rep_idx="+std::to_string(rep_idx);
      
      if ( !safe_load_strv_bgz( fh, buf[rep_idx], rep_cnt, idx_str, 0,
                                false, vb,vt+3,tc ) ) return( false );
    }
    if ( p0 ) std::cerr<<_fo<<"Done.\n"<<std::endl;
    
    return( true );
};

/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *  Class:: Idat
 *   - Stores a signle idat ( one color channel)::
 * 
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
class Idat
{
private:

  // std::map< _S,_S > iparam;  
  // std::vector< _I > tangos;
  // std::vector< _H > signal;
  // std::vector< _d > pvalue;
  
  // Map of Key(field)/Value pairs for scalars::
  // [Key][Val]: Scalar Fields
  std::map<_S, _I > sI_dat;
  std::map<_S, _S > sS_dat;
  
  // Vector of Tango Addresses by class {All,Neg}
  // [Add]: All Tango Addresses
  // [Neg]: Neg Tango Addresses
  std::map<_S, std::vector<_I> > vI_dat;
  
  // Vectors of data[_H] {Sig,Sds,Rep} in order of Tango Address
  // [Sig]: All Signal Intesnsity,
  // [Sds]: All Signal Standard Deviation,
  // [Rep]: All Bead Replicate Counts
  std::map<_S, std::vector<_H> > vH_dat;
  
  // Vectors of data[_d] {Sig,Sds,Rep} in order of Tango Address
  // [Det]: All Detection P-values
  std::map<_S, std::vector<_d> > vd_dat;
  
  // Map of Tango Addresses to Vector Index::
  // [Map][Add]=[Idx]: Tango Address to Array Index
  std::map<_S, std::map<_I,_I> > add_map;
  
  // Matrix[_S][_S] of Run Info::
  std::vector< std::vector<_S> > run_mat;
  
public: 
  
  Idat( const _S idat_path,
        const _B col,
        const std::vector<_I>& pval_vec = {},
        const _Y rm_pval_outliers = false,
        
        const _I vb = 0,
        const _I vt = 1,
        const _I tc = 0,
        const _S ft="Idat" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    
    if ( p0 ) std::cerr
      <<_fo<<"Loading Idat::\n"
      <<_fo<<"\t (vb,vt,tc): ("<<vb<<","<<vt<<","<<tc<<")\n"
      <<_fo<<"\t  Idat_Path: '"<<idat_path<<"'\n"
      <<std::endl;
    
    this->parse_idat( idat_path,col,pval_vec, rm_pval_outliers, vb,vt+1,tc );
    
    if ( p0 ) std::cerr<<_fo<<"Done Loading Idat.\n"<<std::endl;
  };
  
  // Contstructor if Background Addresses are not provided::
  Idat( const _S idat_path,
        const _B col,
        const _Y rm_pval_outliers = false,

        const _I vb = 0,
        const _I vt = 0,
        const _I tc = 0,
        const _S ft = "AS_Idat" ) :
    Idat( idat_path, col,DEFAULT_UINT_VEC,rm_pval_outliers, vb,vt,tc ) {};
  
  ~Idat() {};

  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                      Idat Binary Parser:: Version 3.0
   * 
   * parse_idat()
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  _Y parse_idat( const _S idat_path,
                 const _B col,
                 const std::vector<_I>& pval_vec,
                 const _Y rm_pval_outliers = false,
                 
                 const _I vb = 0,
                 const _I vt = 1,
                 const _I tc = 0,
                 const _S ft="parse_idat" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    // const _Y p1 = vb > vt + 1;
    // const _Y p2 = vb > vt + 2;
    const _Y p3 = vb > vt + 3;
    const _Y p4 = vb > vt + 4;

    _Y success = true;
    
    if ( !legal_file( idat_path ) ) return( false );
    const _Y isgzipped = is_gzipped( idat_path );
    
    if ( p0 ) std::cerr 
      <<_fo<<"Setting Default Idat Storage Params...\n"
      <<_fo<<"\t  (vb,vt,tc): (" <<vb<<","<<vt<<","<<tc<<")\n"
      <<_fo<<"\t      Idat Path = '"<<idat_path<<"\n" 
      <<_fo<<"\t      IsGZipped = '"<<isgzipped<<"'"
      << std::endl;
    
    const _S IDAT_STR("IDAT");
    const _I IDAT_LEN = IDAT_STR.size();
    const _I IDAT_VER = 3;
    const _I HEAD_MIN = 7;
    const _I INFO_LEN = 5;

    /*
     * (Re)-Initialize All Address Space [ Idx,Avg,Sds,Rep ]
     * 
     */
    _I Address_Cnt = 0;
    _I Red_Grn_Int = 0;
    
    _L idat_version;
    _I header_cnt;
    
    _S Sentrix_Idx, Sentrix_Pos, Chip_Format, Sentrix_Man,
    OPA_Value_1, Description, PlateNumber, Well_Number,
    Unknown_001, Unknown_002;
    
    // Clear All just for good measure::
    Sentrix_Idx.clear(); Sentrix_Pos.clear();
    Sentrix_Man.clear(); Chip_Format.clear();
    OPA_Value_1.clear(); Description.clear();
    PlateNumber.clear(); Well_Number.clear();
    Unknown_001.clear(); Unknown_002.clear();
    
    // Address (Tango):: Identifiers
    std::map<_I,_I> Address_Map; // map( Tango Address => Vector Index )
    std::vector<_I> Address_Idx; // Tango Address
    std::vector<_I> MBlocks_Idx; // Some Idat Field, forgot the improtance???
    
    // Address (Tango):: Values
    std::vector<_H> Address_Avg; // Intensity Signal Average
    std::vector<_H> Address_Sds; // Intensity Standard Deviation
    std::vector<_H> Address_Rep; // Bead Replicate Count
    std::vector<_d> Address_Det; // Detection P-values (usually form Neg Ctls)
    
    // Detection P-value Address/Intensity Vector::
    std::vector<_I> Detection_Idx; // Tango Address
    std::vector<_H> Detection_Avg; // Intensity Signal Average
    
    // Run Info:: Matrix
    std::vector< std::vector<_S> > runinfo_mat;
    // Run Info:: Map
    std::map< std::vector<_S>, std::vector<_S> > runinfo_map;
    
    // Clear All just for good measure::
    Address_Idx.clear(); Address_Avg.clear();
    Address_Sds.clear(); Address_Rep.clear();
    Address_Det.clear(); Address_Map.clear();
    MBlocks_Idx.clear(); runinfo_mat.clear();
    Detection_Idx.clear(); Detection_Avg.clear();
    
    if ( isgzipped )
    {
      gzFile gzin = gzopen( idat_path.c_str(), "r" );
      
      // Validate IDAT Header String:: (valid == "IDAT")
      //
      std::vector<char> idat_vec(IDAT_LEN);
      success = safe_load_bgz( gzin, idat_vec[0], IDAT_LEN,
                               "header_str(4c)", 0, false, false, 
                               vb, vt+3, tc+1 );
      _S header_str( idat_vec.begin(), idat_vec.end() );
      if ( p3 ) std::cerr<<_fo<<"Header_Str = '"<<header_str<<"'"<<std::endl;
      
      if ( header_str.compare(IDAT_STR) != 0 ) {
        std::cerr<<_fe<<" FAILED: Header_Str = '"<<header_str<<"' "
                 <<"NOT Supported! Only Value=="<<IDAT_STR<<"!\n"<<std::endl;
        success = false;
        return( success );
      }
      
      // Validate IDAT Version Number:: (valid == 3)
      //
      success = safe_load_bgz( gzin, idat_version, 1, 
                               "Version(L)", 0, true, false, 
                               vb, vt+3, tc+1 );
      if ( p3 ) std::cerr<<_fo<<"Version = '"<<idat_version<<"'"<<std::endl;
      
      if ( idat_version != IDAT_VER ) {
        std::cerr<<_fe<<"FAILED: IDAT Version: '"<<idat_version<< "' "
                 <<"NOT Supported! Only Value=="<<IDAT_VER<<"!\n"<<std::endl;
        return( false );
      }
      
      // Load IDAT Header Count:: (valid >= 7)
      //   TBD:: Update the valid header count to something critical...
      success = safe_load_bgz( gzin, header_cnt, 1, 
                               "Header_Cnt(I)", 0, true, false, 
                               vb, vt+3, tc+1 );
      if ( p3 ) std::cerr<<_fo<<"Header_Cnt = '"<<header_cnt<<"'"<<std::endl;
      
      if ( header_cnt < HEAD_MIN ) {
        std::cerr<<_fe<<"FAILED: IDAT Header_Cnt: '"<<header_cnt<< "' "
                 <<"NOT Supported! Only Value>="<<HEAD_MIN<<"!\n"<<std::endl;
        return( false );
      }
      
      /*
       * Load IDAT Header Values::
       *   - IDAT_CODE
       *   - IDAT_OFFS
       * 
       */
      std::vector<_S> name_vec;
      std::vector<_B> type_vec;
      std::vector<_I> tidx_vec;
      std::vector<_I> size_vec;
      std::vector<_Q> offs_vec;
      
      std::map<_B, _I> tidx_cnts;
      
      // New Mapping Structures::
      std::vector<_S> hidx_to_name;
      std::map<_S, _B> name_to_type;
      std::map<_S, _I> name_to_hidx;
      std::map<_S, _I> name_to_tidx;
      std::map<_S, _I> name_to_size;
      std::map<_S, _Q> name_to_offs;
      
      _Q prev = 0;
      for ( std::size_t hidx=0; hidx<header_cnt; hidx++ ) {
        // Extract IDAT_CODE: Short
        _h code = 0;
        success = safe_load_bgz( gzin, code, 1, "code(h)", 0, true, false, 
                                 vb, vt+3, tc+1 );
        if ( p4 ) std::cerr
          <<_fo<<"\t IDAT Code["<<hidx<<"] = '"<<code<<"'"<<std::endl;
        if ( !success ) return( success );
        
        // Extract IDAT_OFFS:: _Q
        _Q offs = 0;
        success = safe_load_bgz( gzin, offs, 1, "offs(ULL)", 0, true, false, 
                                 vb, vt+3, tc+1 );
        if ( !success ) return( success );
        
        _Q diff = offs - prev;
        prev = offs;
        
        if ( p4 ) std::cerr
          <<_fo<<"\t   Dif_Set["<<hidx<<"] = '"<<diff<<"'"<< "\n"
          <<_fo<<"\t   Off_Set["<<hidx<<"] = '"<<offs<<"'"
          <<std::endl;
        
        // Switch Case (Code) => [ NAME, TYPE ]
        _S name;
        _B type;
        _I size;
        
        switch( code )
        {
        case 1000 :
          name = "Address_Cnt";
          type = 'i';
          size = 1;
          break;
          
        case 102 :
          name = "Address_Idx";
          type = 'I';
          size = 0;
          break;
          
        case 103 :
          name = "Address_Sds";
          type = 'H';
          size = 0;
          break;
          
        case 104 :
          name = "Address_Avg";
          type = 'H';
          size = 0;
          break;
          
        case 107 :
          name = "Address_Rep";
          type = 'H';
          size = 0;
          break;
          
        case 200 :
          name = "MidBlockIdx";
          type = 'S';
          size = 4;
          break;
          
        case 300 :
          name = "RunInfo_Idx";
          type = 'S';
          size = 2;
          break;
          
        case 400 :
          name = "Red_Grn_Int";
          type = 'i';
          size = 1;
          break;
          
        case 401 :
          name = "Sentrix_Man";
          type = 's';
          size = 1;
          break;
          
        case 402 :
          name = "Sentrix_Idx";
          type = 's';
          size = 1;
          break;
          
        case 403 :
          name = "Chip_Format";
          type = 's';
          size = 1;
          break;
          
        case 404 :
          name = "Sentrix_Pos";
          type = 's';
          size = 1;
          break;
          
        case 405 :
          name = "OPA_Value_1";
          type = 's';
          size = 1;
          break;
          
        case 410 :
          name = "Unknown_002";
          type = 's';
          size = 1;
          break;
          
        case 406 :
          name = "Sample_Id_1";
          type = 's';
          size = 1;
          break;
          
        case 407 :
          name = "Description";
          type = 's';
          size = 1;
          break;
          
        case 408 :
          name = "PlateNumber";
          type = 's';
          size = 1;
          break;
          
        case 409 :
          name = "Well_Number";
          type = 's';
          size = 1;
          break;
          
        case 510 :
          name = "Unknown_001";
          type = 's';
          size = 1;
          break;
          
        default :
          std::cerr <<_fe<<"Invalid IDAT Code: '"<<code<<"'\n" 
                    <<std::endl;
        }
        tidx_cnts[ type ]++;
        _I tidx = tidx_cnts[ type ];
        
        if ( p4 ) std::cerr
          <<_fo<<"\t IDAT Name["<<hidx<<"] = '"<<name<<"'\n"
          <<_fo<<"\t IDAT Type["<<hidx<<"] = '"<<type<<"'\n"
          <<_fo<<"\t IDAT Tidx["<<hidx<<"] = '"<<tidx<<"'"
          <<std::endl;
        
        name_vec.push_back( name );
        type_vec.push_back( type );
        tidx_vec.push_back( tidx );
        size_vec.push_back( size );
        offs_vec.push_back( offs );
        
        // New Mapping Structures::
        hidx_to_name.push_back( name );
        name_to_hidx.emplace( name, hidx );
        name_to_tidx.emplace( name, tidx );
        name_to_type.emplace( name, type );
        name_to_size.emplace( name, size );
        name_to_offs.emplace( name, offs );
      }
      
      if ( p4 ) std::cerr<<_fo<<"Header Loaded!!!\n"<<std::endl;
      
      /*
       * Read Address and Count and Allocate Address Space [ Idx,Avg,Sds,Rep ]
       * 
       */
      int tmp_add_cnt = 0;
      success = safe_load_bgz( gzin, tmp_add_cnt, 1, "Address_Cnt", 
                               name_to_offs["Address_Cnt"], false, false, 
                               vb, vt+3, tc+1 );
      Address_Cnt = (_I) tmp_add_cnt;
      
      if ( p4 ) std::cerr
        <<_fo<<"\t Address Cnt["<<name_to_hidx["Address_Cnt"]<<"] = "
        << "'"<<Address_Cnt<<"'\n"<<std::endl;
      
      Address_Idx.resize( Address_Cnt );
      Address_Avg.resize( Address_Cnt );
      Address_Sds.resize( Address_Cnt );
      Address_Rep.resize( Address_Cnt );
      
      /*
       * Load Data::
       * 
       */
      if ( p3 ) std::cerr <<_fo<<"Data Loading..." <<std::endl;
      
      _I data_cnt = offs_vec.size();
      for ( std::size_t d_idx=0; d_idx < data_cnt; d_idx++ )
      {
        // This should already be assigned!!!
        if ( name_vec[d_idx].compare("Address_Cnt")==0 &&
             Address_Cnt != 0 ) continue;
        
        int int_add_cnt = -1;
        int int_red_grn = -1;
        
        switch( type_vec[d_idx] )
        {
        case 'i' :
          if ( name_vec[d_idx].compare("Address_Cnt")==0 ) success =
            safe_load_bgz( gzin, int_add_cnt, 1,
                           name_vec[d_idx], offs_vec[d_idx], false, false, 
                           vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Red_Grn_Int")==0 ) success = 
            safe_load_bgz( gzin, int_red_grn, 1, 
                           name_vec[d_idx], offs_vec[d_idx], false, false, 
                           vb,vt+3,tc);
          break;
          
        case 'I' :
          if ( name_vec[d_idx].compare("Address_Idx")==0 ) success = 
            safe_load_bgz( gzin, Address_Idx[0], Address_Cnt, 
                           name_vec[d_idx], offs_vec[d_idx], false, false, 
                           vb,vt+3,tc);
          break;
          
        case 'H' :
          if ( name_vec[d_idx].compare("Address_Sds")==0 ) success = 
            safe_load_bgz( gzin, Address_Sds[0], Address_Cnt, 
                           name_vec[d_idx], offs_vec[d_idx], false, false, 
                           vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Address_Avg")==0 ) success = 
            safe_load_bgz( gzin, Address_Avg[0], Address_Cnt, 
                           name_vec[d_idx], offs_vec[d_idx], false, false, 
                           vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Address_Rep")==0 ) success = 
            safe_load_bgz( gzin, Address_Rep[0], Address_Cnt, 
                           name_vec[d_idx], offs_vec[d_idx], false, false, 
                           vb,vt+3,tc);
          break;
          
        case 's' :
          if ( name_vec[d_idx].compare("Sentrix_Idx")==0 ) success =
            safe_load_str_bgz( gzin, Sentrix_Idx, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Sentrix_Pos")==0 ) success =
            safe_load_str_bgz( gzin, Sentrix_Pos, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Sentrix_Man")==0 ) success =
            safe_load_str_bgz( gzin, Sentrix_Man, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Chip_Format")==0 ) success =
            safe_load_str_bgz( gzin, Chip_Format, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("OPA_Value_1")==0 ) success =
            safe_load_str_bgz( gzin, OPA_Value_1, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Description")==0 ) success =
            safe_load_str_bgz( gzin, Description, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("PlateNumber")==0 ) success =
            safe_load_str_bgz( gzin, PlateNumber, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Well_Number")==0 ) success =
            safe_load_str_bgz( gzin, Well_Number, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Unknown_001")==0 ) success =
            safe_load_str_bgz( gzin, Unknown_001, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          
          if ( name_vec[d_idx].compare("Unknown_002")==0 ) success =
            safe_load_str_bgz( gzin, Unknown_002, 1,
                               name_vec[d_idx], offs_vec[d_idx], false, true, 
                               vb,vt+3,tc);
          break;
          
        case 'S' :
          
          if ( name_vec[d_idx].compare("MidBlockIdx")==0 ||
               name_vec[d_idx].compare("RunInfo_Idx")==0 )
          {
            _I mag_len = 0;
            success =
              safe_load_bgz( gzin, mag_len, 1, name_vec[d_idx], offs_vec[d_idx],
                             false, false, vb,vt+3,tc );
            if ( !success ) return( success );
            
            const _Q offnxt = offs_vec[d_idx] + 4;
            if ( name_vec[d_idx].compare("MidBlockIdx")==0 )
            {
              MBlocks_Idx.resize( mag_len );
              success = safe_load_bgz( gzin, MBlocks_Idx[0], mag_len,
                                       name_vec[d_idx], offnxt, false, false,
                                       vb,vt+3,tc);
            }
            
            if ( name_vec[d_idx].compare("RunInfo_Idx")==0 ) {
              success = safe_load_strm_bgz( gzin, runinfo_mat, mag_len, INFO_LEN, 
                                            name_vec[d_idx], offnxt, false, 
                                            vb,vt+3,tc );
              
              if ( false ) {
                for ( std::size_t ii=0; ii<runinfo_mat.size(); ii++) {
                  if ( runinfo_mat[ii].size() != 5 ) {
                    std::cerr<<_fe<<"Run_Info["<<ii<<"].size() != 5!\n"<<std::endl;
                    return( false );
                  }
                  // Remove All White Space::
                  for ( std::size_t jj=0; jj<runinfo_mat[ii].size(); jj++)
                    std::replace( runinfo_mat[ii][jj].begin(), runinfo_mat[ii][jj].end(), ' ','_' );
                  
                  // Data Format::
                  //
                  // _S date = runinfo_mat[ii][0]; // '07/10/2022 9:15:37 PM'
                  // _S step = runinfo_mat[ii][1]; // 'Decoding | Scan | Register | Extract'
                  // _S info = runinfo_mat[ii][2]; // 'CallsToUsed=1082126|CallsToUnused=5100|CallsToInvalid=52036'
                  // _S code = runinfo_mat[ii][3]; // 'AutoDecode | iScan_Control_Software'
                  // _S vers = runinfo_mat[ii][4]; // '3.3.28'
                  
                  std::vector<_S> date_key( { "RunInfo",runinfo_mat[ii][1],runinfo_mat[ii][3],"Date" } );
                  std::vector<_S> vers_key( { "RunInfo",runinfo_mat[ii][1],runinfo_mat[ii][3],"Version" } );
                  std::vector<_S> info_key( { "RunInfo",runinfo_mat[ii][1],runinfo_mat[ii][3],"Info" } );
                  
                  runinfo_map[ date_key ].push_back( runinfo_mat[ii][0] );
                  runinfo_map[ vers_key ].push_back( runinfo_mat[ii][4] );
                  runinfo_map[ info_key ].push_back( runinfo_mat[ii][2] );
                }
              }
            }
            // if ( name_vec[d_idx].compare("RunInfo_Idx")==0 )
            // {
            //   for ( std::size_t ii=0; ii < runinfo_mat.size(); ii++) {
            //     for ( std::size_t jj=0; jj < runinfo_mat[ii].size(); jj++) {
            //       std::cerr << "\t run_vec["<<ii<<", "<<jj<<"] = '"
            //                 << runinfo_mat[ii][jj] << "'" << std::endl;
            //     }
            //     std::cerr << std::endl;
            //   }
            // }
            
          }
          break;
          
        default :
          std::cerr <<_fe<<"Invalid IDAT Code: '"<<type_vec[d_idx]<<"'\n" 
                    << std::endl;
        }
        if ( !success ) return( success );
        
        if ( Address_Cnt == 0 && int_add_cnt != -1 )
          Address_Cnt = (_I) int_add_cnt;
        
        if ( int_red_grn != -1 ) Red_Grn_Int = (_I) int_red_grn;
      }
      
    } else {
      /*
       * Only Supporting gzipped Idats for now out of complete lazyineess
       */
      std::cerr<<_fe<<"Idats must be gzipped!: '"<<idat_path<<"' "
               <<"Returning fail code...\n"<<std::endl;
      return( false );
    }
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     * 
     *       Validation/Tango Address Mapping/PnegECDF Detection P-values::
     *                      
     *  1. Quick Validation:: Address_Idx == MBlocks_Idx???
     *  2. Map Tango Addresses to Vector Indexes
     *  3. PnegECDF Background Signal & Detection P-values:: (if addresses are provided)
     *     a. Initialize PnegECDF Background Signal Data
     *     b. Load PnegECDF Background Signal Data 
     *     c. Calculate PnegECDF and Set Detection P-values 
     * 
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    // 3.a Initialize PnegECDF Background Signal Data
    const _Y calc_pval = pval_vec.size() != 0;
    std::map<_I,_I> pval_map;
    if ( calc_pval ) for ( auto& add : pval_vec ) pval_map[ add ] = 0;
    
    for ( std::size_t ii=0; ii<Address_Cnt; ii++ ) {
      // 1. Quick Validation:: Address_Idx == MBlocks_Idx???
      if ( Address_Idx[ii] != MBlocks_Idx[ii] ) {
        std::cerr<<_fe<<"Address["<<ii<<"] != MBlocks["<<ii<<"]:: "
                 <<Address_Idx[ii]<<" != "<<MBlocks_Idx[ii]<<"!"
                 << std::endl;
        return( false );
      }
      
      // 2. Map Tango Addresses to Vector Indexes
      Address_Map[ Address_Idx[ii] ] = ii;
      
      // 3.b Load PnegECDF Background Signal Data 
      if ( calc_pval && pval_map.find(Address_Idx[ii]) != pval_map.end() ) {
        Detection_Idx.push_back( Address_Idx[ii] );
        Detection_Avg.push_back( Address_Avg[ii] );
      }
    }
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     * 
     *                           Set All Variables ::
     * 
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    this->sI_dat.clear();
    this->sS_dat.clear();
    this->vI_dat.clear();
    this->vH_dat.clear();
    this->vd_dat.clear();
    
    this->add_map.clear();
    this->run_mat.clear();
    
    this->sI_dat[ "Address_Cnt" ] = Address_Cnt;
    this->sI_dat[ "Red_Grn_Int" ] = Red_Grn_Int;
    this->sI_dat[ "IdatVersion" ] = (_I) idat_version;
    // this->sI_dat[ "Header_Cnt" ]  = header_cnt;
    
    _S col_str = ""; col_str = col_str + (_c) col;
    this->sS_dat[ "Sentrix_Idx" ] = Sentrix_Idx;
    this->sS_dat[ "Sentrix_Pos" ] = Sentrix_Pos;
    this->sS_dat[ "Col_Channel" ] = col_str;
    this->sS_dat[ "Chip_Format" ] = Chip_Format;
    this->sS_dat[ "PlateNumber" ] = PlateNumber;
    this->sS_dat[ "Well_Number" ] = Well_Number;
    this->sS_dat[ "Idat_Path" ] = idat_path;
    // this->sS_dat[ "Sentrix_Man" ] = Sentrix_Man;
    // this->sS_dat[ "OPA_Value_1" ] = OPA_Value_1;
    // this->sS_dat[ "Description" ] = Description;
    // this->sS_dat[ "Unknown_001" ] = Unknown_001;
    // this->sS_dat[ "Unknown_002" ] = Unknown_002;
    
    // this->vI_dat[ "MBlocks" ] = MBlocks_Idx; // Some Idat Field, forgot the improtance???
    this->vI_dat[ "Add" ] = Address_Idx; // Tango Address
    this->vI_dat[ "Neg" ] = Detection_Idx; // Tango Address for Background (Neg Ctls)
    
    this->vH_dat[ "Sig" ] = Address_Avg; // Intensity (Signal Average)
    this->vH_dat[ "Sds" ] = Address_Sds; // Intensity Standard Deviation
    this->vH_dat[ "Rep" ] = Address_Rep; // Bead Replicate Count
    // this->vH_dat[ "Det_Sig" ] = Detection_Avg; // Intensity (Signal Average) for Background (Neg Ctls)
    
    // Map of Tango Addresses to Vector Index::
    // [Map][Add][Idx]: Tango Address to Array Index
    for ( auto& it : Address_Map ) this->add_map["Map"][ it.first ] = it.second;
    
    
    
    // 3.c Calculate PnegECDF and Set Detection P-values 
    //
    //   - Move outside and after assignment
    //   - Input Address Vector
    if ( calc_pval ) if ( !this->set_pvals( Detection_Avg,Address_Avg,Address_Det, 
                          rm_pval_outliers, vb,vt+1,tc ) ) return( false );
    
    
    // 4. Above is where we can loop over Address_Det[n]...
    
    // Detection P-values (usually form Neg Ctls)
    this->vd_dat[ "Det" ] = Address_Det;
    
    for ( std::size_t ii=0; ii < runinfo_mat.size(); ii++) {
      this->run_mat.push_back( runinfo_mat[ii] );
      // for ( std::size_t jj=0; jj < runinfo_mat[ii].size(); jj++) {
      //   std::cerr << "\t run_vec["<<ii<<", "<<jj<<"] = '"
      //             << runinfo_mat[ii][jj] << "'" << std::endl;
      // }
      // std::cerr << std::endl;
    }
    
    // Run Info:: Matrix Map
    // for ( auto& it : runinfo_map )
    //   this->data_vSvS_as.push_back( it.first, it.second );
    
    if ( p0 ) std::cerr<<_fo<<"Done. Return: '"<<success<<"'\n"<<std::endl;
    
    return( success );
  }; // parse_idat();
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                              Access Methods::
   * 
   * Scalar Methods::
   *   get_Scalar_I()
   *   get_Scalar_S()
   *   
   * Vector Methods::
   *   get_Vector_I()
   *   get_Vector_H()
   *   get_Vector_d()
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  // Get Scalar Fields::
  //
  _I get_Scalar_I( _S field ) { return( this->sI_dat[ field ] ); };
  _S get_Scalar_S( _S field ) { return( this->sS_dat[ field ] ); };
  
  // Get Scalar Data:: Failed Attempt to get one at a time...
  // _I get_Scalar_I( _I field, _I idx ) { return( this->vI_dat[ field ][idx] ); };
  // _H get_Scalar_H( _S field, _I idx ) { return( this->vH_dat[ field ][idx] ); };
  // _d get_Scalar_d( _d field, _I idx ) { return( this->vd_dat[ field ][idx] ); };
  
  // Get Vector Fields::
  //
  std::vector<_I> get_Vector_I( _S field ) { return( this->vI_dat[ field ] ); };
  std::vector<_H> get_Vector_H( _S field ) { return( this->vH_dat[ field ] ); };
  std::vector<_d> get_Vector_d( _S field ) { return( this->vd_dat[ field ] ); };
  
  // Run Info:: Matrix
  // std::vector< std::vector<_S> > get_runinfo_mat() const { return( this->runinfo_mat ); };
  
  // Detection P-value Data (Usually Negative Controls)  
  // std::vector<_I> get_Detection_Idx() const { return( this->Detection_Idx ); };
  // std::vector<_H> get_Detection_Avg() const { return( this->Detection_Avg ); };
  
  // Address Look Up::
  //
  // template <typename T>
  // _i get_Data( _I a, _S& s, T& v ) {
  //   if ( this->add_map["Map"].find(a) == this->add_map["Map"].end() )
  //     return( -1 );
  //   // std::size_t idx = this->add_map["Map"][a];
  //   if ( a == 0 ) return( 1 );
  //   
  //   if ( this->vI_dat.find(s) != this->vI_dat.end() ) {
  //     v = (T) this->vI_dat[ s ][ this->add_map["Map"][ a ] ];
  //     return( 1 );
  //   } else if ( this->vH_dat.find(s) != this->vH_dat.end() ) {
  //     v = (T) this->vH_dat[ s ][ this->add_map["Map"][ a ] ];
  //     return( 1 );
  //   } else if ( this->vd_dat.find(s) != this->vd_dat.end() ) {
  //     v = (T) this->vd_dat[ s ][ this->add_map["Map"][ a ] ];
  //     return( 1 );
  //   } else {
  //     return( -2 );
  //   }
  //   return( -3 );
  // };
  
  template <typename T>
  _Y get_Select_Data_Pair( const _S& field,
                           const std::pair<_I,_I>& pair,
                           std::pair<_S,_S>& keyU,
                           std::pair<_S,_S>& keyM,
                           std::pair< T, T>& data,
                           const _I vb = 0, const _I vt = 1, const _I tc = 0,
                           const _S ft="get_Select_Data_Pair" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    const _Y p8 = vb > vt + 8;
    
    _Y success = true;
    keyU.first = "";    keyU.second = "";
    keyM.first = "";    keyM.second = "";
    data.first = (T) 0; data.second = (T) 0;
    
    if ( this->vI_dat.find(field) == this->vI_dat.end() &&
         this->vH_dat.find(field) == this->vH_dat.end() && 
         this->vd_dat.find(field) == this->vd_dat.end() &&
         this->add_map.find("Map") == this->add_map.end() &&
         this->sS_dat.find("Col_Channel") == this->sS_dat.end() ) {
      std::cerr<<_fe<<"Field = '"<<field<<"' is NOT initialized!"<<std::endl;
      return( false );
    }
    
    // Set Infinium Design and Color Channels::
    keyU.first = field; keyU.second = "U"+this->sS_dat[ "Col_Channel" ];
    keyM.first = field; keyM.second = "M"+this->sS_dat[ "Col_Channel" ];
    std::size_t map_cnt = this->add_map["Map"].size();
    
    if ( p0 ) std::cerr
      <<_fo<<"Starting...\n"
      <<_fo<<"\t    pair = ["<<pair.first<<", "<<pair.second<<"]\n"
      <<_fo<<"\t    keyU = ["<<keyU.first<<", "<<keyU.second<<"]\n"
      <<_fo<<"\t    keyM = ["<<keyM.first<<", "<<keyM.second<<"]\n"
      <<_fo<<"\t map_cnt = '"<<map_cnt<<"'\n"
      <<std::endl;
    
    if ( map_cnt == 0 ) {
      std::cerr<<_fe<<"Invalid Address Map Initialization Sizes!"<<std::endl;
      return( false );
    }
    
    /*
     * Check UnMethylated/Infinium II Address (A)
     * 
     */
    if ( success && pair.first != 0 ) {
      if ( success && this->add_map["Map"].find(pair.first) != this->add_map["Map"].end() ) {
        if ( success && this->vI_dat.find( field ) != this->vI_dat.end() ) {
          data.first = (T) this->vI_dat[ field ][ this->add_map["Map"][ pair.first ] ];
        } else if ( success && this->vH_dat.find( field ) != this->vH_dat.end() ) {
          data.first = (T) this->vH_dat[ field ][ this->add_map["Map"][ pair.first ] ];
        } else if ( success && this->vd_dat.find( field ) != this->vd_dat.end() ) {
          data.first = (T) this->vd_dat[ field ][ this->add_map["Map"][ pair.first ] ];
        } else {
          if ( p8 ) std::cerr
            <<_fe<<"VAL.1["<<pair.first<<"', "<<pair.second<<"]\n"<<std::endl;
          success = false;
        }
      } else {
        if ( p8 ) std::cerr
          <<_fe<<"MAP.1["<<pair.first<<"', "<<pair.second<<"]\n"<<std::endl;
        success = false;
      }
    }
    
    /*
     * Check Methylated Infinium I Address (B)
     * 
     */
    if ( success && pair.second != 0 ) {
      if ( success && this->add_map["Map"].find(pair.second) != this->add_map["Map"].end() ) {
        if ( success && this->vI_dat.find( field ) != this->vI_dat.end() ) {
          data.second = (T) this->vI_dat[ field ][ this->add_map["Map"][ pair.second ] ];
        } else if ( success && this->vH_dat.find( field ) != this->vH_dat.end() ) {
          data.second = (T) this->vH_dat[ field ][ this->add_map["Map"][ pair.second ] ];
        } else if ( success && this->vd_dat.find( field ) != this->vd_dat.end() ) {
          data.second = (T) this->vd_dat[ field ][ this->add_map["Map"][ pair.second ] ];
        } else {
          if ( p8 ) std::cerr
            <<_fe<<"VAL.2["<<pair.first<<"', "<<pair.second<<"]\n"<<std::endl;
          success = false;
        }
      } else {
        if ( p8 ) std::cerr
          <<_fe<<"MAP.2["<<pair.first<<"', "<<pair.second<<"]\n"<<std::endl;
        success = false;
      }
    }
    if ( !success ) {
      keyU.first = ""; keyU.second = ""; data.first  = (T) 0;
      keyM.first = ""; keyM.second = ""; data.second = (T) 0;
    }
    
    return( success );
  }; // get_Select_Data_Pair();
  
  template <typename T>
  _Y get_Select_Data_Pairs( _S& field,
                            std::vector< std::pair<_I,_I> >& par,
                            std::map< std::pair<_S,_S>, std::vector<T> >& dat,
                            const _I vb = 0, const _I vt = 1, const _I tc = 0,
                            const _S ft="get_Select_Data_Pairs" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    
    _Y success = true;
    dat.clear();
    
    // QC:: Validate Query::
    if ( field.compare("Sig") == 0 && field.compare("Sds") == 0 &&
         field.compare("Rep") == 0 && field.compare("Add") == 0 &&
         field.compare("Neg") == 0 && field.compare("Det") == 0 ) {
      std::cerr<<_fe<<"Field = '"<<field<<"' Not support for Idats!"<<std::endl;
      return( false );
    }
    if ( this->vI_dat.find(field) == this->vI_dat.end() &&
         this->vH_dat.find(field) == this->vH_dat.end() &&
         this->vd_dat.find(field) == this->vd_dat.end() ) {
      std::cerr<<_fe<<"Field = '"<<field<<"' is NOT initialized!"<<std::endl;
      return( false );
    }
    
    // Set Infinium Design and Color Channels::
    std::pair<_S,_S> pair_U;
    pair_U.first  = field;
    pair_U.second = "U"+this->sS_dat[ "Col_Channel" ];
    
    std::pair<_S,_S> pair_M;
    pair_M.first  = field;
    pair_M.second = "M"+this->sS_dat[ "Col_Channel" ];
    
    // std::size_t col_vec_cnt = col_vec.size();
    std::size_t par_cnt = par.size();
    std::size_t map_cnt = this->add_map["Map"].size();
    if ( p0 ) std::cerr
      <<_fo<<"Starting...\n"
      <<_fo<<"\t   field = '"<<field<<"'\n"
      <<_fo<<"\t  pair_U = '"<<pair_U.second<<"'\n"
      <<_fo<<"\t  pair_M = '"<<pair_M.second<<"'\n"
      <<_fo<<"\t par_cnt = '"<<par_cnt<<"'\n"
      <<_fo<<"\t map_cnt = '"<<map_cnt<<"'\n"
      <<std::endl;
    
    if ( par_cnt == 0 || map_cnt == 0 ) {
      std::cerr<<_fe<<"Invalid Input/Initialization Sizes!"<<std::endl;
      return( false );
    }
    
    for ( std::size_t idx=0; idx<par_cnt; idx++ ) {
      T dU = (T) 0; T dM = (T) 0;
      // if ( idx < 10 )
      //   std::cerr<<_fe<<"Attempting(First):: idx='"<<idx<<"': first='"<<par[idx].first<<"', second='"<<par[idx].second<<"'"<<std::endl;
      
      /*
       * Check UnMethylated/Infinium II Address (A)
       *
       */
      if ( success && this->add_map["Map"].find(par[idx].first) != this->add_map["Map"].end() ) {
        if ( success && this->vI_dat.find( field ) != this->vI_dat.end() ) {
          dU = (T)  this->vI_dat[ field ][ this->add_map["Map"][ par[idx].first ] ];
        } else if ( success && this->vH_dat.find( field ) != this->vH_dat.end() ) {
          dU = (T)  this->vH_dat[ field ][ this->add_map["Map"][ par[idx].first ] ];
        } else if ( success && this->vd_dat.find( field ) != this->vd_dat.end() ) {
          dU = (T)  this->vd_dat[ field ][ this->add_map["Map"][ par[idx].first ] ];
        } else {
          std::cerr<<_fe<<"Failed to find first VAL:: idx='"<<idx<<"': "
                   <<"first='"<<par[idx].first<<"', "
                   <<"second='"<<par[idx].second<<"'\n"<<std::endl;
          success = false;
        }
      } else {
        std::cerr<<_fe<<"Failed to find first MAP:: idx='"<<idx<<"': "
                 <<"first='"<<par[idx].first<<"', "
                 <<"second='"<<par[idx].second<<"'\n"<<std::endl;
        success = false;
      }
      
      /*
       * Check Methylated Infinium I Address (B)
       *
       */
      if ( par[idx].second != 0 ) {
        // if ( idx < 10 )
        //   std::cerr<<_fe<<"Attempting(Second):: idx='"<<idx<<"': first='"<<par[idx].first<<"', second='"<<par[idx].second<<"'"<<std::endl;
        
        if ( this->add_map["Map"].find(par[idx].second) != this->add_map["Map"].end() ) {
          if ( success && this->vI_dat.find( field ) != this->vI_dat.end() ) {
            dM = (T) this->vI_dat[ field ][ this->add_map["Map"][ par[idx].second ] ];
          } else if ( success && this->vH_dat.find( field ) != this->vH_dat.end() ) {
            dM = (T) this->vH_dat[ field ][ this->add_map["Map"][ par[idx].second ] ];
          } else if ( success && this->vd_dat.find( field ) != this->vd_dat.end() ) {
            dM = (T) this->vd_dat[ field ][ this->add_map["Map"][ par[idx].second ] ];
          } else {
            std::cerr<<_fe<<"Failed to find second VAL:: idx='"<<idx<<"': "
                     <<"first='"<<par[idx].first<<"', "
                     <<"second='"<<par[idx].second<<"'\n"<<std::endl;
            success = false;
          }
        } else {
          std::cerr<<_fe<<"Failed to find second MAP:: idx='"<<idx<<"': "
                   <<"first='"<<par[idx].first<<"', "
                   <<"second='"<<par[idx].second<<"'\n"<<std::endl;
          success = false;
        }
      }
      if ( !success ) break;
      
      dat[pair_U].push_back( dU );
      dat[pair_M].push_back( dM );
    }
    if ( par_cnt != dat[pair_U].size() || par_cnt != dat[pair_U].size() )
      success = false;
    if ( !success ) std::cerr<<_fe<<"Failed to match sizes!"<<std::endl;
    if ( !success ) dat.clear();
    
    return( success );
  }; // get_Select_Data_Pairs();
  
  template <typename T>
  _Y get_Select_Data_Vec( _S& field,
                          std::vector<_I>& add,
                          std::vector<T>&  dat,
                          
                          const _I vb = 0,
                          const _I vt = 1,
                          const _I tc = 0,
                          const _S ft="get_Select_Data_Vec" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    dat.clear();
    
    if ( field.compare("Sig") == 0 && 
         field.compare("Sds") == 0 && 
         field.compare("Rep") == 0 && 
         field.compare("Add") == 0 && 
         field.compare("Neg") == 0 && 
         field.compare("Det") == 0 ) {
      std::cerr<<_fe<<"Field = '"<<field<<"' Not support for Idats!"<<std::endl;
      return(false);
    }
    
    if ( this->vI_dat.find(field) == this->vI_dat.end() &&
         this->vH_dat.find(field) == this->vH_dat.end() && 
         this->vd_dat.find(field) == this->vd_dat.end() ) {
      std::cerr<<_fe<<"Field = '"<<field<<"' is NOT initialized!"<<std::endl;
      return(false);
    }
    
    _I add_map_cnt = this->add_map["Map"].size();
    if ( p0 ) std::cerr
      <<_fo<<"Starting...\n"
      <<_fo<<"\t      field  = '"<<field<<"'\n"
      <<_fo<<"\t add_map_cnt = '"<<add_map_cnt<<"'"
      <<std::endl;
    if ( add_map_cnt == 0 ) return( false );
    
    _I add_size = add.size();
    for ( std::size_t ii=0; ii<add_size; ii++ ) {
      T d = 0;
      if ( add[ii] != 0 && this->add_map["Map"].find(add[ii]) != this->add_map["Map"].end() ) {
        if ( this->vI_dat.find(field) != this->vI_dat.end() ) {
          d = (T) this->vI_dat[ field ][ this->add_map["Map"][ add[ii] ] ];
        } else if ( this->vH_dat.find(field) != this->vH_dat.end() ) {
          d = (T) this->vH_dat[ field ][ this->add_map["Map"][ add[ii] ] ];
        } else if ( this->vd_dat.find(field) != this->vd_dat.end() ) {
          d = (T) this->vd_dat[ field ][ this->add_map["Map"][ add[ii] ] ];
        } else {
          return( false );
        }
      }
      dat.push_back( d );
    }
    if ( add_size != dat.size() ) 
      std::cerr<<_fe<<"Failed to match sizes!"<<std::endl;
    
    if ( add_size == dat.size() ) return( true );
    return( false );
  }; // get_Select_Data_Vec();
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                           Manipulation Methods::
   * 
   *   set_pvals()
   *   reset_pvals()
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  _Y set_pvals( const std::vector<_H>& d_vec,
                const std::vector<_H>& s_vec,
                std::vector<_d>& p_vec,
                const _Y rm_outliers = false,
                
                const _I vb = 0,
                const _I vt = 4,
                const _I tc = 0,
                const _S ft="set_pvals" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    const _Y p4 = vb > vt + 4;
    
    p_vec.clear();
    
    if ( p0 ) std::cerr 
      <<_fo<<"Starting...\n"
      <<_fo<<"\t  (vb,vt,tc): (" <<vb<<","<<vt<<","<<tc<<")\n"
      <<_fo<<"\t d_vec.size() = '"<<d_vec.size()<<"'\n" 
      <<_fo<<"\t s_vec.size() = '"<<s_vec.size()<<"'\n" 
      <<_fo<<"\t p_vec.size() = '"<<p_vec.size()<<"'\n"
      <<_fo<<"\t  rm_outliers = '"<<rm_outliers<<"'"
      << std::endl;
    
    /*
     * Calculate ecdf and calculate p-values::
     *   TBD:: Remove Outliers
     *      - mean + n*std 
     *      - recalculated p-value with its own distribution
     *   TBD:: Generate Plot
     */
    
    std::vector<_H> b_vec( d_vec );
    std::sort( b_vec.begin(), b_vec.end() );
    auto ecdf = empirical_cumulative_distribution_function( std::move(b_vec) );
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                          Remove/Calculate Outliers::
     *                          
     *  NOTE:: So far this is not recomended...
     *  TBD:: Record Number of Outliers as a Class Field. This will help 
     *    determine if background recalculation adds a bunch of junk...
     *    
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    if ( rm_outliers ) {
      _d outlier_cut = 0.2;
      std::vector<_H> n_vec;
      std::vector<_H> o_vec( d_vec );
      std::sort( o_vec.begin(), o_vec.end() );
      
      _I o_cnt = o_vec.size();
      // _H o_med = 0;
      // _d o_std = 0.0;
      
      _H o_med = median_vec( o_vec );
      _d o_avg = mean_vec( o_vec );
      _d o_std = stdev_vec( o_vec );
      // _d o_std = stdev_vec( o_vec, o_avg );
      _d o_cut = o_avg + ( 1*o_std );
      _I o_inc = int( o_vec.size() / 20 );
      
      if ( p4 ) {
        std::cerr<<_fo<<"\n\nBEG OUTLIER TESTING:: "
                 <<"o_cnt="<<o_cnt<<", "
                 <<"o_med="<<o_med<<", "
                 <<"o_avg="<<o_avg<<", "
                 <<"o_std="<<o_std<<", "
                 <<"o_cut="<<o_cut
                 <<std::endl;
        for ( std::size_t ii=0; ii<o_vec.size(); ii += o_inc ) {
          _d det = 1.0 - ((_d) ecdf(o_vec[ii]) );
          _Y o_pass = o_vec[ii] < o_cut; 
          std::cerr<<_fo<<"\t sig["<<ii<<"]="<<o_vec[ii]
                   <<", det="<<det
                   <<", o_pass="<<o_pass<<"!"
                   <<std::endl;
        }
        std::cerr<<_fo<<"END OUTLIER TESTING!\n"<<std::endl;
      }
      for ( std::size_t ii=0; ii<o_vec.size(); ii++ ) {
        _d det = 1.0 - ((_d) ecdf(o_vec[ii]) );
        if ( det >= outlier_cut ) n_vec.push_back( o_vec[ii] );
      }
      ecdf = empirical_cumulative_distribution_function(std::move(n_vec));
      
      this->sI_dat[ "Detection_Cnt0" ] = s_vec.size();
      this->sI_dat[ "Outlier_Cnt0" ] = n_vec.size();
    }
    
    _I s_len = s_vec.size();
    for ( std::size_t ii=0; ii<s_len; ii++ )
      p_vec.push_back( 1.0 - ( (_d) ecdf(s_vec[ii]) ) );
    
    if ( p0 ) std::cerr<<_fo<<"Done. Pvals.size(): '"<<p_vec.size()<<"'\n"<<std::endl;
    
    return( true );
  }; // set_pvals();
  
  _Y reset_pvals( const std::vector<_I>& add_vec,
                  const _Y rm_outliers = false,
                  
                  const _I vb = 0,
                  const _I vt = 4,
                  const _I tc = 0,
                  const _S ft="reset_pvals" )
  {
    const _S tb(tc, '\t');
    const _S _fo("["+ft+"]: "+tb);
    const _S _fe("["+ft+"]: ERROR: ");
    const _S _fw("["+ft+"]: "+tb+"Warning: ");
    
    const _Y p0 = vb > vt + 0;
    const _Y p1 = vb > vt + 1;
    const _Y p4 = vb > vt + 4;
    
    if ( p0 ) std::cerr 
      <<_fo<<"Starting...\n"
      <<_fo<<"\t  (vb,vt,tc): (" <<vb<<","<<vt<<","<<tc<<")\n"
      <<_fo<<"\t add_vec.size() = '"<<add_vec.size()<<"'\n" 
      <<_fo<<"\t add_map.size() = '"<<this->add_map["Map"].size()<<"'\n" 
      <<_fo<<"\t    rm_outliers = '"<<rm_outliers<<"'"
      << std::endl;

    this->vd_dat[ "Det" ].clear(); // p_vec.clear();
    
    // Used to remove and count replicate addresses
    std::map< _I,_I > add_cnt_map;
    std::vector< _H > neg_sig_vec;
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                         Re-extract Raw Intensities::
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    _I add_mat_cnt  = 0;
    _I add_mis_cnt  = 0;
    _I add_dup_cnt  = 0;
    _I add_vec_size = add_vec.size();
    for ( std::size_t ii=0; ii<add_vec_size; ii++ ) {
      add_cnt_map[add_vec[ii]]++;
      
      if ( add_cnt_map[add_vec[ii]] > 1 ) {
        add_dup_cnt++;
      } else {
        if ( this->add_map["Map"].find(add_vec[ii]) == this->add_map["Map"].end() ) {
          // std::cerr 
          // <<_fe<<"Failed to find Address["<<ii<<"]='"<<add_vec[ii]<<"'\n"<< std::endl;
          add_mis_cnt++;
        } else {
          _I add_idx = this->add_map["Map"][add_vec[ii]];
          _H add_sig = this->vH_dat[ "Sig" ][add_idx];
          neg_sig_vec.push_back( add_sig );
          add_mat_cnt++;
        }
      }
    }
    
    if ( p4 && add_mis_cnt > 0 ) std::cerr
      <<_fw<<"Address Miss Count = '"<<add_mis_cnt<<"' > 0!!!\n"<<std::endl;

    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                                Build ECDF::
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    std::vector<_H> b_vec( neg_sig_vec );
    std::sort( b_vec.begin(), b_vec.end() );
    auto ecdf = empirical_cumulative_distribution_function( std::move(b_vec) );

    if ( p1 ) std::cerr 
      <<_fo<<"Intermediate Stats::\n"
      <<_fo<<"\t      add_dup_cnt = '"<<add_dup_cnt<<"'\n" 
      <<_fo<<"\t      add_mat_cnt = '"<<add_mat_cnt<<"'\n" 
      <<_fo<<"\t      add_mis_cnt = '"<<add_mis_cnt<<"'\n" 
      <<_fo<<"\t neg_sig_vec.size = '"<<neg_sig_vec.size()<<"'\n" 
      << std::endl;
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                          Remove/Calculate Outliers::
     *                          
     *  NOTE:: So far this is not recomended...
     *  TBD:: Record Number of Outliers as a Class Field. This will help 
     *    determine if background recalculation adds a bunch of junk...
     *    
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    if ( rm_outliers ) {
      _d n_std = 1.0;
      _d outlier_cut = 0.2;
      std::vector<_H> n_vec;
      std::vector<_H> o_vec( neg_sig_vec );
      std::sort( o_vec.begin(), o_vec.end() );
      
      _I o_cnt = o_vec.size();

      _H o_med = median_vec( o_vec );
      _d o_avg = mean_vec( o_vec );
      _d o_std = stdev_vec( o_vec );
      
      _d o_cut = o_avg + ( n_std * o_std );
      _I o_inc = int( o_vec.size() / 20 );
      
      if ( p4 ) {
        std::cerr<<_fo<<"\n\nBEG OUTLIER TESTING:: "
                 <<"o_cnt="<<o_cnt<<", "
                 <<"o_med="<<o_med<<", "
                 <<"o_avg="<<o_avg<<", "
                 <<"o_std="<<o_std<<", "
                 <<"o_cut="<<o_cut
                 <<std::endl;
        for ( std::size_t ii=0; ii<o_vec.size(); ii += o_inc ) {
          _d det = 1.0 - ((_d) ecdf(o_vec[ii]) );
          _Y o_pass = o_vec[ii] < o_cut;
          
          if ( p4 ) std::cerr
            <<_fo<<"\t sig["<<ii<<"]="<<o_vec[ii]
            <<", det="<<det
            <<", o_pass="<<o_pass<<"!"
            <<std::endl;
        }
        if ( p4 ) std::cerr<<_fo<<"END OUTLIER TESTING!\n"<<std::endl;
      }
      for ( std::size_t ii=0; ii<o_vec.size(); ii++ ) {
        _d det = 1.0 - ((_d) ecdf(o_vec[ii]) );
        if ( det >= outlier_cut ) n_vec.push_back( o_vec[ii] );
      }

      this->sI_dat[ "Detection_Cnt1" ] = neg_sig_vec.size();
      this->sI_dat[ "Outlier_Cnt1" ] = n_vec.size();
      
      if ( p1 ) std::cerr
        <<_fo<<"Outlier Reduction Counts:: "
        <<"n_cnt="<<n_vec.size()<<", "
        <<"o_cnt="<<neg_sig_vec.size()<<", "
        <<std::endl;
      ecdf = empirical_cumulative_distribution_function(std::move(n_vec));
    }
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                          Re-assign New DetP Values::
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    _I sigs_len = this->vH_dat[ "Sig" ].size();
    for ( std::size_t ii=0; ii<sigs_len; ii++ )
      this->vd_dat[ "Det" ].push_back( 1.0 - ( (_d) ecdf(this->vH_dat[ "Sig" ][ii]) ) );

    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     *                                    Done::
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    if ( p0 ) std::cerr
      <<_fo<<"Done. DetP.size(): '"<<this->vd_dat[ "Det" ].size()<<"'\n"
      <<std::endl;
    
    return( true );
  }; // reset_pvals()
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                               Output Methods::
   * 
   *   to_string()
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  _S to_string( const _S format = "tab",
                const _S prefix = "[Idat]: ",
                const _B sep = ',',
                
                const _I vb = 0,
                const _I vt = 1,
                const _I tc = 0,
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
      // Basic Scalar Pringing::
      //   for ( auto& it : sI_dat ) ss <<prefix<<it.first<<"='"<<it.second<<"'\n";
      //   for ( auto& it : sS_dat ) ss <<prefix<<it.first<<"='"<<it.second<<"'\n";
      
      // Map of Key(field)/Value pairs for scalars::
      // [Key][Val]: Scalar Fields
      for ( auto& it : sI_dat )
        ss<<prefix<<it.first<<" = '"<<it.second<<"'"<<std::endl;
      ss<<std::endl;
      
      for ( auto& it : sS_dat )
        ss<<prefix<<it.first<<" = '"<<it.second<<"'"<<std::endl;
      ss<<std::endl;
      
      
      // for ( auto& it : sI_dat )
      //   ss<<prefix<<it.first<<" = '"<<stringify_sum(it.second)<<"'"<<std::endl;
      // ss<<std::endl;
      
      // for ( auto& it : sS_dat )
      //   ss<<prefix<<it.first<<" = '"<<stringify_sum(it.second)<<"'"<<std::endl;
      // ss<<std::endl;
      
      // Vector of Tango Addresses by class {All,Neg}
      // [Add]: All Tango Addresses
      // [Neg]: Neg Tango Addresses
      for ( auto& it : vI_dat )
        ss<<prefix<<it.first<<" = '"<<stringify_sum(it.second)<<"'"<<std::endl;
      ss<<std::endl;
      
      // Vectors of data[_H] {Sig,Sds,Rep} in order of Tango Address
      // [Sig]: All Signal Intesnsity, 
      // [Sds]: All Signal Standard Deviation, 
      // [Rep]: All Bead Replicate Counts
      for ( auto& it : vH_dat )
        ss<<prefix<<it.first<<" = '"<<stringify_sum(it.second)<<"'"<<std::endl;
      ss<<std::endl;
      
      // // Vectors of data[_d] {Sig,Sds,Rep} in order of Tango Address
      // // [Det]: All Detection P-values
      for ( auto& it : vd_dat )
        ss<<prefix<<it.first<<" = '"<<stringify_sum(it.second)<<"'"<<std::endl;
      ss<<std::endl;
      
      // Only Missing::
      // Matrix[_S][_S] of Run Info::
      //  std::vector< std::vector<_S> > run_mat;
      
      //
      // Example of Stringify Methods::
      //
      // if ( f.compare("sum")==0 ) {
      //   for( size_t ii=0; ii<sz; ii++ ) {
      //     ss <<pre
      //        <<"'"<<stringify_sum(idx_to_key[ii])<<"'"
      //        <<"["<<stringify_sum(key_to_idx[idx_to_key[ii]])<<"] = "
      //        <<"'"<<stringify_sum(idx_to_dat[ii])<<"'"
      //        << RET;
      //   }
      // } else if ( f.compare("all")==0 ) {
      //   for( size_t ii=0; ii<sz; ii++ ) {
      //     ss <<pre
      //        <<"'"<<stringify_all(idx_to_key[ii])<<"'"
      //        <<"["<<stringify_all(key_to_idx[idx_to_key[ii]])<<"] = "
      //        <<"'"<<stringify_all(idx_to_dat[ii])<<"'"
      //        << RET;
      //   }
      // }
    }
    
    return( ss.str() );
  };

};

};

#endif /* idat.h */

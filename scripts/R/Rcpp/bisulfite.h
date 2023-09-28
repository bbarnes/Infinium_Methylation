#ifndef __BISULFITE_H__
#define __BISULFITE_H__

#include "stl_templates.h"
#include "fast_ostream.h"

// using namespace std;
using namespace stl_templates;

namespace bisulfite {

// Idat Pair Fields::
_SC _S SIG_STR("Sig");
_SC _S DET_STR("Det");

_SC _S PVAL_STR("Pval");
_SC _S BETA_STR("Beta");

_SC _S UG_STR("UG");
_SC _S UR_STR("UR");
_SC _S MG_STR("MG");
_SC _S MR_STR("MR");

// Sample Manifest Fields::
_SC _S CGN_STR = "CGN";
_SC _S PRB_STR = "Prb";
_SC _S COL_STR = "col";
_SC _S COL_REF = "col_ref";

_SC _S ADD_STR = "Add";
_SC _S MAP_STR = "Map";
_SC _S ALL_STR = "All";
_SC _S WRK_STR = "Workflow";
_SC _S PRG_STR = "Progress";
_SC _S MASK_STR = "mask";

_SC _S MIN_PVAL_STR = "Min_Pval";
_SC _S MIN_PERI_STR = "Min_Perc_Inside";
_SC _S MIN_PERO_STR = "Min_Pval_Outside";

_SC _S MIN_BETA_STR = "Min_Beta";
_SC _S MAX_BETA_STR = "Max_Beta";

// Sample Manifest Stats Fields::
_SC _S TOT_CNT_STR  = "Total_Count";
_SC _S PASS_CNT_STR = "Pass_Count";

_SC _S SENTRIX_NAME = "Sentrix_Name";
_SC _S SENTRIX_ID   = "Sentrix_Idx";
_SC _S SENTRIX_POS  = "Sentrix_Pos";

_SC _S CHIP_FORMAT  = "Chip_Format";
_SC _S ADDRESS_CNT  = "Address_Cnt";
_SC _S IDATVERSION  = "IdatVersion";

_SC _S PROBE_ID   = "Cgn";
_SC _S MANIFEST   = "Man";
_SC _S CONTROLS   = "Ctl";
_SC _S COLOR_STR  = "Col";
_SC _S PROBE_DIN  = "Din";
_SC _S CHROMOSOME = "Chr";
_SC _S COORDINATE = "Pos";
_SC _S ANNOTATION = "Ann";

_SC _S LOCUS_CLUSTER = "Loc";
_SC _S PROBE_CLUSTER = "Prb";

_SC _d dNEG_ONE = -1.0;

// Used for setting default empty vectors in functions
static std::vector<_S> DEFAULT_STRING_VEC;
static std::vector<_I> DEFAULT_UINT_VEC;

enum
{
  BASE_A = 0x0, /* binary: 00 */
  BASE_C = 0x1, /* binary: 01 */
  BASE_G = 0x2, /* binary: 10 */
  BASE_T = 0x3, /* binary: 11 */
};

static _c comp_vec[256];
static _c bscX_vec[256];
static _c bscU_vec[256];

static _c bscM1_vec[256][256][2];
static _c bscM2_vec[256];
static _c bscD1_vec[256][256][2];
static _c bscD2_vec[256];
static _i  topbot_vec[256][256];

static _c iupac_exp_len[256];
static std::map< _c, std::vector<_c> > iupac_exp_nuc;

template <typename T>
struct Static_cast {
  template <typename U>
  T operator () (const U& x) const { return T(x); }
};

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *   Complement and Bisulfite Conversion Matrix Mapping Initializaiton::
 *
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static _Y init_bscX = false;
static _Y init_bscU = false;
static _Y init_bscM = false;
static _Y init_bscD = false;
static _Y init_rvcp = false;
static _Y init_topbot = false;
static _Y init_exp_iupac = false;

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                         Reverse Complement Table::
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static void iupac_exp_init() {
  
  // iupac_exp_len, iupac_exp_nuc, init_exp_iupac
  
  for (_I ii = 0; ii < 256; ii++)
    iupac_exp_len[ii] = 1;
  
  iupac_exp_len[ (_B) 'R' ] = 2;
  iupac_exp_len[ (_B) 'r' ] = 2;
  
  iupac_exp_len[ (_B) 'Y' ] = 2;
  iupac_exp_len[ (_B) 'y' ] = 2;
  
  iupac_exp_len[ (_B) 'S' ] = 2;
  iupac_exp_len[ (_B) 's' ] = 2;
  
  iupac_exp_len[ (_B) 'W' ] = 2;
  iupac_exp_len[ (_B) 'w' ] = 2;
  
  iupac_exp_len[ (_B) 'K' ] = 2;
  iupac_exp_len[ (_B) 'k' ] = 2;
  
  iupac_exp_len[ (_B) 'M' ] = 2;
  iupac_exp_len[ (_B) 'm' ] = 2;
  
  iupac_exp_len[ (_B) 'B' ] = 3;
  iupac_exp_len[ (_B) 'b' ] = 3;
  
  iupac_exp_len[ (_B) 'D' ] = 3;
  iupac_exp_len[ (_B) 'd' ] = 3;
  
  iupac_exp_len[ (_B) 'H' ] = 3;
  iupac_exp_len[ (_B) 'h' ] = 3;
  
  iupac_exp_len[ (_B) 'V' ] = 3;
  iupac_exp_len[ (_B) 'v' ] = 3;
  
  std::vector<_c> c_vec;
  c_vec.push_back('A');
  c_vec.push_back('G');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'R', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'r', c_vec) );
  c_vec.clear();
  
  c_vec.push_back('C');
  c_vec.push_back('T');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'Y', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'y', c_vec) );
  c_vec.clear();
  
  // S
  c_vec.push_back('C');
  c_vec.push_back('G');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'S', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 's', c_vec) );
  c_vec.clear();
  
  // W
  c_vec.push_back('A');
  c_vec.push_back('T');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'W', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'w', c_vec) );
  c_vec.clear();
  
  // K
  c_vec.push_back('G');
  c_vec.push_back('T');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'K', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'k', c_vec) );
  c_vec.clear();
  
  // M
  c_vec.push_back('A');
  c_vec.push_back('C');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'M', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'm', c_vec) );
  c_vec.clear();
  
  // B
  c_vec.push_back('C');
  c_vec.push_back('G');
  c_vec.push_back('T');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'B', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'b', c_vec) );
  c_vec.clear();
  
  // D
  c_vec.push_back('A');
  c_vec.push_back('G');
  c_vec.push_back('T');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'D', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'd', c_vec) );
  c_vec.clear();
  
  // H
  c_vec.push_back('A');
  c_vec.push_back('C');
  c_vec.push_back('T');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'H', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'h', c_vec) );
  c_vec.clear();
  
  // V
  c_vec.push_back('A');
  c_vec.push_back('C');
  c_vec.push_back('G');
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'V', c_vec) );
  iupac_exp_nuc.insert( pair< _c, std::vector<_c> >( 'v', c_vec) );
  c_vec.clear();
  
  init_exp_iupac = true;
}


/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                         Reverse Complement Table::
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static void revcomp_init() {
  
  for (_I ii = 0; ii < 256; ii++)
    comp_vec[ii] = (_B) ii;
  
  comp_vec[ (_B) 'A'] = 'T'; comp_vec[ (_B) 'a'] = 't';
  comp_vec[ (_B) 'C'] = 'G'; comp_vec[ (_B) 'c'] = 'g';
  comp_vec[ (_B) 'G'] = 'C'; comp_vec[ (_B) 'g'] = 'c';
  comp_vec[ (_B) 'T'] = 'A'; comp_vec[ (_B) 't'] = 'a';
  comp_vec[ (_B) 'U'] = 'A'; comp_vec[ (_B) 'u'] = 'a';
  comp_vec[ (_B) 'M'] = 'K'; comp_vec[ (_B) 'm'] = 'k';
  comp_vec[ (_B) 'R'] = 'Y'; comp_vec[ (_B) 'r'] = 'y';
  comp_vec[ (_B) 'W'] = 'W'; comp_vec[ (_B) 'w'] = 'w';
  comp_vec[ (_B) 'S'] = 'S'; comp_vec[ (_B) 's'] = 's';
  comp_vec[ (_B) 'Y'] = 'R'; comp_vec[ (_B) 'y'] = 'r';
  comp_vec[ (_B) 'K'] = 'M'; comp_vec[ (_B) 'k'] = 'm';
  comp_vec[ (_B) 'V'] = 'B'; comp_vec[ (_B) 'v'] = 'b';
  comp_vec[ (_B) 'H'] = 'D'; comp_vec[ (_B) 'h'] = 'd';
  comp_vec[ (_B) 'D'] = 'H'; comp_vec[ (_B) 'd'] = 'h';
  comp_vec[ (_B) 'B'] = 'V'; comp_vec[ (_B) 'b'] = 'v';
  comp_vec[ (_B) 'N'] = 'N'; comp_vec[ (_B) 'n'] = 'n';
  
  // Other Special Characters To Swap Around::
  comp_vec[ (_B) '['] = ']'; comp_vec[ (_B) ']'] = '[';
  
  init_rvcp = true;
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                      De-Methylate Table:: Used for Alignment
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static void bscX_init() {
  
  /*
   * 
   * Should loop through entire ASCII space to map in one-to-one for most
   *   ASCII characters where this is meaningless, or there is no covnersion.
   *   This prevents null pointer exections when non-converted nucleotides
   *   are passed in. After we then assign the true conversions.
   *   
   *           x <- tr(x, 'RYSWKMBDHV', 'ATCATATAAA')
   *   
   */
  
  for (_I ii = 0; ii < 256; ii++)
    bscX_vec[ii] = (_B) ii;
  
  bscX_vec[ (_B) 'R'] = 'a'; bscX_vec[ (_B) 'r'] = 'a';
  bscX_vec[ (_B) 'Y'] = 't'; bscX_vec[ (_B) 'y'] = 't';
  bscX_vec[ (_B) 'S'] = 'c'; bscX_vec[ (_B) 's'] = 'c';
  bscX_vec[ (_B) 'W'] = 'a'; bscX_vec[ (_B) 'w'] = 'a';
  bscX_vec[ (_B) 'K'] = 't'; bscX_vec[ (_B) 'k'] = 't';
  bscX_vec[ (_B) 'M'] = 'a'; bscX_vec[ (_B) 'm'] = 'a';
  bscX_vec[ (_B) 'B'] = 't'; bscX_vec[ (_B) 'b'] = 't';
  bscX_vec[ (_B) 'D'] = 'a'; bscX_vec[ (_B) 'd'] = 'a';
  bscX_vec[ (_B) 'H'] = 'a'; bscX_vec[ (_B) 'h'] = 'a';
  bscX_vec[ (_B) 'V'] = 'a'; bscX_vec[ (_B) 'v'] = 'a';
  bscX_vec[ (_B) 'N'] = 'N'; bscX_vec[ (_B) 'n'] = 'n';
  
  init_bscX = true;
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                       Unmethylated Conversion Table::
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static void bscU_init() {
  
  /*
   * 
   * Should loop through entire ASCII space to map in one-to-one for most
   *   ASCII characters where this is meaningless, or there is no covnersion.
   *   This prevents null pointer exections when non-converted nucleotides
   *   are passed in. After we then assign the true conversions.
   *   
   *           if (uc) x <- tr(x, 'CYSMBHV', 'TTKWKWD')
   *           else    x <- tr(x, 'CYSMBHV', 'ttkwkwd')
   *   
   */
  
  for (_I ii = 0; ii < 256; ii++)
    bscU_vec[ii] = (_B) ii;
  
  bscU_vec[ (_B) 'C'] = 't'; bscU_vec[ (_B) 'c'] = 't';
  bscU_vec[ (_B) 'Y'] = 't'; bscU_vec[ (_B) 'y'] = 't';
  bscU_vec[ (_B) 'S'] = 'k'; bscU_vec[ (_B) 's'] = 'k';
  bscU_vec[ (_B) 'M'] = 'w'; bscU_vec[ (_B) 'm'] = 'w';
  bscU_vec[ (_B) 'B'] = 'k'; bscU_vec[ (_B) 'b'] = 'k';
  bscU_vec[ (_B) 'H'] = 'w'; bscU_vec[ (_B) 'h'] = 'w';
  bscU_vec[ (_B) 'V'] = 'd'; bscU_vec[ (_B) 'v'] = 'd';
  bscU_vec[ (_B) 'N'] = 'N'; bscU_vec[ (_B) 'n'] = 'n';
  
  init_bscU = true;
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                        Methylated Conversion Table::
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static void bscM_init() {
  
  /*
   * Second Pass Conversion Array::
   *   This performs basic bisulfite conversion on nucleotides that were not
   *   protected.
   * 
   *   Should loop through entire ASCII space to map in one-to-one for most
   *   ASCII characters where this is meaningless, or there is no conversion.
   *   This prevents null pointer exections when non-converted nucleotides
   *   are passed in. After we then assign the true conversions.
   *   
   *                  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
   *
   */
  
  for (_I ii = 0; ii < 256; ii++)
    bscM2_vec[ii] = (_B) ii;
  
  bscM2_vec[ (_B) 'C'] = 't';
  bscM2_vec[ (_B) 'Y'] = 't';
  bscM2_vec[ (_B) 'S'] = 'k';
  bscM2_vec[ (_B) 'M'] = 'w';
  bscM2_vec[ (_B) 'B'] = 'k';
  bscM2_vec[ (_B) 'H'] = 'w';
  bscM2_vec[ (_B) 'V'] = 'd';
  bscM2_vec[ (_B) 'N'] = 'N';
  
  /*
   * First Pass Conversion Array:: bscU1_vec()
   * 
   *   Should loop through entire ASCII space and fill in all combinations
   *   then update the actual conversion pairs after...
   *   The default is now in three loops 256 x 256 x 2::
   *   Should loop through entire ASCII space to map in one-to-one for most
   *   ASCII characters where this is meaningless, or there is no conversion.
   *   This prevents null pointer exections when non-converted nucleotides
   *   are passed in. After we then assign the true conversions.
   *   
   */
  
  for (_I ii = 0; ii < 256; ii++)
    for (_I jj = 0; jj < 256; jj++)
      for (_I kk = 0; kk < 2; kk++) {
        if (kk == 0) bscM1_vec[ii][jj][kk] = (_B) ii;
        if (kk == 1) bscM1_vec[ii][jj][kk] = (_B) jj;
      }
      
      bscM1_vec[ (_B) 'C' ][  (_B) 'G' ][0] = 'c';
  bscM1_vec[ (_B) 'C' ][  (_B) 'G' ][1] = 'G';
  bscM1_vec[ (_B) 'C' ][  (_B) 'R' ][0] = 'c';
  bscM1_vec[ (_B) 'C' ][  (_B) 'R' ][1] = 'R';
  bscM1_vec[ (_B) 'C' ][  (_B) 'S' ][0] = 'c';
  bscM1_vec[ (_B) 'C' ][  (_B) 'S' ][1] = 'S';
  bscM1_vec[ (_B) 'C' ][  (_B) 'K' ][0] = 'c';
  bscM1_vec[ (_B) 'C' ][  (_B) 'K' ][1] = 'K';
  bscM1_vec[ (_B) 'C' ][  (_B) 'B' ][0] = 'c';
  bscM1_vec[ (_B) 'C' ][  (_B) 'B' ][1] = 'B';
  bscM1_vec[ (_B) 'C' ][  (_B) 'D' ][0] = 'c';
  bscM1_vec[ (_B) 'C' ][  (_B) 'D' ][1] = 'D';
  bscM1_vec[ (_B) 'C' ][  (_B) 'V' ][0] = 'c';
  bscM1_vec[ (_B) 'C' ][  (_B) 'V' ][1] = 'V';
  
  bscM1_vec[ (_B) 'Y' ][  (_B) 'G' ][0] = 'y';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'G' ][1] = 'G';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'R' ][0] = 'y';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'R' ][1] = 'R';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'S' ][0] = 'y';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'S' ][1] = 'S';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'K' ][0] = 'y';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'K' ][1] = 'K';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'B' ][0] = 'y';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'B' ][1] = 'B';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'D' ][0] = 'y';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'D' ][1] = 'D';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'V' ][0] = 'y';
  bscM1_vec[ (_B) 'Y' ][  (_B) 'V' ][1] = 'V';
  
  bscM1_vec[ (_B) 'S' ][  (_B) 'G' ][0] = 's';
  bscM1_vec[ (_B) 'S' ][  (_B) 'G' ][1] = 'G';
  bscM1_vec[ (_B) 'S' ][  (_B) 'R' ][0] = 's';
  bscM1_vec[ (_B) 'S' ][  (_B) 'R' ][1] = 'R';
  bscM1_vec[ (_B) 'S' ][  (_B) 'S' ][0] = 's';
  bscM1_vec[ (_B) 'S' ][  (_B) 'S' ][1] = 'S';
  bscM1_vec[ (_B) 'S' ][  (_B) 'K' ][0] = 's';
  bscM1_vec[ (_B) 'S' ][  (_B) 'K' ][1] = 'K';
  bscM1_vec[ (_B) 'S' ][  (_B) 'B' ][0] = 's';
  bscM1_vec[ (_B) 'S' ][  (_B) 'B' ][1] = 'B';
  bscM1_vec[ (_B) 'S' ][  (_B) 'D' ][0] = 's';
  bscM1_vec[ (_B) 'S' ][  (_B) 'D' ][1] = 'D';
  bscM1_vec[ (_B) 'S' ][  (_B) 'V' ][0] = 's';
  bscM1_vec[ (_B) 'S' ][  (_B) 'V' ][1] = 'V';
  
  bscM1_vec[ (_B) 'M' ][  (_B) 'G' ][0] = 'm';
  bscM1_vec[ (_B) 'M' ][  (_B) 'G' ][1] = 'G';
  bscM1_vec[ (_B) 'M' ][  (_B) 'R' ][0] = 'm';
  bscM1_vec[ (_B) 'M' ][  (_B) 'R' ][1] = 'R';
  bscM1_vec[ (_B) 'M' ][  (_B) 'S' ][0] = 'm';
  bscM1_vec[ (_B) 'M' ][  (_B) 'S' ][1] = 'S';
  bscM1_vec[ (_B) 'M' ][  (_B) 'K' ][0] = 'm';
  bscM1_vec[ (_B) 'M' ][  (_B) 'K' ][1] = 'K';
  bscM1_vec[ (_B) 'M' ][  (_B) 'B' ][0] = 'm';
  bscM1_vec[ (_B) 'M' ][  (_B) 'B' ][1] = 'B';
  bscM1_vec[ (_B) 'M' ][  (_B) 'D' ][0] = 'm';
  bscM1_vec[ (_B) 'M' ][  (_B) 'D' ][1] = 'D';
  bscM1_vec[ (_B) 'M' ][  (_B) 'V' ][0] = 'm';
  bscM1_vec[ (_B) 'M' ][  (_B) 'V' ][1] = 'V';
  
  bscM1_vec[ (_B) 'B' ][  (_B) 'G' ][0] = 'b';
  bscM1_vec[ (_B) 'B' ][  (_B) 'G' ][1] = 'G';
  bscM1_vec[ (_B) 'B' ][  (_B) 'R' ][0] = 'b';
  bscM1_vec[ (_B) 'B' ][  (_B) 'R' ][1] = 'R';
  bscM1_vec[ (_B) 'B' ][  (_B) 'S' ][0] = 'b';
  bscM1_vec[ (_B) 'B' ][  (_B) 'S' ][1] = 'S';
  bscM1_vec[ (_B) 'B' ][  (_B) 'K' ][0] = 'b';
  bscM1_vec[ (_B) 'B' ][  (_B) 'K' ][1] = 'K';
  bscM1_vec[ (_B) 'B' ][  (_B) 'B' ][0] = 'b';
  bscM1_vec[ (_B) 'B' ][  (_B) 'B' ][1] = 'B';
  bscM1_vec[ (_B) 'B' ][  (_B) 'D' ][0] = 'b';
  bscM1_vec[ (_B) 'B' ][  (_B) 'D' ][1] = 'D';
  bscM1_vec[ (_B) 'B' ][  (_B) 'V' ][0] = 'b';
  bscM1_vec[ (_B) 'B' ][  (_B) 'V' ][1] = 'V';
  
  bscM1_vec[ (_B) 'H' ][  (_B) 'G' ][0] = 'h';
  bscM1_vec[ (_B) 'H' ][  (_B) 'G' ][1] = 'G';
  bscM1_vec[ (_B) 'H' ][  (_B) 'R' ][0] = 'h';
  bscM1_vec[ (_B) 'H' ][  (_B) 'R' ][1] = 'R';
  bscM1_vec[ (_B) 'H' ][  (_B) 'S' ][0] = 'h';
  bscM1_vec[ (_B) 'H' ][  (_B) 'S' ][1] = 'S';
  bscM1_vec[ (_B) 'H' ][  (_B) 'K' ][0] = 'h';
  bscM1_vec[ (_B) 'H' ][  (_B) 'K' ][1] = 'K';
  bscM1_vec[ (_B) 'H' ][  (_B) 'B' ][0] = 'h';
  bscM1_vec[ (_B) 'H' ][  (_B) 'B' ][1] = 'B';
  bscM1_vec[ (_B) 'H' ][  (_B) 'D' ][0] = 'h';
  bscM1_vec[ (_B) 'H' ][  (_B) 'D' ][1] = 'D';
  bscM1_vec[ (_B) 'H' ][  (_B) 'V' ][0] = 'h';
  bscM1_vec[ (_B) 'H' ][  (_B) 'V' ][1] = 'V';
  
  bscM1_vec[ (_B) 'V' ][  (_B) 'G' ][0] = 'v';
  bscM1_vec[ (_B) 'V' ][  (_B) 'G' ][1] = 'G';
  bscM1_vec[ (_B) 'V' ][  (_B) 'R' ][0] = 'v';
  bscM1_vec[ (_B) 'V' ][  (_B) 'R' ][1] = 'R';
  bscM1_vec[ (_B) 'V' ][  (_B) 'S' ][0] = 'v';
  bscM1_vec[ (_B) 'V' ][  (_B) 'S' ][1] = 'S';
  bscM1_vec[ (_B) 'V' ][  (_B) 'K' ][0] = 'v';
  bscM1_vec[ (_B) 'V' ][  (_B) 'K' ][1] = 'K';
  bscM1_vec[ (_B) 'V' ][  (_B) 'B' ][0] = 'v';
  bscM1_vec[ (_B) 'V' ][  (_B) 'B' ][1] = 'B';
  bscM1_vec[ (_B) 'V' ][  (_B) 'D' ][0] = 'v';
  bscM1_vec[ (_B) 'V' ][  (_B) 'D' ][1] = 'D';
  bscM1_vec[ (_B) 'V' ][  (_B) 'V' ][0] = 'v';
  bscM1_vec[ (_B) 'V' ][  (_B) 'V' ][1] = 'V';
  
  init_bscM = true;
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                      Degenerate Conversion Table::
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static void bscD_init() {
  
  /*
   * Second Pass Conversion Array::
   *   This performs basic bisulfite conversion on nucleotides that were not
   *   protected.
   * 
   *   Should loop through entire ASCII space to map in one-to-one for most
   *   ASCII characters where this is meaningless, or there is no covnersion.
   *   This prevents null pointer exections when non-converted nucleotides
   *   are passed in. After we then assign the true conversions.
   *   
   *        stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPD)
   *
   */
  
  //  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  for (_I ii = 0; ii < 256; ii++)
    bscD2_vec[ii] = (_B) ii;
  
  bscD2_vec[ (_B) 'C'] = 't';
  bscD2_vec[ (_B) 'Y'] = 't';
  bscD2_vec[ (_B) 'S'] = 'k';
  bscD2_vec[ (_B) 'M'] = 'w';
  bscD2_vec[ (_B) 'B'] = 'k';
  bscD2_vec[ (_B) 'H'] = 'w';
  bscD2_vec[ (_B) 'V'] = 'd';
  bscD2_vec[ (_B) 'N'] = 'N';
  
  /*
   * First Pass Conversion Array:: bscU1_vec()
   * 
   *   Should loop through entire ASCII space and fill in all combinations
   *   then update the actual conversion pairs after...
   *   The default is now in three loops 256 x 256 x 2::
   *   Should loop through entire ASCII space to map in one-to-one for most
   *   ASCII characters where this is meaningless, or there is no conversion.
   *   This prevents null pointer exections when non-converted nucleotides
   *   are passed in. After we then assign the true conversions.
   *   
   */
  
  for (_I ii = 0; ii < 256; ii++)
    for (_I jj = 0; jj < 256; jj++)
      for (_I kk = 0; kk < 2; kk++) {
        if (kk == 0) bscD1_vec[ii][jj][kk] = (_B) ii;
        if (kk == 1) bscD1_vec[ii][jj][kk] = (_B) jj;
      }
      
      bscD1_vec[ (_B) 'C' ][  (_B) 'G' ][0] = 'y';
  bscD1_vec[ (_B) 'C' ][  (_B) 'G' ][1] = 'G';
  bscD1_vec[ (_B) 'C' ][  (_B) 'R' ][0] = 'y';
  bscD1_vec[ (_B) 'C' ][  (_B) 'R' ][1] = 'R';
  bscD1_vec[ (_B) 'C' ][  (_B) 'S' ][0] = 'y';
  bscD1_vec[ (_B) 'C' ][  (_B) 'S' ][1] = 'S';
  bscD1_vec[ (_B) 'C' ][  (_B) 'K' ][0] = 'y';
  bscD1_vec[ (_B) 'C' ][  (_B) 'K' ][1] = 'K';
  bscD1_vec[ (_B) 'C' ][  (_B) 'B' ][0] = 'y';
  bscD1_vec[ (_B) 'C' ][  (_B) 'B' ][1] = 'B';
  bscD1_vec[ (_B) 'C' ][  (_B) 'D' ][0] = 'y';
  bscD1_vec[ (_B) 'C' ][  (_B) 'D' ][1] = 'D';
  bscD1_vec[ (_B) 'C' ][  (_B) 'V' ][0] = 'y';
  bscD1_vec[ (_B) 'C' ][  (_B) 'V' ][1] = 'V';
  
  bscD1_vec[ (_B) 'Y' ][  (_B) 'G' ][0] = 'y';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'G' ][1] = 'G';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'R' ][0] = 'y';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'R' ][1] = 'R';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'S' ][0] = 'y';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'S' ][1] = 'S';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'K' ][0] = 'y';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'K' ][1] = 'K';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'B' ][0] = 'y';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'B' ][1] = 'B';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'D' ][0] = 'y';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'D' ][1] = 'D';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'V' ][0] = 'y';
  bscD1_vec[ (_B) 'Y' ][  (_B) 'V' ][1] = 'V';
  
  bscD1_vec[ (_B) 'S' ][  (_B) 'G' ][0] = 'b';
  bscD1_vec[ (_B) 'S' ][  (_B) 'G' ][1] = 'G';
  bscD1_vec[ (_B) 'S' ][  (_B) 'R' ][0] = 'b';
  bscD1_vec[ (_B) 'S' ][  (_B) 'R' ][1] = 'R';
  bscD1_vec[ (_B) 'S' ][  (_B) 'S' ][0] = 'b';
  bscD1_vec[ (_B) 'S' ][  (_B) 'S' ][1] = 'S';
  bscD1_vec[ (_B) 'S' ][  (_B) 'K' ][0] = 'b';
  bscD1_vec[ (_B) 'S' ][  (_B) 'K' ][1] = 'K';
  bscD1_vec[ (_B) 'S' ][  (_B) 'B' ][0] = 'b';
  bscD1_vec[ (_B) 'S' ][  (_B) 'B' ][1] = 'B';
  bscD1_vec[ (_B) 'S' ][  (_B) 'D' ][0] = 'b';
  bscD1_vec[ (_B) 'S' ][  (_B) 'D' ][1] = 'D';
  bscD1_vec[ (_B) 'S' ][  (_B) 'V' ][0] = 'b';
  bscD1_vec[ (_B) 'S' ][  (_B) 'V' ][1] = 'V';
  
  bscD1_vec[ (_B) 'M' ][  (_B) 'G' ][0] = 'h';
  bscD1_vec[ (_B) 'M' ][  (_B) 'G' ][1] = 'G';
  bscD1_vec[ (_B) 'M' ][  (_B) 'R' ][0] = 'h';
  bscD1_vec[ (_B) 'M' ][  (_B) 'R' ][1] = 'R';
  bscD1_vec[ (_B) 'M' ][  (_B) 'S' ][0] = 'h';
  bscD1_vec[ (_B) 'M' ][  (_B) 'S' ][1] = 'S';
  bscD1_vec[ (_B) 'M' ][  (_B) 'K' ][0] = 'h';
  bscD1_vec[ (_B) 'M' ][  (_B) 'K' ][1] = 'K';
  bscD1_vec[ (_B) 'M' ][  (_B) 'B' ][0] = 'h';
  bscD1_vec[ (_B) 'M' ][  (_B) 'B' ][1] = 'B';
  bscD1_vec[ (_B) 'M' ][  (_B) 'D' ][0] = 'h';
  bscD1_vec[ (_B) 'M' ][  (_B) 'D' ][1] = 'D';
  bscD1_vec[ (_B) 'M' ][  (_B) 'V' ][0] = 'h';
  bscD1_vec[ (_B) 'M' ][  (_B) 'V' ][1] = 'V';
  
  bscD1_vec[ (_B) 'B' ][  (_B) 'G' ][0] = 'b';
  bscD1_vec[ (_B) 'B' ][  (_B) 'G' ][1] = 'G';
  bscD1_vec[ (_B) 'B' ][  (_B) 'R' ][0] = 'b';
  bscD1_vec[ (_B) 'B' ][  (_B) 'R' ][1] = 'R';
  bscD1_vec[ (_B) 'B' ][  (_B) 'S' ][0] = 'b';
  bscD1_vec[ (_B) 'B' ][  (_B) 'S' ][1] = 'S';
  bscD1_vec[ (_B) 'B' ][  (_B) 'K' ][0] = 'b';
  bscD1_vec[ (_B) 'B' ][  (_B) 'K' ][1] = 'K';
  bscD1_vec[ (_B) 'B' ][  (_B) 'B' ][0] = 'b';
  bscD1_vec[ (_B) 'B' ][  (_B) 'B' ][1] = 'B';
  bscD1_vec[ (_B) 'B' ][  (_B) 'D' ][0] = 'b';
  bscD1_vec[ (_B) 'B' ][  (_B) 'D' ][1] = 'D';
  bscD1_vec[ (_B) 'B' ][  (_B) 'V' ][0] = 'b';
  bscD1_vec[ (_B) 'B' ][  (_B) 'V' ][1] = 'V';
  
  bscD1_vec[ (_B) 'H' ][  (_B) 'G' ][0] = 'h';
  bscD1_vec[ (_B) 'H' ][  (_B) 'G' ][1] = 'G';
  bscD1_vec[ (_B) 'H' ][  (_B) 'R' ][0] = 'h';
  bscD1_vec[ (_B) 'H' ][  (_B) 'R' ][1] = 'R';
  bscD1_vec[ (_B) 'H' ][  (_B) 'S' ][0] = 'h';
  bscD1_vec[ (_B) 'H' ][  (_B) 'S' ][1] = 'S';
  bscD1_vec[ (_B) 'H' ][  (_B) 'K' ][0] = 'h';
  bscD1_vec[ (_B) 'H' ][  (_B) 'K' ][1] = 'K';
  bscD1_vec[ (_B) 'H' ][  (_B) 'B' ][0] = 'h';
  bscD1_vec[ (_B) 'H' ][  (_B) 'B' ][1] = 'B';
  bscD1_vec[ (_B) 'H' ][  (_B) 'D' ][0] = 'h';
  bscD1_vec[ (_B) 'H' ][  (_B) 'D' ][1] = 'D';
  bscD1_vec[ (_B) 'H' ][  (_B) 'V' ][0] = 'h';
  bscD1_vec[ (_B) 'H' ][  (_B) 'V' ][1] = 'V';
  
  bscD1_vec[ (_B) 'V' ][  (_B) 'G' ][0] = 'n';
  bscD1_vec[ (_B) 'V' ][  (_B) 'G' ][1] = 'G';
  bscD1_vec[ (_B) 'V' ][  (_B) 'R' ][0] = 'n';
  bscD1_vec[ (_B) 'V' ][  (_B) 'R' ][1] = 'R';
  bscD1_vec[ (_B) 'V' ][  (_B) 'S' ][0] = 'n';
  bscD1_vec[ (_B) 'V' ][  (_B) 'S' ][1] = 'S';
  bscD1_vec[ (_B) 'V' ][  (_B) 'K' ][0] = 'n';
  bscD1_vec[ (_B) 'V' ][  (_B) 'K' ][1] = 'K';
  bscD1_vec[ (_B) 'V' ][  (_B) 'B' ][0] = 'n';
  bscD1_vec[ (_B) 'V' ][  (_B) 'B' ][1] = 'B';
  bscD1_vec[ (_B) 'V' ][  (_B) 'D' ][0] = 'n';
  bscD1_vec[ (_B) 'V' ][  (_B) 'D' ][1] = 'D';
  bscD1_vec[ (_B) 'V' ][  (_B) 'V' ][0] = 'n';
  bscD1_vec[ (_B) 'V' ][  (_B) 'V' ][1] = 'V';
  
  init_bscD = true;
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                      Degenerate Conversion Table::
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static void topbot_init() {
  
  for (_I ii = 0; ii < 256; ii++)
    for (_I jj = 0; jj < 256; jj++)
      topbot_vec[ii][jj] = 0;
  
  /*
   * 
   *  A/C  => TOP => 1
   *  A/G  => TOP => 1
   *  T/C  => TOP => 1
   *  T/G  => TOP => 1
   *
   *  C/A  => BOT => -1
   *  C/T  => BOT => -1
   *  G/A  => BOT => -1
   *  G/T  => BOT => -1
   *
   *  T == 84
   *  B == 66
   *  
   */
  
  /*
   * 
   * TOP Values::
   * 
   */
  topbot_vec[ (_B) 'A' ][  (_B) 'C' ] =  1;
  topbot_vec[ (_B) 'A' ][  (_B) 'c' ] =  1;
  topbot_vec[ (_B) 'A' ][  (_B) 'G' ] =  1;
  topbot_vec[ (_B) 'A' ][  (_B) 'g' ] =  1;
  topbot_vec[ (_B) 'a' ][  (_B) 'C' ] =  1;
  topbot_vec[ (_B) 'a' ][  (_B) 'c' ] =  1;
  topbot_vec[ (_B) 'a' ][  (_B) 'G' ] =  1;
  topbot_vec[ (_B) 'a' ][  (_B) 'g' ] =  1;
  
  topbot_vec[ (_B) 'T' ][  (_B) 'C' ] =  1;
  topbot_vec[ (_B) 'T' ][  (_B) 'c' ] =  1;
  topbot_vec[ (_B) 'T' ][  (_B) 'G' ] =  1;
  topbot_vec[ (_B) 'T' ][  (_B) 'g' ] =  1;
  topbot_vec[ (_B) 't' ][  (_B) 'C' ] =  1;
  topbot_vec[ (_B) 't' ][  (_B) 'c' ] =  1;
  topbot_vec[ (_B) 't' ][  (_B) 'G' ] =  1;
  topbot_vec[ (_B) 't' ][  (_B) 'g' ] =  1;
  
  /*
   * 
   * BOT Values::
   * 
   */
  topbot_vec[ (_B) 'C' ][  (_B) 'A' ] = -1;
  topbot_vec[ (_B) 'C' ][  (_B) 'a' ] = -1;
  topbot_vec[ (_B) 'C' ][  (_B) 'T' ] = -1;
  topbot_vec[ (_B) 'C' ][  (_B) 't' ] = -1;
  topbot_vec[ (_B) 'c' ][  (_B) 'A' ] = -1;
  topbot_vec[ (_B) 'c' ][  (_B) 'a' ] = -1;
  topbot_vec[ (_B) 'c' ][  (_B) 'T' ] = -1;
  topbot_vec[ (_B) 'c' ][  (_B) 't' ] = -1;
  
  topbot_vec[ (_B) 'G' ][  (_B) 'A' ] = -1;
  topbot_vec[ (_B) 'G' ][  (_B) 'a' ] = -1;
  topbot_vec[ (_B) 'G' ][  (_B) 'T' ] = -1;
  topbot_vec[ (_B) 'G' ][  (_B) 't' ] = -1;
  topbot_vec[ (_B) 'g' ][  (_B) 'A' ] = -1;
  topbot_vec[ (_B) 'g' ][  (_B) 'a' ] = -1;
  topbot_vec[ (_B) 'g' ][  (_B) 'T' ] = -1;
  topbot_vec[ (_B) 'g' ][  (_B) 't' ] = -1;
  
  //
  // BOT Values:: Extended R/Y
  //
  // topbot_vec[ (_B) 'C' ][  (_B) 'Y' ] = -1;
  // topbot_vec[ (_B) 'C' ][  (_B) 'y' ] = -1;
  // topbot_vec[ (_B) 'C' ][  (_B) 'R' ] = -1;
  // topbot_vec[ (_B) 'C' ][  (_B) 'r' ] = -1;
  // topbot_vec[ (_B) 'c' ][  (_B) 'Y' ] = -1;
  // topbot_vec[ (_B) 'c' ][  (_B) 'y' ] = -1;
  // topbot_vec[ (_B) 'c' ][  (_B) 'R' ] = -1;
  // topbot_vec[ (_B) 'c' ][  (_B) 'r' ] = -1;
  // 
  // topbot_vec[ (_B) 'G' ][  (_B) 'Y' ] = -1;
  // topbot_vec[ (_B) 'G' ][  (_B) 'y' ] = -1;
  // topbot_vec[ (_B) 'G' ][  (_B) 'R' ] = -1;
  // topbot_vec[ (_B) 'G' ][  (_B) 'r' ] = -1;
  // topbot_vec[ (_B) 'g' ][  (_B) 'Y' ] = -1;
  // topbot_vec[ (_B) 'g' ][  (_B) 'y' ] = -1;
  // topbot_vec[ (_B) 'g' ][  (_B) 'R' ] = -1;
  // topbot_vec[ (_B) 'g' ][  (_B) 'r' ] = -1;
  
  init_topbot = true;
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                            BSMAP Offset Look Up::
 *                            
 * TBD:: This is a lazy function lookup...
 *   It should be replaced with a table like
 *   the ones above. Simple for now...
 * NOTE:: Could just return a pre-loaded std::map
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

inline _I get_bsmap_offset( _S& pt,
                            _I inf,
                            _S& srd )
{
  _i offset = -999999;
  if ( pt.compare("cg") == 0 ) {
    if ( inf == 1 ) {
      if ( srd.compare("--") == 0 ) {
        offset = 48;
      } else if ( srd.compare("-+") == 0 ) {
        offset = -1;
      } else if ( srd.compare("+-") == 0 ) {
        offset = 0;
      } else if ( srd.compare("++") == 0 ) {
        offset = 49;
      } else {
        // Throw error...
      }
      
    } else if ( inf == 2 ) {
      if ( srd.compare("--") == 0 ) {
        offset = 49;
      } else if ( srd.compare("-+") == 0 ) {
        offset = -2;
      } else if ( srd.compare("+-") == 0 ) {
        offset = -1;
      } else if ( srd.compare("++") == 0 ) {
        offset = 50;
      } else {
        // Throw error...
      }
      
    } else {
      // Throw error...
    }
  }
  
  return( offset );
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                           Basic Math Functions::
 *                              
 *  Internal Functions::
 *   avg_vec()
 *   med_vec()
 *   std_vec()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */


/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *                          Basic Math Functions::
 * 
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

/*
 * Common Mathmatical Operator's Overload:: Vectors
 *  TBD:: Add additional operators/combinations as needed...
 *  
 */

template<typename T>
inline vector<T> operator-( const std::vector<T>& vec, const T c )
{
  vector<T> rvec;
  for ( auto& x : vec ) rvec.push_back( x - c );
  return( rvec );
};

template <typename T>
inline std::vector<T> operator-( const std::vector<T>& x_vec, const std::vector<T>& y_vec )
{
  vector<T> r_vec;
  _t n = x_vec.size();
  if ( y_vec.size() < n ) n = y_vec.size();
  for ( _t ii=0; ii<n; ii++ ) r_vec.push_back( x_vec[ii] - y_vec[ii] );
  return( r_vec );
};

template<typename T>
inline vector<T> operator*( const std::vector<T>& vec, const T c )
{
  vector<T> rvec;
  for ( auto& x : vec ) rvec.push_back( x * c );
  return( rvec );
};

template <typename T>
inline vector<T> operator*( const std::vector<T>& x_vec, const std::vector<T>& y_vec )
{
  vector<T> r_vec;
  _t n = x_vec.size();
  if ( y_vec.size() < n ) n = y_vec.size();
  for ( _t ii=0; ii<n; ii++ ) r_vec.push_back( x_vec[ii] * y_vec[ii] );
  return( r_vec );
};

/*
 * Common Mathmatical Statistics:: Vectors
 *  TBD:: Vectorize (both in place and new vector)
 *  
 */

// Not reals sure what simple_abs() is for...
template<typename T>
inline T simple_abs( T a, T b )
{
  if ( a > b ) return( a - b);
  return( b - a );
};

template<typename T>
inline _d set_precision( T v, _i precision = -1 )
{
  if ( precision >= 0 ) {
    _d pv = pow( 10, precision );
    v = (_d) round( v*pv ) / pv;
  }
  return( v );
};

template<typename T>
inline T sum_vec( const std::vector<T>& v )
{
  T s = 0;
  for ( auto& x : v ) s += x;
  return( s );
};

template<typename T>
inline T median_vec( std::vector<T> vec, _Y is_sorted = false )
{
  _t size = vec.size();
  // The first one returns nan which is what we want, the other two are just
  //  listed for fullness
  if ( size == 0 ) return( sum_vec(vec)/vec.size() );
  // if (size == 0) return(-1.0);
  // if (size == 0) throw domain_error("[median]: median of an empty vector");
  
  if ( !is_sorted ) std::sort( vec.begin(), vec.end() );
  _t mid = int( size/2 );
  
  return( size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2 : vec[mid] );
};

template<typename T>
inline _d mean_vec( const std::vector<T>& vec )
{
  return( sum_vec(vec) / vec.size() );
};

template<typename T>
inline _d sqsum_vec( const std::vector<T>& vec )
{
  _d s = 0;
  for ( auto& x : vec ) s += std::pow(x, 2);
  return( s );
};

template<typename T>
inline _d stdev_vec( const std::vector<T>& vec )
{
  _t size = vec.size();
  return( pow(sqsum_vec(vec) / size - pow(sum_vec(vec) / size, 2), 0.5) );
};

template<typename T>
inline vector<_d> quantile3_vec( std::vector<T> vec, bool is_sorted = false ) {
  _t size = vec.size();
  
  // std::vector<_d> nan_vec( {std::nan,std::nan,std::nan} );
  // if (size == 0) return( nan_vec );
  if ( size == 0 ) throw domain_error("[quantile3_vec]: median of an empty vector");
  if ( !is_sorted ) std::sort( vec.begin(), vec.end() );
  
  _t mid = int( size/2 );
  
  vector<_d> Q1v( vec.cbegin(), vec.cbegin() + mid );
  vector<_d> Q3v( vec.cbegin() + mid+1, vec.cend() );
  
  vector<_d> quantiles;
  quantiles.push_back( median_vec(Q1v, true) );
  quantiles.push_back( median_vec(vec, true) );
  quantiles.push_back( median_vec(Q3v, true) );
  return( quantiles );
};

template<typename T>
inline std::vector<_d> stats5_vec(std::vector<T> vec, _Y is_sorted = false ) {
  _t size = vec.size();
  
  // Should we return a vector of zero and nan?
  // if (size == 0) 
  //   return( std::vector<_d>( {std::nan,std::nan,std::nan,std::nan,std::nan} ) );
  if (size == 0) throw domain_error("[stats5_vec]: median of an empty vector");
  
  vector<_d> stats = quantile3_vec( vec, is_sorted );
  stats.push_back( mean_vec(vec) );
  stats.push_back( stdev_vec(vec) );
  
  return( stats );
};

template<typename T>
inline _d pearsoncoeff_vec( const std::vector<T>& X, const std::vector<T>& Y )
{
  return( sum_vec((X - mean_vec(X))*(Y - mean_vec(Y))) / (X.size()*stdev_vec(X)* stdev_vec(Y)) );
};

// Commonly used for calculating Delta Beta(dB) count::
template <typename T>
inline _I pass_count( const std::vector<T>& X, const std::vector<T>& Y, 
               const _d cut = 0.2, const _i cmp = 0, const _Y is_abs = false )
{
  /*
   *        Int_Value | Comparison | English_Meaning
   *     [ -----------|------------|---------------------------- ]
   * cmp |    -2      |     <      | Less Than
   * cmp |    -1      |     <=     | Less Than or Equal To
   * cmp |     0      |     ==     | Equal To
   * cmp |     1      |     >      | Greater Than or Equal To
   * cmp |     2      |     >=     | Greater Than
   *     [ -----------|------------|---------------------------- ]
   * 
   */
  
  _I p_cnt = 0;
  std::vector<T> D = X - Y;
  if ( is_abs ) for ( _t ii=0; ii<D.size(); ii++ )
    D[ii]=std::abs(D[ii]);
  
  if ( cmp == -2 ) for ( auto& v : D ) if ( v <  cut ) p_cnt++;
  if ( cmp == -1 ) for ( auto& v : D ) if ( v <= cut ) p_cnt++;
  if ( cmp ==  0 ) for ( auto& v : D ) if ( v == cut ) p_cnt++;
  if ( cmp ==  1 ) for ( auto& v : D ) if ( v >= cut ) p_cnt++;
  if ( cmp ==  2 ) for ( auto& v : D ) if ( v >  cut ) p_cnt++;
  return( p_cnt );
};

// Function to find t-test of two set of statistical data.
template <typename T>
inline _d tTest_vec( std::vector<T> X, std::vector<T> Y )
{
  _d sds1 = stdev_vec( X );
  _d sds2 = stdev_vec( Y );
  
  // Formula to find t-test of two set of data.
  _d t_test = (mean_vec(X) - mean_vec(Y)) / std::sqrt( (sds1 * sds1) / X.size() + (sds2 * sds2) / Y.size() );
  
  return( t_test );
}

/*
 * Boost Accumulator Stats::
 *   accumulator_set<double, features<tag::mean, tag::variance>> acc;
 *   
 *   TBD:: Need to find Boost MAD( Median Absolute Deviation )
 */
template <typename T>
inline std::vector<_d> get_vec_stats( const std::vector<T> X,
                               const std::vector<_I>& I,
                               const _I min_size = 0 )
{
  std::vector<_d> r_vec;
  
  if ( (X.size() <= min_size) || 
       (I.size() != 0 && I.size() <= min_size ) ) return( r_vec );
  
  accumulator_set<T, features<tag::count, tag::min, tag::max, 
                              tag::mean,  tag::variance, tag::median > > acc;
  
  if ( I.size() != 0 ) {
    for ( _t ii=0; ii<I.size(); ii++ ) acc( X[ I[ii] ] );
  } else {
    for ( auto& x : X ) acc(x);
  }
  
  r_vec.push_back( (_d) boost::accumulators::count(acc) );
  r_vec.push_back( (_d) boost::accumulators::min(acc) );
  r_vec.push_back( (_d) boost::accumulators::max(acc) );
  r_vec.push_back( (_d) boost::accumulators::mean(acc) );
  r_vec.push_back( (_d) std::sqrt(boost::accumulators::variance(acc)) );
  r_vec.push_back( (_d) boost::accumulators::median(acc) );
  
  return( r_vec );
};

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                              System Functions::
 *                              
 *  Internal Functions::
 *   exec_cmd()
 *   is_gzipped()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

/*
 * exec_cmd:: c++ function to run a command
 * 
 */
static _S exec_cmd( _S cmd,
                    const _I vb = 0,
                    const _I vt = 2 )
{
  const _S func_tag("[exec_cmd]: ");
  const _S func_err("[exec_cmd]: ERROR: ");
  const _S func_wrn("[exec_cmd]: Warning: ");
  
  const _Y p1 = vb >= vt + 1;
  
  if ( p1 ) std::cerr 
    << func_tag << "\tExecuting Cmd = '" << cmd << "'\n" << std::endl;
  
  _c buffer[128];
  _S result = "";
  FILE* pipe = popen(cmd.c_str(), "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  
  try {
    while (fgets(buffer, sizeof buffer, pipe) != NULL) {
      result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  
  return result;
};

/*
 * exec_cmd:: c++ function to check if a file is actually gzipped based on 
 *   contents and NOT file extension.
 * 
 */
static _Y is_gzipped( _S file,
                      
                      const _I vb = 0,
                      const _I vt = 1,
                      const _I tc = 0,
                      const _S ft="is_gzipped" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb > vt + 0;
  
  _c buffer[128];
  _S result = "";
  
  // First determine if the file is gzipped::
  _S cmd = "file " + file;
  
  FILE* pipe = popen(cmd.c_str(), "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  try {
    while (fgets(buffer, sizeof buffer, pipe) != NULL) {
      result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  
  _Y is_zipped = true;
  _i res_idx = result.find("gzip compressed");
  if ( res_idx < 0 ) is_zipped = false;
  
  _S g_mssg = " is NOT gzipped!";
  if ( is_zipped ) g_mssg = " is gzipped!";
  
  if ( p0 ) std::cerr
    <<_fo<< "File = '"<<file<<"':: '"<<g_mssg<<"'\n"
    <<std::endl;
  
  return( is_zipped );
};

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * Internal Static c++ Functions::
 * 
 *                         Basic File IO:: Validation
 *                         
 * Internal::
 *     file_exist()
 *     file_readable()
 *     legal_file()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// static _Y path_exists(const _S &s)
// {
//   struct stat buffer;
//   return (stat (s.c_str(), &buffer) == 0);
// }

/*
 * Basic System IO:: build_path(recurssive(build_directory))
 *  NOTE:: Easiest solution:: SYSTEM CALL!!!
 *  
 *  get_file_name() == fs::path(filename).filename()
 *  
 *  Unfortunately #include <filesystem> doesn't work on macOS
 *  
 static _S get_file_name(const string& s) {
 
 size_t i = s.rfind(SYS_SEP, s.length());
 if (i != string::npos) {
 return( s.substr(i+1, s.length() - i) );
 }
 return( "" );
 }
 */

static _Y path_exists( const _S& s ) {
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}

static _Y build_path( const _S& s ) {
  
  _Y present = path_exists( s );
  // std::cerr << "STATUS: present='"<<present<<"'"<<std::endl;
  if ( present ) return( true );
  
  _S cmd;
  cmd = "mkdir -p "+s;
  // std::cerr << "STATUS: cmd='"<<cmd<<"'"<<std::endl;
  
  _i sys_ret = system( cmd.c_str() );
  // _Y pos_val = path_exists( s );
  // std::cerr << "STATUS: sys_ret='"<<sys_ret<<"', pos_val='"<<pos_val<<"'"<<std::endl;                                                                                                                                              
  if ( sys_ret == 0 ) return( true );
  
  std::cerr << "ERROR: Failed Command: '"<<cmd<<"'\n"
            << "ERROR: System Return = '"<<sys_ret<<"'\n"
            << std::endl;
  // assert( sys_ret == 0 );                                                                                                                                                                               
  return( false );
}

static _Y file_exist( _S path ) {
  FILE *fp = fopen( path.c_str(), "r");
  if ( fp ) {
    fclose(fp);
    return true;
  }
  return errno != ENOENT;
};

static _Y file_readable( _S path ) {
  FILE *fp = fopen( path.c_str(), "r");
  if ( fp ) {
    fclose(fp);
    return true;
  }
  return errno != ENOENT && errno != EPERM;
};

static _Y legal_file( _S path,
                      
                      const _I vb = 0,
                      const _I vt = 4,
                      const _I tc = 0,
                      const _S ft="legal_file" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p8 = vb > vt + 8;
  
  _Y success = true;
  success = file_exist( path );
  if ( p8 && !success ) std::cerr
    <<_fo<<"File Does not exist["<<success<<"] = '"<<path<<"'!"<<std::endl;
  if ( !success ) return( success );
  
  success = file_readable( path );
  if ( p8 && !success ) std::cerr
    <<_fo<<"Not readable file["<<success<<"] = '"<<path<<"'!"<<std::endl;
  if ( !success ) return( success );
  
  return( success );
};

static _Y clean_file( _S path )
{
  _Y clean = file_exist( path );
  if ( !clean ) return( !clean );
  clean = std::remove( path.c_str() );
  return( clean );
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * 
 *                         Old Functions for Sorting a Map::
 *                         
 *                              DELTE THESE FUNCTIONS!!!
 *                         
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// Bullshit this doesnt work...
// template<typename K, typename T>
// _Y cmp_temp_desc( std::pair<K, T>& a,
//                     std::pair<K, T>& b ) {
//   return( a.second < b.second );
// };
// 
// Bullshit this doesnt work...
// template<typename K>
// _Y cmp_I_desc( std::pair< K, _I>& a,
//                  std::pair< K, _I>& b ) {
//   return( a.second < b.second );
// };

// _Y cmp_SI_desc( std::pair<_S, _I>& a,
//                 std::pair<_S, _I>& b ) {
//   return( a.second < b.second );
// };

// template<typename K, typename T>
// _Y cmp_temp_ascc( std::pair<K, T>& a,
//                     std::pair<K, T>& b ) {
//   return( a.second < b.second );
// };

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * Internal Static c++ Functions::
 * 
 *                         Binary DNA Conversion Functions::
 *                         
 * Internal::
 *     bit_to_dna()
 *     dna_to_bit()
 *     init_bit_stack()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// Previously static
inline _c* bit_to_dna( uint8_t* m_data,
                       const _t m_len ) {

  _c* dna_str = new _c[m_len + 1];
  
  /* for each base of the DNA sequence */
  for (_t i = 0; i < m_len; ++i)
  {
    uint8_t shift = 6 - 2 * (i % 4);
    uint8_t mask = BASE_MASK << shift;
    
    /* get the i-th DNA base */
    uint8_t base = (m_data[i / 4] & mask) >> shift;
    
    switch (base)
    {
    case BASE_A:
      dna_str[i] = 'A';
      break;
    case BASE_C:
      dna_str[i] = 'C';
      break;
    case BASE_G:
      dna_str[i] = 'G';
      break;
    case BASE_T:
      dna_str[i] = 'T';
      break;
    default:
      throw std::runtime_error("invalid DNA base");
    }
  }
  dna_str[m_len] = '\0';
  
  return( dna_str );
};

inline uint8_t* dna_to_bit( const _c* dna_str,
                            const _t dna_len ) {

  /* number of bytes necessary to store dna_str as a bitset */
  _t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);
  uint8_t* bit_data = new uint8_t[dna_bytes];
  std::memset(bit_data, 0, dna_bytes);

  /* for each base of the DNA sequence */
  for (_t i = 0; i < dna_len; ++i)
  {
    uint8_t shift = 6 - 2 * (i % 4);
    
    switch (dna_str[i])
    {
    case 'A':
      bit_data[i / 4] |= BASE_A << shift;
      break;
    case 'C':
      bit_data[i / 4] |= BASE_C << shift;
      break;
    case 'G':
      bit_data[i / 4] |= BASE_G << shift;
      break;
    case 'T':
      bit_data[i / 4] |= BASE_T << shift;
      break;
    default:
      throw std::invalid_argument("invalid DNA base");
    }
    shift = (shift == 0) ? 6 : shift - 2;
  }
  
  return( bit_data );
};

std::vector<uint8_t*> init_bit_stack( const _I s,
                                      const _I l,
                                      
                                      const _I vb = 0,
                                      const _I vt = 2,
                                      const _I tc = 0,
                                      const _S ft="init_bit_stack" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0  = vb > vt + 0;
  const _Y p1  = vb > vt + 1;
  // const _Y p8  = vb > vt + 8;
  
  std::vector<uint8_t*> bit_stack( s );
  _t dna_bytes = (l / 4) + (l % 4 != 0);
  
  if ( p0 ) std::cerr 
    <<_fo<< "Allocating Binary Stack(Bit Vector)::\n"
    <<_fo<< "\t stack size = " << s << "\n"
    <<_fo<< "\t seq length = " << l << "\n"
    <<_fo<< "\t  dna_bytes = " << dna_bytes << "\n"
    << std::endl;
  
  for ( _I stack_idx = 0; stack_idx < s; stack_idx++ )
  {
    uint8_t* stack_bit = new uint8_t[dna_bytes];
    std::memset(stack_bit, 0, dna_bytes);
    _i bit_size = sizeof(stack_bit);
    
    bit_stack[stack_idx] = stack_bit;
    _c* bit_seq = bit_to_dna( bit_stack[stack_idx], l );
    uint8_t* seq_bit = dna_to_bit( bit_seq, l );
    
    _Y bit_match = false;
    if ( stack_bit == seq_bit ) bit_match = true;
    
    if ( p1 ) std::cerr
      <<_fo<< "New Bit::\n"
      <<_fo<< "\t   bit_size[" << stack_idx << "] = " << bit_size << "\n"
      <<_fo<< "\t bit_to_dna[" << stack_idx << "]:: " << bit_seq << "\n"
      <<_fo<< "\t  bit match[" << stack_idx << "]:: " << bit_match << "\n"
      << std::endl;
  }
  
  return( bit_stack );
};

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                      Auxillary Manifest Functions::
 * 
 * Helper Functions::
 * - is_cgn()
 * - is_ch()
 * - is_rs()
 * - is_nv()
 * - is_ctl()
 * - is_chr()
 * - is_sex_chr()
 * - get_probe_type()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// Helper functions::
inline _Y is_cgn( const _S& cgn ) {
  if ( cgn.size() >= 3 && cgn[0] == 'c' && cgn[1] == 'g' &&
       isdigit( (_i)cgn[2] ) != 0 ) return( true );
  return( false );
};

inline _Y is_ch( const _S& cgn ) {
  if ( cgn.size() >= 3 && cgn[0] == 'c' && cgn[1] == 'h' &&
       ( isdigit( (_i)cgn[2] ) != 0 || cgn[2] == '.' ) ) return( true );
  return( false );
};

inline _Y is_rs( const _S& cgn ) {
  if ( cgn.size() >= 3 && cgn[0] == 'r' && cgn[1] == 's' &&
       isdigit( (_i)cgn[2] ) != 0 ) return( true );
  return( false );
};

inline _Y is_nv( const _S& cgn ) {
  if ( cgn.size() >= 3 && cgn[0] == 'n' && cgn[1] == 'v' && cgn[2] == '_' )
    return( true );
  return( false );
};

inline _Y is_ctl( const _S& cgn ) {
  if ( cgn.size() >= 4 && 
       cgn[0] == 'c' && cgn[1] == 't' && cgn[2] == 'l' && cgn[3] == '_' ) 
    return( true );
  return( false );
};

inline _Y is_chr( const _S& cgn ) {
  if ( cgn.size() >= 4 && 
       cgn[0] == 'c' && cgn[1] == 'h' && cgn[2] == 'r' && 
       ( isdigit( (_i)cgn[3] ) != 0 || 
       cgn[3] == 'X' || cgn[3] == 'Y' || cgn[3] == 'Z' ) ) return( true );
  return( false );
};

inline _Y is_sex_chr( const _S& cgn ) {
  if ( cgn.size() == 4 && is_chr( cgn ) && 
       ( cgn[3] == 'X' || cgn[3] == 'Y' || cgn[3] == 'Z' ) ) return( true );
  return( false );
};

inline _S get_probe_type( const _S& cgn ) {
  if ( is_cgn( cgn ) ) return( "cg" );
  if ( is_ch( cgn ) )  return( "ch" );
  if ( is_rs( cgn ) )  return( "rs" );
  if ( is_nv( cgn ) )  return( "nv" );
  if ( is_ctl( cgn ) ) return( "ct" );
  return( "uk" );
};

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                              Parsing Functions::
 *                              
 *  Internal Functions::
 *   str_to_upper()
 *   str_to_char_idx()
 *   str_to_str_idx()
 *   vec_to_unique_chars()
 *   chop_valid_str()
 *
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

template <typename T>
inline _t identical_vecs( const std::vector<T> a,
                   const std::vector<T> b )
{
  _t  a_size = a.size();
  if ( a_size != b.size() ) return( a_size+1 );
  
  for ( _t ii=0; ii<a_size; ii++ )
    if ( a[ii] != b[ii] ) return( ii+1 );
    
    return( 0 );
}

void str_to_upper( _S& input )
{
  std::for_each(std::begin(input), std::end(input), [](_c& c) {
    c = static_cast<_c>(std::toupper(static_cast<_B>(c)));
  });
}

static _c str_to_char_idx( const _S& s,
                           const _Y uc = false,
                           const _Y lc = false,
                           const _I idx = 0 )
{
  if ( s.size() == 0 ) return(' ');
  _c ch = s[idx];
  if ( uc ) ch = std::toupper( ch );
  if ( lc ) ch = std::toupper( ch );
  return( ch );
};

static _S str_to_str_idx( const _S& s,
                          const _Y uc = false,
                          const _Y lc = false,
                          const _I idx = 0 )
{
  if ( s.size() == 0 ) return("");
  _c ch = s[idx];
  if ( uc ) ch = std::toupper( ch );
  if ( lc ) ch = std::toupper( ch );
  _S s2;
  s2 = ch;
  return( s2 );
};

std::vector<_c> vec_to_unique_chars( std::vector<_S> v )
{
  
  std::sort(v.begin(), v.end());
  auto v_last = std::unique(v.begin(), v.end());
  v.erase(v_last, v.end());
  auto v_len = v.size();
  
  std::vector<_c> s;
  for ( _t ii = 0; ii < v_len; ii++ ) {
    for ( _c ch : v[ii] ) s.push_back( std::toupper( ch ) );
  }
  
  std::sort(s.begin(), s.end());
  auto s_last = std::unique(s.begin(), s.end());
  s.erase(s_last, s.end());
  
  return( s );
}

_Y chop_valid_str( _S s,
                   std::vector<_S>& r,
                   std::vector<_c> v,
                   const _I vb = 0,
                   const _I vt = 3 )
{
  const _S func_tag("[chop_FR_str]: ");
  const _S func_err("[chop_FR_str]: ERROR: ");
  const _S func_wrn("[chop_FR_str]: Warning: ");
  
  _I s_len = s.size();
  _I v_len = v.size();
  
  if ( s_len == 0 )
  {
    std::cerr << func_err << "Strand String is empty!" << "\n"
              << func_err << "\t     s = '" << s << "'\n"
              << func_err << "\t s_len = '" << s_len << "'\n"
              << std::endl;
    return( false );
  }
  if ( v_len == 0 )
  {
    std::cerr << func_err << "Valid Vector is empty!" << "\n"
              << func_err << "\t v_len = '" << v_len << "'\n"
              << std::endl;
    return( false );
  }
  
  // Unique Vector (sort/unique/erase):: v => w
  std::vector<_c> w;
  for ( _c ch : v ) w.push_back( std::toupper( ch ) );
  
  std::sort(w.begin(), w.end());
  auto w_last = std::unique(w.begin(), w.end());
  w.erase(w_last, w.end());
  auto w_len = w.size();
  
  // Unique Vector (sort/unique/erase):: s => c
  std::vector<_c> c;
  for ( _c ch : s ) c.push_back( std::toupper( ch ) );
  
  // Unique Vector (sort/unique/erase):: v
  std::sort(c.begin(), c.end());
  auto c_last = std::unique(c.begin(), c.end());
  c.erase(c_last, c.end());
  auto c_len = c.size();
  
  // Compare All Against All::
  for ( _t ii = 0; ii < c_len; ii++ ) {
    _Y valid = false;
    for ( _t jj = 0; jj < w_len; jj++ ) { 
      if ( c[ii] == w[jj] ) valid = true;
      
      if ( vb >= vt + 1 ) std::cerr 
        << func_tag << "ii = " << ii << ", valid = " << valid << "\n"
        << func_tag << "\t src["<<ii<<"] = " << c[ii] << "\n"
        << func_tag << "\t vec["<<jj<<"] = " << w[jj] << "\n"
        << std::endl;
    }
    if ( !valid ) {
      std::cerr << func_err << "Illegal Strand Value!" << "\n"
                << func_err << "\t     s = '" << s << "'\n"
                << func_err << "\t c[ii] = '" << c[ii] << "'\n"
                << std::endl;
      r.clear();
      return( false );
    }
    
    _S r_str;
    r_str = c[ii];
    r.push_back( r_str );
    if ( vb >= vt + 1) std::cerr 
      << func_tag << "Adding Valid Substring: '" << c[ii] << "'" << std::endl;
  }
  
  return( true );
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 *                              Parsing Functions::
 *                              
 * Internal Functions::
 *   parse_line()
 *   parse_into_postion()
 *   read_file_fast()
 *   parse_str_buffer()
 *   read_file_zlib()
 *   cgntop_file_to_bgz()
 *   load_cgntop_bgz()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

/*
 * parse_line:: split string by delimiter into string vector
 * 
 */
void parse_line( const string& s, 
                 const _c c,
                 std::vector<_S>& v )
{
  _S::size_type i = 0;
  _S::size_type j = s.find(c);
  
  v.clear();
  while (j != _S::npos) {
    v.push_back(s.substr(i, j-i));
    i = ++j;
    j = s.find(c, j);
    
    if ( j == _S::npos) v.push_back(s.substr(i, s.length() ) );
  }
  
  // if ( j == _S::npos ) v.push_back( s );
  if ( v.size() == 0 ) v.push_back( s );
};

/*
 * parse_into_postion:: Wrapper for spliting strings by delimiter into a 
 *   string vector
 * 
 */
static void parse_into_postion( const _S& s, 
                                const _c c,
                                std::vector<std::vector<_S>>& v )
{
  _S::size_type i = 0;
  _S::size_type j = s.find(c);
  
  _i vec_idx = 0;
  while (j != _S::npos) {
    v[vec_idx].push_back(s.substr(i, j-i));
    i = ++j;
    j = s.find(c, j);
    vec_idx++;
    
    if (j == _S::npos)
      v[vec_idx].push_back(s.substr(i, s.length()));
  }
};

static _S read_file_fast( _S path,
                          const _I vb = 0,
                          const _I vt = 3 )
{
  const _S func_tag("[read_file_fast]: ");
  const _S func_err("[read_file_fast]: ERROR: ");
  const _S func_wrn("[read_file_fast]: Warning: ");
  
  const _Y p1 = vb >= vt + 1;
  
  std::ifstream in( path.c_str() );
  _S c_str;
  
  // Scan Complete File from Begining to End::
  in.seekg( 0,std::ios::end );
  
  if ( p1 ) std::cerr 
    << func_tag << "\tFile Size = " << in.tellg() << std::endl;
  
  // Adjust the size of the resulting string with the file size(tellg function)
  c_str.resize( in.tellg() );
  
  // Return to the beginning of the file::
  in.seekg( 0, std::ios::beg );
  
  // Each value the pointer will copied into contents
  in.read( &c_str[0], c_str.size() );
  
  // Close file handle::
  in.close();
  
  if ( p1 ) std::cerr 
    << func_tag << "\tFile Size = " << c_str.size() << "\n" << std::endl;
  
  return( c_str );
};

/*
 * parse_str_buffer:: split string buffer by delimiter into string vector and
 *   return remaining string
 * 
 */
_S parse_str_buffer( const string& s,
                     std::vector<_S>& n,
                     std::vector<std::vector<_S>>& v,
                     
                     const _c sep = ',',
                     const _c ret = '\n',
                     
                     const _I vb = 0,
                     const _I vt = 1,
                     const _I tc = 0,
                     const _S ft="parse_str_buffer" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  // const _Y p0  = vb > vt + 0;
  const _Y p1  = vb > vt + 1;
  // const _Y p8  = vb > vt + 8;
  _S p = "";
  
  // if ( n_len == 0 ) n_len = n.size();
  _I n_len = n.size();
  
  std::vector<_S> lines_vec;
  parse_line( s, ret, lines_vec );
  _I lines_cnt = lines_vec.size();
  
  _I lines_max = lines_cnt;
  if ( lines_cnt > 1 ) lines_max = lines_cnt-1;
  
  // Ensure v matrix is intialized::
  if ( v.size() == 0 ) v.resize( lines_max );
  
  std::vector<_S> data_vec;
  for ( _I ii=0; ii<lines_max; ii++ ) {
    
    parse_line( lines_vec[ii], sep, data_vec );
    _I data_len = data_vec.size();
    
    if ( n_len == 0 ) {
      n_len = data_len;
      if ( p1 ) std::cerr <<_fo<< "Calculated n_len="<<n_len<<"\n" <<std::endl;
      
      v.resize( n_len );
      n.clear();
      n.resize( n_len );
      for ( _I jj=0; jj<n_len; jj++ ) n[jj] = "V" + std::to_string(jj);
    }
    if ( n_len != data_len ) {
      std::cerr
      <<_fe<< "Failed during processing! n_len("<<n_len<<") != data_len("<<data_len<<")\n"
      <<std::endl;
      exit (EXIT_FAILURE);
      // return( "" );
    }
    
    for ( _I jj=0; jj<data_len; jj++ ) {
      v[jj].push_back( data_vec[jj] );
    }
    // if ( p8 ) std::cerr <<_fo<< "line["<<ii<<"]='"<<lines_vec[ii]<<"'" <<std::endl;
  }
  if ( lines_cnt > 1 ) p = lines_vec[ lines_cnt-1 ];
  
  return( p );
};

/*
 * read_file_zlib:: read file with zlib and return string matrix::
 * 
 */
_Y read_file_zlib( const string& file,
                   std::vector<_S>& name_vec,
                   std::vector<std::vector<_S>>& data_mat,
                   
                   const _c sep = ',',
                   const _c ret = '\n',
                   
                   const _I vb = 0,
                   const _I vt = 1,
                   const _I tc = 0,
                   const _S ft="read_file_zlib" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0  = vb > vt + 0;
  const _Y p1  = vb > vt + 1;
  // const _Y p8  = vb > vt + 8;
  
  _Y success = true;
  
  if ( p0 ) std::cerr
    <<_fo<< "Starting..." << "\n"
    <<_fo<< "\t file = '" << file << "'\n" 
    <<_fo<< "\t  sep = '" << sep << "'\n"
    << std::endl;
  
  // _I buff_cnt = 0; // Currently not used for anything...
  _S line_pre = "";
  _I name_len = name_vec.size();
  if ( p1 ) std::cerr <<_fo<< "Original name_len="<<name_len<<"\n" <<std::endl;
  data_mat.resize( name_len );
  
  gzFile gzin = gzopen( file.c_str() ,"r" );
  if ( !gzin ) {
    std::cerr <<_fe<< "Failed gzopen file='"<<file<<"'\n"<<std::endl;
    exit (EXIT_FAILURE);
    // return( false );
  }
  
  while (1) {
    _i err;                    
    _i bytes_read;
    _B buffer[GZ_LENGTH];
    bytes_read = gzread (gzin, buffer, GZ_LENGTH - 1);
    buffer[bytes_read] = '\0';
    
    _S line_buf( reinterpret_cast<_c*>(buffer) );
    line_buf = line_pre + line_buf;
    
    line_pre = parse_str_buffer( line_buf, name_vec, data_mat, sep, ret, 
                                 vb,vt+1,tc+1 );
    
    if (bytes_read < GZ_LENGTH - 1) {
      if ( gzeof (gzin) ) {
        break;
      } else {
        const _c * error_string;
        error_string = gzerror (gzin, & err);
        if ( err ) {
          std::cerr 
          <<_fe<< "Failed during processing! gzerror='"<<error_string<<"'\n"
          <<std::endl;
          exit (EXIT_FAILURE);
          // return( false );
        }
      }
    }
    // buff_cnt++;
  }
  gzclose (gzin);
  
  return( success );
}

/*
 * cgntop_file_to_bgz:: read file with zlib and return string matrix::
 * 
 */
_Y cgntop_file_to_bgz( const string& file,
                       const string& bgz_path,
                       // std::vector<_S>& name_vec,
                       const _I top_idx = 0,
                       
                       const _c sep = ',',
                       const _c ret = '\n',
                       
                       const _I vb = 0,
                       const _I vt = 1,
                       const _I tc = 0,
                       const _S ft="cgntop_file_to_bgz" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0  = vb > vt + 0;
  const _Y p1  = vb > vt + 1;
  // const _Y p8  = vb > vt + 8;
  
  _Y success = true;
  
  if ( p0 ) std::cerr
    <<_fo<< "Starting..." << "\n"
    <<_fo<< "\t     file = '" << file << "'\n" 
    <<_fo<< "\t bgz_path = '" << bgz_path << "'\n"
    <<_fo<< "\t  top_idx = '" << top_idx << "'\n"
    <<_fo<< "\t      sep = '" << sep << "'\n"
    << std::endl;
  
  gzFile gzin = gzopen( file.c_str() ,"r" );
  if ( !gzin ) {
    std::cerr <<_fe<< "Failed gzopen file='"<<file<<"'\n"<<std::endl;
    exit (EXIT_FAILURE);
  }
  
  _I buf_cnt = 0;
  _I seq_len = 0;
  _Y seq_init = false;
  _S line_pre = "";
  ostringstream os;
  std::vector<uint8_t*> bit_vec;
  while (1) {
    _i err;                    
    _i bytes_read;
    _B buffer[GZ_LENGTH];
    bytes_read = gzread (gzin, buffer, GZ_LENGTH - 1);
    buffer[bytes_read] = '\0';
    
    _S line_buf( reinterpret_cast<_c*>(buffer) );
    line_buf = line_pre + line_buf;
    
    std::vector<_S> name_vec;
    std::vector<std::vector<_S>> data_mat;
    line_pre = parse_str_buffer( line_buf, name_vec, data_mat, sep, ret,
                                 vb,vt+1,tc+1 );
    
    for ( _I ii=0; ii<data_mat[top_idx].size(); ii++ ) {
      if ( p1 ) std::cerr <<_fo<< "Top_Seq["<<top_idx<<"]["<<ii<<"]='"
                          <<data_mat[top_idx][ii]<<"'" <<std::endl;
      
      _I cur_seq_len = data_mat[top_idx][ii].size();
      if ( !seq_init && cur_seq_len != 0 ) seq_len = cur_seq_len;
      
      // Ensure all sequences are same length or zero::
      if ( cur_seq_len > 0 && cur_seq_len != seq_len ) {
        std::cerr 
        <<_fe<< "Failed during processing! Inconsitent Sequence Lengths;"
        <<"cur_seq_len='"<<cur_seq_len<<"' vs. seq_len='"<<seq_len<<"'\n"
        <<std::endl;
        exit (EXIT_FAILURE);
      }
      
      // Convert to 8 bit and push onto bit vector
      uint8_t* seq_bit = dna_to_bit( data_mat[top_idx][ii].c_str(), cur_seq_len );
      bit_vec.push_back(seq_bit);
      
      // if ( p1 ) std::cerr <<_fo<< "bit_size["<<top_idx<<"]["<<ii<<"]='"
      //                     <<sizeof(seq_bit)<<"'\n" <<std::endl;
      
      // Extra QC Check:: Skip 
      if ( false ) {
        // if ( p1 ) std::cerr <<_fo<< "seq_bit["<<top_idx<<"]["<<ii<<"]='"
        //                     <<seq_bit<<"'" <<std::endl;
        _c* bit_seq = bit_to_dna( seq_bit, data_mat[top_idx][ii].size() );
        if ( p1 ) std::cerr <<_fo<< "bit_seq["<<top_idx<<"]["<<ii<<"]='"
                            <<bit_seq<<"'\n" <<std::endl;
      }
    }
    
    if (bytes_read < GZ_LENGTH - 1) {
      if ( gzeof (gzin) ) {
        break;
      } else {
        const _c * error_string;
        error_string = gzerror (gzin, & err);
        if ( err ) {
          std::cerr 
          <<_fe<< "Failed during processing! gzerror["<<buf_cnt<<"]='"<<error_string<<"'\n"
          <<std::endl;
          exit (EXIT_FAILURE);
        }
      }
    }
    buf_cnt++;
  }
  gzclose (gzin);
  
  _I seq_cnt = bit_vec.size();
  _I bit_size = sizeof(bit_vec[0]);
  if ( p1 ) std::cerr <<_fo<< "Writing "<<seq_cnt<<" binary top sequences\n"
                      <<_fo<< "bit_size ='"<<bit_size<<"'\n"
                      <<std::endl;
  
  // Write full output vector
  gzFile gzo = gzopen( bgz_path.c_str(), "w");
  // Write sequence count and sequence length::
  //  TBD:: write dataType as well. This is done somewhere else as well, but
  //    we'll keep this simple for now...
  os.write( reinterpret_cast<_c *>(&seq_cnt), sizeof(unsigned) );
  os.write( reinterpret_cast<_c *>(&seq_len), sizeof(unsigned) );
  os.write( reinterpret_cast<_c *>(&bit_size), sizeof(unsigned) );
  
  // Write all the data::
  os.write( reinterpret_cast<_c*> (&bit_vec[0]),
            seq_cnt * bit_size );
  gzwrite( gzo, os.str().c_str(), os.str().size() );
  gzclose( gzo );
  
  return( success );
}

/*
 * load_cgntop_bgz:: load CG# (cgn by position, topSeq as uint8_t bits)
 *   TBD:: Expand function to template<dataType> instead of just uint8_t
 *   TBD:: Expand to be more general (maybe already done in other functions)
 *   
 */
_Y load_cgntop_bgz( const string& file,
                    std::vector<_S>& v,
                    
                    const _I vb = 0,
                    const _I vt = 1,
                    const _I tc = 0,
                    const _S ft="load_cgntop_bgz" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0  = vb > vt + 0;
  const _Y p1  = vb > vt + 1;
  const _Y p8  = vb > vt + 8;
  
  _Y success = true;
  
  if ( p0 ) std::cerr
    <<_fo<< "Starting..." << "\n"
    <<_fo<< "\t     file = '" << file << "'\n" 
    << std::endl;
  
  gzFile gzin = gzopen(file.c_str(),"r");
  
  _I seq_cnt;
  _i rv_cnt =  gzread( gzin, &seq_cnt, sizeof(unsigned) );
  if ( p1 ) std::cerr <<_fo<< "Sequence Counts='"<<seq_cnt<<"'" <<std::endl;
  
  _I seq_len;
  _i rv_len =  gzread( gzin, &seq_len, sizeof(unsigned) );
  if ( p1 ) std::cerr <<_fo<< "Sequence Length='"<<seq_len<<"'" <<std::endl;
  
  _I bit_size;
  _i rv_bit =  gzread( gzin, &bit_size, sizeof(unsigned) );
  if ( p1 ) std::cerr <<_fo<< "Bit Size='"<<bit_size<<"'" <<std::endl;
  
  std::vector<uint8_t*> bit_vec = init_bit_stack( seq_cnt, seq_len, vb,vt+10 );
  _i rv_dat = gzread( gzin, &bit_vec[0], seq_cnt*ceil(seq_len/bit_size) );
  gzclose(gzin);
  
  if ( p8 ) std::cerr <<_fo<< "rv_cnt='"<<rv_cnt<<"'\n"
                      <<_fo<< "rv_len='"<<rv_len<<"'\n"
                      <<_fo<< "rv_bit='"<<rv_bit<<"'\n"
                      <<_fo<< "rv_dat='"<<rv_dat<<"'\n"
                      <<std::endl;
  
  for( _I ii=0; ii<seq_cnt; ii++)
  {
    // uint8_t* seq_bit = &bit_vec[ii];
    // _c* bit_seq = bit_to_dna( seq_bit, seq_len );
    
    _c* bit_seq = bit_to_dna( bit_vec[ii], seq_len );
    
    if ( p1 ) 
      std::cerr <<_fo<< "bit_seq["<<ii<<"]='"<<bit_seq<<"'" <<std::endl;
  }
  
  return( success );
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * Internal Static c++ Functions::
 * 
 *                           CGN/TOP Stack Functions::
 * 
 * Internal::
 *     write_assoc_stack()
 *     load_assoc_stack()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

/*
 * Write Stack
 * 
 */

static _Y write_assoc_stack( const std::vector<_S> stack,
                             const _I idx_beg,
                             const _I idx_end,
                             const _I idx_cnt,
                             const _I seq_len,
                             
                             const _S file,
                             const _Y all = true,
                             const _c sep = 'b',
                             const _I vb = 0, 
                             const _I vt = 0 )
{
  const _S func_tag ("[write_assoc_stack]: ");
  const _S func_err ("[write_assoc_stack]: ERROR: ");
  const _S func_wrn ("[write_assoc_stack]: Warning: ");
  
  const _Y p1 = vb >= vt + 1;
  const _Y p2 = vb >= vt + 2;
  const _Y p3 = vb >= vt + 3;
  
  _Y success = true;
  const _S idx_file = file + ".idx";
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                             Check Files Exist::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  std::ofstream OUT_FH ( file, std::ios::binary | std::ios::out );
  if ( !OUT_FH.is_open() )
  {
    success = false;
    std::cerr << func_err << "Failed to open file: '" << file << "'\n" 
              << std::endl;
    return( success );
  }
  
  std::ofstream IDX_FH ( idx_file, std::ios::out );
  if ( !IDX_FH.is_open() )
  {
    OUT_FH.close();
    success = false;
    std::cerr << func_err << "Failed to open idx file: '" << idx_file << "'\n" 
              << std::endl;
    return( success );
  }
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                             Write Stack Data::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  _I stack_size = idx_end - idx_beg + 1;
  std::vector<uint8_t*> bit_vec = init_bit_stack( stack_size,
                                                  seq_len,
                                                  vb, vt+4 );
  std::vector<_I> pos_vec( stack_size );
  
  _c sep_val = '\t';
  if ( sep == 'c' ) sep_val = ',';
  
  if ( p1 ) std::cerr 
    << func_tag << "Writing Stack Data::\n"
    << func_tag << "\t    idx_beg = " << idx_beg << "\n"
    << func_tag << "\t    idx_end = " << idx_end << "\n"
    << func_tag << "\t    idx_cnt = " << idx_cnt << "\n"
    << func_tag << "\t    seq_len = " << seq_len << "\n"
    << func_tag << "\t        sep = " << sep << "\n"
    << func_tag << "\t    sep_val = '" << sep_val << "'\n"
    << func_tag << "\t full stack = '" << all << "'\n"
    << func_tag << "\t stack_size = " << stack_size << "\n"
    << func_tag << "\t stack file = " << file << "\n"
    << std::endl;
  
  _I bit_idx = 0;
  for ( _I all_idx = idx_beg; all_idx < stack_size; all_idx++ )
  {
    if ( stack[all_idx].size() != 0 )
    {
      uint8_t* cur_bit = dna_to_bit(
        stack[all_idx].c_str(),
        stack[all_idx].size() );
      // const _i bit_size = sizeof( cur_bit );
      
      if ( sep == 'b' &&  all ) bit_vec[all_idx] = cur_bit;
      if ( sep == 'b' && !all ) bit_vec[bit_idx] = cur_bit;
      if ( sep == 'b' && !all ) pos_vec[bit_idx] = all_idx;
      
      if ( sep != 'b' &&  all ) OUT_FH << stack[all_idx] << std::endl;
      if ( sep != 'b' && !all ) OUT_FH << all_idx << sep_val
                                       << stack[all_idx] << std::endl;
      bit_idx++;
    } else {
      if ( sep != 'b' &&  all ) OUT_FH << std::endl;
    }
  }
  
  if ( sep == 'b' )
  {
    if ( p2 ) std::cerr << func_tag << "\t Writing Binary Data..." << std::endl;
    
    if ( !all && bit_idx < bit_vec.size() ) bit_vec.resize( bit_idx );
    OUT_FH.write( reinterpret_cast<_c*> (&bit_vec[0]),
                  bit_vec.size() * sizeof(bit_vec[0]) );
    
    if ( p3 ) {
      // _I cgn_idx;
      _I bit_size = bit_vec.size();
      std::cerr << "Bit Size = " << bit_size << std::endl;
      for ( _I bit_idx = 0; bit_idx < bit_size; bit_idx++ )
      {
        // cgn_idx = bit_idx;
        // if ( !all ) cgn_idx = pos_vec[bit_idx];
        _I bit_len = sizeof( bit_idx );
        _c* bit_seq = bit_to_dna( bit_vec[bit_idx], seq_len );
        
        std::cerr 
          << func_tag << "\t bit_len[" << bit_idx << "] = " << bit_len << "\n"
          << func_tag << "\t bit_seq[" << bit_idx << "] = " << bit_seq << "\n"
          << std::endl;
      }
    }
    
    if ( !all && bit_idx < pos_vec.size() ) pos_vec.resize( bit_idx );
    if ( !all )
      OUT_FH.write( reinterpret_cast<_c*> (&pos_vec[0]), 
                    pos_vec.size() * sizeof(pos_vec[0]) );
    
    if ( p2 ) std::cerr << func_tag << "\t Wrote Binary Data!" << std::endl;
  }
  OUT_FH.close();
  if ( p1 ) std::cerr << func_tag << "Done writing Stack Data!\n" << std::endl;
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                             Write Index Data::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  if ( p1 ) std::cerr 
    << func_tag << "Writing Stack Index::\n"
    << func_tag << "\t index file = " << idx_file << "\n"
    << std::endl;
  
  IDX_FH << idx_beg << '\t'
         << idx_end << '\t'
         << idx_cnt << '\t'
         << seq_len << '\t'
         << sep << '\t'
         << all << '\t'
         << std::endl;
  IDX_FH.close();
  if ( p1 ) std::cerr << func_tag << "Done Writing Stack Index!\n" << std::endl;
  
  return( success );
};

/*
 * load_stack
 *   - binary + index
 *   - csv/tsv (index stack/full stack)
 * 
 */
static std::vector<_S> load_assoc_stack( const _S file,
                                         const _I vb = 0, 
                                         const _I vt = 0 )
{
  const _S func_tag ("[load_assoc_stack]: ");
  const _S func_err ("[load_assoc_stack]: ERROR: ");
  const _S func_wrn ("[load_assoc_stack]: Warning: ");
  
  const _Y p1 = vb >= vt + 1;
  const _Y p2 = vb >= vt + 2;
  
  const _S idx_file = file + ".idx";
  std::vector<_S> nil_stack;
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                             Check Files Exist::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  std::ifstream DAT_FH ( file, std::ios::binary | std::ios::in );
  if ( !DAT_FH.is_open() )
  {
    std::cerr << func_err << "Failed to open file: '" << file << "'\n" 
              << std::endl;
    return( nil_stack );
  }
  
  std::ifstream IDX_FH ( idx_file, std::ios::in );
  if ( !IDX_FH.is_open() )
  {
    DAT_FH.close();
    std::cerr << func_err << "Failed to open idx file: '" << idx_file << "'\n" 
              << std::endl;
    return( nil_stack );
  }
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                              Load Index Data::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  if ( p1 ) std::cerr 
    << func_tag << "Loading Stack Index::\n"
    << func_tag << "\t index file = " << idx_file << "\n"
    << std::endl;
  
  _S line;
  // std::vector<_S> lines;
  std::vector<_S> idx_vec;
  getline( IDX_FH, line );
  IDX_FH.close();
  
  parse_line( line, '\t', idx_vec );
  const _I idx_beg = stoi( idx_vec[0] );
  const _I idx_end = stoi( idx_vec[1] );
  const _I idx_cnt = stoi( idx_vec[2] );
  const _I seq_len = stoi( idx_vec[3] );
  const _c sep = (_c) idx_vec[4][0];
  _c sep_val = '\t';
  if ( sep == 'c') sep_val = ',';
  const _I all_int = stoi( idx_vec[5] );
  _Y all = false;
  if ( all_int > 0 ) all = true;
  
  _t dna_bytes = (seq_len / 4) + (seq_len % 4 != 0);
  uint8_t* zero_bit = new uint8_t[dna_bytes];
  std::memset(zero_bit, 0, dna_bytes);
  _S zero_seq ( bit_to_dna(zero_bit, seq_len) );
  
  _I stack_size = idx_cnt;
  if ( all ) stack_size = idx_end - idx_beg + 1;
  
  std::vector<uint8_t*> bit_vec = init_bit_stack( stack_size,
                                                  seq_len,
                                                  vb, vt+3 );
  std::vector<_I> pos_vec( stack_size );
  std::vector<_S> cgn_stack( stack_size );
  
  if ( p1 ) std::cerr 
    << func_tag << "Loaded Stack Index::\n"
    << func_tag << "\t    idx_beg = " << idx_beg << "\n"
    << func_tag << "\t    idx_end = " << idx_end << "\n"
    << func_tag << "\t    idx_cnt = " << idx_cnt << "\n"
    << func_tag << "\t    seq_len = " << seq_len << "\n"
    << func_tag << "\t        sep = " << sep << "\n"
    << func_tag << "\t    sep_val = '" << sep_val << "'\n"
    << func_tag << "\t    all_int = '" << all_int << "'\n"
    << func_tag << "\t full stack = '" << all << "'\n"
    << func_tag << "\t stack_size = " << stack_size << "\n"
    << func_tag << "\t stack file = " << file << "\n"
    << std::endl;
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                              Load Stack Data::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  if ( p1 ) std::cerr << func_tag << "Loading Stack Data!\n" << std::endl;
  
  if ( sep == 'b' )
  {
    if ( p2 ) std::cerr << func_tag << "\t Loading Binary Data..." << std::endl;
    
    DAT_FH.read( reinterpret_cast<_c*> (&bit_vec[0]), 
                 stack_size*sizeof(bit_vec[0]) );
    
    if ( !all )
      DAT_FH.read( reinterpret_cast<_c*> (&pos_vec[0]), 
                   stack_size*sizeof(pos_vec[0]) );
    
    if ( p2 ) std::cerr << func_tag << "\t Loaded Binary Data!\n" << std::endl;
    
    /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
     * 
     *                         Converting Binary to DNA::
     * 
     * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
    
    _I cgn_idx = 0;
    _I pos_size = 0;
    _I bit_size = bit_vec.size();
    if ( !all ) pos_size = pos_vec.size();
    
    if ( p2 ) std::cerr 
      << func_tag << "\t Converting Binary Data...\n" 
      << func_tag << "\t\t bit_size = " << bit_size << "\n"
      << func_tag << "\t\t pos_size = " << pos_size << "\n"
      << std::endl;
    
    /*
     * 
     * Write it with or without cgn_idx (all) consideration...
     * 
     */
    for ( _I bit_idx = 0; bit_idx < bit_size; bit_idx++ )
    {
      
      cgn_idx = bit_idx;
      if ( !all ) cgn_idx = pos_vec[bit_idx];
      
      if ( p2 ) std::cerr 
        << func_tag << "\t  all_val[" << bit_idx << "] = " << all << "\n"
        << func_tag << "\t  bit_idx[" << bit_idx << "] = " << bit_idx << "\n"
        << func_tag << "\t  cgn_idx[" << cgn_idx << "] = " << cgn_idx << "\n"
        << std::endl;
      
      if ( !all ) cgn_idx = pos_vec[cgn_idx];
      
      _c* bit_seq = bit_to_dna( bit_vec[bit_idx], seq_len );
      cgn_stack[cgn_idx] = (_S) bit_seq;
      
      // if ( bit_vec[bit_idx] == zero_bit ) cgn_stack[cgn_idx] = "";
      if ( cgn_stack[cgn_idx].compare( zero_seq ) == 0 ) 
        cgn_stack[cgn_idx] = "";
      
      if ( p2 ) std::cerr 
        << func_tag << "\t  bit_seq[" << bit_idx << "] = " << bit_seq << "\n"
        << func_tag << "\t  cgn_seq[" << cgn_idx << "] = " << cgn_stack[cgn_idx] << "\n"
        << func_tag << "\t zero_seq[" << cgn_idx << "] = " << zero_seq << "\n"
        << std::endl;
    }
    
    if ( p2 ) std::cerr 
      << func_tag << "\t Converted Binary Data...\n" << std::endl;
  }
  else
  {
    _Y gzipped = is_gzipped( file );
    
    _S c_str;
    if ( gzipped ) {
      _S gzip_cmd = "gzip -dc " + file;
      c_str = exec_cmd( gzip_cmd, vb, vt+1 );
    } else {
      c_str = read_file_fast( file, vb, vt+1 );
    }
    
    _I contents_size = c_str.size();
    if ( p1 ) std::cerr
      << func_tag << "Contents Size = '" << contents_size << "'\n"
      << func_tag << "         File = '" << file << "'\n"
      << std::endl;
    
    if ( all )
    {
      parse_line( c_str, '\n', cgn_stack );
      
      _I cgn_stack_size = (_I) cgn_stack.size();
      if ( p1 ) std::cerr
        << func_tag << "Cgn Stack Size = " << cgn_stack_size
        << std::endl;
      
      if ( p2 )
      {
        for ( _I row = 0; row < cgn_stack_size; row++ )
        {
          std::cerr
          << func_tag << "\t c_str[" << row << "] = '" 
          << cgn_stack[row]
          << std::endl;
        }
      }
    }
    else
    {
      // Parse contents into Lines::
      std::vector<_S> lines;
      parse_line( c_str, '\n', lines );
      _I lines_size = lines.size();
      
      if ( p1 ) std::cerr
        << func_tag << "Lines Size = " << lines_size << "\n"
        << std::endl;
      
      // Parse Columns into Vectors
      // std::vector<std::vector< _S>> column( 2 );
      // std::vector< std::vector<_S> > two_cols( lines_size );
      std::vector< std::vector<_S> > two_cols( stack_size );
      
      for ( _I row = 0; row < lines_size - 1; row++ )
      {
        parse_line( lines[row], sep_val, two_cols[row] );
        
        _S val0 = two_cols[row][0];
        _S val1 = two_cols[row][1];
        
        _I cgn_idx = (_I) stoi( val0 );
        cgn_stack[cgn_idx] = "";
        cgn_stack[cgn_idx] = val1;
        
        _I cgn_len = cgn_stack.size();
        if ( p2 ) std::cerr
          << func_tag << "\t cgn_seq[" << val0 << ", " << cgn_idx << "] = '" 
          << val1 << "'\t cgn_len = '" << cgn_len << "'; \t cgn_stack = " 
          << cgn_stack[cgn_idx]
          << std::endl;
      }
    }
  }
  DAT_FH.close();
  
  if ( p1 ) std::cerr << func_tag << "Done Loading Stack Data!\n" << std::endl;
  
  if ( p2 )
  {
    _I cgn_size = cgn_stack.size();
    std::cerr << "CGN STACK:: size = " << cgn_size << std::endl;
    for ( _I cgn_idx = 0; cgn_idx < cgn_size; cgn_idx++ )
    {
      std::cerr
      << func_tag << "\t cgn_seq[" << cgn_idx << "] = " << cgn_stack[cgn_idx] << "\n"
      << std::endl;
    }
  }
  
  delete [] zero_bit;
  
  return( cgn_stack );
};

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * Internal c++ Functions::
 * 
 *                         Mutate Sequence Functions::
 *                         
 *                              Reverse/Complement
 *                            Bisulfite Conversion
 *                               De-methylation
 * 
 * Internal::
 *   convert()
 *   str_to_lower()
 *   str_to_upper()
 *   mutate_seq()
 *   is_top_bot()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

static void replace_nuc_at( _S& s,
                            const _I i,
                            _S  r,
                            _S  o = "",
                            const _Y uc = false )
{
  const _S func_tag ("[replace_nuc_at]: ");
  const _S func_err ("[replace_nuc_at]: ERROR: ");
  const _S func_wrn ("[replace_nuc_at]: Warning: ");
  
  const _I seq_len = s.size();
  const _I org_len = o.size();
  
  _c r_nuc = r[0];
  if ( uc ) r_nuc = std::toupper( r_nuc );
  
  if ( i < seq_len ) 
  {
    if ( org_len > 0 ) {
      
      _S seq_nuc( 1, s[i] );
      for (auto & c: o) c = toupper(c);
      for (auto & c: seq_nuc ) c = toupper(c);
      
      // Old Code::
      // ::str_to_upper( o );
      // ::str_to_upper( seq_nuc );
      
      if ( o.compare(seq_nuc) != 0 )
      {
        std::cerr
        << func_wrn << "Original Nucleotide Does Not Match::\n"
        << func_wrn << "\t Index = " << i << "\n"
        << func_wrn << "\t   Seq = " << seq_nuc << "\n"
        << func_wrn << "\t  User = " << o << "\n"
        << func_wrn << "\t  User = " << r << "\n"
        << std::endl;
      }
    }
    s[i] = r_nuc;
  }
}

static void mutate_seq( _S& s,
                        const _S m,
                        const _Y uc = false,
                        const _I vb = 0, 
                        const _I vt = 1 )
{
  
  const _I seq_len = s.size();
  const _I end_idx = seq_len - 1;
  const _I mod_len = m.size();
  
  // Loop through mutations tasks
  for ( _I m_idx = 0; m_idx < mod_len; m_idx++ ) {
    /*
     *  Reverse()
     */
    if ( m[m_idx] == 'R' || m[m_idx] == 'r' ) std::reverse(s.begin(), s.end());
    
    /*
     *  Complement()
     */
    if ( m[m_idx] == 'C' || m[m_idx] == 'c' ) {
      if ( !init_rvcp ) revcomp_init();
      for ( _I ii=0; ii < seq_len; ii++ )
        s[ii] = comp_vec[ (_B) s[ii] ];
    }
    
    /*
     *  De-Methylate()
     */
    if ( m[m_idx] == 'X' || m[m_idx] == 'x' ) {
      if ( !init_bscX ) bscX_init();
      for ( _I ii=0; ii < seq_len; ii++ )
        s[ii] = bscX_vec[ (_I) s[ii] ];
    }
    
    /*
     *  Bisulfite-Convert-Unmethylated()
     */
    if ( m[m_idx] == 'U' || m[m_idx] == 'u' ) {
      if ( !init_bscU ) bscU_init();
      for ( _I ii=0; ii < seq_len; ii++ )
        s[ii] = bscU_vec[ (_I) s[ii] ];
    }
    
    /*
     *  Bisulfite-Convert-Methylated()
     */
    if ( m[m_idx] == 'M' || m[m_idx] == 'm' ) {
      if ( !init_bscM ) bscM_init();
      for ( _I ii=0; ii<end_idx; ii++ ) {
        s[ii+0] = bscM1_vec[ (_I) toupper(s[ii]) ][ (_I) toupper(s[ii+1]) ][0];
        s[ii+1] = bscM1_vec[ (_I) toupper(s[ii]) ][ (_I) toupper(s[ii+1]) ][1];
        s[ii+0] = bscM2_vec[ (_I) s[ii] ];
      }
      s[end_idx] = bscM2_vec[ (_I) s[end_idx] ];
    }
    
    /*
     *  Bisulfite-Convert-Degenerate()
     */
    if (m[m_idx] == 'D' || m[m_idx] == 'd') {
      if ( !init_bscD ) bscD_init();
      for ( _I ii=0; ii < end_idx; ii++ ) {
        s[ii+0] = bscD1_vec[ (_I) toupper(s[ii]) ][ (_I) toupper(s[ii+1]) ][0];
        s[ii+1] = bscD1_vec[ (_I) toupper(s[ii]) ][ (_I) toupper(s[ii+1]) ][1];
        s[ii+0] = bscD2_vec[ (_I) s[ii] ];
      }
      s[end_idx] = bscD2_vec[ (_I) s[end_idx] ];
    }
    
    if ( uc || (_I) m[m_idx] < 91 )
      for ( _I ii=0; ii < seq_len; ii++ ) s[ii] = std::toupper( s[ii] );
  }
  
};

_Y parse_probe_id( _S& s,
                   _S& cgn,
                   _B& tbs,
                   _B& cos,
                   _I& inf,
                   _B& ale,
                   const unsigned sep='_',
                   const _I vb = 0,
                   const _I vt = 1,
                   const _I tc = 0,
                   const _S ft="parse_probe_id" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  // const _Y p0  = vb > vt + 0;
  _Y success = true;
  
  // Clear data::
  cgn = "";
  tbs = '0';
  cos = '0';
  inf = 0;
  ale = '0';
  
  // Split data
  std::vector<_S> prb_vec;
  parse_line( s, sep, prb_vec );
  
  // Validate Split::
  if ( prb_vec.size() < 3 ) {
    std::cerr <<_fe<< "Probe_ID='"<<s<<"' vector length("<<prb_vec.size()
              <<") < 3\n"<<std::endl;
    return( false );
  }
  
  // Validate Strands::
  if ( prb_vec[1].size() < 3 ) {
    std::cerr <<_fe<< "Probe_ID='"<<s<<"' std='"<<prb_vec[1]<<"' length("
              <<prb_vec[1].size()<<") < 3\n"<<std::endl;
    return( false );
  }
  
  // Validate Allele::
  if ( prb_vec[2].size() < 1 ) {
    std::cerr <<_fe<< "Probe_ID='"<<s<<"' allele='"<<prb_vec[2]<<"' length("
              <<prb_vec[2].size()<<") < 1\n"<<std::endl;
    return( false );
  }
  
  // Update Assignments::
  cgn = prb_vec[0];
  
  if ( prb_vec[1][0] == 'T' || prb_vec[1][0] == 'B' ) tbs = prb_vec[1][0];
  // else tb = '0';
  
  if ( prb_vec[1][1] == 'C' || prb_vec[1][1] == 'O' ) cos = prb_vec[1][1];
  // else tb = '0';
  
  if ( prb_vec[1][2] == '1' ) inf = 1;
  else if ( prb_vec[1][2] == '2' ) inf = 2;
  // else inf = 0;
  
  ale = prb_vec[2][0];
  
  return( success );
}

static _i is_top_bot( const _S& s, 
                      const _i var_len = 2,
                      const _Y uc = false,
                      const _I vb = 0,
                      const _I vt = 1,
                      const _I tc = 0,
                      const _S ft="is_top_bot" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0  = vb > vt + 0;
  const _Y p1  = vb > vt + 1;
  const _Y p2  = vb > vt + 2;
  
  if ( !init_topbot ) topbot_init();
  
  _i tb_ret = 0;
  _i seq_len = s.size();
  _i mid_len = (_i) ( seq_len/2 ) - 1;
  _i mid_idx = seq_len % 2;
  
  if ( seq_len - var_len < 3 && p1 ) {
    std::cerr
    <<_fw<< "\nSequence = " << s  << "\n"
    <<_fw<< "\t  Sequence is too short (less than three base pairs"
    <<_fw<< "after var_len removal) to calculate top/bot: seq len = "
    << seq_len << "\n"
    << std::endl;
    
    return( tb_ret );
  }
  
  if ( p0 ) std::cerr
    <<_fo<< "Sequence = " << s << "\n"
    <<_fo<< "\t  len = " << seq_len << "\n"
    <<_fo<< "\t  mid = " << mid_idx << "\n"
    << std::endl;
  
  if ( var_len == 1 ) mid_len += seq_len % 2;
  
  _i ii = mid_len - 1;
  _i jj = mid_len + var_len;
  
  if ( p1 ) {
    std::cerr <<_fo<< "\tMid("<< mid_len <<") = [" << s[mid_len];
    if (var_len == 2)
      std::cerr << s[mid_len+1] << "] (" << mid_len+1 << ")\n" << std::endl;
    else
      std::cerr << "]\n" << std::endl;
  }
  
  while( ii >= 0 && jj < seq_len ) {
    tb_ret = topbot_vec[ (_I) s[ii] ][ (_I) s[jj] ];
    
    if ( p2 ) std::cerr
      <<_fo<<"TopBot = '"<<tb_ret<<"'"<<RET
      <<_fo<<TAB2<<"  ups["<<ii<<"] = '"<<s[ii]<<"'"<<RET
      <<_fo<<TAB2<<"  dns["<<jj<<"] = '"<<s[jj]<<"'"<<RET
      << std::endl;
    
    if ( tb_ret ) break;
    ii--;
    jj++;
  }
  
  return( tb_ret );
};

static void expand_seq_vec( std::vector<_S>& seqs, 
                            const std::vector<_c>& nucs,
                            _t pos_idx )
{
  const _I rep_cnt = nucs.size();
  _S base_seq = seqs[0];
  
  std::vector<_S> seqs_cp = seqs;
  for (_I ii = 0; ii < rep_cnt - 1; ii++)
  {
    seqs.insert( seqs.end(), seqs_cp.begin(), seqs_cp.end() );
  }
  std::sort( seqs.begin(), seqs.end() );
  
  // Now re-loop on the new seqs mutating via modulo::
  const _I new_cnt = seqs.size();
  for (_I ii = 0; ii < new_cnt; ii++)
  {
    _I idx = ii % rep_cnt;
    seqs[ii][pos_idx] = nucs[idx];
  }
};

static std::vector<_S> expand_iupac_seq( const _S& seq,
                                         const _I vb = 0,
                                         const _I vt = 1,
                                         const _I tc = 0,
                                         const _S ft="expand_iupac_seq" ) {
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p2  = vb > vt + 2;
  const _I seq_len = seq.size();
  
  if ( !init_exp_iupac ) iupac_exp_init();
  
  std::vector<_S> exps_vec;
  exps_vec.push_back(seq);
  _I rep_num = 1;
  
  // Return if sequence length is zero
  if ( seq_len == 0 ) return( exps_vec );
  
  for ( _t ii=0; ii<seq_len; ii++) {
    _t exp_len = iupac_exp_len[ (_B) seq[ii] ];
    
    if ( exp_len > 1 ) {
      rep_num = rep_num * exp_len;
      
      // Should really search first with find() and throw error if not found..
      if ( iupac_exp_nuc.find( seq[ii] ) == iupac_exp_nuc.end() ) {
        
      }
      std::vector<_c> rep_nucs = iupac_exp_nuc[ (_c) seq[ii] ];
      expand_seq_vec( exps_vec, rep_nucs, ii );
      
      if ( p2 ) std::cerr
        <<_fo<<"Exp_Len["<<ii<<"] = '"<<exp_len<<"'"<<RET
        <<_fo<<"Rep_Num["<<ii<<"] = '"<<rep_num<<"'"<<RET
        << std::endl;
    }
  }
  // Remove any possible duplicates and ensure/enforce alphabetical order::
  std::sort( exps_vec.begin(), exps_vec.end() );
  
  vector<_S>::iterator ip;
  ip = std::unique( exps_vec.begin(), exps_vec.end() );
  exps_vec.resize( std::distance(exps_vec.begin(), ip) );
  
  return( exps_vec );
};

static _i improbe_seq_umd( std::vector<std::vector<_S>> & vec,
                           _S bsc_seq,
                           const _I vec_idx,
                           _I key_idx,
                           const _i nxb_pos,
                           const _i inf_num,
                           const _c bc_char,
                           const _c co_char,
                           std::map< _S, std::map< _S, std::map< _S, std::vector<_I> > > > & tri_maps,
                           // std::map< _S, _S > & col_map,
                           
                           const _I prb_len = 50,
                           const _I vb = 0,
                           const _I vt = 1,
                           const _I tc = 0 )
{
  const _S tb(tc, '\t');
  const _S func_tag("[improbe_seq_umd]: "+tb);
  const _S func_err("[improbe_seq_umd]: ERROR: ");
  const _S func_wrn("[improbe_seq_umd]: "+tb+"Warning: ");
  
  const _Y p1 = vb >= vt + 1;
  
  /*
   * Sanity Check::
   * 
   */
  const _I min_len = 2 + ( prb_len * 2 );
  const _I seq_len = bsc_seq.size();
  if ( seq_len < min_len ||
       nxb_pos < 0 || 
       ( inf_num != 1 && inf_num != 2 && inf_num != 3) )
  {
    std::cerr
    << func_err << "Invalid Indicie Values! Exiting..." << "\n"
    << func_err << "\t seq_len = " << seq_len << "\n"
    << func_err << "\t min_len = " << min_len << "\n"
    << func_err << "\t nxb_pos = " << nxb_pos << "\n"
    << func_err << "\t inf_num = " << inf_num << "\n"
    << std::endl;
    
    return( -1 );
  }
  
  _c bc_char_lc = std::tolower(bc_char);
  const _S rc_str_lc = "rc";
  _S bc_str_lc;
  bc_str_lc = bc_char_lc;
  
  const _c BC = std::toupper(bc_char);
  const _c CO = std::toupper(co_char);
  
  if ( BC!='N' ) mutate_seq( bsc_seq, bc_str_lc, false, vb, vt+6 );
  if ( CO=='O' ) mutate_seq( bsc_seq, rc_str_lc, false, vb, vt+6 );
  
  _i cpg_pos = -1;
  _i sec_pos = -1;
  _i bod_pos = -1;
  _i end_pos = -1;
  
  cpg_pos = nxb_pos + 1;
  sec_pos = cpg_pos + 1;
  bod_pos = sec_pos + 1;
  end_pos = bod_pos + (prb_len - 2);
  
  /*
   * Extract Subsequences and Make Assignments::
   * 
   */
  _I org_idx = key_idx;
  vec[key_idx++][vec_idx] = bsc_seq;
  vec[key_idx++][vec_idx] = bsc_seq.substr( nxb_pos, 1 );
  vec[key_idx++][vec_idx] = bsc_seq.substr( cpg_pos, 1 );
  vec[key_idx++][vec_idx] = bsc_seq.substr( sec_pos, 1 );
  vec[key_idx++][vec_idx] = bsc_seq.substr( bod_pos, end_pos-bod_pos );
  vec[key_idx++][vec_idx] = bsc_seq.substr( end_pos, 1 );
  
  if ( inf_num==1 ) 
    vec[key_idx][vec_idx] = vec[key_idx-4][vec_idx] + vec[key_idx-3][vec_idx] + vec[key_idx-2][vec_idx];
  if ( inf_num==2)
    vec[key_idx][vec_idx] = vec[key_idx-3][vec_idx] + vec[key_idx-2][vec_idx] + vec[key_idx-1][vec_idx];
  if ( inf_num==3)
    vec[key_idx][vec_idx] = vec[key_idx-5][vec_idx] + vec[key_idx-4][vec_idx] + vec[key_idx-3][vec_idx] + vec[key_idx-2][vec_idx] + vec[key_idx-1][vec_idx];
  
  _S prb_seq = "";
  if ( inf_num==1 ) prb_seq = bsc_seq.substr( cpg_pos, prb_len );
  if ( inf_num==2 ) prb_seq = bsc_seq.substr( sec_pos, prb_len );
  if ( inf_num==3 ) prb_seq = bsc_seq.substr( nxb_pos, prb_len+2 );
  
  if ( vec[key_idx][vec_idx].compare( prb_seq ) != 0 )
  {
    std::cerr
    << func_err << "Failed to match Probe Seuqnces! Exiting..." << "\n"
    << func_err << "\t key_idx = " << key_idx << "\n"
    << func_err << "\t vec_idx = " << vec_idx << "\n"
    << func_err << "\t inf_num = " << inf_num << "\n"
    << func_err << "\t prb_seq = " << prb_seq << "\n"
    << func_err << "\t key_seq = " << vec[key_idx][vec_idx] << "\n"
    << func_err << "\n"
    << func_err << "\t bsc_seq = " << vec[org_idx+0][vec_idx] << "\n"
    << func_err << "\t nxb_seq = " << vec[org_idx+1][vec_idx] << "\n"
    << func_err << "\t cpg_seq = " << vec[org_idx+2][vec_idx] << "\n"
    << func_err << "\t sec_seq = " << vec[org_idx+3][vec_idx] << "\n"
    << func_err << "\t bod_seq = " << vec[org_idx+4][vec_idx] << "\n"
    << func_err << "\t end_seq = " << vec[org_idx+5][vec_idx] << "\n"
    << func_err << "\t key_seq = " << vec[org_idx+6][vec_idx] << "\n"
    << func_err << "\t prb_seq = " << prb_seq << "\n"
    << std::endl;
    
    return( -1 );
  } else {
    if ( p1 ) std::cerr
      << func_tag << "Passed Probe Matching!\n"
      << func_tag << "\t prb_seq = " << prb_seq << "\n"
      << func_tag << "\t key_seq = " << vec[key_idx][vec_idx] << "\n"
      << std::endl;
  }
  key_idx++;
  
  /*
   * Reverse Complement (5'->3') all Sequence Data to Match Probe Order
   * 
   */
  for ( _t ii = org_idx+1; ii < key_idx; ii++ )
  {
    mutate_seq( vec[ii][vec_idx], rc_str_lc, false, vb, vt+6 );
  }
  
  /*
   * Add Tri-Seq and 46-mer to Tri-Map::
   * 
   */
  _S tri_seq = prb_seq.substr( 0,  3 );
  _S n46_seq = prb_seq.substr( 3, 46 );
  
  /*
   const _Y success_insert = insert_tri_map( tri_maps, tri_seq, n46_seq, 
   bc_str_lc, vec_idx, 
   vb+99, vt+1, tc+1 );
   if ( !success_insert )
   {
   std::cerr
   << func_err << "Failed to Insert Data into Tri-Map!!!" << "\n"
   << func_err << "\t tri_seq = " << tri_seq << "\n"
   << func_err << "\t n46_seq = " << n46_seq << "\n"
   << func_err << "\t bsc_str = " << bc_str_lc << "\n"
   << func_err << "\t vec_idx = " << vec_idx << "\n"
   << std::endl;
   
   return( -1 );
   }
   */
  
  return( key_idx );
};

static _Y improbe_seq( std::vector<std::vector<_S>> & vec,
                       const _S ids_str,
                       const _S fwd_seq,
                       const _I vec_idx,
                       std::vector<_S> & names,
                       std::map< _S, std::map< _S, std::map< _S, std::vector<_I> > > > & tri_maps,
                       // std::map< _S, _S > & col_map,
                       
                       _S frs_str,
                       _S cos_str,
                       _S bsc_str,
                       _S din_str,
                       
                       const _I prb_len = 50,
                       const _I vb = 0,
                       const _I vt = 1,
                       const _I tc = 0 )
{
  const _S tb(tc, '\t');
  const _S func_tag("[improbe_seq]: "+tb);
  const _S func_err("[improbe_seq]: ERROR: ");
  const _S func_wrn("[improbe_seq]: "+tb+"Warning: ");
  
  const _Y p1 = vb >= vt + 1;
  const _Y p2 = vb >= vt + 2;
  const _Y p8 = vb >= vt + 8;
  
  const _I ret_size = names.size();
  
  if ( p1 ) std::cerr
    << func_tag << "Starting..." << "\n"
    << func_tag << "\t  vec_idx = " << vec_idx << "\n" 
    << func_tag << "\t  fwd_seq = " << fwd_seq << "\n" 
    << func_tag << "\t  frs_str = " << frs_str << "\n"
    << func_tag << "\t  cos_str = " << cos_str << "\n"
    << func_tag << "\t  bsc_str = " << bsc_str << "\n"
    << func_tag << "\t  din_str = " << din_str << "\n"
    << func_tag << "\t ret_size = " << ret_size << "\n"
    << std::endl;
  
  /*
   * 
   * Update/Modify Input Variables::
   * 
   */
  // Resize/Clean vec
  if ( vec.size() < ret_size ) vec.resize( ret_size );
  // for ( _t ii = 0; ii < ret_size; ii++ ) vec[ii][vec_idx].clear();
  
  const _I min_len = 2 + ( prb_len * 2 );
  const _I fwd_len = fwd_seq.size();
  const _I din_idx = (_i) ( fwd_len / 2 ) - 1;
  
  const _c FR = str_to_char_idx( frs_str, true );
  const _c CO = str_to_char_idx( cos_str, true );
  const _c D0 = str_to_char_idx( din_str, true, false, 0 );
  const _c D1 = str_to_char_idx( din_str, true, false, 1 );
  
  const _c D0_FWD = str_to_char_idx( fwd_seq, true, din_idx );
  const _c D1_FWD = str_to_char_idx( fwd_seq, true, din_idx+1 );
  
  /*
   * 
   * Input Sanity Checks::
   * 
   */
  if ( D0=='C' && D1=='G' && D0_FWD=='C' && D1_FWD=='G' )
  {
    std::cerr
    << func_err << "Failed to match CG! Exiting..." << "\n"
    << func_err << "\t din_fwd = " << D0_FWD << D1_FWD << "\n"
    << func_err << "\t din_str = " << D0 << D1 << "\n"
    << std::endl;
    
    return( false );
  }
  
  if ( fwd_len < min_len )
  {
    std::cerr
    << func_err << "Failed Min Template Sequence Length! Exiting..." << "\n"
    << func_err << "\t fwd_len = " << fwd_len << "\n"
    << func_err << "\t min_len = " << min_len << "\n"
    << std::endl;
    
    return( false );
  }
  
  /*
   * 
   * Perform Strand Flipping & Bisulfite Conversion::
   * 
   */
  _S bsc_seq = fwd_seq;
  if ( FR=='R' ) mutate_seq( bsc_seq, "rc", false, vb, vt+1 );
  
  /*
   * Make Input Data Assignments::
   * 
   */
  _i key_idx = 0;
  vec[key_idx++][vec_idx] = ids_str;
  vec[key_idx++][vec_idx] = din_str;
  vec[key_idx++][vec_idx] = frs_str;
  vec[key_idx++][vec_idx] = cos_str;
  vec[key_idx++][vec_idx] = bsc_str;
  vec[key_idx++][vec_idx] = fwd_seq;
  
  /*
   * 
   * Set Substring Indicies::
   * 
   */
  _i nxb_pos = -1;
  
  // Default Next Base Offset for RS#
  if ( D0=='R' && D1=='S' && FR=='F' && CO=='C' ) nxb_pos = din_idx - 1; // 60;
  if ( D0=='R' && D1=='S' && FR=='R' && CO=='C' ) nxb_pos = din_idx + 0; // 61;
  if ( D0=='R' && D1=='S' && FR=='F' && CO=='O' ) nxb_pos = din_idx + 0; // 61;
  if ( D0=='R' && D1=='S' && FR=='R' && CO=='O' ) nxb_pos = din_idx - 1; // 60;
  
  // Default Next Base Offset for CH#
  if ( D0=='C' && D1=='H' && FR=='F' && CO=='C' ) nxb_pos = din_idx - 1; // 60;
  if ( D0=='C' && D1=='H' && FR=='R' && CO=='C' ) nxb_pos = din_idx + 0; // 61;
  if ( D0=='C' && D1=='H' && FR=='F' && CO=='O' ) nxb_pos = din_idx + 0; // 61;
  if ( D0=='C' && D1=='H' && FR=='R' && CO=='O' ) nxb_pos = din_idx - 1; // 60;
  
  // Default Next Base Offset for CG#
  if ( D0=='C' && D1=='G' && CO=='C') nxb_pos = din_idx - 1; // 60;
  if ( D0=='C' && D1=='G' && CO=='O') nxb_pos = din_idx + 0; // 61;
  
  if ( nxb_pos == -1 )
  {
    std::cerr << func_err << "Unsupported DiNucleotide Class: '" 
              << D0 << D1 << "'\n" << std::endl;
    return( false );
  }
  
  _i pre_idx = key_idx;
  for ( _c BC : bsc_str )
  {
    _i inf_num = 0;
    if ( BC=='U' || BC=='u' ) inf_num = 1;
    if ( BC=='M' || BC=='m' ) inf_num = 1;
    if ( BC=='D' || BC=='d' ) inf_num = 2;
    if ( BC=='N' || BC=='n' ) inf_num = 3;
    
    if ( p8 ) std::cerr
      << func_tag << "Pre improbe_seq_umd()" << "\n"
      << func_tag << "\t      BC = '" << BC << "'\n"
      << func_tag << "\t inf_num = '" << inf_num << "'\n"
      << func_tag << "\t pre_idx = '" << pre_idx << "'\n"
      << std::endl;
    
    key_idx = improbe_seq_umd( vec, bsc_seq, 
                               vec_idx, key_idx, nxb_pos, inf_num, 
                               BC, CO, tri_maps, prb_len,
                               vb, vt+1, tc+1 );
    
    if ( key_idx <= pre_idx )
    {
      std::cerr
      << func_err << "Inconsistent Previous and New Key Index! Exiting..." << "\n"
      << func_err << "\t pre_idx = " << pre_idx << "\n"
      << func_err << "\t key_idx = " << key_idx << "\n"
      << std::endl;
      
      return( false );
    }
    pre_idx = key_idx;
  }
  
  /*
   * 
   * Verbose Output::
   * 
   */
  if ( p2 )
  {
    
    _S pre_str;
    pre_str = vec[0][vec_idx] +
      ": PR=" + vec[1][vec_idx] +
      ": FR=" + vec[2][vec_idx] +
      ", CO=" + vec[3][vec_idx] +
      ", BC=" + vec[4][vec_idx] + ": ";
    
    _I vec_size = pre_idx;
    for ( _t ii = 0; ii < vec_size; ii++ )
    {
      std::cerr << func_tag << pre_str << names[ii] << "["<<ii<<"] = " << vec[ii][vec_idx]
                << std::endl;
    }
    std::cerr << "\nDone.\n" << std::endl;
  }
  
  return( true );
};

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * String/Vector/Matrix Manipulation::
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

std::vector<_i> strVec_to_intVec( const std::vector<_S>& vec )
{
  std::vector<_i> ivec;
  _I v_size = vec.size();
  if ( v_size == 0 ) return ivec;
  
  ivec.clear();
  ivec.resize(v_size);
  for (_I ii=0; ii < v_size; ++ii) ivec[ii] = std::stoi(vec[ii]);
  
  return ivec;
};

std::vector<_I> strVec_to_uintVec( const std::vector<_S>& vec )
{
  std::vector<_I> ivec;
  _I v_size = vec.size();
  if ( v_size == 0 ) return ivec;
  
  ivec.clear();
  ivec.resize(v_size);
  for (_I ii=0; ii < v_size; ++ii) ivec[ii] = std::stoul(vec[ii]);
  
  return ivec;
};

// Parse flat vector( elements = nrows * ncols ) into a 2d matrix[nrows][ncols]
template<typename T>
std::vector< std::vector<T> > vec_to_mat( const std::vector<T>& vec,
                                          _t ncols, 
                                          const _S ft="vec_to_mat" )
{
  const _S _fe("["+ft+"]: ERROR: Incompatible dimensions: ncols="+
               "["+std::to_string(ncols)+" x ]"+vec.size()+RET2 );
  
  if( ncols == 0 || vec.size() % ncols != 0 ) 
    throw std::domain_error( _fe.c_str() );
  
  const auto nrows = vec.size() / ncols;
  
  // declare an empty matrix. eventually, we build this up to an nrows x ncols matrix
  // the final matrix would be a collection of nrows rows
  // exch row in the matrix would a collection (a vector) containing ncols elements
  std::vector< std::vector<T> > mat;
  const auto begin = std::begin( vec );
  
  // add rows one by one to the matrix
  for( _t row = 0 ; row < nrows ; ++row ) // for each row [0,nrows-1] in the matrix
  {
    // add the row (a vector of ncols elements)
    // for example, if ncols = 12,
    // row 0 would contain elements in positions 0, 1, 2, ...10, 11 ie. [0,12)
    // row 1 would contain elements in positions 12, 13, ...23 ie. [12,24)
    // row 2 would contain elements in positions 24, 25, ...35 ie. [24,36)
    // in general, row r would contain elements at positions [ r*12, (r+1)*12 ]
    mat.push_back( { begin + row*ncols, begin + (row+1)*ncols } ) ;
    // the above as akin to:
    // construct the row containing elements at positions as described earlier
    // const std::vector<T> this_row( begin + row*ncols, begin + (row+1)*ncols ) ;
    // mat.push_back(this_row) ; // and add this row to the back of the vector
  }
  return( mat );
}

inline _S strVec_to_str( const std::vector<_S>& v, 
                         const _S sep = COM ) {
  if ( v.size() == 0 ) return( "" );
  
  std::ostringstream oss;
  for ( auto& x : v ) oss<<x<<sep;
  _S str = oss.str();
  str.resize( str.size() - sep.size() );
    
  return( str );
};

_S boolVec_to_str( const std::vector<_Y>& v, 
                   const _S del = "," )
{
  if ( v.size() == 0 ) return( "" );
  
  std::ostringstream oss;
  _I end = v.size()-1;
  for ( _I ii = 0; ii < end; ii++ ) oss << v[ii] << del;
  oss << v[end];
  
  return oss.str();
};

_S intVec_to_str( const std::vector<_i>& v,
                  const _S del = "," )
{
  if ( v.size() == 0 ) return( "" );
  
  std::ostringstream oss;
  _I end = v.size()-1;
  for ( _I ii = 0; ii < end; ii++ ) oss << v[ii] << del;
  oss << v[end];
  
  return oss.str();
};

_S uintVec_to_str( const std::vector<_I>& v,
                   const _S del = "," )
{
  if ( v.size() == 0 ) return( "" );
  
  std::ostringstream oss;
  _I end = v.size()-1;
  for ( _I ii = 0; ii < end; ii++ ) oss << v[ii] << del;
  oss << v[end];
  
  return oss.str();
};

_S ullVec_to_str( const std::vector<_Q>& v,
                  const _S del = "," )
{
  if ( v.size() == 0 ) return( "" );
  
  std::ostringstream oss;
  _I end = v.size()-1;
  for ( _I ii = 0; ii < end; ii++ ) oss << v[ii] << del;
  oss << v[end];
  
  return oss.str();
};

std::vector<_S> doubeVec_to_strVec( const std::vector<double>& vec )
{
  std::vector<_S> tempStr;
  
  for (_I i(0); i < vec.size(); ++i){
    std::ostringstream doubleStr;
    doubleStr << vec[i];    
    tempStr.push_back(doubleStr.str());
  }
  
  return tempStr;
}

_S dblVec_to_str( const std::vector<double>& v, 
                  const _S del ="," )
{
  if ( v.size() == 0 ) return( "" );
  
  std::ostringstream oss;
  std::vector<_S> vec = doubeVec_to_strVec(v);
  _I end = vec.size()-1;
  for ( _I ii = 0; ii < end; ii++ ) oss << vec[ii] << del;
  oss << vec[end];
  
  return oss.str();
}

template <typename T>
_S vec_to_str( const std::vector<T>& v,
               const _S del = ",",
               _I beg = 0,
               _I end = 0 )
{
  _I v_size = v.size();
  if ( v_size == 0 ) return( "" );
  
  if ( end == 0 ) end = v_size;
  if ( beg > v_size ) return( "" ) ;
  if ( end > v_size ) end = v_size;
  if ( beg >= end ) return( "" );
  
  std::ostringstream oss;
  for ( _I ii = beg; ii < end-1; ii++ ) oss << v[ii] << del;
  oss << v[end-1];
  
  return oss.str();
};

/*
 * Safe Seek Template for Binary Gzipped Files::
 * 
 */
template <typename F>
static _Y safe_seek_bgz( F& fh,
                         const _Q off_set = 0,
                         const _Y seek_cur = false,
                         
                         const _I vb = 0,
                         const _I vt = 4,
                         const _I tc = 0,
                         const _S func_tag = "safe_seek_bgz" )
{
  const _S tb(tc, '\t');
  const _S func_out("["+func_tag+"]: "+tb);
  const _S func_err("["+func_tag+"]: ERROR: ");
  const _S func_wrn("["+func_tag+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb >= vt;
  const _Y p1 = vb >= vt + 1;
  
  z_off_t gv;
  _Y success = true;
  
  if ( p0 ) std::cerr 
    << func_out << "Starting...\n"
    << func_out << "\t    off_set = " << off_set << "\n"
    << func_out << "\t   seek_cur = " << seek_cur
    << std::endl;
  
  if ( off_set > 0 ) {
    if ( p1 ) std::cerr 
      << func_out << "\t Attempting to seek to: '" <<off_set<<"'"<< std::endl;
    
    if ( seek_cur ) {
      try { 
        gv = gzseek( fh, off_set, SEEK_CUR );
      } catch (const std::exception& e) {
        std::cerr << func_err << "Caught ERROR: '"<<e.what()<<"'\n"<< std::endl;
        success = false;
      }
    } else {
      try { 
        gv = gzseek( fh, off_set, SEEK_SET );
      } catch (const std::exception& e) {
        std::cerr << func_err << "Caught ERROR: '"<<e.what()<<"'\n"<< std::endl;
        success = false;
      }
    }
    if ( p1 ) std::cerr 
      << func_out << "\t gzseek(gv): '" << gv << "'" << std::endl;
  }
  
  if ( p0 ) std::cerr 
    << func_out << "Done. Return: '" << success << "'\n" << std::endl;
  
  return( success );
}

/*
 * Safe Read Template for Binary Gzipped Files:: Scalar/Vector Values
 * 
 */
template <typename F, class T>
static _Y safe_read_bgz( F& fh,
                         T& buf,
                         const _I rep_cnt = 1,
                         const _S bin_str = "",
                         
                         const _I vb = 0,
                         const _I vt = 4,
                         const _I tc = 0,
                         const _S func_tag = "safe_read_bgz" )
{
  const _S tb(tc, '\t');
  const _S func_out("["+func_tag+"]: "+tb);
  const _S func_err("["+func_tag+"]: ERROR: ");
  const _S func_wrn("["+func_tag+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb >= vt;
  const _Y p1 = vb >= vt + 1;
  
  _i rv = 0;
  _Y success = true;
  
  if ( p0 ) std::cerr 
    << func_out << "Starting...\n" 
    << func_out << "\t    rep_cnt = " << rep_cnt << "\n"
    << func_out << "\t    bin_str = " << bin_str
    << std::endl;
  
  try {
    rv = gzread( fh, &buf, sizeof( buf ) * rep_cnt );
  } catch (const std::exception& e) {
    std::cerr << func_err << "Caught ERROR: '"<< e.what() << "'\n" << std::endl;
    success = false;
  }
  if ( p1 ) std::cerr 
    << func_out << "\t gzread(rv): '" << rv  << "'\n" << std::endl;
  
  if ( p0 ) std::cerr 
    << func_out << "Done. Return: '" << success << "'\n" << std::endl;
  
  return( success );
};

template <typename F, class T>
static _Y safe_load_bgz( F& fh,
                         T& buf,
                         _I rep_cnt = 1,
                         const _S bin_str = "",
                         
                         const _Q off_set = 0,
                         
                         const _Y seek_cur = false,
                         const _Y ilmn_str = false,
                         
                         const _I vb = 0,
                         const _I vt = 4,
                         const _I tc = 0,
                         const _S func_tag = "safe_load_bgz" )
{
  const _S tb(tc, '\t');
  const _S func_out("["+func_tag+"]: "+tb);
  const _S func_err("["+func_tag+"]: ERROR: ");
  const _S func_wrn("["+func_tag+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb >= vt;
  // const _Y p1 = vb >= vt + 1;
  
  _Y success = true;
  
  if ( p0 ) std::cerr 
    << func_out << "Starting...\n"
    << func_out << "\t    rep_cnt = " << rep_cnt << "\n"
    << func_out << "\t    bin_str = " << bin_str << "\n"
    << func_out << "\t    off_set = " << off_set << "\n"
    << func_out << "\t   seek_cur = " << seek_cur << "\n"
    << func_out << "\t   ilmn_str = " << ilmn_str
    << std::endl;
  
  /*
   * 1.) Seek to the position if off_set is non-zero
   * 
   */
  if ( off_set > 0 )
    success = safe_seek_bgz( fh, off_set, seek_cur, vb, vt+1, tc+1 );
  if ( !success ) return( success );
  
  /*
   * 3.) Finally Read and store the data
   * 
   */
  success = safe_read_bgz( fh, buf, rep_cnt, bin_str, vb, vt+1, tc+1 );
  if ( !success ) return( success );
  
  if ( p0 ) std::cerr 
    << func_out << "Done. Return: '" << success << "'\n" << std::endl;
  
  return( success );
};

_Y write_df_bgz( const _S& spec_str,
                 const _S& file,
                 
                 std::vector<_S>& name_vec,
                 std::vector<std::vector<_S>>& str_vecs,
                 std::vector<std::vector<_i>>& int_vecs,
                 std::vector<std::vector<_Y>>& bool_vecs,
                 
                 const _I precision = 1000000,
                 const _Y write_tsv = false,
                 
                 const _I vb = 0,
                 const _I vt = 1,
                 const _I tc = 0,
                 const _S ft="write_df_bgz" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0  = vb > vt + 0;
  const _Y p1  = vb > vt + 1;
  const _Y p10 = vb > vt + 10;
  
  _Y success = true;
  _Y testing = false;
  const _B idx_sep = '\t';
  const _B csv_sep = ',';
  
  _I nrow_cnt  = 0;
  if ( str_vecs.size() != 0 ) nrow_cnt = str_vecs[0].size();
  else if ( int_vecs.size() != 0 ) nrow_cnt = int_vecs[0].size();
  else if ( bool_vecs.size() != 0 ) nrow_cnt = bool_vecs[0].size();
  else {
    std::cerr <<_fe<< "All input vectors are empty! Exiting...\n"<<std::endl;
    return( false );
  }
  
  _I str_size  = str_vecs.size();
  _I int_size  = int_vecs.size();
  _I bool_size = bool_vecs.size();
  
  _I vecs_cnt  = spec_str.size();
  
  const _S bgz_path( file+".bgz" );
  const _S idx_path( file+".idx" );
  const _S tsv_path( file+".tsv" );
  
  const _Y bgz_clean = clean_file( bgz_path );
  const _Y idx_clean = clean_file( idx_path );
  const _Y tsv_clean = clean_file( tsv_path );
  
  std::vector<_Q> offs_vec;
  
  if ( p0 ) std::cerr
    <<_fo<< "Starting..." << "\n"
    <<_fo<< "\t    spec_str = " << spec_str << "\n" 
    <<_fo<< "\t   precision = " << precision << "\n\n"
    <<_fo<< "\t    str_size = " << str_size << "\n" 
    <<_fo<< "\t    int_size = " << int_size << "\n" 
    <<_fo<< "\t   bool_size = " << bool_size << "\n\n"
    <<_fo<< "\t    vecs_cnt = " << vecs_cnt << "\n" 
    <<_fo<< "\t    nrow_cnt = " << nrow_cnt << "\n\n" 
    <<_fo<< "\t   write_tsv = " << write_tsv << "\n" 
    <<_fo<< "\t    bgz_path = " << bgz_path <<"clean="<<bgz_clean<<"\n" 
    <<_fo<< "\t    idx_path = " << idx_path <<"clean="<<idx_clean<<"\n" 
    <<_fo<< "\t    tsv_path = " << tsv_path <<"clean="<<tsv_clean<<"\n"
    << std::endl;
  
  std::ofstream IDX_FH ( idx_path, std::ios::out );
  if ( !IDX_FH.is_open() ) {
    std::cerr <<_fe<< "Failed to open idx file: '" << idx_path << "'\n"
              << std::endl;
    return( false );
  }
  
  std::ofstream TSV_FH ( tsv_path, std::ios::out );
  if ( write_tsv ) {
    if ( !TSV_FH.is_open() ) {
      std::cerr <<_fe<< "Failed to open csv file: '" << tsv_path << "'\n"
                << std::endl;
      return( false );
    }
  }  
  /*
   * Basic Implementation First::
   *   - Write Binary for Address Vectors Only
   *   - Write all other data to Index file...
   */
  ostringstream os;
  gzFile gzo = gzopen( bgz_path.c_str(), "w");
  
  _I cur_len  = 0;
  _I str_idx  = 0;
  _I int_idx  = 0;
  _I bool_idx = 0;
  
  _c char_spec = '0';
  
  _S log_str = "";
  _S any_str = "";
  _S len_str = "";
  _S int_str = "";
  
  /*
   * 
   * Write Data File:: BGZ (Binary Gzipped)
   * 
   * 
   */
  std::vector<_I> len_vec( nrow_cnt );
  for ( _I ii = 0; ii < vecs_cnt; ii++ )
  {
    char_spec = spec_str[ii];
    
    if ( testing || write_tsv ) {
      log_str.clear();
      any_str.clear();
      len_str.clear();
      int_str.clear();
    }
    
    if ( p10 ) std::cerr
      <<_fo<< "Current Vector["<<ii<<"]: name='"<<name_vec[ii]<<"', char_spec='"
      <<char_spec<<"'"<<std::endl;
    
    if ( char_spec == 'c' ) {
      if ( testing || write_tsv ) log_str = strVec_to_str( str_vecs[str_idx] );
      if ( testing ) any_str = vec_to_str( str_vecs[str_idx] );
      
      len_vec.clear();
      // For strings we also need to write a vector of their lengths...
      for ( _I jj = 0; jj < nrow_cnt; jj++ ) {
        _I str_size = (_I)str_vecs[ii][jj].size();
        len_vec.push_back( str_size );
      }
      //  len_vec[jj] = (_I)str_vecs[ii][jj].size();
      
      os.write( reinterpret_cast<_c*> (&len_vec[0]), len_vec.size() * sizeof(len_vec[0]) );
      
      if ( testing || write_tsv ) len_str = vec_to_str( len_vec );
      
      if ( testing ) {
        int_str = uintVec_to_str( len_vec );
        
        for ( _I jj = 0; jj < nrow_cnt; jj++ ) std::cerr
          <<_fo<< "Seq Lens["<<ii<<","<<jj<<"]: name='"<<name_vec[ii]
          <<"', char_spec='"<<char_spec<<"', seq_len'"<<str_vecs[ii][jj].size()<<"'\n"
          <<_fo<<TAB<<" Seq_Lens: len_str = '"<<len_str<<"'\n"
          <<_fo<<TAB<<" Seq_Lens: int_str = '"<<int_str<<"'\n"
          <<std::endl;
      }
      
      cur_len = str_vecs[str_idx].size();
      os.write( reinterpret_cast<_c*> (&str_vecs[str_idx][0]), cur_len * sizeof(str_vecs[str_idx][0]) );
      str_idx++;
    } else if ( char_spec == 'i' || char_spec == 'd' || char_spec == 'n' ) {
      if ( testing || write_tsv ) log_str = intVec_to_str( int_vecs[int_idx] );
      if ( testing ) any_str = vec_to_str( int_vecs[int_idx] );
      
      cur_len = int_vecs[int_idx].size();
      // os.write( reinterpret_cast<_c*> (&int_vecs[int_idx][0]), cur_len * sizeof(int_vecs[int_idx][0]) );
      os.write( reinterpret_cast<_c*> (&int_vecs[int_idx][0]), cur_len * sizeof(_i) );
      int_idx++;
    } else if ( char_spec == 'l' ) {
      if ( testing || write_tsv ) log_str = boolVec_to_str( bool_vecs[bool_idx] );
      if ( testing ) any_str = vec_to_str( bool_vecs[bool_idx] );
      
      cur_len =  bool_vecs[bool_idx].size();
      os.write( reinterpret_cast<_c*> (&bool_vecs[bool_idx]), cur_len );
      bool_idx++;
    } else {
      std::cerr <<_fe<< "Unsupported character spec["<<ii<<"] ='"
                <<char_spec<<"'\n"<<std::endl;
      return( false );
    }
    
    if ( cur_len != nrow_cnt ) {
      std::cerr <<_fe<< "Number of rows do not match: "
                << "cur_len("<<cur_len<<") != "
                <<"nrow_cnt("<<nrow_cnt<<"); "
                <<"character spec["<<ii<<"] ='"<<char_spec<<"', "
                <<"name ='"<<name_vec[ii]<<"'\n"<<std::endl;
      // return( false );
    }
    // Write to buffer and then to gzip::
    offs_vec.push_back( gzoffset(gzo) );
    gzwrite( gzo, os.str().c_str(), os.str().size() );
    
    /*
     * Write Data File:: TSV (tab delimited)
     */
    if ( write_tsv ) TSV_FH 
      << name_vec[ii] << csv_sep 
      << char_spec << csv_sep
      << log_str << idx_sep
      << len_str
      << std::endl;
    
    if ( testing ) std::cerr
      <<_fo<< "Current Vector["<<ii<<"]: name='"<<name_vec[ii]<<"', char_spec='"<<char_spec<<"'\n"
      <<_fo<< "Current Vector["<<ii<<"]: log_str='"<<log_str<<"'\n"
      <<_fo<< "Current Vector["<<ii<<"]: any_str='"<<any_str<<"'\n"
      <<_fo<< "Current Vector["<<ii<<"]: len_str='"<<len_str<<"', size='"<<len_vec.size()<<"'\n"
      <<_fo<< "Current Vector["<<ii<<"]: int_str='"<<int_str<<"'\n"
      <<std::endl;
  }  
  gzclose( gzo );
  if ( write_tsv ) TSV_FH.close();
  if ( p1 ) std::cerr <<_fo<<"Wrote: BGZ='" << bgz_path << "'.\n" << std::endl;
  
  /*
   * 
   * Write Index File:: IDX
   * 
   * 
   */
  IDX_FH << "precision" << idx_sep 
         << "spec_str" << idx_sep
         << "nrow_cnt" << idx_sep
         << "name_str" << idx_sep
         << "off_sets"
         << std::endl;
  
  IDX_FH << precision << idx_sep
         << spec_str << idx_sep
         << nrow_cnt << idx_sep
         << vec_to_str(name_vec) << idx_sep
         << vec_to_str(offs_vec)
         << std::endl;
  // << strVec_to_str(name_vec) << idx_sep
  // << ullVec_to_str(offs_vec)
  
  IDX_FH.close();
  if ( p1 ) std::cerr <<_fo<<"Wrote: IDX='" << idx_path << "'.\n" << std::endl;
  
  return( success );
};

// std::vector<std::vector<double>>& dbl_vecs,

_Y load_df_bgz( const _S& file,
                
                std::vector<_S>& name_vec,
                std::vector<std::vector<_S>>& str_vecs,
                std::vector<std::vector<_i>>& int_vecs,
                std::vector<std::vector<_Y>>& bool_vecs,
                
                const _I vb = 0,
                const _I vt = 1,
                const _I tc = 0,
                const _S ft="load_df_bgz" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0  = vb > vt + 0;
  const _Y p1  = vb > vt + 1;
  const _Y p10 = vb > vt + 10;
  
  _Y success = true;
  // _Y testing = false;
  _Y testing = true;
  
  const _B idx1_sep = '\t';
  const _B idx2_sep = ',';
  
  const _S bgz_path( file+".bgz" );
  const _S idx_path( file+".idx" );
  
  name_vec.clear();
  str_vecs.clear();
  int_vecs.clear();
  bool_vecs.clear();
  
  if ( p0 ) std::cerr
    <<_fo<< "Starting..." << "\n"
    <<_fo<< "\t    bgz_path = " << bgz_path << "\n" 
    <<_fo<< "\t    idx_path = " << idx_path << "\n" 
    << std::endl;
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                             Check Files Exist::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  gzFile gzo = gzopen( bgz_path.c_str(), "r" );
  
  std::ifstream IDX_FH ( idx_path, std::ios::in );
  if ( !IDX_FH.is_open() ) {
    std::cerr <<_fe<<"Failed to open idx file: '"<<idx_path<<"'\n"<< std::endl;
    return( false );
  }
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                              Load Index Data::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  if ( p1 ) std::cerr <<_fo<<"Loading BGZ Index:: '"<<idx_path<<"'"<<std::endl;
  
  _S line;
  std::vector<_S> idx_head_vec;
  std::vector<_S> idx_data_vec;
  
  getline( IDX_FH, line );
  parse_line( line, idx1_sep, idx_head_vec );
  getline( IDX_FH, line );
  parse_line( line, idx1_sep, idx_data_vec );
  IDX_FH.close();
  
  _I precision = 0;
  _S spec_str = "";
  _I nrow_cnt = 0;
  _I ncol_cnt = 0;
  std::vector<_S> strs_vec;
  std::vector<_Q> offs_vec;
  
  for ( _I ii = 0; ii < idx_head_vec.size(); ii++ ) {
    if ( p1 ) std::cerr
      <<_fo<<"Parsing parameter["<<ii<<"]: '"<<idx_head_vec[ii]<<"', "
      << " string value = '"<<idx_data_vec[ii]<<"'"<<std::endl;
    
    if ( idx_head_vec[ii].compare("precision") == 0 ) {
      precision = std::stoul( idx_data_vec[ii] );
    } else if ( idx_head_vec[ii].compare("spec_str") == 0 ) {
      spec_str = idx_data_vec[ii];
    } else if ( idx_head_vec[ii].compare("nrow_cnt") == 0 ) {
      nrow_cnt = std::stoul( idx_data_vec[ii] );
    } else if ( idx_head_vec[ii].compare("name_str") == 0 ) {
      parse_line( idx_data_vec[ii], idx2_sep, name_vec );
    } else if ( idx_head_vec[ii].compare("off_sets") == 0 ) {
      parse_line( idx_data_vec[ii], idx2_sep, strs_vec );
      offs_vec.resize( strs_vec.size() );
      for ( _I jj = 0; jj < strs_vec.size(); jj++ )
        offs_vec[jj] = std::stoull( strs_vec[jj] );
    } else {
      if ( p1 ) std::cerr
        <<_fw<<"Unknown index parameter: '"<<idx_head_vec[ii]<<"' Skipping...\n"
        << std::endl;
    }
  }
  ncol_cnt = spec_str.size();
  
  if ( ncol_cnt == 0 ) {
    if ( p1 ) std::cerr
      <<_fo<<"Failed to load index correctly! ncol_cnt ='"<<ncol_cnt<<"'\n"
      << std::endl;
    return( false );
  }
  
  if ( p1 ) std::cerr
    <<_fo<< "Index Data::" << "\n"
    <<_fo<< "\t    precision = '" << precision << "'\n" 
    <<_fo<< "\t    spec_str  = '" << spec_str << "'\n" 
    <<_fo<< "\t    nrow_cnt  = '" << nrow_cnt << "'\n" 
    <<_fo<< "\t    ncol_cnt  = '" << ncol_cnt << "'\n" 
    <<_fo<< "\t    name_vec  = '" << vec_to_str(name_vec) << "'\n" 
    <<_fo<< "\t    strs_vec  = '" << vec_to_str(strs_vec) << "'\n" 
    << std::endl;
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   * 
   *                              Load Stack Data::
   * 
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  if ( p1 ) std::cerr <<_fo<<"Loading BGZ Data:: '"<<bgz_path<<"'"<<std::endl;
  
  _I str_cnt  = 0;
  _I int_cnt  = 0;
  _I dbl_cnt  = 0;
  _I bool_cnt = 0;
  for ( _I ii = 0; ii < ncol_cnt; ii++ )
  {
    _c char_spec = spec_str[ii];
    
    switch (char_spec)  {
    case 'c':
      str_cnt++;
      break;
      
    case 'i':
      int_cnt++;
      break;
      
    case 'd':
      dbl_cnt++;
      int_cnt++;
      break;
      
    case 'n':
      dbl_cnt++;
      int_cnt++;
      break;
      
    case 'l':
      bool_cnt++;
      break;
      
    default:
      std::cerr <<_fe<< "Unsupported character spec["<<ii<<"] ='"
                <<char_spec<<"'\n"<<std::endl;
      return( false );
    }
  }
  
  if ( p1 ) std::cerr
    <<_fo<< "Spec Type Counts::" << "\n"
    <<_fo<< "\t    str_cnt  = '" << str_cnt << "'\n" 
    <<_fo<< "\t    int_cnt  = '" << int_cnt << "'\n" 
    <<_fo<< "\t    dbl_cnt  = '" << dbl_cnt << "'\n" 
    <<_fo<< "\t    bool_cnt = '" << bool_cnt << "'\n" 
    << std::endl;
  
  // Resize Storage Vectors::
  str_vecs.resize( str_cnt );
  int_vecs.resize( int_cnt );
  bool_vecs.resize( bool_cnt );
  
  _i rv = 0;
  _I cur_len  = 0;
  _I str_idx  = 0;
  _I int_idx  = 0;
  _I bool_idx = 0;
  std::vector<_I> len_vecT( nrow_cnt );
  // std::vector<_i> int_vecT( nrow_cnt );
  
  for ( _I ii = 0; ii < ncol_cnt; ii++ )
  {
    _S log_str = "";
    _S any_str = "";
    _c char_spec = spec_str[ii];
    
    if ( p10 ) std::cerr
      <<_fo<< "Current Vector["<<ii<<"]: name='"<<name_vec[ii]<<"', char_spec='"
      <<char_spec<<"'"<<std::endl;
    
    if ( char_spec == 'c' ) {
      
      str_idx++;
    } else if ( char_spec == 'i' || char_spec == 'd' || char_spec == 'n' ) {
      
      int_vecs[int_idx].clear();
      int_vecs[int_idx].resize( nrow_cnt );
      
      rv = gzread( gzo, &int_vecs[int_idx][0], nrow_cnt*sizeof(_i) );
      std::cerr<<_fo<<"\t val[0]='"<<int_vecs[int_idx][0]<<"'"<<std::endl;
      std::cerr<<_fo<<"\t val[1]='"<<int_vecs[int_idx][1]<<"'"<<std::endl;
      std::cerr<<_fo<<"\t val[2]='"<<int_vecs[int_idx][2]<<"'"<<std::endl;
      std::cerr<<std::endl;
      
      cur_len = int_vecs.size();
      
      //
      // Temp Data Structure for Testing, real code above...
      //
      
      // std::vector<int> tmp_vec;
      // tmp_vec.clear();
      // tmp_vec.resize( nrow_cnt );
      // 
      // rv = gzread( gzo, &tmp_vec[0], nrow_cnt*sizeof(int) );
      // std::cerr <<_fo<< "val[0]='"<< tmp_vec[0] <<"' int_idx='"<<int_idx<<"'"<<std::endl;
      // std::cerr <<_fo<< "val[1]='"<< tmp_vec[1] <<"' int_idx='"<<int_idx<<"'"<<std::endl;
      // std::cerr <<_fo<< "val[2]='"<< tmp_vec[2] <<"' int_idx='"<<int_idx<<"'\n"<<std::endl;
      // 
      // cur_len = tmp_vec.size();
      
      int_idx++;
    } else if ( char_spec == 'l' ) {
      
      bool_idx++;
    } else {
      if ( p10 ) std::cerr <<_fe<< "Unsupported character spec["<<ii<<"] ='"
                           <<char_spec<<"'\n"<<std::endl;
      return( false );
    }
    
    if ( cur_len != nrow_cnt ) {
      std::cerr <<_fe<< "Number of rows do not match: "<<cur_len<<" != "
                <<nrow_cnt<<"; character spec["<<ii<<"] ='"<<char_spec<<"', "
                <<"name ='"<<name_vec[ii]<<"'\n"<<std::endl;
      // return( false );
    }
    
    //
    // NOTE:: May need to use dblVec_to_str() method now...
    //
    if ( testing ) std::cerr
      <<_fo<< "rv="<<rv<<"\n"
      <<_fo<< "Current Vector["<<ii<<"]: str_idx='"<<str_idx<<"'\n"
      <<_fo<< "Current Vector["<<ii<<"]: log_str='"<<log_str<<"'\n"
      <<_fo<< "Current Vector["<<ii<<"]: any_str='"<<any_str<<"'\n"
      <<_fo<< "Current Vector["<<ii<<"]: boolIdx='"<<bool_idx<<"'\n"
      <<std::endl;
  }
  gzclose( gzo );
  
  if ( testing ) {
    for ( _I ii = 0; ii < int_vecs.size(); ii++ ) {
      for ( _I jj = 0; jj < int_vecs[ii].size(); jj++ ) {
        std::cerr <<_fo<< "\t int_vecs["<<ii<<"]["<<jj<<"]='"<< int_vecs[ii][jj] <<"'" <<std::endl;
      }
    }
  }
  
  return( success );
}

// This function is present so I don't see the stupid warning message of unused
//   functions during compile pile time...
//
void unused_mark_test()
{
  _c c = 'c';
  _I idx = 0;
  std::vector<_I> iv;
  std::vector<_S> sv;
  std::vector<std::vector<_S>> sv2;
  _S s = "";
  _S f = "";
  std::map< _S, std::map< _S, std::map< _S, std::vector<_I> > > > xxx;
  
  _S x1 = str_to_str_idx( s );
  parse_into_postion( s, c, sv2 );
  std::vector<_S> lc = load_assoc_stack( f, 0, 0 );
  _i tb = is_top_bot( s );
  std::vector<_S> ex = expand_iupac_seq( s );
  
  _i counter = 0;
  _S file;
  std::vector<_S> stack;
  _Y w_bool = write_assoc_stack( stack, 0, 0, 0, 0, file );
  if ( w_bool ) counter++;
  std::cerr << "counter='"<<counter<<"'\n" <<std::endl;
  
  std::vector<std::vector<_S>> a;
  _I cint = 0; 
  std::vector<_S> nv;
  replace_nuc_at(s, idx, s );
  // tb = false;
  _Y success = improbe_seq( sv2,
                            s, s, idx, 
                            sv, xxx,
                            s, s, s, s,
                            50, 0, 0, 0 );
  std::cerr << "tb='"<<tb<<"'\n" <<std::endl;
  std::cerr << "success='"<<success<<"'\n" <<std::endl;
  std::cerr << "   cint='"<<cint<<"'\n" <<std::endl;
  
  
  success = legal_file( file );
  // success = get_file_name( file );
  success = build_path( file );
  std::cerr << "success='"<<success<<"'\n" <<std::endl;
  
  cint++;
};

};

#endif /* bisulfite.h */

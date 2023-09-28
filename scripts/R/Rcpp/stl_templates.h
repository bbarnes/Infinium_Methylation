#ifndef __STL_TEMPLATES_H__
#define __STL_TEMPLATES_H__

// Core STL Classes::
#include <regex>
#include <string>
#include <iterator>
#include <vector>
#include <map>
#include <cstring>   /* std::memset */
#include <ranges>
#include <utility>
#include <tuple>

// File IO::
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept> /* std::invalid_argument */

// System Stuff:: i.e. piping, etc.
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>

// File IO:: Binary/Gzipped::
#include <zlib.h>

// Math Stuff::
#include <math.h>
#include <cmath>

// Template Type Detection::
#include <experimental/type_traits>
#include <cstddef>

// Standard Stuff::
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <ctype.h>
#include <cctype>
#include <stdint.h>
#include <inttypes.h>
#include <type_traits>
#include <iomanip>
#include <algorithm>
#include <libgen.h>
#include <unistd.h>

// Boost Libaries::
#include <boost/regex.hpp>

#include <boost/math/distributions/empirical_cumulative_distribution_function.hpp>
using boost::math::empirical_cumulative_distribution_function;

#include <boost/math/statistics/t_test.hpp>
using boost::math::statistics::two_sample_t_test;

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/max.hpp>
using namespace boost::accumulators;

#include <iostream>
#include <boost/tuple/tuple.hpp>

using namespace std;

namespace stl_templates {

#define GZ_LENGTH 0x1000
#define BASE_MASK 0x3 /* binary: 11 */
#define _SC static const

// Coding-Lazyness + Binary-Code Reminder Definitions::
typedef char               _c;
typedef signed char        _b;
typedef unsigned char      _B;
typedef short              _h;
typedef unsigned short     _H;
typedef int                _i;
typedef unsigned int       _I;
typedef float              _f;
typedef double             _d;
typedef long               _l;
typedef unsigned long      _L;
typedef long long          _q;
typedef unsigned long long _Q;
typedef std::string        _S;
typedef bool               _Y;
typedef std::size_t        _t;

// Common-Code-Lazyness
_SC _B _RET = '\n';
_SC _B _TAB = '\t';
_SC _B _COM = ',';

_SC _B _RET2 = '\n';
_SC _B _TAB2 = '\t';
_SC _B _COM2 = ',';

_SC _S RET = "\n";
_SC _S TAB = "\t";
_SC _S COM = ",";

_SC _S RET2 = "\n";
_SC _S TAB2 = "\t";
_SC _S COM2 = ",";

_SC _S BRK = "# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #\n";

#ifdef _WIN32
_SC char SYS_SEP = '\\';
#else
_SC char SYS_SEP = '/';
#endif

// Template Type Detection::
typedef _i not_vec;
typedef _i not_lst;

// Vector Growth Rate, but probably not needed...
_SC _I STACK_GROW_RATE = 2;

/*
 * Boost Tuplet Class::
 * 
 * boost::tuple<int, double, std::string> t(1, 3.14, "PI");
 * 
 * template<typename T1, typename T2, typename T3>
 class Triplet
 {
public:
 // boost::tuple<T1, T2, T3> dat;
 // std::tuple<int, int, int> three;
 T1 first;
 T2 middle;
 T3 last;
 
 Triplet(const T1 &x, const T2 &y, const T3 &z)
 {
 std::tuple<T1, T2, T3> t( x,y,z );
 this->first  = std::get<0>( t );
 this->middle = std::get<1>( t );
 this->last   = std::get<2>( t );
 }
 };
 * 
 */

/*
 * #include <iostream>
 * #include <utility>
 * #include <tuple>
 * 
 * template<class T>
 * void test(T t)
 * {
 * int a[std::tuple_size<T>::value]; // can be used at compile time
 * std::cout << std::tuple_size<T>::value << '\n'; // or at run time
 * }
 *  
 * int main()
 * {
 * test(std::make_tuple(1, 2, 3.14));
 * test(std::make_pair(1, 3.14));
 *  }
 * Reports 3 and 2 respectively...
 * 
 */


/*
 * Testing for Template("data structures"), i.e. NOT type(char, _i, double...)
 * 
 */

// Duplicate Call to is_specialization
template<typename Test_Vec, template<typename...> class Ref_Vec> struct is_specialization_vec : std::false_type {};

template<template<typename...> class Ref_Vec, typename... Args> struct is_specialization_vec<Ref_Vec<Args...>, Ref_Vec>: std::true_type {};

// TBD:: Attempt to include STL Pair Class::
template<typename Test_Pair, template<typename...> class Ref_Pair> struct is_specialization_pair : std::false_type {};

template<template<typename...> class Ref_Pair, typename... Args>
struct is_specialization_pair<Ref_Pair<Args...>, Ref_Pair>: std::true_type {};

/*
 * stringify detection definitions::
 * 
 */
// 1- detecting if std::to_string is valid on T
template<typename T>
using std_to_string_expression = decltype(std::to_string(std::declval<T>()));

template<typename T>
constexpr _Y has_std_to_string = std::experimental::is_detected<std_to_string_expression, T>::value;

// 2- detecting if to_string is valid on T
template<typename T>
using to_string_expression = decltype(to_string(std::declval<T>()));

template<typename T>
constexpr _Y has_to_string = std::experimental::is_detected<to_string_expression, T>::value;

// 3- detecting if T can be sent to an ostringstream
template<typename T>
using ostringstream_expression = decltype(std::declval<std::ostringstream&>() << std::declval<T>());

template<typename T>
constexpr _Y has_ostringstream = std::experimental::is_detected<ostringstream_expression, T>::value;

// 4 - detecting if T is a class:: Pair
template<typename T>
using is_template_pair = decltype( is_specialization_pair<T, std::pair>::value );
// using is_template_pair = decltype( is_specialization_vec<T, std::pair>::value );
// using is_template_pair = decltype( is_specialization<T, std::pair>::value );

template<typename T>
constexpr _Y has_tempate_pair = std::experimental::is_detected<is_template_pair, T>::value;

// 5 - detecting if T is a class:: Vector
template<typename T>
using is_template_vector = decltype( is_specialization_vec<T, std::vector>::value );

template<typename T>
constexpr _Y has_tempate_vector = std::experimental::is_detected<is_template_vector, T>::value;

// 6 - detecting if T is a class:: Map
template<typename T>
using is_template_map = decltype( is_specialization_vec<T, std::map>::value );

template<typename T>
constexpr _Y has_tempate_map = std::experimental::is_detected<is_template_map, T>::value;

// template<typename T>
// constexpr _Y is_tempate_class = std::is_class<T>::value;

/*
 * stringify template definitions:: All Elements
 * 
 */
// 1::  std::to_string is valid on T
template<typename T, typename 
  std::enable_if<has_std_to_string<T>, _i>::type = 0>
_S stringify_all(T const& t)
{
  return std::to_string(t);
};

// 2::  std::to_string is not valid on T, but to_string is
template<typename T, typename 
  std::enable_if<
    !has_std_to_string<T> && 
    has_to_string<T>, _i>::type = 0>
    _S stringify_all(T const& t)
    {
      return to_string(t);
    };

// 3::  Neither std::string nor to_string work on T, let's stream it then
template<typename T, typename 
  std::enable_if<
    !has_std_to_string<T> && 
    !has_to_string<T> && 
    has_ostringstream<T>, _i>::type = 0>
    _S stringify_all(T const& t)
    {
      std::ostringstream oss;
      oss << t;
      return oss.str();
    };

// 4:: Neither std::string nor to_string work on T, nor stream it then must be
//     standard template libaray
template<typename T, typename
  std::enable_if<
    !has_std_to_string<T> &&
    !has_to_string<T> &&
    !has_ostringstream<T>, _i>::type = 0 &&
    // is_specialization<T, std::pair>::value &&
    true >
    _S stringify_all(T const& t)
    {
      std::ostringstream oss;
      // Cancled:: Old attampt to use std::pair::
      // Cancled:: oss << stringify_all( t.first ) << COM << stringify_all( t.second );
      for(auto it = std::begin(t); it != std::end(t); ++it) {
        if ( it != std::begin(t) ) oss << COM;
        oss << stringify_all( *it );
      }
      
      return oss.str();
    };

/*
 * stringify template definitions:: First Element
 * 
 */
// 1::  std::to_string is valid on T
template<typename T, typename 
  std::enable_if<has_std_to_string<T>, _i>::type = 0>
_S stringify_beg(T const& t)
{
  return std::to_string(t);
};

// 2::  std::to_string is not valid on T, but to_string is
template<typename T, typename 
  std::enable_if<
    !has_std_to_string<T> && 
    has_to_string<T>, _i>::type = 0>
    _S stringify_beg(T const& t)
    {
      return to_string(t);
    };

// 3::  Neither std::string nor to_string work on T, let's stream it then
template<typename T, typename 
  std::enable_if<
    !has_std_to_string<T> && 
    !has_to_string<T> && 
    has_ostringstream<T>, _i>::type = 0>
    _S stringify_beg(T const& t)
    {
      std::ostringstream oss;
      oss << t;
      return oss.str();
    };

// 4:: Neither std::string nor to_string work on T, nor stream it then must be
//     standard template libaray
template<typename T, typename
  std::enable_if<
    !has_std_to_string<T> &&
    !has_to_string<T> &&
    !has_ostringstream<T>, _i>::type = 0 &&
    true > _S stringify_beg(T const& t)
    {
      std::ostringstream oss;
      // for (auto& it : M) {
      //   A.push_back(it);
      // }
      
      // CURRENT POINT::
      for ( auto& it : t ) oss << stringify(it);
      // Two Previous Attempts Below
      //   oss << t.front();
      //   oss << t[0];
      return oss.str();
    };

/*
 * stringify template definitions:: Last Element
 * 
 */
// 1::  std::to_string is valid on T
template<typename T, typename 
  std::enable_if<has_std_to_string<T>, _i>::type = 0>
_S stringify_end(T const& t)
{
  return std::to_string(t);
};

// 2::  std::to_string is not valid on T, but to_string is
template<typename T, typename 
  std::enable_if<
    !has_std_to_string<T> && 
    has_to_string<T>, _i>::type = 0>
    _S stringify_end(T const& t)
    {
      return to_string(t);
    };

// 3::  Neither std::string nor to_string work on T, let's stream it then
template<typename T, typename 
  std::enable_if<
    !has_std_to_string<T> && 
    !has_to_string<T> && 
    has_ostringstream<T>, _i>::type = 0>
    _S stringify_end(T const& t)
    {
      std::ostringstream oss;
      oss << t;
      return oss.str();
    };

// 4:: Neither std::string nor to_string work on T, nor stream it then must be
//     standard template libaray
template<typename T, typename
  std::enable_if<
    !has_std_to_string<T> &&
    !has_to_string<T> &&
    !has_ostringstream<T>, _i>::type = 0 &&
    true >
    _S stringify_end(T const& t)
    {
      std::ostringstream oss;
      oss << t.back();
      return oss.str();
    };

/*
 * stringify template definitions:: Element Summary
 * 
 */
// 1::  std::to_string is valid on T
template<typename T, typename 
  std::enable_if<has_std_to_string<T>, _i>::type = 0>_S stringify_sum(T const& t)
  {
    return std::to_string(t);
  };

// 2::  std::to_string is not valid on T, but to_string is
template<typename T, typename 
  std::enable_if<
    !has_std_to_string<T> && 
    has_to_string<T>, _i>::type = 0>
    _S stringify_sum(T const& t)
    {
      return to_string(t);
    };

// 3::  Neither std::string nor to_string work on T, let's stream it then
template<typename T, typename 
  std::enable_if<
    !has_std_to_string<T> && 
    !has_to_string<T> && 
    has_ostringstream<T>, _i>::type = 0>
    _S stringify_sum(T const& t)
    {
      std::ostringstream oss;
      oss << t;
      return oss.str();
    };

// 4:: Neither std::string nor to_string work on T, nor stream it then must be
//     standard template libaray
template<typename T, typename
  std::enable_if<
    !has_std_to_string<T> &&
    !has_to_string<T> &&
    !has_ostringstream<T>, _i>::type = 0 &&
    true >
    _S stringify_sum(T const& t)
    {
      const auto [min, max] = std::minmax_element(begin(t), end(t));
      
      std::ostringstream oss;
      // oss << "Index Range: ["<<t.begin()<<":"<<t.end()<<"]; "
      // oss << "Index Range:" ~ 0.0,
      
      oss << "Index Range: ["<<t.front()<<":"<<t.back()<<"]; ";
      // oss << "Index Range: ["<<t.begin()<<":"<<t.end()<<"]; "
      oss << "Index_Range:"<< "Value Range: ["<<*min<<":"<<*max<<"]";
      return oss.str();
    };

// 4.5:: Neither std::string nor to_string work on T, nor stream it then must be
//     standard template libaray
// template<typename T, typename
//   std::enable_if<
//     !has_std_to_string<T> &&
//     !has_to_string<T> &&
//     !has_ostringstream<T>, _i>::type = 0 &&
//     true >
//     _S stringify_sum(T const& t)
//     {
//       const auto [min, max] = std::minmax_element(begin(t), end(t));
//       
//       std::ostringstream oss;
//       oss << "Index Range: ["<<t.begin()<<":"<<t.end()<<"]; "
//       // oss << "Index Range: ["<<t.front()<<":"<<t.back()<<"]; "
//       oss << "Index Range:" ~ 0.0,
//       // oss << "Index_Range:"<< "Value Range: ["<<*min<<":"<<*max<<"]";
//       return oss.str();
//     };

// 5:: Neither std::string nor to_string work on T, nor stream it then must be
//     standard template libaray
// template<typename T, typename
//   std::enable_if<
//     !has_std_to_string<T> &&
//     !has_to_string<T> &&
//     !has_ostringstream<T>, std::vector<_i> >::type = 0 &&
//     true >
//     _S stringify_sum(T const& t)
//     {
//       const auto [min, max] = std::minmax_element(begin(t), end(t));
//       
//       std::ostringstream oss;
//       oss<<" Index Range: ["<<t.front()<<":"<<t.back()<<"]; "
//          <<"Value Range: ["<<*min<<":"<<*max<<"]";
//       return oss.str();
//     };

/*
 * System Independet Handling in Binary (1 byte = 8 bits): Reading
 *   https://www.cplusplus.com/articles/DzywvCM9/
 * 
 */
uint8_t ReadU8( istream& fh)
{
  uint8_t val;
  uint8_t bytes[1];
  
  fh.read( (char*)bytes, 1 );  // read 2 bytes from the file handle
  val = bytes[0];  // construct the 16-bit value from those bytes
  
  // Simplify... ??? Maybe not...
  // fh.read( (char*)val, 1 );
  
  return val;
};

uint16_t ReadU16( istream& fh)
{
  uint16_t val;
  uint8_t bytes[2];
  
  fh.read( (char*)bytes, 2 );  // read 2 bytes from the file handle
  val = bytes[0] | (bytes[1] << 8);  // construct the 16-bit value from those bytes
  
  return val;
};

uint32_t ReadU32( istream& fh)
{
  uint32_t val;
  uint8_t bytes[4];
  
  fh.read( (char*)bytes, 4 );  // read 2 bytes from the file handle
  val = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);  // construct the 32-bit value from those bytes
  
  return val;
};

_S ReadString(istream& fh)
{
  uint32_t len = ReadU32(fh);
  
  char* buffer = new char[len];
  fh.read(buffer, len);
  
  _S str( buffer, len );
  delete[] buffer;
  
  return str;
};

uint16_t ReadU16_BGZ( gzFile& fh)
{
  uint16_t val;
  uint8_t bytes[2];
  
  gzread( fh, (char*)bytes, 2 );  // read 2 bytes from the gzFile handle
  val = bytes[0] | (bytes[1] << 8);  // construct the 16-bit value from those bytes
  
  return val;
};

/*
 * System Independet Handling in Binary (1 byte = 8 bits): Writing
 *   https://www.cplusplus.com/articles/DzywvCM9/
 * 
 */
void WriteU8( ostream& fh, uint8_t val)
{
  uint8_t bytes[1];
  
  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;  // low byte
  // bytes[1] = (val >> 8) & 0xFF;  // high byte
  
  // write those bytes to the file handle
  fh.write( (char*)bytes, 1 );
};

void WriteU16( ostream& fh, uint16_t val)
{
  uint8_t bytes[2];
  
  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;  // low byte
  bytes[1] = (val >> 8) & 0xFF;  // high byte
  
  // write those bytes to the file handle
  fh.write( (char*)bytes, 2 );
};

void WriteU32( ostream& fh, uint32_t val)
{
  uint8_t bytes[4];
  
  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;  // low byte
  bytes[1] = (val >> 8) & 0xFF;  // high byte
  bytes[2] = (val >> 16) & 0xFF;  // high byte
  bytes[3] = (val >> 24) & 0xFF;  // high byte
  
  // write those bytes to the file handle
  fh.write( (char*)bytes, 2 );
};

void WriteString( ostream& fh, _S str)
{
  uint32_t len = str.length();
  
  WriteU32(fh, len);
  fh.write( str.c_str(), len );
};

/*
 * Odds are WriteUN_BGZ is going to be horrible at compression...
 * 
 */
void WriteU8_BGZ( gzFile& fh, uint8_t val)
{
  uint8_t bytes[1];
  
  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;  // low byte
  // bytes[1] = (val >> 8) & 0xFF;  // high byte
  
  ostringstream os;
  // os.write( reinterpret_cast<char*> (&this->Address_Idx[0]), this->Address_Cnt * sizeof(this->Address_Idx[0]) );
  // os.write( reinterpret_cast<char*> (&bytes[0]), this->Address_Cnt * sizeof(this->Address_Idx[0]) );
  os.write( reinterpret_cast<char*> (bytes), 1 );
  
  // write those bytes to the gzFile handle
  // gzwrite( fh, (char*)bytes, 2 );
  gzwrite( fh, os.str().c_str(), os.str().size() );
};

void WriteU16_BGZ( gzFile& fh, uint16_t val)
{
  uint8_t bytes[2];
  
  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;  // low byte
  bytes[1] = (val >> 8) & 0xFF;  // high byte
  
  ostringstream os;
  // os.write( reinterpret_cast<char*> (&this->Address_Idx[0]), this->Address_Cnt * sizeof(this->Address_Idx[0]) );
  // os.write( reinterpret_cast<char*> (&bytes[0]), this->Address_Cnt * sizeof(this->Address_Idx[0]) );
  os.write( reinterpret_cast<char*> (bytes), 2 );
  
  // write those bytes to the gzFile handle
  // gzwrite( fh, (char*)bytes, 2 );
  gzwrite( fh, os.str().c_str(), os.str().size() );
};



struct assoc_stack_type {
  _B dtype;
  _B struc;
  _I depth;
};

/*
 * Determine Data Type::
 * 
 */
template<class T>
assoc_stack_type get_assoc_stack_type(T test){ 
  _B dtype = 'x';
  _B struc = 'x';
  _I depth =  0;
  
  if constexpr (std::is_same<T, _c>::value) {
    dtype = 'c';
    struc = 't';
    depth =  0;
    // return( "c0" );
  } else if constexpr (std::is_same<T, _S>::value) {
    dtype = 'S';
    struc = 't';
    depth =  0;
    // return( "S0" );
  } else if constexpr (std::is_same<T, _b>::value) {
    dtype = 'b';
    struc = 't';
    depth =  0;
    // return( "b0" );
  } else if constexpr (std::is_same<T, _B>::value) {
    dtype = 'B';
    struc = 't';
    depth =  0;
    // return( "B0" );
  } else if constexpr (std::is_same<T, _h>::value) {
    dtype = 'h';
    struc = 't';
    depth =  0;
    // return( "h0" );
  } else if constexpr (std::is_same<T, _H>::value) {
    dtype = 'H';
    struc = 't';
    depth =  0;
    // return( "H0" );
  } else if constexpr (std::is_same<T, _i>::value) {
    dtype = 'i';
    struc = 't';
    depth =  0;
    // return( "i0" );
  } else if constexpr (std::is_same<T, _I>::value) {
    dtype = 'I';
    struc = 't';
    depth =  0;
    // return( "I0" );
  } else if constexpr (std::is_same<T, _l>::value) {
    dtype = 'l';
    struc = 't';
    depth =  0;
    // return( "l0" );
  } else if constexpr (std::is_same<T, _L>::value) {
    dtype = 'L';
    struc = 't';
    depth =  0;
    // return( "L0" );
  } else if constexpr (std::is_same<T, _q>::value) {
    dtype = 'q';
    struc = 't';
    depth =  0;
    // return( "q0" );
  } else if constexpr (std::is_same<T, _Q>::value) {
    dtype = 'Q';
    struc = 't';
    depth =  0;
    // return( "Q0" );
    
    /*
     * Vector Class Matching::
     * 
     */
  } else if constexpr (std::is_same<T, std::vector<_c>>::value) {
    dtype = 'c';
    struc = 'v';
    depth =  1;
    // return( "cV1" );
  } else if constexpr (std::is_same<T, std::vector<_S>>::value) {
    dtype = 'S';
    struc = 'v';
    depth =  1;
    // return( "SV1" );
  } else if constexpr (std::is_same<T, std::vector<_b>>::value) {
    dtype = 'b';
    struc = 'v';
    depth =  1;
    // return( "bV1" );
  } else if constexpr (std::is_same<T, std::vector<_B>>::value) {
    dtype = 'B';
    struc = 'v';
    depth =  1;
    // return( "BV1" );
  } else if constexpr (std::is_same<T, std::vector<_h>>::value) {
    dtype = 'h';
    struc = 'v';
    depth =  1;
    // return( "hV1" );
  } else if constexpr (std::is_same<T, std::vector<_H>>::value) {
    dtype = 'H';
    struc = 'v';
    depth =  1;
    // return( "HV1" );
  } else if constexpr (std::is_same<T, std::vector<_i>>::value) {
    dtype = 'i';
    struc = 'v';
    depth =  1;
    // return( "iV1" );
  } else if constexpr (std::is_same<T, std::vector<_I>>::value) {
    dtype = 'I';
    struc = 'v';
    depth =  1;
    // return( "IV1" );
  } else if constexpr (std::is_same<T, std::vector<_l>>::value) {
    dtype = 'l';
    struc = 'v';
    depth =  1;
    // return( "lV1" );
  } else if constexpr (std::is_same<T, std::vector<_L>>::value) {
    dtype = 'L';
    struc = 'v';
    depth =  1;
    // return( "LV1" );
  } else if constexpr (std::is_same<T, std::vector<_q>>::value) {
    dtype = 'q';
    struc = 'v';
    depth =  1;
    // return( "qV1" );
  } else if constexpr (std::is_same<T, std::vector<_Q>>::value) {
    dtype = 'Q';
    struc = 'v';
    depth =  1;
    // return( "QV1" );
    
    /*
     * Map Class Match::
     * 
     */
  } else if constexpr (std::is_same<T, std::map<_c,_I>>::value) {
    dtype = 'c';
    struc = 'm';
    depth =  1;
    // return( "cM1" );
  } else if constexpr (std::is_same<T, std::map<_S,_I>>::value) {
    dtype = 'S';
    struc = 'm';
    depth =  1;
    // return( "SM1" );
  } else if constexpr (std::is_same<T, std::map<_b,_I>>::value) {
    dtype = 'b';
    struc = 'm';
    depth =  1;
    // return( "bM1" );
  } else if constexpr (std::is_same<T, std::map<_B,_I>>::value) {
    dtype = 'B';
    struc = 'm';
    depth =  1;
    // return( "BM1" );
  } else if constexpr (std::is_same<T, std::map<_h,_I>>::value) {
    dtype = 'h';
    struc = 'm';
    depth =  1;
    // return( "hM1" );
  } else if constexpr (std::is_same<T, std::map<_H,_I>>::value) {
    dtype = 'H';
    struc = 'm';
    depth =  1;
    // return( "HM1" );
  } else if constexpr (std::is_same<T, std::map<_i,_I>>::value) {
    dtype = 'i';
    struc = 'm';
    depth =  1;
    // return( "iM1" );
  } else if constexpr (std::is_same<T, std::map<_I,_I>>::value) {
    dtype = 'I';
    struc = 'm';
    depth =  1;
    // return( "IM1" );
  } else if constexpr (std::is_same<T, std::map<_l,_I>>::value) {
    dtype = 'l';
    struc = 'm';
    depth =  1;
    // return( "lM1" );
  } else if constexpr (std::is_same<T, std::map<_L,_I>>::value) {
    dtype = 'L';
    struc = 'm';
    depth =  1;
    // return( "LM1" );
  } else if constexpr (std::is_same<T, std::map<_q,_I>>::value) {
    dtype = 'q';
    struc = 'm';
    depth =  1;
    // return( "qM1" );
  } else if constexpr (std::is_same<T, std::map<_Q,_I>>::value) {
    dtype = 'Q';
    struc = 'm';
    depth =  1;
    // return( "QM1" );
    
    /*
     * Vector<Vector> Class Match::
     * 
     */
  } else if constexpr (std::is_same<T, std::vector<std::vector<_c>>>::value) {
    dtype = 'c';
    struc = 'v';
    depth =  2;
    // return( "cV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_S>>>::value) {
    dtype = 'S';
    struc = 'v';
    depth =  2;
    // return( "SV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_b>>>::value) {
    dtype = 'b';
    struc = 'v';
    depth =  2;
    // return( "bV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_B>>>::value) {
    dtype = 'B';
    struc = 'v';
    depth =  2;
    // return( "BV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_h>>>::value) {
    dtype = 'h';
    struc = 'v';
    depth =  2;
    // return( "hV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_H>>>::value) {
    dtype = 'H';
    struc = 'v';
    depth =  2;
    // return( "HV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_i>>>::value) {
    dtype = 'i';
    struc = 'v';
    depth =  2;
    // return( "iV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_I>>>::value) {
    dtype = 'I';
    struc = 'v';
    depth =  2;
    // return( "IV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_l>>>::value) {
    dtype = 'l';
    struc = 'v';
    depth =  2;
    // return( "lV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_L>>>::value) {
    dtype = 'L';
    struc = 'v';
    depth =  2;
    // return( "LV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_q>>>::value) {
    dtype = 'q';
    struc = 'v';
    depth =  2;
    // return( "qV2" );
  } else if constexpr (std::is_same<T, std::vector<std::vector<_Q>>>::value) {
    dtype = 'Q';
    struc = 'v';
    depth =  2;
    // return( "QV2" );
    
    /*
     * Map<Vector> Class Match::
     * 
     */
  } else if constexpr (std::is_same<T, std::map<std::vector<_c>,_I>>::value) {
    dtype = 'c';
    struc = 'm';
    depth =  2;
    // return( "cMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_S>,_I>>::value) {
    dtype = 'S';
    struc = 'm';
    depth =  2;
    // return( "SMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_b>,_I>>::value) {
    dtype = 'b';
    struc = 'm';
    depth =  2;
    // return( "bMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_B>,_I>>::value) {
    dtype = 'B';
    struc = 'm';
    depth =  2;
    // return( "BMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_h>,_I>>::value) {
    dtype = 'h';
    struc = 'm';
    depth =  2;
    // return( "hMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_H>,_I>>::value) {
    dtype = 'H';
    struc = 'm';
    depth =  2;
    // return( "HMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_i>,_I>>::value) {
    dtype = 'i';
    struc = 'm';
    depth =  2;
    // return( "iMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_I>,_I>>::value) {
    dtype = 'I';
    struc = 'm';
    depth =  2;
    // return( "IMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_l>,_I>>::value) {
    dtype = 'l';
    struc = 'm';
    depth =  2;
    // return( "lMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_L>,_I>>::value) {
    dtype = 'L';
    struc = 'm';
    depth =  2;
    // return( "LMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_q>,_I>>::value) {
    dtype = 'q';
    struc = 'm';
    depth =  2;
    // return( "qMV" );
  } else if constexpr (std::is_same<T, std::map<std::vector<_Q>,_I>>::value) {
    dtype = 'Q';
    struc = 'm';
    depth =  2;
    // return( "QMV" );
  };
  
  assoc_stack_type x = { dtype, struc, depth };
  
  return( x );
};

};

#endif /* assoc_stack.h */

#ifndef __FAST_OSTREAM_H__
#define __FAST_OSTREAM_H__

#include <iomanip>
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <cstring>

namespace idat {

template<unsigned short STACK_SIZE>
class stack_buf
{
public:
  static const unsigned short stack_size = STACK_SIZE;
  stack_buf() :_v(), _stack_size(0) {}
  ~stack_buf() = default;
  stack_buf(const stack_buf& other):stack_buf(other, delegate_copy_move {})
  {}
  
  stack_buf(stack_buf&& other):stack_buf(other, delegate_copy_move {})
  {
    other.clear();
  }
  template<class T1>
  stack_buf& operator=(T1&& other)
  {
    _stack_size = other._stack_size;
    if (other.vector_used())
      _v = std::forward<T1>(other)._v;
    else
      std::copy_n(other._stack_array.begin(), other._stack_size, _stack_array.begin());
    return *this;
  }
  
  void append(const char* buf, std::size_t buf_size)
  {
    //If we are aleady using _v, forget about the stack
    if (vector_used())
    {
      _v.insert(_v.end(), buf, buf + buf_size);
    }
    //Try use the stack
    else
    {
      if (_stack_size + buf_size <= STACK_SIZE)
      {
        std::memcpy(&_stack_array[_stack_size], buf, buf_size);
        _stack_size += buf_size;
      }
      //Not enough stack space. Copy all to _v
      else
      {
        _v.reserve(_stack_size + buf_size);
        _v.insert(_v.end(), _stack_array.begin(), _stack_array.begin() + _stack_size);
        _v.insert(_v.end(), buf, buf + buf_size);
      }
    }
  }       
  void clear()
  {
    _stack_size = 0;
    _v.clear();
  }
  
  const char* data() const
  {
    if (vector_used())
      return _v.data();
    else
      return _stack_array.data();
  }
  
  std::size_t size() const
  {
    if (vector_used())
      return _v.size();
    else
      return _stack_size;
  }
  
private:
  struct delegate_copy_move {};
  template<class T1>
  stack_buf(T1&& other, delegate_copy_move)
  {
    _stack_size = other._stack_size;
    if (other.vector_used())
      _v = std::forward<T1>(other)._v;
    else
      std::copy_n(other._stack_array.begin(), other._stack_size, _stack_array.begin());
  }
  
  inline bool vector_used() const
  {
    return !_v.empty();
  }
  
  std::vector<char> _v;
  std::array<char, STACK_SIZE> _stack_array;
  std::size_t _stack_size;
};

//
// Custom std::streambuf
//
class stack_devicebuf :public std::streambuf
{
public:
  static const unsigned short stack_size = 256;
  using stackbuf_t = stack_buf<stack_size>;
  
  stack_devicebuf() = default;
  ~stack_devicebuf() = default;
  
  stack_devicebuf(const stack_devicebuf& other) :std::basic_streambuf<char>(), _stackbuf(other._stackbuf)
  {}
  
  stack_devicebuf(stack_devicebuf&& other):
    std::basic_streambuf<char>(),
    _stackbuf(std::move(other._stackbuf))
    {
      other.clear();
    }
  
  stack_devicebuf& operator=(stack_devicebuf other)
  {
    std::swap(_stackbuf, other._stackbuf);
    return *this;
  }
  
  const stackbuf_t& buf() const
  {
    return _stackbuf;
  }
  std::size_t size() const
  {
    return _stackbuf.size();
  }
  
  void clear()
  {
    _stackbuf.clear();
  }
  
protected:
  // copy the given buffer into the accumulated fast buffer
  std::streamsize xsputn(const char_type* s, std::streamsize count) override
  {
    _stackbuf.append(s, static_cast<unsigned int>(count));
    return count;
  }
  
  int_type overflow(int_type ch) override
  {
    if (traits_type::not_eof(ch))
    {
      char c = traits_type::to_char_type(ch);
      xsputn(&c, 1);
    }
    return ch;
  }
private:
  stackbuf_t _stackbuf;
};

//
// fast ostringstream
//
class fast_oss :public std::ostream
{
public:
  fast_oss() :std::ostream(&_dev) {}
  ~fast_oss() = default;
  
  fast_oss(const fast_oss& other) :std::basic_ios<char>(), std::ostream(&_dev), _dev(other._dev)
  {}
  
  fast_oss(fast_oss&& other) :std::basic_ios<char>(), std::ostream(&_dev), _dev(std::move(other._dev))
  {
    other.clear();
  }
  
  
  fast_oss& operator=(fast_oss other)
  {
    swap(*this, other);
    return *this;
  }
  
  void swap(fast_oss& first, fast_oss& second) // nothrow
  {
    using std::swap;
    swap(first._dev, second._dev);
  }
  
  std::string str()
  {
    auto& buffer = _dev.buf();
    const char*data = buffer.data();
    return std::string(data, data+buffer.size());
  }
  
  const stack_devicebuf::stackbuf_t& buf() const
  {
    return _dev.buf();
  }
  
  std::size_t size() const
  {
    return _dev.size();
  }    
  void clear()
  {
    _dev.clear();
  }
  
private:
  stack_devicebuf _dev;
};

};

#endif /* fast_ostream.h */

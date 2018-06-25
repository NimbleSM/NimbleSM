#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace nimble
{
namespace quanta
{
namespace out
{
namespace internal
{
constexpr size_t max_uchar_digits     = (int)(sizeof(char) * 2.4082399653118) + 1;
constexpr size_t max_ushort_digits    = (int)(sizeof(short) * 2.4082399653118) + 1;
constexpr size_t max_uint_digits      = (int)(sizeof(int) * 2.4082399653118) + 1;
constexpr size_t max_ulong_digits     = (int)(sizeof(long) * 2.4082399653118) + 1;
constexpr size_t max_ulonglong_digits = (int)(sizeof(long long) * 2.4082399653118) + 1;

constexpr size_t max_schar_digits     = max_uchar_digits + 1;
constexpr size_t max_sshort_digits    = max_ushort_digits + 1;
constexpr size_t max_sint_digits      = max_uint_digits + 1;
constexpr size_t max_slong_digits     = max_ulong_digits + 1;
constexpr size_t max_slonglong_digits = max_ulong_digits + 1;

template<class int_t>
struct num_digits;
template<>
struct num_digits<char>
{
  static constexpr bool is_signed = true;
  static constexpr size_t digits  = max_schar_digits;
};
template<>
struct num_digits<signed char>
{
  static constexpr bool is_signed = true;
  static constexpr size_t digits  = max_schar_digits;
};
template<>
struct num_digits<short>
{
  static constexpr bool is_signed = true;
  static constexpr size_t digits  = max_sshort_digits;
};
template<>
struct num_digits<int>
{
  static constexpr bool is_signed = true;
  static constexpr size_t digits  = max_sint_digits;
};
template<>
struct num_digits<long>
{
  static constexpr bool is_signed = true;
  static constexpr size_t digits  = max_slong_digits;
};
template<>
struct num_digits<long long>
{
  static constexpr bool is_signed = true;
  static constexpr size_t digits  = max_slonglong_digits;
};

template<>
struct num_digits<unsigned char>
{
  static constexpr bool is_signed = false;
  static constexpr size_t digits  = max_uchar_digits;
};
template<>
struct num_digits<unsigned short>
{
  static constexpr bool is_signed = false;
  static constexpr size_t digits  = max_ushort_digits;
};
template<>
struct num_digits<unsigned int>
{
  static constexpr bool is_signed = false;
  static constexpr size_t digits  = max_uint_digits;
};
template<>
struct num_digits<unsigned long>
{
  static constexpr bool is_signed = false;
  static constexpr size_t digits  = max_ulong_digits;
};
template<>
struct num_digits<unsigned long long>
{
  static constexpr bool is_signed = false;
  static constexpr size_t digits  = max_ulonglong_digits;
};

template<class int_t>
std::string& append_primitive_unsigned_int_to(std::string& s, int_t value)
{
  std::array<char, num_digits<int_t>::digits> buffer;
  auto buffer_iter = buffer.begin() + buffer.size();
  do
  {
    int_t old_value = value;
    value /= 10;
    int digit = old_value - (value * 10);
    --buffer_iter;
    *buffer_iter = '0' + digit;
  } while (value);
  s.append(buffer_iter, buffer.end());
  return s;
}
template<class int_t>
std::string& append_primitive_signed_int_to(std::string& s, int_t value)
{
  std::array<char, num_digits<int_t>::digits> buffer;
  auto buffer_iter = buffer.begin() + buffer.size();
  if (value >= 0)
  {
    do
    {
      int_t old_value = value;
      value /= 10;
      int digit = old_value - (value * 10);
      --buffer_iter;
      *buffer_iter = '0' + digit;
    } while (value);
  }
  else
  {
    do
    {
      int_t old_value = value;
      value /= 10;
      int digit = (value * 10) - old_value;
      --buffer_iter;
      *buffer_iter = '0' + digit;
    } while (value);
    --buffer_iter;
    *buffer_iter = '-';
  }
  s.append(buffer_iter, buffer.end());
  return s;
}
}   // namespace internal
std::string& append_to(std::string& s, const std::string& val);
std::string& append_to(std::string& s, char value);
std::string& append_to(std::string& s, int value);
std::string& append_to(std::string& s, long value);
std::string& append_to(std::string& s, long long value);
std::string& append_to(std::string& s, unsigned short value);
std::string& append_to(std::string& s, unsigned int value);
std::string& append_to(std::string& s, unsigned long value);
std::string& append_to(std::string& s, unsigned long long value);
std::string& append_to(std::string& s, float value);
std::string& append_to(std::string& s, double value);
std::string& append_to(std::string& s, long double value);
template<class A, class B>
std::string& append_to(std::string& s, const std::pair<A, B>& _pair);

template<class elem_t, class alloc_t>
std::string& append_to(std::string& s, const std::vector<elem_t, alloc_t>& vector);
template<class elem_t, class alloc_t>
std::string& append_to(std::string& s,
                       const std::vector<elem_t, alloc_t>& vector,
                       const std::string& _start,
                       const std::string& _end,
                       const std::string& separator);
template<class elem_t, class alloc_t, class append_method_t>
std::string& append_to(std::string& s,
                       const std::vector<elem_t, alloc_t>& vector,
                       const std::string& _start,
                       const std::string& _end,
                       const std::string& separator,
                       append_method_t&& custom_append_to);

template<class elem_t, size_t N>
std::string& append_to(std::string& s, const std::array<elem_t, N>& vector);
template<class elem_t, size_t N>
std::string& append_to(std::string& s,
                       const std::array<elem_t, N>& vector,
                       const std::string& _start,
                       const std::string& _end,
                       const std::string& separator);
template<class elem_t, size_t N, class append_method_t>
std::string& append_to(std::string& s,
                       const std::array<elem_t, N>& vector,
                       const std::string& _start,
                       const std::string& _end,
                       const std::string& separator,
                       append_method_t&& custom_append_to);

template<class iter_t, class sentry_t>
std::string& append_range_to(std::string& s,
                             iter_t _start,
                             sentry_t _end,
                             const std::string& separator);
template<class iter_t, class sentry_t, class append_method_t>
std::string& append_range_to(std::string& s,
                             iter_t _start,
                             sentry_t _end,
                             const std::string& separator,
                             append_method_t&& custom_append_to);

std::string& append_to(std::string& s, const std::string& val) { return s.append(val); }
std::string& append_to(std::string& s, char value)
{
  s.push_back(value);
  return s;
}
std::string& append_to(std::string& s, short value)
{
  return internal::append_primitive_signed_int_to(s, value);
}
std::string& append_to(std::string& s, int value)
{
  return internal::append_primitive_signed_int_to(s, value);
}
std::string& append_to(std::string& s, long value)
{
  return internal::append_primitive_signed_int_to(s, value);
}
std::string& append_to(std::string& s, long long value)
{
  return internal::append_primitive_signed_int_to(s, value);
}
std::string& append_to(std::string& s, unsigned short value)
{
  return internal::append_primitive_unsigned_int_to(s, value);
}
std::string& append_to(std::string& s, unsigned int value)
{
  return internal::append_primitive_unsigned_int_to(s, value);
}
std::string& append_to(std::string& s, unsigned long value)
{
  return internal::append_primitive_unsigned_int_to(s, value);
}
std::string& append_to(std::string& s, unsigned long long value)
{
  return internal::append_primitive_unsigned_int_to(s, value);
}
std::string& append_to(std::string& s, float value) { return s += std::to_string(value); }
std::string& append_to(std::string& s, double value) { return s += std::to_string(value); }
std::string& append_to(std::string& s, long double value) { return s += std::to_string(value); }

template<class A, class B>
std::string& append_to(std::string& s, const std::pair<A, B>& _pair) {
  s.push_back('[');
  append_to(s, _pair.first);
  s.push_back(',');
  s.push_back(' ');
  append_to(s, _pair.second);
  s.push_back(']');
  return s;
}

template<class elem_t, class alloc_t>
std::string& append_to(std::string& s, const std::vector<elem_t, alloc_t>& vector)
{
  s.push_back('[');
  append_range_to(s, vector.begin(), vector.end(), ", ");
  s.push_back(']');
  return s;
}
template<class elem_t, class alloc_t>
std::string& append_to(std::string& s,
                       const std::vector<elem_t, alloc_t>& vector,
                       const std::string& _start,
                       const std::string& _end,
                       const std::string& separator)
{
  s.append(_start);
  append_range_to(s, vector.begin(), vector.end(), separator);
  return s.append(_end);
}
template<class elem_t, class alloc_t, class append_method_t>
std::string& append_to(std::string& s,
                       const std::vector<elem_t, alloc_t>& vector,
                       const std::string& _start,
                       const std::string& _end,
                       const std::string& separator,
                       append_method_t&& custom_append_to)
{
  s.append(_start);
  append_range_to(s, vector.begin(), vector.end(), separator, custom_append_to);
  return s.append(_end);
}
template<class elem_t, size_t N>
std::string& append_to(std::string& s, const std::array<elem_t, N>& vector)
{
  s.push_back('[');
  append_range_to(s, vector.begin(), vector.end(), ", ");
  s.push_back(']');
  return s;
}
template<class elem_t, size_t N>
std::string& append_to(std::string& s,
                       const std::array<elem_t, N>& vector,
                       const std::string& _start,
                       const std::string& _end,
                       const std::string& separator)
{
  s.append(_start);
  append_range_to(s, vector.begin(), vector.end(), separator);
  return s.append(_end);
}
template<class elem_t, size_t N, class append_method_t>
std::string& append_to(std::string& s,
                       const std::array<elem_t, N>& vector,
                       const std::string& _start,
                       const std::string& _end,
                       const std::string& separator,
                       append_method_t&& custom_append_to)
{
  s.append(_start);
  append_range_to(s, vector.begin(), vector.end(), separator, custom_append_to);
  return s.append(_end);
}

template<class iter_t, class sentry_t>
std::string& append_range_to(std::string& s,
                             iter_t _start,
                             sentry_t _end,
                             const std::string& separator)
{
  if (_start == _end)
  {
    return s;
  }
  for (;;)
  {
    append_to(s, *_start);
    ++_start;
    if (_start == _end)
      break;
    append_to(s, separator);
  }
  return s;
}
template<class iter_t, class sentry_t, class append_method_t>
std::string& append_range_to(std::string& s,
                             iter_t _start,
                             sentry_t _end,
                             const std::string& separator,
                             append_method_t&& custom_append_to)
{
  if (_start == _end)
  {
    return s;
  }
  for (;;)
  {
    custom_append_to(s, *_start);
    ++_start;
    if (_start == _end)
      break;
    append_to(s, separator);
  }
  return s;
}
}   // namespace out
}   // namespace quanta
}   // namespace nimble
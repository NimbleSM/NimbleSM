#pragma once 

#if __cpp_constexpr >= 201304
#define CONSTEXPR_CPP14 constexpr /* marks a function as constexpr iff c++14 (or above) constexpr is supported */
#else 
#define CONSTEXPR_CPP14 /* marks a function as constexpr iff c++14 (or above) constexpr is supported */
#endif

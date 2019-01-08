#pragma once
template <class First, class... Rest>
using FirstOf = First;

template<class T> T iof();
template<class Range>
using ElemType = typename std::decay<decltype(std::begin(iof<Range>()))>::type;

template<class Func, class... Input>
using OutputType = decltype(iof<Func>()(iof<Input>()...));

template<class Func, class... Input>
using DecayedOutputType = typename std::decay<decltype(std::begin(iof<Range>()))>::type;
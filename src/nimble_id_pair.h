/*
//@HEADER
// ************************************************************************
//
//                                NimbleSM
//                             Copyright 2018
//   National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
// retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
// NO EVENT SHALL NTESS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact David Littlewood (djlittl@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef NIMBLE_ID_PAIR_H
#define NIMBLE_ID_PAIR_H

#include <algorithm>

#ifndef NIMBLE_HAVE_DARMA
  #include <vector>
#endif

struct id_pair;
struct id_pair_minor;
struct id_pair {
	int major;
	int minor;
	id_pair() = default;
	id_pair(const id_pair&) = default;
	id_pair(id_pair&&) = default;
	id_pair(int majorid, int minorid) : major {majorid}, minor { minorid } {};
	//minor is ignored for comparison purposes.
	bool operator >(const id_pair& other) const {
		return major > other.major;
	}
	bool operator >=(const id_pair& other) const {
		return major >= other.major;
	}
	bool operator <(const id_pair& other) const {
		return major < other.major;
	}
	bool operator <=(const id_pair& other) const {
		return major <= other.major;
	}
	bool operator ==(const id_pair& other) const {
		return major == other.major;
	}
	bool operator !=(const id_pair& other) const {
		return major != other.major;
	}
	//minor is ignored for comparison purposes.
	bool operator >(int _major) const {
		return major > _major;
	}
	bool operator >=(int _major) const {
		return major >= _major;
	}
	bool operator <(int _major) const {
		return major < _major;
	}
	bool operator <=(int _major) const {
		return major <= _major;
	}
	bool operator ==(int _major) const {
		return major == _major;
	}
	bool operator !=(int _major) const {
		return major != _major;
	}
	id_pair& operator = (const id_pair& other) {
		major = other.major;
		minor = other.minor;
		return *this;
	}
	id_pair& operator = (id_pair&& other) {
		major = other.major;
		minor = other.minor;
		return *this;
	}
	template<class iteratorT>
	static std::vector<int> GetMinorIDs(iteratorT begin, iteratorT end) {
		std::vector<int> vect { std::distance(begin, end) };
		for(int& i : vect) {
			i = (*begin).minor;
			++begin;
		}
		return vect;
	}
	template<class beginT>
	static std::vector<int> GetMinorIDs(beginT begin, size_t size) {
		std::vector<int> vect = std::vector<int>(size);
		for(int& i : vect) {
			i = (*begin).minor;
			++begin;
		}
		return vect;
	}
	template<class iteratorT>
	static std::vector<int> GetMajorIDs(iteratorT begin, iteratorT end) {
		std::vector<int> vect { std::distance(begin, end) };
		for(int& i : vect) {
			i = (*begin).major;
			++begin;
		}
		return vect;
	}
	template<class beginT>
	static std::vector<int> GetMajorIDs(beginT begin, size_t size) {
		std::vector<int> vect = std::vector<int>(size);
		for(int& i : vect) {
			i = (*begin).major;
			++begin;
		}
		return vect;
	}
	static void SortByGlobalID(id_pair* begin, id_pair* end) {
		std::sort(begin, end);
	}
	static void SortByLocalID(id_pair* begin, id_pair* end) {
		//Does a reinterpret cast to being locally preferenced id pairs.
		std::sort((id_pair_minor*)begin, (id_pair_minor*)end);
	}
	static void StableSortByGlobalID(id_pair* begin, id_pair* end) {
		std::stable_sort(begin, end);
	}
	static void StableSortByLocalID(id_pair* begin, id_pair* end) {
		//Does a reinterpret cast to being locally preferenced id pairs.
		std::stable_sort((id_pair_minor*)begin, (id_pair_minor*)end);
	}
	template<class listT> static void SortByGlobalID(listT& list) {
		SortByGlobalID(&list.front(), &list.back() + 1);
	}
	template<class listT> static void StableSortByGlobalID(listT& list) {
		StableSortByGlobalID(&list.front(), &list.back() + 1);
	}
	template<class listT> static void SortByLocalID(listT& list) {
		SortByLocalID(&list.front(), &list.back() + 1);
	}
	template<class listT> static void StableSortByLocalID(listT& list) {
		StableSortByLocalID(&list.front(), &list.back() + 1);
	}
};
struct id_pair_minor : id_pair {
	using id_pair::id_pair;
	id_pair_minor(const id_pair_minor&) = default;
	id_pair_minor(id_pair_minor&&) = default;
	bool operator >(const id_pair_minor& other) const {
		return minor > other.minor;
	}
	bool operator >=(const id_pair_minor& other) const {
		return minor >= other.minor;
	}
	bool operator <(const id_pair_minor& other) const {
		return minor < other.minor;
	}
	bool operator <=(const id_pair_minor& other) const {
		return minor <= other.minor;
	}
	bool operator ==(const id_pair_minor& other) const {
		return minor == other.minor;
	}
	bool operator !=(const id_pair_minor& other) const {
		return minor != other.minor;
	}
	bool operator >(int _minor) const {
		return minor > _minor;
	}
	bool operator >=(int _minor) const {
		return minor >= _minor;
	}
	bool operator <(int _minor) const {
		return minor < _minor;
	}
	bool operator <=(int _minor) const {
		return minor <= _minor;
	}
	bool operator ==(int _minor) const {
		return minor == _minor;
	}
	bool operator !=(int _minor) const {
		return minor != _minor;
	}
	id_pair_minor& operator = (const id_pair_minor& other) {
		major = other.major;
		minor = other.minor;
		return *this;
	}
	id_pair_minor& operator = (id_pair_minor&& other) {
		major = other.major;
		minor = other.minor;
		return *this;
	}
};

#endif // NIMBLE_ID_PAIR_H

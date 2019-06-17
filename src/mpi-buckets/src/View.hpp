#pragma once
#include <iterator>
#include <type_traits>
#include "defines.hpp"
#include "meta.hpp"


template <class Iter>
struct View {
   private:
    constexpr static auto getIter() noexcept -> Iter;
    constexpr static auto getrIter() noexcept -> Iter &;
    constexpr static auto getcrIter() noexcept -> Iter const &;

   public:
    Iter scan;
    Iter ending;
    View()             = default;
    View(View const &) = default;
    View(View &&)      = default;

    using iterator_tag = typename std::iterator_traits<Iter>::iterator_category;

    // using value_type = decltype(head());
    constexpr static bool isRandomAccess =
        std::is_same<std::random_access_iterator_tag, iterator_tag>::value;
    constexpr static bool isDerefNX   = noexcept(*getIter());
    constexpr static bool isCopyNX    = noexcept(Iter{getcrIter()});
    constexpr static bool isTailNX    = isCopyNX && noexcept(Iter{getcrIter() + 1});
    constexpr static bool isHasNX     = noexcept(getcrIter() != getcrIter());
    constexpr static bool isEmptyNX   = noexcept(getcrIter() == getcrIter());
    constexpr static bool isAdvanceNX = noexcept(++getrIter());
    constexpr static bool isAtNX      = noexcept(getcrIter()[0]);

    /**
     * @brief For some reason std::distance isn't noexcept, so this returns true
     * if std::distance is random-access
     *
     */
    constexpr static bool isSizeNX = isRandomAccess;

    using index_t = decltype(scan - ending);

    /* Always mutable */
    /**
     * @brief Advances scan by 1; basically the same as doing *this = tail
     *
     * @return decltype(++scan) returns a reference to the new head
     */
    CONSTEXPR_CPP14 auto advance() noexcept(noexcept(++scan)) -> decltype(++scan) {
        return ++scan;
    }

    /**
     * @brief indexes into the view; only enabled if Iter is random-access.
     *
     * @param index index to access
     * @return std::enable_if<isRandomAccess, decltype(scan[index])>
     */
    CONSTEXPR_CPP14 auto operator[](index_t index) noexcept(noexcept(scan[index]))
        -> std::enable_if<isRandomAccess, decltype(scan[index])> {
        return scan[index];
    }
    constexpr auto operator[](index_t index) const noexcept(noexcept(scan[index]))
        -> std::enable_if<isRandomAccess, decltype(scan[index])> {
        return scan[index];
    }
    CONSTEXPR_CPP14 auto head() noexcept(noexcept(*scan)) -> decltype(*scan) {
        return *scan;
    }
    constexpr auto head() const noexcept(noexcept(*scan)) -> decltype(*scan) {
        return *scan;
    }
    /**
     * @brief Returns the beginning of the view
     *
     * @return decltype(Iter{scan})
     */
    CONSTEXPR_CPP14 auto begin() noexcept(isCopyNX) -> decltype(Iter{scan}) {
        return scan;
    }
    /**
     * @brief Returns the end of the view
     *
     * @return decltype(Iter{scan})
     */
    CONSTEXPR_CPP14 auto end() noexcept(isCopyNX) -> decltype(Iter{scan}) {
        return ending;
    }
    /**
     * @brief Returns the beginning of the view
     *
     * @return decltype(Iter{scan})
     */
    constexpr auto begin() const noexcept(isCopyNX) -> decltype(Iter{scan}) {
        return scan;
    }
    /**
     * @brief Returns the end of the view
     *
     * @return decltype(Iter{scan})
     */
    constexpr auto end() const noexcept(isCopyNX) -> decltype(Iter{scan}) {
        return ending;
    }

    /**
     * @brief returns a View of everything after the head
     *
     * @return View
     */
    constexpr auto tail() const noexcept(isTailNX) -> View {
        return View{scan + 1, ending};
    }
    /**
     * @brief checks if there are any elements in the view
     *
     * @return decltype(scan != ending)
     */
    constexpr auto has() const noexcept(isHasNX) -> decltype(scan != ending) {
        return scan != ending;
    }
    /**
     * @brief checks if the view is empty
     *
     * @return decltype(scan == ending)
     */
    constexpr auto empty() const noexcept(isHasNX) -> decltype(scan == ending) {
        return scan == ending;
    }
    /**
     * @brief Returns the size
     *
     * @return decltype(std::distance(scan, ending))
     */
    constexpr auto size() const noexcept(isSizeNX)
        -> decltype(std::distance(scan, ending)) {
        return std::distance(scan, ending);
    }
};
/* 
Template deduction guides are a c++17 feature; this enables them
only if the code is being compiled with a version greater than c++17
*/
#if __cplusplus >= 201703
/**
 * @brief Specifies how to deduce the template parameters of View given some inputs
 * to the constructor
 * 
 * @tparam Iterator
 */
template <class Iterator>
View(Iterator, Iterator)->View<Iterator>;
#endif

namespace ViewDeets {
/**
 * @brief Creates a view of the given iterable list
 * 
 * @tparam List 
 * @param list 
 * @return View<decltype(std::begin(list))> 
 */
template <class List>
constexpr auto of(List && list) 
    noexcept(
        noexcept(
            View<decltype(std::begin(list))>{std::begin(list), std::end(list)}
        )
    ) -> View<decltype(std::begin(list))> {
    return View<decltype(std::begin(list))>{std::begin(list), std::end(list)};
}
}  // namespace ViewDeets

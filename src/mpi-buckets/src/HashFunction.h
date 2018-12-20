#pragma once
#include <cstddef>
#include <cstdint>
#include <random>
#include <type_traits>
#include "logic.h"
/*Supported by:
icc 18.0.0 or later
gcc 4.7.1 or later
clang 3.5 or later
msvc 19.14 or later
*/

/**
 * @brief Implements a hash function that can hash 1 or more integer-valued inputs
 * https://en.wikipedia.org/wiki/Universal_hashing#Avoiding_modular_arithmetic_2
 *
 * Multiple values are hashed as
 * h(v1) for 1 value
 * h(h(v1) ^ v2) for 2 values
 * h(h(h(v1) ^ v2) ^ v3) for 3 values
 * h(h(h(h(v1) ^ v2) ^ v3) ^ v4) for 4 values, and so on
 * where h is the hash function.
 *
 * this is done so that hash(a, b) != hash(b, a)
 */
class HashFunction {
   private:
    /**
     * @brief handles the trivial case for the hash helper (don't do anything, return the
     * accumulator)
     *
     * @param accumulator keeps track of the hash so far
     * @return uint64_t hashed value
     */

    constexpr auto hash_helper(uint64_t accumulator) const noexcept -> uint64_t {
        return accumulator;
    }
    /**
     * @brief handles hashing a single input and an accumulator
     *
     * @param accumulator keeps track of the hash so far
     * @param input value to hash with the accumulator
     * @return uint64_t hashed value
     */
    constexpr auto hash_helper(uint64_t accumulator, uint64_t input) const noexcept
        -> uint64_t {
        return hash(accumulator ^ input);
    }
    /**
     * @brief Handles the general case hashing of multiple integral inputs
     *
     * @tparam In Types of the other integral inputs
     * @param accumulator keeps track of the hash so far
     * @param input input to hash this go-around
     * @param rest rest of the inputs to hash
     * @return uint64_t hashed value
     */
    template <class... In>
    constexpr auto hash_helper(uint64_t accumulator, uint64_t input, In... rest) const
        noexcept -> uint64_t {
        return hash_helper(hash(accumulator ^ input), rest...);
    }
    /**
     * @brief Mashes the higher and lower bits - used to prevent hashes that
     * follow a linear sequence. If the mashing step doesn't occur,
     * then (int)hash(1) - (int)hash(0) == (int)hash(2) - (int)hash(1)
     * which is bad
     *
     * @param input
     * @return uint64_t
     */
    constexpr auto mash(uint64_t input) const noexcept -> uint64_t {
        return (input ^ (input << 32)) >> 32;
        //return input >> 32;
    }

   public:
    // default_a and default_b were chosen as the first and second values
    // produced by the default_constructed std::mt19937_64.
    // default_a should be odd and default_b's input doesn't matter
    constexpr static uint64_t default_a = 14514284786278117030ull | 1;
    constexpr static uint64_t default_b = 4620546740167642908ull;

    uint64_t a = default_a;
    uint64_t b = default_b;

    constexpr HashFunction()                    = default;
    constexpr HashFunction(HashFunction const&) = default;
    constexpr HashFunction(HashFunction&&)      = default;

    auto operator=(HashFunction const&) -> HashFunction& = default;
    auto operator=(HashFunction &&) -> HashFunction& = default;

    ~HashFunction() = default;
    /**
     * @brief Construct a new Hash Function object
     *
     * @param a - first parameter of hash function. Must be an odd number; constructor
     * automatically automatically converts into an odd number using bitwise or
     * @param b - second parameter of hash function.
     */
    constexpr HashFunction(uint64_t a, uint64_t b) noexcept : a(a | 1), b(b) {}
    /**
     * @brief Construct a new Hash Function object using a random number generator
     *
     * @tparam RNG type representing the random number generator
     * @param gen the random number generator
     */
    template <class RNG>
    HashFunction(RNG&& gen) noexcept(noexcept(gen())) {
        constexpr static uint64_t min = 1;
        constexpr static uint64_t max = std::numeric_limits<uint64_t>::max();
        std::uniform_int_distribution<uint64_t> sampleFrom{min, max};
        a = sampleFrom(gen) | 1;
        b = sampleFrom(gen);
    }

    /**
     * @brief Hashes the input according to the parameters of the hash function
     *
     * @param input the value to hash
     * @return uint64_t the hashed value, an integer in the range [0, 2^32)
     */
    constexpr auto hash(uint64_t input) const noexcept -> uint64_t {
        return mash(((uint64_t)a * (uint64_t)input + b));
    }
    /**
     * @brief Hashes multiple inputs into a single output according to the parameters of
     * the hash function
     *
     * @tparam In A type implicitly convertible to uint64_t
     * @param input a list of inputs
     * @param rest
     * @return uint64_t
     */
    template <class... In>
    constexpr auto hash(uint64_t input, In... rest) const noexcept -> uint64_t {
        static_assert(all_true(std::is_integral<In>::value...),
                      "Inputs to hash function must all be integral types");
        return hash_helper(0, input, (uint64_t)rest...);
    }
    /**
     * @brief Hashes the input according to the parameters of the hash function
     *
     * @param input the value to hash
     * @return uint64_t the hashed value, an integer in the range [0, 2^32)
     */
    constexpr auto operator()(uint64_t input) const noexcept -> uint64_t {
        return hash(input);
    }
    /**
     * @brief Hashes multiple inputs into a single output according to the parameters of
     * the hash function
     *
     * @tparam In A type implicitly convertible to uint64_t
     * @param input a list of inputs
     * @param rest
     * @return uint64_t
     */
    template <class... In>
    constexpr auto operator()(uint64_t first, In... rest) const noexcept -> uint64_t {
        static_assert(all_true(std::is_integral<In>::value...),
                      "Inputs to hash function must all be integral types");
        return hash_helper(0, first, (uint64_t)rest...);
    }
};

/**
 * @brief This block contains a set of checks to verify that the hash function
 * implememtation has nice properties. It basically contains compile-time
 * unit tests :)
 * 
 */
namespace {

constexpr static HashFunction DefaultHash = HashFunction{};

constexpr auto hdiff(uint64_t val, uint64_t skip) -> int64_t {
    return (int64_t)DefaultHash(val + skip) - (int64_t)DefaultHash(val); 
}

constexpr auto checkThatHashingIsNonlinear(uint64_t max, uint64_t skip) -> bool {
    return max == 0 || skip == 0
        ? true
        : hdiff(0, skip) != hdiff(max, skip) 
            && checkThatHashingIsNonlinear(max - 1, skip);
}

//Hashing should be non-commutative - in other words, it should depend on the
//order of the inputs. If it doesn't satisfy this, there will be wierd edge
//cases where the hashing fails
static_assert(
    DefaultHash(0, 1) != DefaultHash(1, 0),
    "Hashes with the same inputs in a different order should be distinct"
);

static_assert(
    DefaultHash(0) != 0, 
    "0 should not hash to 0 under the default hash"
);

//Ensures that hash(n + k) - hash(n) != hash(k) - hash(0)
//From n=1 to 50 and from k = 1 to 4
static_assert(
    checkThatHashingIsNonlinear(50, 1) &&
    checkThatHashingIsNonlinear(50, 2) &&
    checkThatHashingIsNonlinear(50, 3) && 
    checkThatHashingIsNonlinear(50, 4), 
    "hash function implementation produces linear results"
);
}  // namespace

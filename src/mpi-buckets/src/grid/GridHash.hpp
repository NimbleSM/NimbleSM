#pragma once
#include <cstdint>
#include "HashFunction.hpp"

/**
 * @brief Represents a hash function that will hash a coordinate (x, y, z)
 * by first finding the corresponding grid cell, and hashing the indicies of
 * the grid cell
 *
 */
class GridHash {
   public:
    /**
     * @brief returns the floor of the input as an integer
     *
     * @param val input to take the floor of
     * @return int64_t
     */
    constexpr static auto ifloor(double val) -> int64_t {
        /*
            NB: This function was originally written as:

                return (int64_t)floor(val);

            However that code failed to compile on clang.
            After doing research I found out that <cmath> functions can't be
            used in constexpr contexts, so I've switched to the below
            implementation.

            It returns the exact value as an int for inputs that are exact
            integers. Otherwise, it returns

                (int64_t)val if val is positive, and
                (int64_t)(val - 1) if val is negative
        */
        return (double)(int64_t)val == val
                   ? (int64_t)val
                   : val < 0 ? (int64_t)(val - 1.0) : (int64_t)val;
    }

    /**
     * @brief Used to compute a complete hash
     *
     */
    HashFunction hashCellBounds;
    /**
     * @brief The scale is the inverse of the cell size
     *
     */
    double scale = 1.0;

    GridHash()                = default;
    GridHash(GridHash const&) = default;
    GridHash(GridHash&&)      = default;

    /**
     * @brief Constructs a grid hash
     *
     */
    constexpr GridHash(HashFunction const& hashCellBounds, double cell_size) noexcept
        : hashCellBounds(hashCellBounds), scale(1. / cell_size) {}

    auto operator=(GridHash const&) -> GridHash& = default;
    auto operator=(GridHash &&) -> GridHash& = default;

    /**
     * @brief returns the amount coordinates are scaled to get the scale
     *
     * @return double
     */
    constexpr auto getScale() const noexcept -> double { return scale; }
    /**
     * @brief Get the size of grid cells
     *
     * @return double
     */
    constexpr auto getCellSize() const noexcept -> double { return 1.0 / scale; }
    /**
     * @brief Get the Hash Function used to do the hashing
     *
     * @return HashFunction
     */
    constexpr auto getHashFunction() const noexcept -> HashFunction {
        return hashCellBounds;
    }
    /**
     * @brief Returns the cell coordinate of the spacial coordinate
     * Equivilant to ifloor(coordinate / getCellSize())
     *
     * @param coordinate coordinate to find the index of
     * @return constexpr int64_t cell coordinate
     */
    constexpr int64_t index(double coordinate) const noexcept {
        return ifloor(coordinate * scale);
    }

    /**
     * @brief Hashes the grid cell given by x, y, and z
     *
     * @param cell_x x coordinate of cell
     * @param cell_y y coordinate of cell
     * @param cell_z z coordinate of cell
     * @return uint64_t
     */
    constexpr auto hash(int64_t cell_x, int64_t cell_y, int64_t cell_z) const noexcept -> uint32_t {
        return hashCellBounds((uint64_t)cell_x, (uint64_t)cell_y, (uint64_t)cell_z);
    }
    /**
     * @brief Hashes the grid cell given by x, y, and z
     *
     * @param cell_x x coordinate of cell
     * @param cell_y y coordinate of cell
     * @param cell_z z coordinate of cell
     * @return uint64_t hashed integer
     */
    constexpr auto operator()(int64_t cell_x, int64_t cell_y, int64_t cell_z) const noexcept -> uint32_t {
        return hashCellBounds((uint64_t)cell_x, (uint64_t)cell_y, (uint64_t)cell_z);
    }
    /**
     * @brief Finds the grid cell containing x, y, and z, and returns a hash
     * of the mininum corner of the grid cell
     *
     * @param x x coordinate of point
     * @param y y coordinate of point
     * @param z z coordinate of point
     * @return uint64_t
     */
    constexpr auto hash(double x, double y, double z) const noexcept -> uint32_t {
        return hash(ifloor(x * scale), ifloor(y * scale), ifloor(z * scale));
    }
    /**
     * @brief Same as calling hash member function
     *
     * @param x x coordinate of point
     * @param y y coordinate of point
     * @param z z coordinate of point
     * @return uint64_t hashed integer
     */
    constexpr auto operator()(double x, double y, double z) const noexcept -> uint32_t {
        return hash(ifloor(x * scale), ifloor(y * scale), ifloor(z * scale));
    }
};

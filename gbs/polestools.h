#pragma once
#include <vector>
#include <array>


namespace gbs {
    /**
     * Transposes a matrix of poles.
     * 
     * @param poles A 2D vector representing a matrix of poles.
     * @return A 2D vector representing the transposed matrix of poles.
     */
    template <typename T, size_t dim>
    std::vector<std::vector<std::array<T, dim>>> transpose_poles(const std::vector<std::vector<std::array<T, dim>>>& poles) {
        if (poles.empty()) return {};

        size_t nu = poles.size();
        size_t nv = poles[0].size();

        std::vector<std::vector<std::array<T, dim>>> poles_t(nv, std::vector<std::array<T, dim>>(nu));

        for (size_t i = 0; i < nu; ++i) {
            for (size_t j = 0; j < nv; ++j) {
                poles_t[j][i] = poles[i][j];
            }
        }

        return poles_t;
    }

    /**
     * Transposes a block of a matrix of poles.
     * 
     * @param input The input matrix of poles to transpose.
     * @param output The output matrix where the transposed block will be stored.
     * @param startRow The starting row index for the block.
     * @param startCol The starting column index for the block.
     * @param blockSize The size of the block to be transposed.
     */
    template <typename T, size_t dim>
    // Function to transpose a small block of the matrix
    void transpose_block(const std::vector<std::vector<std::array<T, dim>>>& input,
                        std::vector<std::vector<std::array<T, dim>>>& output,
                        size_t startRow, size_t startCol, size_t blockSize) {
        size_t endRow = std::min(startRow + blockSize, input.size());
        size_t endCol = std::min(startCol + blockSize, input[0].size());

        for (size_t i = startRow; i < endRow; ++i) {
            for (size_t j = startCol; j < endCol; ++j) {
                output[j][i] = input[i][j];
            }
        }
    }

    /**
     * Transposes a matrix of poles using block-based transposition for efficiency.
     * 
     * @param poles The matrix of poles to be transposed.
     * @param blockSize The size of each block used in the transposition.
     * @return The transposed matrix of poles.
     */
    template <typename T, size_t dim>
    auto transpose_poles(const std::vector<std::vector<std::array<T, dim>>>& poles, size_t blockSize) {

        size_t nu = poles.size();
        size_t nv = poles[0].size();

        std::vector<std::vector<std::array<T, dim>>> poles_t(nv, std::vector<std::array<T, dim>>(nu));

        for (size_t i = 0; i < nu; i += blockSize) {
            for (size_t j = 0; j < nv; j += blockSize) {
                transpose_block(poles, poles_t, i, j, blockSize);
            }
        }

        return poles_t;
    }

    /**
     * Flattens a 2D vector of poles into a 1D vector.
     * 
     * @param poles_curves The 2D vector of poles to be flattened.
     * @return A 1D vector containing all the poles.
     */
    template <typename T, size_t dim>
    auto flatten_poles(const std::vector<std::vector<std::array<T, dim>>>& poles_curves) {
        std::vector<std::array<T, dim>> flattened;

        // Calculate total size for memory reservation
        size_t totalSize = std::accumulate(poles_curves.begin(), poles_curves.end(), size_t(0),
            [](size_t sum, const auto& row) { return sum + row.size(); });
        flattened.reserve(totalSize);

        // Flatten the vector
        for (const auto& row : poles_curves) {
            std::ranges::copy(row, std::back_inserter(flattened));
        }

        return flattened;
    }
}
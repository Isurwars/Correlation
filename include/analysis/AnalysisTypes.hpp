#pragma once

#include "math/Precision.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <span>
#include <vector>

namespace correlation::analysis {

struct MaxFrames {
  int value;
};

struct StartFrame {
  size_t value;
};

struct EndFrame {
  size_t value;
};

/**
 * @struct BondCutoffRange
 * @brief Represents the squared minimum and maximum bond distance bounds for an element pair.
 */
struct alignas(64) BondCutoffRange {
  real_t min_sq{static_cast<real_t>(0.0)}; ///< Squared minimum bond distance (Å²).
  real_t max_sq{static_cast<real_t>(0.0)}; ///< Squared maximum bond distance (Å²).

  [[nodiscard]] constexpr bool isValid() const noexcept {
    return min_sq >= 0.0 && max_sq > 0.0 && min_sq <= max_sq;
  }

  constexpr bool operator==(const BondCutoffRange &other) const noexcept = default;
};

/**
 * @class BondCutoffMatrix
 * @brief Contiguous 2D matrix of element-pair bond cutoff ranges.
 *
 * Stores cutoff ranges in a single cache-line-aligned flat buffer while
 * providing RowProxy and std::span row views for seamless operator[i][j] indexing
 * with zero heap fragmentation.
 */
class BondCutoffMatrix {
public:
  /**
   * @class RowProxy
   * @brief Mutable row proxy providing element access, span conversion, and compatibility resize.
   */
  class RowProxy {
  public:
    RowProxy() = default;
    RowProxy(BondCutoffMatrix &matrix, size_t row) noexcept : matrix_(&matrix), row_(row) {}

    /** @return Reference to element at index @p col. */
    [[nodiscard]] BondCutoffRange &operator[](size_t col) const noexcept;

    /** @return Number of elements in this row. */
    [[nodiscard]] size_t size() const noexcept;

    /** @return True if row has zero elements. */
    [[nodiscard]] bool empty() const noexcept;

    /** @brief Compatibility resize ensuring matrix column count matches @p new_cols. */
    void resize(size_t new_cols);

    /** @brief Advances row offset by @p n rows. */
    void advance(std::ptrdiff_t n) noexcept {
      row_ = static_cast<size_t>(static_cast<std::ptrdiff_t>(row_) + n);
    }

    [[nodiscard]] size_t row() const noexcept { return row_; }
    [[nodiscard]] const BondCutoffMatrix *matrix() const noexcept { return matrix_; }

    explicit operator std::span<BondCutoffRange>() const noexcept;
    explicit operator std::span<const BondCutoffRange>() const noexcept;

    [[nodiscard]] BondCutoffRange *begin() const noexcept;
    [[nodiscard]] BondCutoffRange *end() const noexcept;

  private:
    BondCutoffMatrix *matrix_{nullptr};
    size_t row_{0};
  };

  /**
   * @class ConstRowProxy
   * @brief Immutable row proxy providing element access and span conversion.
   */
  class ConstRowProxy {
  public:
    ConstRowProxy() = default;
    ConstRowProxy(const BondCutoffMatrix &matrix, size_t row) noexcept
        : matrix_(&matrix), row_(row) {}

    /** @return Const reference to element at index @p col. */
    [[nodiscard]] const BondCutoffRange &operator[](size_t col) const noexcept;

    /** @return Number of elements in this row. */
    [[nodiscard]] size_t size() const noexcept;

    /** @return True if row has zero elements. */
    [[nodiscard]] bool empty() const noexcept;

    /** @brief Advances row offset by @p n rows. */
    void advance(std::ptrdiff_t n) noexcept {
      row_ = static_cast<size_t>(static_cast<std::ptrdiff_t>(row_) + n);
    }

    [[nodiscard]] size_t row() const noexcept { return row_; }
    [[nodiscard]] const BondCutoffMatrix *matrix() const noexcept { return matrix_; }

    explicit operator std::span<const BondCutoffRange>() const noexcept;

    [[nodiscard]] const BondCutoffRange *begin() const noexcept;
    [[nodiscard]] const BondCutoffRange *end() const noexcept;

  private:
    const BondCutoffMatrix *matrix_{nullptr};
    size_t row_{0};
  };

  /** @brief Default constructor creating an empty matrix. */
  BondCutoffMatrix() = default;

  /**
   * @brief Constructs an n x n square matrix.
   * @param[in] n Row and column count.
   */
  explicit BondCutoffMatrix(size_t n) : rows_(n), cols_(n), data_(n * n) {}

  /**
   * @brief Constructs a matrix with specified dimensions and default value.
   * @param[in] rows Number of rows.
   * @param[in] cols Number of columns.
   * @param[in] default_val Initial value for all matrix entries.
   */
  BondCutoffMatrix(size_t rows, size_t cols, BondCutoffRange default_val = {})
      : rows_(rows), cols_(cols), data_(rows * cols, default_val) {}

  /**
   * @brief Compatibility constructor replicating an initial row across all rows.
   * @param[in] rows Number of rows.
   * @param[in] initial_row Initial row vector to copy into every row.
   */
  BondCutoffMatrix(size_t rows, const std::vector<BondCutoffRange> &initial_row)
      : rows_(rows), cols_(initial_row.size()), data_(rows * initial_row.size()) {
    for (size_t row = 0; row < rows; ++row) {
      std::ranges::copy(initial_row, data_.begin() + static_cast<std::ptrdiff_t>(row * cols_));
    }
  }

  /**
   * @brief Initializer-list constructor for 2D row initializers.
   * @param[in] list Initializer list of vector rows.
   */
  BondCutoffMatrix(std::initializer_list<std::vector<BondCutoffRange>> list)
      : rows_(list.size()), cols_(rows_ > 0 ? list.begin()->size() : 0), data_(rows_ * cols_) {
    size_t r = 0;
    for (const auto &row : list) {
      std::ranges::copy(row, data_.begin() + static_cast<std::ptrdiff_t>(r * cols_));
      ++r;
    }
  }

  /** @return Number of rows in the matrix. */
  [[nodiscard]] size_t size() const noexcept { return rows_; }

  /** @return Number of rows in the matrix. */
  [[nodiscard]] size_t rows() const noexcept { return rows_; }

  /** @return Number of columns in the matrix. */
  [[nodiscard]] size_t cols() const noexcept { return cols_; }

  /** @return True if matrix has zero elements. */
  [[nodiscard]] bool empty() const noexcept { return data_.empty(); }

  /** @brief Clears all elements, resetting matrix dimensions to zero. */
  void clear() noexcept {
    data_.clear();
    rows_ = 0;
    cols_ = 0;
  }

  /**
   * @brief Resizes matrix rows, defaulting columns to rows if uninitialized.
   * @param[in] rows New row count.
   */
  void resize(size_t rows) {
    if (cols_ == 0) {
      cols_ = rows;
    }
    rows_ = rows;
    data_.resize(rows_ * cols_);
  }

  /**
   * @brief Compatibility resize replicating @p initial_row across all rows.
   * @param[in] rows New row count.
   * @param[in] initial_row Row to replicate.
   */
  void resize(size_t rows, const std::vector<BondCutoffRange> &initial_row) {
    cols_ = initial_row.size();
    rows_ = rows;
    data_.resize(rows_ * cols_);
    for (size_t row = 0; row < rows_; ++row) {
      std::ranges::copy(initial_row, data_.begin() + static_cast<std::ptrdiff_t>(row * cols_));
    }
  }

  /**
   * @brief Resizes matrix to specified dimensions.
   * @param[in] rows New row count.
   * @param[in] cols New column count.
   * @param[in] default_val Default value for populated cells.
   */
  void resize(size_t rows, size_t cols, BondCutoffRange default_val = {}) {
    rows_ = rows;
    cols_ = cols;
    data_.assign(rows * cols, default_val);
  }

  /**
   * @brief Ensures matrix has at least @p new_cols columns.
   * @param[in] new_cols Column count to set.
   */
  void ensureCols(size_t new_cols) {
    if (cols_ == new_cols) {
      return;
    }
    cols_ = new_cols;
    data_.resize(rows_ * cols_);
  }

  /**
   * @brief Accesses row via mutable RowProxy.
   * @param[in] row Zero-based row index.
   * @return RowProxy over the row's elements.
   */
  [[nodiscard]] RowProxy operator[](size_t row) noexcept { return {*this, row}; }

  /**
   * @brief Accesses row via immutable ConstRowProxy.
   * @param[in] row Zero-based row index.
   * @return ConstRowProxy over the row's elements.
   */
  [[nodiscard]] ConstRowProxy operator[](size_t row) const noexcept { return {*this, row}; }

  /**
   * @brief Direct 2D mutable element accessor.
   * @param[in] row Zero-based row index.
   * @param[in] col Zero-based column index.
   * @return Reference to the cutoff range.
   */
  [[nodiscard]] BondCutoffRange &operator()(size_t row, size_t col) noexcept {
    return data_[row * cols_ + col];
  }

  /**
   * @brief Direct 2D immutable element accessor.
   * @param[in] row Zero-based row index.
   * @param[in] col Zero-based column index.
   * @return Const reference to the cutoff range.
   */
  [[nodiscard]] const BondCutoffRange &operator()(size_t row, size_t col) const noexcept {
    return data_[row * cols_ + col];
  }

  /** @return Pointer to start of row @p row. */
  [[nodiscard]] BondCutoffRange *rowData(size_t row) noexcept { return data_.data() + row * cols_; }

  /** @return Const pointer to start of row @p row. */
  [[nodiscard]] const BondCutoffRange *rowData(size_t row) const noexcept {
    return data_.data() + row * cols_;
  }

  /** @return Contiguous span of all matrix elements. */
  [[nodiscard]] std::span<const BondCutoffRange> flatSpan() const noexcept { return data_; }

  /** @return Mutable contiguous span of all matrix elements. */
  [[nodiscard]] std::span<BondCutoffRange> flatSpan() noexcept { return data_; }

  /**
   * @class RowIterator
   * @brief Random-access iterator yielding RowProxy views for each row.
   */
  template <bool IsConst> class RowIterator {
  public:
    using difference_type = std::ptrdiff_t;
    using value_type = std::conditional_t<IsConst, ConstRowProxy, RowProxy>;
    using reference = value_type &;
    using pointer = value_type *;
    using iterator_category = std::random_access_iterator_tag;

    RowIterator() = default;
    RowIterator(std::conditional_t<IsConst, const BondCutoffMatrix *, BondCutoffMatrix *> matrix,
                size_t row) noexcept
        : proxy_(*matrix, row) {}

    reference operator*() noexcept { return proxy_; }
    const value_type &operator*() const noexcept { return proxy_; }
    pointer operator->() noexcept { return &proxy_; }
    const value_type *operator->() const noexcept { return &proxy_; }

    RowIterator &operator++() noexcept {
      proxy_.advance(1);
      return *this;
    }
    RowIterator operator++(int) noexcept {
      auto tmp = *this;
      proxy_.advance(1);
      return tmp;
    }
    RowIterator &operator--() noexcept {
      proxy_.advance(-1);
      return *this;
    }
    RowIterator operator--(int) noexcept {
      auto tmp = *this;
      proxy_.advance(-1);
      return tmp;
    }
    RowIterator &operator+=(difference_type n) noexcept {
      proxy_.advance(n);
      return *this;
    }
    RowIterator &operator-=(difference_type n) noexcept {
      proxy_.advance(-n);
      return *this;
    }
    bool operator==(const RowIterator &other) const noexcept {
      return proxy_.matrix() == other.proxy_.matrix() && proxy_.row() == other.proxy_.row();
    }
    auto operator<=>(const RowIterator &other) const noexcept {
      return proxy_.row() <=> other.proxy_.row();
    }

  private:
    value_type proxy_{};
  };

  using iterator = RowIterator<false>;
  using const_iterator = RowIterator<true>;

  [[nodiscard]] iterator begin() noexcept { return {this, 0}; }
  [[nodiscard]] iterator end() noexcept { return {this, rows_}; }
  [[nodiscard]] const_iterator begin() const noexcept { return {this, 0}; }
  [[nodiscard]] const_iterator end() const noexcept { return {this, rows_}; }
  [[nodiscard]] const_iterator cbegin() const noexcept { return begin(); }
  [[nodiscard]] const_iterator cend() const noexcept { return end(); }

  /**
   * @brief Equality operator for BondCutoffMatrix.
   * @param[in] other The other matrix to compare against.
   * @return True if dimensions and contents match.
   */
  bool operator==(const BondCutoffMatrix &other) const noexcept {
    return rows_ == other.rows_ && cols_ == other.cols_ && data_ == other.data_;
  }

private:
  size_t rows_{0};
  size_t cols_{0};
  std::vector<BondCutoffRange> data_;
};

inline BondCutoffRange &BondCutoffMatrix::RowProxy::operator[](size_t col) const noexcept {
  return (*matrix_)(row_, col);
}

inline size_t BondCutoffMatrix::RowProxy::size() const noexcept { return matrix_->cols(); }

inline bool BondCutoffMatrix::RowProxy::empty() const noexcept { return matrix_->cols() == 0; }

inline void BondCutoffMatrix::RowProxy::resize(size_t new_cols) {
  if (matrix_->cols() == 0 && new_cols > 0) {
    matrix_->ensureCols(new_cols);
  }
}

inline BondCutoffMatrix::RowProxy::operator std::span<BondCutoffRange>() const noexcept {
  return {matrix_->rowData(row_), matrix_->cols()};
}

inline BondCutoffMatrix::RowProxy::operator std::span<const BondCutoffRange>() const noexcept {
  return {matrix_->rowData(row_), matrix_->cols()};
}

inline BondCutoffRange *BondCutoffMatrix::RowProxy::begin() const noexcept {
  return matrix_->rowData(row_);
}

inline BondCutoffRange *BondCutoffMatrix::RowProxy::end() const noexcept {
  return matrix_->rowData(row_) + matrix_->cols();
}

inline const BondCutoffRange &
BondCutoffMatrix::ConstRowProxy::operator[](size_t col) const noexcept {
  return (*matrix_)(row_, col);
}

inline size_t BondCutoffMatrix::ConstRowProxy::size() const noexcept { return matrix_->cols(); }

inline bool BondCutoffMatrix::ConstRowProxy::empty() const noexcept { return matrix_->cols() == 0; }

inline BondCutoffMatrix::ConstRowProxy::operator std::span<const BondCutoffRange>() const noexcept {
  return {matrix_->rowData(row_), matrix_->cols()};
}

inline const BondCutoffRange *BondCutoffMatrix::ConstRowProxy::begin() const noexcept {
  return matrix_->rowData(row_);
}

inline const BondCutoffRange *BondCutoffMatrix::ConstRowProxy::end() const noexcept {
  return matrix_->rowData(row_) + matrix_->cols();
}

/**
 * @enum RingType
 * @brief Ring definition used for cycle topology search.
 */
enum class RingType : std::uint8_t {
  King,      ///< King ring: geodesic distance in G equals perimeter distance along C.
  Franzblau, ///< Franzblau shortest-path primitive ring: sub-paths <= floor(|C|/2) are unique
             ///< geodesics.
};

/**
 * @enum RingProjectionMode
 * @brief Projection or filtering strategy for network ring tracing.
 */
enum class RingProjectionMode : std::uint8_t {
  Direct,             ///< Direct search on raw atomic neighbor graph.
  BridgedProjection,  ///< Contract bridging atoms to project network-former polyhedra connectivity
                      ///< (e.g. Si-O-Si -> Si-Si).
  AlternatingTracing, ///< Trace bipartite cycles strictly alternating former and bridging elements
                      ///< (e.g. Si-O-Si-O).
};

} // namespace correlation::analysis

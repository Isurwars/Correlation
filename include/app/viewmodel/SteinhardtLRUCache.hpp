/**
 * @file SteinhardtLRUCache.hpp
 * @brief Thread-safe LRU cache for Steinhardt bond-order parameters across trajectory frames.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/order/SteinhardtCalculator.hpp"

#include <algorithm>
#include <cstddef>
#include <list>
#include <memory>
#include <shared_mutex>
#include <unordered_map>

namespace correlation::analysis {
class StructureAnalyzer;
} // namespace correlation::analysis

namespace correlation::app {

/**
 * @struct FramePhaseDescriptors
 * @brief Unified container holding computed per-atom phase order parameters for a frame.
 */
struct FramePhaseDescriptors {
  std::shared_ptr<const analysis::StructureAnalyzer> analyzer;
  std::shared_ptr<const calculators::SteinhardtCalculator::SteinhardtParams> steinhardt;
  std::vector<real_t> local_entropy;
  std::vector<real_t> twf_nconn;
  std::vector<real_t> q_tet;
  std::vector<real_t> displacement;
  std::vector<real_t> d2_min;
  std::vector<real_t> fermi_ldos;
};

namespace detail {

[[nodiscard]] inline std::size_t
estimateEntryBytes(const calculators::SteinhardtCalculator::SteinhardtParams &params) noexcept {
  const std::size_t elem_count = params.Q4.size() + params.Q6.size() + params.W4_hat.size() +
                                 params.W6_hat.size() + params.Q4_bar.size() + params.Q6_bar.size();
  std::size_t bytes = sizeof(params) + (elem_count * sizeof(real_t));
  if (params.q6m_normalized.has_value()) {
    bytes += params.q6m_normalized->size() * sizeof(std::complex<real_t>);
  }
  return bytes;
}

[[nodiscard]] inline std::size_t
estimateEntryBytes(const FramePhaseDescriptors &desc) noexcept {
  std::size_t bytes = sizeof(desc);
  if (desc.steinhardt != nullptr) {
    bytes += estimateEntryBytes(*desc.steinhardt);
  }
  const std::size_t scalar_count = desc.local_entropy.size() + desc.twf_nconn.size() +
                                   desc.q_tet.size() + desc.displacement.size() +
                                   desc.d2_min.size() + desc.fermi_ldos.size();
  bytes += scalar_count * sizeof(real_t);
  return bytes;
}

template <typename T> [[nodiscard]] inline std::size_t estimateEntryBytes(const T &val) noexcept {
  return sizeof(val);
}

} // namespace detail

/**
 * @class ThreadSafeLRUCache
 * @brief Thread-safe LRU cache storing per-atom parameters or descriptors for up to N trajectory
 * frames, bounded by both frame count and a configurable memory byte budget.
 *
 * Utilizes std::shared_mutex to permit concurrent read queries while synchronizing LRU mutations
 * and evictions under write locks.
 */
template <typename ValueType = calculators::SteinhardtCalculator::SteinhardtParams>
class ThreadSafeLRUCache {
public:
  using ParamsPtr = std::shared_ptr<const ValueType>;
  static constexpr std::size_t kDefaultMaxBytes = 256UL * 1024UL * 1024UL; ///< 256 MiB cap.

  /**
   * @brief Constructs cache with specified maximum frame capacity and memory byte budget.
   * @param[in] max_capacity Maximum number of frames cached simultaneously (default 512).
   * @param[in] max_bytes Maximum resident byte budget (default 256 MiB; 0 disables byte cap).
   */
  explicit ThreadSafeLRUCache(std::size_t max_capacity = 512,
                              std::size_t max_bytes = kDefaultMaxBytes) noexcept
      : capacity_(std::max<std::size_t>(1, max_capacity)), max_bytes_(max_bytes) {}

  ~ThreadSafeLRUCache() = default;
  ThreadSafeLRUCache(const ThreadSafeLRUCache &) = delete;
  ThreadSafeLRUCache &operator=(const ThreadSafeLRUCache &) = delete;
  ThreadSafeLRUCache(ThreadSafeLRUCache &&) noexcept = delete;
  ThreadSafeLRUCache &operator=(ThreadSafeLRUCache &&) noexcept = delete;

  /**
   * @brief Retrieves cached parameters and marks frame as most recently used.
   * @param[in] frame_index Trajectory frame index.
   * @return Shared pointer to cached parameters, or nullptr if not cached.
   */
  [[nodiscard]] ParamsPtr get(std::size_t frame_index) {
    const std::unique_lock lock(mutex_);
    const auto iter = map_.find(frame_index);
    if (iter == map_.end()) {
      return nullptr;
    }
    lru_order_.erase(iter->second.order_it);
    lru_order_.push_front(frame_index);
    iter->second.order_it = lru_order_.begin();
    return iter->second.params;
  }

  /**
   * @brief Inspects cache for frame without updating LRU order.
   * @param[in] frame_index Trajectory frame index.
   * @return Shared pointer to cached parameters, or nullptr if not cached.
   */
  [[nodiscard]] ParamsPtr peek(std::size_t frame_index) const {
    const std::shared_lock lock(mutex_);
    const auto iter = map_.find(frame_index);
    if (iter == map_.end()) {
      return nullptr;
    }
    return iter->second.params;
  }

  /**
   * @brief Inserts or updates cached parameters for a frame.
   * @param[in] frame_index Trajectory frame index.
   * @param[in] params Computed parameters.
   */
  void put(std::size_t frame_index, ParamsPtr params) {
    if (!params) {
      return;
    }
    const std::size_t entry_bytes = detail::estimateEntryBytes(*params);
    const std::unique_lock lock(mutex_);
    const auto iter = map_.find(frame_index);
    if (iter != map_.end()) {
      current_bytes_ = (current_bytes_ >= iter->second.bytes) ? (current_bytes_ - iter->second.bytes)
                                                              : 0;
      lru_order_.erase(iter->second.order_it);
      lru_order_.push_front(frame_index);
      iter->second.order_it = lru_order_.begin();
      iter->second.params = std::move(params);
      iter->second.bytes = entry_bytes;
      current_bytes_ += entry_bytes;
      evictIfNeededLocked();
      return;
    }

    lru_order_.push_front(frame_index);
    map_[frame_index] = CacheEntry{
        .params = std::move(params),
        .order_it = lru_order_.begin(),
        .bytes = entry_bytes,
    };
    current_bytes_ += entry_bytes;
    evictIfNeededLocked();
  }

  /**
   * @brief Checks if a given frame index resides in the cache.
   * @param[in] frame_index Trajectory frame index.
   * @return True if frame is cached.
   */
  [[nodiscard]] bool contains(std::size_t frame_index) const {
    const std::shared_lock lock(mutex_);
    return map_.contains(frame_index);
  }

  /**
   * @brief Gets current count of cached frames.
   * @return Number of entries in cache.
   */
  [[nodiscard]] std::size_t size() const {
    const std::shared_lock lock(mutex_);
    return map_.size();
  }

  /**
   * @brief Gets maximum frame capacity of the cache.
   * @return Cache capacity limit.
   */
  [[nodiscard]] std::size_t capacity() const noexcept {
    const std::shared_lock lock(mutex_);
    return capacity_;
  }

  /**
   * @brief Updates maximum frame capacity, evicting oldest entries if shrinking.
   * @param[in] new_capacity New maximum frame count (clamped to >= 1).
   */
  void setCapacity(std::size_t new_capacity) {
    const std::unique_lock lock(mutex_);
    capacity_ = std::max<std::size_t>(1, new_capacity);
    evictIfNeededLocked();
  }

  /**
   * @brief Gets current estimated memory footprint of cached entries in bytes.
   * @return Estimated resident bytes.
   */
  [[nodiscard]] std::size_t memoryBytes() const noexcept {
    const std::shared_lock lock(mutex_);
    return current_bytes_;
  }

  /**
   * @brief Gets configured maximum byte budget (0 indicates unlimited).
   * @return Maximum byte budget.
   */
  [[nodiscard]] std::size_t maxBytes() const noexcept {
    const std::shared_lock lock(mutex_);
    return max_bytes_;
  }

  /**
   * @brief Updates the byte budget, evicting oldest entries if over budget.
   * @param[in] new_max_bytes New byte budget (0 disables byte cap).
   */
  void setMaxBytes(std::size_t new_max_bytes) {
    const std::unique_lock lock(mutex_);
    max_bytes_ = new_max_bytes;
    evictIfNeededLocked();
  }

  /**
   * @brief Counts cached entries satisfying a predicate.
   * @tparam Pred Callable taking `const ValueType &` and returning bool.
   * @param[in] pred Predicate evaluated on each non-null cached value.
   * @return Number of matching cached entries.
   */
  template <typename Pred> [[nodiscard]] std::size_t countMatching(Pred &&pred) const {
    const std::shared_lock lock(mutex_);
    std::size_t count = 0;
    for (const auto &[idx, entry] : map_) {
      if (entry.params && pred(*entry.params)) {
        ++count;
      }
    }
    return count;
  }

  /**
   * @brief Purges all entries from the cache.
   */
  void clear() {
    const std::unique_lock lock(mutex_);
    map_.clear();
    lru_order_.clear();
    current_bytes_ = 0;
  }

private:
  struct CacheEntry {
    ParamsPtr params;
    std::list<std::size_t>::iterator order_it;
    std::size_t bytes{0};
  };

  void evictIfNeededLocked() {
    // Always keep at least 1 entry if present (the most recently inserted frame)
    while (map_.size() > 1 && !lru_order_.empty() &&
           (map_.size() > capacity_ || (max_bytes_ > 0 && current_bytes_ > max_bytes_))) {
      const auto oldest_idx = lru_order_.back();
      lru_order_.pop_back();
      if (const auto iter = map_.find(oldest_idx); iter != map_.end()) {
        current_bytes_ = (current_bytes_ >= iter->second.bytes)
                             ? (current_bytes_ - iter->second.bytes)
                             : 0;
        map_.erase(iter);
      }
    }
    while (map_.size() > capacity_ && !lru_order_.empty()) {
      const auto oldest_idx = lru_order_.back();
      lru_order_.pop_back();
      if (const auto iter = map_.find(oldest_idx); iter != map_.end()) {
        current_bytes_ = (current_bytes_ >= iter->second.bytes)
                             ? (current_bytes_ - iter->second.bytes)
                             : 0;
        map_.erase(iter);
      }
    }
  }

  mutable std::shared_mutex mutex_;
  std::size_t capacity_{512};
  std::size_t max_bytes_{kDefaultMaxBytes};
  std::size_t current_bytes_{0};
  std::list<std::size_t> lru_order_;
  std::unordered_map<std::size_t, CacheEntry> map_;
};

using SteinhardtLRUCache = ThreadSafeLRUCache<calculators::SteinhardtCalculator::SteinhardtParams>;
using PhaseDescriptorsLRUCache = ThreadSafeLRUCache<FramePhaseDescriptors>;

} // namespace correlation::app

/**
 * @file SteinhardtLRUCache.hpp
 * @brief Thread-safe LRU cache for Steinhardt bond-order parameters across trajectory frames.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/order/SteinhardtCalculator.hpp"

#include <cstddef>
#include <list>
#include <memory>
#include <shared_mutex>
#include <unordered_map>

namespace correlation::app {

/**
 * @struct FramePhaseDescriptors
 * @brief Unified container holding computed per-atom phase order parameters for a frame.
 */
struct FramePhaseDescriptors {
  std::shared_ptr<const calculators::SteinhardtCalculator::SteinhardtParams> steinhardt;
  std::vector<real_t> local_entropy;
};

/**
 * @class ThreadSafeLRUCache
 * @brief Thread-safe LRU cache storing per-atom parameters or descriptors for up to N trajectory
 * frames.
 *
 * Utilizes std::shared_mutex to permit concurrent read queries while synchronizing LRU mutations
 * and evictions under write locks.
 */
template <typename ValueType = calculators::SteinhardtCalculator::SteinhardtParams>
class ThreadSafeLRUCache {
public:
  using ParamsPtr = std::shared_ptr<const ValueType>;

  /**
   * @brief Constructs cache with specified maximum capacity.
   * @param[in] max_capacity Maximum number of frames cached simultaneously (default 32).
   */
  explicit ThreadSafeLRUCache(std::size_t max_capacity = 32) noexcept : capacity_(max_capacity) {}

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
    const std::unique_lock lock(mutex_);
    const auto iter = map_.find(frame_index);
    if (iter != map_.end()) {
      lru_order_.erase(iter->second.order_it);
      lru_order_.push_front(frame_index);
      iter->second.order_it = lru_order_.begin();
      iter->second.params = std::move(params);
      return;
    }

    if (map_.size() >= capacity_ && !lru_order_.empty()) {
      const auto oldest_idx = lru_order_.back();
      lru_order_.pop_back();
      map_.erase(oldest_idx);
    }

    lru_order_.push_front(frame_index);
    map_[frame_index] = CacheEntry{
        .params = std::move(params),
        .order_it = lru_order_.begin(),
    };
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
  [[nodiscard]] std::size_t capacity() const noexcept { return capacity_; }

  /**
   * @brief Purges all entries from the cache.
   */
  void clear() {
    const std::unique_lock lock(mutex_);
    map_.clear();
    lru_order_.clear();
  }

private:
  struct CacheEntry {
    ParamsPtr params;
    std::list<std::size_t>::iterator order_it;
  };

  mutable std::shared_mutex mutex_;
  std::size_t capacity_{32};
  std::list<std::size_t> lru_order_;
  std::unordered_map<std::size_t, CacheEntry> map_;
};

using SteinhardtLRUCache = ThreadSafeLRUCache<calculators::SteinhardtCalculator::SteinhardtParams>;
using PhaseDescriptorsLRUCache = ThreadSafeLRUCache<FramePhaseDescriptors>;

} // namespace correlation::app

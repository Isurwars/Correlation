/**
 * @file OrbV3GNNModel.cpp
 * @brief Native C++23 implementation of ORB-v3 GNN + LDOSHead inference engine.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "mlip/OrbV3GNNModel.hpp"
#include "core/MappedFile.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string_view>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <vector>

namespace correlation::mlip {

namespace {

constexpr std::string_view K_EXPECTED_MAGIC = "CORRORB3";
constexpr std::uint32_t K_EXPECTED_VERSION = 1;
constexpr size_t K_ATOM_EMB_ROWS = 119;
constexpr size_t K_REF_ENERGY_COUNT = 118;
constexpr size_t K_LDOS_RES_BLOCKS = 3;
constexpr float K_RMSNORM_EPS = 1.1920929e-7F;
constexpr float K_LAYERNORM_EPS = 1.0e-5F;

struct MLP3RMSNormWeights {
  std::span<const float> w0;
  std::span<const float> b0;
  std::span<const float> w1;
  std::span<const float> b1;
  std::span<const float> w2;
  std::span<const float> b2;
  std::span<const float> rms_w;
  size_t in_dim{0};
  size_t hidden_dim{0};
  size_t out_dim{0};
};

struct AttentionStackWeights {
  MLP3RMSNormWeights edge_mlp;
  MLP3RMSNormWeights node_mlp;
  std::span<const float> recv_attn_w;
  std::span<const float> recv_attn_b;
  std::span<const float> send_attn_w;
  std::span<const float> send_attn_b;
};

struct EnergyHeadWeights {
  float norm_mean{0.0F};
  float norm_std{1.0F};
  std::span<const float> w0;
  std::span<const float> b0;
  std::span<const float> w1;
  std::span<const float> b1;
  std::span<const float> ref_energies;
};

struct ForcesHeadWeights {
  float norm_mean{0.0F};
  float norm_std{1.0F};
  std::span<const float> w0;
  std::span<const float> b0;
  std::span<const float> w1;
  std::span<const float> b1;
};

struct StressHeadWeights {
  float diag_mean{0.0F};
  float diag_std{1.0F};
  float off_mean{0.0F};
  float off_std{1.0F};
  std::span<const float> w0;
  std::span<const float> b0;
  std::span<const float> w1;
  std::span<const float> b1;
};

struct LDOSResBlockWeights {
  std::span<const float> ln_w;
  std::span<const float> ln_b;
  std::span<const float> fc_w;
  std::span<const float> fc_b;
};

struct LDOSHeadWeights {
  std::span<const float> node_norm_w;
  std::span<const float> node_norm_b;
  std::span<const float> fc0_w;
  std::span<const float> fc0_b;
  std::array<LDOSResBlockWeights, K_LDOS_RES_BLOCKS> blocks{};
  std::span<const float> out_norm_w;
  std::span<const float> out_norm_b;
  std::span<const float> fc_out_w;
  std::span<const float> fc_out_b;
};

class WeightCursor {
public:
  WeightCursor(const float *base, size_t total_floats) noexcept
      : base_(base), total_(total_floats) {}

  [[nodiscard]] std::span<const float> take(size_t count) {
    if (offset_ + count > total_) {
      throw std::runtime_error("OrbV3GNNModel: weight file truncated (buffer underrun).");
    }
    const std::span<const float> slice(base_ + offset_, count);
    offset_ += count;
    return slice;
  }

  [[nodiscard]] size_t remaining() const noexcept { return total_ - offset_; }

private:
  const float *base_{nullptr};
  size_t total_{0};
  size_t offset_{0};
};

[[nodiscard]] MLP3RMSNormWeights takeMLP3RMSNorm(WeightCursor &cursor, size_t in_dim,
                                                 size_t hidden_dim, size_t out_dim) {
  MLP3RMSNormWeights mlp;
  mlp.in_dim = in_dim;
  mlp.hidden_dim = hidden_dim;
  mlp.out_dim = out_dim;
  mlp.w0 = cursor.take(hidden_dim * in_dim);
  mlp.b0 = cursor.take(hidden_dim);
  mlp.w1 = cursor.take(hidden_dim * hidden_dim);
  mlp.b1 = cursor.take(hidden_dim);
  mlp.w2 = cursor.take(out_dim * hidden_dim);
  mlp.b2 = cursor.take(out_dim);
  mlp.rms_w = cursor.take(out_dim);
  return mlp;
}

[[nodiscard]] inline float silu(float val) noexcept { return val / (1.0F + std::exp(-val)); }

[[nodiscard]] inline float sigmoid(float val) noexcept { return 1.0F / (1.0F + std::exp(-val)); }

[[nodiscard]] inline float softplusBeta2(float val) noexcept {
  const float scaled = 2.0F * val;
  if (scaled > 20.0F) {
    return val;
  }
  if (scaled < -20.0F) {
    return 0.5F * std::exp(scaled);
  }
  return 0.5F * std::log1p(std::exp(scaled));
}

[[nodiscard]] inline float dotProductF32(const float *lhs, const float *rhs,
                                         size_t count) noexcept {
  float acc = 0.0F;
  for (size_t k = 0; k < count; ++k) {
    acc += lhs[k] * rhs[k];
  }
  return acc;
}

void computeLinearRow(const float *in_row, size_t in_dim, const float *weight, const float *bias,
                      size_t out_dim, bool apply_silu, float *out_row) noexcept {
  for (size_t col = 0; col < out_dim; ++col) {
    const float sum = dotProductF32(in_row, weight + col * in_dim, in_dim) + bias[col];
    out_row[col] = apply_silu ? silu(sum) : sum;
  }
}

void linearForward(std::span<const float> input, size_t num_rows, size_t in_dim,
                   std::span<const float> weight, std::span<const float> bias, size_t out_dim,
                   bool apply_silu, std::span<float> output) {
  if (num_rows == 0) {
    return;
  }
  const float *in_ptr = input.data();
  const float *w_ptr = weight.data();
  const float *b_ptr = bias.data();
  float *out_ptr = output.data();

  if (num_rows * out_dim >= 2048) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_rows),
                      [&](const tbb::blocked_range<size_t> &range) {
                        for (size_t row = range.begin(); row != range.end(); ++row) {
                          computeLinearRow(in_ptr + row * in_dim, in_dim, w_ptr, b_ptr, out_dim,
                                           apply_silu, out_ptr + row * out_dim);
                        }
                      });
  } else {
    for (size_t row = 0; row < num_rows; ++row) {
      computeLinearRow(in_ptr + row * in_dim, in_dim, w_ptr, b_ptr, out_dim, apply_silu,
                       out_ptr + row * out_dim);
    }
  }
}

void applyRMSNormInPlace(std::span<float> data, size_t num_rows, size_t dim,
                         std::span<const float> gamma) noexcept {
  if (num_rows == 0 || dim == 0) {
    return;
  }
  const float inv_dim = 1.0F / static_cast<float>(dim);
  const float *g_ptr = gamma.data();

  for (size_t row = 0; row < num_rows; ++row) {
    float *row_ptr = data.data() + row * dim;
    const float sum_sq = dotProductF32(row_ptr, row_ptr, dim);
    const float inv_rms = 1.0F / std::sqrt(sum_sq * inv_dim + K_RMSNORM_EPS);
    for (size_t col = 0; col < dim; ++col) {
      row_ptr[col] = row_ptr[col] * inv_rms * g_ptr[col];
    }
  }
}

void applyLayerNorm(std::span<const float> input, size_t num_rows, size_t dim,
                    std::span<const float> gamma, std::span<const float> beta,
                    std::span<float> output) noexcept {
  if (num_rows == 0 || dim == 0) {
    return;
  }
  const float inv_dim = 1.0F / static_cast<float>(dim);
  const float *g_ptr = gamma.data();
  const float *b_ptr = beta.data();

  for (size_t row = 0; row < num_rows; ++row) {
    const float *in_row = input.data() + row * dim;
    float *out_row = output.data() + row * dim;

    float sum = 0.0F;
    for (size_t col = 0; col < dim; ++col) {
      sum += in_row[col];
    }
    const float mean = sum * inv_dim;

    float var_sum = 0.0F;
    for (size_t col = 0; col < dim; ++col) {
      const float diff = in_row[col] - mean;
      var_sum += diff * diff;
    }
    const float inv_std = 1.0F / std::sqrt(var_sum * inv_dim + K_LAYERNORM_EPS);

    for (size_t col = 0; col < dim; ++col) {
      out_row[col] = (in_row[col] - mean) * inv_std * g_ptr[col] + b_ptr[col];
    }
  }
}

void evaluateMLP3RMSNorm(std::span<const float> input, size_t num_rows,
                         const MLP3RMSNormWeights &mlp, std::vector<float> &scratch_a,
                         std::vector<float> &scratch_b, std::span<float> output) {
  if (num_rows == 0) {
    return;
  }
  scratch_a.resize(num_rows * mlp.hidden_dim);
  scratch_b.resize(num_rows * mlp.hidden_dim);

  linearForward(input, num_rows, mlp.in_dim, mlp.w0, mlp.b0, mlp.hidden_dim, true, scratch_a);
  linearForward(scratch_a, num_rows, mlp.hidden_dim, mlp.w1, mlp.b1, mlp.hidden_dim, true,
                scratch_b);
  linearForward(scratch_b, num_rows, mlp.hidden_dim, mlp.w2, mlp.b2, mlp.out_dim, false, output);
  applyRMSNormInPlace(output, num_rows, mlp.out_dim, mlp.rms_w);
}

void addResidualInPlace(std::span<float> dst, std::span<const float> delta) noexcept {
  const size_t total = dst.size();
  float *d_ptr = dst.data();
  const float *s_ptr = delta.data();
  for (size_t i = 0; i < total; ++i) {
    d_ptr[i] += s_ptr[i];
  }
}

} // namespace

struct OrbV3GNNModel::Impl {
  OrbV3ModelHeader header{};
  std::unique_ptr<core::MappedFile> mapped_file;
  std::span<const float> atom_emb;
  MLP3RMSNormWeights encoder_node;
  MLP3RMSNormWeights encoder_edge;
  std::vector<AttentionStackWeights> stacks;
  EnergyHeadWeights energy_head;
  ForcesHeadWeights forces_head;
  StressHeadWeights stress_head;
  LDOSHeadWeights ldos_head;
  bool loaded{false};

  void parseMappedWeights() {
    if (!mapped_file || mapped_file->size() < sizeof(OrbV3ModelHeader)) {
      throw std::runtime_error("OrbV3GNNModel: file is smaller than 64-byte header.");
    }
    std::memcpy(&header, mapped_file->data(), sizeof(OrbV3ModelHeader));

    const std::string_view magic_view(header.magic.data(), header.magic.size());
    if (magic_view != K_EXPECTED_MAGIC || header.version != K_EXPECTED_VERSION) {
      throw std::runtime_error("OrbV3GNNModel: invalid magic or version in weight file.");
    }

    const size_t payload_bytes = mapped_file->size() - sizeof(OrbV3ModelHeader);
    if (payload_bytes % sizeof(float) != 0) {
      throw std::runtime_error("OrbV3GNNModel: payload byte size is not a multiple of 4.");
    }

    const auto *float_base =
        reinterpret_cast<const float *>(mapped_file->data() + sizeof(OrbV3ModelHeader));
    WeightCursor cursor(float_base, payload_bytes / sizeof(float));

    const size_t d_lat = header.latent_dim;
    const size_t d_mlp = header.mlp_hidden_dim;
    const size_t d_head = header.head_hidden_dim;
    const size_t d_ldos = header.ldos_hidden_dim;
    const size_t n_sh = static_cast<size_t>(header.l_max + 1) * (header.l_max + 1);
    const size_t d_edge_in = static_cast<size_t>(header.num_rbf) * n_sh;

    atom_emb = cursor.take(K_ATOM_EMB_ROWS * d_lat);
    encoder_node = takeMLP3RMSNorm(cursor, d_lat, d_mlp, d_lat);
    encoder_edge = takeMLP3RMSNorm(cursor, d_edge_in, d_mlp, d_lat);

    stacks.resize(header.num_stacks);
    for (size_t idx = 0; idx < header.num_stacks; ++idx) {
      stacks[idx].edge_mlp = takeMLP3RMSNorm(cursor, d_lat * 3, d_mlp, d_lat);
      stacks[idx].node_mlp = takeMLP3RMSNorm(cursor, d_lat * 3, d_mlp, d_lat);
      stacks[idx].recv_attn_w = cursor.take(d_lat);
      stacks[idx].recv_attn_b = cursor.take(1);
      stacks[idx].send_attn_w = cursor.take(d_lat);
      stacks[idx].send_attn_b = cursor.take(1);
    }

    const auto e_norm = cursor.take(2);
    energy_head.norm_mean = e_norm[0];
    energy_head.norm_std = e_norm[1];
    energy_head.w0 = cursor.take(d_head * d_lat);
    energy_head.b0 = cursor.take(d_head);
    energy_head.w1 = cursor.take(d_head);
    energy_head.b1 = cursor.take(1);
    energy_head.ref_energies = cursor.take(K_REF_ENERGY_COUNT);

    const auto f_norm = cursor.take(2);
    forces_head.norm_mean = f_norm[0];
    forces_head.norm_std = f_norm[1];
    forces_head.w0 = cursor.take(d_head * d_lat);
    forces_head.b0 = cursor.take(d_head);
    forces_head.w1 = cursor.take(3 * d_head);
    forces_head.b1 = cursor.take(3);

    const auto s_norm = cursor.take(4);
    stress_head.diag_mean = s_norm[0];
    stress_head.diag_std = s_norm[1];
    stress_head.off_mean = s_norm[2];
    stress_head.off_std = s_norm[3];
    stress_head.w0 = cursor.take(d_head * d_lat);
    stress_head.b0 = cursor.take(d_head);
    stress_head.w1 = cursor.take(6 * d_head);
    stress_head.b1 = cursor.take(6);

    ldos_head.node_norm_w = cursor.take(d_lat);
    ldos_head.node_norm_b = cursor.take(d_lat);
    ldos_head.fc0_w = cursor.take(d_ldos * d_lat);
    ldos_head.fc0_b = cursor.take(d_ldos);
    for (auto &block : ldos_head.blocks) {
      block.ln_w = cursor.take(d_ldos);
      block.ln_b = cursor.take(d_ldos);
      block.fc_w = cursor.take(d_ldos * d_ldos);
      block.fc_b = cursor.take(d_ldos);
    }
    ldos_head.out_norm_w = cursor.take(d_ldos);
    ldos_head.out_norm_b = cursor.take(d_ldos);
    ldos_head.fc_out_w = cursor.take(static_cast<size_t>(header.ldos_bins) * d_ldos);
    ldos_head.fc_out_b = cursor.take(header.ldos_bins);

    if (cursor.remaining() != 0) {
      throw std::runtime_error("OrbV3GNNModel: unexpected trailing bytes in weight file.");
    }
    loaded = true;
  }

  void runInteractionStack(const AttentionStackWeights &stack, size_t num_atoms, size_t num_edges,
                           std::span<const int64_t> senders, std::span<const int64_t> receivers,
                           std::span<const real_t> cutoffs, std::vector<float> &nodes,
                           std::vector<float> &edges, std::vector<float> &scratch_a,
                           std::vector<float> &scratch_b) const {
    const size_t d_lat = header.latent_dim;
    std::vector<float> recv_gates(num_edges, 0.0F);
    std::vector<float> send_gates(num_edges, 0.0F);
    std::vector<float> edge_concat(num_edges * d_lat * 3, 0.0F);

    for (size_t idx_e = 0; idx_e < num_edges; ++idx_e) {
      const float *edge_vec = edges.data() + idx_e * d_lat;
      const auto cutoff_f = static_cast<float>(cutoffs[idx_e]);
      const float r_logit =
          dotProductF32(edge_vec, stack.recv_attn_w.data(), d_lat) + stack.recv_attn_b[0];
      const float s_logit =
          dotProductF32(edge_vec, stack.send_attn_w.data(), d_lat) + stack.send_attn_b[0];
      recv_gates[idx_e] = sigmoid(r_logit) * cutoff_f;
      send_gates[idx_e] = sigmoid(s_logit) * cutoff_f;

      const auto src_idx = static_cast<size_t>(senders[idx_e]);
      const auto dst_idx = static_cast<size_t>(receivers[idx_e]);
      float *dst_concat = edge_concat.data() + idx_e * d_lat * 3;
      std::memcpy(dst_concat, edge_vec, d_lat * sizeof(float));
      std::memcpy(dst_concat + d_lat, nodes.data() + src_idx * d_lat, d_lat * sizeof(float));
      std::memcpy(dst_concat + 2 * d_lat, nodes.data() + dst_idx * d_lat, d_lat * sizeof(float));
    }

    std::vector<float> updated_edges(num_edges * d_lat, 0.0F);
    evaluateMLP3RMSNorm(edge_concat, num_edges, stack.edge_mlp, scratch_a, scratch_b,
                        updated_edges);

    std::vector<float> sent_msg(num_atoms * d_lat, 0.0F);
    std::vector<float> recv_msg(num_atoms * d_lat, 0.0F);
    for (size_t idx_e = 0; idx_e < num_edges; ++idx_e) {
      const auto src_idx = static_cast<size_t>(senders[idx_e]);
      const auto dst_idx = static_cast<size_t>(receivers[idx_e]);
      const float s_gate = send_gates[idx_e];
      const float r_gate = recv_gates[idx_e];
      const float *u_edge = updated_edges.data() + idx_e * d_lat;
      float *s_dst = sent_msg.data() + src_idx * d_lat;
      float *r_dst = recv_msg.data() + dst_idx * d_lat;
      for (size_t k = 0; k < d_lat; ++k) {
        s_dst[k] += u_edge[k] * s_gate;
        r_dst[k] += u_edge[k] * r_gate;
      }
    }

    std::vector<float> node_concat(num_atoms * d_lat * 3, 0.0F);
    for (size_t idx_n = 0; idx_n < num_atoms; ++idx_n) {
      float *dst_concat = node_concat.data() + idx_n * d_lat * 3;
      std::memcpy(dst_concat, nodes.data() + idx_n * d_lat, d_lat * sizeof(float));
      std::memcpy(dst_concat + d_lat, recv_msg.data() + idx_n * d_lat, d_lat * sizeof(float));
      std::memcpy(dst_concat + 2 * d_lat, sent_msg.data() + idx_n * d_lat, d_lat * sizeof(float));
    }

    std::vector<float> updated_nodes(num_atoms * d_lat, 0.0F);
    evaluateMLP3RMSNorm(node_concat, num_atoms, stack.node_mlp, scratch_a, scratch_b,
                        updated_nodes);

    addResidualInPlace(nodes, updated_nodes);
    addResidualInPlace(edges, updated_edges);
  }

  void evaluateLDOSHead(std::span<const float> nodes, size_t num_atoms, MLIPOutput &output) const {
    const size_t d_lat = header.latent_dim;
    const size_t d_ldos = header.ldos_hidden_dim;
    const size_t n_bins = header.ldos_bins;

    std::vector<float> normed_nodes(num_atoms * d_lat);
    applyLayerNorm(nodes, num_atoms, d_lat, ldos_head.node_norm_w, ldos_head.node_norm_b,
                   normed_nodes);

    std::vector<float> hidden(num_atoms * d_ldos);
    linearForward(normed_nodes, num_atoms, d_lat, ldos_head.fc0_w, ldos_head.fc0_b, d_ldos, true,
                  hidden);

    std::vector<float> block_normed(num_atoms * d_ldos);
    std::vector<float> block_out(num_atoms * d_ldos);
    for (const auto &block : ldos_head.blocks) {
      applyLayerNorm(hidden, num_atoms, d_ldos, block.ln_w, block.ln_b, block_normed);
      linearForward(block_normed, num_atoms, d_ldos, block.fc_w, block.fc_b, d_ldos, true,
                    block_out);
      addResidualInPlace(hidden, block_out);
    }

    applyLayerNorm(hidden, num_atoms, d_ldos, ldos_head.out_norm_w, ldos_head.out_norm_b,
                   block_normed);

    std::vector<float> logits(num_atoms * n_bins);
    linearForward(block_normed, num_atoms, d_ldos, ldos_head.fc_out_w, ldos_head.fc_out_b, n_bins,
                  false, logits);

    output.ldos_bins = n_bins;
    output.ldos_e_min = static_cast<real_t>(header.ldos_emin);
    output.ldos_e_max = static_cast<real_t>(header.ldos_emax);
    output.ldos_sigma = static_cast<real_t>(header.ldos_sigma);
    output.ldos.resize(num_atoms * n_bins);

    for (size_t idx = 0; idx < logits.size(); ++idx) {
      output.ldos[idx] = static_cast<real_t>(softplusBeta2(logits[idx]));
    }
  }

  void evaluatePhysicalHeads(std::span<const float> nodes, std::span<const int64_t> atomic_numbers,
                             size_t num_atoms, MLIPOutput &output) const {
    const size_t d_lat = header.latent_dim;
    const size_t d_head = header.head_hidden_dim;
    const float inv_n = 1.0F / static_cast<float>(num_atoms);

    // Mean-pooled node representation [1, d_lat]
    std::vector<float> mean_node(d_lat, 0.0F);
    for (size_t i = 0; i < num_atoms; ++i) {
      const float *n_ptr = nodes.data() + i * d_lat;
      for (size_t k = 0; k < d_lat; ++k) {
        mean_node[k] += n_ptr[k] * inv_n;
      }
    }

    // 1. Energy Head (system total + per-atom decomposition)
    std::vector<float> e_hid(d_head);
    float e_raw = 0.0F;
    linearForward(mean_node, 1, d_lat, energy_head.w0, energy_head.b0, d_head, true, e_hid);
    linearForward(e_hid, 1, d_head, energy_head.w1, energy_head.b1, 1, false,
                  std::span<float>(&e_raw, 1));

    const double interaction_e =
        static_cast<double>(e_raw * energy_head.norm_std + energy_head.norm_mean) *
        static_cast<double>(num_atoms);
    double ref_sum = 0.0;
    for (size_t i = 0; i < num_atoms; ++i) {
      const auto z_val = static_cast<size_t>(atomic_numbers[i]);
      if (z_val < K_REF_ENERGY_COUNT) {
        ref_sum += static_cast<double>(energy_head.ref_energies[z_val]);
      }
    }
    output.total_energy = static_cast<real_t>(interaction_e + ref_sum);

    // Site-resolved per-atom energy evaluation
    std::vector<float> atom_e_hid(num_atoms * d_head);
    std::vector<float> atom_e_raw(num_atoms);
    linearForward(nodes, num_atoms, d_lat, energy_head.w0, energy_head.b0, d_head, true,
                  atom_e_hid);
    linearForward(atom_e_hid, num_atoms, d_head, energy_head.w1, energy_head.b1, 1, false,
                  atom_e_raw);
    output.per_atom_energy.resize(num_atoms);
    for (size_t i = 0; i < num_atoms; ++i) {
      const auto z_val = static_cast<size_t>(atomic_numbers[i]);
      const float ref_i = (z_val < K_REF_ENERGY_COUNT) ? energy_head.ref_energies[z_val] : 0.0F;
      const float e_i = (atom_e_raw[i] * energy_head.norm_std + energy_head.norm_mean) + ref_i;
      output.per_atom_energy[i] = static_cast<real_t>(e_i);
    }

    // 2. Forces Head (with net-force zero-mean removal)
    std::vector<float> f_hid(num_atoms * d_head);
    std::vector<float> f_raw(num_atoms * 3);
    linearForward(nodes, num_atoms, d_lat, forces_head.w0, forces_head.b0, d_head, true, f_hid);
    linearForward(f_hid, num_atoms, d_head, forces_head.w1, forces_head.b1, 3, false, f_raw);

    double f_mean_x = 0.0;
    double f_mean_y = 0.0;
    double f_mean_z = 0.0;
    for (size_t i = 0; i < num_atoms; ++i) {
      const float fx_denorm = f_raw[i * 3 + 0] * forces_head.norm_std + forces_head.norm_mean;
      const float fy_denorm = f_raw[i * 3 + 1] * forces_head.norm_std + forces_head.norm_mean;
      const float fz_denorm = f_raw[i * 3 + 2] * forces_head.norm_std + forces_head.norm_mean;
      f_raw[i * 3 + 0] = fx_denorm;
      f_raw[i * 3 + 1] = fy_denorm;
      f_raw[i * 3 + 2] = fz_denorm;
      f_mean_x += static_cast<double>(fx_denorm);
      f_mean_y += static_cast<double>(fy_denorm);
      f_mean_z += static_cast<double>(fz_denorm);
    }
    const auto inv_atoms_d = 1.0 / static_cast<double>(num_atoms);
    f_mean_x *= inv_atoms_d;
    f_mean_y *= inv_atoms_d;
    f_mean_z *= inv_atoms_d;

    output.forces.resize(num_atoms);
    for (size_t i = 0; i < num_atoms; ++i) {
      output.forces[i] = correlation::math::Vector3<real_t>{
          static_cast<real_t>(static_cast<double>(f_raw[i * 3 + 0]) - f_mean_x),
          static_cast<real_t>(static_cast<double>(f_raw[i * 3 + 1]) - f_mean_y),
          static_cast<real_t>(static_cast<double>(f_raw[i * 3 + 2]) - f_mean_z),
      };
    }

    // 3. Stress Head (Voigt [xx, yy, zz, yz, xz, xy] -> 3x3 symmetric tensor)
    std::vector<float> s_hid(d_head);
    std::array<float, 6> s_raw{};
    linearForward(mean_node, 1, d_lat, stress_head.w0, stress_head.b0, d_head, true, s_hid);
    linearForward(s_hid, 1, d_head, stress_head.w1, stress_head.b1, 6, false, s_raw);

    const auto s_xx = static_cast<real_t>(s_raw[0] * stress_head.diag_std + stress_head.diag_mean);
    const auto s_yy = static_cast<real_t>(s_raw[1] * stress_head.diag_std + stress_head.diag_mean);
    const auto s_zz = static_cast<real_t>(s_raw[2] * stress_head.diag_std + stress_head.diag_mean);
    const auto s_yz = static_cast<real_t>(s_raw[3] * stress_head.off_std + stress_head.off_mean);
    const auto s_xz = static_cast<real_t>(s_raw[4] * stress_head.off_std + stress_head.off_mean);
    const auto s_xy = static_cast<real_t>(s_raw[5] * stress_head.off_std + stress_head.off_mean);

    output.stress =
        correlation::math::Matrix3<real_t>(correlation::math::Vector3<real_t>{s_xx, s_xy, s_xz},
                                           correlation::math::Vector3<real_t>{s_xy, s_yy, s_yz},
                                           correlation::math::Vector3<real_t>{s_xz, s_yz, s_zz});
  }
};

std::string OrbV3GNNModel::findDefaultModelPath(const std::string &preferred_path) {
  if (!preferred_path.empty() && std::filesystem::is_regular_file(preferred_path)) {
    return preferred_path;
  }

  std::vector<std::filesystem::path> candidates = {
      "models/cellulose_ldos.bin",
      "../models/cellulose_ldos.bin",
      "../../models/cellulose_ldos.bin",
      "/home/isurwars/Projects/Cellulose/exported_models/cellulose_ldos.bin",
  };

  if (const char *home_env = std::getenv("HOME"); home_env != nullptr) {
    const std::filesystem::path home_dir(home_env);
    candidates.push_back(home_dir / ".local/share/correlation/models/cellulose_ldos.bin");
    candidates.push_back(home_dir / "Projects/Correlation/models/cellulose_ldos.bin");
    candidates.push_back(home_dir / "Projects/Cellulose/exported_models/cellulose_ldos.bin");
  }

  for (const auto &candidate : candidates) {
    std::error_code err_code;
    if (std::filesystem::is_regular_file(candidate, err_code)) {
      return candidate.string();
    }
  }
  return {};
}

OrbV3GNNModel::OrbV3GNNModel(std::string model_path)
    : impl_(std::make_unique<Impl>()), model_path_(findDefaultModelPath(model_path)) {
  if (model_path_.empty() && !model_path.empty()) {
    model_path_ = std::move(model_path);
  }
  if (model_path_.empty()) {
    return;
  }
  try {
    impl_->mapped_file = std::make_unique<core::MappedFile>(model_path_);
    impl_->parseMappedWeights();
  } catch (const std::exception &e) {
    std::cerr << "[OrbV3GNNModel] Failed to load model '" << model_path_ << "': " << e.what()
              << '\n';
    impl_->loaded = false;
  }
}

OrbV3GNNModel::~OrbV3GNNModel() = default;
OrbV3GNNModel::OrbV3GNNModel(OrbV3GNNModel &&) noexcept = default;
OrbV3GNNModel &OrbV3GNNModel::operator=(OrbV3GNNModel &&) noexcept = default;

std::string OrbV3GNNModel::getModelName() const { return "ORB-v3-LDOS [" + model_path_ + "]"; }

bool OrbV3GNNModel::isLoaded() const noexcept { return impl_ && impl_->loaded; }

const OrbV3ModelHeader &OrbV3GNNModel::header() const noexcept { return impl_->header; }

MLIPOutput OrbV3GNNModel::evaluate(const correlation::core::Cell &cell) const {
  MLIPOutput output;
  const size_t num_atoms = cell.atoms().size();
  if (num_atoms == 0) {
    return output;
  }

  if (!isLoaded()) {
    throw std::runtime_error("OrbV3GNNModel: model weights are not loaded ('" + model_path_ +
                             "').");
  }

  const auto compat = PeriodicGraphBuilder::validateElementCompatibility(cell);
  if (compat.has_unmapped_elements) {
    std::string symbols_str;
    for (const auto &sym : compat.unmapped_symbols) {
      if (!symbols_str.empty()) {
        symbols_str += ", ";
      }
      symbols_str += "'" + sym + "'";
    }
    throw std::runtime_error(
        "OrbV3GNNModel: unmapped atomic species (" + symbols_str +
        ") detected (Z=0). Remap numeric LAMMPS atom types to chemical elements before running "
        "eDoS/MLIP.");
  }

  if (compat.has_out_of_domain_elements) {
    std::cerr << "[OrbV3GNNModel] Warning: cell contains elements outside Cellulose fine-tuning "
                 "set {H, C, O}; evaluating with pretrained OMat24 ORB-v3 backbone embeddings.\n";
  }

  const OrbDescriptorConfig graph_cfg{
      .r_max = static_cast<real_t>(impl_->header.r_max),
      .num_rbf = impl_->header.num_rbf,
      .l_max = impl_->header.l_max,
      .max_num_neighbors = impl_->header.max_neighbors,
      .include_self_loops = false,
      .compute_orb_features = true,
  };
  const auto graph = PeriodicGraphBuilder::buildOrbGraph(cell, graph_cfg);
  const size_t num_edges = graph.edge_count;
  const size_t d_lat = impl_->header.latent_dim;
  const size_t d_edge_in = impl_->encoder_edge.in_dim;

  // 1. Lookup initial atom embeddings [N, d_lat]
  std::vector<float> raw_nodes(num_atoms * d_lat, 0.0F);
  for (size_t i = 0; i < num_atoms; ++i) {
    const auto z_val = static_cast<size_t>(std::clamp<int64_t>(graph.atomic_numbers[i], 0, 118));
    std::memcpy(raw_nodes.data() + i * d_lat, impl_->atom_emb.data() + z_val * d_lat,
                d_lat * sizeof(float));
  }

  // Convert edge_orb_features_flat to float32 if real_t is double
  std::vector<float> raw_edges(num_edges * d_edge_in, 0.0F);
  for (size_t idx = 0; idx < raw_edges.size(); ++idx) {
    raw_edges[idx] = static_cast<float>(graph.edge_orb_features_flat[idx]);
  }

  // 2. Encode nodes and edges
  std::vector<float> nodes(num_atoms * d_lat, 0.0F);
  std::vector<float> edges(num_edges * d_lat, 0.0F);
  std::vector<float> scratch_a;
  std::vector<float> scratch_b;

  evaluateMLP3RMSNorm(raw_nodes, num_atoms, impl_->encoder_node, scratch_a, scratch_b, nodes);
  evaluateMLP3RMSNorm(raw_edges, num_edges, impl_->encoder_edge, scratch_a, scratch_b, edges);

  // 3. Message-passing interaction stacks
  const std::span<const int64_t> senders(graph.edge_index_flat.data(), num_edges);
  const std::span<const int64_t> receivers(graph.edge_index_flat.data() + num_edges, num_edges);
  for (const auto &stack : impl_->stacks) {
    impl_->runInteractionStack(stack, num_atoms, num_edges, senders, receivers,
                               graph.edge_cutoff_envelope_flat, nodes, edges, scratch_a, scratch_b);
  }

  // 4. Evaluate LDOS and physical heads
  impl_->evaluateLDOSHead(nodes, num_atoms, output);
  impl_->evaluatePhysicalHeads(nodes, graph.atomic_numbers, num_atoms, output);

  return output;
}

} // namespace correlation::mlip

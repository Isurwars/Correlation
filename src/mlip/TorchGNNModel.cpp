/**
 * @file TorchGNNModel.cpp
 * @brief Implementation of TorchGNNModel for LibTorch GNN evaluation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "mlip/TorchGNNModel.hpp"
#include <iostream>

#ifdef CORRELATION_HAS_LIBTORCH
#include <torch/script.h>
#include <torch/torch.h>
#endif

namespace correlation::mlip {

struct TorchGNNModel::Impl {
#ifdef CORRELATION_HAS_LIBTORCH
  torch::jit::script::Module module;
  torch::Device device{torch::kCPU};
  bool loaded{false};
#else
  bool loaded{false};
#endif
};

TorchGNNModel::TorchGNNModel(std::string model_path, std::string device, real_t cutoff)
    : impl_(std::make_unique<Impl>()), model_path_(std::move(model_path)), device_str_(std::move(device)),
      cutoff_(cutoff) {
#ifdef CORRELATION_HAS_LIBTORCH
  try {
    if (device_str_ == "cuda" || device_str_.rfind("cuda:", 0) == 0) {
      if (torch::cuda::is_available()) {
        impl_->device = torch::Device(device_str_);
      } else {
        std::cerr << "[TorchGNNModel] CUDA requested but not available. Falling back to CPU.\n";
        impl_->device = torch::Device(torch::kCPU);
      }
    } else {
      impl_->device = torch::Device(torch::kCPU);
    }

    impl_->module = torch::jit::load(model_path_, impl_->device);
    impl_->module.eval();
    impl_->loaded = true;
  } catch (const c10::Error &e) {
    std::cerr << "[TorchGNNModel] Failed to load TorchScript model: " << e.what() << "\n";
    impl_->loaded = false;
  } catch (const std::exception &e) {
    std::cerr << "[TorchGNNModel] Error initializing model: " << e.what() << "\n";
    impl_->loaded = false;
  }
#else
  std::cerr << "[TorchGNNModel] LibTorch support was not compiled in (-DCORRELATION_ENABLE_LIBTORCH=OFF).\n";
#endif
}

TorchGNNModel::~TorchGNNModel() = default;
TorchGNNModel::TorchGNNModel(TorchGNNModel &&) noexcept = default;
TorchGNNModel &TorchGNNModel::operator=(TorchGNNModel &&) noexcept = default;

std::string TorchGNNModel::getModelName() const { return "TorchGNN [" + model_path_ + "]"; }

bool TorchGNNModel::isLoaded() const noexcept { return impl_ && impl_->loaded; }

MLIPOutput TorchGNNModel::evaluate(const correlation::core::Cell &cell) const {
  MLIPOutput output;
  const size_t number_atoms = cell.atoms().size();
  if (number_atoms == 0) {
    return output;
  }

#ifdef CORRELATION_HAS_LIBTORCH
  if (!isLoaded()) {
    throw std::runtime_error("[TorchGNNModel] Model is not loaded or failed initialization.");
  }

  // 1. Build periodic neighbor graph with PBC shifts
  const auto graph = PeriodicGraphBuilder::buildGraph(cell, cutoff_, false);
  const size_t E = graph.edge_count;

  // 2. Select tensor precision based on real_t
  const torch::ScalarType tensor_dtype = (sizeof(real_t) == sizeof(double)) ? torch::kFloat64 : torch::kFloat32;

  // 3. Construct input tensors (zero-copy from flat buffers into torch tensors)
  torch::Tensor pos = torch::from_blob(const_cast<real_t *>(graph.positions_flat.data()),
                                       {static_cast<int64_t>(number_atoms), 3}, tensor_dtype)
                          .to(impl_->device);

  torch::Tensor atomic_numbers = torch::from_blob(const_cast<int64_t *>(graph.atomic_numbers.data()),
                                                  {static_cast<int64_t>(number_atoms)}, torch::kInt64)
                                     .to(impl_->device);

  torch::Tensor edge_index;
  if (E > 0) {
    edge_index = torch::from_blob(const_cast<int64_t *>(graph.edge_index_flat.data()), {2, static_cast<int64_t>(E)},
                                  torch::kInt64)
                     .to(impl_->device);
  } else {
    edge_index = torch::empty({2, 0}, torch::TensorOptions().dtype(torch::kInt64).device(impl_->device));
  }

  torch::Tensor edge_shift;
  if (E > 0) {
    edge_shift = torch::from_blob(const_cast<real_t *>(graph.edge_shifts_flat.data()), {static_cast<int64_t>(E), 3},
                                  tensor_dtype)
                     .to(impl_->device);
  } else {
    edge_shift = torch::empty({0, 3}, torch::TensorOptions().dtype(tensor_dtype).device(impl_->device));
  }

  torch::Tensor cell_t =
      torch::from_blob(const_cast<real_t *>(graph.cell_flat.data()), {3, 3}, tensor_dtype).to(impl_->device);

  // 4. Run forward pass
  std::vector<torch::jit::IValue> inputs;
  inputs.emplace_back(pos);
  inputs.emplace_back(atomic_numbers);
  inputs.emplace_back(edge_index);
  inputs.emplace_back(edge_shift);
  inputs.emplace_back(cell_t);

  torch::NoGradGuard no_grad;
  auto result_ivalue = impl_->module.forward(inputs);

  // 5. Unpack outputs
  if (result_ivalue.isGenericDict()) {
    auto dict = result_ivalue.toGenericDict();

    // Extract Forces
    if (dict.contains("forces")) {
      auto forces_t = dict.at("forces").toTensor().to(torch::kCPU).to(tensor_dtype).contiguous();
      output.forces.resize(number_atoms);
      const auto *f_ptr = forces_t.data_ptr<real_t>();
      for (size_t i = 0; i < number_atoms; ++i) {
        output.forces[i] = correlation::math::Vector3<real_t>{f_ptr[i * 3 + 0], f_ptr[i * 3 + 1], f_ptr[i * 3 + 2]};
      }
    }

    // Extract LDOS
    if (dict.contains("ldos")) {
      auto ldos_t = dict.at("ldos").toTensor().to(torch::kCPU).to(tensor_dtype).contiguous();
      const auto sizes = ldos_t.sizes();
      if (sizes.size() == 2 && sizes[0] == static_cast<int64_t>(number_atoms)) {
        const size_t n_bins = static_cast<size_t>(sizes[1]);
        output.ldos_bins = n_bins;
        output.ldos.resize(number_atoms, std::vector<real_t>(n_bins, static_cast<real_t>(0.0)));
        const auto *ldos_ptr = ldos_t.data_ptr<real_t>();
        for (size_t i = 0; i < number_atoms; ++i) {
          for (size_t b = 0; b < n_bins; ++b) {
            output.ldos[i][b] = ldos_ptr[i * n_bins + b];
          }
        }
      }
    }

    // Extract Total Energy
    if (dict.contains("energy")) {
      auto energy_t = dict.at("energy").toTensor().to(torch::kCPU).to(tensor_dtype);
      output.total_energy = energy_t.item<real_t>();
    }
  }

#else
  // Fallback if LibTorch was not enabled at compile-time
  output.total_energy = static_cast<real_t>(0.0);
  output.per_atom_energy.resize(number_atoms, static_cast<real_t>(0.0));
  output.forces.resize(number_atoms, correlation::math::Vector3<real_t>{0.0, 0.0, 0.0});
  output.stress = correlation::math::Matrix3<real_t>();
#endif

  return output;
}

} // namespace correlation::mlip

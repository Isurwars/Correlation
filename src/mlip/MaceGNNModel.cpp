/**
 * @file MaceGNNModel.cpp
 * @brief Implementation of MaceGNNModel for LibTorch MACE evaluation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "mlip/MaceGNNModel.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <iostream>

#ifdef CORRELATION_HAS_LIBTORCH
#include <torch/script.h>
#include <torch/torch.h>
#endif

namespace correlation::mlip {

struct MaceGNNModel::Impl {
#ifdef CORRELATION_HAS_LIBTORCH
  torch::jit::script::Module module;
  torch::Device device{torch::kCPU};
  bool loaded{false};
#else
  bool loaded{false};
#endif
};

MaceGNNModel::MaceGNNModel(std::string model_path, std::string device, real_t cutoff)
    : impl_(std::make_unique<Impl>()), model_path_(std::move(model_path)),
      device_str_(std::move(device)), cutoff_(cutoff) {
#ifdef CORRELATION_HAS_LIBTORCH
  try {
    if (device_str_ == "cuda" || device_str_.rfind("cuda:", 0) == 0) {
      if (torch::cuda::is_available()) {
        impl_->device = torch::Device(device_str_);
      } else {
        std::cerr << "[MaceGNNModel] CUDA requested but not available. Falling back to CPU.\n";
        impl_->device = torch::Device(torch::kCPU);
      }
    } else {
      impl_->device = torch::Device(torch::kCPU);
    }

    impl_->module = torch::jit::load(model_path_, impl_->device);
    impl_->module.eval();
    impl_->loaded = true;
  } catch (const c10::Error &e) {
    std::cerr << "[MaceGNNModel] Failed to load TorchScript model: " << e.what() << "\n";
    impl_->loaded = false;
  } catch (const std::exception &e) {
    std::cerr << "[MaceGNNModel] Error initializing model: " << e.what() << "\n";
    impl_->loaded = false;
  }
#else
  std::cerr
      << "[MaceGNNModel] LibTorch support was not compiled in (-DCORRELATION_ENABLE_LIBTORCH=OFF).\n";
#endif
}

MaceGNNModel::~MaceGNNModel() = default;
MaceGNNModel::MaceGNNModel(MaceGNNModel &&) noexcept = default;
MaceGNNModel &MaceGNNModel::operator=(MaceGNNModel &&) noexcept = default;

std::string MaceGNNModel::getModelName() const { return "MACE [" + model_path_ + "]"; }

bool MaceGNNModel::isLoaded() const noexcept { return impl_ && impl_->loaded; }

#ifdef CORRELATION_HAS_LIBTORCH
namespace {

void extractMaceForces(const c10::Dict<std::string, torch::Tensor> &dict,
                       torch::ScalarType tensor_dtype, size_t number_atoms, MLIPOutput &output) {
  if (!dict.contains("forces")) {
    return;
  }
  const auto forces_t = dict.at("forces").to(torch::kCPU).to(tensor_dtype).contiguous();
  output.forces.resize(number_atoms);
  const auto *f_ptr = forces_t.data_ptr<real_t>();
  for (size_t i = 0; i < number_atoms; ++i) {
    output.forces[i] =
        correlation::math::Vector3<real_t>{f_ptr[i * 3 + 0], f_ptr[i * 3 + 1], f_ptr[i * 3 + 2]};
  }
}

void extractMaceEnergy(const c10::Dict<std::string, torch::Tensor> &dict,
                       torch::ScalarType tensor_dtype, MLIPOutput &output) {
  if (dict.contains("energy")) {
    const auto energy_t = dict.at("energy").to(torch::kCPU).to(tensor_dtype).contiguous();
    output.total_energy = *energy_t.data_ptr<real_t>();
  }

  if (dict.contains("node_energy")) {
    const auto node_e_t = dict.at("node_energy").to(torch::kCPU).to(tensor_dtype).contiguous();
    const size_t n_atoms = static_cast<size_t>(node_e_t.numel());
    output.per_atom_energy.resize(n_atoms);
    const auto *e_ptr = node_e_t.data_ptr<real_t>();
    std::copy_n(e_ptr, n_atoms, output.per_atom_energy.begin());
  }
}

void extractMaceStress(const c10::Dict<std::string, torch::Tensor> &dict,
                       torch::ScalarType tensor_dtype, MLIPOutput &output) {
  std::string stress_key;
  if (dict.contains("stress")) {
    stress_key = "stress";
  } else if (dict.contains("virial")) {
    stress_key = "virial";
  }

  if (stress_key.empty()) {
    return;
  }

  const auto stress_t = dict.at(stress_key).to(torch::kCPU).to(tensor_dtype).contiguous();
  const auto *s_ptr = stress_t.data_ptr<real_t>();
  if (stress_t.numel() == 9) {
    for (size_t r = 0; r < 3; ++r) {
      for (size_t c = 0; c < 3; ++c) {
        output.stress(r, c) = s_ptr[r * 3 + c];
      }
    }
  }
}

} // namespace
#endif

MLIPOutput MaceGNNModel::evaluate(const correlation::core::Cell &cell) const {
  MLIPOutput output;
  const size_t number_atoms = cell.atoms().size();
  if (number_atoms == 0) {
    return output;
  }

#ifdef CORRELATION_HAS_LIBTORCH
  if (!isLoaded()) {
    throw std::runtime_error("[MaceGNNModel] Model is not loaded or failed initialization.");
  }

  const auto graph = PeriodicGraphBuilder::buildGraph(cell, cutoff_, false);
  const size_t E = graph.edge_count;

  const torch::ScalarType tensor_dtype =
      (sizeof(real_t) == sizeof(double)) ? torch::kFloat64 : torch::kFloat32;

  auto positions_buf = graph.positions_flat;
  torch::Tensor pos =
      torch::from_blob(positions_buf.data(), {static_cast<int64_t>(number_atoms), 3}, tensor_dtype)
          .to(impl_->device)
          .clone();

  auto atomic_numbers_buf = graph.atomic_numbers;
  torch::Tensor atomic_numbers =
      torch::from_blob(atomic_numbers_buf.data(), {static_cast<int64_t>(number_atoms)},
                       torch::kInt64)
          .to(impl_->device)
          .clone();

  torch::Tensor edge_index;
  if (E > 0) {
    auto edge_index_buf = graph.edge_index_flat;
    edge_index =
        torch::from_blob(edge_index_buf.data(), {2, static_cast<int64_t>(E)}, torch::kInt64)
            .to(impl_->device)
            .clone();
  } else {
    edge_index =
        torch::empty({2, 0}, torch::TensorOptions().dtype(torch::kInt64).device(impl_->device));
  }

  torch::Tensor unit_shifts;
  if (E > 0) {
    auto edge_shifts_buf = graph.edge_shifts_flat;
    unit_shifts =
        torch::from_blob(edge_shifts_buf.data(), {static_cast<int64_t>(E), 3}, tensor_dtype)
            .to(impl_->device)
            .clone();
  } else {
    unit_shifts =
        torch::empty({0, 3}, torch::TensorOptions().dtype(tensor_dtype).device(impl_->device));
  }

  auto cell_buf = graph.cell_flat;
  torch::Tensor cell_t =
      torch::from_blob(cell_buf.data(), {3, 3}, tensor_dtype).to(impl_->device).clone();

  c10::Dict<std::string, torch::Tensor> input_dict;
  input_dict.insert("positions", pos);
  input_dict.insert("node_attrs", atomic_numbers);
  input_dict.insert("edge_index", edge_index);
  input_dict.insert("unit_shifts", unit_shifts);
  input_dict.insert("shifts", unit_shifts);
  input_dict.insert("cell", cell_t);

  torch::NoGradGuard no_grad;
  std::vector<torch::jit::IValue> inputs{input_dict};
  auto result_ivalue = impl_->module.forward(inputs);

  if (result_ivalue.isGenericDict()) {
    auto dict = result_ivalue.toGenericDict();
    c10::Dict<std::string, torch::Tensor> tensor_dict;
    for (const auto &entry : dict) {
      if (entry.key().isString() && entry.value().isTensor()) {
        tensor_dict.insert(entry.key().toStringRef(), entry.value().toTensor());
      }
    }
    extractMaceForces(tensor_dict, tensor_dtype, number_atoms, output);
    extractMaceEnergy(tensor_dict, tensor_dtype, output);
    extractMaceStress(tensor_dict, tensor_dtype, output);
  }
#else
  // Fallback when compiled without LibTorch
  output.forces.resize(number_atoms, correlation::math::Vector3<real_t>{});
#endif

  return output;
}

} // namespace correlation::mlip

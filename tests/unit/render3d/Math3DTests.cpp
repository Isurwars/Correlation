/**
 * @file Math3DTests.cpp
 * @brief Unit tests for 3D transformation mathematics (Vector4, Matrix4, Quaternion).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/Math3D.hpp"
#include <gtest/gtest.h>
#include <numbers>

namespace correlation::render3d {

TEST(Math3DTests, Vector4BasicOperations) {
  const Vector4<float> vec1(1.0F, 2.0F, 3.0F, 1.0F);
  const Vector4<float> vec2(4.0F, 5.0F, 6.0F, 1.0F);

  const auto sum = vec1 + vec2;
  EXPECT_FLOAT_EQ(sum.x(), 5.0F);
  EXPECT_FLOAT_EQ(sum.y(), 7.0F);
  EXPECT_FLOAT_EQ(sum.z(), 9.0F);
  EXPECT_FLOAT_EQ(sum.w(), 2.0F);

  const auto diff = vec2 - vec1;
  EXPECT_FLOAT_EQ(diff.x(), 3.0F);
  EXPECT_FLOAT_EQ(diff.y(), 3.0F);
  EXPECT_FLOAT_EQ(diff.z(), 3.0F);
  EXPECT_FLOAT_EQ(diff.w(), 0.0F);

  const auto scaled = vec1 * 2.0F;
  EXPECT_FLOAT_EQ(scaled.x(), 2.0F);
  EXPECT_FLOAT_EQ(scaled.y(), 4.0F);
  EXPECT_FLOAT_EQ(scaled.z(), 6.0F);
  EXPECT_FLOAT_EQ(scaled.w(), 2.0F);

  EXPECT_THROW([[maybe_unused]] auto div_zero = vec1 / 0.0F, std::runtime_error);

  const auto point = vec1.toVector3();
  EXPECT_FLOAT_EQ(point.x(), 1.0F);
  EXPECT_FLOAT_EQ(point.y(), 2.0F);
  EXPECT_FLOAT_EQ(point.z(), 3.0F);
}

TEST(Math3DTests, Matrix4IdentityAndMultiplication) {
  const auto ident = Matrix4<float>::identity();
  const Vector4<float> vec(3.0F, 4.0F, 5.0F, 1.0F);

  const auto transformed = ident * vec;
  EXPECT_FLOAT_EQ(transformed.x(), 3.0F);
  EXPECT_FLOAT_EQ(transformed.y(), 4.0F);
  EXPECT_FLOAT_EQ(transformed.z(), 5.0F);
  EXPECT_FLOAT_EQ(transformed.w(), 1.0F);

  const auto ident2 = ident * ident;
  for (std::size_t row_idx = 0; row_idx < 4; ++row_idx) {
    for (std::size_t col_idx = 0; col_idx < 4; ++col_idx) {
      EXPECT_FLOAT_EQ(ident(row_idx, col_idx), ident2(row_idx, col_idx));
    }
  }
}

TEST(Math3DTests, Matrix4TranslationAndScaling) {
  const auto trans_mat = Matrix4<float>::translation(10.0F, -5.0F, 2.0F);
  const math::Vector3<float> point(1.0F, 2.0F, 3.0F);
  const auto transformed_pt = trans_mat.transformPoint(point);

  EXPECT_FLOAT_EQ(transformed_pt.x(), 11.0F);
  EXPECT_FLOAT_EQ(transformed_pt.y(), -3.0F);
  EXPECT_FLOAT_EQ(transformed_pt.z(), 5.0F);

  const auto scale_mat = Matrix4<float>::scaling(2.0F, 3.0F, 4.0F);
  const auto scaled_pt = scale_mat.transformPoint(point);
  EXPECT_FLOAT_EQ(scaled_pt.x(), 2.0F);
  EXPECT_FLOAT_EQ(scaled_pt.y(), 6.0F);
  EXPECT_FLOAT_EQ(scaled_pt.z(), 12.0F);
}

TEST(Math3DTests, Matrix4LookAtAndViewMatrix) {
  const math::Vector3<float> eye(0.0F, 0.0F, 10.0F);
  const math::Vector3<float> target(0.0F, 0.0F, 0.0F);
  const math::Vector3<float> up_direction(0.0F, 1.0F, 0.0F);

  const auto view = Matrix4<float>::lookAt(eye, target, up_direction);
  const auto eye_transformed = view.transformPoint(eye);

  EXPECT_NEAR(eye_transformed.x(), 0.0F, 1e-5F);
  EXPECT_NEAR(eye_transformed.y(), 0.0F, 1e-5F);
  EXPECT_NEAR(eye_transformed.z(), 0.0F, 1e-5F);

  const auto target_transformed = view.transformPoint(target);
  EXPECT_NEAR(target_transformed.x(), 0.0F, 1e-5F);
  EXPECT_NEAR(target_transformed.y(), 0.0F, 1e-5F);
  EXPECT_NEAR(target_transformed.z(), -10.0F, 1e-5F);
}

TEST(Math3DTests, Matrix4PerspectiveProjection) {
  const float fov = 45.0F * (std::numbers::pi_v<float> / 180.0F);
  const float aspect = 1.3333F;
  const float near_plane = 0.1F;
  const float far_plane = 100.0F;

  const auto proj = Matrix4<float>::perspective(fov, aspect, near_plane, far_plane);
  EXPECT_THROW([[maybe_unused]] const auto invalid =
                   Matrix4<float>::perspective(fov, 0.0F, near_plane, far_plane),
               std::invalid_argument);
  EXPECT_THROW([[maybe_unused]] const auto invalid2 =
                   Matrix4<float>::perspective(fov, aspect, -1.0F, far_plane),
               std::invalid_argument);

  const math::Vector3<float> pt_near(0.0F, 0.0F, -near_plane);
  const auto proj_near = proj.transformPoint(pt_near);
  EXPECT_NEAR(proj_near.z(), 0.0F, 1e-4F);
}

TEST(Math3DTests, QuaternionRotation) {
  const math::Vector3<float> axis(0.0F, 1.0F, 0.0F);
  const float angle = 90.0F * (std::numbers::pi_v<float> / 180.0F);
  const auto quat = Quaternion<float>::fromAxisAngle(axis, angle);

  const math::Vector3<float> test_vec(1.0F, 0.0F, 0.0F);
  const auto rotated = quat.rotate(test_vec);

  EXPECT_NEAR(rotated.x(), 0.0F, 1e-5F);
  EXPECT_NEAR(rotated.y(), 0.0F, 1e-5F);
  EXPECT_NEAR(rotated.z(), -1.0F, 1e-5F);

  const auto rot_matrix = quat.toMatrix4();
  const auto matrix_rotated = rot_matrix.transformPoint(test_vec);
  EXPECT_NEAR(matrix_rotated.x(), 0.0F, 1e-5F);
  EXPECT_NEAR(matrix_rotated.y(), 0.0F, 1e-5F);
  EXPECT_NEAR(matrix_rotated.z(), -1.0F, 1e-5F);
}

} // namespace correlation::render3d

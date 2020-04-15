#pragma once

template<typename T>
inline
Eigen::Matrix<T, Eigen::Dynamic, 1>
csr_log_sum_exp(int m, int n,
                const std::vector<int>& v,
                const std::vector<int>& u,
                const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                std::ostream* pstream__) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> result(m);
  result.setZero();

  for (int row = 0; row < m; ++row) {
    int idx = csr_u_to_z(u, row); // u[row + 1] - u[row]
    int row_end_in_w = (u[row] - stan::error_index::value) + idx;

    Eigen::Matrix<T, Eigen::Dynamic, 1> b_sub(idx);
    b_sub.setZero();

    int i = 0;

    for (int nze = u[row] - stan::error_index::value; nze < row_end_in_w; ++nze, ++i) {
      b_sub.coeffRef(i) = b.coeffRef(v[nze] - stan::error_index::value);
    }

    const double max = b_sub.maxCoeff();

    result.coeffRef(row) = max + std::log((b_sub.array() - max).exp().sum());
  }

  return result;
}

// template<typename T_MATRIX, typename T_VECTOR>
template<typename T>
inline
// Eigen::Matrix<stan::return_type_t<T_MATRIX, T_VECTOR>, Eigen::Dynamic, 1>
Eigen::Matrix<T, Eigen::Dynamic, 1>
csr_log_sum_exp2(int m, int n,
                 const std::vector<int>& v,
                 const std::vector<int>& u,
                 const Eigen::Matrix<T, Eigen::Dynamic, 1>& w,
                 const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                 std::ostream* pstream__) {
  // using result_t = stan::return_type_t<T_MATRIX, T_VECTOR>;

  Eigen::Matrix<T, Eigen::Dynamic, 1> result(m);
  result.setZero();

  for (int row = 0; row < m; ++row) {
    int idx = csr_u_to_z(u, row); // u[row + 1] - u[row]
    int row_end_in_w = (u[row] - stan::error_index::value) + idx;

    Eigen::Matrix<T, Eigen::Dynamic, 1> b_sub(idx);
    b_sub.setZero();

    int i = 0;

    for (int nze = u[row] - stan::error_index::value; nze < row_end_in_w; ++nze, ++i) {
      b_sub.coeffRef(i) = b.coeffRef(v[nze] - stan::error_index::value);
    }

    const double max = b_sub.maxCoeff();

    Eigen::Matrix<T, Eigen::Dynamic, 1> w_sub(w.segment(u[row] - stan::error_index::value, idx));

    result.coeffRef(row) = stan::math::dot_product(w_sub, max + std::log((b_sub.array() - max).exp().sum()));
  }

  return result;
}

/*
template <typename T1, typename T2>
inline Eigen::Matrix<return_type_t<T1, T2>, Eigen::Dynamic, 1>
my_csr_matrix_times_vector(int m, int n,
                        const Eigen::Matrix<T1, Eigen::Dynamic, 1>& w,
                        const std::vector<int>& v, const std::vector<int>& u,
                        const Eigen::Matrix<T2, Eigen::Dynamic, 1>& b) {
  using result_t = return_type_t<T1, T2>;

  check_positive("csr_matrix_times_vector", "m", m);
  check_positive("csr_matrix_times_vector", "n", n);
  check_size_match("csr_matrix_times_vector", "n", n, "b", b.size());
  check_size_match("csr_matrix_times_vector", "m", m, "u", u.size() - 1);
  check_size_match("csr_matrix_times_vector", "w", w.size(), "v", v.size());
  check_size_match("csr_matrix_times_vector", "u/z",
                   u[m - 1] + csr_u_to_z(u, m - 1) - 1, "v", v.size());
  for (int i : v) {
    check_range("csr_matrix_times_vector", "v[]", n, i);
  }

  Eigen::Matrix<result_t, Eigen::Dynamic, 1> result(m);
  result.setZero();
  for (int row = 0; row < m; ++row) {
    int idx = csr_u_to_z(u, row);
    int row_end_in_w = (u[row] - stan::error_index::value) + idx;
    int i = 0;
    Eigen::Matrix<result_t, Eigen::Dynamic, 1> b_sub(idx);
    b_sub.setZero();
    for (int nze = u[row] - stan::error_index::value; nze < row_end_in_w;
         ++nze, ++i) {
      check_range("csr_matrix_times_vector", "j", n, v[nze]);
      b_sub.coeffRef(i) = b.coeffRef(v[nze] - stan::error_index::value);
    }
    Eigen::Matrix<T1, Eigen::Dynamic, 1> w_sub(
      w.segment(u[row] - stan::error_index::value, idx));
    result.coeffRef(row) = dot_product(w_sub, b_sub);
  }
  return result;
}
*/

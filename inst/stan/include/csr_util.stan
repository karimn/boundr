// vector csr_log_sum_exp(int m, int n, int[] v, int[] u, vector b);
// vector csr_log_sum_exp2(int m, int n, vector w, int[] v, int[] u, vector b);

vector csr_log_sum_exp(int m, int n, int[] v, int[] u, vector b) {
  vector[m] result = rep_vector(0, m);

  for (row_index in 1:m) {
    int row_size = u[row_index + 1] - u[row_index];

    if (row_size > 0) {
      int row_pos = u[row_index];
      int row_end = row_pos + row_size - 1;

      vector[row_size] sub_b = b[v[row_pos:row_end]];

      result[row_index] = log_sum_exp(sub_b);
    }
  }

  return result;
}

vector csr_log_sum_exp2(int m, int n, vector w, int[] v, int[] u, vector b) {
  vector[m] result = rep_vector(0, m);

  for (row_index in 1:m) {
    int row_size = u[row_index + 1] - u[row_index];

    if (row_size > 0) {
      int row_pos = u[row_index];
      int row_end = row_pos + row_size - 1;

      vector[row_size] sub_w = w[row_pos:row_end];
      vector[row_size] sub_b = b[v[row_pos:row_end]];

      real max_sub_b = max(sub_b);

      result[row_index] = max_sub_b + log(dot_product(sub_w, exp(sub_b - max_sub_b)));
    }
  }

  return result;
}

vector csr_diff_exp(int m, int n, int[] v, vector b) {
  vector[m] result = rep_vector(0, m);
  int row_pos = 1;

  for (row_index in 1:m) {
    int row_end = row_pos + 1;

    // vector[row_size] sub_w = w[row_pos:row_end];
    vector[2] sub_b = b[v[row_pos:row_end]];

    // real max_sub_b = max(sub_b);

    result[row_index] = exp(sub_b[1]) - exp(sub_b[2]);

    row_pos += 2;
  }

  return result;
}

int[] csr_shift_expand_v(int[] indices, int shift_size, int num_shifts) {
  int num_indices = num_elements(indices);
  int new_csr_v[num_indices * num_shifts];

  int id_pos = 1;

  for (shift_index in 1:num_shifts) {
    int id_end = id_pos + num_indices - 1;

    new_csr_v[id_pos:id_end] = array_add(indices, { (shift_index - 1) * shift_size });

    id_pos = id_end + 1;
  }

  return new_csr_v;
}

int[] csr_shift_expand_u(int[] row_sizes, int num_shifts) {
  int num_rows = num_elements(row_sizes);
  int new_csr_u[num_rows * num_shifts + 1]; // The last element is for padding

  int row_pos_pos = 1;
  int last_size;

  for (shift_index in 1:num_shifts) {
    for (row_index in 1:num_rows) {
      new_csr_u[row_pos_pos] = row_pos_pos > 1 ? new_csr_u[row_pos_pos - 1] + last_size : 1;
      last_size = row_sizes[row_index];

      row_pos_pos += 1;
    }
  }

  new_csr_u[row_pos_pos] = new_csr_u[row_pos_pos - 1] + last_size;

  return new_csr_u;
}

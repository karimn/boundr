functions {
  // These functions were in the multilvlr submodule

  int num_test(int[] to_test, int[] target_val, int test_equality) {
    int num_to_test = num_elements(to_test);
    int num_targets = num_elements(target_val);
    int result = 0;

    int sorted_to_test[num_to_test] = sort_asc(to_test);

    for (to_test_index in 1:num_to_test) {
      int found = 0;

      for (target_index in 1:num_targets) {
        if (sorted_to_test[to_test_index] == target_val[target_index]) {
          if (test_equality) {
            result += 1;
          }

          found = 1;
          break;
        }
      }

      if (!found && (1 - test_equality)) {
        result += 1;
      }
    }

    return(result);
  }

  int num_equals(int[] to_test, int[] target_val) {
    return(num_test(to_test, target_val, 1));
  }

  int[] unique(int[] find_in) {
    int num_to_tally = num_elements(find_in);
    int sorted_find_in[num_to_tally] = sort_asc(find_in);
    int unique_count = 1;
    int unique_found[num_to_tally];

    unique_found[1] = sorted_find_in[1];

    for (tally_index in 2:num_to_tally) {
      if (sorted_find_in[tally_index] != unique_found[unique_count]) {
        unique_count += 1;
        unique_found[unique_count] = sorted_find_in[tally_index];
      }
    }

    return(unique_found[1:unique_count]);
  }

  int num_unique(int[] find_in) {
    return(num_elements(unique(find_in)));
  }

  int[] count(int count_size, int[] find_in) {
    int count_array[count_size] = rep_array(0, count_size);

    for (count_index in 1:count_size) {
      count_array[count_index] = num_equals(find_in, { count_index });
    }

    return(count_array);
  }

  int[] array_add(int[] left, int[] right) {
    int array_size = num_elements(left);
    int array_sum[array_size];
    int right_array_size = num_elements(right);

    if (right_array_size != array_size && right_array_size != 1) {
      reject("Incompatible array sizes.");
    }

    for (array_index in 1:array_size) {
      array_sum[array_index] = left[array_index] + right[right_array_size > 1 ? array_index : 1];
    }

    return(array_sum);
  }

  int[] extract_group_pos_end(int[] group_sizes, int group_to_extract) {
    int group_pos = group_to_extract > 1 ? sum(group_sizes[1:(group_to_extract - 1)]) + 1 : 1;
    int group_end = group_pos + group_sizes[group_to_extract] - 1;

    return({ group_pos, group_end });
  }

  int[] array_extract_group_values(int[] all_values, int[] group_sizes, int[] groups_to_extract) {
    int num_groups_to_extract = num_elements(groups_to_extract);
    int num_values[num_groups_to_extract] = group_sizes[groups_to_extract];
    int extracted[sum(num_values)];

    int extracted_pos = 1;

    for (group_index in 1:num_groups_to_extract) {
      int extracted_end = extracted_pos + num_values[group_index] - 1;
      int actual_curr_group_index = groups_to_extract[group_index];

      int group_pos_end[2] = extract_group_pos_end(group_sizes, actual_curr_group_index);

      extracted[extracted_pos:extracted_end] = all_values[group_pos_end[1]:group_pos_end[2]];

      extracted_pos = extracted_end + 1;
    }

    return(extracted);
  }

  int[] array_product(int[] left, int[] right) {
    int array_size = num_elements(left);
    int array_prod[array_size];
    int right_array_size = num_elements(right);

    if (right_array_size != array_size && right_array_size != 1) {
      reject("Incompatible array sizes.");
    }

    for (array_index in 1:array_size) {
      array_prod[array_index] = left[array_index] * right[right_array_size > 1 ? array_index : 1];
    }

    return(array_prod);
  }

  int[] seq(int from, int to, int by) {
    int reverse = from > to;
    int seq_len = (((1 - 2 * reverse) * (to - from)) / by) + 1;
    int result_seq[seq_len];

    for (seq_index in 1:seq_len) {
      result_seq[seq_index] =  from + (1 - 2 * reverse) * ((seq_index - 1) * by);
    }

    return(result_seq);
  }

  // Below functions only defined in this file

  int num_gt_zero(vector v) {
    int num_found = 0;

    for (i in 1:num_elements(v)) {
      if (v[i] > 0) {
        num_found += 1;
      }
    }

    return num_found;
  }

  int[] which_compare_zero(vector v, int num_found, int gt) {
    int found_indices[num_found];
    int found_pos = 1;

    for (i in 1:num_elements(v)) {
      if (gt * (v[i] > 0) + (1 - gt) * (v[i] <= 0)) {
        found_indices[found_pos] = i;
        found_pos += 1;
      }
    }

    return found_indices;
  }

  int[] calculate_level_size(int num_levels, int[,] unique_entity_ids) {
    // int num_levels = dims(unique_entity_ids)[2];
    int level_size[num_levels];

    for (level_index in 1:num_levels) {
      level_size[level_index] = num_unique(unique_entity_ids[, level_index]);
    }

    return level_size;
  }

  int[] calculate_num_in_level_entities(int[,] unique_entity_ids, int[] level_size) {
    int num_levels = num_elements(level_size);
    int num_unique_entities_in_level_entities[sum(level_size)];

    int level_entity_pos = 1;

    for (level_index in 1:num_levels) {
      int level_entity_end = level_entity_pos + level_size[level_index] - 1;

      num_unique_entities_in_level_entities[level_entity_pos:level_entity_end] = count(level_size[level_index], unique_entity_ids[, level_index]);

      level_entity_pos = level_entity_end + 1;
    }

    return num_unique_entities_in_level_entities;
  }

  int[] calculate_entity_total_num_candidates(int[] num_unique_entity_candidate_groups, int[] unique_entity_candidate_groups, int[] candidate_group_size) {
    int num_unique_entities = num_elements(num_unique_entity_candidate_groups);
    int entity_total_num_candidates[num_unique_entities];

    int entity_candidate_group_pos = 1;

    for (entity_index in 1:num_unique_entities) {
      int num_entity_candidate_groups = num_unique_entity_candidate_groups[entity_index];
      int entity_candidate_group_end = entity_candidate_group_pos + num_entity_candidate_groups - 1;

      int curr_candidate_groups[num_entity_candidate_groups] = unique_entity_candidate_groups[entity_candidate_group_pos:entity_candidate_group_end];

      entity_total_num_candidates[entity_index] = sum(candidate_group_size[curr_candidate_groups]);

      entity_candidate_group_pos = entity_candidate_group_end + 1;
    }

    return entity_total_num_candidates;
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

  matrix rep_corr_estimand_rng(matrix r_prob_mat, int[] rep_corr_outcomes, int[] rep_corr_cond, int[,] experiment_assign_entity, matrix[] type_response_value) {
    int num_r_types = rows(r_prob_mat);
    int num_rep_corr_estimands = num_elements(rep_corr_outcomes) / 2;
    int num_experiment_entities = size(experiment_assign_entity);

    matrix[num_rep_corr_estimands, 2] rep_corr_estimand;

    // https://en.wikipedia.org/wiki/Correlation_and_dependence

    matrix[num_rep_corr_estimands, 2] rep_outcome_n = rep_matrix(0, num_rep_corr_estimands, 2); // n
    matrix[num_rep_corr_estimands, 2] rep_outcome_product_sum = rep_matrix(0, num_rep_corr_estimands, 2); // \sum_i x_i y_i
    matrix[num_rep_corr_estimands, 2] rep_outcome1_sum = rep_matrix(0, num_rep_corr_estimands, 2); // \sum_i x_i
    matrix[num_rep_corr_estimands, 2] rep_outcome2_sum = rep_matrix(0, num_rep_corr_estimands, 2); // \sum_i y_i
    matrix[num_rep_corr_estimands, 2] rep_outcome1_squared_sum = rep_matrix(0, num_rep_corr_estimands, 2); // \sum_i x_i^2
    matrix[num_rep_corr_estimands, 2] rep_outcome2_squared_sum = rep_matrix(0, num_rep_corr_estimands, 2); // \sum_i y_i^2

    for (entity_exp_index in 1:num_experiment_entities) {
      int entity_exp_id = experiment_assign_entity[entity_exp_index, 2];
      int entity_exp_size = experiment_assign_entity[entity_exp_index, 3];

      vector[num_r_types] rep_sample_types = to_vector(multinomial_rng(r_prob_mat[, experiment_assign_entity[entity_exp_index, 1]], entity_exp_size));

      for (rep_corr_index in 1:num_rep_corr_estimands) {
        int rep_corr_pos = 2 * (rep_corr_index - 1) + 1;
        int outcome_id[2] = rep_corr_outcomes[rep_corr_pos:(rep_corr_pos + 1)];

        row_vector[num_r_types] outcome1 = type_response_value[outcome_id[1], entity_exp_id];
        row_vector[num_r_types] outcome2 = type_response_value[outcome_id[2], entity_exp_id];
        row_vector[num_r_types] cond = type_response_value[rep_corr_cond[rep_corr_index], entity_exp_id];

        rep_outcome_n[rep_corr_index, 1] += (1 - cond) * rep_sample_types;
        rep_outcome_n[rep_corr_index, 2] += cond * rep_sample_types;

        rep_outcome_product_sum[rep_corr_index, 1] += (outcome1 .* outcome2 .* (1 - cond)) * rep_sample_types;
        rep_outcome_product_sum[rep_corr_index, 2] += (outcome1 .* outcome2 .* cond) * rep_sample_types;

        rep_outcome1_sum[rep_corr_index, 1] += (outcome1 .* (1 - cond)) * rep_sample_types;
        rep_outcome1_sum[rep_corr_index, 2] += (outcome1 .* cond) * rep_sample_types;

        rep_outcome2_sum[rep_corr_index, 1] += (outcome2 .* (1 - cond)) * rep_sample_types;
        rep_outcome2_sum[rep_corr_index, 2] += (outcome2 .* cond) * rep_sample_types;

        rep_outcome1_squared_sum[rep_corr_index, 1] += square(outcome1 .* (1 - cond)) * rep_sample_types;
        rep_outcome1_squared_sum[rep_corr_index, 2] += square(outcome1 .* cond) * rep_sample_types;

        rep_outcome2_squared_sum[rep_corr_index, 1] += square(outcome2 .* (1 - cond)) * rep_sample_types;
        rep_outcome2_squared_sum[rep_corr_index, 2] += square(outcome2 .* cond) * rep_sample_types;
      }
    }

    for (rep_corr_index in 1:num_rep_corr_estimands) {
      rep_corr_estimand[rep_corr_index, 1] =
        ((rep_outcome_n[rep_corr_index, 1] * rep_outcome_product_sum[rep_corr_index, 1]) - (rep_outcome1_sum[rep_corr_index, 1] * rep_outcome2_sum[rep_corr_index, 1])) /
        (sqrt(rep_outcome_n[rep_corr_index, 1] * rep_outcome1_squared_sum[rep_corr_index, 1] - square(rep_outcome1_sum[rep_corr_index, 1])) *
         sqrt(rep_outcome_n[rep_corr_index, 1] * rep_outcome2_squared_sum[rep_corr_index, 1] - square(rep_outcome2_sum[rep_corr_index, 1])));

      if (is_nan(rep_corr_estimand[rep_corr_index, 1])) {
        rep_corr_estimand[rep_corr_index, 1] = 0;
      }

      rep_corr_estimand[rep_corr_index, 2] =
        ((rep_outcome_n[rep_corr_index, 2] * rep_outcome_product_sum[rep_corr_index, 2]) - (rep_outcome1_sum[rep_corr_index, 2] * rep_outcome2_sum[rep_corr_index, 2])) /
        (sqrt(rep_outcome_n[rep_corr_index, 2] * rep_outcome1_squared_sum[rep_corr_index, 2] - square(rep_outcome1_sum[rep_corr_index, 2])) *
         sqrt(rep_outcome_n[rep_corr_index, 2] * rep_outcome2_squared_sum[rep_corr_index, 2] - square(rep_outcome2_sum[rep_corr_index, 2])));

      if (is_nan(rep_corr_estimand[rep_corr_index, 2])) {
        rep_corr_estimand[rep_corr_index, 2] = 0;
      }
    }

    return rep_corr_estimand;
  }

 vector calculate_r_type_joint_prob(int num_r_types, int num_discrete_r_types, int num_discretized_r_types, int discrete_group_size,
                                    int[] num_compatible_discretized_r_types, int[] compatible_discretized_r_types, int[,] compatible_discretized_pair_ids,
                                    vector discrete_prob, matrix[] beta, int entity_index) {
   int num_discretized_variables = size(beta);
   vector[num_r_types] joint_prob;

   for (discrete_index in 1:num_discrete_r_types) {
     int discretized_pos = 1 + (discrete_index - 1) * num_discretized_r_types;
     int discretized_end = discrete_index * num_discretized_r_types;
     int discrete_group_pos = 1 + (discrete_index - 1) * discrete_group_size;
     int discrete_group_end = discrete_index * discrete_group_size;
     int r_prob_pos = discrete_group_pos;
     int r_prob_end = discrete_group_end;

     vector[num_discretized_r_types] first_discretized_prob = softmax(beta[1, discretized_pos:discretized_end, entity_index]);

     joint_prob[r_prob_pos:r_prob_end] = discrete_prob[discrete_index] * first_discretized_prob[compatible_discretized_pair_ids[1, discrete_group_pos:discrete_group_end]];

     for (discretized_var_index in 2:num_discretized_variables) {
       vector[sum(num_compatible_discretized_r_types)] curr_cond_discretized_prob;
       int cond_prob_pos = 1;

       for (discretized_type_index in 1:num_discretized_r_types) {
         int cond_prob_end = cond_prob_pos + num_compatible_discretized_r_types[discretized_type_index] - 1;
         int compatible_ids[num_compatible_discretized_r_types[discretized_type_index]] = compatible_discretized_r_types[cond_prob_pos:cond_prob_end];

         curr_cond_discretized_prob[cond_prob_pos:cond_prob_end] =  softmax(beta[discretized_var_index, discretized_pos:discretized_end, entity_index][compatible_ids]);

         cond_prob_pos = cond_prob_end + 1;
       }

       joint_prob[r_prob_pos:r_prob_end] =
         joint_prob[r_prob_pos:r_prob_end] .* curr_cond_discretized_prob[compatible_discretized_pair_ids[discretized_var_index, discrete_group_pos:discrete_group_end]];
     }
   }

   return joint_prob;
 }
}

data {
  int<lower = 0> num_obs;

  int<lower = 1> num_r_types;
  int<lower = 1> num_discrete_r_types;
  int<lower = 1, upper = num_discrete_r_types> discrete_r_type_id[num_r_types];

  // Discretized Variables

  int<lower = 0> num_cutpoints;
  vector[num_cutpoints] cutpoints;

  int<lower = 0> num_discretized_r_types;
  int<lower = 1, upper = num_discretized_r_types> num_compatible_discretized_r_types[num_discretized_r_types];
  int<lower = 1, upper = num_discretized_r_types> compatible_discretized_r_types[sum(num_compatible_discretized_r_types)];

  int<lower = 1, upper = sum(num_compatible_discretized_r_types)> compatible_discretized_pair_ids[max(0, num_cutpoints - 2), num_r_types];

  // Experimental Assignment

  int<lower = 1> num_experiment_types;
  simplex[num_experiment_types] experiment_types_prob;

  // Outcomes

  int<lower = 1> num_responses;
  // matrix[num_experiment_types, num_r_types] type_response_value[num_responses];

  // Background Variable Types

  int<lower = 1> num_bg_variables;
  int<lower = 1> num_bg_variable_types[num_bg_variables];
  int<lower = 1, upper = num_r_types> num_bg_variable_type_combo_members[sum(num_bg_variable_types)];
  int<lower = 1, upper = num_r_types> bg_variable_type_combo_members[sum(num_bg_variable_type_combo_members)];

  // Levels

  int<lower = 0> num_levels;

  int<lower = 1> num_unique_entities;
  int<lower = 1> unique_entity_ids[num_unique_entities, num_levels];
  int<lower = 1, upper = num_unique_entities> obs_unique_entity_id[num_obs];

  // For replication purposes, we need to know the various experimental assignments of the entity members
  // int<lower = num_experiment_types> num_experiment_entities;
  // int experiment_assign_entity[num_experiment_entities, 3]; // [unique entity, experiment assignment, # obs]

  // Candidate Groups

  int<lower = 1> num_candidate_groups;
  int<lower = 1, upper = num_r_types> candidate_group_size[num_candidate_groups];
  int<lower = 1, upper = num_r_types> candidate_group_ids[sum(candidate_group_size)];

  int<lower = 1, upper = num_candidate_groups> num_unique_entity_candidate_groups[num_obs > 0 ? num_unique_entities : 0];
  int<lower = 1, upper = num_candidate_groups> unique_entity_candidate_groups[sum(num_unique_entity_candidate_groups)];
  row_vector<lower = 1>[sum(num_unique_entity_candidate_groups)] num_unique_entity_in_candidate_groups;

  int<lower = 1, upper = num_candidate_groups> obs_candidate_group[num_obs];

  int<lower = 1, upper = num_obs> obs_in_unique_entity_in_candidate_groups[num_obs];

  // Estimation

  int<lower = 0> num_discrete_estimands;
  int<lower = 0, upper = num_discrete_estimands> num_atom_estimands;
  int<lower = 0, upper = num_discrete_estimands> num_diff_estimands;
  int<lower = 0, upper = num_atom_estimands> num_discretized_groups;
  int<lower = 0> num_mean_diff_estimands;
  int<lower = 0> num_utility_diff_estimands;

  int<lower = 1, upper = num_atom_estimands> diff_estimand_atoms[num_diff_estimands * 2];
  int<lower = 1, upper = num_discretized_groups> mean_diff_estimand_atoms[num_mean_diff_estimands * 2];
  int<lower = 1, upper = num_discretized_groups> utility_diff_estimand_atoms[num_utility_diff_estimands * 2];

  int<lower = 1, upper = num_atom_estimands> discretized_group_ids[num_discretized_groups * (num_cutpoints - 2)];

  int<lower = 0> num_abducted_estimands;
  int<lower = 1, upper = num_atom_estimands> abducted_estimand_ids[num_abducted_estimands];

  int<lower = 0, upper = num_levels> num_estimand_levels;
  int<lower = 1, upper = num_levels> estimand_levels[num_estimand_levels];

  int<lower = 0, upper = num_levels> num_between_entity_diff_levels;
  int<lower = 1, upper = num_levels> between_entity_diff_levels[num_between_entity_diff_levels];

  int<lower = 0> abducted_prob_size[num_abducted_estimands];
  int<lower = 1, upper = num_r_types * num_experiment_types> abducted_prob_index[sum(abducted_prob_size)];

  int<lower = 0> est_prob_size[num_atom_estimands];
  int<lower = 1, upper = num_r_types * num_experiment_types> est_prob_index[sum(est_prob_size)];

  int<lower = 0> num_discrete_utility_values;
  vector[num_discrete_utility_values] utility;

  // Replication

  // int<lower = 0> num_rep_corr_estimands;
  // int<lower = 1, upper = num_responses> rep_corr_outcomes[num_rep_corr_estimands * 2];
  // int<lower = 1, upper = num_responses> rep_corr_cond[num_rep_corr_estimands];

  // Cross Validation

  int<lower = -1, upper = num_levels> log_lik_level; /* -1 : no cross validation
                                                         0 : cross validation at the observation level
                                                        >0 : cross validation at specified level */

  // Priors

  real<lower = 0> discrete_beta_hyper_sd;
  real<lower = 0> discretized_beta_hyper_sd;

  vector<lower = 0>[num_levels] tau_level_sigma;

  // Configuration

  int<lower = 1, upper = 2> run_type;
  int<lower = 0, upper = 1> use_random_binpoint;
  int<lower = 0, upper = 1> generate_rep;
  int<lower = 0, upper = 1> calculate_marginal_prob;
}

transformed data {
  int RUN_TYPE_PRIOR_PREDICT = 1;
  int RUN_TYPE_FIT = 2;

  int num_all_estimands = num_discrete_estimands + (1 + (num_discrete_utility_values > 0)) * num_discretized_groups + num_mean_diff_estimands + num_utility_diff_estimands;

  int num_discretized_variables = max(num_cutpoints - 2, 0);

  int discrete_group_size = num_r_types / num_discrete_r_types;

  vector[max(num_cutpoints - 1, 0)] cutpoint_midpoints;
  vector[max(num_cutpoints - 1, 0)] discretized_binwidth;

  real discretize_bin_alpha[max(num_cutpoints - 1, 0)];
  real discretize_bin_beta[max(num_cutpoints - 1, 0)];

  real entity_discretize_bin_alpha[num_unique_entities * (num_cutpoints - 1) * num_discretized_groups];
  real entity_discretize_bin_beta[num_unique_entities * (num_cutpoints - 1) * num_discretized_groups];

  int<lower = num_r_types> num_r_types_full = num_r_types * num_experiment_types; // Unnested, with experimental variables uncollapsed
  int<lower = 1, upper = num_r_types> experiment_r_type_index[num_r_types_full];

  vector<lower = 0, upper = 1>[num_r_types_full] full_experiment_types_prob = to_vector(rep_matrix(experiment_types_prob, num_r_types));

  int<lower = 1> level_size[num_levels] = calculate_level_size(num_levels, unique_entity_ids);
  vector<lower = 0, upper = num_obs>[num_unique_entities] num_obs_in_unique_entity;
  simplex[max(1, num_unique_entities)] unique_entity_prop;

  int<lower = 1, upper = num_unique_entities> num_unique_entities_in_estimand_level_entities[sum(level_size[estimand_levels])] =
    calculate_num_in_level_entities(unique_entity_ids[, estimand_levels], level_size[estimand_levels]);

  int<lower = 1, upper = num_unique_entities> unique_entities_in_level_entities[num_unique_entities, num_levels];

  int<lower = 1, upper = num_obs> num_obs_in_log_lik_level_entities[log_lik_level > 0 ? level_size[log_lik_level] : 0];
  int<lower = 1, upper = num_obs> obs_in_log_lik_level_entities[log_lik_level > 0 ? num_obs : 0];

  vector<lower = -1, upper = 1>[num_diff_estimands * 2 * num_unique_entities] vec_diff = to_vector(rep_matrix([1, -1]', num_diff_estimands * num_unique_entities));
  vector<lower = -1, upper = 1>[num_mean_diff_estimands * 2 * num_unique_entities] vec_mean_diff = to_vector(rep_matrix([1, -1]', num_mean_diff_estimands * num_unique_entities));
  vector<lower = -1, upper = 1>[num_utility_diff_estimands * 2 * num_unique_entities] vec_utility_diff = to_vector(rep_matrix([1, -1]', num_utility_diff_estimands * num_unique_entities));

  int vec_1_size = max({ sum(candidate_group_size[unique_entity_candidate_groups]),
                         sum(candidate_group_size[obs_candidate_group]),
                         num_unique_entities * sum(abducted_prob_size), num_unique_entities * sum(est_prob_size),
                         num_unique_entities * num_r_types });

  vector<lower = 1, upper = 1>[vec_1_size] vec_1 = rep_vector(1, vec_1_size);

  int<lower = 1> entity_total_num_candidates[num_obs > 0 ? num_unique_entities : 0] =
    calculate_entity_total_num_candidates(num_unique_entity_candidate_groups, unique_entity_candidate_groups, candidate_group_size);

  // int<lower = 1, upper = num_r_types * num_unique_entities> entity_prob_ids[num_r_types * num_unique_entities];
  // int<lower = 1, upper = num_r_types * num_unique_entities + 1> entity_prob_csr_row_pos[num_unique_entities + 1];
  //
  // int<lower = 1, upper = num_unique_entities> prob_repeater_ids[num_r_types * num_unique_entities];
  // int<lower = 1, upper = num_r_types * num_unique_entities + 1> prob_repeater_csr_row_pos[num_r_types * num_unique_entities + 1];

  int<lower = 1, upper = num_r_types * num_unique_entities> entity_candidate_group_ids[sum(entity_total_num_candidates)];
  int<lower = 1, upper = sum(candidate_group_size[unique_entity_candidate_groups]) + 1> entity_candidate_group_csr_row_pos[sum(num_unique_entity_candidate_groups) + 1];

  int<lower = 1, upper = num_r_types * num_unique_entities> obs_candidate_group_ids[sum(candidate_group_size[obs_candidate_group])];
  int<lower = 1, upper = sum(candidate_group_size[obs_candidate_group]) + 1> obs_candidate_group_csr_row_pos[num_obs + 1];

  int<lower = 1, upper = num_r_types_full * num_unique_entities> entity_abducted_prob_ids[sum(abducted_prob_size) * num_unique_entities];
  int<lower = 1> entity_abducted_prob_csr_row_pos[num_abducted_estimands ? num_unique_entities * num_abducted_estimands + 1 : 0];

  int<lower = 1> long_entity_abducted_index[num_unique_entities * num_abducted_estimands];

  int<lower = 1, upper = num_r_types_full * num_unique_entities> entity_est_prob_ids[sum(est_prob_size) * num_unique_entities];
  int<lower = 1> entity_est_prob_csr_row_pos[num_atom_estimands > 0 ? num_unique_entities * num_atom_estimands + 1 : 0];

  int<lower = 1, upper = num_atom_estimands * num_unique_entities> entity_diff_estimand_ids[num_diff_estimands * 2 * num_unique_entities];
  int<lower = 1> entity_diff_estimand_csr_row_pos[num_diff_estimands > 0 ? num_unique_entities * num_diff_estimands + 1 : 0];

  int<lower = 1, upper = num_discretized_groups * num_unique_entities> entity_mean_diff_estimand_ids[num_mean_diff_estimands * 2 * num_unique_entities];
  int<lower = 1> entity_mean_diff_estimand_csr_row_pos[num_mean_diff_estimands > 0 ? num_unique_entities * num_mean_diff_estimands + 1 : 0];

  int<lower = 1, upper = num_discretized_groups * num_unique_entities> entity_utility_diff_estimand_ids[num_utility_diff_estimands * 2 * num_unique_entities];
  int<lower = 1> entity_utility_diff_estimand_csr_row_pos[num_utility_diff_estimands > 0 ? num_unique_entities * num_utility_diff_estimands + 1 : 0];

  int<lower = 1> total_num_bg_variable_types = sum(num_bg_variable_types);

  vector<lower = 0, upper = 1>[sum(num_bg_variable_type_combo_members) * num_unique_entities] marginal_prob_csr_vec;
  int<lower = 1, upper = num_r_types * num_unique_entities> entity_marginal_prob_ids[sum(num_bg_variable_type_combo_members) * num_unique_entities];
  int<lower = 1> entity_marginal_prob_csr_row_pos[total_num_bg_variable_types + 1];

  vector<lower = 0, upper = 1>[num_all_estimands * num_unique_entities * num_estimand_levels] level_estimands_csr_vec;
  int<lower = 1, upper = num_all_estimands * num_unique_entities> entity_estimand_ids[num_all_estimands * num_unique_entities * num_estimand_levels];
  int<lower = 1> entity_estimand_csr_row_pos[num_estimand_levels > 0 ? num_all_estimands * sum(level_size[estimand_levels]) + 1 : 0];

  vector<lower = -1, upper = 1>[num_discretized_groups * (2 * num_discretized_variables + 1) * num_unique_entities] entity_histogram_vec;
  int<lower = 1, upper = num_atom_estimands * num_unique_entities + 1> entity_histogram_ids[num_discretized_groups * (2 * num_discretized_variables + 1) * num_unique_entities];
  int<lower = 1> entity_histogram_csr_row_pos[num_discretized_groups > 0 ? num_discretized_groups * (num_cutpoints - 1) * num_unique_entities + 1 : 0];

  vector<lower = min(cutpoint_midpoints), upper = max(cutpoint_midpoints)>[num_discretized_groups * (num_cutpoints - 1) * num_unique_entities] entity_midpoint_vec;
  int<lower = 1, upper = num_discretized_groups * (num_cutpoints - 1) * num_unique_entities> entity_midpoint_ids[num_discretized_groups * (num_cutpoints - 1) * num_unique_entities];
  int<lower = 1> entity_midpoint_csr_row_pos[num_discretized_groups > 0 ? num_discretized_groups * num_unique_entities + 1 : 0];

  vector<lower = min(utility), upper = max(utility)>[num_discrete_utility_values > 0 ? num_discretized_groups * (num_cutpoints - 1) * num_unique_entities : 0] entity_utility_vec;

  vector<lower = -1, upper = 1>[2 * num_atom_estimands * (sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels)] between_entity_diff_csr_vec =
    to_vector(rep_matrix([-1, 1]', num_atom_estimands * (sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels))); // [-1, 1, -1, 1, ...]'
  int<lower = 1, upper = num_all_estimands * sum(level_size)> between_entity_diff_csr_ids[2 * num_atom_estimands * (sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels)];
  int<lower = 1> between_entity_diff_csr_row_pos[num_between_entity_diff_levels > 0 ? num_atom_estimands * (sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels) + 1 : 0];

  int<lower = 2> nonzero_beta_offsets[max(num_discretized_r_types - 1, 0) * num_discrete_r_types];

  if (num_discretized_r_types > 0) {
    int index_pos = 1;
    for (offset_index in 2:(num_discretized_r_types * num_discrete_r_types)) {
      if (offset_index % num_discretized_r_types != 1) {
        nonzero_beta_offsets[index_pos] = offset_index;
        index_pos += 1;
      }
    }
  }

  // TODO calculate the expected among and check it
  // if (num_discretized_r_types > 0 && num_r_types != num_discrete_r_types * num_discretized_r_types) {
  //   reject("Error in specified r_type sizes.")
  // }

  for (cutpoint_index in 2:num_cutpoints) {
    cutpoint_midpoints[cutpoint_index - 1] = mean(cutpoints[(cutpoint_index - 1):cutpoint_index]);

    discretized_binwidth[cutpoint_index - 1] = cutpoints[cutpoint_index] - cutpoints[cutpoint_index - 1];

    discretize_bin_alpha[cutpoint_index - 1] = cutpoints[cutpoint_index - 1];
    discretize_bin_beta[cutpoint_index - 1] = cutpoints[cutpoint_index];
  }

  if (log_lik_level > 0) {
    num_obs_in_log_lik_level_entities = calculate_num_in_level_entities(unique_entity_ids[obs_unique_entity_id, log_lik_level:log_lik_level], { level_size[log_lik_level] });
    obs_in_log_lik_level_entities = sort_indices_asc(unique_entity_ids[obs_unique_entity_id, log_lik_level]);
  }

  for (level_index in 1:num_levels) {
    unique_entities_in_level_entities[, level_index] = sort_indices_asc(unique_entity_ids[, level_index]);
  }

  // entity_prob_ids = csr_shift_expand_v(seq(1, num_r_types, 1), num_r_types, num_unique_entities);
  // entity_prob_csr_row_pos = csr_shift_expand_u( { num_r_types }, num_unique_entities);
  //
  // prob_repeater_ids = csr_shift_expand_v(rep_array(1, num_r_types), 1, num_unique_entities);
  // prob_repeater_csr_row_pos =  csr_shift_expand_u(rep_array(1, num_r_types), num_unique_entities);

  if (num_abducted_estimands > 0) {
    entity_abducted_prob_ids = csr_shift_expand_v(abducted_prob_index, num_r_types_full, num_unique_entities);
    entity_abducted_prob_csr_row_pos = csr_shift_expand_u(abducted_prob_size, num_unique_entities);
  }

  if (num_atom_estimands > 0) {
    entity_est_prob_ids = csr_shift_expand_v(est_prob_index, num_r_types_full, num_unique_entities);
    entity_est_prob_csr_row_pos = csr_shift_expand_u(est_prob_size, num_unique_entities);
  }

  if (num_diff_estimands > 0) {
    entity_diff_estimand_ids = csr_shift_expand_v(diff_estimand_atoms, num_atom_estimands, num_unique_entities);
    entity_diff_estimand_csr_row_pos = csr_shift_expand_u(rep_array(2, num_diff_estimands), num_unique_entities);
  }

  if (num_mean_diff_estimands > 0) {
    entity_mean_diff_estimand_ids = csr_shift_expand_v(mean_diff_estimand_atoms, num_discretized_groups, num_unique_entities);
    entity_mean_diff_estimand_csr_row_pos = csr_shift_expand_u(rep_array(2, num_mean_diff_estimands), num_unique_entities);
  }

  if (num_utility_diff_estimands > 0) {
    entity_utility_diff_estimand_ids = csr_shift_expand_v(utility_diff_estimand_atoms, num_discretized_groups, num_unique_entities);
    entity_utility_diff_estimand_csr_row_pos = csr_shift_expand_u(rep_array(2, num_utility_diff_estimands), num_unique_entities);
  }

  if (num_obs > 0) {
    {
      int entity_candidate_group_pos = 1;
      int entity_candidate_pos = 1;
      int entity_candidate_group_csr_row_pos_pos = 1;
      int last_candidate_group_size;

      int long_entity_abducted_index_pos = 1;

      for (entity_index in 1:num_unique_entities) {
        int num_entity_candidate_groups = num_unique_entity_candidate_groups[entity_index];
        int entity_candidate_group_end = entity_candidate_group_pos + num_entity_candidate_groups - 1;

        int curr_candidate_groups[num_entity_candidate_groups] = unique_entity_candidate_groups[entity_candidate_group_pos:entity_candidate_group_end];

        num_obs_in_unique_entity[entity_index] = num_equals(obs_unique_entity_id, { entity_index });

        // entity_total_num_candidates[entity_index] = sum(candidate_group_size[curr_candidate_groups]);

        for (candidate_group_index in 1:num_entity_candidate_groups) {
          int curr_candidate_group = curr_candidate_groups[candidate_group_index];
          int candidate_pos = curr_candidate_group > 1 ? sum(candidate_group_size[1:(curr_candidate_group - 1)]) + 1 : 1;
          int candidate_end = sum(candidate_group_size[1:curr_candidate_group]);

          int entity_candidate_end = entity_candidate_pos + candidate_group_size[curr_candidate_group] - 1;

          entity_candidate_group_ids[entity_candidate_pos:entity_candidate_end] = array_add(candidate_group_ids[candidate_pos:candidate_end], { (entity_index - 1) * num_r_types });

          entity_candidate_group_csr_row_pos[entity_candidate_group_csr_row_pos_pos] =
            entity_candidate_group_csr_row_pos_pos > 1 ?
            entity_candidate_group_csr_row_pos[entity_candidate_group_csr_row_pos_pos - 1] + last_candidate_group_size : 1;

          last_candidate_group_size = candidate_group_size[curr_candidate_groups][candidate_group_index];

          entity_candidate_pos = entity_candidate_end + 1;
          entity_candidate_group_csr_row_pos_pos += 1;
        }

        entity_candidate_group_csr_row_pos[entity_candidate_group_csr_row_pos_pos] = entity_candidate_group_csr_row_pos[entity_candidate_group_csr_row_pos_pos - 1] + last_candidate_group_size;

        entity_candidate_group_pos = entity_candidate_group_end + 1;

        if (num_abducted_estimands > 0) {
          int long_entity_abducted_index_end = long_entity_abducted_index_pos + num_abducted_estimands - 1;

          long_entity_abducted_index[long_entity_abducted_index_pos:long_entity_abducted_index_end] = array_add(abducted_estimand_ids, { (entity_index - 1) * num_atom_estimands });

          long_entity_abducted_index_pos = long_entity_abducted_index_end + 1;
        }
      }
    }

    {
      int obs_candidate_pos = 1;
      int last_candidate_group_size;

      for (obs_index in 1:num_obs) {
        int entity_index = obs_unique_entity_id[obs_index];
        int obs_candidate_end = obs_candidate_pos + candidate_group_size[obs_candidate_group[obs_index]] - 1;

        int curr_candidate_group = obs_candidate_group[obs_index];

        int candidate_pos = curr_candidate_group > 1 ? sum(candidate_group_size[1:(curr_candidate_group - 1)]) + 1 : 1;
        int candidate_end = sum(candidate_group_size[1:curr_candidate_group]);

        obs_candidate_group_ids[obs_candidate_pos:obs_candidate_end] = array_add(candidate_group_ids[candidate_pos:candidate_end], { (entity_index - 1) * num_r_types });

        obs_candidate_group_csr_row_pos[obs_index] = obs_index > 1 ? obs_candidate_group_csr_row_pos[obs_index - 1] + last_candidate_group_size : 1;

        last_candidate_group_size = candidate_group_size[curr_candidate_group];

        obs_candidate_pos = obs_candidate_end + 1;
      }

      obs_candidate_group_csr_row_pos[num_obs + 1] = obs_candidate_group_csr_row_pos[num_obs] + last_candidate_group_size;

      unique_entity_prop = num_obs_in_unique_entity / num_obs;
    }
  } else {
    obs_candidate_group_csr_row_pos[1] = 1;

    unique_entity_prop = rep_vector(1.0 / num_unique_entities, num_unique_entities);

    entity_candidate_group_csr_row_pos[1] = 1;

    num_obs_in_unique_entity = rep_vector(0, num_unique_entities);
  }


  {
    int latent_type_marginal_members_pos = 1;
    int entity_marginal_prob_pos = 1;
    int entity_marginal_prob_csr_row_pos_pos = 1;
    int last_entity_marginal_prob_size;

    for (latent_type_index in 1:total_num_bg_variable_types) {
      int curr_var_size = num_bg_variable_type_combo_members[latent_type_index];
      int latent_type_marginal_members_end = latent_type_marginal_members_pos + curr_var_size - 1;

      entity_marginal_prob_csr_row_pos[entity_marginal_prob_csr_row_pos_pos] =
        entity_marginal_prob_csr_row_pos_pos > 1 ?
        entity_marginal_prob_csr_row_pos[entity_marginal_prob_csr_row_pos_pos - 1] + last_entity_marginal_prob_size : 1;

      last_entity_marginal_prob_size = curr_var_size * num_unique_entities;

      for (entity_index in 1:num_unique_entities) {
        real entity_prop = num_obs_in_unique_entity[entity_index] / max(num_obs, 1);
        int entity_marginal_prob_end = entity_marginal_prob_pos + curr_var_size - 1;

        marginal_prob_csr_vec[entity_marginal_prob_pos:entity_marginal_prob_end] = rep_vector(entity_prop, curr_var_size);

        entity_marginal_prob_ids[entity_marginal_prob_pos:entity_marginal_prob_end] =
          array_add(bg_variable_type_combo_members[latent_type_marginal_members_pos:latent_type_marginal_members_end], { (entity_index - 1) * num_r_types });

        entity_marginal_prob_pos = entity_marginal_prob_end + 1;
      }

      latent_type_marginal_members_pos = latent_type_marginal_members_end + 1;
      entity_marginal_prob_csr_row_pos_pos += 1;
    }

    entity_marginal_prob_csr_row_pos[entity_marginal_prob_csr_row_pos_pos] = entity_marginal_prob_csr_row_pos[entity_marginal_prob_csr_row_pos_pos - 1] + last_entity_marginal_prob_size;
  }

  if (num_estimand_levels > 0) {
    int estimand_entity_pos = 1;
    int level_estimands_csr_vec_pos = 1;
    int entity_estimand_csr_row_pos_pos = 1;
    int last_num_unique_entities;

    for (estimand_level_index in 1:num_estimand_levels) {
      int curr_estimand_level = estimand_levels[estimand_level_index];
      int curr_estimand_level_size = level_size[curr_estimand_level];
      int curr_num_unique_entities[curr_estimand_level_size] = array_extract_group_values(num_unique_entities_in_estimand_level_entities, level_size[estimand_levels], { estimand_level_index });

      for (estimand_level_entity_index in 1:curr_estimand_level_size) {
        int unique_entities_in_curr_level_entity[curr_num_unique_entities[estimand_level_entity_index]] =
          array_extract_group_values(unique_entities_in_level_entities[, curr_estimand_level], curr_num_unique_entities, { estimand_level_entity_index });

        int unique_entities_estimand_offset[curr_num_unique_entities[estimand_level_entity_index]] =
          array_product(array_add(unique_entities_in_curr_level_entity, { -1 }), { num_all_estimands });

        vector[curr_num_unique_entities[estimand_level_entity_index]] curr_num_obs_in_unique_entity = num_obs_in_unique_entity[unique_entities_in_curr_level_entity];

        for (estimand_index in 1:num_all_estimands) {
          int level_estimands_csr_vec_end = level_estimands_csr_vec_pos + curr_num_unique_entities[estimand_level_entity_index] - 1;

          entity_estimand_ids[level_estimands_csr_vec_pos:level_estimands_csr_vec_end] = array_add(unique_entities_estimand_offset, { estimand_index });
          level_estimands_csr_vec[level_estimands_csr_vec_pos:level_estimands_csr_vec_end] = curr_num_obs_in_unique_entity / sum(curr_num_obs_in_unique_entity);

          entity_estimand_csr_row_pos[entity_estimand_csr_row_pos_pos] =
            entity_estimand_csr_row_pos_pos > 1 ?
            entity_estimand_csr_row_pos[entity_estimand_csr_row_pos_pos - 1] + last_num_unique_entities : 1;

          last_num_unique_entities = curr_num_unique_entities[estimand_level_entity_index];

          entity_estimand_csr_row_pos_pos += 1;
          level_estimands_csr_vec_pos = level_estimands_csr_vec_end + 1;
        }

        estimand_entity_pos += 1;
      }
    }

    entity_estimand_csr_row_pos[entity_estimand_csr_row_pos_pos] = entity_estimand_csr_row_pos[entity_estimand_csr_row_pos_pos - 1] + last_num_unique_entities;
  }

  if (num_discretized_groups > 0) {
    int discretized_group_pos = 1;
    int histogram_vec_pos = 1;
    int histogram_ids_pos = 1;

    int midpoint_vec_pos = 1;
    int midpoint_ids_pos = 1;

    int histogram_ids_size = num_discretized_groups * (2 * num_discretized_variables + 1);
    int histogram_last_id_pos[num_discretized_groups];
    int last_histogram_one_pos = (num_unique_entities * num_atom_estimands) + 1;

    int histogram_ids[histogram_ids_size];

    int midpoint_ids[num_discretized_groups * (num_cutpoints - 1)];

    // 1,2,2,2,2,....
    int histogram_row_sizes[(num_cutpoints - 1) * num_discretized_groups] = to_array_1d(rep_array(append_array({ 1 }, rep_array(2, num_discretized_variables)), num_discretized_groups));

    entity_histogram_csr_row_pos = csr_shift_expand_u(histogram_row_sizes, num_unique_entities);
    entity_midpoint_csr_row_pos = csr_shift_expand_u(rep_array(num_cutpoints - 1, num_discretized_groups), num_unique_entities);

    for (discretized_group_index in 1:num_discretized_groups) {
      int histogram_ids_end = histogram_ids_pos + (2 * num_discretized_variables + 1) - 1;
      int discretized_group_end = discretized_group_pos + num_discretized_variables - 1;
      int midpoint_ids_end = midpoint_ids_pos + (num_cutpoints - 1) - 1;

      int curr_group_members[num_discretized_variables] = discretized_group_ids[discretized_group_pos:discretized_group_end];

      histogram_last_id_pos[discretized_group_index] = histogram_ids_end;

      for (entity_index in 1:num_unique_entities) {
        int histogram_vec_end = histogram_vec_pos + (2 * num_discretized_variables + 1) - 1;
        int midpoint_vec_end = midpoint_vec_pos + (num_cutpoints - 1) - 1;

        entity_histogram_vec[histogram_vec_pos:histogram_vec_end] = append_row(to_vector(rep_matrix([1, -1]', num_discretized_variables)), 1);
        entity_midpoint_vec[midpoint_vec_pos:midpoint_vec_end] = cutpoint_midpoints;

        if (num_discrete_utility_values > 0) {
          entity_utility_vec[midpoint_vec_pos:midpoint_vec_end] = utility;
        }

        entity_discretize_bin_alpha[midpoint_vec_pos:midpoint_vec_end] = discretize_bin_alpha;
        entity_discretize_bin_beta[midpoint_vec_pos:midpoint_vec_end] = discretize_bin_beta;

        histogram_vec_pos = histogram_vec_end + 1;
        midpoint_vec_pos = midpoint_vec_end + 1;
      }

      for (group_member_index in 1:num_discretized_variables) {
        int offset = 2 * (group_member_index - 1);

        histogram_ids[(histogram_ids_pos + offset):(histogram_ids_pos + offset + 1)] = rep_array(curr_group_members[group_member_index], 2);
      }

      histogram_ids[(histogram_ids_end - 1):histogram_ids_end] = { curr_group_members[num_discretized_variables], 0 };

      midpoint_ids[midpoint_ids_pos:midpoint_ids_end] = seq(1 + (discretized_group_index - 1) * (num_cutpoints - 1), discretized_group_index * (num_cutpoints - 1), 1);

      histogram_ids_pos = histogram_ids_end + 1;
      discretized_group_pos = discretized_group_end + 1;
      midpoint_ids_pos = midpoint_ids_end + 1;
    }

    entity_histogram_ids = csr_shift_expand_v(histogram_ids, num_atom_estimands, num_unique_entities);

    for (entity_index in 1:num_unique_entities) {
      entity_histogram_ids[array_add(histogram_last_id_pos, { (entity_index - 1) * histogram_ids_size })] = rep_array(last_histogram_one_pos, num_discretized_groups);
    }

    entity_midpoint_ids = csr_shift_expand_v(midpoint_ids, num_discretized_groups * (num_cutpoints - 1), num_unique_entities);
  }

  if (num_between_entity_diff_levels > 0) {
    int between_csr_ids_pos = 1;

    between_entity_diff_csr_row_pos = seq(1, 2 * num_atom_estimands * (sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels) + 1, 2);

    for (between_diff_level_index_index in 1:num_between_entity_diff_levels) {
      int between_diff_level_index = between_entity_diff_levels[between_diff_level_index_index];
      int between_diff_level_offset = between_diff_level_index > 1 ? sum(level_size[:(between_diff_level_index - 1)]) * num_all_estimands : 0;

      for (between_diff_entity_index in 1:(level_size[between_diff_level_index] - 1)) {
        int diff_pair[2] = { between_diff_level_offset, between_diff_level_offset + num_all_estimands * between_diff_entity_index };

        for (estimand_index in 1:num_atom_estimands) {
          between_entity_diff_csr_ids[between_csr_ids_pos:(between_csr_ids_pos + 1)] = array_add(diff_pair, { estimand_index });
          between_csr_ids_pos += 2;
        }
      }
    }
  }

  {
    int exp_r_pos = 1;

    for (r_type_index in 1:num_r_types) {
      int exp_r_end = exp_r_pos + num_experiment_types - 1;

      experiment_r_type_index[exp_r_pos:exp_r_end] = rep_array(r_type_index, num_experiment_types);

      exp_r_pos = exp_r_end + 1;
    }
  }
}

parameters {
  vector[num_discrete_r_types - 1] toplevel_discrete_beta;

  matrix[num_discrete_r_types - 1, sum(level_size)] discrete_level_beta_raw;
  matrix<lower = 0>[num_discrete_r_types - 1, num_levels] discrete_level_beta_sigma;

  matrix[max(0, num_discretized_r_types - 1), num_discrete_r_types] toplevel_discretized_beta[num_discretized_variables];

  matrix[max(num_discretized_r_types - 1, 0) * num_discrete_r_types, sum(level_size)] discretized_level_beta_raw[num_discretized_variables];
  matrix<lower = 0>[max(num_discretized_r_types - 1, 0), num_levels] discretized_level_beta_sigma[num_discretized_variables];
}

transformed parameters {
  vector<lower = 0, upper = 1>[num_r_types * num_unique_entities] r_prob; // Entities first, Types second

  vector[run_type == RUN_TYPE_FIT ? sum(num_unique_entity_candidate_groups) : 0] entity_candidates_group_logp;

  matrix[num_discrete_r_types, num_unique_entities] discrete_beta = rep_matrix(append_row(0, toplevel_discrete_beta), num_unique_entities);

  matrix[num_discretized_r_types * num_discrete_r_types, num_unique_entities] discretized_beta[num_discretized_variables];

  for (discretized_var_index in 1:num_discretized_variables) {
    matrix[num_discretized_r_types, num_discrete_r_types] curr_discretized_beta =
      append_row(rep_row_vector(0, num_discrete_r_types), toplevel_discretized_beta[discretized_var_index]);

    discretized_beta[discretized_var_index] = rep_matrix(to_vector(curr_discretized_beta), num_unique_entities);
  }

  if (num_levels > 0) {
    int level_entity_pos = 1;

    for (level_index in 1:num_levels) {
      int level_entity_end = level_entity_pos + level_size[level_index] - 1;

      matrix[num_discrete_r_types - 1, level_size[level_index]] curr_discrete_level_beta =
        discrete_level_beta_raw[, level_entity_pos:level_entity_end] .* rep_matrix(discrete_level_beta_sigma[, level_index], level_size[level_index]);

      discrete_beta[2:] += curr_discrete_level_beta[, unique_entity_ids[, level_index]];

      for (discretized_var_index in 1:num_discretized_variables) {
        matrix[max(num_discretized_r_types - 1, 0) * num_discrete_r_types, level_size[level_index]] curr_discretized_level_beta =
          discretized_level_beta_raw[discretized_var_index, , level_entity_pos:level_entity_end] .*
          rep_matrix(to_vector(rep_matrix(discretized_level_beta_sigma[discretized_var_index, , level_index], num_discrete_r_types)), level_size[level_index]);

        discretized_beta[discretized_var_index, nonzero_beta_offsets] += curr_discretized_level_beta[, unique_entity_ids[, level_index]];
      }

      level_entity_pos = level_entity_end + 1;
    }
  }

  for (entity_index in 1:num_unique_entities) {
    vector[num_discrete_r_types] curr_discrete_prob = softmax(discrete_beta[, entity_index]);

    if (num_discretized_r_types > 0) {
      int r_prob_pos = (entity_index - 1) * num_r_types + 1;
      int r_prob_end = entity_index * num_r_types;

      r_prob[r_prob_pos:r_prob_end] = calculate_r_type_joint_prob(num_r_types, num_discrete_r_types, num_discretized_r_types, discrete_group_size,
                                                                  num_compatible_discretized_r_types, compatible_discretized_r_types, compatible_discretized_pair_ids,
                                                                  curr_discrete_prob, discretized_beta, entity_index);
      // for (discrete_index in 1:num_discrete_r_types) {
      //   int discretized_pos = 1 + (discrete_index - 1) * num_discretized_r_types;
      //   int discretized_end = discrete_index * num_discretized_r_types;
      //   int discrete_group_pos = 1 + (discrete_index - 1) * discrete_group_size;
      //   int discrete_group_end = discrete_index * discrete_group_size;
      //   int r_prob_pos = (entity_index - 1) * num_r_types + discrete_group_pos;
      //   int r_prob_end = (entity_index - 1) * num_r_types + discrete_group_end;
      //
      //   vector[num_discretized_r_types] first_discretized_prob = softmax(discretized_beta[1, discretized_pos:discretized_end, entity_index]);
      //
      //   r_prob[r_prob_pos:r_prob_end] = curr_discrete_prob[discrete_index] * first_discretized_prob[compatible_discretized_pair_ids[1, discrete_group_pos:discrete_group_end]];
      //
      //   for (discretized_var_index in 2:num_discretized_variables) {
      //     vector[sum(num_compatible_discretized_r_types)] curr_cond_discretized_prob;
      //     int cond_prob_pos = 1;
      //
      //     for (discretized_type_index in 1:num_discretized_r_types) {
      //       int cond_prob_end = cond_prob_pos + num_compatible_discretized_r_types[discretized_type_index] - 1;
      //       int compatible_ids[num_compatible_discretized_r_types[discretized_type_index]] = compatible_discretized_r_types[cond_prob_pos:cond_prob_end, 1];
      //
      //       curr_cond_discretized_prob[cond_prob_pos:cond_prob_end] =
      //         softmax(discretized_beta[discretized_var_index, discretized_pos:discretized_end, entity_index][compatible_ids]);
      //
      //       cond_prob_pos = cond_prob_end + 1;
      //     }
      //
      //     r_prob[r_prob_pos:r_prob_end] =
      //       r_prob[r_prob_pos:r_prob_end] .* curr_cond_discretized_prob[compatible_discretized_pair_ids[discretized_var_index, discrete_group_pos:discrete_group_end]];
      //   }
      // }
    } else {
      int r_prob_pos = (entity_index - 1) * num_r_types + 1;
      int r_prob_end = (entity_index) * num_r_types;

      r_prob[r_prob_pos:r_prob_end] = curr_discrete_prob[discrete_r_type_id];
    }
  }
  // } else {
  //   vector[num_discrete_r_types] curr_discrete_prob = softmax(discrete_beta[, 1]);
  //
  //   if (num_discretized_r_types > 0) {
  //     r_prob = calculate_r_type_joint_prob(num_r_types, num_discrete_r_types, num_discretized_r_types, discrete_group_size,
  //                                          num_compatible_discretized_r_types, compatible_discretized_r_types, compatible_discretized_pair_ids,
  //                                          curr_discrete_prob, discretized_beta, 1);
  //     // for (discrete_index in 1:num_discrete_r_types) {
  //     //   int discretized_pos = 1 + (discrete_index - 1) * num_discretized_r_types;
  //     //   int discretized_end = discrete_index * num_discretized_r_types;
  //     //   int discrete_group_pos = 1 + (discrete_index - 1) * discrete_group_size;
  //     //   int discrete_group_end = discrete_index * discrete_group_size;
  //     //   int r_prob_pos = discrete_group_pos;
  //     //   int r_prob_end = discrete_group_end;
  //     //
  //     //   vector[num_discretized_r_types] first_discretized_prob = softmax(discretized_beta[1, discrete_index, discretized_pos:discretized_end]');
  //     //
  //     //   r_prob[r_prob_pos:r_prob_end] = discrete_prob[discrete_index] * first_discretized_prob[compatible_discretized_pair_ids[1, discrete_group_pos:discrete_group_end]];
  //     //
  //     //   for (discretized_var_index in 2:num_discretized_variables) {
  //     //     vector[sum(num_compatible_discretized_r_types)] curr_cond_discretized_prob;
  //     //     int cond_prob_pos = 1;
  //     //
  //     //     for (discretized_type_index in 1:num_discretized_r_types) {
  //     //       int cond_prob_end = cond_prob_pos + num_compatible_discretized_r_types[discretized_type_index] - 1;
  //     //       int compatible_ids[num_compatible_discretized_r_types[discretized_type_index]] = compatible_discretized_r_types[cond_prob_pos:cond_prob_end, 1];
  //     //
  //     //       curr_cond_discretized_prob[cond_prob_pos:cond_prob_end] =
  //     //         softmax(discretized_beta[discretized_var_index, discrete_index, discretized_pos:discretized_end][compatible_ids]');
  //     //
  //     //       cond_prob_pos = cond_prob_end + 1;
  //     //     }
  //     //
  //     //     r_prob[r_prob_pos:r_prob_end] =
  //     //       r_prob[r_prob_pos:r_prob_end] .* curr_cond_discretized_prob[compatible_discretized_pair_ids[discretized_var_index, discrete_group_pos:discrete_group_end]];
  //     //   }
  //     // }
  //   } else {
  //     r_prob = curr_discrete_prob;
  //   }
  // }

  if (run_type == RUN_TYPE_FIT) {
    entity_candidates_group_logp = log(csr_matrix_times_vector(sum(num_unique_entity_candidate_groups),
                                                               num_r_types * num_unique_entities,
                                                               vec_1[1:sum(candidate_group_size[unique_entity_candidate_groups])],
                                                               entity_candidate_group_ids,
                                                               entity_candidate_group_csr_row_pos,
                                                               r_prob));
  }
}

model {
  toplevel_discrete_beta ~ normal(0, discrete_beta_hyper_sd);

  if (num_levels > 0) {
    to_vector(discrete_level_beta_raw) ~ std_normal();

    for (level_index in 1:num_levels) {
      discrete_level_beta_sigma[, level_index] ~ normal(0, tau_level_sigma[level_index]);
    }
  }

  for (discretized_var_index in 1:num_discretized_variables) {
    to_vector(toplevel_discretized_beta[discretized_var_index]) ~ normal(0, discretized_beta_hyper_sd);

    if (num_levels > 0) {
      to_vector(discretized_level_beta_raw[discretized_var_index]) ~ std_normal();

      for (level_index in 1:num_levels) {
        discretized_level_beta_sigma[discretized_var_index, , level_index] ~ normal(0, tau_level_sigma[level_index]);
      }
    }
  }

  if (run_type == RUN_TYPE_FIT) {
    target += num_unique_entity_in_candidate_groups * entity_candidates_group_logp;
  }
}

generated quantities {
  vector[run_type == RUN_TYPE_FIT && log_lik_level > -1 ? (log_lik_level > 0 ? level_size[log_lik_level] : num_obs) : 0] log_lik;

  vector<lower = 0, upper = 1>[calculate_marginal_prob ? total_num_bg_variable_types : 0] marginal_p_r;
  // matrix<lower = 0, upper = 1>[calculate_marginal_prob ? total_num_bg_variable_types : 0, sum(level_size)] level_marginal_p_r;

  vector<lower = 0, upper = 1.01>[num_abducted_estimands * num_unique_entities] total_abducted_prob;

  matrix<lower = 0, upper = 1.01>[num_atom_estimands, num_unique_entities] iter_atom_estimand;
  matrix<lower = -1.01, upper = 1.01>[num_diff_estimands, num_unique_entities] iter_diff_estimand;

  matrix[num_all_estimands, num_unique_entities] iter_entity_estimand;
  vector[num_all_estimands] iter_estimand;
  matrix[num_all_estimands, sum(level_size[estimand_levels])] iter_level_entity_estimand;
  matrix[num_all_estimands, run_type == RUN_TYPE_PRIOR_PREDICT ? num_estimand_levels : 0] iter_level_entity_estimand_sd;
  matrix[num_atom_estimands, sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels] iter_between_level_entity_diff_estimand;

  vector<lower = 0, upper = 1>[num_discretized_groups * (num_cutpoints - 1) * num_unique_entities] iter_entity_discretized_histogram_vec;
  matrix[num_discretized_groups, num_unique_entities] iter_entity_discretized_mean;
  matrix<lower = min(utility), upper = max(utility)>[num_discrete_utility_values > 0 ? num_discretized_groups : 0, num_unique_entities] iter_entity_discretized_utility;
  matrix[num_mean_diff_estimands, num_unique_entities] iter_mean_diff_estimand;
  matrix[num_mean_diff_estimands, num_unique_entities] iter_utility_diff_estimand;

  // matrix<lower = -1, upper = 1>[generate_rep && num_rep_corr_estimands > 0 && run_type == RUN_TYPE_FIT ? num_rep_corr_estimands : 0, 2] rep_corr_estimand;

  if (run_type == RUN_TYPE_FIT && log_lik_level > -1) {
    log_lik = log(csr_matrix_times_vector(num_obs,
                                          num_r_types * num_unique_entities,
                                          vec_1[1:sum(candidate_group_size[obs_candidate_group])],
                                          obs_candidate_group_ids,
                                          obs_candidate_group_csr_row_pos,
                                          r_prob));

    // for (obs_index in 1:num_obs) {
    // }
  }

  // if (log_lik_level > 0) {
  //   int level_entity_obs_pos = 1;
  //
  //   for (level_entity_index in 1:level_size[log_lik_level]) {
  //     int num_obs_in_level_entity = num_obs_in_log_lik_level_entities[level_entity_index];
  //     int level_entity_obs_end = level_entity_obs_pos + num_obs_in_level_entity - 1;
  //
  //     log_lik[level_entity_index] = sum(obs_log_lik[obs_in_log_lik_level_entities[level_entity_obs_pos:level_entity_obs_end]]);
  //
  //     level_entity_obs_pos = level_entity_obs_end + 1;
  //   }
  // } else if (log_lik_level == 0) {
  //   // log_lik = obs_log_lik;
  // }

    // for (obs_index in 1:num_obs) {
    //   int candidate_end = candidate_pos + num_r_candidates[obs_index] - 1;
    //
    //   int curr_r_candidates[num_r_candidates[obs_index]] = r_candidates[candidate_pos:candidate_end];
    //   vector[num_r_candidates[obs_index]] curr_candidate_r_prob = r_prob[curr_r_candidates, obs_unique_entity_id[obs_index]];
    //
    //   obs_log_lik[obs_index] = log(sum(curr_candidate_r_prob)); // + curr_cont_outcome_log_lik;
    //
    //   candidate_pos = candidate_end + 1;
    // }

  if (num_discrete_estimands > 0) {
    matrix[num_experiment_types, num_r_types * num_unique_entities] full_r_prob = rep_matrix(r_prob', num_experiment_types) .* rep_matrix(experiment_types_prob, num_r_types * num_unique_entities);
    vector[num_r_types_full * num_unique_entities] full_r_prob_vec = to_vector(full_r_prob);

    vector[num_atom_estimands * num_unique_entities] iter_atom_estimand_vec =
      csr_matrix_times_vector(num_atom_estimands * num_unique_entities,
                              num_r_types_full * num_unique_entities,
                              vec_1[1:(num_unique_entities * sum(est_prob_size))],
                              entity_est_prob_ids,
                              entity_est_prob_csr_row_pos,
                              full_r_prob_vec);

    vector[num_diff_estimands * num_unique_entities] iter_diff_estimand_vec;

    if (num_abducted_estimands > 0) {
      total_abducted_prob = csr_matrix_times_vector(
        num_abducted_estimands * num_unique_entities,
        num_r_types_full * num_unique_entities,
        vec_1[1:(num_unique_entities * sum(abducted_prob_size))],
        entity_abducted_prob_ids,
        entity_abducted_prob_csr_row_pos,
        full_r_prob_vec);

      iter_atom_estimand_vec[long_entity_abducted_index] = iter_atom_estimand_vec[long_entity_abducted_index] ./ total_abducted_prob;
    }

    iter_atom_estimand = to_matrix(iter_atom_estimand_vec, num_atom_estimands, num_unique_entities);

    if (num_diff_estimands > 0) {
      iter_diff_estimand_vec = csr_matrix_times_vector(
        num_diff_estimands * num_unique_entities,
        num_atom_estimands * num_unique_entities,
        vec_diff,
        entity_diff_estimand_ids,
        entity_diff_estimand_csr_row_pos,
        iter_atom_estimand_vec);

      iter_diff_estimand = to_matrix(iter_diff_estimand_vec, num_diff_estimands, num_unique_entities);
    }

    if (num_discretized_groups > 0) {
     vector[num_discretized_groups * num_unique_entities] iter_entity_discretized_mean_vec;

      iter_entity_discretized_histogram_vec =
        csr_matrix_times_vector(
          num_discretized_groups * (num_cutpoints - 1) * num_unique_entities,
          num_atom_estimands * num_unique_entities + 1,
          entity_histogram_vec,
          entity_histogram_ids,
          entity_histogram_csr_row_pos,
          append_row(iter_atom_estimand_vec, 1)
        );

      iter_entity_discretized_mean_vec = csr_matrix_times_vector(
        num_discretized_groups * num_unique_entities,
        num_discretized_groups * (num_cutpoints - 1) * num_unique_entities,
        use_random_binpoint ? to_vector(uniform_rng(entity_discretize_bin_alpha, entity_discretize_bin_beta)) : entity_midpoint_vec,
        entity_midpoint_ids,
        entity_midpoint_csr_row_pos,
        iter_entity_discretized_histogram_vec
      );

      if (num_discrete_utility_values > 0) {
        vector[num_discretized_groups * num_unique_entities] iter_entity_discretized_utility_vec =
          csr_matrix_times_vector(
            num_discretized_groups * num_unique_entities,
            num_discretized_groups * (num_cutpoints - 1) * num_unique_entities,
            entity_utility_vec,
            entity_midpoint_ids,
            entity_midpoint_csr_row_pos,
            iter_entity_discretized_histogram_vec
          );

        iter_entity_discretized_utility = to_matrix(iter_entity_discretized_utility_vec, num_discretized_groups, num_unique_entities);

        if (num_utility_diff_estimands > 0) {
          vector[num_utility_diff_estimands * num_unique_entities] iter_utility_diff_estimand_vec = csr_matrix_times_vector(
            num_utility_diff_estimands * num_unique_entities,
            num_discretized_groups * num_unique_entities,
            vec_utility_diff,
            entity_utility_diff_estimand_ids,
            entity_utility_diff_estimand_csr_row_pos,
            iter_entity_discretized_utility_vec);

          iter_utility_diff_estimand = to_matrix(iter_utility_diff_estimand_vec, num_utility_diff_estimands, num_unique_entities);
        }
     }

     iter_entity_discretized_mean = to_matrix(iter_entity_discretized_mean_vec, num_discretized_groups, num_unique_entities);

     if (num_mean_diff_estimands > 0) {
       vector[num_mean_diff_estimands * num_unique_entities] iter_mean_diff_estimand_vec = csr_matrix_times_vector(
         num_mean_diff_estimands * num_unique_entities,
         num_discretized_groups * num_unique_entities,
         vec_mean_diff,
         entity_mean_diff_estimand_ids,
         entity_mean_diff_estimand_csr_row_pos,
         iter_entity_discretized_mean_vec);

        iter_mean_diff_estimand = to_matrix(iter_mean_diff_estimand_vec, num_mean_diff_estimands, num_unique_entities);
     }

    }

    iter_entity_estimand =
      append_row(
        append_row(
          iter_atom_estimand,
          append_row(
            iter_diff_estimand,
            append_row(
              iter_entity_discretized_mean,
              iter_entity_discretized_utility
            )
          )
        ),
      append_row(iter_mean_diff_estimand, iter_utility_diff_estimand));

    iter_estimand = iter_entity_estimand * unique_entity_prop;

    if (num_estimand_levels > 0) {
      vector[num_all_estimands * sum(level_size[estimand_levels])] iter_level_entity_estimand_vec =
        csr_matrix_times_vector(
          num_all_estimands * sum(level_size[estimand_levels]),
          num_all_estimands * num_unique_entities,
          level_estimands_csr_vec,
          entity_estimand_ids,
          entity_estimand_csr_row_pos,
          to_vector(iter_entity_estimand)
        );

      iter_level_entity_estimand = to_matrix(iter_level_entity_estimand_vec, num_all_estimands, sum(level_size[estimand_levels]));

      if (run_type == RUN_TYPE_PRIOR_PREDICT) {
        int level_entity_pos = 1;

        for (est_level_index_index in 1:num_estimand_levels) {
          int est_level_index = estimand_levels[est_level_index_index];
          int level_entity_end = level_entity_pos + level_size[est_level_index] - 1;

          vector[num_all_estimands] level_entity_means = iter_level_entity_estimand[, level_entity_pos:level_entity_end] * rep_vector(1.0 / level_size[est_level_index], level_size[est_level_index]);

          iter_level_entity_estimand_sd[, est_level_index] =
            sqrt(
              square(iter_level_entity_estimand[, level_entity_pos:level_entity_end] - rep_matrix(level_entity_means, level_size[est_level_index]))
              * rep_vector(1.0 / (level_size[est_level_index] - 1), level_size[est_level_index])
            );

          level_entity_pos = level_entity_end + 1;
        }
      }

      if (num_between_entity_diff_levels > 0) {
        vector[num_atom_estimands * (sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels)] iter_between_level_entity_diff_estimand_vec =
          csr_matrix_times_vector(
            num_atom_estimands * (sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels),
            num_all_estimands * sum(level_size[estimand_levels]),
            between_entity_diff_csr_vec,
            between_entity_diff_csr_ids,
            between_entity_diff_csr_row_pos,
            iter_level_entity_estimand_vec
          );

        iter_between_level_entity_diff_estimand = to_matrix(iter_between_level_entity_diff_estimand_vec, num_atom_estimands, sum(level_size[between_entity_diff_levels]) - num_between_entity_diff_levels);
      }
    }
  }

  if (calculate_marginal_prob) {
    // vector[total_num_bg_variable_types * sum(level_size)] level_marginal_p_r_vec =
    //   csr_matrix_times_vector(total_num_bg_variable_types * sum(level_size),
    //                           num_r_types * num_unique_entities,
    //                           // marginal_prob_csr_vec,
    //                           // entity_marginal_prob_ids,
    //                           // entity_marginal_prob_csr_row_pos,
    //                           r_prob);
    //
    // level_marginal_p_r = to_matrix(level_marginal_p_r_vec);

    marginal_p_r = csr_matrix_times_vector(total_num_bg_variable_types,
                                           num_r_types * num_unique_entities,
                                           marginal_prob_csr_vec,
                                           entity_marginal_prob_ids,
                                           entity_marginal_prob_csr_row_pos,
                                           r_prob);
  }

  // if (generate_rep && num_rep_corr_estimands > 0 && run_type == RUN_TYPE_FIT) {
  //   matrix[num_r_types, num_unique_entities] r_prob_mat = to_matrix(r_prob, num_r_types, num_unique_entities);
  //
  //   rep_corr_estimand = rep_corr_estimand_rng(r_prob_mat, rep_corr_outcomes, rep_corr_cond, experiment_assign_entity, type_response_value);
  // }
}

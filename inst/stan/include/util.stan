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

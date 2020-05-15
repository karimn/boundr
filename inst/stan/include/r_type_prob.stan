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

vector calculate_r_type_joint_log_prob(int num_r_types, int num_discrete_r_types, int num_discretized_r_types, int discrete_group_size,
                                       int[] num_discretized_types_conditional, int[] discretized_types_conditional,
                                       int[] num_compatible_discretized_r_types, int[] compatible_discretized_r_types, int[,] compatible_discretized_pair_ids,
                                       vector discrete_log_prob, matrix[] beta, int entity_index) {
 int num_discretized_variables = size(beta);
 vector[num_r_types] joint_log_prob;

 int discretized_cond_pos = 1;

 for (discrete_index in 1:num_discrete_r_types) {
    int discretized_cond_end = discretized_cond_pos + num_discretized_types_conditional[discrete_index] - 1;

   // int discretized_pos = 1 + (discrete_index - 1) * num_discretized_r_types;
   // int discretized_end = discrete_index * num_discretized_r_types;

   int discrete_group_pos = 1 + (discrete_index - 1) * discrete_group_size;
   int discrete_group_end = discrete_index * discrete_group_size;
   int r_prob_pos = discrete_group_pos;
   int r_prob_end = discrete_group_end;

   // vector[num_discretized_r_types] first_discretized_log_prob = log_softmax(beta[1, discretized_pos:discretized_end, entity_index]);
   vector[num_discretized_r_types] first_discretized_log_prob = rep_vector(negative_infinity(), num_discretized_r_types);

   first_discretized_log_prob[discretized_types_conditional[discretized_cond_pos:discretized_cond_end]] =
      log_softmax(beta[1, discretized_cond_pos:discretized_cond_end, entity_index]);

   joint_log_prob[r_prob_pos:r_prob_end] = discrete_log_prob[discrete_index] + first_discretized_log_prob[compatible_discretized_pair_ids[1, discrete_group_pos:discrete_group_end]];

   for (discretized_var_index in 2:num_discretized_variables) {
     vector[sum(num_compatible_discretized_r_types)] curr_cond_discretized_log_prob;
     int cond_prob_pos = 1;

     vector[num_discretized_r_types] curr_beta = rep_vector(negative_infinity(), num_discretized_r_types);

     curr_beta[discretized_types_conditional[discretized_cond_pos:discretized_cond_end]] =
       beta[discretized_var_index, discretized_cond_pos:discretized_cond_end, entity_index];

     for (discretized_type_index in 1:num_discretized_r_types) {
       int cond_prob_end = cond_prob_pos + num_compatible_discretized_r_types[discretized_type_index] - 1;
       int compatible_ids[num_compatible_discretized_r_types[discretized_type_index]] = compatible_discretized_r_types[cond_prob_pos:cond_prob_end];

       // curr_cond_discretized_log_prob[cond_prob_pos:cond_prob_end] = log_softmax(beta[discretized_var_index, discretized_pos:discretized_end, entity_index][compatible_ids]);
       curr_cond_discretized_log_prob[cond_prob_pos:cond_prob_end] = log_softmax(curr_beta[compatible_ids]);

       cond_prob_pos = cond_prob_end + 1;
     }

     joint_log_prob[r_prob_pos:r_prob_end] += curr_cond_discretized_log_prob[compatible_discretized_pair_ids[discretized_var_index, discrete_group_pos:discrete_group_end]];
   }

    discretized_cond_pos = discretized_cond_end + 1;
 }

 return joint_log_prob;
}



data_dir = '/Users/alex/data_to_remove/aortic_13673972_384_9c6b2e5_true_bicuspid_nostab_circ_3pt9_rad_1pt7_new_initial_cond/field_data'

file_name = 'aortic_no_partition_3840945.csv'

load aortic_no_partition_384_final_data.mat

indices_global = params.layer_indices(1).indices_global; 

valve_with_reference.leaflets(1) = import_valve_from_csv(valve_with_reference.leaflets(1), file_name, indices_global, data_dir)


INSTALL(FILES
	DataAllozymes(diploids_and_tetraploids).txt
	DataSSR(with_spatial_distance_matrix).txt
	DataAFLP(estimation_of_marker_heritability).txt
	DataHaplo(with_distance_matrix_between_haplotypes).txt
	DataAFLP(kinship_estimation_corrected_by_heritability).txt 
	DataSSR(fine-scale_spatial_structure).txt
	
DESTINATION "${CMAKE_DATA_DIR}/example_data")


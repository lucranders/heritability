"""
This is a boilerplate pipeline 'calculateZZt'
generated using Kedro 0.17.6
"""

from kedro.pipeline import pipeline, node,Pipeline
from .nodes import calculate_K_matrix, make_K_matrix_non_negative

def create_pipeline(**kwargs):
	template_calculate_matrix = Pipeline([
	node(
			func=calculate_K_matrix,
			inputs=['params:params_run',"input_data"],
			outputs="specified_matrices",
			name="calculate_matrices",
			tags=['Genotype','Matrix','Calculation']
	),
	node(
			func=make_K_matrix_non_negative,
			inputs=['params:params_run',"specified_matrices"],
			outputs="non_negative_specified_matrices",
			name="correction_to_non_negative_specified_matrices",
			tags=['Genotype','Matrix','Correction']
	),
	])

	# template_GCTA = Pipeline([
	# node(
	# 		func=calculate_GCTA,
	# 		inputs=['params:overrideSetsChrs',"selected_Sample","params:gcta","bed_Files"],
	# 		outputs="gcta_matrix",
	# 		name="calculate_gcta_matrix",
	# 		tags=['Genotype','Matrix','Calculation']
	# ),
	# node(
	# 		func=correct_GCTA,
	# 		inputs=['params:overrideSetsChrs',"selected_Sample","gcta_matrix"],
	# 		outputs="corrected_gcta_matrix",
	# 		name="correction_gcta_matrix",
	# 		tags=['Genotype','Matrix','Correction']
	# ),
	# ])


	calculate_matrices_pipe = pipeline(
	pipe=template_calculate_matrix,
	inputs={'input_data':'_3:bed_file_creation.bed_files'},
	namespace="_4:matrices_calculation"
	)
	# fullSetChromosomesGCTA = pipeline(
	# pipe=template_GCTA,
	# parameters={"params:overrideSetsChrs": "params:listChrsFull", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
	# namespace="herit_" + setName_
	# )
	# fullSetComplete = fullSetChromosomes + fullSetChromosomesGCTA

	# setName_ = 'TwoSets'
	# twoSetsChromosomes = pipeline(
	# pipe=templateK_C_alphaCalculus,
	# parameters={"params:overrideSetsChrs": "params:listChrsPartitioned2", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
	# namespace="herit_" + setName_
	# )

	# twoSetsChromosomesGCTA = pipeline(
	# pipe=template_GCTA,
	# parameters={"params:overrideSetsChrs": "params:listChrsPartitioned2", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
	# namespace="herit_" + setName_
	# )

	# twoSetsComplete = twoSetsChromosomes + twoSetsChromosomesGCTA

#     setName_ = 'TwentyTwoSets'
#     twentyTwoSetsChromosomes = pipeline(
#         pipe=templateK_C_alphaCalculus,
#         parameters={"params:overrideSetsChrs": "params:listChrsPartitioned22", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
#         namespace="herit_" + setName_
#     )

#     entirepipeline = fullSetChromosomes + twoSetsChromosomes + twentyTwoSetsChromosomes
	# entirepipeline = fullSetComplete + twoSetsComplete
	return calculate_matrices_pipe
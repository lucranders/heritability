"""
This is a boilerplate pipeline 'calculateZZt'
generated using Kedro 0.17.6
"""

from kedro.pipeline import pipeline, node,Pipeline
from .nodes import calculate_K_C_alpha, correct_K_C_alpha, calculate_GCTA, correct_GCTA

def create_pipeline(**kwargs):
	templateK_C_alphaCalculus = Pipeline([
	node(
			func=calculate_K_C_alpha,
			inputs=['params:overrideSetsChrs',"selected_Sample","params:alphas","bed_Files"],
			outputs="k_c_alpha_matrix",
			name="calculate_k_c_alpha_matrix",
			tags=['Genotype','Matrix','Calculation']
	),
	node(
			func=correct_K_C_alpha,
			inputs=['params:overrideSetsChrs',"selected_Sample","params:alphas","k_c_alpha_matrix"],
			outputs="corrected_k_c_alpha_matrix",
			name="correction_k_c_alpha_matrix",
			tags=['Genotype','Matrix','Correction']
	),
	])

	template_GCTA = Pipeline([
	node(
			func=calculate_GCTA,
			inputs=['params:overrideSetsChrs',"selected_Sample","params:gcta","bed_Files"],
			outputs="gcta_matrix",
			name="calculate_gcta_matrix",
			tags=['Genotype','Matrix','Calculation']
	),
	node(
			func=correct_GCTA,
			inputs=['params:overrideSetsChrs',"selected_Sample","gcta_matrix"],
			outputs="corrected_gcta_matrix",
			name="correction_gcta_matrix",
			tags=['Genotype','Matrix','Correction']
	),
	])


	setName_ = 'FullSet'
	fullSetChromosomes = pipeline(
	pipe=templateK_C_alphaCalculus,
	parameters={"params:overrideSetsChrs": "params:listChrsFull", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
	namespace="herit_" + setName_
	)
	fullSetChromosomesGCTA = pipeline(
	pipe=template_GCTA,
	parameters={"params:overrideSetsChrs": "params:listChrsFull", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
	namespace="herit_" + setName_
	)
	fullSetComplete = fullSetChromosomes + fullSetChromosomesGCTA

	setName_ = 'TwoSets'
	twoSetsChromosomes = pipeline(
	pipe=templateK_C_alphaCalculus,
	parameters={"params:overrideSetsChrs": "params:listChrsPartitioned2", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
	namespace="herit_" + setName_
	)

	twoSetsChromosomesGCTA = pipeline(
	pipe=template_GCTA,
	parameters={"params:overrideSetsChrs": "params:listChrsPartitioned2", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
	namespace="herit_" + setName_
	)

	twoSetsComplete = twoSetsChromosomes + twoSetsChromosomesGCTA

#     setName_ = 'TwentyTwoSets'
#     twentyTwoSetsChromosomes = pipeline(
#         pipe=templateK_C_alphaCalculus,
#         parameters={"params:overrideSetsChrs": "params:listChrsPartitioned22", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
#         namespace="herit_" + setName_
#     )

#     entirepipeline = fullSetChromosomes + twoSetsChromosomes + twentyTwoSetsChromosomes
	entirepipeline = fullSetComplete + twoSetsComplete
	return entirepipeline
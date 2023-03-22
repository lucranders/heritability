from kedro.pipeline import Pipeline, node, pipeline
from .nodes import select_sample ,filter_snps, remove_bad_snps, create_bed_files

def create_pipeline(**kwargs):
	first_pipe_temp = Pipeline([
		node(
			select_sample,
			["dataset1" , "dataset2" , "params:params_run"],
			outputs="selected_sample",
			name="sample_selection",
		),
		])
	second_pipe_temp = Pipeline([
		node(
			filter_snps,
			["params:params_run","params:override"],
			outputs="initial_list_snps",
			name="build_first_list_snps",
		),
		# node(
		# 	remove_bad_snps,
		# 	["initial_list_snps"],
		# 	outputs="final_list_snps",
		# 	name="build_second_list_snps",
		# ),
	])
	third_pipe_temp = Pipeline([
		node(
			create_bed_files,
			["params:override"],
			outputs="bed_files",
			name="build_bed_files",
		),
		])
	first_pipe = pipeline(
	pipe=first_pipe_temp,
	inputs={"dataset1": "phenotypeData","dataset2": "genotypedData"},
	namespace="_1:select_sample"
	)
	second_pipe = pipeline(
	pipe=second_pipe_temp,
    parameters={"params:override": "_1:select_sample.selected_sample"},
	namespace="_2:filter_snps"
	)
	third_pipe = pipeline(
	pipe=third_pipe_temp,
    parameters={"params:override": "_2:filter_snps.initial_list_snps"},
	namespace="_3:bed_file_creation"
	)
	final_pipe = first_pipe + second_pipe + third_pipe
	return final_pipe

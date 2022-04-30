from kedro.pipeline import Pipeline, node
from .nodes import createTempFolder, selectSample ,createBedFiles

def create_pipeline(**kwargs):
    return Pipeline(
        [node(
            createTempFolder,
            ['params:snpsParams','params:sampParams'],
            outputs='params:pathTempFiles',
            name="create_temporary_folder_for_further_analysis",
        ),
        node(
            selectSample,
            ["phenotypeData" , "genotypedData" , "params:sampParams","params:pathTempFiles"],
            outputs="selected_Sample",
            name="sample_selection",
        ),
        node(
            createBedFiles,
            ["params:pathTempFiles","params:snpsParams","selected_Sample"],
            outputs="bed_Files",
            name="filter_desired_snps_based_on_parameters",
        ),
        ]
    )

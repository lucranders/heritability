from kedro.pipeline import Pipeline, node
from .nodes import createTempFolder, selectSample ,createBedFiles

def create_pipeline(**kwargs):
    return Pipeline(
        [node(
            selectSample,
            ["phenotypeData" , "genotypedData" , "params:sampParams","params:snpsParams"],
            outputs="selected_Sample",
            name="sample_selection",
        ),
        node(
            createBedFiles,
            ["params:snpsParams","selected_Sample"],
            outputs="bed_Files",
            name="filter_desired_snps_based_on_parameters",
        ),
        ]
    )

from kedro.pipeline import Pipeline, node
from .nodes import selectSample,createTmpFolder,createBedFiles

def create_pipeline(**kwargs):
    return Pipeline(
        [node(
            createTmpFolder,
            ["params:pathTempFiles" , "params:snpsParams"],
            outputs="params:tempFolderPath",
            name="temporaryFolderCreation",
        ),
        node(
            selectSample,
            ["phenotypeData" , "params:sampParams","params:tempFolderPath"],
            outputs="selectedSample",
            name="sampleSelection",
        ),
        node(
            createBedFiles,
            ["params:plink","params:tempFolderPath","params:pathVcfFiles","params:snpsParams"],
            outputs=None,
            name="filterSnps",
        ),
        ]
    )

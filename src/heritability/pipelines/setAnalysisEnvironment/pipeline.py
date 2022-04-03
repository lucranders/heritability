from kedro.pipeline import Pipeline, node
from .nodes import createTempFolder, selectSample ,createBedFiles

def create_pipeline(**kwargs):
    return Pipeline(
        [node(
            createTempFolder,
            ['params:snpsParams', 'params:pathTemp'],
            outputs='params:pathTempFiles',
            name="createTemporaryFolder",
        ),
        node(
            selectSample,
            ["phenotypeData" , "params:sampParams","params:pathTempFiles"],
            outputs="selectedSample",
            name="sampleSelection",
        ),
        node(
            createBedFiles,
            ["params:plink","params:tempFolderPath","params:pathVcfFiles","params:snpsParams","selectedSample"],
            outputs="bedFiles",
            name="filterSnps",
        ),
        ]
    )

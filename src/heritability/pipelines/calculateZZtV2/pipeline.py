"""
This is a boilerplate pipeline 'calculateZZtV2'
generated using Kedro 0.17.6
"""

from kedro.pipeline import Pipeline, node
from heritability.pipelines.calculateZZt.nodes import calculateGCTA,correctGRM

def create_pipeline(**kwargs):
    return Pipeline([
        node(
            calculateGCTA,
            ["params:nameFilePt1",'params:listChrsPt1','params:nameMatrixPt1',"params:tempFolderPath","params:gcta","params:numThreads","bed_Files"],
            outputs='calculated_ZZt_without_chromosome_6',
            name="calculate_ZZt_without_chromosome_6",
        ),
        node(
            calculateGCTA,
            ["params:nameFilePt2",'params:listChrsPt2','params:nameMatrixPt2',"params:tempFolderPath","params:gcta","params:numThreads","bed_Files"],
            outputs='calculated_ZZt_only_chromosome_6',
            name="calculate_ZZt_only_chromosome_6",
        ),
        node(
            correctGRM,
            ["params:tempFolderPath",'params:nameMatrixPt1',"calculated_ZZt_without_chromosome_6"],
            outputs='corrected_ZZt_without_chromosome_6',
            name="ZZt_without_chromosome_6_corrected",
        ),
        node(
            correctGRM,
            ["params:tempFolderPath",'params:nameMatrixPt2',"calculated_ZZt_only_chromosome_6"],
            outputs='corrected_ZZt_only_chromosome_6',
            name="ZZt_only_chromosome_6_corrected",
        ),
    ])

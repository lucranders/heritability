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
            ['params:listChrsPartitioned2',"selected_Sample","params:numThreads","bed_Files"],
            outputs='calculated_ZZt_partitioned_2_components',
            name="calculate_ZZt_partitioned_2_components",
        ),
        node(
            correctGRM,
            ["selected_Sample",'params:listChrsPartitioned2',"calculated_ZZt_partitioned_2_components"],
            outputs='corrected_ZZt_partitioned_2_components',
            name="ZZt_partitioned_2_components",
        ),
        node(
            calculateGCTA,
            ['params:listChrsPartitioned22',"selected_Sample","params:numThreads","bed_Files"],
            outputs='calculated_ZZt_partitioned_22_components',
            name="calculate_ZZt_partitioned_22_components",
        ),
        node(
            correctGRM,
            ["selected_Sample",'params:listChrsPartitioned22',"calculated_ZZt_partitioned_22_components"],
            outputs='corrected_ZZt_partitioned_22_components',
            name="ZZt_partitioned_22_components",
        ),
    ])

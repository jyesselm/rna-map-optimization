"""RNA-MAP Optimization Toolkit.

Parameter optimization tools for RNA-MAP workflows.
"""

__version__ = "0.1.0"
__author__ = "Joe Yesselman"
__email__ = "jyesselm@unl.edu"

from rna_map_optimization.bowtie2 import (
    build_bowtie2_index,
    generate_bit_vectors_and_analyze,
    params_to_bowtie2_args,
    parse_bowtie2_stats,
    parse_sam_quality,
    run_bowtie2_alignment,
)
from rna_map_optimization.optimizer import run_optimization

__all__ = [
    "build_bowtie2_index",
    "generate_bit_vectors_and_analyze",
    "params_to_bowtie2_args",
    "parse_bowtie2_stats",
    "parse_sam_quality",
    "run_bowtie2_alignment",
    "run_optimization",
    "cli",
    "optimizer",
    "bowtie2",
    "utils",
]


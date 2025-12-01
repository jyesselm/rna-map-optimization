"""Bowtie2 alignment utilities.

This module provides common functions for running Bowtie2 alignments,
parsing results, and converting parameters to command-line arguments.
"""

import subprocess
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from rna_map_mini.logger import get_logger
except ImportError:
    try:
        from rna_map_mini.core.logger import get_logger
    except ImportError:
        import logging

        logging.basicConfig(level=logging.INFO)

        def get_logger(name: str):
            """Create a logger instance."""
            logger = logging.getLogger(name)
            if not logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter(
                    "%(name)s - %(levelname)s - %(message)s"
                )
                handler.setFormatter(formatter)
                logger.addHandler(handler)
                logger.setLevel(logging.INFO)
            return logger

try:
    from rna_map_mini.analysis.bit_vector_iterator import BitVectorIterator
    from rna_map_mini.analysis.mutation_histogram import MutationHistogram
except ImportError:
    try:
        from rna_map_mini.core.bit_vector_iterator import BitVectorIterator
        from rna_map_mini.core.mutation_histogram import MutationHistogram
    except ImportError:
        try:
            from rna_map_mini.bit_vector_iterator import BitVectorIterator
            from rna_map_mini.mutation_histogram import MutationHistogram
        except ImportError:
            raise ImportError(
                "Could not import BitVectorIterator or MutationHistogram from rna_map_mini. "
                "Please ensure rna-map-mini is installed: "
                "pip install git+https://github.com/jyesselm/rna-map-mini.git"
            )

log = get_logger("BOWTIE2_UTILS")


def parse_bowtie2_stats(stderr_file: Path) -> Dict[str, int | float | bool]:
    """Parse Bowtie2 alignment statistics from stderr output.

    Args:
        stderr_file: Path to Bowtie2 stderr output

    Returns:
        Dictionary with alignment statistics
    """
    stats: Dict[str, int | float | bool] = {
        "total_reads": 0,
        "aligned_0_times": 0,
        "aligned_exactly_1_time": 0,
        "aligned_more_than_1_time": 0,
        "overall_alignment_rate": 0.0,
        "concordant_0_times": 0,
        "concordant_exactly_1_time": 0,
        "concordant_more_than_1_time": 0,
        "is_paired_end": False,
    }

    if not stderr_file.exists():
        return stats

    with open(stderr_file) as f:
        content = f.read()
        lines = content.split("\n")

        for line in lines:
            line = line.strip()
            if "reads; of these:" in line:
                stats["total_reads"] = int(line.split()[0])
            elif "were paired" in line:
                stats["is_paired_end"] = True
            elif "aligned 0 times" in line and "concordantly" not in line:
                stats["aligned_0_times"] = int(line.split()[0])
            elif "aligned exactly 1 time" in line and "concordantly" not in line:
                stats["aligned_exactly_1_time"] = int(line.split()[0])
            elif "aligned >1 times" in line and "concordantly" not in line:
                stats["aligned_more_than_1_time"] = int(line.split()[0])
            elif "aligned concordantly 0 times" in line:
                stats["concordant_0_times"] = int(line.split()[0])
            elif "aligned concordantly exactly 1 time" in line:
                stats["concordant_exactly_1_time"] = int(line.split()[0])
            elif "aligned concordantly >1 times" in line:
                stats["concordant_more_than_1_time"] = int(line.split()[0])
            elif "overall alignment rate" in line:
                rate_str = line.split("%")[0].split()[-1]
                stats["overall_alignment_rate"] = float(rate_str)

    return stats


def parse_sam_quality(
    sam_file: Path, mapq_cutoff: int = 20
) -> Dict[str, int | float | Dict[int, int]]:
    """Parse SAM file to calculate quality metrics.

    Args:
        sam_file: Path to SAM file
        mapq_cutoff: Minimum MAPQ score for high-quality alignments

    Returns:
        Dictionary with quality metrics
    """
    metrics: Dict[str, int | float | Dict[int, int]] = {
        "total_alignments": 0,
        "high_quality_alignments": 0,
        "avg_mapq": 0.0,
        "median_mapq": 0.0,
        "mapq_distribution": {},
    }

    if not sam_file.exists():
        return metrics

    mapq_scores: List[int] = []
    mapq_counts: Dict[int, int] = {}

    with open(sam_file) as f:
        for line in f:
            if line.startswith("@"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue

            metrics["total_alignments"] = int(metrics["total_alignments"]) + 1
            mapq = int(fields[4])
            mapq_scores.append(mapq)
            mapq_counts[mapq] = mapq_counts.get(mapq, 0) + 1

            if mapq >= mapq_cutoff:
                metrics["high_quality_alignments"] = (
                    int(metrics["high_quality_alignments"]) + 1
                )

    if mapq_scores:
        metrics["avg_mapq"] = sum(mapq_scores) / len(mapq_scores)
        sorted_scores = sorted(mapq_scores)
        mid = len(sorted_scores) // 2
        metrics["median_mapq"] = (
            sorted_scores[mid]
            if len(sorted_scores) % 2 == 1
            else (sorted_scores[mid - 1] + sorted_scores[mid]) / 2
        )
        metrics["mapq_distribution"] = mapq_counts

    return metrics


def build_bowtie2_index(fasta: Path, index_dir: Path) -> Path:
    """Build Bowtie2 index from FASTA file.

    Args:
        fasta: Path to FASTA file
        index_dir: Directory to store index

    Returns:
        Path to index basename
    """
    index_dir.mkdir(parents=True, exist_ok=True)
    index_name = index_dir / fasta.stem

    if (index_name.with_suffix(".1.bt2")).exists():
        log.info(f"Index already exists: {index_name}")
        return index_name

    log.info(f"Building Bowtie2 index: {index_name}")
    cmd = ["bowtie2-build", str(fasta), str(index_name)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Failed to build index: {result.stderr}")

    return index_name


def run_bowtie2_alignment(
    index: Path,
    fastq1: Path,
    fastq2: Optional[Path],
    output_sam: Path,
    bt2_args: List[str],
    threads: int = 4,
) -> Tuple[Path, float]:
    """Run Bowtie2 alignment with given parameters.

    Args:
        index: Path to Bowtie2 index
        fastq1: Path to first FASTQ file
        fastq2: Path to second FASTQ file (optional)
        output_sam: Path to output SAM file
        bt2_args: List of Bowtie2 arguments
        threads: Number of threads

    Returns:
        Tuple of (stderr_file, elapsed_time)
    """
    stderr_file = output_sam.parent / "bowtie2_stderr.txt"

    cmd = ["bowtie2"] + bt2_args + ["-x", str(index), "-S", str(output_sam)]

    if fastq2 and fastq2.exists():
        cmd.extend(["-1", str(fastq1), "-2", str(fastq2)])
    else:
        cmd.extend(["-U", str(fastq1)])

    if not any(arg.startswith("-p") for arg in bt2_args):
        cmd.extend(["-p", str(threads)])

    log.info(f"Running: {' '.join(cmd)}")

    start_time = time.time()
    with open(stderr_file, "w") as stderr:
        result = subprocess.run(
            cmd, stderr=stderr, stdout=subprocess.PIPE, text=True
        )

    elapsed_time = time.time() - start_time

    if result.returncode != 0:
        log.warning(f"Bowtie2 returned non-zero exit code: {result.returncode}")

    return stderr_file, elapsed_time


def generate_bit_vectors_and_analyze(
    sam_file: Path,
    ref_seqs: Dict[str, str],
    paired: bool,
    qscore_cutoff: int = 25,
    mapq_cutoff: int = 20,
) -> Dict:
    """Generate bit vectors and extract quality metrics including signal-to-noise.

    Args:
        sam_file: Path to SAM file
        ref_seqs: Dictionary of reference sequences
        paired: Whether reads are paired-end
        qscore_cutoff: Quality score cutoff for bit vector generation
        mapq_cutoff: MAPQ cutoff for filtering

    Returns:
        Dictionary with bit vector metrics including signal-to-noise ratio
        and per-construct statistics
    """
    metrics: Dict = {
        "total_bit_vectors": 0,
        "accepted_bit_vectors": 0,
        "rejected_low_mapq": 0,
        "mutation_distribution": defaultdict(int),
        "avg_mutations_per_read": 0.0,
        "reads_with_mutations": 0,
        "reads_without_mutations": 0,
        "total_mutations": 0,
        "coverage_positions": 0,
        "avg_coverage": 0.0,
        "signal_to_noise": 0.0,
        "constructs": {},
    }

    if not sam_file.exists():
        return metrics

    try:
        mut_histos: Dict[str, MutationHistogram] = {}
        for ref_name, ref_seq in ref_seqs.items():
            mut_histos[ref_name] = MutationHistogram(
                ref_name, ref_seq, "DMS", 1, len(ref_seq)
            )

        iterator = BitVectorIterator(
            sam_file, ref_seqs, paired=paired, use_pysam=False, qscore_cutoff=qscore_cutoff
        )

        for bit_vector in iterator:
            metrics["total_bit_vectors"] += 1

            if not bit_vector.reads:
                continue

            if any(read.mapq < mapq_cutoff for read in bit_vector.reads):
                metrics["rejected_low_mapq"] += 1
                continue

            if not bit_vector.data:
                continue

            metrics["accepted_bit_vectors"] += 1

            ref_name = bit_vector.reads[0].rname
            if ref_name in mut_histos:
                mh = mut_histos[ref_name]
                mh.num_reads += 1
                mh.num_aligned += 1

                total_muts = 0
                for pos in mh.get_nuc_coords():
                    if pos not in bit_vector.data:
                        continue
                    read_bit = bit_vector.data[pos]
                    mh.cov_bases[pos] += 1

                    if read_bit in ["A", "C", "G", "T"]:
                        mh.mut_bases[pos] += 1
                        mh.mod_bases[read_bit][pos] += 1
                        total_muts += 1
                    elif read_bit == "1":
                        mh.del_bases[pos] += 1
                    elif read_bit == "?":
                        mh.info_bases[pos] += 1

                mh.num_of_mutations[total_muts] += 1

            mut_count = sum(
                1 for v in bit_vector.data.values() if v in ["A", "C", "G", "T"]
            )
            metrics["mutation_distribution"][mut_count] += 1
            metrics["total_mutations"] += mut_count

            if mut_count > 0:
                metrics["reads_with_mutations"] += 1
            else:
                metrics["reads_without_mutations"] += 1

            metrics["coverage_positions"] += len(bit_vector.data)

        total_snr = 0.0
        total_weight = 0
        for ref_name, mh in mut_histos.items():
            if mh.num_aligned > 0:
                snr = mh.get_signal_to_noise()
                total_snr += snr * mh.num_aligned
                total_weight += mh.num_aligned

                metrics["constructs"][ref_name] = {
                    "aligned_reads": mh.num_aligned,
                    "signal_to_noise": snr,
                    "sequence_length": len(mh.ref_seq),
                    "total_reads": mh.num_reads,
                }

        if total_weight > 0:
            metrics["signal_to_noise"] = total_snr / total_weight

        if metrics["accepted_bit_vectors"] > 0:
            metrics["avg_mutations_per_read"] = (
                metrics["total_mutations"] / metrics["accepted_bit_vectors"]
            )
            metrics["avg_coverage"] = (
                metrics["coverage_positions"] / metrics["accepted_bit_vectors"]
            )

        metrics["mutation_distribution"] = dict(metrics["mutation_distribution"])

    except Exception as e:
        log.warning(f"Error generating bit vectors: {e}")

    return metrics


def params_to_bowtie2_args(params: Dict) -> List[str]:
    """Convert parameter dictionary to Bowtie2 argument list.

    Args:
        params: Parameter dictionary

    Returns:
        List of Bowtie2 arguments
    """
    args: List[str] = []

    flag_map: Dict[str, str] = {
        "local": "--local",
        "end_to_end": "--end-to-end",
        "no_unal": "--no-unal",
        "no_discordant": "--no-discordant",
        "no_mixed": "--no-mixed",
        "seed_length": "-L",
        "seed_mismatches": "-N",
        "maxins": "-X",
        "minins": "-I",
        "score_min": "--score-min",
        "mismatch_penalty": "--mp",
        "gap_penalty_read": "--rdg",
        "gap_penalty_ref": "--rfg",
        "seed_interval": "-i",
        "np_penalty": "--np",
        "n_ceil": "--n-ceil",
        "gbar": "--gbar",
        "match_bonus": "--ma",
        "extension_effort": "-D",
        "repetitive_effort": "-R",
        "very_fast_local": "--very-fast-local",
        "fast_local": "--fast-local",
        "sensitive_local": "--sensitive-local",
        "very_sensitive_local": "--very-sensitive-local",
    }

    for key, value in params.items():
        if key in flag_map and value:
            if isinstance(value, bool) and value:
                args.append(flag_map[key])
            elif isinstance(value, (int, str)):
                args.extend([flag_map[key], str(value)])

    return args


#!/usr/bin/env python3
"""Verify alignments using secondary alignment tools.

This script uses alternative alignment methods to verify Bowtie2 alignments
for small subsamples. Useful for:
- Validating that aligned reads should actually align
- Checking if unaligned reads should have aligned
- Comparing alignment quality between methods
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from Bio import pairwise2
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("Warning: Biopython not available. Install with: pip install biopython")

try:
    import parasail
    PARASAIL_AVAILABLE = True
except ImportError:
    PARASAIL_AVAILABLE = False
    print("Warning: parasail not available. Install with: pip install parasail")


def check_minimap2() -> bool:
    """Check if minimap2 is available."""
    try:
        result = subprocess.run(
            ["minimap2", "--version"],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_blast() -> bool:
    """Check if BLAST is available."""
    try:
        result = subprocess.run(
            ["blastn", "-version"],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def read_fasta(fasta_path: Path) -> Dict[str, str]:
    """Read FASTA file."""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            sequences[current_name] = "".join(current_seq)
    
    return sequences


def read_fastq(fastq_path: Path, max_reads: Optional[int] = None) -> List[Tuple[str, str, str]]:
    """Read FASTQ file and return list of (name, sequence, quality)."""
    reads = []
    with open(fastq_path) as f:
        while True:
            name_line = f.readline().strip()
            if not name_line:
                break
            if name_line.startswith("@"):
                name = name_line[1:].split()[0]
                seq = f.readline().strip()
                f.readline()  # Skip +
                qual = f.readline().strip()
                reads.append((name, seq, qual))
                if max_reads and len(reads) >= max_reads:
                    break
    return reads


def parse_sam_alignment(sam_line: str) -> Optional[Dict]:
    """Parse a single SAM line."""
    fields = sam_line.strip().split("\t")
    if len(fields) < 11 or fields[0].startswith("@"):
        return None
    
    return {
        "qname": fields[0],
        "flag": int(fields[1]),
        "rname": fields[2],
        "pos": int(fields[3]) if fields[3] != "0" else None,
        "mapq": int(fields[4]),
        "cigar": fields[5],
        "seq": fields[9],
        "qual": fields[10],
    }


def biopython_align(seq1: str, seq2: str, method: str = "global") -> Dict:
    """Perform pairwise alignment using Biopython.
    
    Args:
        seq1: First sequence (read)
        seq2: Second sequence (reference)
        method: 'global' or 'local'
        
    Returns:
        Dictionary with alignment metrics
    """
    if not BIOPYTHON_AVAILABLE:
        return {"error": "Biopython not available"}
    
    # Scoring: match=2, mismatch=-1, gap_open=-2, gap_extend=-1
    if method == "global":
        alignments = pairwise2.align.globalms(
            seq1, seq2, 2, -1, -2, -1
        )
    else:
        alignments = pairwise2.align.localms(
            seq1, seq2, 2, -1, -2, -1
        )
    
    if not alignments:
        return {"aligned": False, "score": 0, "identity": 0.0}
    
    best = alignments[0]
    aligned_seq1, aligned_seq2, score, start, end = best
    
    # Calculate identity
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != "-")
    total = len([x for x in aligned_seq1 if x != "-"])
    identity = matches / total if total > 0 else 0.0
    
    # Count mismatches and gaps
    mismatches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != "-" and b != "-")
    gaps = aligned_seq1.count("-") + aligned_seq2.count("-")
    
    return {
        "aligned": True,
        "score": score,
        "identity": identity,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "length": len(aligned_seq1),
        "alignment": (aligned_seq1, aligned_seq2),
    }


def parasail_align(seq1: str, seq2: str, method: str = "local") -> Dict:
    """Perform alignment using parasail (SIMD-accelerated).
    
    Args:
        seq1: First sequence (read)
        seq2: Second sequence (reference)
        method: 'local' or 'global'
        
    Returns:
        Dictionary with alignment metrics
    """
    if not PARASAIL_AVAILABLE:
        return {"error": "parasail not available"}
    
    # Scoring matrix: match=2, mismatch=-1
    matrix = parasail.matrix_create("ACGT", 2, -1)
    
    if method == "local":
        result = parasail.sw(seq1, seq2, 2, -1, matrix)
    else:
        result = parasail.nw(seq1, seq2, 2, -1, matrix)
    
    if result.score < 0:
        return {"aligned": False, "score": result.score}
    
    # Calculate identity from traceback
    traceback = parasail.trace(result, seq1, seq2, matrix)
    aligned_seq1 = traceback.traceback.query
    aligned_seq2 = traceback.traceback.ref
    
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != "-")
    total = len([x for x in aligned_seq1 if x != "-"])
    identity = matches / total if total > 0 else 0.0
    
    mismatches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != "-" and b != "-")
    gaps = aligned_seq1.count("-") + aligned_seq2.count("-")
    
    return {
        "aligned": True,
        "score": result.score,
        "identity": identity,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "length": len(aligned_seq1),
    }


def minimap2_align(
    fastq_path: Path,
    fasta_path: Path,
    output_sam: Path,
    preset: str = "sr",
) -> Dict:
    """Run minimap2 alignment.
    
    Args:
        fastq_path: Path to FASTQ file
        fasta_path: Path to reference FASTA
        output_sam: Path to output SAM file
        preset: Minimap2 preset (sr=short reads, map-ont, etc.)
        
    Returns:
        Dictionary with alignment statistics
    """
    if not check_minimap2():
        return {"error": "minimap2 not available"}
    
    cmd = [
        "minimap2",
        f"-x{preset}",
        "-a",  # SAM output
        str(fasta_path),
        str(fastq_path),
    ]
    
    with open(output_sam, "w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        return {"error": f"minimap2 failed: {result.stderr}"}
    
    # Parse SAM to get statistics
    stats = {
        "total_reads": 0,
        "aligned": 0,
        "unaligned": 0,
        "avg_identity": 0.0,
    }
    
    identities = []
    with open(output_sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            stats["total_reads"] += 1
            fields = line.strip().split("\t")
            if len(fields) >= 11 and fields[2] != "*":
                stats["aligned"] += 1
                # Try to extract identity from NM tag
                for field in fields[11:]:
                    if field.startswith("NM:i:"):
                        nm = int(field.split(":")[2])
                        read_len = len(fields[9])
                        identity = 1.0 - (nm / read_len) if read_len > 0 else 0.0
                        identities.append(identity)
                        break
            else:
                stats["unaligned"] += 1
    
    if identities:
        stats["avg_identity"] = sum(identities) / len(identities)
    
    return stats


def verify_bowtie2_alignment(
    sam_file: Path,
    fasta: Path,
    fastq: Path,
    max_reads: int = 100,
    method: str = "biopython",
) -> Dict:
    """Verify Bowtie2 alignments using secondary method.
    
    Args:
        sam_file: Path to Bowtie2 SAM output
        fasta: Path to reference FASTA
        fastq: Path to original FASTQ file
        max_reads: Maximum number of reads to verify
        method: Verification method ('biopython', 'parasail', 'minimap2')
        
    Returns:
        Dictionary with verification results
    """
    ref_seqs = read_fasta(fasta)
    reads_dict = {name: (seq, qual) for name, seq, qual in read_fastq(fastq, max_reads)}
    
    results = {
        "total_checked": 0,
        "aligned_in_bowtie2": 0,
        "aligned_in_verify": 0,
        "agreement": 0,
        "disagreement": 0,
        "bowtie2_only": 0,
        "verify_only": 0,
        "identity_stats": [],
        "disagreements": [],
    }
    
    # Read SAM file
    bowtie2_alignments = {}
    with open(sam_file) as f:
        for line in f:
            if line.startswith("@"):
                continue
            aln = parse_sam_alignment(line)
            if aln and aln["pos"] is not None:
                bowtie2_alignments[aln["qname"]] = aln
    
    # Verify each read
    for read_name, (read_seq, read_qual) in reads_dict.items():
        if read_name not in bowtie2_alignments:
            # Read didn't align in Bowtie2 - check if it should
            results["total_checked"] += 1
            continue
        
        results["aligned_in_bowtie2"] += 1
        results["total_checked"] += 1
        
        bowtie2_aln = bowtie2_alignments[read_name]
        ref_name = bowtie2_aln["rname"]
        
        if ref_name not in ref_seqs:
            continue
        
        ref_seq = ref_seqs[ref_name]
        
        # Perform verification alignment
        if method == "biopython":
            verify_result = biopython_align(read_seq, ref_seq, method="local")
        elif method == "parasail":
            verify_result = parasail_align(read_seq, ref_seq, method="local")
        else:
            continue
        
        if verify_result.get("aligned", False):
            results["aligned_in_verify"] += 1
            identity = verify_result.get("identity", 0.0)
            results["identity_stats"].append(identity)
            
            # Check agreement
            if identity > 0.7:  # Reasonable threshold
                results["agreement"] += 1
            else:
                results["disagreement"] += 1
                results["disagreements"].append({
                    "read": read_name,
                    "bowtie2_mapq": bowtie2_aln["mapq"],
                    "verify_identity": identity,
                    "verify_score": verify_result.get("score", 0),
                })
        else:
            results["bowtie2_only"] += 1
            results["disagreements"].append({
                "read": read_name,
                "bowtie2_mapq": bowtie2_aln["mapq"],
                "issue": "Aligned in Bowtie2 but not in verification",
            })
    
    # Calculate statistics
    if results["identity_stats"]:
        results["avg_identity"] = sum(results["identity_stats"]) / len(results["identity_stats"])
        results["min_identity"] = min(results["identity_stats"])
        results["max_identity"] = max(results["identity_stats"])
    
    return results


def main():
    """Main verification function."""
    parser = argparse.ArgumentParser(
        description="Verify Bowtie2 alignments using secondary methods"
    )
    
    parser.add_argument(
        "--sam",
        type=Path,
        required=True,
        help="Path to Bowtie2 SAM output file",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Path to reference FASTA file",
    )
    parser.add_argument(
        "--fastq",
        type=Path,
        required=True,
        help="Path to original FASTQ file",
    )
    parser.add_argument(
        "--method",
        choices=["biopython", "parasail", "minimap2"],
        default="biopython",
        help="Verification method (default: biopython)",
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=100,
        help="Maximum number of reads to verify (default: 100)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("verification_results.json"),
        help="Output JSON file for results",
    )
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("Alignment Verification")
    print("=" * 80)
    print(f"SAM file: {args.sam}")
    print(f"Reference: {args.fasta}")
    print(f"Reads: {args.fastq}")
    print(f"Method: {args.method}")
    print(f"Max reads: {args.max_reads}")
    print()
    
    # Check method availability
    if args.method == "biopython" and not BIOPYTHON_AVAILABLE:
        print("ERROR: Biopython not available. Install with: pip install biopython")
        return 1
    
    if args.method == "parasail" and not PARASAIL_AVAILABLE:
        print("ERROR: parasail not available. Install with: pip install parasail")
        return 1
    
    if args.method == "minimap2" and not check_minimap2():
        print("ERROR: minimap2 not available. Install with: conda install -c bioconda minimap2")
        return 1
    
    # Run verification
    if args.method == "minimap2":
        # For minimap2, run full alignment
        output_sam = args.output.parent / "minimap2_verification.sam"
        results = minimap2_align(args.fastq, args.fasta, output_sam)
        print("Minimap2 alignment results:")
        print(f"  Total reads: {results.get('total_reads', 0)}")
        print(f"  Aligned: {results.get('aligned', 0)}")
        print(f"  Unaligned: {results.get('unaligned', 0)}")
        print(f"  Avg identity: {results.get('avg_identity', 0.0):.2%}")
    else:
        # For biopython/parasail, verify individual alignments
        results = verify_bowtie2_alignment(
            args.sam, args.fasta, args.fastq, args.max_reads, args.method
        )
        
        print("Verification Results:")
        print(f"  Total reads checked: {results['total_checked']}")
        print(f"  Aligned in Bowtie2: {results['aligned_in_bowtie2']}")
        print(f"  Aligned in {args.method}: {results['aligned_in_verify']}")
        print(f"  Agreement: {results['agreement']}")
        print(f"  Disagreement: {results['disagreement']}")
        print(f"  Bowtie2 only: {results['bowtie2_only']}")
        
        if results.get("avg_identity"):
            print(f"  Avg identity: {results['avg_identity']:.2%}")
            print(f"  Min identity: {results['min_identity']:.2%}")
            print(f"  Max identity: {results['max_identity']:.2%}")
        
        if results["disagreements"]:
            print(f"\n  Disagreements found: {len(results['disagreements'])}")
            print("  First 5 disagreements:")
            for d in results["disagreements"][:5]:
                print(f"    - {d['read']}: MAPQ={d.get('bowtie2_mapq', 'N/A')}, "
                      f"Identity={d.get('verify_identity', 0.0):.2%}")
    
    # Save results
    import json
    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to: {args.output}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


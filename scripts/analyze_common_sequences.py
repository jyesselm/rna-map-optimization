#!/usr/bin/env python3
"""Analyze and optionally remove common sequences from multiple reference sequences.

This script helps identify common regions between reference sequences and can
create a new FASTA file with common sequences removed, keeping only unique regions.
"""

import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict


def fasta_to_dict(fasta_path: Path) -> Dict[str, str]:
    """Read FASTA file and return dictionary of sequence names to sequences."""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences


def find_common_substrings(
    sequences: Dict[str, str],
    min_length: int = 20
) -> List[Tuple[str, int, int, Set[str]]]:
    """Find common substrings across sequences.
    
    Args:
        sequences: Dictionary mapping sequence names to sequences
        min_length: Minimum length of common substring to report
        
    Returns:
        List of tuples: (substring, start_pos, end_pos, set_of_sequence_names)
    """
    seq_list = list(sequences.items())
    common_regions = []
    
    if len(seq_list) < 2:
        return common_regions
    
    # Compare all pairs
    for i, (name1, seq1) in enumerate(seq_list):
        for j, (name2, seq2) in enumerate(seq_list[i+1:], start=i+1):
            # Find longest common substring
            lcs = longest_common_substring(seq1, seq2, min_length)
            if lcs:
                substring, pos1, pos2 = lcs
                common_regions.append((substring, pos1, pos2, {name1, name2}))
    
    # Merge overlapping regions
    merged = merge_overlapping_regions(common_regions, sequences)
    return merged


def longest_common_substring(
    seq1: str,
    seq2: str,
    min_length: int = 20
) -> Tuple[str, int, int] | None:
    """Find longest common substring between two sequences.
    
    Returns:
        Tuple of (substring, position_in_seq1, position_in_seq2) or None
    """
    m, n = len(seq1), len(seq2)
    max_len = 0
    end_pos1 = 0
    end_pos2 = 0
    
    # Dynamic programming table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
                if dp[i][j] > max_len:
                    max_len = dp[i][j]
                    end_pos1 = i
                    end_pos2 = j
            else:
                dp[i][j] = 0
    
    if max_len >= min_length:
        start_pos1 = end_pos1 - max_len
        start_pos2 = end_pos2 - max_len
        substring = seq1[start_pos1:end_pos1]
        return (substring, start_pos1, start_pos2)
    
    return None


def merge_overlapping_regions(
    regions: List[Tuple[str, int, int, Set[str]]],
    sequences: Dict[str, str]
) -> List[Tuple[str, int, int, Set[str]]]:
    """Merge overlapping common regions."""
    if not regions:
        return []
    
    # Group by sequence pairs
    merged = []
    seen = set()
    
    for substring, pos1, pos2, seq_names in regions:
        key = (pos1, pos2, tuple(sorted(seq_names)))
        if key in seen:
            continue
        seen.add(key)
        merged.append((substring, pos1, pos2, seq_names))
    
    return merged


def find_unique_regions(
    sequences: Dict[str, str],
    common_regions: List[Tuple[str, int, int, Set[str]]],
    min_unique_length: int = 50
) -> Dict[str, List[Tuple[int, int]]]:
    """Find unique regions for each sequence (regions not in common).
    
    Args:
        sequences: Dictionary of sequences
        common_regions: List of common regions
        min_unique_length: Minimum length of unique region to keep
        
    Returns:
        Dictionary mapping sequence names to lists of (start, end) tuples
    """
    unique_regions = {}
    
    for seq_name, seq in sequences.items():
        # Find which common regions affect this sequence
        affected_regions = []
        for substring, pos1, pos2, seq_names in common_regions:
            if seq_name in seq_names:
                # Find position in this sequence
                if seq_name == list(seq_names)[0]:
                    affected_regions.append((pos1, pos1 + len(substring)))
                else:
                    affected_regions.append((pos2, pos2 + len(substring)))
        
        # Sort by start position
        affected_regions.sort()
        
        # Find unique regions (gaps between common regions)
        unique = []
        start = 0
        
        for common_start, common_end in affected_regions:
            if common_start > start:
                unique_len = common_start - start
                if unique_len >= min_unique_length:
                    unique.append((start, common_start))
            start = max(start, common_end)
        
        # Add final region if exists
        if start < len(seq):
            unique_len = len(seq) - start
            if unique_len >= min_unique_length:
                unique.append((start, len(seq)))
        
        unique_regions[seq_name] = unique
    
    return unique_regions


def remove_common_sequences(
    sequences: Dict[str, str],
    min_unique_length: int = 50
) -> Dict[str, str]:
    """Remove common sequences, keeping only unique regions.
    
    Args:
        sequences: Dictionary of sequences
        min_unique_length: Minimum length of unique region to keep
        
    Returns:
        Dictionary of sequences with common regions removed
    """
    common_regions = find_common_substrings(sequences, min_length=20)
    unique_regions = find_unique_regions(sequences, common_regions, min_unique_length)
    
    result = {}
    for seq_name, seq in sequences.items():
        if seq_name in unique_regions and unique_regions[seq_name]:
            # Concatenate unique regions with N's as separators
            unique_parts = [seq[start:end] for start, end in unique_regions[seq_name]]
            result[seq_name] = 'N' * 10 + 'N'.join(unique_parts) + 'N' * 10
        else:
            # No unique regions found, keep original
            result[seq_name] = seq
    
    return result


def write_fasta(sequences: Dict[str, str], output_path: Path):
    """Write sequences to FASTA file."""
    with open(output_path, 'w') as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze common sequences in multiple reference FASTA files"
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Input FASTA file with multiple sequences"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output FASTA file with common sequences removed (optional)"
    )
    parser.add_argument(
        "--min-common-length",
        type=int,
        default=20,
        help="Minimum length of common sequence to detect (default: 20)"
    )
    parser.add_argument(
        "--min-unique-length",
        type=int,
        default=50,
        help="Minimum length of unique region to keep (default: 50)"
    )
    parser.add_argument(
        "--analyze-only",
        action="store_true",
        help="Only analyze common sequences, don't create output file"
    )
    
    args = parser.parse_args()
    
    # Read sequences
    print(f"Reading sequences from {args.fasta}")
    sequences = fasta_to_dict(args.fasta)
    print(f"Found {len(sequences)} sequences")
    
    for name, seq in sequences.items():
        print(f"  - {name}: {len(seq)} bp")
    
    if len(sequences) < 2:
        print("\nWarning: Need at least 2 sequences to find common regions")
        return
    
    # Find common regions
    print(f"\nAnalyzing common sequences (min length: {args.min_common_length} bp)...")
    common_regions = find_common_substrings(sequences, min_length=args.min_common_length)
    
    if not common_regions:
        print("No common sequences found!")
        return
    
    print(f"\nFound {len(common_regions)} common region(s):")
    for i, (substring, pos1, pos2, seq_names) in enumerate(common_regions, 1):
        print(f"\n  Common region {i}:")
        print(f"    Length: {len(substring)} bp")
        print(f"    Sequences: {', '.join(sorted(seq_names))}")
        print(f"    Preview: {substring[:50]}..." if len(substring) > 50 else f"    Sequence: {substring}")
    
    # Calculate statistics
    total_common = sum(len(substring) for substring, _, _, _ in common_regions)
    total_length = sum(len(seq) for seq in sequences.values())
    avg_length = total_length / len(sequences)
    common_percent = (total_common / avg_length) * 100 if avg_length > 0 else 0
    
    print(f"\nStatistics:")
    print(f"  Total common sequence length: {total_common} bp")
    print(f"  Average sequence length: {avg_length:.1f} bp")
    print(f"  Common sequence percentage: {common_percent:.1f}%")
    
    # Find unique regions
    print(f"\nFinding unique regions (min length: {args.min_unique_length} bp)...")
    unique_regions = find_unique_regions(sequences, common_regions, args.min_unique_length)
    
    print("\nUnique regions per sequence:")
    for seq_name, regions in unique_regions.items():
        if regions:
            total_unique = sum(end - start for start, end in regions)
            print(f"  {seq_name}: {len(regions)} region(s), {total_unique} bp total")
            for start, end in regions:
                print(f"    [{start}:{end}] ({end-start} bp)")
        else:
            print(f"  {seq_name}: No unique regions found (all common)")
    
    # Create output with common sequences removed
    if not args.analyze_only and args.output:
        print(f"\nCreating output file with common sequences removed: {args.output}")
        unique_sequences = remove_common_sequences(sequences, args.min_unique_length)
        write_fasta(unique_sequences, args.output)
        
        print("\nOutput sequences:")
        for name, seq in unique_sequences.items():
            print(f"  - {name}: {len(seq)} bp (original: {len(sequences[name])} bp)")
        
        print(f"\n✓ Saved to {args.output}")
        print("\nRecommendation: Test alignment with this file and compare:")
        print("  - MAPQ scores (should be higher)")
        print("  - Multi-mapping rate (should be lower)")
        print("  - Alignment rate (may be slightly lower if reads span common regions)")


if __name__ == "__main__":
    main()


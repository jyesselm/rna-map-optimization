#!/usr/bin/env python3
"""View alignments in a human-readable format.

This script displays SAM alignments in an easy-to-read format showing:
- Reference sequence
- Read sequence aligned to reference
- Matches, mismatches, insertions, deletions
- Quality scores and MAPQ
- Position information
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.text import Text
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False
    print("Note: Install 'rich' for colored output: pip install rich")


def parse_cigar(cigar: str, read_seq: str, ref_seq: str, ref_start: int) -> Tuple[str, str, str]:
    """Parse CIGAR string and create aligned sequences.
    
    Returns:
        (aligned_ref, aligned_read, match_string)
    """
    if cigar == "*":
        return "", "", ""
    
    # Parse CIGAR
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    
    ref_pos = ref_start - 1  # Convert to 0-based
    read_pos = 0
    aligned_ref = []
    aligned_read = []
    match_string = []
    
    for length, op in cigar_ops:
        length = int(length)
        
        if op == 'M':  # Match or mismatch
            for i in range(length):
                if ref_pos < len(ref_seq) and read_pos < len(read_seq):
                    ref_base = ref_seq[ref_pos]
                    read_base = read_seq[read_pos]
                    aligned_ref.append(ref_base)
                    aligned_read.append(read_base)
                    if ref_base == read_base:
                        match_string.append('|')
                    else:
                        match_string.append('X')
                    ref_pos += 1
                    read_pos += 1
                else:
                    break
                    
        elif op == 'I':  # Insertion in read
            for i in range(length):
                if read_pos < len(read_seq):
                    aligned_ref.append('-')
                    aligned_read.append(read_seq[read_pos])
                    match_string.append(' ')
                    read_pos += 1
                    
        elif op == 'D':  # Deletion in read
            for i in range(length):
                if ref_pos < len(ref_seq):
                    aligned_ref.append(ref_seq[ref_pos])
                    aligned_read.append('-')
                    match_string.append(' ')
                    ref_pos += 1
                    
        elif op == 'N':  # Skipped region (intron)
            for i in range(length):
                aligned_ref.append('N')
                aligned_read.append('-')
                match_string.append(' ')
                ref_pos += 1
                
        elif op == 'S':  # Soft clip (part of read not aligned)
            for i in range(length):
                if read_pos < len(read_seq):
                    aligned_ref.append(' ')
                    aligned_read.append(read_seq[read_pos].lower())
                    match_string.append(' ')
                    read_pos += 1
                    
        elif op == 'H':  # Hard clip (not in read sequence)
            pass  # Already clipped from sequence
            
        elif op == '=':  # Match
            for i in range(length):
                if ref_pos < len(ref_seq) and read_pos < len(read_seq):
                    ref_base = ref_seq[ref_pos]
                    read_base = read_seq[read_pos]
                    aligned_ref.append(ref_base)
                    aligned_read.append(read_base)
                    match_string.append('|')
                    ref_pos += 1
                    read_pos += 1
                    
        elif op == 'X':  # Mismatch
            for i in range(length):
                if ref_pos < len(ref_seq) and read_pos < len(read_seq):
                    ref_base = ref_seq[ref_pos]
                    read_base = read_seq[read_pos]
                    aligned_ref.append(ref_base)
                    aligned_read.append(read_base)
                    match_string.append('X')
                    ref_pos += 1
                    read_pos += 1
    
    return ''.join(aligned_ref), ''.join(aligned_read), ''.join(match_string)


def parse_sam_line(line: str) -> Optional[Dict]:
    """Parse a SAM line."""
    fields = line.strip().split("\t")
    if len(fields) < 11 or fields[0].startswith("@"):
        return None
    
    # Extract optional tags
    tags = {}
    for field in fields[11:]:
        if ':' in field:
            tag, tag_type, value = field.split(':', 2)
            tags[tag] = value
    
    return {
        "qname": fields[0],
        "flag": int(fields[1]),
        "rname": fields[2],
        "pos": int(fields[3]) if fields[3] != "0" else None,
        "mapq": int(fields[4]),
        "cigar": fields[5],
        "rnext": fields[6],
        "pnext": int(fields[7]) if fields[7] != "0" else None,
        "tlen": int(fields[8]) if fields[8] != "0" else None,
        "seq": fields[9],
        "qual": fields[10],
        "tags": tags,
    }


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


def format_sequence(seq: str, match_str: str, use_color: bool = True) -> str:
    """Format sequence with colors/symbols for matches/mismatches."""
    if not use_color or not RICH_AVAILABLE:
        return seq
    
    from rich.text import Text
    text = Text()
    for i, (base, match) in enumerate(zip(seq, match_str)):
        if match == '|':
            text.append(base, style="green")
        elif match == 'X':
            text.append(base, style="red")
        elif base == '-':
            text.append(base, style="yellow")
        elif base.islower():
            text.append(base, style="dim")
        else:
            text.append(base)
    return text


def display_alignment_rich(
    aln: Dict,
    ref_seq: str,
    console: Console,
    max_width: int = 80,
    show_quality: bool = True,
):
    """Display alignment using rich formatting."""
    from rich.panel import Panel
    from rich.text import Text
    
    # Parse CIGAR and create alignment
    ref_start = aln["pos"] - 1 if aln["pos"] else 0
    aligned_ref, aligned_read, match_str = parse_cigar(
        aln["cigar"], aln["seq"], ref_seq, ref_start
    )
    
    # Calculate statistics
    matches = match_str.count('|')
    mismatches = match_str.count('X')
    insertions = aligned_read.count('-')
    deletions = aligned_ref.count('-')
    total = len([x for x in match_str if x in '|X'])
    identity = (matches / total * 100) if total > 0 else 0.0
    
    # Create header
    header = Text()
    header.append(f"Read: {aln['qname']}", style="bold cyan")
    header.append(f" | Reference: {aln['rname']}", style="bold blue")
    header.append(f" | Position: {aln['pos']}", style="bold")
    header.append(f" | MAPQ: {aln['mapq']}", style="bold yellow")
    
    # Create alignment display
    lines = []
    
    # Reference line with position numbers
    ref_line = Text()
    ref_line.append(f"Ref  {aln['pos']:>6} ", style="dim")
    for i, (base, match) in enumerate(zip(aligned_ref, match_str)):
        if i > 0 and i % 10 == 0:
            ref_line.append(f" {i+aln['pos']:<6}", style="dim")
        if match == '|':
            ref_line.append(base, style="green")
        elif match == 'X':
            ref_line.append(base, style="red")
        elif base == '-':
            ref_line.append(base, style="yellow")
        else:
            ref_line.append(base)
    lines.append(ref_line)
    
    # Match line
    match_line = Text()
    match_line.append("      " + " " * len(str(aln['pos'])))
    for i, match in enumerate(match_str):
        if i > 0 and i % 10 == 0:
            match_line.append(" " * 7)
        if match == '|':
            match_line.append("|", style="green")
        elif match == 'X':
            match_line.append("X", style="red")
        else:
            match_line.append(" ")
    lines.append(match_line)
    
    # Read line
    read_line = Text()
    read_line.append("Read     ")
    for i, (base, match) in enumerate(zip(aligned_read, match_str)):
        if i > 0 and i % 10 == 0:
            read_line.append(" " * 7)
        if match == '|':
            read_line.append(base, style="green")
        elif match == 'X':
            read_line.append(base, style="red")
        elif base == '-':
            read_line.append(base, style="yellow")
        elif base.islower():
            read_line.append(base, style="dim")
        else:
            read_line.append(base)
    lines.append(read_line)
    
    # Quality line (if requested)
    if show_quality and aln["qual"]:
        qual_line = Text()
        qual_line.append("Qual    ")
        read_pos = 0
        for i, match in enumerate(match_str):
            if i > 0 and i % 10 == 0:
                qual_line.append(" " * 7)
            if match in '|X ' and read_pos < len(aln["qual"]):
                q = ord(aln["qual"][read_pos]) - 33
                if q >= 30:
                    qual_line.append(str(q % 10), style="green")
                elif q >= 20:
                    qual_line.append(str(q % 10), style="yellow")
                else:
                    qual_line.append(str(q % 10), style="red")
                read_pos += 1
            elif aligned_read[i] != '-':
                read_pos += 1
            else:
                qual_line.append(" ")
        lines.append(qual_line)
    
    # Statistics
    stats = Text()
    stats.append(f"Identity: {identity:.1f}%", style="bold")
    stats.append(f" | Matches: {matches}", style="green")
    stats.append(f" | Mismatches: {mismatches}", style="red")
    if insertions > 0:
        stats.append(f" | Insertions: {insertions}", style="yellow")
    if deletions > 0:
        stats.append(f" | Deletions: {deletions}", style="yellow")
    
    # Display
    console.print()
    console.print(Panel(header, border_style="cyan"))
    for line in lines:
        console.print(line)
    console.print(stats)
    console.print()


def create_aligned_read_string(
    aln: Dict,
    ref_seq: str,
    ref_start_pos: int = 1,
) -> Tuple[str, str, Dict]:
    """Create aligned read string that matches reference positions.
    
    Returns:
        (aligned_read_str, highlight_str, stats)
        aligned_read_str: Read sequence aligned to reference positions (full length of ref_seq)
        highlight_str: String indicating matches (space) vs mismatches (!) vs deletions (-)
    """
    # Parse CIGAR and create alignment
    ref_start = aln["pos"] - 1 if aln["pos"] else 0
    aligned_ref, aligned_read, match_str = parse_cigar(
        aln["cigar"], aln["seq"], ref_seq, ref_start
    )
    
    # Create a full-length string matching reference positions
    aligned_read_str = [' '] * len(ref_seq)  # Initialize with spaces
    highlight_str = [' '] * len(ref_seq)  # Initialize with spaces
    
    # Map aligned positions to reference positions
    ref_pos_in_aligned = 0
    for i, (ref_base, read_base, match) in enumerate(zip(aligned_ref, aligned_read, match_str)):
        if ref_base != '-':  # This position exists in reference
            ref_idx = ref_start + ref_pos_in_aligned
            if ref_idx < len(ref_seq):
                if read_base == '-':
                    # Deletion: show dash
                    aligned_read_str[ref_idx] = '-'
                    highlight_str[ref_idx] = '-'  # Deletion marker
                else:
                    # Match or mismatch
                    aligned_read_str[ref_idx] = read_base
                    if match == '|':
                        highlight_str[ref_idx] = ' '  # Match (no highlight)
                    elif match == 'X':
                        highlight_str[ref_idx] = '!'  # Mismatch
                ref_pos_in_aligned += 1
        elif read_base != '-' and read_base != ' ':  # Insertion (not in reference)
            # For insertions, we could show them but they don't map to ref positions
            # Skip for now - they're shown as lowercase in the read
            pass
    
    # Calculate statistics
    matches = match_str.count('|')
    mismatches = match_str.count('X')
    insertions = sum(1 for r in aligned_read if r == '-')
    deletions = sum(1 for r in aligned_read_str if r == '-')
    total = len([x for x in match_str if x in '|X'])
    identity = (matches / total * 100) if total > 0 else 0.0
    
    stats = {
        "matches": matches,
        "mismatches": mismatches,
        "insertions": insertions,
        "deletions": deletions,
        "identity": identity,
    }
    
    return ''.join(aligned_read_str), ''.join(highlight_str), stats


def display_alignment_compact_rich(
    ref_seq: str,
    alignments: List[Dict],
    console,
    ref_name: str = "",
    show_quality: bool = False,
):
    """Display reference as one line, then each read aligned below with highlights."""
    from rich.text import Text
    from rich.panel import Panel
    
    # Display reference sequence
    ref_text = Text()
    ref_text.append("Ref: ", style="bold blue")
    for i, base in enumerate(ref_seq):
        if i > 0 and i % 10 == 0:
            ref_text.append(f" {i:<3}", style="dim")
        ref_text.append(base, style="bold")
    console.print(ref_text)
    console.print()
    
    # Display each read
    for aln in alignments:
        # Create aligned read string
        aligned_read_str, highlight_str, stats = create_aligned_read_string(
            aln, ref_seq, aln["pos"] if aln["pos"] else 1
        )
        
        # Create read display
        read_text = Text()
        read_text.append(f"{aln['qname'][:30]:<30} ", style="dim")
        read_text.append(f"pos:{aln['pos']:<4} ", style="dim")
        read_text.append(f"MAPQ:{aln['mapq']:<3} ", style="dim")
        
        # Add aligned sequence with highlights (already aligned to ref positions)
        for i in range(len(ref_seq)):
            base = aligned_read_str[i] if i < len(aligned_read_str) else ' '
            highlight = highlight_str[i] if i < len(highlight_str) else ' '
            
            if base == ' ':
                read_text.append(" ", style="dim")
            elif highlight == '!':
                read_text.append(base, style="bold red")  # Mismatch
            elif highlight == '-':
                read_text.append("-", style="bold yellow")  # Deletion
            elif base.islower():
                read_text.append(base, style="dim")  # Soft clip
            else:
                read_text.append(base, style="green")  # Match
        
        # Add statistics
        read_text.append(f" | id:{stats['identity']:.1f}%", style="dim")
        
        console.print(read_text)
    
    console.print()


def display_alignment_compact_plain(
    ref_seq: str,
    alignments: List[Dict],
    ref_name: str = "",
    show_quality: bool = False,
):
    """Display reference as one line, then each read aligned below with highlights."""
    # Display reference sequence
    print("Ref: ", end="")
    for i, base in enumerate(ref_seq):
        if i > 0 and i % 10 == 0:
            print(f" {i:<3}", end="")
        print(base, end="")
    print()
    print()
    
    # Display each read
    for aln in alignments:
        # Create aligned read string
        aligned_read_str, highlight_str, stats = create_aligned_read_string(
            aln, ref_seq, aln["pos"] if aln["pos"] else 1
        )
        
        # Print read header
        print(f"{aln['qname'][:30]:<30} pos:{aln['pos']:<4} MAPQ:{aln['mapq']:<3} ", end="")
        
        # Print aligned sequence with highlights (already aligned to ref positions)
        for i in range(len(ref_seq)):
            base = aligned_read_str[i] if i < len(aligned_read_str) else ' '
            highlight = highlight_str[i] if i < len(highlight_str) else ' '
            
            if base == ' ':
                print(" ", end="")
            elif highlight == '!':
                print(f"\033[91m{base}\033[0m", end="")  # Red for mismatch
            elif highlight == '-':
                print(f"\033[93m-\033[0m", end="")  # Yellow for deletion
            elif base.islower():
                print(f"\033[2m{base}\033[0m", end="")  # Dim for soft clip
            else:
                print(f"\033[92m{base}\033[0m", end="")  # Green for match
        
        # Add statistics
        print(f" | id:{stats['identity']:.1f}%")
    
    print()


def display_alignment_plain(
    aln: Dict,
    ref_seq: str,
    max_width: int = 80,
    show_quality: bool = True,
):
    """Display alignment in plain text format (legacy - kept for compatibility)."""
    # This is the old format - we'll use compact format by default now
    pass


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="View alignments in human-readable format"
    )
    
    parser.add_argument(
        "--sam",
        type=Path,
        required=True,
        help="Path to SAM file",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Path to reference FASTA file",
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=10,
        help="Maximum number of alignments to display (default: 10)",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ to display (default: 0)",
    )
    parser.add_argument(
        "--read-name",
        type=str,
        help="Display specific read by name",
    )
    parser.add_argument(
        "--no-quality",
        action="store_true",
        help="Don't show quality scores",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Save output to file instead of printing",
    )
    
    args = parser.parse_args()
    
    # Read reference sequences
    ref_seqs = read_fasta(args.fasta)
    
    # Setup console
    if RICH_AVAILABLE:
        console = Console(record=args.output is not None, width=120)
    else:
        console = None
    
    # Collect alignments grouped by reference
    alignments_by_ref = {}
    
    with open(args.sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            
            aln = parse_sam_line(line)
            if not aln:
                continue
            
            # Filter by MAPQ
            if aln["mapq"] < args.min_mapq:
                continue
            
            # Filter by read name
            if args.read_name and aln["qname"] != args.read_name:
                continue
            
            # Check if reference exists
            if aln["rname"] not in ref_seqs:
                continue
            
            if aln["pos"] is None:
                continue
            
            # Group by reference
            ref_name = aln["rname"]
            if ref_name not in alignments_by_ref:
                alignments_by_ref[ref_name] = []
            
            alignments_by_ref[ref_name].append(aln)
    
    if not alignments_by_ref:
        print("No alignments found matching criteria.")
        return 1
    
    # Display each reference with its alignments
    for ref_name, alignments in alignments_by_ref.items():
        # Limit number of reads
        alignments = alignments[:args.max_reads]
        ref_seq = ref_seqs[ref_name]
        
        if RICH_AVAILABLE and console:
            display_alignment_compact_rich(
                ref_seq,
                alignments,
                console,
                ref_name=ref_name,
                show_quality=not args.no_quality,
            )
        else:
            display_alignment_compact_plain(
                ref_seq,
                alignments,
                ref_name=ref_name,
                show_quality=not args.no_quality,
            )
    
    # Save if requested
    if args.output and RICH_AVAILABLE and console:
        if str(args.output).endswith('.html'):
            console.save_html(str(args.output))
        else:
            # Save as text
            with open(args.output, 'w') as f:
                f.write(console.export_text())
        print(f"\nOutput saved to {args.output}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


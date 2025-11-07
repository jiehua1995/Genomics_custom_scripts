#!/usr/bin/env python3
"""
quick_qc.py - Quality control and visualization for FASTA/FASTQ datasets

This standalone script accepts a FASTA/FASTQ file (optionally gzipped) or a
folder of such files and produces a self-contained interactive HTML report
using Plotly.

Features:
- Auto-detect FASTA vs FASTQ and nucleotide vs protein where possible
- Streaming/chunked parsing to handle large files
- Multi-file support with per-file coloring
- Multi-threaded parsing across files for speed
- Outputs a single self-contained HTML report: {input_name}_QC_report.html

Dependencies: plotly, pandas, biopython, tqdm

All logs and messages are printed in English.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import os
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple

import plotly.graph_objects as go
import plotly.offline as pyo
import plotly.io as pio
from Bio import SeqIO
from tqdm import tqdm
import pandas as pd

# Add colorama for colored progress bars
try:
    from colorama import init as colorama_init, Fore, Back, Style
    colorama_init()
    HAS_COLORAMA = True
except ImportError:
    HAS_COLORAMA = False
    # Define dummy constants if colorama not available
    class Fore:
        GREEN = RED = YELLOW = CYAN = RESET = ""
    class Back:
        GREEN = RED = YELLOW = CYAN = RESET = ""
    class Style:
        BRIGHT = RESET_ALL = ""
import plotly.io as pio
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Constants
MAX_POS_QUALITY_PLOT = 500  # max per-base position to collect quality for plotting


def open_maybe_gzip(path: str):
    """Open a file with transparent gzip support."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "r", errors="ignore")


def detect_format_and_type(path: str, explicit_type: str = None) -> Tuple[str, str]:
    """Detect file format (fasta/fastq) and sequence type (dna/rna/protein).

    Returns (format, seq_type). format in {"fasta","fastq"}. seq_type in {"dna","rna","protein","unknown"}.
    """
    # Try to infer from extension
    lname = path.lower()
    fmt = None
    if lname.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        fmt = 'fasta'
    if lname.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
        fmt = 'fastq'

    # Peek a few records
    try:
        with open_maybe_gzip(path) as fh:
            # choose parser based on extension guess first
            if fmt == 'fastq':
                recs = list(SeqIO.parse(fh, 'fastq'))[:50]
            elif fmt == 'fasta':
                recs = list(SeqIO.parse(fh, 'fasta'))[:50]
            else:
                # try both, prefer fastq if it yields records
                fh.seek(0)
                recs = list(SeqIO.parse(fh, 'fastq'))[:50]
                if not recs:
                    fh.seek(0)
                    recs = list(SeqIO.parse(fh, 'fasta'))[:50]
            if not recs:
                return (fmt or 'fasta', 'unknown')

            # detect sequence alphabet
            seqs = [str(r.seq).upper() for r in recs]
            letters = Counter(''.join(seqs))
            # simple heuristics
            nuc_score = sum(letters.get(x, 0) for x in 'ATGCU')
            aa_score = sum(letters.get(x, 0) for x in 'EFILPQZ')  # Z uncommon in nuc
            seq_type = 'protein' if aa_score > nuc_score else 'dna'
            if explicit_type:
                seq_type = explicit_type
            # determine format more robustly: if any record has .letter_annotations for phred
            if fmt is None:
                # reuse parser detection
                # if recs came from fastq parse earlier this would be fine
                # we conservatively set fasta if uncertain
                fmt = 'fastq' if isinstance(recs[0].letter_annotations, dict) and recs[0].letter_annotations else 'fasta'
            return (fmt, seq_type)
    except Exception:
        # fallback
        return (fmt or 'fasta', explicit_type or 'unknown')


def process_file(path: str, seq_type_arg: str = None) -> Dict:
    """Stream a single file, compute statistics and minimal collections for plotting.

    Returns a dict of summary statistics and small arrays for plotting.
    """
    fmt, detected_type = detect_format_and_type(path, seq_type_arg)
    results = {
        'path': path,
        'format': fmt,
        'type': detected_type,
        'num_seqs': 0,
        'total_bases': 0,
        'lengths': [],
        'gc_pct': [],
        'aa_counts': Counter(),
        'qual_pos_sums': defaultdict(int),
        'qual_pos_counts': defaultdict(int),
        'qual_scores': [],
    }

    # Stream parsing
    open_func = open_maybe_gzip
    with open_func(path) as fh:
        parser = 'fastq' if fmt == 'fastq' else 'fasta'
        for rec in SeqIO.parse(fh, parser):
            seq = str(rec.seq)
            seqlen = len(seq)
            results['num_seqs'] += 1
            results['total_bases'] += seqlen
            results['lengths'].append(seqlen)
            if detected_type in ('dna', 'rna', 'unknown'):
                # compute GC%
                seq_upper = seq.upper()
                g = seq_upper.count('G')
                c = seq_upper.count('C')
                gc = (g + c) / seqlen * 100 if seqlen > 0 else 0
                results['gc_pct'].append(gc)
            else:
                # amino acid composition
                results['aa_counts'].update(list(seq.upper()))

            # Qualities for FASTQ
            if parser == 'fastq' and hasattr(rec, 'letter_annotations'):
                quals = rec.letter_annotations.get('phred_quality', [])
                if quals:
                    results['qual_scores'].extend(quals)
                    for i, q in enumerate(quals[:MAX_POS_QUALITY_PLOT]):
                        results['qual_pos_sums'][i] += q
                        results['qual_pos_counts'][i] += 1

    # Post-process aa composition to fraction
    if results['aa_counts']:
        total_aa = sum(results['aa_counts'].values())
        results['aa_comp'] = {k: v / total_aa * 100 for k, v in results['aa_counts'].items()}
    else:
        results['aa_comp'] = {}

    return results


def process_file_with_progress(path: str, seq_type_arg: str = None, progress_bar=None) -> Dict:
    """Stream a single file with progress tracking, compute statistics and minimal collections for plotting.

    Returns a dict of summary statistics and small arrays for plotting.
    """
    fmt, detected_type = detect_format_and_type(path, seq_type_arg)
    results = {
        'path': path,
        'format': fmt,
        'type': detected_type,
        'num_seqs': 0,
        'total_bases': 0,
        'lengths': [],
        'gc_pct': [],
        'aa_counts': Counter(),
        'qual_pos_sums': defaultdict(int),
        'qual_pos_counts': defaultdict(int),
        'qual_scores': [],
    }

    # Get file size to estimate progress
    try:
        file_size = os.path.getsize(path)
    except:
        file_size = 0
    
    bytes_read = 0
    update_interval = max(1, file_size // 100)  # Update progress every 1% of file
    last_update = 0

    # Stream parsing with progress tracking
    open_func = open_maybe_gzip
    with open_func(path) as fh:
        parser = 'fastq' if fmt == 'fastq' else 'fasta'
        
        for rec in SeqIO.parse(fh, parser):
            seq = str(rec.seq)
            seqlen = len(seq)
            results['num_seqs'] += 1
            results['total_bases'] += seqlen
            results['lengths'].append(seqlen)
            
            if detected_type in ('dna', 'rna', 'unknown'):
                # compute GC%
                seq_upper = seq.upper()
                g = seq_upper.count('G')
                c = seq_upper.count('C')
                gc = (g + c) / seqlen * 100 if seqlen > 0 else 0
                results['gc_pct'].append(gc)
            else:
                # amino acid composition
                results['aa_counts'].update(list(seq.upper()))

            # Qualities for FASTQ
            if parser == 'fastq' and hasattr(rec, 'letter_annotations'):
                quals = rec.letter_annotations.get('phred_quality', [])
                if quals:
                    results['qual_scores'].extend(quals)
                    for i, q in enumerate(quals[:MAX_POS_QUALITY_PLOT]):
                        results['qual_pos_sums'][i] += q
                        results['qual_pos_counts'][i] += 1

            # Update progress based on estimated file position
            if progress_bar and file_size > 0:
                try:
                    current_pos = fh.tell()
                    if current_pos - last_update > update_interval:
                        progress_pct = min(100, int((current_pos / file_size) * 100))
                        progress_bar.n = progress_pct
                        progress_bar.total = 100
                        progress_bar.refresh()
                        last_update = current_pos
                except:
                    pass  # Some file objects don't support tell()

    # Final progress update
    if progress_bar:
        progress_bar.n = 100
        progress_bar.total = 100
        progress_bar.refresh()

    # Post-process aa composition to fraction
    if results['aa_counts']:
        total_aa = sum(results['aa_counts'].values())
        results['aa_comp'] = {k: v / total_aa * 100 for k, v in results['aa_counts'].items()}
    else:
        results['aa_comp'] = {}

    return results


def make_report(results_list: List[Dict], out_html: str):
    """Generate an interactive HTML report with Plotly and save to out_html.

    The report is self-contained (includes plotly.js) and can be opened offline.
    """
    # Build summary table
    rows = []
    for r in results_list:
        n50 = calculate_nx(sorted(r['lengths'], reverse=True), 50) if r['lengths'] else 0
        rows.append({
            'file': os.path.basename(r['path']),
            'format': r['format'],
            'type': r['type'],
            'num_seqs': r['num_seqs'],
            'total_bases': r['total_bases'],
            'mean_length': int(pd.Series(r['lengths']).mean()) if r['lengths'] else 0,
            'N50': int(n50),
            'mean_GC': float(pd.Series(r['gc_pct']).mean()) if r['gc_pct'] else None,
        })
    df = pd.DataFrame(rows)

    # Create figures
    figs = []

    # Length distribution overlayed with improved binning
    fig_len = go.Figure()
    colors = px_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
    
    # Calculate intelligent bins for length distribution
    all_lengths = []
    for r in results_list:
        all_lengths.extend(r['lengths'])
    
    if all_lengths:
        min_len, max_len = min(all_lengths), max(all_lengths)
        # Use more bins for better resolution
        if max_len - min_len > 1000:
            # For large ranges, use sqrt of data points or 100, whichever is larger
            nbins = max(100, int(len(all_lengths) ** 0.5))
        else:
            # For smaller ranges, use 50 bins minimum
            nbins = max(50, (max_len - min_len) // 10)
        
        # Limit bins to reasonable maximum
        nbins = min(nbins, 200)
    else:
        nbins = 50
        
    for i, r in enumerate(results_list):
        if not r['lengths']:
            continue
        fig_len.add_trace(go.Histogram(
            x=r['lengths'], 
            name=os.path.basename(r['path']), 
            opacity=0.75,
            nbinsx=nbins
        ))
    fig_len.update_layout(
        title='Sequence/Read Length Distribution', 
        barmode='overlay', 
        xaxis_title='Length (bp)', 
        yaxis_title='Count'
    )
    figs.append(('Length distribution', fig_len))

    # GC% or AA composition
    # If nucleotide datasets
    if any(r['gc_pct'] for r in results_list):
        fig_gc = go.Figure()
        
        # Calculate intelligent bins for GC distribution (0-100%)
        all_gc = []
        for r in results_list:
            all_gc.extend(r['gc_pct'])
        
        # Use 100 bins for GC% to get 1% resolution, or adaptive based on data spread
        if all_gc:
            gc_range = max(all_gc) - min(all_gc)
            if gc_range > 50:
                nbins_gc = 100  # 1% resolution for wide distributions
            elif gc_range > 20:
                nbins_gc = 50   # 0.5% resolution for medium distributions  
            else:
                nbins_gc = max(20, int(gc_range * 2))  # Higher resolution for narrow distributions
        else:
            nbins_gc = 50
            
        for r in results_list:
            if not r['gc_pct']:
                continue
            fig_gc.add_trace(go.Histogram(
                x=r['gc_pct'], 
                name=os.path.basename(r['path']), 
                opacity=0.75,
                nbinsx=nbins_gc
            ))
        fig_gc.update_layout(
            title='GC% Distribution', 
            barmode='overlay', 
            xaxis_title='GC%',
            yaxis_title='Count'
        )
        figs.append(('GC content', fig_gc))
    else:
        # amino acid composition per file as stacked bars
        # collect common amino acids
        aa_set = set()
        for r in results_list:
            aa_set |= set(r['aa_comp'].keys())
        aa_list = sorted(aa_set)
        fig_aa = go.Figure()
        for r in results_list:
            vals = [r['aa_comp'].get(aa, 0) for aa in aa_list]
            fig_aa.add_trace(go.Bar(name=os.path.basename(r['path']), x=aa_list, y=vals))
        fig_aa.update_layout(title='Amino Acid Composition (%)', barmode='group', xaxis_title='Amino Acid')
        figs.append(('AA composition', fig_aa))

    # FASTQ quality plots
    if any(r['format'] == 'fastq' for r in results_list):
        # overall quality score histogram per file with improved binning
        fig_q = go.Figure()
        
        # Calculate intelligent bins for quality scores (typically 0-40+ range)
        all_quals = []
        for r in results_list:
            all_quals.extend(r['qual_scores'])
        
        if all_quals:
            min_q, max_q = min(all_quals), max(all_quals)
            # Quality scores are typically integers, so use 1-unit bins if range is reasonable
            if max_q - min_q <= 50:
                # Use 1-unit bins for quality scores
                nbins_q = max(max_q - min_q + 1, 20)
            else:
                # For very wide ranges, use 50 bins
                nbins_q = 50
        else:
            nbins_q = 40  # Default for typical quality range 0-40
        
        for r in results_list:
            if not r['qual_scores']:
                continue
            fig_q.add_trace(go.Histogram(
                x=r['qual_scores'], 
                name=os.path.basename(r['path']), 
                opacity=0.75,
                nbinsx=nbins_q
            ))
        fig_q.update_layout(
            title='Per-base Quality Score Distribution (all bases)', 
            barmode='overlay', 
            xaxis_title='Phred quality',
            yaxis_title='Count'
        )
        figs.append(('Quality distribution', fig_q))

        # mean quality by position (limited to MAX_POS_QUALITY_PLOT)
        fig_qpos = go.Figure()
        for r in results_list:
            if not r['qual_pos_counts']:
                continue
            positions = sorted(r['qual_pos_counts'].keys())
            means = [r['qual_pos_sums'][p] / r['qual_pos_counts'][p] for p in positions]
            fig_qpos.add_trace(go.Scatter(x=positions, y=means, mode='lines+markers', name=os.path.basename(r['path'])))
        fig_qpos.update_layout(title='Mean Quality by Read Position', xaxis_title='Position', yaxis_title='Mean Phred quality')
        figs.append(('Mean quality by position', fig_qpos))

    # Compose modern HTML with enhanced styling
    modern_css = """
    <style>
        /* Modern CSS styling inspired by MultiQC but with contemporary design */
        :root {
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --accent-color: #e74c3c;
            --success-color: #27ae60;
            --warning-color: #f39c12;
            --background-color: #f8f9fa;
            --card-background: #ffffff;
            --text-color: #2c3e50;
            --border-color: #dee2e6;
            --shadow: 0 0.125rem 0.25rem rgba(0, 0, 0, 0.075);
            --shadow-lg: 0 0.5rem 1rem rgba(0, 0, 0, 0.15);
        }

        * {
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Roboto', 'Oxygen', 'Ubuntu', 'Cantarell', 'Open Sans', 'Helvetica Neue', sans-serif;
            line-height: 1.6;
            color: var(--text-color);
            background-color: var(--background-color);
            margin: 0;
            padding: 0;
            font-size: 14px;
        }

        .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }

        .header {
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
            color: white;
            padding: 30px 0;
            margin-bottom: 30px;
            box-shadow: var(--shadow-lg);
        }

        .header-content {
            max-width: 1200px;
            margin: 0 auto;
            padding: 0 20px;
        }

        .header h1 {
            margin: 0 0 10px 0;
            font-size: 2.5rem;
            font-weight: 300;
            text-shadow: 0 2px 4px rgba(0,0,0,0.3);
        }

        .header .subtitle {
            margin: 0;
            opacity: 0.9;
            font-size: 1.1rem;
            font-weight: 400;
        }

        .summary-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }

        .card {
            background: var(--card-background);
            border-radius: 8px;
            padding: 24px;
            box-shadow: var(--shadow);
            border: 1px solid var(--border-color);
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }

        .card:hover {
            transform: translateY(-2px);
            box-shadow: var(--shadow-lg);
        }

        .card-title {
            font-size: 0.875rem;
            color: #6c757d;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            margin: 0 0 8px 0;
            font-weight: 600;
        }

        .card-value {
            font-size: 2rem;
            font-weight: 700;
            color: var(--primary-color);
            margin: 0;
            line-height: 1.2;
        }

        .card-subtitle {
            font-size: 0.875rem;
            color: #6c757d;
            margin: 4px 0 0 0;
        }

        .files-table-container {
            background: var(--card-background);
            border-radius: 8px;
            padding: 24px;
            margin-bottom: 30px;
            box-shadow: var(--shadow);
            border: 1px solid var(--border-color);
        }

        .section-title {
            font-size: 1.5rem;
            font-weight: 600;
            color: var(--primary-color);
            margin: 0 0 20px 0;
            padding-bottom: 10px;
            border-bottom: 2px solid var(--secondary-color);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.875rem;
        }

        th {
            background: var(--background-color);
            font-weight: 600;
            text-align: left;
            padding: 12px 16px;
            color: var(--text-color);
            border-bottom: 2px solid var(--border-color);
        }

        td {
            padding: 12px 16px;
            border-bottom: 1px solid var(--border-color);
        }

        tr:hover {
            background: #f8f9fa;
        }

        .plot-container {
            background: var(--card-background);
            border-radius: 8px;
            padding: 24px;
            margin-bottom: 30px;
            box-shadow: var(--shadow);
            border: 1px solid var(--border-color);
        }

        .plot-title {
            font-size: 1.25rem;
            font-weight: 600;
            color: var(--primary-color);
            margin: 0 0 20px 0;
            display: flex;
            align-items: center;
        }

        .plot-title::before {
            content: '';
            width: 4px;
            height: 20px;
            background: var(--secondary-color);
            margin-right: 12px;
            border-radius: 2px;
        }

        .badge {
            display: inline-block;
            padding: 4px 8px;
            border-radius: 12px;
            font-size: 0.75rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.3px;
        }

        .badge-dna {
            background: rgba(46, 204, 113, 0.1);
            color: var(--success-color);
        }

        .badge-rna {
            background: rgba(52, 152, 219, 0.1);
            color: var(--secondary-color);
        }

        .badge-protein {
            background: rgba(231, 76, 60, 0.1);
            color: var(--accent-color);
        }

        .badge-fastq {
            background: rgba(243, 156, 18, 0.1);
            color: var(--warning-color);
        }

        .badge-fasta {
            background: rgba(155, 89, 182, 0.1);
            color: #9b59b6;
        }

        .navigation {
            position: sticky;
            top: 0;
            background: var(--card-background);
            border-bottom: 1px solid var(--border-color);
            padding: 16px 0;
            margin-bottom: 30px;
            z-index: 100;
        }

        .nav-links {
            max-width: 1200px;
            margin: 0 auto;
            padding: 0 20px;
            display: flex;
            gap: 24px;
            flex-wrap: wrap;
        }

        .nav-link {
            color: var(--text-color);
            text-decoration: none;
            font-weight: 500;
            padding: 8px 16px;
            border-radius: 6px;
            transition: all 0.2s ease;
        }

        .nav-link:hover {
            background: var(--background-color);
            color: var(--secondary-color);
        }

        @media (max-width: 768px) {
            .container {
                padding: 10px;
            }
            
            .header h1 {
                font-size: 2rem;
            }
            
            .summary-cards {
                grid-template-columns: 1fr;
            }
            
            .nav-links {
                flex-direction: column;
                gap: 8px;
            }
        }

        /* Plotly plot styling */
        .plotly-graph-div {
            margin: 0 !important;
        }

        /* Loading animation */
        .loading {
            display: inline-block;
            width: 20px;
            height: 20px;
            border: 3px solid rgba(52, 152, 219, 0.3);
            border-radius: 50%;
            border-top-color: var(--secondary-color);
            animation: spin 1s ease-in-out infinite;
        }

        @keyframes spin {
            to { transform: rotate(360deg); }
        }
    </style>
    """

    # Calculate summary statistics
    total_files = len(results_list)
    total_sequences = sum(r['num_seqs'] for r in results_list)
    total_bases = sum(r['total_bases'] for r in results_list)
    file_types = [r['format'].upper() for r in results_list]
    seq_types = [r['type'].upper() for r in results_list]
    
    # Create navigation links
    nav_sections = []
    nav_sections.append('<a href="#summary" class="nav-link">📊 Summary</a>')
    nav_sections.append('<a href="#files" class="nav-link">📁 Files</a>')
    nav_sections.append('<a href="#length" class="nav-link">📏 Length Distribution</a>')
    
    if any(r['gc_pct'] for r in results_list):
        nav_sections.append('<a href="#gc" class="nav-link">🧬 GC Content</a>')
    else:
        nav_sections.append('<a href="#aa" class="nav-link">🧪 Amino Acids</a>')
        
    if any(r['format'] == 'fastq' for r in results_list):
        nav_sections.append('<a href="#quality" class="nav-link">⭐ Quality</a>')
        nav_sections.append('<a href="#position" class="nav-link">📍 Position Quality</a>')

    # Create modern HTML structure
    html_parts = [
        modern_css,
        '<div class="header">',
        '  <div class="header-content">',
        '    <h1>🧬 Genomics QC Report</h1>',
        f'    <p class="subtitle">Analysis of {total_files} file(s) • Generated by quick_qc</p>',
        '  </div>',
        '</div>',
        
        '<div class="navigation">',
        '  <div class="nav-links">',
        '    ' + '\n    '.join(nav_sections),
        '  </div>',
        '</div>',
        
        '<div class="container">',
        
        # Summary cards
        '  <div id="summary" class="summary-cards">',
        '    <div class="card">',
        '      <div class="card-title">Total Files</div>',
        f'      <div class="card-value">{total_files:,}</div>',
        f'      <div class="card-subtitle">{", ".join(set(file_types))}</div>',
        '    </div>',
        '    <div class="card">',
        '      <div class="card-title">Total Sequences</div>',
        f'      <div class="card-value">{total_sequences:,}</div>',
        f'      <div class="card-subtitle">{", ".join(set(seq_types))}</div>',
        '    </div>',
        '    <div class="card">',
        '      <div class="card-title">Total Bases</div>',
        f'      <div class="card-value">{total_bases:,}</div>',
        f'      <div class="card-subtitle">{total_bases / 1e6:.2f} Mbp</div>' if total_bases > 1e6 else f'      <div class="card-subtitle">{total_bases / 1e3:.2f} Kbp</div>',
        '    </div>',
        '  </div>',
        
        # Files table
        '  <div id="files" class="files-table-container">',
        '    <h2 class="section-title">📁 File Details</h2>',
    ]
    
    # Enhanced table with badges
    enhanced_df = df.copy()
    enhanced_df['File'] = enhanced_df['file'].apply(lambda x: f'<strong>{x}</strong>')
    
    # Add badges for file types and sequence types
    badge_html = []
    for r in results_list:
        fmt_badge = f'<span class="badge badge-{r["format"].lower()}">{r["format"].upper()}</span>'
        type_badge = f'<span class="badge badge-{r["type"].lower()}">{r["type"].upper()}</span>'
        badge_html.append(f'{fmt_badge} {type_badge}')
    
    enhanced_df['Type'] = badge_html
    # Remove the original 'file' column and reorder
    cols = ['File', 'Type'] + [col for col in enhanced_df.columns if col not in ['File', 'Type', 'file']]
    enhanced_df = enhanced_df[cols]
    
    html_parts.append(enhanced_df.to_html(index=False, escape=False, table_id='files-table'))
    html_parts.append('  </div>')

    # Add plots with modern styling
    plot_id = 0
    for title, fig in figs:
        plot_id += 1
        
        # Determine section ID and icon
        section_id = 'length' if 'Length' in title else ('gc' if 'GC' in title else ('aa' if 'Amino' in title else ('quality' if 'Quality' in title and 'Position' not in title else 'position')))
        
        # Add modern plot container
        html_parts.extend([
            f'  <div id="{section_id}" class="plot-container">',
            f'    <h2 class="plot-title">{title}</h2>',
        ])
        
        # Update figure styling for modern look
        fig.update_layout(
            font=dict(family="system-ui, -apple-system, sans-serif", size=12),
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            margin=dict(l=60, r=60, t=40, b=60),
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1,
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="rgba(0,0,0,0.1)",
                borderwidth=1
            )
        )
        
        # Add grid and improve axis styling
        fig.update_xaxes(
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(0,0,0,0.1)',
            showline=True,
            linewidth=1,
            linecolor='rgba(0,0,0,0.2)',
            mirror=True
        )
        fig.update_yaxes(
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(0,0,0,0.1)',
            showline=True,
            linewidth=1,
            linecolor='rgba(0,0,0,0.2)',
            mirror=True
        )
        
        # Generate plot HTML with improved settings
        plot_html = pyo.plot(
            fig, 
            output_type='div', 
            include_plotlyjs=False,
            config={
                'displayModeBar': True,
                'modeBarButtonsToRemove': ['lasso2d', 'select2d'],
                'displaylogo': False,
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': f'qc_plot_{plot_id}',
                    'height': 600,
                    'width': 1000,
                    'scale': 2
                }
            }
        )
        
        html_parts.append(f'    {plot_html}')
        html_parts.append('  </div>')

    # Close container
    html_parts.append('</div>')

    # Footer with metadata
    import datetime
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_parts.extend([
        '<div style="background: var(--background-color); padding: 20px; margin-top: 30px; border-top: 1px solid var(--border-color); text-align: center; color: #6c757d; font-size: 0.875rem;">',
        f'  <p>Generated by <strong>quick_qc</strong> on {current_time}</p>',
        '  <p>Report powered by <a href="https://plotly.com/" style="color: var(--secondary-color); text-decoration: none;">Plotly.js</a></p>',
        '</div>'
    ])

    # Generate final HTML with embedded plotly.js for offline viewing
    body = '\n'.join(html_parts)
    
    # Use plotly's built-in method to get the plotly.js library content
    plotly_lib = pyo.get_plotlyjs()
    
    # Create modern, self-contained HTML
    final_html = f"""<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <title>🧬 Genomics QC Report</title>
    <script type='text/javascript'>{plotly_lib}</script>
</head>
<body>
{body}
</body>
</html>"""

    with open(out_html, 'w', encoding='utf-8') as fh:
        fh.write(final_html)


def calculate_nx(lengths: List[int], x: int) -> int:
    """Calculate Nx statistic.
    lengths should be sorted in descending order.
    """
    if not lengths:
        return 0
    threshold = sum(lengths) * (x / 100.0)
    cumulative = 0
    for l in lengths:
        cumulative += l
        if cumulative >= threshold:
            return l
    return 0


def collect_inputs(path_or_dir: str) -> List[str]:
    """Collect input files: either single file or all supported files in a directory."""
    if os.path.isdir(path_or_dir):
        exts = ('.fa', '.fasta', '.fq', '.fastq', '.fa.gz', '.fasta.gz', '.fq.gz', '.fastq.gz')
        files = [os.path.join(path_or_dir, f) for f in os.listdir(path_or_dir) if f.lower().endswith(exts)]
        return sorted(files)
    if os.path.isfile(path_or_dir):
        return [path_or_dir]
    raise FileNotFoundError(path_or_dir)


def parse_args():
    parser = argparse.ArgumentParser(description='Generate QC HTML report for FASTA/FASTQ files')
    parser.add_argument('input', help='[REQUIRED] Input FASTA/FASTQ file or directory containing such files')
    parser.add_argument('--type', choices=['dna', 'rna', 'protein'], help='[OPTIONAL] Specify sequence type when auto-detection is ambiguous')
    parser.add_argument('-o', '--output', help='[OPTIONAL] Output HTML file (default: {input_name}_QC_report.html)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='[OPTIONAL] Number of threads to parse files in parallel (default: 4)')
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        files = collect_inputs(args.input)
    except FileNotFoundError:
        print('Input file or directory not found', file=sys.stderr)
        sys.exit(1)
    if not files:
        print('No supported FASTA/FASTQ files found in directory', file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} file(s) to process")
    
    # First, detect file types and display them
    if HAS_COLORAMA:
        print(f"\n{Fore.CYAN}Detecting file types...{Style.RESET_ALL}")
    else:
        print("\nDetecting file types...")
    file_types = {}
    for f in files:
        try:
            fmt, seq_type = detect_format_and_type(f, args.type)
            file_types[f] = (fmt, seq_type)
            if HAS_COLORAMA:
                print(f"  {Fore.GREEN}{os.path.basename(f)}{Style.RESET_ALL}: {Fore.YELLOW}{fmt.upper()}{Style.RESET_ALL} ({Fore.CYAN}{seq_type}{Style.RESET_ALL})")
            else:
                print(f"  {os.path.basename(f)}: {fmt.upper()} ({seq_type})")
        except Exception as e:
            file_types[f] = ('unknown', 'unknown')
            if HAS_COLORAMA:
                print(f"  {Fore.RED}{os.path.basename(f)}{Style.RESET_ALL}: Error detecting type - {e}", file=sys.stderr)
            else:
                print(f"  {os.path.basename(f)}: Error detecting type - {e}", file=sys.stderr)

    if HAS_COLORAMA:
        print(f"\n{Fore.YELLOW}Processing {len(files)} files...{Style.RESET_ALL}")
    else:
        print(f"\nProcessing {len(files)} files...")
    
    # Process files with detailed progress tracking
    results = []
    if HAS_COLORAMA:
        total_progress = tqdm(total=len(files), desc=f"{Fore.BLUE}Overall Progress{Style.RESET_ALL}", position=0, 
                             bar_format="{desc}: {n}/{total} files ({percentage:3.0f}%)")
    else:
        total_progress = tqdm(total=len(files), desc="Overall Progress", position=0, 
                             bar_format="{desc}: {n}/{total} files ({percentage:3.0f}%)")
    
    for i, file_path in enumerate(files):
        filename = os.path.basename(file_path)
        fmt, seq_type = file_types[file_path]
        
        # Create a nested progress bar for current file with correct total and colors
        if HAS_COLORAMA:
            file_desc = f"{Fore.GREEN}Processing {filename}{Style.RESET_ALL}"
        else:
            file_desc = f"Processing {filename}"
        file_progress = tqdm(total=100, desc=file_desc, position=1, leave=False,
                           bar_format="{desc}: {percentage:3.0f}%|{bar}|")
        
        try:
            # Process file with progress callback
            result = process_file_with_progress(file_path, args.type, file_progress)
            results.append(result)
            file_progress.close()
            if HAS_COLORAMA:
                total_progress.set_description(f"{Fore.BLUE}Overall Progress{Style.RESET_ALL} ({Fore.GREEN}{filename} ✓{Style.RESET_ALL})")
            else:
                total_progress.set_description(f"Overall Progress ({filename} ✓)")
            total_progress.update(1)
        except Exception as e:
            file_progress.close()
            if HAS_COLORAMA:
                total_progress.set_description(f"{Fore.BLUE}Overall Progress{Style.RESET_ALL} ({Fore.RED}{filename} ✗{Style.RESET_ALL})")
            else:
                total_progress.set_description(f"Overall Progress ({filename} ✗)")
            total_progress.update(1)
            print(f'Error processing {filename}: {e}', file=sys.stderr)

    total_progress.close()
    if HAS_COLORAMA:
        print(f"\n{Fore.GREEN}Successfully processed {len(results)} files{Style.RESET_ALL}")
    else:
        print(f"\nSuccessfully processed {len(results)} files")

    # Determine output filename
    base = os.path.basename(args.input.rstrip('/\\'))
    out_html = args.output or f"{base}_QC_report.html"
    make_report(results, out_html)
    print(f'Report written to {out_html}')


if __name__ == '__main__':
    main()

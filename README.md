# fasta_compare
A Python script for analysing DNA barcode sequences from multiple FASTA files, evaluating sequence quality based on multiple metrics, and selecting the best sequences according to configurable BOLD BIN quality criteria.

## Overview
This tool can process multiple FASTA files containing DNA sequences for barcode markers (cox1/COI, rbcl, matk), perform quality assessment, and output filtered sequences to separate FASTA files along with detailed analysis reports. It's specifically designed for barcode sequence analysis and BOLD (Barcode of Life Data) BIN compliance.

## Features
- **Quality assessment**: Evaluate sequences based on gaps, ambiguous bases, and continuous stretches (without gaps or ambiguous bases).
- **Flexible ranking systems**: 'Standard' and 'relaxed' barcode quality criteria.
- **Target-specific processing**: Handles and extracts different barcode sequences using known barcoding regions.
- **Comprehensive output and reporting**: Produces detailed CSV output with the following metrics:

| Column | Description |
|--------|-------------|
| file | Source FASTA file path |
| process_id | Extracted process identifier (e.g. BSNHM001-24) |
| parameters | Extracted BGEE run parameters (e.g. r_1_s_100) |
| seq_id | Full sequence identifier/found from fasta header (e.g. BSNHM001-24_r_1_s_100_BSNHM001-24) |
| length | Full sequence length (including gaps ('-') and ambiguous bases ('N') |
| leading_gaps | Count of leading gap ('-') characters |
| trailing_gaps | Count of trailing gap ('-') characters |
| internal_gaps | Count of internal gap ('-') characters |
| ambiguous_bases | Count of 'N' characters |
| longest_stretch | Longest continuous sequence without gaps or ambiguous bases |
| barcode_length | Length of barcode region (fixed by `--target`) |
| barcode_ambiguous_bases | Count of N in barcode region |
| barcode_longest_stretch | Longest continuous barcode sequence without gaps or ambiguous bases |
| barcode_rank | Barcode quality rank (1-6) |
| full_rank | Full sequence quality rank (1-3) |
| best_sequence | Whether this is the best sequence for the process_id |
| selected_full_fasta | Whether this sequence was output to full FASTA |
| selected_barcode_fasta | Whether this sequence was output to barcode FASTA |

## Installation
### Requirements/Dependencies
```bash
pip install biopython
# + Standard Python libraries (csv, re, os, argparse, logging, datetime)
```

## Usage
### Basic Usage
```bash
python fasta_compare.py \
    --output-csv results.csv \
    --output-fasta best_sequences.fasta \
    --output-barcode best_barcodes.fasta \
    --input file1.fasta file2.fasta \
    --target cox1
```
### Advanced Usage
```bash
# With custom rank threshold and relaxed criteria
python fasta_compare.py \
    --output-csv results.csv \
    --output-fasta best_sequences.fasta \
    --output-barcode best_barcodes.fasta \
    --input *.fasta \
    --target cox1 \
    --rank 2 \
    --relaxed \
    --verbose
```

## Command Line Arguments
### Required Arguments
| Argument | Short | Description |
|----------|-------|-------------|
| `--output-csv` | `-o` | Path to output CSV analysis file |
| `--output-fasta` | `-of` | Path to output FASTA file for best full sequences |
| `--output-barcode` | `-ob` | Path to output FASTA file for best barcode regions |
| `--input` | `-i` | Input FASTA files to analyze (space-separated) |
| `--target` | `-t` | Target genetic marker (cox1/COI/CO1, rbcl/RBCL, or matk/MATK) |

### Optional Arguments
| Argument | Description | Default |
|----------|-------------|---------|
| `--rank` | Maximum acceptable barcode rank for selection (see below for rank explanation) | 3 |
| `--relaxed` | Use relaxed barcode ranking criteria (see below for explanation) | False |
| `--log-file` | Custom path for log file | Auto-generated timestamp |
| `--verbose` `-v` | Enable detailed debug logging | False |

## Barcode Regions by Target
| Target | Barcode Region | Description |
|--------|----------------|-------------|
| cox1/COI | 40-700 | Cytochrome c oxidase subunit I |
| rbcl | 1-700 | RuBisCO large subunit |
| matk | 1-900 | Maturase K |
* ITS(2) region to be added in the near future

## Ranking Systems
### Standard Barcode Ranking (1-6, lower is better)
| Rank | Criteria |
|------|----------|
| 1 | No ambiguous bases, longest stretch ≥ 650 |
| 2 | No ambiguous bases, longest stretch ≥ 500 |
| 3 | No ambiguous bases, 300 ≤ longest stretch ≤ 499 |
| 4 | No ambiguous bases, 1 ≤ longest stretch ≤ 299 |
| 5 | Has ambiguous bases |
| 6 | Other cases |

### Relaxed Barcode Ranking (1-6, lower is better)
| Rank | Criteria | Notes |
|------|----------|-------|
| 1 | No ambiguous bases, longest stretch ≥ 500 | BIN compliant |
| 2 | No ambiguous bases, 200 ≤ longest stretch ≤ 499 | Good but no BIN |
| 3 | No ambiguous bases, 100 ≤ longest stretch ≤ 199 | Not great, better than nothing |
| 4 | <6 ambiguous bases, longest stretch ≥ 500 | |
| 5 | <6 ambiguous bases, longest stretch ≥ 300 | |
| 6 | Other cases | |

### Full Sequence Ranking (1-3, lower is better)
| Rank | Criteria |
|------|----------|
| 1 | No ambiguous bases |
| 2 | Has ambiguous bases |
| 3 | Other cases |

## Processing Logic
### 1. Sequence Analysis Loop
For each FASTA file:
```
For each sequence record:
├── Parse sequence ID and extract process_id/parameters
├── Calculate basic metrics:
│   ├── Length, leading/trailing gaps, internal gaps
│   ├── Ambiguous bases (N count)
│   └── Longest continuous stretch without gaps
├── Extract barcode region based on target marker
├── Calculate barcode-specific metrics:
│   ├── Barcode length, ambiguous bases
│   └── Longest barcode stretch without gaps
├── Determine rankings:
│   ├── Barcode rank (standard or relaxed)
│   └── Full sequence rank
└── Store all metrics and sequence record
```

### 2. Best Sequence Selection
For each unique process_id:
```
Selection criteria (in order of priority):
├── Barcode rank (lower is better)
├── Full sequence rank (lower is better)  
├── Barcode longest stretch (higher is better)
├── Full sequence longest stretch (higher is better)
├── Internal gaps (lower is better)
├── Ambiguous bases (lower is better)
└── Unique identifier (for deterministic selection)
Result = One "best" sequence per process_id
```

### 3. Output Sequence Selection
The script uses a two-tier selection process:
```
For each process_id:
├── First attempt: Find sequences where:
│   ├── full_rank == 1 (no ambiguous bases)
│   └── barcode_rank <= rank_threshold
├── If found: Use for BOTH full sequence and barcode output
├── Second attempt: Find sequences where:
│   ├── full_rank == 2 (has ambiguous bases)
│   └── barcode_rank <= rank_threshold
├── If found: Use ONLY for barcode output (no full sequence)
└── If none found: No output for this process_id
```

### 4. Sequence Formatting
#### Full Sequences
```
Original (full) sequence
├── Trim leading and trailing gaps ('-', '~')
├── Remove internal '~' characters (stitch sequence)
├── Replace internal '-' characters with 'N'
└── Trim leading and trailing 'N' characters
```
#### Barcode Sequences
```
Extract barcode region
├── Remove internal '~' characters (stitch sequence)
├── For gaps marked with '-':
│   ├── If ≤ 6 bases: fill with 'N'
│   └── If > 6 bases: keep longest fragment
└── Trim leading and trailing 'N' characters
```

### 5. Output Generation
```
Generate three output files:
├── CSV report: All sequences with complete metrics
├── Full sequences FASTA: Best full sequences (one per process_id)
└── Barcode FASTA: Best barcode regions (one per process_id)
```

## Valid FASTA header formats
The script expects sequence IDs to be in a particular format to correctly parse process_id and parameters, e.g.:
```
For BOLD123_r_1_s_100_BOLD123:
  - process_id = BOLD123
  - parameters = r_1_s_100_BOLD123

For SAMPLE_A_r_trial1:
  - process_id = SAMPLE_A
  - parameters = r_trial1

For complex_name_here_r_N_s_N_other_stuff:
  - process_id = complex_name_here
  - parameters = r_N_s_N_other_stuff

For Consensus_BOLD123_r_param1_param2:
  - process_id = `BOLD123`
  - parameters = `r_param1_param2`
```

## Gap Character Handling
| Character | Type | Processing |
|-----------|------|------------|
| `-` | Standard gap | Converted to N (if ≤6) or causes fragmentation (if >6) |
| `~` | Sequencing gap | Removed, sequences stitched together |
| `N` | Ambiguous base | Counted in quality metrics, trimmed from ends |


## Performance Considerations
- **Memory usage**: Sequences are held in memory during processing
- **Progress reporting**: Logs progress every 100 sequences for large files
- **File handling**: Processes files sequentially, not in parallel
- **Output formatting**: FASTA sequences written without line breaks for consistency

## Troubleshooting
### Common Issues
1. **"Input file does not exist"**: Check file paths and permissions
2. **"Invalid target marker"**: Use cox1/COI/CO1, rbcl/RBCL, or matk/MATK
3. **"Invalid rank value"**: Rank must be between 1-6
4. **Empty output files**: No sequences met the quality criteria; try `--relaxed` or a higher `--rank`

## Authors
Created by Ben Price & Dan Parsons @ Natural History Museum UK (NHMUK)

## Contributing
- Please suggest any improvements in 'Issues' - contributions welcome!

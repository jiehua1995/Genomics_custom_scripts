#!/usr/bin/env python3
"""
jiehua_custom - Main CLI entry point

This script provides a unified command-line interface to access all tools
in the jiehua_custom package. It can display help information about all
available tools and their functionality.

Usage:
    jiehua_custom --help       # Show this help and list all tools
    jiehua_custom --list       # List all available tools with descriptions
    jiehua_custom --version    # Show package version
"""

import argparse
import sys
from . import __version__, __author__, __email__


def get_tools_info():
    """
    Return information about all available tools in the package.
    
    Returns:
        List of tuples: (command_name, description, usage_example)
    """
    tools = [
        (
            "gff2bed",
            "Convert GFF files to BED format with coordinate adjustment",
            "gff2bed -i annotations.gff -o annotations.bed"
        ),
        (
            "file-checksum", 
            "Calculate file hashes for data integrity verification",
            "file-checksum -p /path/to/data -o checksums.txt -a sha256"
        ),
        (
            "fasta-stats",
            "Analyze FASTA files for sequence length and GC content statistics", 
            "fasta-stats genome.fasta -summarize -o stats.txt"
        )
    ]
    return tools

def print_main_help():
    """Print comprehensive help information for the jiehua_custom package."""
    print(f"""
jiehua_custom v{__version__}
A collection of bioinformatics command-line utilities

Author: {__author__} <{__email__}>

AVAILABLE TOOLS:
""")
    
    tools = get_tools_info()
    for cmd, desc, example in tools:
        print(f"  {cmd:<15} - {desc}")
        print(f"  {' ' * 15}   Example: {example}")
        print()
    
    print("""USAGE:
  # Run individual tools directly:
  gff2bed --help              # Get help for gff2bed tool
  file-checksum --help        # Get help for file-checksum tool  
  fasta-stats --help          # Get help for fasta-stats tool
  
  # Package information:
  jiehua_custom --help        # Show this help message
  jiehua_custom --list        # List tools only
  jiehua_custom --version     # Show version

INSTALLATION:
  # From GitHub (recommended)
  git clone https://github.com/jiehua1995/Genomics_custom_scripts.git
  cd Genomics_custom_scripts
  pip install -r requirements.txt
  pip install -e .
  
  # Alternative: direct pip install from GitHub  
  pip install git+https://github.com/jiehua1995/Genomics_custom_scripts.git

MORE INFO:
  Repository: https://github.com/jiehua1995/Genomics_custom_scripts
  Documentation: See README.md or individual tool help pages
""")

def print_tools_list():
    """Print a simple list of available tools."""
    print("Available tools in jiehua_custom:")
    tools = get_tools_info()
    for cmd, desc, _ in tools:
        print(f"  {cmd:<15} - {desc}")

def print_version():
    """Print version information."""
    print(f"jiehua_custom version {__version__}")
    print(f"Author: {__author__} <{__email__}>")

def main():
    """
    Main entry point for the jiehua_custom command.
    
    This function handles the main CLI interface and routes to appropriate
    help/information functions based on user input.
    """
    parser = argparse.ArgumentParser(
        description="jiehua_custom - Bioinformatics command-line utilities",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
For detailed help on individual tools, use:
  gff2bed --help
  file-checksum --help  
  fasta-stats --help

Examples:
  jiehua_custom --help       # Show comprehensive help
  jiehua_custom --list       # List available tools
  jiehua_custom --version    # Show version info
        """,
        add_help=False  # We'll handle --help manually
    )
    
    parser.add_argument("--help", "-h", action="store_true", 
                       help="Show help message and list all tools")
    parser.add_argument("--list", "-l", action="store_true",
                       help="List all available tools")  
    parser.add_argument("--version", "-v", action="store_true",
                       help="Show version information")
    
    # Parse args, but allow unknown args (for potential future extensions)
    args, unknown = parser.parse_known_args()
    
    # Handle different command options
    if args.version:
        print_version()
    elif args.list:
        print_tools_list()
    elif args.help or len(sys.argv) == 1:
        # Show help by default if no args, or if --help specified
        print_main_help()
    else:
        # If unknown arguments, suggest correct usage
        if unknown:
            print(f"Unknown arguments: {' '.join(unknown)}")
            print("Use 'jiehua_custom --help' for usage information.")
            sys.exit(1)
        else:
            print_main_help()

if __name__ == "__main__":
    main()
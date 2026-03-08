#!/usr/bin/env python3

import argparse
from pathlib import Path

from GeneTicker.core import run_codon_analysis


def main():
    """
    Main function for command-line execution: parses arguments and calls the
    core analysis function.
    """
    parser = argparse.ArgumentParser(
        description="Calculate codon frequencies from CDS features in a GBFF file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        type=Path,  # Use pathlib.Path directly for argument type
        help="The path to the input GBFF file (.gbff, .gbk, etc.)",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Optional: The path and filename for the output file. Format is inferred from the extension (.csv, .xlsx, .json, etc.). If not provided, results are printed to the console.",
    )

    # Flag to include non-canonical codon-AA pairs
    parser.add_argument(
        "--include-non-canonical",
        action="store_true",
        help="Include non-canonical codon-amino acid pairs in the output. By default, only canonical pairs (according to the detected translation table) are included.",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Indicate whether to print detailed processing information to the console.",
    )

    args = parser.parse_args()

    # Call the main analysis function with parsed arguments
    run_codon_analysis(
        input_file=args.input_file,
        output=args.output,
        include_non_canonical=args.include_non_canonical,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()

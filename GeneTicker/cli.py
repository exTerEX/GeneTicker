#!/usr/bin/env python3

import argparse
from pathlib import Path

from GeneTicker.core import DEFAULT_MIN_FREQ, run_codon_analysis


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

    # New flag to disable the frequency filter
    parser.add_argument(
        "--no-freq-filter",
        action="store_true",
        help=f"Do not filter out low-frequency codons. By default, codons with < {DEFAULT_MIN_FREQ:.2f}%% frequency are removed.",
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
        no_freq_filter=args.no_freq_filter,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()

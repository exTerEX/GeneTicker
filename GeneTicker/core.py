"""
Calculates the codon frequency from all CDS (Coding DNA Sequence) features
in a given GenBank Flat File (GBFF). Results are either written to a file
(format determined by the --output extension), or printed to the console.

This script uses Biopython to parse the file, iterates through all features,
extracts the DNA sequence for each CDS, and counts the occurrences of
each 64-possible codon, mapping them to their corresponding amino acid.
"""

from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from GeneTicker.utils import export_results_to_file, print_results

# Define a type alias for our codon counter for clarity
CodonAACounter = Counter[Tuple[str, str, str]]

# Define the set of valid, unambiguous DNA bases
VALID_BASES = set("ATCG")

# Default minimum frequency to include a codon in the final report
DEFAULT_MIN_FREQ = 0.05


def get_codon_table(table_id: int = 1) -> Dict[str, str]:
    """
    Fetch a codon table from BioPython's NCBI translation tables.

    :param table_id: NCBI translation table ID (default 1 = Standard)
    :type table_id: int
    :return: Dictionary mapping codons to amino acids
    :rtype: Dict[str, str]
    """
    try:
        table = CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError:
        print(f"Warning: Unknown translation table {table_id}, falling back to standard (1)")
        table = CodonTable.unambiguous_dna_by_id[1]

    codon_dict: Dict[str, str] = {}

    # Add regular codon -> AA mappings
    for codon, aa in table.forward_table.items():
        codon_dict[codon] = aa

    # Add stop codons
    for codon in table.stop_codons:
        codon_dict[codon] = "*"

    return codon_dict


def calculate_codon_frequency_to_df(
    gbff_file: Path, include_non_canonical: bool = False, verbose: bool = False
) -> Optional[pd.DataFrame]:
    """
    Parses a GenBank file, calculates codon frequencies from all CDS features,
    and returns the results as a Pandas DataFrame.

    The first codon of each CDS is tracked as 'START', and stop codons as 'STOP'.
    The translation table is automatically detected from /transl_table qualifiers.

    :param gbff_file: The file path to the input .gbff file
    :type gbff_file: Path
    :param include_non_canonical: Include non-canonical codon-AA pairs, defaults to False
    :type include_non_canonical: bool, optional
    :param verbose: If True, prints detailed progress information, defaults to False
    :type verbose: bool, optional
    :return: A Pandas DataFrame with codon frequency data, or None if no data is produced.
    :rtype: Optional[pd.DataFrame]"""
    if verbose:
        print(f"--- Starting analysis of {gbff_file.name} ---")

    codon_aa_counts: CodonAACounter = Counter()
    total_cds_processed = 0
    transl_table_counts: Counter[int] = Counter()

    try:
        for record in SeqIO.parse(gbff_file, "genbank"):
            record: SeqRecord
            if verbose:
                print(f"Processing record: {record.id}")

            for feature in record.features:
                feature: SeqFeature

                if feature.type == "CDS":
                    locus_tag = feature.qualifiers.get("locus_tag", ["?"])[0]
                    translation_str = feature.qualifiers.get("translation", [""])[0]

                    # Track translation table usage
                    transl_table_str = feature.qualifiers.get("transl_table", ["1"])[0]
                    try:
                        transl_table_id = int(transl_table_str)
                        transl_table_counts[transl_table_id] += 1
                    except ValueError:
                        transl_table_counts[1] += 1  # Default to standard

                    if not translation_str:
                        if verbose:
                            print(f"  > Warning: CDS {locus_tag} is missing /translation qualifier. Skipping.")
                        continue

                    if "pseudo" in feature.qualifiers:
                        if verbose:
                            print(f"  > Skipping pseudogene: {locus_tag}")
                        continue

                    total_cds_processed += 1

                    cds_seq: Seq = feature.extract(record.seq)

                    cds_seq_str = str(cds_seq).upper()
                    if not set(cds_seq_str).issubset(VALID_BASES):
                        if verbose:
                            print(
                                f"  > Warning: CDS {locus_tag} (length {len(cds_seq)}) "
                                "contains non-ATCG bases. Skipping."
                            )
                        continue

                    cds_len = len(cds_seq)
                    truncated_len = (cds_len // 3) * 3

                    if cds_len != truncated_len and verbose:
                        print(
                            f"  > Warning: CDS {locus_tag} "
                            f"(length {cds_len}) is not a multiple of 3. "
                            f"Truncating to {truncated_len} bases for counting."
                        )

                    # Iterate over the sequence in 3-base-pair steps (codons)
                    for i in range(0, truncated_len, 3):
                        codon = cds_seq_str[i : i + 3]
                        aa_index = i // 3

                        aa_code: str
                        special_type: str = ""

                        # Use the translation string for the AA code for all non-stop codons.
                        if aa_index < len(translation_str):
                            aa_code = translation_str[aa_index]

                            if i == 0:
                                special_type = "START"
                        else:
                            # Handle STOP codon
                            aa_code = "*"
                            special_type = "STOP"

                        # Use a tuple as the key to avoid string parsing later
                        codon_aa_counts[(codon, aa_code, special_type)] += 1

    except FileNotFoundError:
        print(f"Error: File not found at {gbff_file}")
        return None
    except Exception as e:
        print(f"An error occurred during parsing: {e}")
        return None

    if verbose:
        print(f"--- Analysis complete. Processed {total_cds_processed} CDS features. ---")

    total_codons = sum(codon_aa_counts.values())
    if total_codons == 0:
        return None

    # Determine the most common translation table used in this file
    if transl_table_counts:
        detected_table_id = transl_table_counts.most_common(1)[0][0]
    else:
        detected_table_id = 1  # Default to standard

    if verbose:
        print(f"Using NCBI translation table {detected_table_id}")

    codon_table = get_codon_table(detected_table_id)

    data: List[dict] = []
    # Unpack the tuple key
    for (codon, aa_code, special_type), count in sorted(codon_aa_counts.items()):
        frequency_percent = (count / total_codons) * 100

        # Determine if this codon-AA pair is canonical
        canonical_aa = codon_table.get(codon, "")
        is_canonical = (aa_code == canonical_aa) or (special_type == "STOP" and canonical_aa == "*")

        # Build special type string, adding Sec/Pyl markers
        special_types = []
        if special_type:  # START or STOP
            special_types.append(special_type)
        if aa_code == "U":
            special_types.append("Sec")  # Selenocysteine
        elif aa_code == "O":
            special_types.append("Pyl")  # Pyrrolysine
        special_type_str = ",".join(special_types)

        data.append(
            {
                "Amino Acid": "*" if special_type == "STOP" or aa_code == "*" else aa_code,
                "Codon": codon,
                "Count": count,
                "Freq. (%)": frequency_percent,
                "Special Type": special_type_str,
                "Canonical": is_canonical,
            }
        )
    df = pd.DataFrame(data)

    if df.empty:
        return None

    df = df[["Amino Acid", "Codon", "Count", "Freq. (%)", "Special Type", "Canonical"]]

    if not include_non_canonical:
        if verbose:
            print("Filtering out non-canonical codon-amino acid pairs...")
        df = df[df["Canonical"] == True].reset_index(drop=True)
        if verbose:
            print(f"Remaining codons after filtering: {len(df)}")

    return df


def run_codon_analysis(
    input_file: Path, output: Optional[Path] = None, include_non_canonical: bool = False, verbose: bool = False
) -> Optional[pd.DataFrame]:
    """
    Runs the full codon analysis pipeline.

    This function is callable by other Python scripts.

    :param input_file: The path to the input GBFF file
    :type input_file: Path
    :param output: Path for output file, defaults to None
    :type output: Optional[Path], optional
    :param include_non_canonical: Include non-canonical codon-AA pairs, defaults to False
    :type include_non_canonical: bool, optional
    :param verbose: If True, prints detailed progress information, defaults to False
    :type verbose: bool, optional
    :return: The resulting Pandas DataFrame if analysis is successful, otherwise None.
    :rtype: Optional[pd.DataFrame]
    """
    if not input_file.exists() or not input_file.is_file():
        print(f"Error: Input file not found or is not a file: {input_file}")
        return None

    df = calculate_codon_frequency_to_df(input_file, include_non_canonical, verbose=verbose)

    if df is None:
        if verbose:
            print("No valid CDS features were processed or no data remained after filtering.")
        return None

    if output:
        # If the output path is provided, try to export to file
        export_success = export_results_to_file(df, output)

        if not export_success:
            # If export failed due to unsupported extension or error,
            # fall back to printing the result as requested.
            print("\n--- Falling back to console output due to export failure ---")
            print_results(df, input_file.name)
    else:
        # If no output path is provided, print to console
        print_results(df, input_file.name)

    return df

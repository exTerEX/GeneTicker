from pathlib import Path

import pandas as pd


def export_results_to_file(df: pd.DataFrame, output_path: Path, verbose: bool = False) -> bool:
    """
    Writes a Pandas DataFrame to a file, inferring the format from the extension.
    Returns True on success, False on failure (e.g., unsupported extension).

    :param df: Dataframe to export
    :type df: pd.DataFrame
    :param output_path: Ouput file path. Provide extension to infer format.
    :type output_path: Path
    :param verbose: Show additional information for debugging purposes, defaults to False
    :type verbose: bool, optional
    :return: Indicate success or failure of the export operation
    :rtype: bool
    """
    if verbose:
        print(f"\n--- Creating Report at {output_path} ({output_path.suffix.lower()}) ---")

    ext = output_path.suffix.lower()

    try:
        if ext == ".xlsx":
            df.to_excel(output_path, index=False, sheet_name="Codon Frequencies")
        elif ext == ".csv":
            df.to_csv(output_path, index=False)
        elif ext == ".tsv":
            df.to_csv(output_path, index=False, sep="\t")
        elif ext == ".json":
            df.to_json(output_path, orient="records", indent=4)
        elif ext == ".pickle":
            df.to_pickle(output_path)
        elif ext == ".xml":
            # Requires lxml and possibly openpyxl for XML functionality
            df.to_xml(output_path, index=False)
        elif ext in (".feather", ".parquet", ".orc"):
            # Requires pyarrow or fastparquet engine
            if ext == ".feather":
                df.to_feather(output_path)
            elif ext == ".parquet":
                df.to_parquet(output_path)
            elif ext == ".orc":
                df.to_orc(output_path)
        else:
            # Handle unsupported extension: print error and return False
            print(
                f"Error: Unsupported output file extension: {ext}. "
                "Please use a supported extension (.csv, .xlsx, .tsv, .json, .pickle, .xml, .feather, .parquet, or .orc)."
            )
            return False  # Indicate failure

        if verbose:
            print(f"Successfully exported {len(df)} codon entries to {output_path}")
        return True  # Indicate success

    except ImportError as e:
        print(
            f"Dependency Error for {ext} export: {e}. Please ensure the necessary libraries (e.g., 'pyarrow' for feather/parquet/orc, 'lxml' for xml) are installed."
        )
        return False
    except Exception as e:
        print(f"An unexpected error occurred during file export to {output_path}: {e}")
        return False


def print_results(df: pd.DataFrame, input_file_name: str):
    """
    Prints the codon-to-amino-acid frequency results in a formatted table to the console.
    """
    total_codons = df["Count"].sum()

    print(f"\n=== Codon Frequency Report for {input_file_name} ===")
    print(f"Total Valid Codons Counted: {total_codons}\n")

    # Custom formatting for console output - adjusted width for new column
    print(f"{'AA':<4} {'Codon':<6} {'Count':<10} {'Freq. (%)':>10} {'Special Type':<12}")
    print("-" * 49)

    # Iterate over the DataFrame rows to print
    for index, row in df.iterrows():
        print(
            f"{row['Amino Acid']:<4} "
            f"{row['Codon']:<6} "
            f"{row['Count']:<10} "
            f"{row['Freq. (%)']:>10.4f} "
            f"{row['Special Type']:<12}"
        )
    print("=" * 49)

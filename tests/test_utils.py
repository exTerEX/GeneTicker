import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pandas.testing as pdt

from GeneTicker.utils import export_results_to_file, print_results


def _sample_df():
    return pd.DataFrame(
        [
            {"Amino Acid": "Ala", "Codon": "GCT", "Count": 10, "Freq. (%)": 50.0, "Special Type": "none"},
            {"Amino Acid": "Gly", "Codon": "GGT", "Count": 5, "Freq. (%)": 25.0, "Special Type": "none"},
            {"Amino Acid": "Ser", "Codon": "TCT", "Count": 5, "Freq. (%)": 25.0, "Special Type": "start"},
        ]
    )


def test_export_csv_success(tmp_path: Path):
    df = _sample_df()
    out = tmp_path / "out.csv"
    ok = export_results_to_file(df, out)
    assert ok is True
    assert out.exists()
    # Read back and compare
    read = pd.read_csv(out)
    # Ensure Count column types comparable
    read["Count"] = read["Count"].astype(int)
    pdt.assert_frame_equal(df.reset_index(drop=True), read.reset_index(drop=True))


def test_export_tsv_success(tmp_path: Path):
    df = _sample_df()
    out = tmp_path / "out.tsv"
    ok = export_results_to_file(df, out)
    assert ok is True
    assert out.exists()
    read = pd.read_csv(out, sep="\t")
    read["Count"] = read["Count"].astype(int)
    pdt.assert_frame_equal(df.reset_index(drop=True), read.reset_index(drop=True))


def test_export_json_and_pickle_success(tmp_path: Path):
    df = _sample_df()
    out_json = tmp_path / "out.json"
    out_pickle = tmp_path / "out.pickle"

    assert export_results_to_file(df, out_json) is True
    assert out_json.exists()
    read_json = pd.read_json(out_json, orient="records")
    # Ensure Freq. (%) preserves float dtype when read back from JSON
    read_json["Freq. (%)"] = read_json["Freq. (%)"].astype(float)
    # JSON may cast ints as int64; ensure equality ignoring dtypes for simplicity
    pdt.assert_frame_equal(df.reset_index(drop=True), read_json.reset_index(drop=True))

    assert export_results_to_file(df, out_pickle) is True
    assert out_pickle.exists()
    read_pickle = pd.read_pickle(out_pickle)
    pdt.assert_frame_equal(df.reset_index(drop=True), read_pickle.reset_index(drop=True))


def test_export_unsupported_extension(tmp_path: Path, capsys):
    df = _sample_df()
    out = tmp_path / "out.unsupported"
    ok = export_results_to_file(df, out)
    captured = capsys.readouterr()
    assert ok is False
    assert "Unsupported output file extension" in captured.out
    assert not out.exists()


def test_export_parquet_raises_importerror(monkeypatch, tmp_path: Path):
    """
    Simulate a missing dependency for parquet by forcing DataFrame.to_parquet to raise ImportError.
    The function should catch ImportError and return False.
    """
    df = _sample_df()
    out = tmp_path / "out.parquet"

    def _raise_importerror(self, path, *args, **kwargs):
        raise ImportError("pyarrow is not installed")

    monkeypatch.setattr(pd.DataFrame, "to_parquet", _raise_importerror)
    ok = export_results_to_file(df, out)
    assert ok is False
    assert not out.exists()


def test_print_results_formatting(capsys):
    df = _sample_df()
    # Intentionally change order to ensure printing iterates rows
    df = df.loc[[0, 1, 2]].reset_index(drop=True)
    print_results(df, "sample_input.fasta")
    captured = capsys.readouterr()
    out = captured.out

    # Check header and totals
    assert "Codon Frequency Report for sample_input.fasta" in out
    assert "Total Valid Codons Counted: 20" in out

    # Check table header and separators
    assert "AA" in out and "Codon" in out and "Count" in out and "Freq. (%)" in out
    assert "-" * 49 in out
    # Check that each codon line appears with its codon and count
    assert "GCT" in out and "10" in out
    assert "GGT" in out and "5" in out
    assert "TCT" in out and "5" in out

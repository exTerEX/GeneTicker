import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pathlib import Path as P
import textwrap

from GeneTicker.core import calculate_codon_frequency_to_df, run_codon_analysis


def _write_simple_gbff(path: P):
    """Write a minimal GenBank (GBFF) file with a single CDS feature.

    The sequence contains three codons: ATG GCT GGT which translate to M G G
    (translation string "MGG"). This ensures predictable counts of 1 per codon.
    """
    content = textwrap.dedent(
        """
        LOCUS       TEST01                9 bp    DNA     linear   01-JAN-1980
        DEFINITION  .
        ACCESSION   TEST
        VERSION     TEST.1
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
        FEATURES             Location/Qualifiers
             CDS             1..9
                             /locus_tag="TEST1"
                             /translation="MGG"
        ORIGIN
                1 atggctggt
        //
        """
    )
    path.write_text(content)


def test_calculate_codon_frequency_basic(tmp_path: Path):
    gbff = tmp_path / "sample.gbff"
    _write_simple_gbff(gbff)

    # Use 0.0 threshold to include all codons
    df = calculate_codon_frequency_to_df(gbff, 0.0, verbose=False)
    assert df is not None

    # Expect three codons each with count 1
    assert len(df) == 3
    counts = dict(zip(df['Codon'].tolist(), df['Count'].tolist()))
    assert counts.get('ATG') == 1
    assert counts.get('GCT') == 1
    assert counts.get('GGT') == 1


def test_calculate_codon_frequency_missing_file(tmp_path: Path):
    missing = tmp_path / "no_such.gbff"
    df = calculate_codon_frequency_to_df(missing, 0.0, verbose=False)
    assert df is None


def test_run_codon_analysis_fallback_on_export_failure(monkeypatch, tmp_path: Path, capsys):
    gbff = tmp_path / "sample2.gbff"
    _write_simple_gbff(gbff)

    # Import the module to monkeypatch the bound names inside it
    import GeneTicker.core as core

    # Force export to fail so run_codon_analysis prints fallback message
    monkeypatch.setattr(core, "export_results_to_file", lambda df, out: False)
    # Prevent print_results from producing noisy output during test
    monkeypatch.setattr(core, "print_results", lambda df, name: None)

    outpath = tmp_path / "out.csv"
    df = run_codon_analysis(gbff, output=outpath, no_freq_filter=True, verbose=True)
    assert df is not None

    captured = capsys.readouterr()
    assert "Falling back to console output" in captured.out

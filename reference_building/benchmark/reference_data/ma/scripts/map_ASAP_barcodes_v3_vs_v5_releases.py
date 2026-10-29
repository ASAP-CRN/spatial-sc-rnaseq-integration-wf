#!/usr/bin/env python3
"""
Map ASAP cohort cell barcodes between the v3.1.0 and v5.0.0 release structures.

The two releases encode barcodes differently:
  v3.1.0: TTGAACGGTGCCTACG-1_ASAP_PMBDS_000038_s001_Rep1.cleaned_unfiltered.h5ad
  v5.0.0: CTAGTGACAATCTACG-1_DS_PMDBS_0001.ASAP_PMDBS_000002_s002_BiolRep1_0

Both encode the same barcode sequence ("TTGAACGGTGCCTACG-1") and ASAP
sample ID, so matching on (sample + barcode sequence) -- after normalizing
v3.1.0's "PM{B<>D}S" typo -- maps the two releases 1:1 for most cells; 
the replicate/suffix portion is dropped, since it isn't directly comparable.
"""

import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))

from utils.logging_extra import log_run_command

_V3_SEED = "reference_building/reference_h5ad/deconvolution_reference_300k_c2lUda_Kamath200Rand.h5ad"
_V5_FINAL_CSV = (
    "reference_building/da_neurons_identification_analysis/reference_data/asap_cohort/"
    "pmdbs_sc_rnaseq/cohort_analysis/5.0.0/2026-10-19T14-59-48Z/asap-cohort.final_metadata.csv.gz"
)
_V5_UNFILTERED_H5AD = (
    "reference_building/da_neurons_identification_analysis/reference_data/asap_cohort/"
    "pmdbs_sc_rnaseq/cohort_analysis/5.0.0/2026-10-19T14-59-48Z/"
    "asap-cohort.merged_cleaned_unfiltered.h5ad"
)
_V3_FINAL_CSV = (
    "reference_building/da_neurons_identification_analysis/reference_data/asap_cohort/"
    "pmdbs_sc_rnaseq/cohort_analysis/3.1.0/2025-12-01T22-36-30Z/asap-cohort.final_metadata.csv.gz"
)
_OUTDIR = "reference_building/benchmark/asap_barcode_mapping_v3_vs_v5"
_STEM = "asap_barcode_mapping_v3_vs_v5"

_V3_BARCODE_RE = re.compile(
    r"^(?P<seq>[ACGT]+-\d+)_(?P<sample>ASAP_PM(?:BD|DB)S_\d+_s\d+)_(?P<rep>.+?)\.cleaned_unfiltered\.h5ad$"
)
_V5_BARCODE_RE = re.compile(
    r"^(?P<seq>[ACGT]+-\d+)_(?P<dataset>DS_PMDBS_\d+)\.(?P<sample>ASAP_PM(?:BD|DB)S_\d+_s\d+)_(?P<rep>.+)$"
)


def _normalize_sample(sample: str) -> str:
    """Correct the ASAP_PMBDS_ typo to ASAP_PMDBS_, matching sc-rnaseq-wf's own convention."""
    return sample.replace("PMBDS", "PMDBS")


def load_barcode_index(path: str) -> pd.Index:
    """
    Load raw barcode strings from an H5AD or a CSV(.gz) file, by extension.

    Parameters
    ----------
    path : str
        Path to an `.h5ad` file (barcodes read from `obs_names`) or a
        `.csv`/`.csv.gz` file (barcodes read from its first column).

    Returns
    -------
    pandas.Index
        Raw barcode strings, not yet parsed into (seq, sample).
    """
    if path.endswith(".h5ad"):
        adata = ad.read_h5ad(path, backed="r")
        return adata.obs_names
    header = pd.read_csv(path, nrows=0)
    idx_col = header.columns[0]
    return pd.Index(pd.read_csv(path, usecols=[idx_col], dtype=str)[idx_col])


def parse_v3_barcodes(barcodes) -> pd.DataFrame:
    """
    Parse v3.1.0-style barcodes into (barcode, seq, sample) rows.

    Parameters
    ----------
    barcodes : iterable of str
        Raw v3.1.0 barcodes (from `load_barcode_index`), e.g.
        "TTGAACGGTGCCTACG-1_ASAP_PMBDS_000038_s001_Rep1.cleaned_unfiltered.h5ad".
        A barcode not matching the expected pattern is dropped, with a warning.

    Returns
    -------
    pandas.DataFrame
        barcode, seq, sample (sample is typo-normalized to ASAP_PMDBS_...).
    """
    rows = []
    n_unmatched = 0
    for bc in barcodes:
        m = _V3_BARCODE_RE.match(bc)
        if not m:
            n_unmatched += 1
            continue
        rows.append({"barcode": bc, "seq": m["seq"], "sample": _normalize_sample(m["sample"])})
    if n_unmatched:
        print(f"[WARNING] {n_unmatched:,} v3.1.0 barcode(s) did not match the expected pattern; skipped.",
              file=sys.stderr)
    return pd.DataFrame(rows, columns=["barcode", "seq", "sample"])


def parse_v5_barcodes(barcodes) -> pd.DataFrame:
    """
    Parse v5.0.0-style barcodes into (barcode, seq, sample) rows.

    Parameters
    ----------
    barcodes : iterable of str
        Raw v5.0.0 obs_names/index, e.g.
        "CTAGTGACAATCTACG-1_DS_PMDBS_0001.ASAP_PMDBS_000002_s002_BiolRep1_0".
        A barcode not matching the expected pattern (Kamath/Siletti spike-ins,
        which use a completely different convention) is silently excluded --
        out of scope for this mapping, not a parsing failure.

    Returns
    -------
    pandas.DataFrame
        barcode, seq, sample.
    """
    rows = []
    for bc in barcodes:
        m = _V5_BARCODE_RE.match(bc)
        if not m:
            continue
        rows.append({"barcode": bc, "seq": m["seq"], "sample": _normalize_sample(m["sample"])})
    return pd.DataFrame(rows, columns=["barcode", "seq", "sample"])


def build_mapping(v3_df: pd.DataFrame, v5_df: pd.DataFrame) -> pd.DataFrame:
    """
    Join v3.1.0 and v5.0.0 barcodes on (sample, seq), one row per v3.1.0 barcode.

    Status accounts for ambiguity on *either* side: if more than one v3.1.0
    barcode or more than one v5.0.0 barcode shares a (sample, seq) key (e.g.
    the same physical barcode sequenced across two replicate runs), every
    v3.1.0 barcode in that group is marked "ambiguous_in_v5.0.0_final" rather than guessing
    a pairing.

    Parameters
    ----------
    v3_df : pandas.DataFrame
        From `parse_v3_barcodes`.
    v5_df : pandas.DataFrame
        From `parse_v5_barcodes`.

    Returns
    -------
    pandas.DataFrame
        v3.1.0_barcode, v5.0.0_barcode, status:
          "1:1_vs_v5.0.0_final"       : exactly one v5.0.0 barcode shares this
                                        (sample, seq) key.
          "ambiguous_in_v5.0.0_final" : more than one v3.1.0 and/or v5.0.0
                                        barcode shares it; v5.0.0_barcode
                                        lists every candidate, ";"-separated.
          "not_mappable"              : no v5.0.0 barcode shares it (may later
                                        be reclassified as
                                        "only_in_cleaned_unfiltered_in_v5.0.0"
                                        by `main`).
    """
    v5_grouped = v5_df.groupby(["sample", "seq"])["barcode"].apply(list)
    v3_grouped = v3_df.groupby(["sample", "seq"])["barcode"].apply(list)

    rows = []
    for key, v3_barcodes in v3_grouped.items():
        v5_barcodes = v5_grouped.get(key, [])
        if not v5_barcodes:
            for v3_bc in v3_barcodes:
                rows.append({"v3.1.0_barcode": v3_bc, "v5.0.0_barcode": "", "status": "not_mappable"})
        elif len(v3_barcodes) == 1 and len(v5_barcodes) == 1:
            rows.append({
                "v3.1.0_barcode": v3_barcodes[0],
                "v5.0.0_barcode": v5_barcodes[0],
                "status": "1:1_vs_v5.0.0_final",
            })
        else:
            for v3_bc in v3_barcodes:
                rows.append({
                    "v3.1.0_barcode": v3_bc,
                    "v5.0.0_barcode": ";".join(v5_barcodes),
                    "status": "ambiguous_in_v5.0.0_final",
                })
    return pd.DataFrame(rows, columns=["v3.1.0_barcode", "v5.0.0_barcode", "status"])


def reclassify_filtered_cells(
    mapping: pd.DataFrame, v3_df: pd.DataFrame, v5_unfiltered_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Relabel "not_mappable" rows found in v5.0.0's pre-QC-filter H5AD.

    A v3.1.0 barcode whose (sample, seq) key is present in v5.0.0's
    merged_cleaned_unfiltered H5AD but absent from the final (QC-filtered)
    metadata was dropped by v5.0.0's QC filtering, not lost at merge -- these
    are relabelled "only_in_cleaned_unfiltered_in_v5.0.0" in place, with
    "v5.0.0_barcode" populated from the pre-filter H5AD (";"-joined if more
    than one candidate shares the key).

    Parameters
    ----------
    mapping : pandas.DataFrame
        From `build_mapping`.
    v3_df : pandas.DataFrame
        From `parse_v3_barcodes` -- the same frame passed to `build_mapping`.
    v5_unfiltered_df : pandas.DataFrame
        From `parse_v5_barcodes`, applied to the pre-filter merged H5AD.

    Returns
    -------
    pandas.DataFrame
        `mapping`, with eligible "not_mappable" rows relabelled and their
        "v5.0.0_barcode" filled in.
    """
    unfiltered_grouped = v5_unfiltered_df.groupby(["sample", "seq"])["barcode"].apply(list)
    key_by_v3_barcode = dict(zip(v3_df["barcode"], zip(v3_df["sample"], v3_df["seq"])))

    is_not_mappable = mapping["status"] == "not_mappable"
    keys = mapping.loc[is_not_mappable, "v3.1.0_barcode"].map(key_by_v3_barcode)
    matched_barcodes = keys.map(lambda key: unfiltered_grouped.get(key, []))

    found_in_unfiltered = matched_barcodes.map(len) > 0
    reclassify_idx = found_in_unfiltered[found_in_unfiltered].index
    mapping.loc[reclassify_idx, "status"] = "only_in_cleaned_unfiltered_in_v5.0.0"
    mapping.loc[reclassify_idx, "v5.0.0_barcode"] = matched_barcodes.loc[reclassify_idx].map(";".join)
    return mapping


def reclassify_missing_samples(
    mapping: pd.DataFrame, v3_df: pd.DataFrame, v5_unfiltered_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Relabel "not_mappable" rows whose ASAP sample is absent from v5.0.0 entirely.

    Distinguishes a whole sample never carried into v5.0.0 (this label) from a
    single cell within an otherwise-present sample that simply has no v5.0.0
    counterpart (remains "not_mappable").

    Parameters
    ----------
    mapping : pandas.DataFrame
        From `build_mapping`, already passed through `reclassify_filtered_cells`.
    v3_df : pandas.DataFrame
        From `parse_v3_barcodes` -- the same frame passed to `build_mapping`.
    v5_unfiltered_df : pandas.DataFrame
        From `parse_v5_barcodes`, applied to the pre-filter merged H5AD.

    Returns
    -------
    pandas.DataFrame
        `mapping`, with eligible "not_mappable" rows relabelled
        "sample_id_not_found_in_v5.0.0".
    """
    v5_samples = set(v5_unfiltered_df["sample"])
    sample_by_v3_barcode = dict(zip(v3_df["barcode"], v3_df["sample"]))

    is_not_mappable = mapping["status"] == "not_mappable"
    sample_missing = mapping.loc[is_not_mappable, "v3.1.0_barcode"].map(
        lambda bc: sample_by_v3_barcode[bc] not in v5_samples
    )
    reclassify_idx = sample_missing[sample_missing].index
    mapping.loc[reclassify_idx, "status"] = "sample_id_not_found_in_v5.0.0"
    return mapping


def add_in_v3_final_column(mapping: pd.DataFrame, v3_final_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add an "in_v3.1.0_final" boolean column to `mapping`.

    True when the v3.1.0 barcode also appears in v3.1.0's own final
    (QC-filtered) metadata -- i.e. it would have been included whether
    Strategy 1's reference was built from v3.1.0's final metadata or its
    pre-QC merged H5AD; False means it only ever existed pre-QC in v3.1.0.

    Parameters
    ----------
    mapping : pandas.DataFrame
        From `build_mapping`, keyed on the "v3.1.0_barcode" column.
    v3_final_df : pandas.DataFrame
        From `parse_v3_barcodes`, applied to v3.1.0's own final metadata CSV.

    Returns
    -------
    pandas.DataFrame
        `mapping` with an added "in_v3.1.0_final" boolean column.
    """
    v3_final_barcodes = set(v3_final_df["barcode"])
    mapping["in_v3.1.0_final"] = mapping["v3.1.0_barcode"].isin(v3_final_barcodes)
    return mapping


def build_summary_table(mapping: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize `mapping` by (in_v3.1.0_final, status) counts and percentages.

    Parameters
    ----------
    mapping : pandas.DataFrame
        From `build_mapping`, after `add_in_v3_final_column`.

    Returns
    -------
    pandas.DataFrame
        in_v3.1.0_final, status, count, "% of v3.1.0" (rounded to 1 decimal),
        sorted by in_v3.1.0_final (descending) then count (descending).
    """
    total = len(mapping)
    summary = (
        mapping.groupby(["in_v3.1.0_final", "status"])
        .size()
        .reset_index(name="count")
    )
    summary["% of v3.1.0"] = (summary["count"] / total * 100).round(1)
    return summary.sort_values(
        ["in_v3.1.0_final", "count"], ascending=[False, False]
    ).reset_index(drop=True)


def summary_table_to_markdown(summary: pd.DataFrame) -> str:
    """
    Render `summary` as a GitHub-flavoured Markdown table.

    Parameters
    ----------
    summary : pandas.DataFrame
        From `build_summary_table`.

    Returns
    -------
    str
        Markdown table text, newline-terminated.
    """
    columns = list(summary.columns)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in summary.itertuples(index=False)
    ]
    return "\n".join([header, separator, *rows]) + "\n"


def main(
    v3_seed: str,
    v5_final_csv: str,
    v5_unfiltered_h5ad: str,
    v3_final_csv: str,
    outdir: str,
) -> None:
    """
    Build and write the v3.1.0-to-v5.0.0 barcode mapping TSV and MD summary.

    Parameters
    ----------
    v3_seed : str
        Path to the v3.1.0 barcode set to map (`.h5ad` or `.csv`/`.csv.gz`):
        either the Strategy 1 deconvolution reference H5AD (ASAP barcodes
        plus Siletti/Kamath cells, which are skipped as unmatched) or
        v3.1.0's own full final metadata CSV (same file as `v3_final_csv`).
    v5_final_csv : str
        Path to the v5.0.0 final (QC-filtered) metadata CSV (only the index
        is read).
    v5_unfiltered_h5ad : str
        Path to the v5.0.0 pre-QC-filter merged_cleaned_unfiltered H5AD, used
        to tell "dropped by QC filtering" apart from "never merged".
    v3_final_csv : str
        Path to v3.1.0's own final (QC-filtered) metadata CSV (only the index
        is read), used to populate the "in_v3.1.0_final" column.
    outdir : str
        Destination directory. Outputs are written as `{_STEM}.tsv.gz` and
        `{_STEM}_summary.md`.

    Returns
    -------
    None
    """
    print(f"[INFO] Loading v3.1.0 barcodes from {v3_seed}")
    v3_df = parse_v3_barcodes(load_barcode_index(v3_seed))
    print(f"  {len(v3_df):,} v3.1.0 barcode(s) parsed")

    print(f"[INFO] Loading v5.0.0 barcodes from {v5_final_csv}")
    v5_index = load_barcode_index(v5_final_csv)
    v5_df = parse_v5_barcodes(v5_index)
    print(f"  {len(v5_df):,} ASAP v5.0.0 barcode(s) parsed "
          f"({len(v5_index) - len(v5_df):,} non-ASAP spike-in barcode(s) excluded)")

    print("[INFO] Building mapping ...")
    mapping = build_mapping(v3_df, v5_df)

    print(f"[INFO] Loading v5.0.0 pre-QC-filter barcodes from {v5_unfiltered_h5ad}")
    v5_unfiltered_df = parse_v5_barcodes(load_barcode_index(v5_unfiltered_h5ad))
    print(f"  {len(v5_unfiltered_df):,} ASAP v5.0.0 pre-QC-filter barcode(s) parsed")

    print("[INFO] Reclassifying not_mappable barcodes found in the pre-QC-filter H5AD ...")
    mapping = reclassify_filtered_cells(mapping, v3_df, v5_unfiltered_df)

    print("[INFO] Reclassifying not_mappable barcodes whose sample is absent from v5.0.0 ...")
    mapping = reclassify_missing_samples(mapping, v3_df, v5_unfiltered_df)

    print(f"[INFO] Loading v3.1.0 final-metadata barcodes from {v3_final_csv}")
    v3_final_df = parse_v3_barcodes(load_barcode_index(v3_final_csv))
    print(f"  {len(v3_final_df):,} v3.1.0 final-metadata barcode(s) parsed")
    mapping = add_in_v3_final_column(mapping, v3_final_df)

    counts = mapping["status"].value_counts()
    print(f"[INFO] {len(mapping):,} v3.1.0 barcode(s) mapped:")
    for status in (
        "1:1_vs_v5.0.0_final",
        "ambiguous_in_v5.0.0_final",
        "only_in_cleaned_unfiltered_in_v5.0.0",
        "sample_id_not_found_in_v5.0.0",
        "not_mappable",
    ):
        print(f"  {status}: {counts.get(status, 0):,}")
    print(f"  in_v3.1.0_final: {mapping['in_v3.1.0_final'].sum():,} / {len(mapping):,}")

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    out_path = outdir_path / f"{_STEM}.tsv.gz"
    mapping.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"[DONE] Written {out_path}")

    summary = build_summary_table(mapping)
    summary_path = outdir_path / f"{_STEM}_summary.md"
    summary_path.write_text(f"Barcodes from the v3.1.0 file: `{v3_seed}` mapped vs. the v5.0.0 barcodes\n\n")
    summary_path.write_text(summary_path.read_text() + summary_table_to_markdown(summary))
    print(f"[DONE] Written {summary_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--v3-seed",
        default=_V3_SEED,
        help=(
            "Path to the v3.1.0 barcode set to map (.h5ad or .csv/.csv.gz).\n"
            "Either the Strategy 1 deconvolution reference H5AD (ASAP barcodes\n"
            "plus Siletti/Kamath cells, skipped as unmatched), or the full\n"
            "--v3-final-csv, to check whether the ~20%% not_mappable rate seen\n"
            "for the reference also holds for the entire v3.1.0 cohort.\n"
            f"[Default: {_V3_SEED}]"
        ),
    )
    parser.add_argument(
        "--v5-final-csv",
        default=_V5_FINAL_CSV,
        help=f"Path to the v5.0.0 final (QC-filtered) metadata CSV.\n[Default: {_V5_FINAL_CSV}]",
    )
    parser.add_argument(
        "--v5-unfiltered-h5ad",
        default=_V5_UNFILTERED_H5AD,
        help=(
            "Path to the v5.0.0 pre-QC-filter merged_cleaned_unfiltered H5AD.\n"
            f"[Default: {_V5_UNFILTERED_H5AD}]"
        ),
    )
    parser.add_argument(
        "--v3-final-csv",
        default=_V3_FINAL_CSV,
        help=(
            "Path to v3.1.0's own final (QC-filtered) metadata CSV, used to\n"
            f"populate the in_v3.1.0_final column.\n[Default: {_V3_FINAL_CSV}]"
        ),
    )
    parser.add_argument(
        "--outdir",
        default=_OUTDIR,
        help=(
            f"Destination directory. Outputs are written as {_STEM}.tsv.gz,\n"
            f"{_STEM}_summary.md, and {_STEM}.log.\n[Default: {_OUTDIR}]"
        ),
    )
    args = parser.parse_args()

    log_run_command(Path(args.outdir) / f"{_STEM}.log")
    main(
        args.v3_seed,
        args.v5_final_csv,
        args.v5_unfiltered_h5ad,
        args.v3_final_csv,
        args.outdir,
    )

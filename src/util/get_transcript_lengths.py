#!/usr/bin/env python3

import argparse
import gzip
import io
import os
import re
import sys
from collections.abc import Iterable
from typing import TextIO, cast


ATTRIBUTE_KV_PATTERN: re.Pattern[str] = re.compile(r"(\S+)\s+\"([^\"]+)\"")


def open_maybe_gzip(path: str) -> TextIO:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r", encoding="utf-8")


def parse_gtf_attributes(attribute_field: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    matches: list[tuple[str, str]] = ATTRIBUTE_KV_PATTERN.findall(attribute_field)
    for key, value in matches:
        attributes[key] = value
    return attributes


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals_sorted = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged: list[tuple[int, int]] = []
    cur_start, cur_end = intervals_sorted[0]
    for start, end in intervals_sorted[1:]:
        if start <= cur_end + 1:
            if end > cur_end:
                cur_end = end
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    return merged


def compute_transcript_exonic_lengths(gtf_path: str) -> dict[str, int]:
    transcript_to_intervals: dict[str, list[tuple[int, int]]] = {}

    with open_maybe_gzip(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != "exon":
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            attributes = parse_gtf_attributes(parts[8])
            transcript_id = attributes.get("transcript_id")
            if not transcript_id:
                continue
            if start > end:
                start, end = end, start
            transcript_to_intervals.setdefault(transcript_id, []).append((start, end))

    transcript_to_length: dict[str, int] = {}
    for transcript_id, intervals in transcript_to_intervals.items():
        merged = merge_intervals(intervals)
        # GTF coordinates are 1-based inclusive
        length = sum(end - start + 1 for start, end in merged)
        transcript_to_length[transcript_id] = length

    return transcript_to_length


def default_output_path(input_path: str) -> str:
    base, _ = os.path.splitext(input_path)
    # handle .gtf.gz
    if base.endswith(".gtf"):
        base = base[:-4]
    return f"{base}.transcript_lengths.tsv"


def main(argv: Iterable[str]) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compute transcript lengths (sum of exon lengths, no introns) "
            "from a GTF file and write two columns: transcript_id and length."
        )
    )
    _ = parser.add_argument(
        "--gtf",
        required=False,
        default="/mnt/c/data/nih/mm39.ncbiRefSeq_SIRV.gtf",
        help="Path to input GTF (optionally .gz)",
    )
    _ = parser.add_argument(
        "--out",
        required=False,
        help="Output TSV path (defaults to <gtf without extensions>.transcript_lengths.tsv)",
    )

    args = parser.parse_args(list(argv))
    gtf_path = cast(str, args.gtf)
    out_opt = cast(str | None, args.out)
    out_path: str = out_opt or default_output_path(gtf_path)

    transcript_lengths = compute_transcript_exonic_lengths(gtf_path)

    # Write without header: two tab-separated columns
    with open(out_path, "w", encoding="utf-8") as out_fh:
        for transcript_id in sorted(transcript_lengths):
            _ = out_fh.write(f"{transcript_id}\t{transcript_lengths[transcript_id]}\n")

    print(
        f"Wrote {len(transcript_lengths)} transcripts to '{out_path}'.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))



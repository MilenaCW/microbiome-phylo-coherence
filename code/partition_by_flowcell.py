import argparse
import csv
import gzip
import re
import shutil
from collections import Counter
from pathlib import Path

FASTQ_EXTENSIONS = ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq")
UNKNOWN = "UNKNOWN"


def parse_flowcell_id(header_line: str) -> str:
    header = header_line.strip()
    if header.startswith("@"):
        header = header[1:]
    header = header.split()[0]
    parts = header.split(":")
    if len(parts) >= 3:
        return parts[2]
    return UNKNOWN


def flowcell_from_fastq(path: Path) -> str:
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            first_header = handle.readline()
    except (OSError, gzip.BadGzipFile):
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            first_header = handle.readline()
    if not first_header:
        return UNKNOWN
    return parse_flowcell_id(first_header)


def apply_flowcell_regex(flowcell_id: str, pattern: re.Pattern | None) -> str:
    if pattern is None or flowcell_id == UNKNOWN:
        return flowcell_id
    m = pattern.search(flowcell_id)
    if m and m.lastindex:
        return m.group(1)
    return flowcell_id


def collect_fastq_files(input_dir: Path) -> list[Path]:
    files: list[Path] = []
    seen: set[Path] = set()
    for pattern in FASTQ_EXTENSIONS:
        for f in sorted(input_dir.glob(pattern)):
            if f not in seen:
                files.append(f)
                seen.add(f)
    return sorted(files)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Sort FASTQ files into per-flowcell subdirectories. "
            "Flowcell IDs are read from each file's first header line "
            "(Illumina format: @instrument:run:FLOWCELL:...). "
            "Files whose headers cannot be parsed go into an 'UNKNOWN' directory."
        )
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing flat FASTQ files to sort.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Parent directory where flowcell subdirectories will be created.",
    )
    parser.add_argument(
        "--flowcell-regex",
        default=None,
        metavar="PATTERN",
        help=(
            "Regex with one capture group to extract a substring from the raw "
            "flowcell ID. If the pattern matches, group 1 is used as the directory "
            "name; otherwise the full raw ID is kept. "
            "Example for soil data (ID format '000000000-ATY9B'): '[^-]+-(.+)'"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned moves without moving any files.",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    output_dir = Path(args.output_dir)
    if not args.dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    flowcell_regex = re.compile(args.flowcell_regex) if args.flowcell_regex else None
    if flowcell_regex and not flowcell_regex.groups:
        parser.error("--flowcell-regex must contain exactly one capture group, e.g. '[^-]+-(.+)'")

    fastq_files = collect_fastq_files(input_dir)
    if not fastq_files:
        raise FileNotFoundError(
            f"No FASTQ files found in {input_dir} "
            f"(searched for: {', '.join(FASTQ_EXTENSIONS)})"
        )

    assignments: list[tuple[Path, str]] = []
    for fastq in fastq_files:
        raw_id = flowcell_from_fastq(fastq)
        flowcell_id = apply_flowcell_regex(raw_id, flowcell_regex)
        assignments.append((fastq, flowcell_id))

    counts: Counter[str] = Counter(fc for _, fc in assignments)

    print(f"Found {len(counts)} flowcell(s):")
    for flowcell, n in sorted(counts.items()):
        print(f"  {flowcell:<20} {n} file(s)")

    if args.dry_run:
        print()
        for fastq, flowcell_id in assignments:
            dest = output_dir / flowcell_id / fastq.name
            print(f"DRY RUN: {fastq} -> {dest}")
        return

    for fastq, flowcell_id in assignments:
        dest_dir = output_dir / flowcell_id
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.move(str(fastq), str(dest_dir / fastq.name))
        print(f"Moved {fastq.name} -> {flowcell_id}/")

    summary_path = output_dir / "flowcells.tsv"
    with summary_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["flowcell", "n_files"])
        for flowcell, n in sorted(counts.items()):
            writer.writerow([flowcell, n])
    print(f"\nWrote flowcell summary to {summary_path}")


if __name__ == "__main__":
    main()

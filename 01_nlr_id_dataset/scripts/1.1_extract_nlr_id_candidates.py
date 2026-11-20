#!/usr/bin/env python3
"""Step 1.1 script: extract CC-NLRs that carry an extra domain.

Reads the per-assembly NLRtracker output directories of an annotated proteome
collection and retains every protein that NLRtracker classifies with the given
status and simplified domain architecture (by default `NLR` and `CNLO`, that
is, a coiled-coil NLR carrying at least one domain in addition to CC, NB-ARC
and LRR). For each retained protein the full-length sequence, the NB-ARC
sequence, the InterProScan signatures and the NLRtracker domain table rows are
written out, with the assembly accession prefixed to every sequence identifier.

Assemblies and proteins are processed in sorted order, so the output is
reproducible regardless of the order in which the filesystem lists the input
directories.
"""

import argparse
import re
import sys
from pathlib import Path

ASSEMBLY_RE = re.compile(r"^nlrtracker_(?P<accession>GC[AF]_\d+\.\d+)_")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nlrtracker-dir", required=True, type=Path,
                   help="Directory holding the per-assembly NLRtracker output "
                        "directories, each named nlrtracker_<accession>_<suffix>.")
    p.add_argument("--output-dir", required=True, type=Path,
                   help="Directory to write the extracted candidate files into.")
    p.add_argument("--report-tsv", required=True, type=Path,
                   help="Per-assembly report of NLR and candidate counts.")
    p.add_argument("--output-prefix", default="1.1_nlr_id_candidates",
                   help="Basename shared by the output files "
                        "(default: %(default)s).")
    p.add_argument("--status", default="NLR",
                   help="Required value of the NLRtracker Status column "
                        "(default: %(default)s).")
    p.add_argument("--architecture", default="CNLO",
                   help="Required value of the NLRtracker simplified "
                        "domain-architecture column (default: %(default)s).")
    return p.parse_args()


def fail(message):
    sys.exit(f"error: {message}")


def read_fasta(path):
    """Yield (identifier, sequence) pairs. The identifier is the first word."""
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def find_assemblies(root):
    """Return {accession: directory}, sorted by accession."""
    found = {}
    for entry in sorted(root.iterdir()):
        if not entry.is_dir():
            continue
        if entry.name.startswith("._"):
            continue
        match = ASSEMBLY_RE.match(entry.name)
        if match is None:
            continue
        accession = match.group("accession")
        if accession in found:
            fail(f"assembly {accession} appears in more than one directory under {root}")
        found[accession] = entry
    if not found:
        fail(f"no nlrtracker_<accession>_* directories found under {root}")
    return found


def one_file(directory, suffix):
    # Names beginning with "._" are macOS AppleDouble sidecar files that appear
    # on FAT and exFAT volumes; they are not NLRtracker output.
    matches = sorted(path for path in directory.glob(f"*{suffix}")
                     if not path.name.startswith("._"))
    if len(matches) != 1:
        fail(f"expected exactly one *{suffix} in {directory}, found {len(matches)}")
    return matches[0]


def read_nlrtracker_table(path, status, architecture):
    """Return (n_nlr, [protein_id, ...]) for the rows matching the selection."""
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for required in ("seqname", "Status", "Domain architecture simplified"):
            if required not in header:
                fail(f"{path} has no '{required}' column")
        i_name = header.index("seqname")
        i_status = header.index("Status")
        i_arch = header.index("Domain architecture simplified")
        n_nlr, selected = 0, []
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(i_name, i_status, i_arch):
                fail(f"{path} line {line_number} has too few columns")
            if fields[i_status] != status:
                continue
            n_nlr += 1
            if fields[i_arch] == architecture:
                selected.append(fields[i_name])
    if len(set(selected)) != len(selected):
        fail(f"{path} lists the same seqname more than once among the selected rows")
    return n_nlr, sorted(selected)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_tsv.parent.mkdir(parents=True, exist_ok=True)

    prefix = args.output_prefix
    paths = {
        "protein": args.output_dir / f"{prefix}.fasta",
        "nbarc": args.output_dir / f"{prefix}_nbarc.fasta",
        "gff": args.output_dir / f"{prefix}.gff",
        "domains": args.output_dir / f"{prefix}_domains.tsv",
    }

    assemblies = find_assemblies(args.nlrtracker_dir)
    totals = {"nlr": 0, "candidates": 0, "nbarc": 0}
    domains_header_written = {"done": False}
    report = []

    with open(paths["protein"], "w") as protein_out, \
         open(paths["nbarc"], "w") as nbarc_out, \
         open(paths["gff"], "w") as gff_out, \
         open(paths["domains"], "w") as domains_out, \
         open(args.report_tsv, "w") as report_out:

        report_out.write("assembly_accession\tn_nlr\tn_candidates\n")

        for accession, directory in assemblies.items():
            table = one_file(directory, "_NLRtracker.tsv")
            n_nlr, selected = read_nlrtracker_table(table, args.status, args.architecture)
            wanted = set(selected)
            totals["nlr"] += n_nlr
            totals["candidates"] += len(selected)
            report_out.write(f"{accession}\t{n_nlr}\t{len(selected)}\n")
            if not selected:
                continue

            # Full-length protein sequences.
            sequences = {name: seq for name, seq in read_fasta(one_file(directory, "_NLR.fasta"))
                         if name in wanted}
            missing = wanted - set(sequences)
            if missing:
                fail(f"{accession}: no protein sequence for {sorted(missing)[:5]}")
            for name in selected:
                protein_out.write(f">{accession}|{name}\n{sequences[name]}\n")

            # NB-ARC sequences. A record identifier may carry a coordinate
            # suffix such as XP_000000.1_(1-300); a protein may in principle
            # have more than one NB-ARC region, and every one is kept.
            nbarc = {}
            for name, seq in read_fasta(one_file(directory, "_NBARC.fasta")):
                protein = name.split("_(")[0]
                if protein in wanted:
                    nbarc.setdefault(protein, []).append((name, seq))
            missing = wanted - set(nbarc)
            if missing:
                fail(f"{accession}: no NB-ARC sequence for {sorted(missing)[:5]}")
            for name in selected:
                for record_name, seq in sorted(nbarc[name]):
                    nbarc_out.write(f">{accession}|{record_name}\n{seq}\n")
                    totals["nbarc"] += 1

            # InterProScan signatures. The GFF seqid becomes
            # <accession>|<protein>, so every feature stays traceable to its
            # genome even if a protein identifier recurs in another assembly.
            gff_rows = {}
            with open(one_file(directory, "_NLR.gff3")) as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    name = line.split("\t", 1)[0]
                    if name in wanted:
                        gff_rows.setdefault(name, []).append(line.rstrip("\n"))
            for name in selected:
                for row in gff_rows.get(name, []):
                    gff_out.write(f"{accession}|{row}\n")

            # NLRtracker domain table. The assembly accession is added as a
            # separate first column; the file keeps NLRtracker's own header.
            domain_rows = {}
            with open(one_file(directory, "_Domains.tsv")) as handle:
                header = handle.readline().rstrip("\n")
                if not domains_header_written["done"]:
                    domains_out.write(f"assembly_accession\t{header}\n")
                    domains_header_written["done"] = True
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    name = line.split("\t", 1)[0]
                    if name in wanted:
                        domain_rows.setdefault(name, []).append(line.rstrip("\n"))
            for name in selected:
                for row in domain_rows.get(name, []):
                    domains_out.write(f"{accession}\t{row}\n")

    print(f"Step 1.1 extract_nlr_id_candidates completed: "
          f"assemblies={len(assemblies)}, nlr={totals['nlr']}, "
          f"candidates={totals['candidates']}, nbarc_records={totals['nbarc']}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

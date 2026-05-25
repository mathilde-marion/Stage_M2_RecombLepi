#!/usr/bin/env python3
# Standard library imports
import argparse
import csv
import sys
from collections import defaultdict


# ------------------------------------------------------------------
# Genetic code
# ------------------------------------------------------------------
# This dictionary maps each codon to its translated amino acid.
# "*" denotes a stop codon.
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

# Allowed unambiguous DNA bases
VALID_BASES = set("ACGT")

# For each reference base, these are the 3 possible alternative bases.
# This is used to test all possible single-nucleotide substitutions.
ALT_BASES = {
    "A": ["C", "G", "T"],
    "C": ["A", "G", "T"],
    "G": ["A", "C", "T"],
    "T": ["A", "C", "G"],
}


# ------------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------------
def parse_args():
    """
    Read arguments from the command line.

    Required:
        --fasta : genome FASTA
        --gff3  : annotation in GFF3 format

    Optional:
        --mode  : either:
            - transcript   -> output one line per transcript
            - gene_longest -> output one line per gene, using the longest CDS isoform
        --output : output TSV file path, or "-" for stdout
    """
    parser = argparse.ArgumentParser(
        description="Compute synonymous and nonsynonymous sites from genome FASTA + GFF3 CDS annotations."
    )
    parser.add_argument("--fasta", required=True, help="Genome FASTA")
    parser.add_argument("--gff3", required=True, help="Genome annotation in GFF3")
    parser.add_argument(
        "--mode",
        choices=["transcript", "gene_longest"],
        default="gene_longest",
        help="Output per transcript, or per gene using the longest CDS isoform (default: gene_longest)"
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Output TSV file (default: stdout)"
    )
    return parser.parse_args()


# ------------------------------------------------------------------
# FASTA reader
# ------------------------------------------------------------------
def read_fasta(path):
    """
    Read a FASTA file into a dictionary:
        {sequence_name: sequence}

    Only the first word after ">" is used as the sequence ID.
    Sequences are converted to uppercase.
    """
    seqs = {}
    name = None
    chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            # New FASTA header
            if line.startswith(">"):
                # Save previous sequence before starting a new one
                if name is not None:
                    seqs[name] = "".join(chunks).upper()

                name = line[1:].split()[0]
                chunks = []

            # Sequence line
            else:
                chunks.append(line)

    # Save last sequence at end of file
    if name is not None:
        seqs[name] = "".join(chunks).upper()

    return seqs


# ------------------------------------------------------------------
# GFF3 attribute parser
# ------------------------------------------------------------------
def parse_attributes(attr_str):
    """
    Parse the 9th column of a GFF3 line into a dictionary.

    Example:
        "ID=tx1;Parent=gene1"
    becomes:
        {"ID": "tx1", "Parent": "gene1"}
    """
    attrs = {}
    for field in attr_str.split(";"):
        field = field.strip()
        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key] = value
        else:
            # Rare malformed case: keep the key with empty value
            attrs[field] = ""

    return attrs


# ------------------------------------------------------------------
# Reverse complement
# ------------------------------------------------------------------
def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    Used for CDS on the minus strand.
    """
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1].upper()


# ------------------------------------------------------------------
# GFF3 parser
# ------------------------------------------------------------------
def parse_gff3(path):
    """
    Read the GFF3 and extract:
        1. transcript_to_gene:
            {transcript_id: gene_id}
        2. gene_ids:
            set of known gene IDs
        3. cds_by_parent:
            {parent_id: [list of CDS features]}

    Each CDS feature is stored as a dictionary with:
        seqid, start, end, strand, phase, attrs
    """
    transcript_to_gene = {}
    gene_ids = set()
    cds_by_parent = defaultdict(list)

    with open(path) as fh:
        for line in fh:
            # Skip comments and empty lines
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
#            if len(fields) != 9:
#                continue
            fields=fields[0:9]

            seqid, source, feature_type, start, end, score, strand, phase, attrs_str = fields
            start = int(start)
            end = int(end)
            attrs = parse_attributes(attrs_str)

            feature_type_lower = feature_type.lower()
            feature_id = attrs.get("ID")
            parents = attrs.get("Parent", "")
            parent_list = [x for x in parents.split(",") if x]

            # Record gene IDs
            if feature_type_lower == "gene":
                if feature_id:
                    gene_ids.add(feature_id)

            # Record transcript -> gene mapping
            if feature_type_lower in {"mrna", "transcript"}:
                if feature_id and parent_list:
                    transcript_to_gene[feature_id] = parent_list[0]

            # Record CDS features by parent transcript (or gene in some GFF3 files)
            if feature_type == "CDS":
                if not parent_list:
                    continue

                try:
                    phase_int = int(phase)
                except ValueError:
                    phase_int = 0

                for parent in parent_list:
                    cds_by_parent[parent].append({
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "phase": phase_int,
                        "attrs": attrs,
                    })

    return transcript_to_gene, gene_ids, cds_by_parent


# ------------------------------------------------------------------
# Reconstruct a CDS sequence from its CDS features
# ------------------------------------------------------------------
def build_cds_sequence(cds_features, genome):
    """
    Rebuild the full CDS sequence from a list of CDS features.

    Important points:
    - CDS features are ordered according to strand.
    - The GFF3 phase is used as a consistency check on frame continuity.
    - Bases at the beginning of a CDS fragment are NOT discarded.
    - Codons may span exon boundaries.
    - If the CDS is on the minus strand, fragments are reverse-complemented.
    - The final CDS length must be a multiple of 3, otherwise an error is raised.
    """
    if not cds_features:
        return None, "no_cds"

    strands = {x["strand"] for x in cds_features}
    seqids = {x["seqid"] for x in cds_features}

    if len(strands) != 1:
        return None, "mixed_strand"
    if len(seqids) != 1:
        return None, "multi_seqid"

    strand = next(iter(strands))
    seqid = next(iter(seqids))

    if seqid not in genome:
        return None, "missing_seqid_in_fasta"

    if strand == "+":
        ordered = sorted(cds_features, key=lambda x: (x["start"], x["end"]))
    elif strand == "-":
        ordered = sorted(cds_features, key=lambda x: (x["end"], x["start"]), reverse=True)
    else:
        return None, "invalid_strand"

    pieces = []
    cum_len = 0

    for i, feat in enumerate(ordered):
        start0 = feat["start"] - 1
        end1 = feat["end"]
        frag = genome[seqid][start0:end1]

        if strand == "-":
            frag = reverse_complement(frag)
        else:
            frag = frag.upper()

        phase = feat["phase"]
        if phase not in (0, 1, 2):
            return None, "invalid_phase"

        # Phase consistency check based on cumulative CDS length before this fragment
        expected_phase = (3 - (cum_len % 3)) % 3

        if i == 0:
            if phase != 0:
                sys.stderr.write(
                    "Warning: first CDS fragment has phase {} on {}; proceeding anyway.\n".format(
                        phase, seqid
                    )
                )
        else:
            if phase != expected_phase:
                return None, "phase_inconsistency_expected{}_got{}".format(
                    expected_phase, phase
                )

        pieces.append(frag)
        cum_len += len(frag)

    cds_seq = "".join(pieces).upper()

    if len(cds_seq) % 3 != 0:
        return None, "cds_length_not_multiple_of_3_len{}".format(len(cds_seq))

    return {
        "seqid": seqid,
        "strand": strand,
        "cds_seq": cds_seq,
    }, None

# ------------------------------------------------------------------
# Count synonymous and nonsynonymous sites in a CDS
# ------------------------------------------------------------------
def count_sites_in_cds(cds_seq):
    """
    Compute synonymous and nonsynonymous sites in a CDS.

    Method:
    - Split CDS into codons
    - Remove terminal stop codon if present
    - For each codon position:
        - test the 3 possible single-nucleotide substitutions
        - count how many are synonymous
        - synonymous sites at that position = (# synonymous changes) / 3
        - nonsynonymous sites at that position = (# nonsynonymous changes) / 3

    Returns a dictionary with:
        - codon_count
        - syn_sites
        - nonsyn_sites
        - ignored_ambiguous_codons
        - has_terminal_stop
        - usable_cds_length
    """
    ignored_ambiguous_codons = 0
    has_terminal_stop = False

    # Too short to contain a codon
    if len(cds_seq) < 3:
        return {
            "codon_count": 0,
            "syn_sites": 0.0,
            "nonsyn_sites": 0.0,
            "ignored_ambiguous_codons": 0,
            "has_terminal_stop": False,
            "usable_cds_length": 0,
        }

    # At this stage, CDS length should already have been validated
    if len(cds_seq) % 3 != 0:
        raise ValueError(f"CDS length is not a multiple of 3 (length={len(cds_seq)})")

    # Split into codons
    codons = [cds_seq[i:i+3] for i in range(0, len(cds_seq), 3)]

    # Remove terminal stop codon from the denominator
    if codons and codons[-1] in {"TAA", "TAG", "TGA"}:
        has_terminal_stop = True
        codons = codons[:-1]

    syn_sites = 0.0
    nonsyn_sites = 0.0
    usable_codons = 0

    for codon in codons:
        # Ignore codons containing ambiguous bases
        if len(codon) != 3 or any(base not in VALID_BASES for base in codon):
            ignored_ambiguous_codons += 1
            continue

        if codon not in CODON_TABLE:
            ignored_ambiguous_codons += 1
            continue

        aa_ref = CODON_TABLE[codon]

        # Ignore internal stop codons
        if aa_ref == "*":
            ignored_ambiguous_codons += 1
            continue

        usable_codons += 1

        # Examine each codon position independently
        for pos in range(3):
            ref_base = codon[pos]
            syn_changes = 0

            # Test all 3 possible alternative nucleotides
            for alt in ALT_BASES[ref_base]:
                mut = codon[:pos] + alt + codon[pos+1:]
                aa_alt = CODON_TABLE.get(mut, None)
                if aa_alt is None:
                    continue

                if aa_alt == aa_ref:
                    syn_changes += 1

            # Fraction of synonymous changes at this position
            syn_sites += syn_changes / 3.0

            # The rest are nonsynonymous
            nonsyn_sites += (3 - syn_changes) / 3.0

    return {
        "codon_count": usable_codons,
        "syn_sites": syn_sites,
        "nonsyn_sites": nonsyn_sites,
        "ignored_ambiguous_codons": ignored_ambiguous_codons,
        "has_terminal_stop": has_terminal_stop,
        "usable_cds_length": usable_codons * 3,
    }


# ------------------------------------------------------------------
# Choose the longest transcript for each gene
# ------------------------------------------------------------------
def choose_longest_transcript_per_gene(rows):
    """
    When several transcript isoforms exist for the same gene,
    keep only the one with the longest usable CDS.

    Input:
        rows = list of output dictionaries

    Output:
        list of one row per gene
    """
    best = {}

    for row in rows:
        gene_id = row["gene_id"]
        if gene_id is None:
            continue

        current = best.get(gene_id)

        if current is None or row["usable_cds_length"] > current["usable_cds_length"]:
            best[gene_id] = row

    return list(best.values())


# ------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------
def main():
    """
    Main workflow:
    1. Read arguments
    2. Load genome FASTA
    3. Parse GFF3
    4. Rebuild CDS for each transcript
    5. Count synonymous and nonsynonymous sites
    6. Optionally collapse to one transcript per gene
    7. Write TSV output
    """
    args = parse_args()

    # Read genome sequences
    genome = read_fasta(args.fasta)

    # Parse annotation
    transcript_to_gene, gene_ids, cds_by_parent = parse_gff3(args.gff3)

    rows = []

    # Process each transcript (or CDS parent)
    for parent_id, cds_features in cds_by_parent.items():

        # Infer gene ID and transcript ID
        if parent_id in transcript_to_gene:
            gene_id = transcript_to_gene[parent_id]
            transcript_id = parent_id

        # Some GFF3 files attach CDS directly to gene
        elif parent_id in gene_ids:
            gene_id = parent_id
            transcript_id = parent_id

        # Otherwise keep parent as transcript ID, gene unknown
        else:
            gene_id = None
            transcript_id = parent_id

        # Rebuild CDS sequence from genome + coordinates
        built, err = build_cds_sequence(cds_features, genome)
        if err is not None:
            sys.stderr.write(f"Skipping {parent_id}: {err}\n")
            continue

        # Count synonymous/nonsynonymous sites
        site_stats = count_sites_in_cds(built["cds_seq"])

        # Store output row
        row = {
            "transcript_id": transcript_id,
            "gene_id": gene_id,
            "seqid": built["seqid"],
            "strand": built["strand"],
            "usable_cds_length": site_stats["usable_cds_length"],
            "codon_count": site_stats["codon_count"],
            "syn_sites": site_stats["syn_sites"],
            "nonsyn_sites": site_stats["nonsyn_sites"],
            "total_sites": site_stats["syn_sites"] + site_stats["nonsyn_sites"],
            "ignored_ambiguous_codons": site_stats["ignored_ambiguous_codons"],
            "has_terminal_stop": int(site_stats["has_terminal_stop"]),
        }
        rows.append(row)

    # If requested, collapse to one row per gene using the longest CDS isoform
    if args.mode == "gene_longest":
        rows = choose_longest_transcript_per_gene(rows)

    # Write output
    out_fh = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    writer = csv.writer(out_fh, delimiter="\t")

    # Header
    writer.writerow([
        "gene_id",
        "transcript_id",
        "seqid",
        "strand",
        "usable_cds_length",
        "codon_count",
        "syn_sites",
        "nonsyn_sites",
        "total_sites",
        "ignored_ambiguous_codons",
        "has_terminal_stop",
    ])

    # One line per CDS / gene
    for row in rows:
        writer.writerow([
            row["gene_id"] if row["gene_id"] is not None else ".",
            row["transcript_id"],
            row["seqid"],
            row["strand"],
            row["usable_cds_length"],
            row["codon_count"],
            f"{row['syn_sites']:.6f}",
            f"{row['nonsyn_sites']:.6f}",
            f"{row['total_sites']:.6f}",
            row["ignored_ambiguous_codons"],
            row["has_terminal_stop"],
        ])

    if out_fh is not sys.stdout:
        out_fh.close()


# ------------------------------------------------------------------
# Script entry point
# ------------------------------------------------------------------
if __name__ == "__main__":
    main()

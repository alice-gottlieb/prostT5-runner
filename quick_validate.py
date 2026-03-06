from Bio import SeqIO
from collections import Counter
aa = {r.id: len(r.seq) for r in SeqIO.parse("results_fs_only_full_fs_compare/combined_aa.fasta", "fasta")}
tdi = {r.id: len(r.seq) for r in SeqIO.parse("results_fs_only_full_fs_compare/foldseek_db/all_sequences_3di.fasta", "fasta")}
for sid in aa:
    if sid not in tdi:
        print(f"MISSING: {sid}")
    elif aa[sid] != tdi[sid]:
        print(f"LENGTH MISMATCH: {sid} aa={aa[sid]} 3di={tdi[sid]}")

for r in SeqIO.parse("results_fs_only_full_fs_compare/foldseek_db/all_sequences_3di.fasta", "fasta"):
    counts = Counter(str(r.seq))
    # Flag if >90% of residues are a single character
    most_common_char, most_common_count = counts.most_common(1)[0]
    if most_common_count / len(r.seq) > 0.9:
        print(f"SUSPICIOUS: {r.id} is {most_common_count/len(r.seq):.0%} '{most_common_char}'")
    # Flag unexpected characters
    unexpected = set(str(r.seq)) - set("ACDEFGHIKLMNPQRSTVWY")
    if unexpected:
        print(f"BAD CHARS: {r.id} has {unexpected}")

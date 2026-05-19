# 3Di merged-DB output directory

ProstT5-predicted 3Di codes for all bacterial reference genomes.

## Directory Sturcture


| File | Description |
| --- | --- |
| `all_sequences_3di.fasta` | Every 3Di code in a single file |
| `validation.log` | List of proteins in `all_sequences_3di.fasta` and confirming equal length with protein FAA files from NCBI (some NCBI references genomes are missing due to lack of .faa files present on NCBI) |
| `mergedDB` (+ `.index`, `.dbtype`) | Main amino-acid sub-DB. |
| `mergedDB_h` (+ `.index`, `.dbtype`) | Headers (original FASTA sequence IDs). |
| `mergedDB_ss` (+ `.index`, `.dbtype`) | Predicted 3Di structural sequences. |
| `mergedDB_ss_h` (+ `.index`, `.dbtype`) | Header link for the 3Di sub-DB. Hard/symbolic link to `mergedDB_h` created by `foldseek lndb`. Required for `convert2fasta`. |

## Usage

**Read the 3Di sequences** — `all_sequences_3di.fasta` is a plain FASTA file, so any FASTA-aware tool works. E.g. To count entries:

```bash
grep '^>' all_sequences_3di.fasta | wc -l
``` 
or 
```python
from Bio import SeqIO
recs = list(SeqIO.parse("all_sequences_3di.fasta", "fasta"))
```

**Search the DB against itself (all-vs-all):**
```bash
foldseek search mergedDB mergedDB result_prefix tmp_search -e 0.001 --threads 16
foldseek convertalis mergedDB mergedDB result_prefix results.tsv --threads 16
```

**Search against another DB:**
```bash
foldseek search mergedDB target_db result_prefix tmp_search --threads 16
```

**Cluster:**
```bash
foldseek cluster mergedDB cluster_result tmp_cluster --threads 16
```

## Notes

- Built by `foldseek concatdbs` (no `--preserve-keys`) over many smaller GPU runs of `foldseek createdb` on H200 and RTX2080Ti cards. Then `foldseek lndb` to link `_ss_h` and `_h`. 
- Sequences came from ProstT5 prediction. `--alignment-type 1` (TMalign) is **not** supported, and TM-score / LDDT output is unavailable. The default 3Di Gotoh-Smith-Waterman alignment is fine.

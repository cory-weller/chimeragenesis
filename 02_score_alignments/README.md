# 02_score_alignments

## Score all alignments
```
python score_alignments.py
# Runs for < 1 minute
```

Create `tar` archive of alignments to save disk space
```
cd protein_alignments && \
ls | tar -czv -f alignments.tar.gz --files-from - && \
rm *.pep.aln && \
rm *.pep.aln.2
```

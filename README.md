# tmer
Clusters transcripts by shared fixed-length kmers. This tool was originally devloped to detect redundant isoforms in transcriptome assemblies.

# Usage example
To cluster transcripts in `transcripts.fasta` and write the annotated transcripts to `transcripts_clustered.fasta`, use:

`python tmer.py transcripts.fasta transcripts_clustered.fasta`

Clustered transcript sequences will be named according to `>C<cluster_index>_<isoform_index>`. 

# Notes

This tool is intended to work with transcripts in the 10kb-100kb range.

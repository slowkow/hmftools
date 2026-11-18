# Copy to env.local.sh (gitignored) and edit. Nothing patient-specific belongs in the repo.
export SAGE_TUMOR=TUMOR_ID
export SAGE_NORMAL=NORMAL_ID
# Directory holding <sample>.redux.bam(.bai), .redux.bqr.tsv, .redux.jitter_params.tsv, .redux.ms_table.tsv.gz
export SAGE_BAM_DIR=/path/to/bams
export SAGE_REF_GENOME=/path/to/GRCh38.fasta          # needs .fai and .dict alongside
export SAGE_RESOURCES=/path/to/hmf_pipeline_resources.38
export SAGE_OUT_DIR=/path/outside/the/repo/baseline   # outputs contain patient variants

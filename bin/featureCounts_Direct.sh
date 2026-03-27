BAM_DIR="/media/baseripper/volume_28TB_1/metalOutput/run_P11065/star"
OUT_DIR="/media/baseripper/volume_28TB_1/metalOutput/run_P11065/featurecounts_from_bam"
GTF="/media/baseripper/2TB_Baseripper/Stef/V0.1.1_MetaLL/MetaLL_0.0.1/assets/GRCh38_indexed/gencode.v44.annotation.gtf"
REPO="/media/baseripper/2TB_Baseripper/Stef/MetaLL-dev/MetaLL_0.0.1"

mkdir -p "$OUT_DIR" /tmp/fc_tmp
chmod 777 /tmp/fc_tmp

for bam in "$BAM_DIR"/*.Aligned.sortedByCoord.out.bam; do
  s=$(basename "$bam" .Aligned.sortedByCoord.out.bam)
  singularity exec --env TMPDIR=/tmp/fc_tmp docker://quay.io/biocontainers/subread:2.0.6--he4a0461_0 \
    featureCounts -T 8 -p -s 0 --tmpDir /tmp/fc_tmp -a "$GTF" -o "$OUT_DIR/${s}.featureCounts.txt" "$bam"
done

python3 "$REPO/bin/merge_counts_enhanced.py" \
  --featurecounts "$OUT_DIR"/*.featureCounts.txt \
  --output "/media/baseripper/volume_28TB_1/metalOutput/run_P11065/counts/featurecounts_gene_counts_matrix.tsv"
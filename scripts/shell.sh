### Importing using Manifest

qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/ms_manifest.tsv \
  --output-path demux_seqs.qza

### Denoising the Data

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 151 \
  --o-representative-sequences MS_rep-seqs.qza \
  --o-table MS_table.qza \
  --o-denoising-stats MS_stats.qza

qiime metadata tabulate \
  --m-input-file MS_stats.qza \
  --o-visualization MS_stats.qzv

qiime feature-table summarize \
  --i-table MS_table.qza \
  --o-visualization MS_table.qzv \
  --m-sample-metadata-file /work/group_assignment_1/data/MS/corrected_ms_metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data MS_rep-seqs.qza \
  --o-visualization MS_rep-seqs.qzv 
  

### Building the Classifier

qiime feature-classifier extract-reads \
  --i-sequences /datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --o-reads ref-seqs-trimmed.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed.qza \
  --i-reference-taxonomy /datasets/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads MS_rep-seqs.qza \
  --o-classification MS_taxonomy.qza

qiime metadata tabulate \
  --m-input-file MS_taxonomy.qza \
  --o-visualization MS_taxonomy.qzv

qiime taxa barplot \
  --i-table MS_table.qza \
  --i-taxonomy MS_taxonomy.qza \
  --m-metadata-file /work/group_assignment_1/data/MS/corrected_ms_metadata.tsv \
  --o-visualization MS_taxa-bar-plots.qzv


### Alpha Rarefaction

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences MS_rep-seqs.qza \
  --o-alignment MS_aligned-rep-seqs.qza \
  --o-masked-alignment MS_masked-aligned-rep-seqs.qza \
  --o-tree MS_unrooted-tree.qza \
  --o-rooted-tree MS_rooted-tree.qza 

qiime diversity alpha-rarefaction \
  --i-table MS_table.qza \
  --i-phylogeny MS_rooted-tree.qza \
  --p-max-depth 27052 \
  --m-metadata-file /work/group_assignment_1/data/MS/corrected_ms_metadata.tsv \
  --o-visualization MS_alpha-rarefaction.qzv 


### PICRUSt2 Export

qiime tools export \
   --input-path MS_table.qza \
   --output-path picrust

qiime tools export \
   --input-path MS_rep-seqs.qza \
   --output-path picrust

conda deactivate
conda activate picrust2

picrust2_pipeline.py \
   -s picrust/dna-sequences.fasta \
   -i picrust/feature-table.biom \
   -o picrust_out

conda deactivate
conda activate qiime2-amplicon-2025.4

qiime tools export \
--input-path MS_taxonomy.qza \
--output-path exports 

qiime tools export \
--input-path MS_rooted-tree.qza \
--output-path exports

biom convert \
-i picrust/feature-table.biom \
--to-tsv \
-o exports/feature-table.txt

gunzip picrust_out/combined_EC_predicted.tsv.gz
gunzip picrust_out/combined_KO_predicted.tsv.gz
gunzip picrust_out/pathways_out/path_abun_unstrat.tsv.gz

mv picrust_out/*.tsv exports
mv picrust_out/pathways_out/*.tsv exports


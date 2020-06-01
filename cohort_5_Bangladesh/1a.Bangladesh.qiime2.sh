# Install QIIME2
# conda create -n qiime2-2019.7 qiime=2019.7.0
conda activate qiime2-2019.7

# https://docs.qiime2.org/2018.6/tutorials/overview/
# https://chmi-sops.github.io/mydoc_qiime2.html

# Note: All QIIME 2 visualizers will generate a .qzv file. You can view these files with qiime tools view at https://view.qiime2.org/

# Prepare the metadata and tabulate it (to veryfy that all is fine)
# often qiime complains the file has a wrong encoding. In that case run this: iconv -t UTF-8 -f ISO-8859-1 in.tsv > out.tsv
qiime metadata tabulate \
 --m-input-file Bangladesh_metadata.tsv \
 --o-visualization Bangladesh_metadata.qzv


# These files are too large for GitHub:
# qiime tools import \
 # --type SampleData[JoinedSequencesWithQuality] \
 # --input-path Bangladesh_fastq_manifest_no-C9-C10.csv \
 # --output-path Bangladesh_reads.qza \
 # --input-format SingleEndFastqManifestPhred33

# qiime demux summarize \
 # --i-data Bangladesh_reads.qza \
 # --o-visualization Bangladesh_reads.qzv

# Since the reads are already merged, at this stage you can choose to proceed using Deblur for additional quality control, or you can dereplicate sequences and optionally cluster them into OTUs with q2-vsearch. If you try this option, we strongly encourage you to call qiime quality-filter q-score-joined with a higher min-quality threshold - possibly \
 --p-min-quality 20 or  --p-min-quality 30.
# qiime quality-filter q-score-joined \
 # --p-min-quality 30 \
 # --i-demux Bangladesh_reads.qza \
 # --o-filtered-sequences Bangl.QF30.qza \
 # --o-filter-stats Bangl.QF30.stats.qza



#######################################  deblur start #################################################

# You should pass the sequence length value you selected from the quality score plots for \
 --p-trim-length. This will trim all sequences to this length, and discard any sequences which are not at least this long. Use a trim length of based on the quality score plots
# Andrej: this plot is not helpful. Will get the reads length distribution differently:
# zcat *SeqPrep-merged.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > Bangladesh.readlenght.txt

Rscript -e 'reads<-read.csv(file="Bangladesh.readlenght.txt", sep="", header=FALSE); jpeg("Bangladesh.readlength.jpg" , width=1000, height=1000); plot (reads$V1,reads$V2,type="l",xlab="read length",ylab="occurences",col="blue"); dev.off()'

# 95% of reads are at least 344 nt long.

# Note that deblur (and also vsearch dereplicate-sequences) should be preceded by basic quality-score-based filtering, but this is unnecessary for dada2. Both Deblur and DADA2 contain internal chimera checking methods and abundance filtering, so additional filtering should not be necessary following these methods.
# try deblur on filtered data (qual 30 and 20):

qiime deblur denoise-16S \
 --i-demultiplexed-seqs Bangl.QF30.qza \
 --p-trim-length 344 \
 --p-sample-stats \
 --o-representative-sequences Bangl.QF30.deblur_rep-seqs.qza \
 --o-table Bangl.QF30.deblur_tab.qza \
 --o-stats Bangl.QF30.deblur_stats.qza
qiime feature-table summarize \
 --i-table Bangl.QF30.deblur_tab.qza \
 --o-visualization Bangl.QF30.deblur_tab.qzv \
 --m-sample-metadata-file Bangladesh_metadata.tsv 
qiime feature-table tabulate-seqs \
  --i-data Bangl.QF30.deblur_rep-seqs.qza \
  --o-visualization Bangl.QF30.deblur_rep-seqs.qzv

qiime deblur denoise-16S \
 --i-demultiplexed-seqs Bangl.QF20.qza \
 --p-trim-length 344 \
 --p-sample-stats \
 --o-representative-sequences Bangl.QF20.deblur_rep-seqs.qza \
 --o-table Bangl.QF20.deblur_tab.qza \
 --o-stats Bangl.QF20.deblur_stats.qza
qiime feature-table summarize \
 --i-table Bangl.QF20.deblur_tab.qza \
 --o-visualization Bangl.QF20.deblur_tab.qzv \
 --m-sample-metadata-file Bangladesh_metadata.tsv 
qiime feature-table tabulate-seqs \
  --i-data Bangl.QF20.deblur_rep-seqs.qza \
  --o-visualization Bangl.QF20.deblur_rep-seqs.qzv


#######################################  deblur end #################################################


#######################################  qiime feature-classifier start #######################################

# input is from deblur denoise QualFilt-p30
qiime feature-classifier classify-sklearn \
 --p-n-jobs 2 \
 --i-classifier ../Greengenes_QIIME2/gg_13_8.99_otus.extract-from-primers_ref-seqs.fit-classifier-naive-bayes.qza \
 --i-reads Bangl.QF30.deblur_rep-seqs.qza \
 --o-classification Bangl.QF30.deblur_rep-seqs.sklearn.qza
qiime metadata tabulate \
 --m-input-file Bangl.QF30.deblur_rep-seqs.sklearn.qza \
 --o-visualization Bangl.QF30.deblur_rep-seqs.sklearn.qzv

# input is from deblur denoise QualFilt-p20
qiime feature-classifier classify-sklearn \
 --p-n-jobs 2 \
 --i-classifier ../Greengenes_QIIME2/gg_13_8.99_otus.extract-from-primers_ref-seqs.fit-classifier-naive-bayes.qza \
 --i-reads Bangl.QF20.deblur_rep-seqs.qza \
 --o-classification Bangl.QF20.deblur_rep-seqs.sklearn.qza
qiime metadata tabulate \
 --m-input-file Bangl.QF20.deblur_rep-seqs.sklearn.qza \
 --o-visualization Bangl.QF20.deblur_rep-seqs.sklearn.qzv

#######################################  qiime feature-classifier end #######################################



#######################################  taxonomic analyses start #######################################  

# Taxa collapse (will need it downstream).
# https://docs.qiime2.org/2018.8/plugins/available/taxa/collapse/
# Collapse groups of features that have the same taxonomic assignment through the specified level. The frequencies of all features will be summed when they are collapsed.
# Optional. You can visualize all levels in qiime taxa barplot on the whole dataset.

# QualFilt-p30
qiime taxa collapse \
 --i-table Bangl.QF30.deblur_tab.qza\
 --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 2 \
 --o-collapsed-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-2.qza

qiime taxa collapse \
 --i-table Bangl.QF30.deblur_tab.qza\
 --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 3 \
 --o-collapsed-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.qza

qiime taxa collapse \
 --i-table Bangl.QF30.deblur_tab.qza\
 --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 4 \
 --o-collapsed-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.qza

qiime taxa collapse \
 --i-table Bangl.QF30.deblur_tab.qza\
 --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 5 \
 --o-collapsed-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.qza

qiime taxa collapse \
 --i-table Bangl.QF30.deblur_tab.qza\
 --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 6 \
 --o-collapsed-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.qza

qiime taxa collapse \
 --i-table Bangl.QF30.deblur_tab.qza\
 --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 7 \
 --o-collapsed-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.qza


# barplot (https://docs.qiime2.org/2018.8/plugins/available/taxa/barplot/)
qiime taxa barplot \
 --i-table Bangl.QF30.deblur_tab.qza \
 --i-taxonomy  Bangl.QF30.deblur_rep-seqs.sklearn.qza \
 --m-metadata-file Bangladesh_metadata.tsv \
 --o-visualization Bangl.QF30.deblur_rep-seqs.sklearn.taxa-barplot.qzv
 

# Heatmap (https://docs.qiime2.org/2018.8/plugins/available/feature-table/heatmap/?highlight=taxa%20collapse)
# QualFilt-p30
qiime feature-table heatmap \
 --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-2.qza \
 --m-metadata-file Bangladesh_metadata.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Bangl.QF30.deblur_tab.sklearn.collapsed-lev-2.heatmap.qzv

qiime feature-table heatmap \
 --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.qza \
 --m-metadata-file Bangladesh_metadata.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.heatmap.qzv

qiime feature-table heatmap \
 --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.qza \
 --m-metadata-file Bangladesh_metadata.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.heatmap.qzv

qiime feature-table heatmap \
 --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.qza \
 --m-metadata-file Bangladesh_metadata.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.heatmap.qzv

qiime feature-table heatmap \
 --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.qza \
 --m-metadata-file Bangladesh_metadata.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.heatmap.qzv

qiime feature-table heatmap \
 --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.qza \
 --m-metadata-file Bangladesh_metadata.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.heatmap.qzv



#######################################  Sequence alignment and phylogeny building  #######################################
# If you are sequencing phylogenetic markers (e.g., 16S rRNA genes), you can align these sequences to assess the phylogenetic relationship between each of your features. This phylogeny can then be used by other downstream analyses, such as UniFrac distance analyses.
# https://forum.qiime2.org/t/q2-phylogeny-community-tutorial/4455

# Sequence Alignment
qiime alignment mafft \
 --i-sequences Bangl.QF30.deblur_rep-seqs.qza \
 --o-alignment Bangl.QF30.deblur_rep-seqs.aligned.qza

qiime alignment mask \
 --i-alignment Bangl.QF30.deblur_rep-seqs.aligned.qza \
 --o-masked-alignment Bangl.QF30.deblur_rep-seqs.aligned.masked.qza

# fragment insertion analysis (https://github.com/biocore/q2-fragment-insertion)
# QIIME 2's "Moving Pictures" tutorial suggests constructing a de-novo phylogeny for the fragments, i.e FeatureData[Sequence], to obtain a Phylogeny[Rooted] that can be used for phylogenetic diversity computation. "Fragment insertion" provides an alternative way to acquire the Phylogeny[Rooted] by inserting sequences of FeatureData[Sequence] into a high quality reference phylogeny and thus provides multiple advantages over de-novo phylogenies, e.g. accurate branch lengths, multi-study meta-analyses, mixed region meta-analyses (e.g. V4 and V2).
# Fragment insertion avoids artificially long outgroup branches that would lead to exaggerated separation in beta diversity.
# Fragment insertion enables meta-analyses across different variable 16S regions and fragment length.
# Default reference (phylogeny and matching alignment) is Greengenes 13_8 at 99%.

qiime fragment-insertion sepp \
  --p-threads 4 \
  --i-representative-sequences Bangl.QF30.deblur_rep-seqs.qza \
  --o-tree Bangl.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --o-placements Bangl.QF30.deblur_rep-seqs.fragment-insertion-placements.qza


# After you created the insertion tree, which can be used for phylogenetic diversity computation, e.g. Faith's PD or UniFrac, a typical next step is to filter a feature-table such that it only contains fragments that are in the insertion tree. This becomes necessary, since SEPP might reject insertion of fragments that are too remotely related to everything in the reference alignment/phylogeny. Those rows in your feature-table for fragments not in the phylogeny, will cause diversity computation to fail, since branch lengths cannot be determined.
# In my case, nothing got removed. Anyway, this step may be useful for troubleshooting.

qiime fragment-insertion filter-features \
  --i-table Bangl.QF30.deblur_tab.qza \
  --i-tree Bangl.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --o-filtered-table Bangl.QF30.deblur_tab.fragment-insertion_filtered_table.qza \
  --o-removed-table Bangl.QF30.deblur_tab.fragment-insertion_removed_table.qza

qiime feature-table summarize \
 --i-table Bangl.QF30.deblur_tab.fragment-insertion_filtered_table.qza \
 --o-visualization Bangl.QF30.deblur_tab.fragment-insertion_filtered_table.qzv

qiime feature-table summarize \
 --i-table Bangl.QF30.deblur_tab.fragment-insertion_removed_table.qza \
 --o-visualization Bangl.QF30.deblur_tab.fragment-insertion_removed_table.qzv

#######################################  Sequence alignment and phylogeny building  end #######################################


##### Alpha and beta diversity analysis

# Alpha rarefaction plotting
# Note: The value that you provide for --p-max-depth should be determined by reviewing the “Frequency per sample” information presented in the table.qzv file that was created above. In general, choosing a value that is somewhere around the median frequency seems to work well, but you may want to increase that value if the lines in the resulting rarefaction plot don’t appear to be leveling out, or decrease that value if you seem to be losing many of your samples due to low total frequencies closer to the minimum sampling depth than the maximum sampling depth.

qiime diversity alpha-rarefaction \
  --i-table Bangl.QF30.deblur_tab.fragment-insertion_filtered_table.qza \
  --i-phylogeny Bangl.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --p-max-depth 10000 \
  --p-steps 30 \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization Bangl.QF30.deblur_tab.fragment-insertion_filtered_table.alpha-rarefaction.qzv

# Results:
# Based on the rarefaction curves, I can use all the samples. Perhaps, I could remove only one sample, P4-1 that contains 1,700 (p30) reads, so the second most abundant sample is C6 with 2,313 (p30) reads.
# I will use all samples, it is better to have all the replicates than coverage.


qiime diversity core-metrics-phylogenetic \
  --i-phylogeny Bangl.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --i-table Bangl.QF30.deblur_tab.fragment-insertion_filtered_table.qza \
  --p-sampling-depth 1700 \
  --m-metadata-file Bangladesh_metadata.tsv \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv


# Check if continuous sample metadata columns are correlated with alpha diversity:
qiime diversity alpha-correlation \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-correlation-spearman.qzv

qiime diversity alpha-correlation \
  --p-method pearson \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-correlation-pearson.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization core-metrics-results/shannon-correlation-spearman.qzv

qiime diversity alpha-correlation \
  --p-method pearson \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization core-metrics-results/shannon-correlation-pearson.qzv





# The following commands will test whether distances between samples within a group, such as samples from the same body site, are more similar to each other then they are to samples from the other groups.
# --p-pairwise parameter will also perform pairwise tests that will allow you to determine which specific pairs of groups differ from one another, if any.
# This command can be slow to run, so we’ll run this on specific columns of metadata that we’re interested in exploring, rather than all metadata columns that it’s applicable to.
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization core-metrics-results/unweighted-unifrac_GroupMonth_significance.qzv \
  --p-pairwise


qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/unweighted-unifrac_Group_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/bray_curtis_Group_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization core-metrics-results/bray_curtis_GroupMonth_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/jaccard_Group_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization core-metrics-results/jaccard_GroupMonth_significance.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/weighted-unifrac_Group_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization core-metrics-results/weighted-unifrac_GroupMonth_significance.qzv \
  --p-pairwise


# While our core-metrics-phylogenetic command did already generate some Emperor plots, we want to pass an optional parameter, --p-custom-axes, which is very useful for exploring time series data
# We will generate Emperor plots for unweighted UniFrac and Bray-Curtis so that the resulting plot will contain axes for principal coordinate 1, principal coordinate 2, and days since the experiment start. We will use that last axis to explore how these samples changed over time.
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-custom-axes Month_num \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor_Month.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-custom-axes Month_num \
  --o-visualization core-metrics-results/weighted-unifrac-emperor_Month.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-custom-axes Month_num \
  --o-visualization core-metrics-results/bray-curtis-emperor-Month.qzv

##### Differential taxonomical abundance testing
mkdir Differential_abundance

# Differential abundance testing in microbiome analysis is an active area of research. There are two QIIME 2 plugins that can be used for this: q2-gneiss and q2-composition.
# ANCOM can be applied to identify features that are differentially abundant (i.e. present in different abundances) across sample groups. ANCOM is implemented in the q2-composition plugin.
# ANCOM operates on a FeatureTable[Composition] QIIME 2 artifact, which is based on frequencies of features on a per-sample basis, but cannot tolerate frequencies of zero.
# ANCOM results: https://forum.qiime2.org/t/how-to-interpret-ancom-results/1958/

# Filtering
# check here https://forum.qiime2.org/t/gneiss-zero-balance-error/1857/5
# I should filter the deblur table for features that appear in say 5 or more samples, as singletons or doubletons can affect differential taxonomical abundance testing, including ANCOM (I tried ANCOM with the original deblur table, I did not get any meaningful results).
# Remove all features with a total abundance (summed across all samples) of less than 10 as follows: --p-min-frequency 10
# Remove features that are not present in at least 5 samples: --p-min-samples 5

qiime feature-table filter-features \
  --i-table Bangl.QF30.deblur_tab.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza

mkdir Differential_abundance/ANCOM

qiime composition add-pseudocount \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.qza

qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv

# We’re also often interested in performing a differential abundance test at a specific taxonomic level. 

qiime feature-table filter-features \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv

qiime feature-table filter-features \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv

qiime feature-table filter-features \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv

qiime feature-table filter-features \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv

qiime feature-table filter-features \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv



# gneiss (https://docs.qiime2.org/2018.8/tutorials/gneiss/)
# make sure to use the filtered table else the results will be biased, or you'll run into errors (this happened on deblur unfiltered table: https://forum.qiime2.org/t/gneiss-zero-balance-error/1857/5)

# gneiss Correlation clustering
mkdir Differential_abundance/gneiss

# gneiss correlation-clustering
qiime gneiss correlation-clustering \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --o-clustering Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza

qiime gneiss ilr-hierarchical \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --o-balances Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza

qiime gneiss ols-regression \
  --p-formula "Subject+Group+GroupMonth" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Subject+Group+GroupMonth.qzv
  
qiime gneiss dendrogram-heatmap \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_heatmap_GroupMonth.qzv

qiime gneiss dendrogram-heatmap \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Subject \
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_heatmap_Subject.qzv

qiime gneiss balance-taxonomy \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 2 \
  --p-balance-name 'y0' \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y0_taxa-level-2_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y0' \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y0_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y1' \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y1_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y2' \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y2_taxa-level-5_summary.qzv

# Try different formulas in the ols-regression
qiime gneiss ols-regression \
  --p-formula "GroupMonth" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_GroupMonth.qzv

qiime gneiss ols-regression \
  --p-formula "Group*Month_cat" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Group-x-Month.qzv

qiime gneiss ols-regression \
  --p-formula "Group+Month_cat" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Group+Month.qzv

qiime gneiss ols-regression \
  --p-formula "Subject" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Subject.qzv

qiime gneiss ols-regression \
  --p-formula "Subject+Month_cat" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Subject+Month.qzv

# Gneiss Gradient clustering
# Use only patients (it seems that it's better that the metadata too contains only patients)
qiime gneiss gradient-clustering \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --m-gradient-file Bangladesh_metadata_patients.tsv \
  --m-gradient-column Month_num \
  --o-clustering Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza

qiime gneiss ilr-hierarchical \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --o-balances Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_balances.qza

qiime gneiss ols-regression \
  --p-formula "Subject+Month_cat" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata_patients.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_regression_summary_Subject+Month.qzv

qiime gneiss ols-regression \
  --p-formula "Month_cat" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata_patients.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_regression_summary_Month.qzv

qiime gneiss ols-regression \
  --p-formula "Subject" \
  --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata_patients.tsv \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_regression_summary_Subject.qzv

qiime gneiss dendrogram-heatmap \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Bangladesh_metadata_patients.tsv \
  --m-metadata-column Month_cat \
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_heatmap_Month.qzv

qiime gneiss balance-taxonomy \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 2 \
  --p-balance-name 'y1' \
  --m-metadata-file Bangladesh_metadata_patients.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_y1_taxa-level-2_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 4 \
  --p-balance-name 'y1' \
  --m-metadata-file Bangladesh_metadata_patients.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_y1_taxa-level-4_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y1' \
  --m-metadata-file Bangladesh_metadata_patients.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_y1_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 6 \
  --p-balance-name 'y1' \
  --m-metadata-file Bangladesh_metadata_patients.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.patients.gneiss_gradient_y1_taxa-level-6_summary.qzv


# gneiss Phylogenetic analysis
# Results are similar as with correlation-clustering, but some steps don't work. No need to pursue it further.

# qiime gneiss assign-ids \
  # --i-input-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  # --i-input-tree Bangl.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  # --o-output-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids.qza \
  # --o-output-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_hierarchy.qza

# qiime gneiss ilr-hierarchical \
  # --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  # --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_hierarchy.qza \
  # --o-balances Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_balances.qza

# qiime gneiss ols-regression \
  # --p-formula "Subject+Group+GroupMonth" \
  # --i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_balances.qza \
  # --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_hierarchy.qza \
  # --m-metadata-file Bangladesh_metadata.tsv \
  # --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_regression_summary_Subject+Group+GroupMonth.qzv

# qiime gneiss dendrogram-heatmap \
  # --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  # --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_hierarchy.qza \
  # --m-metadata-file Bangladesh_metadata.tsv \
  # --m-metadata-column GroupMonth \
  # --p-color-map seismic \
  # --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_heatmap_GroupMonth.qzv

# qiime gneiss balance-taxonomy \
  # --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  # --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_hierarchy.qza \
  # --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  # --p-taxa-level 3 \
  # --p-balance-name '210L-6b39ef51-8216-4f46-a41a-4914810d0a67' \
  # --m-metadata-file Bangladesh_metadata.tsv \
  # --m-metadata-column GroupMonth \
  # --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_0210L-6b39ef51-8216-4f46-a41a-4914810d0a67_taxa-level-3_summary.qzv

# qiime gneiss balance-taxonomy \
  # --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  # --i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_hierarchy.qza \
  # --i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  # --p-taxa-level 5 \
  # --p-balance-name '210L-6b39ef51-8216-4f46-a41a-4914810d0a67' \
  # --m-metadata-file Bangladesh_metadata.tsv \
  # --m-metadata-column GroupMonth \
  # --o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_assign-ids_0210L-6b39ef51-8216-4f46-a41a-4914810d0a67_taxa-level-5_summary.qzv

# In theory we don't need assign-ids above, as the new gneiss has the ilr-phylogenetic option to deal with the rooted tree. But with this option I get an error at the ols-regression step :(
# The heatmaps look the same, so I guess it makes no difference.
#qiime gneiss ilr-phylogenetic \
  #--i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  #--i-tree Bangl.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  #--o-balances Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_balances.qza \
  #--o-hierarchy Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_hierarchy.qza

## I'm getting this error here: 'the label [y105] is not in the [index]'
#qiime gneiss ols-regression \
  #--p-formula "Subject+Group+GroupMonth" \
  #--i-table Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_balances.qza \
  #--i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_hierarchy.qza \
  #--m-metadata-file Bangladesh_metadata.tsv \
  #--o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_regression_summary_Subject+Group+GroupMonth.qzv

#qiime gneiss dendrogram-heatmap \
  #--i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  #--i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_hierarchy.qza \
  #--m-metadata-file Bangladesh_metadata.tsv \
  #--m-metadata-column GroupMonth \
  #--p-color-map seismic \
  #--o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_heatmap_GroupMonth.qzv

#qiime gneiss balance-taxonomy \
  #--i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  #--i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_hierarchy.qza \
  #--i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  #--p-taxa-level 3 \
  #--p-balance-name 'y154' \
  #--m-metadata-file Bangladesh_metadata.tsv \
  #--m-metadata-column GroupMonth \
  #--o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_y154_taxa-level-3_summary.qzv

#qiime gneiss balance-taxonomy \
  #--i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  #--i-tree Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_hierarchy.qza \
  #--i-taxonomy Bangl.QF30.deblur_rep-seqs.sklearn.qza \
  #--p-taxa-level 5 \
  #--p-balance-name 'y154' \
  #--m-metadata-file Bangladesh_metadata.tsv \
  #--m-metadata-column GroupMonth \
  #--o-visualization Differential_abundance/gneiss/Bangl.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_y154_taxa-level-5_summary.qzv

##################################
# Longitudinal analysis (https://docs.qiime2.org/2018.8/tutorials/longitudinal/)
mkdir Differential_abundance/longitudinal

# Volatility analysis
qiime longitudinal volatility \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-file core-metrics-results/shannon_vector.qza \
  --m-metadata-file core-metrics-results/evenness_vector.qza \
  --m-metadata-file core-metrics-results/faith_pd_vector.qza \
  --p-default-metric shannon \
  --p-default-group-column Group \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --o-visualization Differential_abundance/longitudinal/volatility.qzv

# Feature volatility analysis
# Using the option --p-where "Group='patient'" resulted in error, and single quoting 'Group' resulted in an empty table.
# "group" is a reserved keyword in SQL. Apparently, single quotes don't help in qiime2, but square bracket work, like this: --p-where "[Group]='patient'"
qiime feature-table filter-samples \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-where "[Group]='patient'" \
  --o-filtered-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza
qiime longitudinal feature-volatility \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_deblur_tab.freq10-minSamp5.patients


qiime feature-table filter-samples \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-where "[Group]='patient'" \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.patients.qza
qiime longitudinal feature-volatility \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_collapsed-lev-4.freq10-minSamp5.patients


qiime feature-table filter-samples \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-where "[Group]='patient'" \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.patients.qza
qiime longitudinal feature-volatility \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_collapsed-lev-5.freq10-minSamp5.patients


qiime feature-table filter-samples \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-where "[Group]='patient'" \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.qza
qiime longitudinal feature-volatility \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_collapsed-lev-6.freq10-minSamp5.patients

qiime feature-table filter-samples \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-where "[Group]='patient'" \
  --o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.patients.qza
qiime longitudinal feature-volatility \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_collapsed-lev-7.freq10-minSamp5.patients



# Non-parametric microbial interdependence test (NMIT)
############ Not working because I don't have more than 1 longitudinal group to compare (I only have "patients").
# run it on the filtered table, but remove controls from it:
# It's recommended to use a collapsed table:
# Using the option --p-where "Group='patient'" resulted in error, and single quoting 'Group' resulted in an empty table.
# "group" is a reserved keyword in SQL. Apparently, single quotes don't help in qiime2, but square bracket work, like this: --p-where "[Group]='patient'"
#qiime feature-table filter-samples \
  #--i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.qza \
  #--m-metadata-file Bangladesh_metadata.tsv \
  #--p-where "[Group]='patient'" \
  #--o-filtered-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.qza

##nmit needs relative frequency table as input:
#qiime feature-table relative-frequency \
  #--i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.qza \
  #--o-relative-frequency-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.rel-freq.qza

#qiime longitudinal nmit \
  #--i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.rel-freq.qza \
  #--m-metadata-file Bangladesh_metadata.tsv \
  #--p-individual-id-column Subject \
  #--p-corr-method pearson \
  #--o-distance-matrix Differential_abundance/longitudinal/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.rel-freq.nmit-dm.qza

## perform PERMANOVA tests to evaluate whether between-group distances are larger than within-group distance.
## It doesn't work: All values in the grouping vector are unique. This method cannot operate on a grouping vector with only unique values (e.g., there are no 'within' distances because each group of objects contains only a single object).
#qiime diversity beta-group-significance \
  #--i-distance-matrix Differential_abundance/longitudinal/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.rel-freq.nmit-dm.qza \
  #--m-metadata-file Bangladesh_metadata.tsv \
  #--m-metadata-column Subject \
  #--o-visualization Differential_abundance/longitudinal/Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.rel-freq.nmit-dm.beta-group-significance.qzv

#Predicting continuous (i.e., numerical) sample data (https://docs.qiime2.org/2018.8/tutorials/sample-classifier/)
qiime sample-classifier regress-samples \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --output-dir Differential_abundance/sample-classifier/regressor_deblur_tab.freq10-minSamp5.patients
qiime sample-classifier regress-samples-ncv \
  --i-table Bangl.QF30.deblur_tab.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --o-predictions Differential_abundance/sample-classifier/regressor_deblur_tab.freq10-minSamp5.patients/predictions-ncv.qza \
  --o-feature-importance Differential_abundance/sample-classifier/regressor_deblur_tab.freq10-minSamp5.patients/importance-ncv.qza
qiime sample-classifier scatterplot \
  --i-predictions Differential_abundance/sample-classifier/regressor_deblur_tab.freq10-minSamp5.patients/predictions-ncv.qza \
  --m-truth-file Bangladesh_metadata.tsv \
  --m-truth-column Month_num \
  --o-visualization Differential_abundance/sample-classifier/regressor_deblur_tab.freq10-minSamp5.patients/scatter.qzv

qiime sample-classifier regress-samples \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --output-dir Differential_abundance/sample-classifier/regressor_sklearn.collapsed-4.freq10-minSamp5.patients
qiime sample-classifier regress-samples-ncv \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --o-predictions Differential_abundance/sample-classifier/regressor_sklearn.collapsed-4.freq10-minSamp5.patients/predictions-ncv.qza \
  --o-feature-importance Differential_abundance/sample-classifier/regressor_sklearn.collapsed-4.freq10-minSamp5.patients/importance-ncv.qza
qiime sample-classifier scatterplot \
  --i-predictions Differential_abundance/sample-classifier/regressor_sklearn.collapsed-4.freq10-minSamp5.patients/predictions-ncv.qza \
  --m-truth-file Bangladesh_metadata.tsv \
  --m-truth-column Month_num \
  --o-visualization Differential_abundance/sample-classifier/regressor_sklearn.collapsed-4.freq10-minSamp5.patients/scatter.qzv


qiime sample-classifier regress-samples \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --output-dir Differential_abundance/sample-classifier/regressor_sklearn.collapsed-5.freq10-minSamp5.patients
qiime sample-classifier regress-samples-ncv \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --o-predictions Differential_abundance/sample-classifier/regressor_sklearn.collapsed-5.freq10-minSamp5.patients/predictions-ncv.qza \
  --o-feature-importance Differential_abundance/sample-classifier/regressor_sklearn.collapsed-5.freq10-minSamp5.patients/importance-ncv.qza
qiime sample-classifier scatterplot \
  --i-predictions Differential_abundance/sample-classifier/regressor_sklearn.collapsed-5.freq10-minSamp5.patients/predictions-ncv.qza \
  --m-truth-file Bangladesh_metadata.tsv \
  --m-truth-column Month_num \
  --o-visualization Differential_abundance/sample-classifier/regressor_sklearn.collapsed-5.freq10-minSamp5.patients/scatter.qzv


qiime sample-classifier regress-samples \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --output-dir Differential_abundance/sample-classifier/regressor_sklearn.collapsed-6.freq10-minSamp5.patients
qiime sample-classifier regress-samples-ncv \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --o-predictions Differential_abundance/sample-classifier/regressor_sklearn.collapsed-6.freq10-minSamp5.patients/predictions-ncv.qza \
  --o-feature-importance Differential_abundance/sample-classifier/regressor_sklearn.collapsed-6.freq10-minSamp5.patients/importance-ncv.qza
qiime sample-classifier scatterplot \
  --i-predictions Differential_abundance/sample-classifier/regressor_sklearn.collapsed-6.freq10-minSamp5.patients/predictions-ncv.qza \
  --m-truth-file Bangladesh_metadata.tsv \
  --m-truth-column Month_num \
  --o-visualization Differential_abundance/sample-classifier/regressor_sklearn.collapsed-6.freq10-minSamp5.patients/scatter.qzv


qiime sample-classifier regress-samples \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --output-dir Differential_abundance/sample-classifier/regressor_sklearn.collapsed-7.freq10-minSamp5.patients
qiime sample-classifier regress-samples-ncv \
  --i-table Bangl.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.patients.qza \
  --m-metadata-file Bangladesh_metadata.tsv \
  --m-metadata-column Month_num \
  --p-estimator RandomForestRegressor \
  --o-predictions Differential_abundance/sample-classifier/regressor_sklearn.collapsed-7.freq10-minSamp5.patients/predictions-ncv.qza \
  --o-feature-importance Differential_abundance/sample-classifier/regressor_sklearn.collapsed-7.freq10-minSamp5.patients/importance-ncv.qza
qiime sample-classifier scatterplot \
  --i-predictions Differential_abundance/sample-classifier/regressor_sklearn.collapsed-7.freq10-minSamp5.patients/predictions-ncv.qza \
  --m-truth-file Bangladesh_metadata.tsv \
  --m-truth-column Month_num \
  --o-visualization Differential_abundance/sample-classifier/regressor_sklearn.collapsed-7.freq10-minSamp5.patients/scatter.qzv


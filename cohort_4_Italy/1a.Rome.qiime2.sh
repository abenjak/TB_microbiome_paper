# Install QIIME2
# conda create -n qiime2-2019.7 qiime=2019.7.0
conda activate qiime2-2019.7

# https://docs.qiime2.org/2018.6/tutorials/overview/
# https://chmi-sops.github.io/mydoc_qiime2.html

# Note: All QIIME 2 visualizers will generate a .qzv file. You can view these files with qiime tools view at https://view.qiime2.org/

# Prepare the metadata and tabulate it (to veryfy that all is fine)
# Sometimes qiime complains the file has a wrong encoding. In that case run this: iconv -t UTF-8 -f ISO-8859-1 in.tsv > out.tsv
qiime metadata tabulate \
 --m-input-file Rome_metadata.v3.tsv \
 --o-visualization Rome_metadata.v3.qzv

# These files are too large for GitHub:
# qiime tools import \
 # --type SampleData[JoinedSequencesWithQuality] \
 # --input-path Rome_fastq_manifest.csv \
 # --output-path Rome_reads.qza \
 # --input-format SingleEndFastqManifestPhred33

# qiime demux summarize \
 # --i-data Rome_reads.qza \
 # --o-visualization Rome_reads.qzv

# Since the reads are already merged, at this stage you can choose to proceed using Deblur for additional quality control, or you can dereplicate sequences and optionally cluster them into OTUs with q2-vsearch. If you try this option, we strongly encourage you to call qiime quality-filter q-score-joined with a higher min-quality threshold - possibly \
 --p-min-quality 20 or  --p-min-quality 30.
# qiime quality-filter q-score-joined \
 # --p-min-quality 30 \
 # --i-demux Rome_reads.qza \
 # --o-filtered-sequences Rome.QF30.qza \
 # --o-filter-stats Rome.QF30.stats.qza


#######################################  deblur start #################################################

# You should pass the sequence length value you selected from the quality score plots for \
 --p-trim-length. This will trim all sequences to this length, and discard any sequences which are not at least this long. Use a trim length of based on the quality score plots
# Andrej: this plot is not helpful. Will get the reads length distribution differently:
# zcat *SeqPrep-merged.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > Rome.readlenght.txt

Rscript -e 'reads<-read.csv(file="Rome.readlength.txt", sep="", header=FALSE); jpeg("Rome.readlength.jpg" , width=1000, height=1000); plot (reads$V1,reads$V2,type="l",xlab="read length",ylab="occurences",col="blue"); dev.off()'

# 99% of reads are at least 333 nt long.

# Note that deblur (and also vsearch dereplicate-sequences) should be preceded by basic quality-score-based filtering, but this is unnecessary for dada2. Both Deblur and DADA2 contain internal chimera checking methods and abundance filtering, so additional filtering should not be necessary following these methods.
# try deblur on filtered data (qual 30 and 20):

# qiime deblur denoise-16S \
 # --i-demultiplexed-seqs Rome.QF30.qza \
 # --p-trim-length 333 \
 # --p-sample-stats \
 # --p-jobs-to-start 4 \
 # --o-representative-sequences Rome.QF30.deblur_rep-seqs.qza \
 # --o-table Rome.QF30.deblur_tab.qza \
 # --o-stats Rome.QF30.deblur_stats.qza
qiime feature-table summarize \
 --i-table Rome.QF30.deblur_tab.qza \
 --o-visualization Rome.QF30.deblur_tab.qzv \
 --m-sample-metadata-file Rome_metadata.v3.tsv 
qiime feature-table tabulate-seqs \
  --i-data Rome.QF30.deblur_rep-seqs.qza \
  --o-visualization Rome.QF30.deblur_rep-seqs.qzv

#Result: RP14-3 has fewest features (1,558), RP43-6 the most (8,183). Overall quite uniform.

#######################################  deblur end #################################################



#######################################  qiime feature-classifier start #######################################

qiime feature-classifier classify-sklearn \
 --p-n-jobs 2 \
 --i-classifier ../Greengenes_QIIME2/gg_13_8.99_otus.extract-from-primers_ref-seqs.fit-classifier-naive-bayes.qza \
 --i-reads Rome.QF30.deblur_rep-seqs.qza \
 --o-classification Rome.QF30.deblur_rep-seqs.sklearn.qza
qiime metadata tabulate \
 --m-input-file Rome.QF30.deblur_rep-seqs.sklearn.qza \
 --o-visualization Rome.QF30.deblur_rep-seqs.sklearn.qzv

#######################################  qiime feature-classifier end #######################################



#######################################  taxonomic analyses start #######################################  
# Taxa collapse (will need it downstream).
# https://docs.qiime2.org/2018.8/plugins/available/taxa/collapse/
# Collapse groups of features that have the same taxonomic assignment through the specified level. The frequencies of all features will be summed when they are collapsed.
# Optional. You can visualize all levels in qiime taxa barplot on the whole dataset.

qiime taxa collapse \
 --i-table Rome.QF30.deblur_tab.qza\
 --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 2 \
 --o-collapsed-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-2.qza

qiime taxa collapse \
 --i-table Rome.QF30.deblur_tab.qza\
 --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 3 \
 --o-collapsed-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.qza

qiime taxa collapse \
 --i-table Rome.QF30.deblur_tab.qza\
 --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 4 \
 --o-collapsed-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.qza

qiime taxa collapse \
 --i-table Rome.QF30.deblur_tab.qza\
 --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 5 \
 --o-collapsed-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.qza

qiime taxa collapse \
 --i-table Rome.QF30.deblur_tab.qza\
 --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 6 \
 --o-collapsed-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.qza

qiime taxa collapse \
 --i-table Rome.QF30.deblur_tab.qza\
 --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
 --p-level 7 \
 --o-collapsed-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.qza

# barplot (https://docs.qiime2.org/2018.8/plugins/available/taxa/barplot/)
qiime taxa barplot \
 --i-table Rome.QF30.deblur_tab.qza \
 --i-taxonomy  Rome.QF30.deblur_rep-seqs.sklearn.qza \
 --m-metadata-file Rome_metadata.v3.tsv \
 --o-visualization Rome.QF30.deblur_rep-seqs.sklearn.taxa-barplot.qzv

# Heatmap (https://docs.qiime2.org/2018.8/plugins/available/feature-table/heatmap/?highlight=taxa%20collapse)
qiime feature-table heatmap \
 --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-2.qza \
 --m-metadata-file Rome_metadata.v3.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Rome.QF30.deblur_tab.sklearn.collapsed-lev-2.heatmap.qzv

qiime feature-table heatmap \
 --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.qza \
 --m-metadata-file Rome_metadata.v3.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.heatmap.qzv

qiime feature-table heatmap \
 --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.qza \
 --m-metadata-file Rome_metadata.v3.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.heatmap.qzv

qiime feature-table heatmap \
 --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.qza \
 --m-metadata-file Rome_metadata.v3.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.heatmap.qzv

qiime feature-table heatmap \
 --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.qza \
 --m-metadata-file Rome_metadata.v3.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.heatmap.qzv

qiime feature-table heatmap \
 --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.qza \
 --m-metadata-file Rome_metadata.v3.tsv \
 --m-metadata-column GroupMonth \
 --o-visualization Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.heatmap.qzv


#######################################  Sequence alignment and phylogeny building  start #######################################
# If you are sequencing phylogenetic markers (e.g., 16S rRNA genes), you can align these sequences to assess the phylogenetic relationship between each of your features. This phylogeny can then be used by other downstream analyses, such as UniFrac distance analyses.
# https://forum.qiime2.org/t/q2-phylogeny-community-tutorial/4455

# Sequence Alignment
qiime alignment mafft \
 --i-sequences Rome.QF30.deblur_rep-seqs.qza \
 --o-alignment Rome.QF30.deblur_rep-seqs.aligned.qza

qiime alignment mask \
 --i-alignment Rome.QF30.deblur_rep-seqs.aligned.qza \
 --o-masked-alignment Rome.QF30.deblur_rep-seqs.aligned.masked.qza


# fragment insertion analysis (https://github.com/biocore/q2-fragment-insertion)
# QIIME 2's "Moving Pictures" tutorial suggests constructing a de-novo phylogeny for the fragments, i.e FeatureData[Sequence], to obtain a Phylogeny[Rooted] that can be used for phylogenetic diversity computation. "Fragment insertion" provides an alternative way to acquire the Phylogeny[Rooted] by inserting sequences of FeatureData[Sequence] into a high quality reference phylogeny and thus provides multiple advantages over de-novo phylogenies, e.g. accurate branch lengths, multi-study meta-analyses, mixed region meta-analyses (e.g. V4 and V2).
# Fragment insertion avoids artificially long outgroup branches that would lead to exaggerated separation in beta diversity.
# Fragment insertion enables meta-analyses across different variable 16S regions and fragment length.
# Default reference (phylogeny and matching alignment) is Greengenes 13_8 at 99%.

qiime fragment-insertion sepp \
  --p-threads 4 \
  --i-representative-sequences Rome.QF30.deblur_rep-seqs.qza \
  --o-tree Rome.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --o-placements Rome.QF30.deblur_rep-seqs.fragment-insertion-placements.qza

# After you created the insertion tree, which can be used for phylogenetic diversity computation, e.g. Faith's PD or UniFrac, a typical next step is to filter a feature-table such that it only contains fragments that are in the insertion tree. This becomes necessary, since SEPP might reject insertion of fragments that are too remotely related to everything in the reference alignment/phylogeny. Those rows in your feature-table for fragments not in the phylogeny, will cause diversity computation to fail, since branch lengths cannot be determined.
# In my case, nothing got removed. Anyway, this step may be useful for troubleshooting.

qiime fragment-insertion filter-features \
  --i-table Rome.QF30.deblur_tab.qza \
  --i-tree Rome.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --o-filtered-table Rome.QF30.deblur_tab.fragment-insertion_filtered_table.qza \
  --o-removed-table Rome.QF30.deblur_tab.fragment-insertion_removed_table.qza

qiime feature-table summarize \
 --i-table Rome.QF30.deblur_tab.fragment-insertion_filtered_table.qza \
 --o-visualization Rome.QF30.deblur_tab.fragment-insertion_filtered_table.qzv

qiime feature-table summarize \
 --i-table Rome.QF30.deblur_tab.fragment-insertion_removed_table.qza \
 --o-visualization Rome.QF30.deblur_tab.fragment-insertion_removed_table.qzv


#######################################  Sequence alignment and phylogeny building  end #######################################

##### Alpha and beta diversity analysis

# Alpha rarefaction plotting
# Note: The value that you provide for --p-max-depth should be determined by reviewing the “Frequency per sample” information presented in the table.qzv file that was created above. In general, choosing a value that is somewhere around the median frequency seems to work well, but you may want to increase that value if the lines in the resulting rarefaction plot don’t appear to be leveling out, or decrease that value if you seem to be losing many of your samples due to low total frequencies closer to the minimum sampling depth than the maximum sampling depth.
# In my case, the maximum sample total frequency of the feature_table is 8183
qiime diversity alpha-rarefaction \
  --i-table Rome.QF30.deblur_tab.fragment-insertion_filtered_table.qza \
  --i-phylogeny Rome.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --p-max-depth 7000 \
  --p-steps 30 \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Rome.QF30.deblur_tab.fragment-insertion_filtered_table.alpha-rarefaction.qzv

# Results:
# Based on the rarefaction curves, I can use all the samples, since the weakest sample RP14-3 (1,558 features) is passed the rarefaction plateau.

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny Rome.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --i-table Rome.QF30.deblur_tab.fragment-insertion_filtered_table.qza \
  --p-sampling-depth 1558 \
  --m-metadata-file Rome_metadata.v3.tsv \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv


# Check if continuous sample metadata columns are correlated with alpha diversity:
qiime diversity alpha-correlation \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization core-metrics-results/faith-pd-correlation-spearman.qzv

qiime diversity alpha-correlation \
  --p-method pearson \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization core-metrics-results/faith-pd-correlation-pearson.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization core-metrics-results/shannon-correlation-spearman.qzv

qiime diversity alpha-correlation \
  --p-method pearson \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization core-metrics-results/shannon-correlation-pearson.qzv


# The following commands will test whether distances between samples within a group, such as samples from the same body site, are more similar to each other then they are to samples from the other groups.
# --p-pairwise parameter will also perform pairwise tests that will allow you to determine which specific pairs of groups differ from one another, if any.
# This command can be slow to run, so we’ll run this on specific columns of metadata that we’re interested in exploring, rather than all metadata columns that it’s applicable to.
# Cannot use the GroupMonth because some groups are represented by only 1 sample
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/unweighted-unifrac_Group_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization core-metrics-results/unweighted-unifrac_GroupMonthSimple_significance.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/bray_curtis_Group_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization core-metrics-results/bray_curtis_GroupMonthSimple_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/jaccard_Group_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization core-metrics-results/jaccard_GroupMonthSimple_significance.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/weighted-unifrac_Group_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization core-metrics-results/weighted-unifrac_GroupMonthSimple_significance.qzv \
  --p-pairwise

# While our core-metrics-phylogenetic command did already generate some Emperor plots, we want to pass an optional parameter, --p-custom-axes, which is very useful for exploring time series data
# We will generate Emperor plots for unweighted UniFrac and Bray-Curtis so that the resulting plot will contain axes for principal coordinate 1, principal coordinate 2, and days since the experiment start. We will use that last axis to explore how these samples changed over time.
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --p-custom-axes Month_num \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor_Month.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --p-custom-axes Month_num \
  --o-visualization core-metrics-results/weighted-unifrac-emperor_Month.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
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
  --i-table Rome.QF30.deblur_tab.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Rome.QF30.deblur_tab.freq10-minSamp5.qza

mkdir Differential_abundance/ANCOM

qiime composition add-pseudocount \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.qza

qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv

qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.ancom-GroupMonthSimple.qzv

qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Group \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.freq10-minSamp5.added-pseudocount.ancom-Group.qzv

# We’re also often interested in performing a differential abundance test at a specific taxonomic level. 

qiime feature-table filter-features \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-3.freq10-minSamp5.added-pseudocount.ancom-GroupMonthSimple.qzv


qiime feature-table filter-features \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.added-pseudocount.ancom-GroupMonthSimple.qzv


qiime feature-table filter-features \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.added-pseudocount.ancom-GroupMonthSimple.qzv


qiime feature-table filter-features \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.added-pseudocount.ancom-GroupMonthSimple.qzv


qiime feature-table filter-features \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.qza \
  --p-min-frequency 10 \
  --p-min-samples 5 \
  --o-filtered-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.qza
qiime composition add-pseudocount \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.qza \
  --o-composition-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.added-pseudocount.qza
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonth \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.added-pseudocount.ancom-GroupMonth.qzv
qiime composition ancom \
  --i-table Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.added-pseudocount.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/ANCOM/Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.added-pseudocount.ancom-GroupMonthSimple.qzv

# gneiss (https://docs.qiime2.org/2018.8/tutorials/gneiss/)
# make sure to use the filtered table else the results will be biased, or you'll run into errors (this happened on deblur unfiltered table: https://forum.qiime2.org/t/gneiss-zero-balance-error/1857/5)

# gneiss Correlation-clustering
mkdir Differential_abundance/gneiss

qiime gneiss correlation-clustering \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --o-clustering Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza

qiime gneiss ilr-hierarchical \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --o-balances Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza

qiime gneiss ols-regression \
  --p-formula "Subject+Group+GroupMonth" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Subject+Group+GroupMonth.qzv

qiime gneiss dendrogram-heatmap \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_heatmap_GroupMonthSimple.qzv

qiime gneiss dendrogram-heatmap \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Subject \
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_heatmap_Subject.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 2 \
  --p-balance-name 'y0' \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y0_taxa-level-2_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 6 \
  --p-balance-name 'y0' \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y0_taxa-level-6_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y0' \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y0_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y1' \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y1_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y2' \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_y2_taxa-level-5_summary.qzv


# Try different formulas in the ols-regression
qiime gneiss ols-regression \
  --p-formula "GroupMonth" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_GroupMonth.qzv

qiime gneiss ols-regression \
  --p-formula "GroupMonthSimple" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_GroupMonthSimple.qzv

qiime gneiss ols-regression \
  --p-formula "Group*Month_cat" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Group-x-Month.qzv

qiime gneiss ols-regression \
  --p-formula "Group+Month_cat" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Group+Month.qzv

qiime gneiss ols-regression \
  --p-formula "Subject" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Subject.qzv

qiime gneiss ols-regression \
  --p-formula "Subject+Month_cat" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Subject+Month.qzv

qiime gneiss ols-regression \
  --p-formula "Group" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Group.qzv

qiime gneiss ols-regression \
  --p-formula "Month_cat" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_cor-clust_regression_summary_Month.qzv


# Gneiss Gradient-clustering
qiime gneiss gradient-clustering \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --m-gradient-file Rome_metadata.v3.tsv \
  --m-gradient-column Month_num \
  --o-clustering Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza

qiime gneiss ilr-hierarchical \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --o-balances Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_balances.qza

qiime gneiss ols-regression \
  --p-formula "Subject+Group+Month_cat" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_regression_summary_Subject+Group+Month.qzv

qiime gneiss ols-regression \
  --p-formula "Subject" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_regression_summary_Subject.qzv

qiime gneiss ols-regression \
  --p-formula "Group+Month_cat" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_regression_summary_Group+Month.qzv

qiime gneiss ols-regression \
  --p-formula "Group" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_regression_summary_Group.qzv

qiime gneiss ols-regression \
  --p-formula "Month_cat" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_regression_summary_Month.qzv


qiime gneiss dendrogram-heatmap \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Month_cat \
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_heatmap_Month.qzv

qiime gneiss dendrogram-heatmap \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Group\
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_heatmap_Group.qzv

qiime gneiss dendrogram-heatmap \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple\
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_heatmap_GroupMonthSimple.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 2 \
  --p-balance-name 'y6' \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column GroupMonthSimple \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_gradient_y6_taxa-level-2_summary.qzv

# Gneiss Gradient-clustering (only with TB patients)

qiime feature-table filter-samples \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --p-where "[Group]='TB'" \
  --o-filtered-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza

Gneiss seems to work better if only TB metadata is given.
qiime gneiss gradient-clustering \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --m-gradient-file Rome_metadata_TB.tsv \
  --m-gradient-column Month_num \
  --o-clustering Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza

qiime gneiss ilr-hierarchical \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --o-balances Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_balances.qza

qiime gneiss ols-regression \
  --p-formula "Subject+Month_cat" \
  --i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_balances.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata_TB.tsv \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_regression_summary_Subject+Month.qzv

qiime gneiss dendrogram-heatmap \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --p-color-map seismic \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_heatmap_Month.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 2 \
  --p-balance-name 'y8' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y8_taxa-level-2_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 4 \
  --p-balance-name 'y8' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y8_taxa-level-4_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y8' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y8_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 6 \
  --p-balance-name 'y8' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y8_taxa-level-6_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 2 \
  --p-balance-name 'y5' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y5_taxa-level-2_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 4 \
  --p-balance-name 'y5' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y5_taxa-level-4_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y5' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y5_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 6 \
  --p-balance-name 'y5' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y5_taxa-level-6_summary.qzv


qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 4 \
  --p-balance-name 'y0' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y0_taxa-level-4_summary.qzv


qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y0' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y0_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 6 \
  --p-balance-name 'y0' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y0_taxa-level-6_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y1' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y1_taxa-level-5_summary.qzv

qiime gneiss balance-taxonomy \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.TB.qza \
  --i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_hierarchy.qza \
  --i-taxonomy Rome.QF30.deblur_rep-seqs.sklearn.qza \
  --p-taxa-level 5 \
  --p-balance-name 'y1' \
  --m-metadata-file Rome_metadata_TB.tsv \
  --m-metadata-column Month_cat \
  --o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.TB.gneiss_gradient_y1_taxa-level-5_summary.qzv


# gneiss Phylogenetic analysis
# Results are similar as with correlation-clustering, but some steps don't work. No need to pursue it further.
qiime gneiss ilr-phylogenetic \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --i-tree Rome.QF30.deblur_rep-seqs.fragment-insertion-tree.qza \
  --o-balances Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_balances.qza \
  --o-hierarchy Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_hierarchy.qza

## I'm getting this error here: 'the label [y105] is not in the [index]'
#qiime gneiss ols-regression \
  #--p-formula "Subject+Group+GroupMonth" \
  #--i-table Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_balances.qza \
  #--i-tree Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_hierarchy.qza \
  #--m-metadata-file Rome_metadata.v3.tsv \
  #--o-visualization Differential_abundance/gneiss/Rome.QF30.deblur_tab.freq10-minSamp5.gneiss_phylogeny_regression_summary_Subject+Group+GroupMonth.qzv

##################################
# Longitudinal analysis (https://docs.qiime2.org/2018.8/tutorials/longitudinal/)
mkdir Differential_abundance/longitudinal

# Linear mixed effects (LME) models test the relationship between a single response variable and one or more independent variables, where observations are made across dependent samples, e.g., in repeated-measures sampling experiments. 
qiime longitudinal linear-mixed-effects \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-file core-metrics-results/shannon_vector.qza \
  --p-metric shannon \
  --p-group-columns Group \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --o-visualization Differential_abundance/longitudinal/linear-mixed-effects.qzv

# Volatility analysis
# The volatility visualizer generates interactive line plots that allow us to assess how volatile a dependent variable is over a continuous, independent variable (e.g., time) in one or more groups.
qiime longitudinal volatility \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-file core-metrics-results/shannon_vector.qza \
  --p-default-metric shannon \
  --m-metadata-file core-metrics-results/evenness_vector.qza \
  --m-metadata-file core-metrics-results/faith_pd_vector.qza \
  --p-default-group-column Group \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --o-visualization Differential_abundance/longitudinal/volatility.qzv

# first-distances identifies the beta diversity distances between successive samples from the same subject.
qiime longitudinal first-distances \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --p-replicate-handling random \
  --o-first-distances Differential_abundance/longitudinal/first-distances.qza

# the output of first-distances is particularly empowering, though, because it allows us to analyze longitudinal changes in beta diversity using actions that cannot operate directly on a distance matrix, such as linear-mixed-effects
qiime longitudinal linear-mixed-effects \
  --m-metadata-file Differential_abundance/longitudinal/first-distances.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --p-metric Distance \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --p-group-columns Group \
  --o-visualization Differential_abundance/longitudinal/first-distances-LME.qzv

# Feature volatility analysis
# Doesnt work on the entire dataset because of "metadata column that contains one or more values that match only one sample"
# Try it on a subset of TB patients (no RP6-2, RP25-5, RP43-4 and RP43-5, based on a new metadata file)
qiime longitudinal feature-volatility \
  --i-table Rome.QF30.deblur_tab.freq10-minSamp5.qza \
  --m-metadata-file Rome_metadata_TB_no-singles.tsv \
  --p-missing-samples ignore \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_deblur_tab.freq10-minSamp5

qiime longitudinal feature-volatility \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-4.freq10-minSamp5.qza \
  --m-metadata-file Rome_metadata_TB_no-singles.tsv \
  --p-missing-samples ignore \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_collapsed-lev-4.freq10-minSamp5

qiime longitudinal feature-volatility \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-5.freq10-minSamp5.qza \
  --m-metadata-file Rome_metadata_TB_no-singles.tsv \
  --p-missing-samples ignore \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_collapsed-lev-5.freq10-minSamp5

qiime longitudinal feature-volatility \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.qza \
  --m-metadata-file Rome_metadata_TB_no-singles.tsv \
  --p-missing-samples ignore \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_collapsed-lev-6.freq10-minSamp5

qiime longitudinal feature-volatility \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-7.freq10-minSamp5.qza \
  --m-metadata-file Rome_metadata_TB_no-singles.tsv \
  --p-missing-samples ignore \
  --p-state-column Month_num \
  --p-individual-id-column Subject \
  --output-dir Differential_abundance/longitudinal/feat-volatility_collapsed-lev-7.freq10-minSamp5


# Non-parametric microbial interdependence test (NMIT)
# run it on the filtered table, but remove controls from it:
# NMIT needs relative frequency table as input. It's recommended to use a collapsed table.
qiime feature-table relative-frequency \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.qza \
  --o-relative-frequency-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.qza

qiime longitudinal nmit \
  --i-table Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --p-individual-id-column Subject \
  --p-corr-method pearson \
  --o-distance-matrix Differential_abundance/longitudinal/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.nmit-dm.qza

## perform PERMANOVA tests to evaluate whether between-group distances are larger than within-group distance.
qiime diversity beta-group-significance \
  --i-distance-matrix Differential_abundance/longitudinal/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.nmit-dm.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-column Group \
  --o-visualization Differential_abundance/longitudinal/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.nmit-dm.beta-group-significance.qzv

# Finally, we can compute principal coordinates and use Emperor to visualize similarities among subjects:
qiime diversity pcoa \
  --i-distance-matrix Differential_abundance/longitudinal/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.nmit-dm.qza \
  --o-pcoa Differential_abundance/longitudinal/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.nmit-pc.qza

qiime emperor plot \
  --i-pcoa Differential_abundance/longitudinal/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.nmit-pc.qza \
  --m-metadata-file Rome_metadata.v3.tsv \
  --o-visualization Differential_abundance/longitudinal/Rome.QF30.deblur_tab.sklearn.collapsed-lev-6.freq10-minSamp5.rel-freq.nmit-emperor.qzv


# Pairwise difference comparisons
# The nonTB timepoints mostly don't match with TB
qiime longitudinal pairwise-differences \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-file core-metrics-results/shannon_vector.qza \
  --p-metric shannon \
  --p-group-column Group \
  --p-state-column MonthSimple_num \
  --p-state-1 0 \
  --p-state-2 2.5 \
  --p-individual-id-column Subject \
  --o-visualization core-metrics-results/pairwise-differences.shannon.time_0-2.5.qzv

qiime longitudinal pairwise-differences \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-file core-metrics-results/shannon_vector.qza \
  --p-metric shannon \
  --p-group-column Group \
  --p-state-column MonthSimple_num \
  --p-state-1 0 \
  --p-state-2 5 \
  --p-individual-id-column Subject \
  --o-visualization core-metrics-results/pairwise-differences.shannon.time_0-5.qzv

qiime longitudinal pairwise-differences \
  --m-metadata-file Rome_metadata.v3.tsv \
  --m-metadata-file core-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Group \
  --p-state-column MonthSimple_num \
  --p-state-1 0 \
  --p-state-2 5 \
  --p-individual-id-column Subject \
  --o-visualization core-metrics-results/pairwise-differences.faith_pd.time_0-5.qzv

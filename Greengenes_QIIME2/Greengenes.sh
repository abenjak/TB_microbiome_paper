############# train the classifier
# gg_13_8.99_otus should have been installed with your QIIME2. If not, you can fetch if from Greengenes.
# Adjust your full paths accordingly:

qiime tools import \
 --type 'FeatureData[Sequence]' \
 --input-path 'Greengenes/gg_13_8_otus/rep_set/99_otus.fasta' \
 --output-path  gg_13_8.99_otus.qza
 
qiime tools import \
 --type 'FeatureData[Taxonomy]' \
 --source-format HeaderlessTSVTaxonomyFormat \
 --input-path Greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
 --output-path gg_13_8.99_ref-taxonomy.qza

# If you want to use Naive Bayes classifier, it has been shown that taxonomic classification accuracy improves when a Naive Bayes classifier is trained on only the region of the target sequences that was sequenced
qiime feature-classifier extract-reads \
 --i-sequences gg_13_8.99_otus.qza \
 --p-f-primer AGRGTTYGATYMTGGCTCAG \
 --p-r-primer TGCTGCCTCCCGTAGGAGT \
 --o-reads gg_13_8.99_otus.extract-from-primers_ref-seqs.qza

# We can now train a Naive Bayes classifier as follows, using the reference reads and taxonomy that we just created.
qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads gg_13_8.99_otus.extract-from-primers_ref-seqs.qza \
 --i-reference-taxonomy gg_13_8.99_ref-taxonomy.qza \
 --o-classifier gg_13_8.99_otus.extract-from-primers_ref-seqs.fit-classifier-naive-bayes.qza


# nf-ncov-voc

# Genomic Analysis

# Functional Annotation

In this module, the genomic data for each lineage is converted to a GVF (Genomic Variant Format) file and annotated with functional information.  GVF is a variant of GFF3 that is standardized for describing genomic mutations; it is used here because it can describe mutations across multiple rows, and because the "#attributes" column can store information in custom key-value pairs.  The key-value pairs added at this stage include for each mutation: VOC/VOI status, clade-defining status (for reference lineages), and functional annotations parsed from Paul Gordon's Pokay repository.  The annotated GVFs are the input to the visualization module.

# Surveillance Report

A TSV file can be produced that summarizes mutation information for SARS-CoV-2 variants.  This file contains the same information found in the GVF files in a more human-readable format, with all reference lineages per variant merged into one concise TSV.

# Visualization



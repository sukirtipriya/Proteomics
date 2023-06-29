# Proteomics

DEqMS is a statistical method for identification of differentially expressed proteins in MS-data, that takes into acount the dependence of variance on the number of PSMs or peptides used for protein quantification.

Mass spectrometry (MS)-based proteomics is widely used for identification and quantification of proteins from complex biological samples. In a typical proteomics experiment, proteins are digested into peptides using a proteolytic enzyme, commonly trypsin, prior to MS-analysis. Peptides are usually fragmented multiple times, thus generating multiple peptide spectrum matches (PSMs) for the same peptide sequence. Peptide identification and quantification is subsequently performed in data analysis workflows using either label-free or labeled approaches. Both quantitative approaches rely on a hierarchical data structure: PSMs are nested into peptides which are then nested into proteins. Protein level quantification is subsequently generated through summarizing peptide or PSM level information.


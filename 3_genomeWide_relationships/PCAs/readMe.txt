This folder contains the analysis, plotting code used to perform the PCAs shown in Fig. S6. Here we describe each file.

PCA_analysis_code.txt:
This file contains all of the code needed to convert the genetic data to an appropriate format, prune to reduce linkage disequilibrium of associated sites, subset the individuals, and perform the PCA.

keep_samples:
This file contains the IDs of individuals to be retained in the second PCAs that only included Northern saxatilis and arcana. This file is referred to in the file PCA_analysis _code.txt.

eigenvalues.xlsx:
This file contains the eigenvalues for the four separate PCAs shown in Fig. S6. The raw eigenvalues and percent variation of each PC are given.

pca_out:
This file contains the raw eigenvectors from the four PCAs. Columns D-M are for the full PCA, N-W are for the full-LD pruned PCA, AA-AJ are for the downsampled PCA, and AK to AT are for the downsampled-LD pruned PCA. The ‘mode’ and ‘species’ columsn give the IDs for the Full PCA. Species_down and mode_down are the same but for the downsampled dataset.

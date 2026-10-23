# Plant CpG Island (CGI) Comparative Analysis

A comparative bioinformatics workflow evaluating upstream CpG island (CGI) sequence characteristics and $k$-mer composition ($k=3$) across three *Brassicaceae* species (*A. thaliana*, *B. rapa*, *B. juncea*). 

This pilot study compares upstream CpG island (CGI) sequence features across three Brassicaceae species to investigate whether these features exhibit species-associated patterns or similar distributions across diploid and polyploid lineages.

## Environment

* Language: Python 3.12
* Libraries: Pandas, NumPy, Biopython, pybedtools, Scikit-learn, SciPy, kneed
* Visualization: Plotly, Matplotlib, Seaborn

## Getting Started (How to Run)
1. Install dependencies: `pip install -r requirements.txt`
2. Install `bedtools` CLI (required by `pybedtools`): e.g. `conda install -c bioconda bedtools`
3. Download raw genome files and place them under `dataset/raw/` (see `dataset/raw/README.md` for sources and exact filenames)
4. Run notebooks in numbered order: `notebook/01_downloader.ipynb` → `notebook/02_preprocess.ipynb` → `notebook/03_feature_eda.ipynb` → `notebook/04_clustering.ipynb` → `notebook/05_decision_tree.ipynb` → `notebook/06_feature_analysis.ipynb`

## Data

|Species|Data Type|Source|Version|Format|
|-|-|-|-|-|
|*Arabidopsis thaliana*|Upstream 1000bp|[TAIR](https://www.arabidopsis.org)|TAIR10, upstream file dated 2010-11-04|FASTA, GFF|
|*Brassica juncea* (var. tumida)|Chromosome-level|[BRAD](http://brassicadb.cn/) (inferred from file naming)|Braju_tum_V2.0|FASTA, GFF|
|*Brassica rapa* (cv. Chiifu)|Chromosome-level|[BRAD](http://brassicadb.cn/)|Brara_Chiifu_V4.0|FASTA, GFF|

*B. juncea* is an amphidiploid (AB genome, formed via hybridization of *B. rapa*-like and *B. nigra*-like ancestors), while *A. thaliana* and *B. rapa* are diploid — a structural difference relevant to interpreting species-level differences observed later.

## Workflow

### 1. `01_downloader.ipynb` — Upstream Sequence Extraction

Extracts upstream sequences for three *Brassicaceae* species. For *B. rapa* and *B. juncea*, the upstream 1,000-bp regions are extracted using mRNA coordinates from GFF annotations; *A. thaliana* uses pre-extracted TAIR10 upstream sequences.

* *B. rapa*, *B. juncea*: computed directly from chromosome FASTA + GFF mRNA coordinates (`+` strand: `[start-1000, start)`; `-` strand: `(end, end+1000]`, reverse-complemented).
* *A. thaliana*: uses TAIR10's pre-extracted upstream FASTA file.
	- **Note**: TAIR10 uses pre-extracted upstream sequences, whereas *Brassica* regions were extracted from mRNA coordinates in GFF files. Differences in upstream-region definitions may affect cross-species comparisons (see Limitations).
* For genes with multiple isoforms, only the first-encountered transcript per base gene ID is kept (a simplification from an earlier draft that selected the longest isoform).

### 2. `02_preprocess.ipynb` — CpG Island Detection

Applies a 200-bp sliding window to each upstream sequence and incrementally updates C, G, and CpG counts as the window moves.

A window qualifies as CGI if:

1. Region length > 200bp (final island length)
2. GC content > 50%
3. Observed/Expected CpG ratio > 0.6

(criteria per Gardiner-Garden & Frommer, 1987).

Consecutive qualifying windows are merged into a single island, tolerating gaps of up to 100bp.

### 3. `03_feature_eda.ipynb` — Exploratory Analysis

* Merges all three species' CGI datasets, tagging each row with species of origin.
* Visualizes GC content / obs-exp ratio / length in 3D to check visual separability.
* Extracts 3-mer (trinucleotide) composition per island, normalized by sequence length. **k=3 was an exploratory choice, not motivated by a specific biological hypothesis** (e.g. codon structure); other k values were not compared.
* Reduces 3-mer features to 11 principal components via PCA (elbow method: no clear knee point found; n=11 selected as the smallest n reaching ≥70% cumulative explained variance).
**Note:** Visualizations in this notebook may not render correctly in GitHub's preview. For the complete notebook outputs, view the notebook in [nbviewer](https://nbviewer.org/github/imbesky/CpG-island-finding/blob/main/notebook/03_feature_eda.ipynb)

### 4. `04_clustering.ipynb` — DBSCAN

* Standardized CGI physical features (GC content, length, obs/exp) before clustering, since DBSCAN relies on Euclidean distance and unscaled 'length' (range ~200-1000+) would dominate over the other two features (range ~0-2).
* `eps` selected via k-distance elbow method;
  * `min_samples` candidates were multiples of feature dimensionality for the CGI feature set (2n, 3n, 4n, 5n where n=3).
  * The 3-mer feature set's min_samples candidates ([10, 18, 27]) do not follow this same clean rule and the original derivation was not preserved — a reminder of why logging parameter choices during exploratory work matters.
* Best (eps, min_samples) selected by silhouette score.

### 5. `05_decision_tree.ipynb` — Species Classification

Trains a decision tree (GridSearchCV to find optimal parameters, 80/20 train/test split) to predict species from 
- (a) CGI physical features
- (b) raw 3-mer composition
Examines which features contribute to the models' species classification predictions.

- **Note**: Scoring is standard accuracy.
- **Note**: The class distribution is imbalanced, and the current analysis does not fully disentangle the effects of class imbalance from genuine biological separability.

### 6. `06_feature_analysis.ipynb` — Statistical Validation

* **Radar chart**: per-species median of each CGI feature, z-score standardized for comparison on one chart. Standardization here is computed across only 3 data points (one median per species) — best read qualitatively, not as a rigorous statistical test.
* **ANOVA + box plot**: one-way ANOVA per feature (much larger sample size — every individual CGI, not species medians) to test whether mean CGI feature values differ among species.
* **KDE plot**: per-species density distributions, independently normalized (`common_norm=False`) so shape is comparable regardless of sample size differences between species.
* **3-mer frequency profile**: compares average proportions of all 64 trinucleotides across species, with an additional focus on GC-rich trinucleotides explored in the decision-tree feature-importance analysis.

## Key Findings & Interpretation

1. Limited Species Separation in CGI Physical Features — or Insufficient Discriminative Signal
	DBSCAN unsupervised clustering placed the vast majority of CGI physical features (GC content, length, and observed/expected CpG ratio) in a single cluster, with some observations classified as noise, across all three species.
This may reflect overlap in the measured feature distributions across the three species, limited discriminative power of the selected features, or sensitivity to the chosen clustering parameters.
The current analysis cannot distinguish among these explanations or establish evolutionary conservation from the clustering results alone.

2. Species-Associated Differences in *B. juncea*
	One-way **ANOVA** detected statistically significant differences among the three species in all three physical features: CGI length, GC content, and observed/expected CpG ratio ($p < 0.001$ for each feature).
Non-parametric **KDE plots** also showed differences in CGI length and GC content distributions, particularly for *B. juncea*.
These findings describe species-associated patterns in the analyzed datasets.

3. Classification Performance and Feature Importance
	Supervised decision tree models achieved approximately 62% accuracy in species classification across the observed runs, compared with 33.3% accuracy expected under uniform random guessing.
Because the class distribution was imbalanced, this comparison alone does not establish strong classification performance.
Feature-importance analysis highlighted **GC content** among the physical features and **GC-rich trinucleotides** among the k-mer features as influential predictors in the observed models.
However, the ranking of individual trinucleotides varied between runs, so no single 3-mer can be considered a stable, dominant species discriminator based on the current results.

## Methodological Evaluation & Limitations & Future Research Proposal

This exploratory pipeline provided an initial framework for comparing upstream CGI features across the three species. Its limitations identify several priorities for further validation and refinement.

### 1. Biological and Data-Definition Limitations & Proposed Solutions
- **Upstream Region Definition Differences:** *A. thaliana* uses pre-extracted TAIR10 upstream sequences, whereas *Brassica* regions were extracted from mRNA coordinates using GFF files.
	- *Next Step:* Re-extract upstream regions using a consistent transcript-centric rule (e.g., representative transcript TSS) to reduce potential definition bias.
- **Window Boundary Effects:** Some of the *B. juncea* CGI lengths cluster near the 1000bp extraction boundary, suggesting some islands may extend beyond this window and remain undetected.
- **Species-Level Aggregation vs. Functional Stratification:** Aggregating all analyzed upstream CGIs may obscure patterns associated with particular gene functions.
	- *Hypothesis:* Housekeeping genes may show more conserved CGI features, whereas **stress-response genes may exhibit greater divergence**. Stratifying genes by Gene Ontology (GO) terms could help test whether CGI patterns differ across functional categories and species.

### 2. Methodological & Machine Learning Refinements
- **Feature Representation ($k$-mer Spectrum):** $k=3$ was used as an exploratory representation but alternative values were not compared. Higher-order k-mers ($k=4..6$) or position-specific representations may capture additional sequence patterns, although their utility remains to be tested.
- **Small-Sample Standardization:** The radar chart's z-scores were computed across only 3 data points (species medians), limiting it to a qualitative visualization rather than a statistically rigorous comparison.
- **Class Imbalance & Metrics:** Accuracy (~62%) may be influenced by class imbalance. Future evaluation should compare accuracy with a majority-class baseline and include balanced accuracy. Class weighting or resampling could also be evaluated, with resampling restricted to training folds to avoid data leakage.
- **Hyperparameter Auditability:** The reason for some parameter choices, such as the candidate min_samples values for 3-mer DBSCAN, was not fully documented. Recording parameter-selection procedures and software versions would improve reproducibility.

## Repository Structure

```
CpG_island_finding
├── notebook
|   ├── 01_downloader.ipynb        # 1. Upstream sequence extraction
|   ├── 02_preprocess.ipynb        # 2. CpG island detection
|   ├── 03_feature_eda.ipynb       # 3. Data merge, EDA, 3-mer extraction, PCA
|   ├── 04_clustering.ipynb        # 4. DBSCAN
|   ├── 05_decision_tree.ipynb     # 5. Species classification
|   └── 06_feature_analysis.ipynb  # 6. Statistical validation & visualization
└── dataset/
    ├── raw/                # Original FASTA/GFF files
    └── processed/          # Workflow intermediate/output CSVs
```

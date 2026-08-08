# Plant CpG Island (CGI) Comparative Analysis

A comparative bioinformatics pipeline evaluating upstream CpG island (CGI) sequence
characteristics and $k$-mer composition ($k=3$) across three *Brassicaceae* species
(*A. thaliana*, *B. rapa*, *B. juncea*). 

This pilot study serves as an empirical benchmark to test whether genome-wide upstream
CGI sequence features retain species-specific divergence signals or reflect conserved
regulatory architecture across diploid and polyploid lineages.

## Environment

* Language: Python 3.12
* Libraries: Pandas, NumPy, Biopython, pybedtools, Scikit-learn, SciPy, kneed
* Visualization: Plotly, Matplotlib, Seaborn

## Getting Started (How to Run)
1. Install dependencies: `pip install -r requirements.txt`
2. Install `bedtools` CLI (required by `pybedtools`): e.g. `conda install -c bioconda bedtools`
3. Download raw genome files and place them under `dataset/raw/`
   (see `dataset/raw/README.md` for sources and exact filenames)
4. Run notebooks in numbered order:
   `01_downloader.ipynb` → `02_preprocess.ipynb` → `03_feature_eda.ipynb`
   → `04_clustering.ipynb` → `05_decision_tree.ipynb` → `06_feature_analysis.ipynb`

## Data

|Species|Data Type|Source|Version|Format|
|-|-|-|-|-|
|*Arabidopsis thaliana*|Upstream 1000bp|[TAIR](https://www.arabidopsis.org)|TAIR10, upstream file dated 2010-11-04|FASTA, GFF|
|*Brassica juncea* (var. tumida)|Chromosome-level|[BRAD](http://brassicadb.cn/) (inferred from file naming)|Braju_tum_V2.0|FASTA, GFF|
|*Brassica rapa* (cv. Chiifu)|Chromosome-level|[BRAD](http://brassicadb.cn/)|Brara_Chiifu_V4.0|FASTA, GFF|

*B. juncea* is an amphidiploid (AB genome, formed via hybridization of
*B. rapa*-like and *B. nigra*-like ancestors), while *A. thaliana* and
*B. rapa* are diploid — a structural difference relevant to interpreting
species-level differences observed later.

## Pipeline

### 1. `01_downloader.ipynb` — Upstream Sequence Extraction

Extracts each gene's upstream 1,000bp region relative to its transcription
start site (TSS), since CGIs predominantly occur in promoter/5'UTR regions.

* *B. rapa*, *B. juncea*: computed directly from chromosome FASTA + GFF
mRNA coordinates (`+` strand: `[start-1000, start)`; `-` strand: `(end, end+1000]`, reverse-complemented).
* *A. thaliana*: uses TAIR10's pre-extracted upstream FASTA file.
**Note**: TAIR defines this upstream region relative to the *maximum extent across all gene models at a locus*, not a single mRNA's TSS — a definition mismatch from the *Brassica* extraction method (see Limitations).
* For genes with multiple isoforms, only the first-encountered transcript per base gene ID is kept (a simplification from an earlier draft that selected the longest isoform).

### 2. `02_preprocess.ipynb` — CpG Island Detection

Applies a 200bp sliding window over each upstream sequence, incrementally
updating C/G/CpG counts per step (O(1) per window vs. O(window size)).

A window qualifies as CGI if:

1. Region length > 200bp
2. GC content > 50%
3. Observed/Expected CpG ratio > 0.6

(criteria per Gardiner-Garden & Frommer, 1987). Consecutive qualifying
windows are merged into a single island, tolerating gaps of up to 100bp.

### 3. `03_feature_eda.ipynb` — Exploratory Analysis

* Merges all three species' CGI datasets, tagging each row with species of origin.
* Visualizes GC content / obs-exp ratio / length in 3D to check visual separability.
* Extracts 3-mer (trinucleotide) composition per island, normalized by sequence length. **k=3 was an exploratory choice, not motivated by a specific biological hypothesis** (e.g. codon structure); other k values were not compared.
* Reduces 3-mer features to 11 principal components via PCA (elbow method: no clear knee point found; n=11 selected as the smallest n reaching ≥70% cumulative explained variance).

### 4. `04_clustering.ipynb` — DBSCAN

* Standardizes CGI physical features (GC content, length, obs/exp) before clustering, since DBSCAN relies on Euclidean distance and unscaled  'length' (range ~200-1000+) would dominate over the other two features (range ~0-2).
* `eps` selected via k-distance elbow method;

  * `min_samples` candidates were multiples of feature dimensionality for the CGI feature set (2n, 3n, 4n, 5n where n=3).
  * The 3-mer feature set's min_samples candidates ([10, 18, 27]) do not follow this same clean rule and the original derivation was not preserved — a reminder of why logging parameter choices during exploratory work matters.
* Best (eps, min_samples) selected by silhouette score.

### 5. `05_decision_tree.ipynb` — Species Classification

Trains a decision tree (GridSearchCV, 80/20 train/test split) to predict
species from (a) CGI physical features and (b) raw 3-mer composition,
identifying which individual features drive whatever separability exists.

**Note**: scoring is standard accuracy. Given *B. juncea*'s amphidiploid
genome likely yields more genes (and thus more CGIs) than the two diploid
species, class imbalance may partly explain the model's tendency to
classify most samples as mustard — a factor not fully disentangled from
genuine biological separability in the current analysis.

### 6. `06_feature_analysis.ipynb` — Statistical Validation

* **Radar chart**: per-species median of each CGI feature, z-score standardized for comparison on one chart. Standardization here is computed across only 3 data points (one median per species) — best read qualitatively, not as a rigorous statistical test.
* **ANOVA + box plot**: one-way ANOVA per feature (much larger sample size — every individual CGI, not species medians) to test whether species-level means differ significantly.
* **KDE plot**: per-species density distributions, independently normalized (`common_norm=False`) so shape is comparable regardless of sample size differences between species.
* **3-mer frequency profile**: full 64-mer average proportion per species, and a focused comparison of GC-centric 3-mers (motivated by the decision tree's feature importance results, where 'GCG' was the dominant predictor).

## Key Findings & Interpretation

1. High Sequence Feature Conservation — or Insufficient Discriminative Signal
DBSCAN unsupervised clustering collapsed the vast majority of CGIs into a
single dense cluster across all three species. This could reflect strong
**evolutionary conservation** of promoter-proximal CGI architecture across
the *Brassicaceae* family — though it may equally indicate that the three
physical features used (GC content, length, obs/exp) lack sufficient
discriminative power to resolve species-level differences. The two
explanations cannot be disentangled with the current feature set alone.

2. Polyploid-Specific Feature Shift in *B. juncea*
One-way ANOVA revealed statistically significant overall differences ($p < 0.001$), driven primarily by *B. juncea* (amphidiploid). Non-parametric KDE plots confirmed a distinct distribution shift in CGI length and GC content for *B. juncea* compared to the two diploid species, highlighting potential structural genomic changes following hybridization.

3. Discriminative Capacity and Feature Dominance
Supervised decision tree models achieved ~62% accuracy (baseline random chance: 33.3%). Feature importance analysis revealed that **GC content** (physical features) and the **'GCG' trinucleotide** (accounting for 83% of feature importance in 3-mer models) served as primary discriminators, indicating that specific GC-rich micro-motifs carry the bulk of species-differentiating signals.

## Critical Methodological Evaluation & Future Research Proposal

While this pipeline established a robust data-processing baseline, critical evaluation
reveals key methodological bottlenecks that inform necessary next steps:

### 1. Biological Bottlenecks & Proposed Solutions
- **Locus Definition Artifacts:** TAIR10's upstream regions rely on locus-level maximum bounds, whereas *Brassica* extractions used transcript TSS.
	- *Next Step:* Re-extract all promoter regions under a single transcript-centric rule (e.g., primary/representative transcript TSS) to eliminate definition bias.
- **Window Boundary Effects:** *B. juncea* CGI lengths cluster near the 1000bp extraction boundary, suggesting some islands may extend beyond this arbitrary window and remain undetected.
- **Genome-Wide Heterogeneity vs. Functional Stratification:** Analyzing all promoter CGIs in aggregate masks sub-population signals.
	- *Hypothesis:* Housekeeping genes undergo strong purifying selection (conserved CGIs), whereas stress-response genes exhibit adaptive divergence. Stratifying CGIs by **GO (Gene Ontology) terms** is expected to resolve species-specific patterns.

### 2. Methodological & Machine Learning Refinements
- **Feature Representation ($k$-mer Spectrum):** $k=3$ was chosen for computational feasibility, but higher-order spatial motifs ($k=4..6$) or positional weight matrices (PWMs) may capture transcription factor binding sites more effectively.
- **Small-Sample Standardization:** The radar chart's z-scores were computed across only 3 data points (species medians), limiting it to a qualitative visualization rather than a statistically rigorous comparison.
- **Class Imbalance & Metrics:** Accuracy (~62%) was susceptible to majority-class bias (*B. juncea* gene abundance). Future models wil incorporate `Balanced Accuracy`, `Macro-F1`, and class-weighted sampling (SMOTE).
- **Hyperparameter Auditability:** Documenting parameter derivations (e.g., $3$-mer DBSCAN `min_samples`) highlighted the critical need for automated experiment tracking (e.g., MLflow/Hydra) in exploratory bioinformatics workflows.

## Repository Structure

```
├── 01_downloader.ipynb        # 1. Upstream sequence extraction
├── 02_preprocess.ipynb        # 2. CpG island detection
├── 03_feature_eda.ipynb       # 3. Data merge, EDA, 3-mer extraction, PCA
├── 04_clustering.ipynb        # 4. DBSCAN
├── 05_decision_tree.ipynb     # 5. Species classification
├── 06_feature_analysis.ipynb  # 6. Statistical validation & visualization
├── dataset/
│   ├── raw/                # Original FASTA/GFF files
│   └── processed/          # Pipeline intermediate/output CSVs
```
This repository provides the source code for the SteMClass app detailed in the study https://doi.org/10.1101/2025.09.02.673063. 
The app itself is hosted at https://stemclass.charite.de/. Users should **use the hosted app** but we are sharing the code for transparency and reproducibility.

## SteMClass

SteMClass is a DNA methylation-based classifier for the characterization of iPSC-derived cell differentiation state.

-   **Classification** using a pre-trained random forest with classifiction cut-off at probability **≥ 0.5** (below that, reported as _Not Classifiable_).
    
-   **Visualizations:** probability bar chart, UMAP projection against the training cohort, and heatmaps:
    
    -   **Global (cell-state marker) heatmap** using the top 10k CpGs that distinguish the selected state from iPSC.
        
    -   **Gene-level heatmap** for CpGs annotated to a user-specified gene.
       
-   **CNV profiling** from the uploaded IDAT pair (450K, EPICv1 or EPICv2), run independently of the classifier:
    -   **Output:** genome-wide CNV plot (labelled and unlabelled) and a table with one row per alteration, including its state (gain/loss/amplification/deletion), log2 ratio, coordinates (hg19 for 450K / EPICv1, hg38 for EPICv2). The stem cell hotspot column lists  known iPSC hotspots hit by each alteration, which are highlighted in green in the plot. Cancer gene census genes are labeled in black. The plots and the table, with hotspot citations, can be downloaded.
        

All model artifacts and reference data are **streamed at runtime** from GitHub Releases, no large local bundles are included in the repo.

## How to use the app (~10 minutes per sample)

1.  Open the **hosted SteMClass app** (https://stemclass.charite.de/).
    
2.  On the **Classification** tab:
    
    -   Upload **exactly two IDAT files** for a single sample (one _Red_ and one _Grn_).
        
    -   Enter a **Sample ID** and select the **array version** (450K, EPICv1, or EPICv2).
        
    -   Click **Preprocess Data**, then **Run Classification**.
        
3.  Review the **Prediction**, **Score**, and **Class Probabilities**.
    
4.  (Optional) **Download the HTML report**.
    
5.  Explore **Sample Visualization** (UMAP), **Cell State Global Beta Values** (marker heatmap), **Cell State Gene Level Beta Values** (gene heatmap) and **CNV Profiling** (DNA methylation array inferred CNV plot and alterations table).
    
6.  Use the **Download** buttons on the plot and table tabs to save.

7. To start over with a new sample, click **New Analysis**.    


## Outputs

-   **Classification result** with class probabilities (classification cut-off 0.5).
    
-   **UMAP** of test sample with the training cohort.
    
-   **Marker heatmap** (top 10k CpGs for selected state vs iPSC).
    
-   **Gene-level heatmap** (CpGs annotated to the given gene).
    
-   **HTML report** reproducing the classification-associated plots.
  
-   **CNV Plot** 
    

----------

## Methods (brief)

-   **Preprocessing:** `minfi::preprocessNoob`, genome mapping, and array compatibility handling. EPICv1/EPICv2 are projected onto the **450K feature space** for compatibility with the trained model.
    
-   **Model:** pre-fitted **random forest**; features aligned to the training set; **median imputation** via a `recipes` pipeline prepped on the training matrix.
    
-   **UMAP:** `uwot::umap` with `n_neighbors = 15`, `min_dist = 0.1`, Euclidean distance, seeded with `set.seed(123)`.

-    **CNV Profiling** SeSAMe (`QCDPB` preprocessing) intensities fitted with conumee2 (https://github.com/hovestadtlab/conumee2) against an array-specific control reference (female and male controls pooled). Copy number is reported for the autosomes (chr1–22). Sex is inferred with SeSAMe `inferSex`. Segments with a gain (log2 ≥ 0.15 ), loss (log2 ≤ -0.17 ), amplifications (AMP, log2 ≥ 0.6 ) and homozygous deletions (HOMDEL, log2 < -0.64 and 0.3 below the background segment). The stem cell hotspot column lists known iPSC hotspots (doi:10.64898/2025.12.01.691615) hit by each alteration, which are highlighted in green in the plot. Cancer gene census genes are labeled in black.
-   **Filtering:** segments covered by fewer than 150 probes and AMP / HOMDEL genes covered by fewer than 10 probes are not reported.
    
----------

## Data & model (fetched at runtime)

At runtime the app streams:

-   Trained model (`final_rf_fit_no_cal.rds`)
    
-   Training matrix (`final_BIH_train_data.txt`) and sample annotations (`final_BIH_train_targets.txt`)
    
-   CpG annotation with top-marker sets (`CpG_450k_annotation_with_top10k_marker.txt`)

-   CNV cohorts of 450K / EPICv1 / EPICv2 specific normal reference (`ref_Conumee.rds`)

----------

## Scope & limitations

-   Designed for **one sample per run** via the UI.

-   Currently the model characterizes 8 distinct differentiation states (see https://doi.org/10.1101/2025.09.02.673063 for details).
    
-   The **probability threshold** (0.5) governs the reject decision.
    
-   Gene-level heatmaps rely on **Illumina 450K gene annotations**.
    
-   Intended as a **research tool**; not for clinical use.
    

----------

## Reproducibility & transparency

-   We publish the exact Shiny code used in the hosted app, including the UI theme and plotting code paths used for the downloadable report.
    
-   All model and reference objects are versioned on GitHub Releases.
    
-   Random elements (UMAP) are seeded to provide consistent embeddings across runs.
    

----------

## Privacy

-   IDATs are read for the active session only and not stored. Outputs are generated on demand.
    
-   When using the hosted app, please follow your institution’s data governance and de-identification policies.

----------

## Issue reporting

If you encounter a **bug affecting the hosted app**, please open an issue with:

-   A short description of the problem and the tab you were using.
    
-   Array type, approximate IDAT size, and the exact error message (if shown).
    
-   Browser/OS and timestamp (UTC) of the attempt.

----------

## Citation

If SteMClass is useful in your work, please cite the associated article (https://doi.org/10.1101/2025.09.02.673063) and acknowledge the app.

----------

## License

This repository (including the Shiny app code and trained model) is released under the
**Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)**.

- You may use, share, and adapt the materials for **non-commercial academic research**.
- **Commercial use is not permitted** without explicit written permission from the authors.
- Full license text: [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

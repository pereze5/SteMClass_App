This repository provides the source code for the SteMClass Shiny app detailed in the study https://doi.org/10.1101/2025.09.02.673063. 
The app itself is hosted at https://pereze5-stemclass.share.connect.posit.cloud/; readers should **use the hosted app** but we are sharing the code for transparency and reproducibility.

## SteMClass

SteMClass is a DNA methylation-based classifier for the characterisation of iPSC-derived cell differentiation state.

-   **Classification** using a pre-trained random forest with classifiction cut-off at probability **≥ 0.6** (below that, reported as _Not Classifiable_).
    
-   **Visualizations:** probability bar chart, UMAP projection against the training cohort, and heatmaps:
    
    -   **Global (cell-state marker) heatmap** using the top 10k CpGs that distinguish the selected state from iPSC.
        
    -   **Gene-level heatmap** for CpGs annotated to a user-specified gene.
        

All model artifacts and reference data are **streamed at runtime** from GitHub Releases; no large local bundles are included in the repo.

## How to use the app (~10 minutes per sample)

1.  Open the **hosted SteMClass app** (https://pereze5-stemclass.share.connect.posit.cloud/).
    
2.  On the **Classification** tab:
    
    -   Upload **exactly two IDAT files** for a single sample (one _Red_ and one _Grn_).
        
    -   Enter a **Sample ID** and select the **array version** (450K, EPICv1, or EPICv2).
        
    -   Click **Preprocess Data**, then **Run Classification**.
        
3.  Review the **Prediction**, **Score**, and **Class Probabilities**.
    
4.  (Optional) **Download the HTML report**.
    
5.  Explore **Sample Visualization** (UMAP), **Cell State Global Beta Values** (marker heatmap), and **Cell State Gene Level Beta Values** (gene heatmap).
    
6.  Use the **Download** buttons on the heatmap tabs to save as PNGs.

7. To start over with a new sample, click **New Analysis**.    

**Upload size:** The app currently accepts uploads up to ~30 MB per session.


## Outputs

-   **Classification result** with class probabilities (reject below 0.6).
    
-   **Probability bar chart** (cut-off line at 0.6).
    
-   **UMAP** of test sample with the training cohort.
    
-   **Marker heatmap** (top 10k CpGs for selected state vs iPSC).
    
-   **Gene-level heatmap** (CpGs annotated to the given gene).
    
-   **HTML report** reproducing the key plots.
    

----------

## Methods (brief)

-   **Preprocessing:** `minfi::preprocessNoob`, genome mapping, and array compatibility handling. EPICv1/EPICv2 are projected onto the **450K feature space** for compatibility with the trained model.
    
-   **Model:** pre-fitted **random forest**; features aligned to the training set; **median imputation** via a `recipes` pipeline prepped on the training matrix.
    
-   **UMAP:** `uwot::umap` with `n_neighbors = 15`, `min_dist = 0.1`, Euclidean distance, seeded with `set.seed(123)`.
    
----------

## Data & model (fetched at runtime)

At runtime the app streams:

-   Trained model (`final_rf_fit_no_cal.rds`)
    
-   Training matrix (`final_BIH_train_data.txt`) and sample annotations (`final_BIH_train_targets.txt`)
    
-   CpG annotation with top-marker sets (`CpG_450k_annotation_with_top10k_marker.txt`)

----------

## Scope & limitations

-   Designed for **one sample per run** via the UI.

-   Currently the model characterizes 8 distinct differentiation states (see https://doi.org/10.1101/2025.09.02.673063 for details).
    
-   The **probability threshold** (0.6) governs the reject decision.
    
-   Gene-level heatmaps rely on **Illumina 450K gene annotations**.
    
-   Intended as a **research tool**; not for clinical use.
    

----------

## Reproducibility & transparency

-   We publish the exact Shiny code used in the hosted app, including the UI theme and plotting code paths used for the downloadable report.
    
-   All model and reference objects are versioned on GitHub Releases.
    
-   Random elements (UMAP) are seeded to provide consistent embeddings across runs.
    

----------

## Privacy

-   IDATs are read for the active session only and not stored; outputs are generated on demand.
    
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

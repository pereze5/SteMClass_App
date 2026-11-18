options(shiny.maxRequestSize = 80*1024^2)

library(shinyjs)
library(thematic)     # For automatic theming of ggplot2 plots
library(bslib)        # For Bootstrap theming
library(shinycssloaders)
library(RColorBrewer)
library(circlize)
library(dplyr)
library(tidyr)

# Load the trained model, training means, and CpG annotation table

# ==== Streamed file URLs ====
urls <- list(
  model        = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/final_rf_fit_prob.rds",
  train_anno   = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/final_train_targets.txt",
  marker_data   = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/final_BIH_marker_data.txt",
  cpg_anno     = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/CpG_450k_annotation_with_top10k_marker.txt",
  train_data = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/final_train_data.txt",
  cal_model = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/calibration_model_deploy.rds",
  prepped_rec = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/prepped_recipe.rds",
  train_gmset = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/gmset_train_2025_flt.rda"
)

# 9‐color Blues palette; interpolate to 255

blues_base <- brewer.pal(9, "Blues")
blues_pal  <- colorRampPalette(blues_base)(255)

heatmap_scale_blues <- colorRamp2(
  seq(0, 1, length.out = 255),
  blues_pal
)

REJECT_LABEL <- "reject"


example_samples <- list(
  list(
    label         = "iPSC Test Sample (GSM9238338)",
    accession     = "GSM9238338 iPSC",
    array_version = "EPICv1",
    grn_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/206687210042_R03C01_Grn.idat",
    red_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/206687210042_R03C01_Red.idat"
  ),
  list(
    label         = "Ectoderm Test Sample (GSM9238401)",
    accession     = "GSM9238401 Ectoderm",
    array_version = "EPICv1",
    grn_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/204964480134_R08C01_Grn.idat",
    red_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/204964480134_R08C01_Red.idat"
  ),
  list(
    label         = "Endoderm Test Sample (GSM9238413)",
    accession     = "GSM9238413 Endoderm",
    array_version = "EPICv1",
    grn_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/207705780010_R07C01_Grn.idat",
    red_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/207705780010_R07C01_Red.idat"
  ),
  list(
    label         = "Mesoderm Test Sample (GSM9238472)",
    accession     = "GSM9238472 Mesoderm",
    array_version = "EPICv1",
    grn_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/207705780010_R06C01_Grn.idat",
    red_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/207705780010_R06C01_Red.idat"
  ),
  list(
    label         = "Endothelial Test Sample (GSM9238440)",
    accession     = "GSM9238440 Endothelial",
    array_version = "EPICv1",
    grn_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/206702500016_R07C01_Grn.idat",
    red_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/206702500016_R07C01_Red.idat"
  ),
  list(
    label         = "Lung Test Sample (GSM9238447)",
    accession     = "GSM9238447 Lung",
    array_version = "EPICv2",
    grn_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/208010430001_R03C01_Grn.idat",
    red_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/208010430001_R03C01_Red.idat"
  ),
  list(
    label         = "NSC Test Sample (GSM9238490)",
    accession     = "GSM9238490 NSC",
    array_version = "EPICv1",
    grn_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/205038940051_R02C01_Grn.idat",
    red_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/205038940051_R02C01_Red.idat"
  ),
  list(
    label         = "Astrocyte Test Sample (GSM9238374)",
    accession     = "GSM9238374 Astrocyte",
    array_version = "EPICv2",
    grn_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/208010440130_R04C01_Grn.idat",
    red_url       = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/208010440130_R04C01_Red.idat"
  )
)


# Named vector for the selectInput choices (name = label, value = accession)
example_choices <- c(
  "- select a sample -" = "",
  setNames(
    sapply(example_samples, `[[`, "accession"),
    sapply(example_samples, `[[`, "label")
  )
)


# ---- helper: stable softmax (in case calibrator returns logits)
softmax_mat <- function(eta) {
  eta <- as.matrix(eta)
  eta <- eta - apply(eta, 1, max)            # stability
  exp_eta <- exp(eta)
  exp_eta / rowSums(exp_eta)
}

# ---- helper: apply calibrator (glmnet OR tidymodels workflow) ----
apply_calibrator <- function(cal_obj, raw_prob_df, class_names) {
  # raw_prob_df must have columns == class_names in that order
  stopifnot(all(class_names %in% colnames(raw_prob_df)))
  
  X <- as.matrix(raw_prob_df[, class_names, drop = FALSE])
  
  # Case A: tidymodels workflow (or parsnip fit) that supports predict(type="prob")
  if (inherits(cal_obj, "workflow")) {
    # workflow predict() usually returns tibble of class probs
    p <- predict(cal_obj, new_data = as.data.frame(X), type = "prob")
    # Ensure column names are exactly class_names
    colnames(p) <- gsub("^\\.pred_", "", colnames(p))
    p <- p[, class_names, drop = FALSE]
    return(as_tibble(p))
  }
  
  # Case B: glmnet model (cv.glmnet or glmnet)
  if (inherits(cal_obj, "cv.glmnet") || inherits(cal_obj, "glmnet")) {
    # Use lambda.1se if available (per your description); else lambda.min; else s not needed
    s_val <- NULL
    if (!is.null(cal_obj$lambda.1se)) s_val <- "lambda.1se"
    if (is.null(s_val) && !is.null(cal_obj$lambda.min)) s_val <- "lambda.min"
    
    # Prefer type="response" for multinomial, which returns probabilities
    p_arr <- if (!is.null(s_val)) {
      predict(cal_obj, newx = X, s = s_val, type = "response")
    } else {
      predict(cal_obj, newx = X, type = "response")
    }
    
    # glmnet multinomial predict(type="response") returns an array: n x K x 1
    if (length(dim(p_arr)) == 3) p_mat <- p_arr[, , 1, drop = FALSE] else p_mat <- p_arr
    
    # Set class names (glmnet often carries dimnames; but enforce)
    colnames(p_mat) <- class_names
    return(as_tibble(p_mat))
  }
  
  stop("Unsupported calibration object type: ", paste(class(cal_obj), collapse = ", "))
}

predict_calibrated_all <- function(
    beta_mat_cpg_by_sample,
    train_data,
    prepped_rec,
    rf_fit_path,
    cal_path,
    class_names = c("Astrocyte","Ectoderm","Endoderm","Endothelial","iPSC","Lung","Mesoderm","NSC"),
    threshold = 0.5,
    apply_reject = TRUE
) {
  if (!requireNamespace("workflows", quietly = TRUE)) stop("Package 'workflows' is required.")
  if (!requireNamespace("parsnip", quietly = TRUE))   stop("Package 'parsnip' is required.")
  
  # 1) Expected CpGs from training
  # After prep(rec) exists:
  expected_cpgs <- recipes::juice(prepped_rec, all_predictors()) %>% colnames()
  
  # 2) Convert to samples x CpGs and align/insert NA for missing CpGs
  X <- t(beta_mat_cpg_by_sample) %>% as.data.frame()  # samples x CpGs
  missing_cpgs <- setdiff(expected_cpgs, colnames(X))
  if (length(missing_cpgs) > 0) {
    X[, missing_cpgs] <- NA_real_
  }
  X <- X[, expected_cpgs, drop = FALSE]
  
  # 3) Prepare new_data for bake: predictors only, numeric, aligned
  X_tbl <- X %>%
    tibble::rownames_to_column("Sample_accession") %>%
    tibble::as_tibble()
  
  # Force all CpG predictors to numeric (protects step_impute_median)
  # This will turn "NA", "", etc. into NA_real_ cleanly
  X_tbl[expected_cpgs] <- lapply(X_tbl[expected_cpgs], function(z) as.numeric(z))
  
  X_tbl <- X_tbl %>% tibble::column_to_rownames("Sample_accession")
  
  # Bake (no outcome column required)
  baked <- recipes::bake(prepped_rec, new_data = X_tbl)
  
  
  # 4) Load models from the paths/connections supplied by the caller
  rf_bundle  <- readRDS(rf_fit_path)
  cal_deploy <- tryCatch(
    readRDS(cal_path),
    error = function(e) stop("Failed to load calibration model: ", e$message)
  )
  
  rf_fit <- if (is.list(rf_bundle) && "model_fit" %in% names(rf_bundle)) rf_bundle$model_fit else rf_bundle
  
  # 5) RF raw probs (workflow-safe)
  raw_probs <- predict(rf_fit, new_data = baked, type = "prob") %>% tibble::as_tibble()
  colnames(raw_probs) <- gsub("^\\.pred_", "", colnames(raw_probs))
  
  missing <- setdiff(class_names, colnames(raw_probs))
  if (length(missing) > 0) stop("RF probability output missing classes: ", paste(missing, collapse = ", "))
  raw_probs <- raw_probs[, class_names, drop = FALSE]
  
  # 6) Calibrated probs
  cal_obj <- cal_deploy
  if (is.list(cal_obj) && "cv" %in% names(cal_obj)) cal_obj <- cal_obj$cv
  cal_probs <- apply_calibrator(cal_obj, raw_probs, class_names)
  
  # 7) Assemble output (one row per sample)
  # After apply_calibrator(), force bare class names
  colnames(cal_probs) <- class_names
  
  out <- dplyr::bind_cols(
    tibble::tibble(Sample_accession = rownames(baked)),
    raw_probs %>% dplyr::rename_with(~ paste0("raw_", .x)),
    cal_probs %>% dplyr::rename_with(~ paste0("cal_", .x))
  )
  
  if (apply_reject) {
    cal_mat <- as.matrix(cal_probs)
    max_prob <- apply(cal_mat, 1, max)
    argmax   <- class_names[max.col(cal_mat, ties.method = "first")]
    pred_class <- ifelse(max_prob < threshold, REJECT_LABEL, argmax)
    
    out <- dplyr::mutate(
      out,
      cal_max_prob = max_prob,
      cal_pred = factor(pred_class, levels = c(class_names, REJECT_LABEL))
    )
  }
  
  out
}

# ==============================================================================
# CNV PROFILING CONFIGURATION + HELPERS
# ==============================================================================
# The CNV workflow uses SeSAMe + conumee2 on the original IDAT intensities.
# Reference/annotation paths can be overridden with environment variables.
# Defaults below preserve the paths used by the supplied standalone CNV script.

CNV_RELEASE_BASE <- "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data"

cnv_config <- list(
  ensemble_hg19 = paste0(CNV_RELEASE_BASE, "/ensemble_gr.rds"),
  ensemble_hg38 = paste0(CNV_RELEASE_BASE, "/ensemble_gr_hg38.rds"),
  cancer_gene_census = paste0(
    CNV_RELEASE_BASE,
    "/CancerGeneCensus_allWed_May_4-08_41_12_2022.csv"
  ),
  detailplot_override = paste0(CNV_RELEASE_BASE, "/cnv_detailplot_override.R"),
  
  # StemCNV-check curated stem cell aberration hotspots (cytobands + genes),
  # shown in the detail plot annotation together with the CGC hit genes.
  stemCNV_hotspots = paste0(CNV_RELEASE_BASE, "/CNV-hotspots.csv"),
  
  # UCSC cytoband tables, needed to turn hotspot cytobands (e.g. "20q11.21")
  # into genomic coordinates for each genome build. For a fully self-contained
  # app, upload both files to the release and point these URLs there.
  cytoband_hg19 = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
  cytoband_hg38 = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
  
  # Existing sex-separated conumee references in the release. For each array
  # the female and male references are pooled at run time into ONE mixed-sex
  # reference: sex does not affect autosomal copy number, so pooling maximises
  # the number of controls. chrX/chrY are not analysed.
  ref_450k_female   = paste0(CNV_RELEASE_BASE, "/ref450K_female_Conumee.rds"),
  ref_450k_male     = paste0(CNV_RELEASE_BASE, "/ref450K_male_Conumee.rds"),
  ref_epicv1_female = paste0(CNV_RELEASE_BASE, "/refEPIC_female_Conumee.rds"),
  ref_epicv1_male   = paste0(CNV_RELEASE_BASE, "/refEPIC_male_Conumee.rds"),
  ref_epicv2_female = paste0(CNV_RELEASE_BASE, "/refEPICv2_female_Conumee.rds"),
  ref_epicv2_male   = paste0(CNV_RELEASE_BASE, "/refEPICv2_male_Conumee.rds")
)

# All CNV analysis uses SeSAMe with the QCDPB preprocessing string, matching the
# release references:
#   Q = quality masking, C = Infinium-I channel inference, D = non-linear dye-bias
#   correction, P = pOOBAH detection masking, B = noob background subtraction.
# Masked probes (Q, P) return NA total intensity and are dropped from the fit.
CNV_SESAME_PREP <- "QCDPB"

# Fixed CNV reference sex (no user selection).
CNV_REFERENCE_SEX <- "mixed"

cnv_reference_label <- function(platform) paste("Mixed-sex", platform, "reference")

# Autosomes only: against a mixed-sex reference chrX/chrY ratios are not
# interpretable, so they are excluded from fitting, segmentation and plots.
CNV_PLOT_CHRS <- paste0("chr", 1:22)

# Focal CNV call thresholds (focal log2 ratio).
CNV_AMP_THRESHOLD    <- 0.6     # focal >= this -> AMP
CNV_HOMDEL_THRESHOLD <- -0.64   # focal <  this -> HOMDEL

# Segment states (same cut-offs as GAIN / LOSS in the focal-call logic). Every
# autosomal segment beyond these is reported in the alterations table; a
# StemCNV-check hotspot is "hit" when such a segment overlaps it in the curated
# direction, or when a high-confidence focal call overlaps it.
CNV_HOTSPOT_GAIN <- 0.15    # segment mean >= this -> GAIN (>= CNV_AMP_THRESHOLD -> AMP)
CNV_HOTSPOT_LOSS <- -0.17   # segment mean <= this -> LOSS (<  CNV_HOMDEL_THRESHOLD -> HOMDEL)

# Sex is predicted with sesame::inferSex() on the same QCDPB-processed data used
# for CNV. If that fails, minfi::getSex() on the classifier's own-array Noob
# data is used as a fallback: yMed - xMed < -2 -> Female (calls within the
# margin below of -2 are marked "borderline").
CNV_SEX_CUTOFF <- -2
CNV_SEX_BORDERLINE_MARGIN <- 0.5

# Genome-plot colour scale (identical to the standalone script).
# y-range of the genome plot: the profile uses +/-1.25 (conumee default); the
# extra space holds the non-overlapping label rows.
CNV_PLOT_YLIM <- c(-2.1, 2.1)

CNV_PLOT_COLS <- c(
  "#39379F", "#39379F", "#39379F", "#39379F",
  "#ECF0F1",
  "#E17E49", "#E17E49", "#E17E49", "#E17E49"
)

# Cache remote CNV assets once per R process.
cnv_asset_cache_dir <- file.path(tempdir(), "stemclass_cnv_assets")
dir.create(cnv_asset_cache_dir, recursive = TRUE, showWarnings = FALSE)

cnv_is_url <- function(x) {
  is.character(x) &&
    length(x) == 1L &&
    nzchar(x) &&
    grepl("^https?://", x, ignore.case = TRUE)
}

cnv_materialize_asset <- function(path_or_url, label = "CNV asset") {
  if (!is.character(path_or_url) || length(path_or_url) != 1L || !nzchar(path_or_url)) {
    stop(label, " is not configured.", call. = FALSE)
  }
  
  if (!cnv_is_url(path_or_url)) {
    if (!file.exists(path_or_url)) {
      stop(label, " was not found: ", path_or_url, call. = FALSE)
    }
    return(path_or_url)
  }
  
  url_no_query <- sub("\\?.*$", "", path_or_url)
  filename <- basename(url_no_query)
  if (!nzchar(filename) || identical(filename, "/")) {
    filename <- paste0("cnv_asset_", abs(sum(utf8ToInt(path_or_url))))
  }
  
  dest <- file.path(cnv_asset_cache_dir, filename)
  
  if (!file.exists(dest) || is.na(file.info(dest)$size) || file.info(dest)$size == 0) {
    tmp <- paste0(dest, ".part")
    if (file.exists(tmp)) unlink(tmp)
    
    message("Downloading ", label, ": ", filename)
    
    status <- tryCatch(
      utils::download.file(
        url = path_or_url,
        destfile = tmp,
        mode = "wb",
        quiet = TRUE
      ),
      error = function(e) e
    )
    
    if (inherits(status, "error") || !file.exists(tmp) ||
        is.na(file.info(tmp)$size) || file.info(tmp)$size == 0) {
      if (file.exists(tmp)) unlink(tmp)
      msg <- if (inherits(status, "error")) conditionMessage(status) else "empty download"
      stop("Failed to download ", label, " from ", path_or_url, ": ", msg, call. = FALSE)
    }
    
    if (!file.rename(tmp, dest)) {
      ok <- file.copy(tmp, dest, overwrite = TRUE)
      unlink(tmp)
      if (!ok) stop("Could not cache ", label, " at ", dest, call. = FALSE)
    }
  }
  
  dest
}

# Source the custom conumee detail-plot override from the GitHub release.
try({
  cnv_detailplot_file <- cnv_materialize_asset(
    cnv_config$detailplot_override,
    "CNV detail-plot override"
  )
  source(cnv_detailplot_file, local = FALSE)
}, silent = TRUE)

cnv_require_packages <- function() {
  required <- c(
    "sesame", "sesameData", "minfi", "conumee2", "GenomicRanges", "IRanges",
    "S4Vectors", "GenomeInfoDb", "DT"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "CNV profiling requires these packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

# Array-specific analysis on each platform's own probe set.
cnv_platform_config <- function(platform) {
  switch(
    platform,
    "450K"   = list(sesame = "HM450",  conumee = "450k",   genome = "hg19", use_exclude_regions = TRUE,  source_label = "450K"),
    "EPICv1" = list(sesame = "EPIC",   conumee = "EPIC",   genome = "hg19", use_exclude_regions = TRUE,  source_label = "EPICv1"),
    "EPICv2" = list(sesame = "EPICv2", conumee = "EPICv2", genome = "hg38", use_exclude_regions = FALSE, source_label = "EPICv2"),
    stop("Unsupported CNV platform: ", platform, call. = FALSE)
  )
}

prepare_ensembl_gr <- function(ensembl_gr) {
  mcols_names <- colnames(S4Vectors::mcols(ensembl_gr))
  gene_col <- dplyr::case_when(
    "name"          %in% mcols_names ~ "name",
    "gene"          %in% mcols_names ~ "gene",
    "gene_name"     %in% mcols_names ~ "gene_name",
    "symbol"        %in% mcols_names ~ "symbol",
    "external_name" %in% mcols_names ~ "external_name",
    TRUE ~ NA_character_
  )
  if (is.na(gene_col)) {
    stop(
      "Could not find a gene-symbol column in the Ensembl GRanges. Available mcols: ",
      paste(mcols_names, collapse = ", "),
      call. = FALSE
    )
  }
  S4Vectors::mcols(ensembl_gr)$name <- as.character(
    S4Vectors::mcols(ensembl_gr)[[gene_col]]
  )
  ensembl_gr
}

make_cgc_detail_regions <- function(ensembl_gr, cancer_gene_census, ref_genome) {
  ensembl_gr <- prepare_ensembl_gr(ensembl_gr)
  cgc_symbols <- unique(as.character(cancer_gene_census$Gene.Symbol))
  cgc_symbols <- setdiff(cgc_symbols, c("CDKN2A", "CDKN2B"))
  
  cgc_gr <- ensembl_gr[
    S4Vectors::mcols(ensembl_gr)$name %in% cgc_symbols
  ]
  
  if (identical(ref_genome, "hg19")) {
    cdkn2ab_range <- GenomicRanges::GRanges(
      seqnames = "chr9",
      ranges = IRanges::IRanges(start = 21967751, end = 22009304),
      strand = "*"
    )
  } else {
    cdkn2ab_range <- GenomicRanges::GRanges(
      seqnames = "chr9",
      ranges = IRanges::IRanges(start = 21967752, end = 21995301),
      strand = "*"
    )
  }
  
  S4Vectors::mcols(cdkn2ab_range)$name <- "CDKN2A/B"
  detail_gr <- c(cgc_gr, cdkn2ab_range)
  detail_gr <- GenomeInfoDb::keepStandardChromosomes(detail_gr, pruning.mode = "coarse")
  GenomeInfoDb::seqlevelsStyle(detail_gr) <- "UCSC"
  detail_gr[
    as.character(GenomicRanges::seqnames(detail_gr)) %in% paste0("chr", 1:22)
  ]
}

strip_conumee_sample_prefix <- function(df) {
  wanted_suffixes <- c(
    "ID", "chrom", "loc.start", "loc.end", "num.mark", "bstat", "pval",
    "seg.mean", "seg.median", "Chromosome", "Start", "End", "Name", "Value"
  )
  # Escape dots so "seg.mean" etc. match literally.
  suffix_pattern <- paste(gsub(".", "\\.", wanted_suffixes, fixed = TRUE), collapse = "|")
  # conumee prefixes columns with the sample name, e.g. "X2066..._R03C01.seg.mean"
  # or "GSM9238338 iPSC.seg.mean". Strip ANY prefix (spaces/dots allowed) that
  # precedes a known suffix; unprefixed names such as "seg.mean" are untouched.
  names(df) <- sub(
    paste0("^.+?\\.(", suffix_pattern, ")$"),
    "\\1",
    names(df)
  )
  df
}

cnv_reference_cache <- new.env(parent = emptyenv())
cnv_resource_cache  <- new.env(parent = emptyenv())

# ---- Progress helpers --------------------------------------------------------
# `progress` is a callback function(value, detail) where value is the absolute
# fraction complete (0-1). The default does nothing so the CNV functions also
# work outside Shiny.
cnv_no_progress <- function(value, detail) invisible(NULL)

# TRUE when an asset is local or already downloaded into the session cache.
cnv_asset_is_cached <- function(path_or_url) {
  if (!cnv_is_url(path_or_url)) return(TRUE)
  dest <- file.path(cnv_asset_cache_dir, basename(sub("\\?.*$", "", path_or_url)))
  file.exists(dest) && isTRUE(file.info(dest)$size > 0)
}

# Report "Downloading ..." or "Loading ... (cached)" and then fetch the asset.
cnv_fetch_with_progress <- function(path_or_url, label, progress, value) {
  fname <- basename(sub("\\?.*$", "", path_or_url))
  detail <- if (cnv_asset_is_cached(path_or_url)) {
    paste0("Loading ", label, " (cached): ", fname)
  } else {
    paste0("Downloading ", label, " (first run only): ", fname)
  }
  progress(value, detail)
  cnv_materialize_asset(path_or_url, label)
}

cnv_load_resources <- function(platform, progress = cnv_no_progress) {
  cfg <- cnv_platform_config(platform)
  key <- cfg$genome
  if (exists(key, envir = cnv_resource_cache, inherits = FALSE)) {
    progress(0.20, paste0("Gene annotation (", cfg$genome, ") already loaded"))
    return(get(key, envir = cnv_resource_cache, inherits = FALSE))
  }
  
  ensembl_src <- if (cfg$genome == "hg19") cnv_config$ensemble_hg19 else cnv_config$ensemble_hg38
  
  # Download (or locate) the annotation files; errors clearly if unavailable.
  ensembl_path <- cnv_fetch_with_progress(
    ensembl_src, paste0("Ensembl gene annotation (", cfg$genome, ")"), progress, 0.08
  )
  cgc_path <- cnv_fetch_with_progress(
    cnv_config$cancer_gene_census, "Cancer Gene Census table", progress, 0.14
  )
  
  hotspot_path <- cnv_fetch_with_progress(
    cnv_config$stemCNV_hotspots, "StemCNV-check hotspot list", progress, 0.15
  )
  cytoband_path <- cnv_fetch_with_progress(
    if (cfg$genome == "hg19") cnv_config$cytoband_hg19 else cnv_config$cytoband_hg38,
    paste0("UCSC cytoband table (", cfg$genome, ")"), progress, 0.16
  )
  
  out <- list(
    ensembl = readRDS(ensembl_path),
    cgc = read.csv(cgc_path, stringsAsFactors = FALSE),
    genome = cfg$genome
  )
  out$cytobands <- read_ucsc_cytobands(cytoband_path)
  out$hotspots <- make_stemcnv_hotspot_regions(
    hotspots = utils::read.csv(hotspot_path, stringsAsFactors = FALSE, check.names = FALSE),
    cytobands = out$cytobands,
    ensembl_gr = out$ensembl,
    genome = cfg$genome
  )
  assign(key, out, envir = cnv_resource_cache)
  out
}

# ---- StemCNV-check hotspot formatting ------------------------------------------
# UCSC cytoBand.txt(.gz): chrom, chromStart (0-based), chromEnd, band name, stain.
read_ucsc_cytobands <- function(path) {
  cb <- utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE,
                          col.names = c("chrom", "start", "end", "band", "stain"))
  cb$start <- cb$start + 1L                     # to 1-based, as GRanges
  cb$arm_band <- paste0(sub("^chr", "", cb$chrom), cb$band)   # e.g. "20q11.21"
  cb
}

# Cytobands covered by a hotspot label. Coarse labels include all sub-bands:
#   "1q25"    -> q25.1, q25.2, q25.3
#   "20q11.2" -> q11.21, q11.22, q11.23
#   "20q11.21"-> q11.21 only
cytoband_rows_for <- function(label, cytobands) {
  exact  <- cytobands$arm_band == label
  child  <- startsWith(cytobands$arm_band, paste0(label, "."))
  deeper <- grepl("\\.", label) & startsWith(cytobands$arm_band, label) &
    grepl("^[0-9]", substring(cytobands$arm_band, nchar(label) + 1))
  which(exact | child | deeper)
}

# Turn the curated CSV into one GRanges region per hotspot (genome-specific):
#   1. collapse duplicate rows (one row per citation) to one row per hotspot,
#      keeping call type(s), max check score and the list of citations;
#   2. map "gband" hotspots to cytoband coordinates, "gene_name" hotspots to
#      Ensembl gene coordinates;
#   3. keep autosomes only (chrX/chrY are not analysed) and UCSC naming.
make_stemcnv_hotspot_regions <- function(hotspots, cytobands, ensembl_gr, genome) {
  req_cols <- c("hotspot", "mapping", "call_type", "check_score", "citation")
  miss <- setdiff(req_cols, names(hotspots))
  if (length(miss)) {
    warning("[CNV] Hotspot list is missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
    return(GenomicRanges::GRanges())
  }
  
  hs <- hotspots
  hs$hotspot <- trimws(hs$hotspot)
  hs <- hs[nzchar(hs$hotspot), , drop = FALSE]
  
  # 1. one row per hotspot
  keys <- unique(hs[, c("hotspot", "mapping")])
  collapsed <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    d <- hs[hs$hotspot == keys$hotspot[i] & hs$mapping == keys$mapping[i], , drop = FALSE]
    data.frame(
      hotspot     = keys$hotspot[i],
      mapping     = keys$mapping[i],
      call_type   = paste(sort(unique(d$call_type)), collapse = "/"),
      check_score = suppressWarnings(max(as.numeric(d$check_score), na.rm = TRUE)),
      citations   = paste(unique(d$citation), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
  
  # 2a. cytoband hotspots
  gb <- collapsed[collapsed$mapping == "gband", , drop = FALSE]
  gb_gr <- lapply(seq_len(nrow(gb)), function(i) {
    rows <- cytoband_rows_for(gb$hotspot[i], cytobands)
    if (!length(rows)) return(NULL)
    chr <- unique(cytobands$chrom[rows])
    if (length(chr) != 1) return(NULL)
    GenomicRanges::GRanges(chr, IRanges::IRanges(min(cytobands$start[rows]), max(cytobands$end[rows])))
  })
  gb_ok <- !vapply(gb_gr, is.null, logical(1))
  
  # 2b. gene hotspots
  ens <- prepare_ensembl_gr(ensembl_gr)
  ens <- GenomeInfoDb::keepStandardChromosomes(ens, pruning.mode = "coarse")
  GenomeInfoDb::seqlevelsStyle(ens) <- "UCSC"
  gn <- collapsed[collapsed$mapping == "gene_name", , drop = FALSE]
  gn_gr <- lapply(gn$hotspot, function(g) {
    hit <- ens[S4Vectors::mcols(ens)$name == g]
    if (!length(hit)) return(NULL)
    hit <- hit[as.character(GenomicRanges::seqnames(hit)) == as.character(GenomicRanges::seqnames(hit))[1]]
    GenomicRanges::GRanges(as.character(GenomicRanges::seqnames(hit))[1],
                           IRanges::IRanges(min(GenomicRanges::start(hit)), max(GenomicRanges::end(hit))))
  })
  gn_ok <- !vapply(gn_gr, is.null, logical(1))
  
  unmapped <- c(gb$hotspot[!gb_ok], gn$hotspot[!gn_ok])
  if (length(unmapped)) {
    message("[CNV] Hotspots not found in the ", genome, " annotation (skipped): ",
            paste(unmapped, collapse = ", "))
  }
  
  mapped_rows <- rbind(gb[gb_ok, , drop = FALSE], gn[gn_ok, , drop = FALSE])
  if (!nrow(mapped_rows)) return(GenomicRanges::GRanges())
  gr <- suppressWarnings(do.call(c, c(gb_gr[gb_ok], gn_gr[gn_ok])))
  S4Vectors::mcols(gr) <- S4DataFrame_from(mapped_rows)
  S4Vectors::mcols(gr)$name <- mapped_rows$hotspot
  
  # 3. autosomes only
  gr <- gr[as.character(GenomicRanges::seqnames(gr)) %in% paste0("chr", 1:22)]
  message(sprintf("[CNV] StemCNV-check hotspots mapped to %s: %d regions (%d cytobands, %d genes)",
                  genome, length(gr), sum(S4Vectors::mcols(gr)$mapping == "gband"),
                  sum(S4Vectors::mcols(gr)$mapping == "gene_name")))
  gr
}

# Hotspots hit in this sample. Returns the hotspot GRanges subset with an
# added "evidence" column, e.g. "segment gain (0.27)" or "focal call: KRAS AMP".
find_hotspot_hits <- function(hotspots, segments, focal_regions = NULL) {
  if (!length(hotspots)) return(hotspots)
  ct <- S4Vectors::mcols(hotspots)$call_type
  want_gain <- grepl("gain", ct)
  want_loss <- grepl("loss", ct)
  evidence <- rep(NA_character_, length(hotspots))
  hit_state <- rep(NA_character_, length(hotspots))
  hit_value <- rep(NA_real_, length(hotspots))
  
  # 1. Segments overlapping the hotspot, in the curated direction.
  seg <- segments[is.finite(as.numeric(segments$seg.mean)), , drop = FALSE]
  if (nrow(seg)) {
    seg_gr <- GenomicRanges::GRanges(
      as.character(seg$chrom),
      IRanges::IRanges(as.integer(seg$loc.start), as.integer(seg$loc.end)),
      seg.mean = as.numeric(seg$seg.mean)
    )
    ov <- GenomicRanges::findOverlaps(hotspots, seg_gr)
    for (i in unique(S4Vectors::queryHits(ov))) {
      m <- S4Vectors::mcols(seg_gr)$seg.mean[S4Vectors::subjectHits(ov)[S4Vectors::queryHits(ov) == i]]
      if (want_gain[i] && any(m >= CNV_HOTSPOT_GAIN)) {
        evidence[i] <- sprintf("segment gain (%.2f)", max(m))
        hit_state[i] <- "GAIN"; hit_value[i] <- max(m)
      } else if (want_loss[i] && any(m <= CNV_HOTSPOT_LOSS)) {
        evidence[i] <- sprintf("segment loss (%.2f)", min(m))
        hit_state[i] <- "LOSS"; hit_value[i] <- min(m)
      }
    }
  }
  
  # 2. High-confidence focal calls overlapping the hotspot (any direction).
  if (!is.null(focal_regions) && length(focal_regions)) {
    ov <- GenomicRanges::findOverlaps(hotspots, focal_regions)
    for (i in unique(S4Vectors::queryHits(ov))) {
      calls <- S4Vectors::mcols(focal_regions)$label[S4Vectors::subjectHits(ov)[S4Vectors::queryHits(ov) == i]]
      txt <- paste0("focal call: ", paste(unique(calls), collapse = ", "))
      evidence[i] <- if (is.na(evidence[i])) txt else paste(evidence[i], txt, sep = "; ")
      if (is.na(hit_state[i])) hit_state[i] <- sub("^.* ", "", calls[1])
    }
  }
  
  keep <- !is.na(evidence)
  hit <- hotspots[keep]
  S4Vectors::mcols(hit)$evidence <- evidence[keep]
  S4Vectors::mcols(hit)$hit_state <- hit_state[keep]
  S4Vectors::mcols(hit)$hit_value <- hit_value[keep]
  hit
}

# Numeric chromosome order (chr1 < chr2 < ... < chr22 < chrX < chrY), so tables
# and downloads are never sorted alphabetically (chr1, chr10, chr11, ...).
chr_rank <- function(chr) {
  x <- sub("^chr", "", as.character(chr))
  x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x[x %in% c("M", "MT")] <- "25"
  suppressWarnings(as.integer(x))
}

order_by_genome <- function(df, chr_col, start_col) {
  if (!nrow(df)) return(df)
  df[order(chr_rank(df[[chr_col]]), as.numeric(df[[start_col]]), na.last = TRUE), , drop = FALSE]
}

# Cytoband span covered by a region, e.g. "20q11.1-q11.23" or "7q35".
segment_band_label <- function(chr, start, end, cytobands) {
  if (is.null(cytobands) || is.na(chr)) return(NA_character_)
  b <- cytobands[cytobands$chrom == chr & cytobands$end >= start & cytobands$start <= end, , drop = FALSE]
  if (!nrow(b)) return(NA_character_)
  b <- b[order(b$start), , drop = FALSE]
  first <- b$band[1]; last <- b$band[nrow(b)]
  paste0(sub("^chr", "", chr), if (first == last) first else paste0(first, "-", last))
}

# One harmonised table of all detected alterations:
#   - every autosomal SEGMENT beyond the GAIN / LOSS cut-offs (AMP / HOMDEL if
#     beyond the focal cut-offs), whether or not it overlaps a hotspot;
#   - every high-confidence FOCAL CALL (AMP / HOMDEL) in CGC genes and CDKN2A/B.
# The "Stem cell hotspot" column lists StemCNV-check hotspots hit by the row:
# for segments, hotspots overlapped in the curated direction; for focal calls,
# hotspots overlapping the gene.
build_alteration_table <- function(diagnostic, gene_regions, hotspots_all, segments, cytobands) {
  cols <- c("Alteration", "Source", "CNV state", "log2 ratio", "Chr", "Start", "End",
            "Size (kb)", "Stem cell hotspot", "Evidence", "Hotspot citations")
  empty <- stats::setNames(
    data.frame(character(), character(), character(), numeric(), character(), integer(),
               integer(), numeric(), character(), character(), character(),
               stringsAsFactors = FALSE),
    cols
  )
  hs_m <- if (length(hotspots_all)) S4Vectors::mcols(hotspots_all) else NULL
  
  hotspot_info <- function(idx) {
    if (!length(idx)) return(list(label = "No", cites = NA_character_))
    list(
      label = paste0("Yes: ", paste(paste0(hs_m$name[idx], " (", hs_m$call_type[idx], ")"), collapse = ", ")),
      cites = paste(unique(unlist(strsplit(hs_m$citations[idx], "; "))), collapse = "; ")
    )
  }
  make_row <- function(alt, src, state, val, chr, st, en, hs, ev) {
    stats::setNames(data.frame(
      alt, src, state, round(val, 3), chr, as.integer(st), as.integer(en),
      round((en - st + 1) / 1000, 1), hs$label, ev, hs$cites,
      stringsAsFactors = FALSE
    ), cols)
  }
  
  rows <- list()
  
  # --- segments ---------------------------------------------------------------
  seg <- segments
  seg$seg.mean <- as.numeric(seg$seg.mean)
  seg <- seg[seg$chrom %in% paste0("chr", 1:22) & is.finite(seg$seg.mean) &
               (seg$seg.mean >= CNV_HOTSPOT_GAIN | seg$seg.mean <= CNV_HOTSPOT_LOSS), , drop = FALSE]
  for (k in seq_len(nrow(seg))) {
    v <- seg$seg.mean[k]
    state <- if (v >= CNV_AMP_THRESHOLD) "AMP" else if (v >= CNV_HOTSPOT_GAIN) "GAIN"
    else if (v < CNV_HOMDEL_THRESHOLD) "HOMDEL" else "LOSS"
    gain <- v > 0
    chr <- as.character(seg$chrom[k]); st <- as.integer(seg$loc.start[k]); en <- as.integer(seg$loc.end[k])
    idx <- integer()
    if (length(hotspots_all)) {
      ov <- S4Vectors::queryHits(GenomicRanges::findOverlaps(
        hotspots_all, GenomicRanges::GRanges(chr, IRanges::IRanges(st, en))))
      dir_ok <- grepl(if (gain) "gain" else "loss", hs_m$call_type[ov])
      idx <- unique(ov[dir_ok])
    }
    n_probes <- if ("num.mark" %in% names(seg)) seg$num.mark[k] else NA
    rows[[length(rows) + 1]] <- make_row(
      alt = segment_band_label(chr, st, en, cytobands),
      src = "Segment", state = state, val = v, chr = chr, st = st, en = en,
      hs = hotspot_info(idx),
      ev = if (is.na(n_probes)) sprintf("segment mean %.2f", v)
      else sprintf("segment mean %.2f; %s probes", v, format(n_probes, big.mark = ","))
    )
  }
  
  # --- focal calls --------------------------------------------------------------
  for (i in seq_len(nrow(diagnostic))) {
    g <- diagnostic$V10[i]
    reg <- gene_regions[S4Vectors::mcols(gene_regions)$name == g]
    if (length(reg)) {
      chr <- as.character(GenomicRanges::seqnames(reg))[1]
      reg <- reg[as.character(GenomicRanges::seqnames(reg)) == chr]
      st <- min(GenomicRanges::start(reg)); en <- max(GenomicRanges::end(reg))
      idx <- if (length(hotspots_all)) unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(
        hotspots_all, GenomicRanges::GRanges(chr, IRanges::IRanges(st, en))))) else integer()
    } else {
      chr <- NA_character_; st <- NA_integer_; en <- NA_integer_; idx <- integer()
    }
    rows[[length(rows) + 1]] <- make_row(
      alt = g, src = "Focal call (CGC)", state = diagnostic$CNV_state[i],
      val = as.numeric(diagnostic$focal[i]), chr = chr, st = st, en = en,
      hs = hotspot_info(idx),
      ev = sprintf("focal log2 %.2f; background %.2f",
                   as.numeric(diagnostic$focal[i]), as.numeric(diagnostic$nonfocal[i]))
    )
  }
  
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  chr_num <- chr_rank(out$Chr)
  out[order(chr_num, out$Start, out$Source != "Segment", na.last = TRUE), , drop = FALSE]
}

S4DataFrame_from <- function(df) {
  S4Vectors::DataFrame(
    mapping = df$mapping, call_type = df$call_type,
    check_score = df$check_score, citations = df$citations
  )
}

# Array-specific mixed-sex reference, cached per platform for the R process.
cnv_reference_paths <- function(platform) {
  switch(
    platform,
    "450K"   = c(female = cnv_config$ref_450k_female,   male = cnv_config$ref_450k_male),
    "EPICv1" = c(female = cnv_config$ref_epicv1_female, male = cnv_config$ref_epicv1_male),
    "EPICv2" = c(female = cnv_config$ref_epicv2_female, male = cnv_config$ref_epicv2_male),
    stop("Unsupported CNV platform: ", platform, call. = FALSE)
  )
}

# Extract the intensity matrix from a saved reference WITHOUT touching its S4
# class. The release .rds files were created with the original 'conumee'
# package; calling is() / @ on them makes R try to load 'conumee'. S4 slots are
# stored as attributes, so the intensity table is read directly with attr().
# Only conumee2 is used for all analysis (the pooled reference is rebuilt with
# conumee2::CNV.load()).
cnv_extract_intensity <- function(obj, label = "reference") {
  if (is.matrix(obj) || is.data.frame(obj)) {
    m <- as.matrix(obj)
  } else {
    int <- attr(obj, "intensity", exact = TRUE)
    if (is.null(int)) {
      stop("Could not find an 'intensity' table in the ", label,
           " (class: ", paste(unclass(class(obj)), collapse = ", "), ").", call. = FALSE)
    }
    m <- as.matrix(as.data.frame(int, check.names = FALSE))
  }
  storage.mode(m) <- "double"
  m
}

# Load the female and male references for one array and pool them into one
# mixed-sex conumee reference (probes present in both; EPICv2 replicate IDs
# collapsed to cg/ch IDs). Cached per platform for the R process.
cnv_load_reference <- function(platform, progress = cnv_no_progress) {
  label <- cnv_reference_label(platform)
  cache_key <- paste(platform, "mixed", sep = "__")
  if (exists(cache_key, envir = cnv_reference_cache, inherits = FALSE)) {
    progress(0.25, paste(label, "already loaded"))
    return(get(cache_key, envir = cnv_reference_cache, inherits = FALSE))
  }
  
  paths <- cnv_reference_paths(platform)
  mats <- list()
  for (sx in names(paths)) {
    f <- cnv_fetch_with_progress(paths[[sx]], paste(sx, platform, "CNV reference"), progress,
                                 if (sx == "female") 0.18 else 0.21)
    obj <- readRDS(f)
    m <- cnv_extract_intensity(obj, paste(sx, platform, "reference"))
    storage.mode(m) <- "double"
    m <- collapse_v2_matrix(m)           # no-op for 450K / EPICv1 IDs
    colnames(m) <- paste0(sx, "_", colnames(m))
    mats[[sx]] <- m
  }
  
  # Probes with no measured value in any control carry no information.
  for (sx in names(mats)) {
    keep_rows <- rowSums(is.finite(mats[[sx]])) > 0
    mats[[sx]] <- mats[[sx]][keep_rows, , drop = FALSE]
  }
  n_f <- nrow(mats$female); n_m <- nrow(mats$male)
  shared <- intersect(rownames(mats$female), rownames(mats$male))
  message(sprintf(
    "[CNV] %s references: female %s probes x %d controls (e.g. %s); male %s probes x %d controls (e.g. %s); %s shared",
    platform,
    format(n_f, big.mark = ","), ncol(mats$female), paste(utils::head(rownames(mats$female), 2), collapse = ", "),
    format(n_m, big.mark = ","), ncol(mats$male),   paste(utils::head(rownames(mats$male), 2), collapse = ", "),
    format(length(shared), big.mark = ",")))
  
  # Pool only if the two references cover (nearly) the same probes. Compare with
  # the LARGER reference so a truncated / differently built file cannot shrink
  # the pooled probe set. Otherwise use the larger reference alone and say why.
  excluded_note <- NULL
  if (length(shared) >= 0.9 * max(n_f, n_m)) {
    pooled <- cbind(mats$female[shared, , drop = FALSE], mats$male[shared, , drop = FALSE])
    n_female <- ncol(mats$female); n_male <- ncol(mats$male)
  } else {
    use <- if (n_f >= n_m) "female" else "male"
    drop <- setdiff(c("female", "male"), use)
    excluded_note <- sprintf(
      "%s %s reference excluded: it shares only %s of %s probes with the %s reference (built on a different probe set)",
      drop, platform, format(length(shared), big.mark = ","),
      format(max(n_f, n_m), big.mark = ","), use)
    warning("[CNV] ", excluded_note, call. = FALSE)
    message("[CNV] ", excluded_note)
    pooled <- mats[[use]]
    n_female <- if (use == "female") ncol(mats$female) else 0L
    n_male   <- if (use == "male")   ncol(mats$male)   else 0L
  }
  
  ref <- conumee2::CNV.load(as.data.frame(pooled))
  attr(ref, "n_female") <- n_female
  attr(ref, "n_male") <- n_male
  attr(ref, "excluded_note") <- excluded_note
  
  progress(0.25, sprintf("%s ready: %d female + %d male controls%s", label, n_female, n_male,
                         if (!is.null(excluded_note)) paste0(" (", excluded_note, ")") else ""))
  assign(cache_key, ref, envir = cnv_reference_cache)
  ref
}


# EPICv2 replicate probe IDs (cg00000029_TC21) -> cg/ch IDs. Applied to the
# query, the reference (already collapsed by the builder) and the annotation,
# so all three always share plain cg/ch IDs.
strip_v2_suffix <- function(x) sub("_[A-Z]{2}[0-9]+$", "", x)

collapse_v2_matrix <- function(mat) {
  mat <- as.matrix(mat)
  ids <- strip_v2_suffix(rownames(mat))
  if (identical(ids, rownames(mat))) return(mat)
  sums <- rowsum(mat, ids, na.rm = TRUE)
  ns   <- rowsum((!is.na(mat)) * 1, ids)
  out  <- sums / ns
  out[ns == 0] <- NA
  out
}

harmonise_anno_ids <- function(anno) {
  ids <- strip_v2_suffix(names(anno@probes))
  if (identical(ids, names(anno@probes))) return(anno)
  names(anno@probes) <- ids
  anno@probes <- anno@probes[!duplicated(ids)]
  anno
}

# Sex from the QCDPB-processed SeSAMe data (sesame::inferSex); minfi fallback.
cnv_predict_sex_sesame <- function(sdf, mset = NULL) {
  sx <- tryCatch(as.character(sesame::inferSex(sdf))[1], error = function(e) NA_character_)
  if (!is.na(sx) && toupper(sx) %in% c("FEMALE", "MALE", "F", "M")) {
    sex <- if (toupper(sx) %in% c("FEMALE", "F")) "Female" else "Male"
    return(list(sex = sex, diff = NA_real_,
                note = paste0("sesame inferSex (", CNV_SESAME_PREP, "-processed data)")))
  }
  if (!is.null(mset)) {
    res <- cnv_predict_sex(mset)
    res$note <- paste0(res$note)
    return(res)
  }
  list(sex = "Indeterminate", diff = NA_real_, note = "sex prediction failed")
}

# Sex from the sample's own-array Noob MethylSet (minfi::getSex rule).
cnv_predict_sex <- function(mset) {
  out <- tryCatch({
    s <- minfi::getSex(minfi::mapToGenome(mset), cutoff = CNV_SEX_CUTOFF)
    list(xMed = as.numeric(s$xMed[1]), yMed = as.numeric(s$yMed[1]))
  }, error = function(e) {
    # Fallback: same rule computed directly (chromosomes from the array annotation).
    tryCatch({
      ann <- minfi::getAnnotation(mset)
      cn  <- log2(pmax(minfi::getMeth(mset)[, 1] + minfi::getUnmeth(mset)[, 1], 1))
      chr <- as.character(ann[rownames(mset), "chr"])
      list(xMed = stats::median(cn[chr == "chrX"], na.rm = TRUE),
           yMed = stats::median(cn[chr == "chrY"], na.rm = TRUE))
    }, error = function(e2) NULL)
  })
  if (is.null(out) || !is.finite(out$xMed) || !is.finite(out$yMed)) {
    return(list(sex = "Indeterminate", diff = NA_real_,
                note = "sex prediction failed (chrX/chrY probes unavailable)"))
  }
  diff <- out$yMed - out$xMed
  sex <- if (diff < CNV_SEX_CUTOFF) "Female" else "Male"
  if (abs(diff - CNV_SEX_CUTOFF) < CNV_SEX_BORDERLINE_MARGIN) sex <- paste(sex, "(borderline)")
  list(
    sex = sex,
    diff = diff,
    note = sprintf("minfi getSex: yMed - xMed = %.2f (female < %.0f)", diff, CNV_SEX_CUTOFF)
  )
}

# Keep only annotation probes present in `ids` (reference probes, or probes with
# a measured value in the query). Bins left with too few probes are removed so
# CNV.bin / CNV.segment never see empty bins. Stops if coverage is below 50%.
restrict_anno_to_probes <- function(anno, ids, what, bin_minprobes = 15) {
  keep <- names(anno@probes) %in% ids
  frac <- mean(keep)
  message(sprintf("[CNV] %s covers %.1f%% of annotation probes (%s of %s)",
                  what, 100 * frac, format(sum(keep), big.mark = ","),
                  format(length(keep), big.mark = ",")))
  if (frac < 0.5) {
    stop(sprintf(paste0(
      "The %s only covers %.1f%% of the annotation probes (example IDs: %s). ",
      "Check that the release references for this array use the same probe IDs as the sample."),
      what, 100 * frac, paste(utils::head(ids, 3), collapse = ", ")), call. = FALSE)
  }
  if (all(keep)) return(anno)
  anno@probes <- anno@probes[keep]
  n_per_bin <- GenomicRanges::countOverlaps(anno@bins, anno@probes)
  anno@bins <- anno@bins[n_per_bin >= bin_minprobes]
  anno
}

restrict_anno_to_ref <- function(anno, ref, bin_minprobes = 15) {
  restrict_anno_to_probes(anno, rownames(ref@intensity), "CNV reference", bin_minprobes)
}

# Detail regions used for plotting: CGC genes (without individual CDKN2A/B)
# restricted to the genes that were called, as in the standalone script.
collapse_cdkn2ab_hits <- function(x) {
  x <- unique(stats::na.omit(as.character(x)))
  if ("CDKN2A/B" %in% x || any(x %in% c("CDKN2A", "CDKN2B"))) {
    x <- setdiff(x, c("CDKN2A", "CDKN2B"))
    x <- unique(c(x, "CDKN2A/B"))
  }
  x
}

cnv_build_annotation <- function(platform, resources, bin_size = 50000L,
                                 reference = NULL, detail_regions = NULL) {
  cfg <- cnv_platform_config(platform)
  if (is.null(detail_regions)) {
    detail_regions <- make_cgc_detail_regions(
      resources$ensembl,
      resources$cgc,
      cfg$genome
    )
  }
  
  excluded <- NULL
  if (isTRUE(cfg$use_exclude_regions)) {
    env <- new.env(parent = emptyenv())
    utils::data("exclude_regions", package = "conumee2", envir = env)
    if (exists("exclude_regions", envir = env, inherits = FALSE)) {
      excluded <- get("exclude_regions", envir = env, inherits = FALSE)
    }
  }
  
  anno <- conumee2::CNV.create_anno(
    bin_minprobes = 15,
    exclude_regions = excluded,
    array_type = cfg$conumee,
    genome = cfg$genome,
    detail_regions = detail_regions,
    bin_minsize = as.integer(bin_size),
    chrXY = FALSE   # autosomes only (mixed-sex reference)
  )
  anno <- harmonise_anno_ids(anno)
  
  # Ensure every annotation probe exists in the array-matched reference.
  if (!is.null(reference)) {
    anno <- restrict_anno_to_ref(anno, reference)
  }
  anno
}

cnv_read_written_table <- function(path) {
  read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

cnv_harmonise_segments <- function(segment_df) {
  segment_df <- strip_conumee_sample_prefix(segment_df)
  
  if (!"seg.mean" %in% names(segment_df)) {
    candidate <- intersect(
      c("seg.mean", "seg_mean", "mean", "Mean", "value", "Value", "log2", "log2ratio"),
      names(segment_df)
    )[1]
    if (is.na(candidate)) stop("Could not identify the conumee segment-mean column.")
    segment_df$seg.mean <- as.numeric(segment_df[[candidate]])
  }
  if (!"chrom" %in% names(segment_df)) {
    candidate <- intersect(c("chrom", "Chromosome", "chr", "seqnames"), names(segment_df))[1]
    if (is.na(candidate)) stop("Could not identify the conumee chromosome column.")
    segment_df$chrom <- as.character(segment_df[[candidate]])
  }
  if (!"loc.start" %in% names(segment_df)) {
    candidate <- intersect(c("loc.start", "Start", "start", "loc_start"), names(segment_df))[1]
    if (is.na(candidate)) stop("Could not identify the conumee segment-start column.")
    segment_df$loc.start <- as.integer(segment_df[[candidate]])
  }
  if (!"loc.end" %in% names(segment_df)) {
    candidate <- intersect(c("loc.end", "End", "end", "loc_end"), names(segment_df))[1]
    if (is.na(candidate)) stop("Could not identify the conumee segment-end column.")
    segment_df$loc.end <- as.integer(segment_df[[candidate]])
  }
  segment_df
}

# Genes used to annotate segments (as ampl_del_gr_725 in the standalone script):
# ALL Cancer Gene Census genes plus the combined CDKN2A/B region.
make_cgc_segment_genes <- function(ensembl_gr, cancer_gene_census, ref_genome) {
  ensembl_gr <- prepare_ensembl_gr(ensembl_gr)
  cgc_gr <- ensembl_gr[
    S4Vectors::mcols(ensembl_gr)$name %in% cancer_gene_census$Gene.Symbol
  ]
  cdkn2ab_range <- GenomicRanges::GRanges(
    seqnames = "chr9",
    ranges = if (identical(ref_genome, "hg19")) {
      IRanges::IRanges(start = 21967751, end = 22009304)
    } else {
      IRanges::IRanges(start = 21967752, end = 21995301)
    },
    strand = "*"
  )
  S4Vectors::mcols(cdkn2ab_range)$name <- "CDKN2A/B"
  S4Vectors::mcols(cgc_gr) <- S4Vectors::mcols(cgc_gr)[, "name", drop = FALSE]
  cgc_gr <- GenomeInfoDb::keepStandardChromosomes(cgc_gr, pruning.mode = "coarse")
  GenomeInfoDb::seqlevelsStyle(cgc_gr) <- "UCSC"   # match conumee's "chr1" naming
  suppressWarnings(c(cgc_gr, cdkn2ab_range))
}

cnv_call_focal_events <- function(segment_df, detail_df, resources, sample_id) {
  cgc <- resources$cgc
  detail_regions <- make_cgc_segment_genes(resources$ensembl, cgc, resources$genome)
  
  segment_df <- cnv_harmonise_segments(segment_df)
  # Autosomal calls only.
  segment_df <- segment_df[segment_df$chrom %in% paste0("chr", 1:22), , drop = FALSE]
  segment_df$ID <- sample_id
  segment_df$V10 <- NA_character_
  
  for (i in seq_len(nrow(segment_df))) {
    event <- GenomicRanges::GRanges(
      seqnames = as.character(segment_df$chrom[i]),
      ranges = IRanges::IRanges(
        start = as.integer(segment_df$loc.start[i]),
        end = as.integer(segment_df$loc.end[i])
      )
    )
    hits <- GenomicRanges::findOverlaps(detail_regions, event)
    hit_gr <- detail_regions[S4Vectors::queryHits(hits)]
    nm <- S4Vectors::mcols(hit_gr)$name
    nm <- nm[nm %in% c(cgc$Gene.Symbol, "CDKN2A/B")]
    segment_df$V10[i] <- paste(unique(nm), collapse = "; ")
  }
  
  segment_anno <- segment_df %>%
    dplyr::mutate(V10 = strsplit(as.character(V10), ";\\s*")) %>%
    tidyr::unnest(cols = "V10") %>%
    dplyr::filter(!is.na(V10), V10 != "", !is.na(seg.mean), seg.mean != 0) %>%
    dplyr::mutate(
      CNV = ifelse(loc.end - loc.start < 4000000, "focal", "nonfocal")
    )
  
  detail_df <- strip_conumee_sample_prefix(detail_df)
  required_detail <- c("Chromosome", "Start", "End", "Name", "Value")
  if (!all(required_detail %in% names(detail_df))) {
    stop(
      "Conumee detail output is missing required columns: ",
      paste(setdiff(required_detail, names(detail_df)), collapse = ", "),
      call. = FALSE
    )
  }
  
  detail_as_segments <- detail_df %>%
    dplyr::transmute(
      chrom = as.character(Chromosome),
      loc.start = as.integer(Start),
      loc.end = as.integer(End),
      V10 = as.character(Name),
      ID = sample_id,
      num.mark = if ("num.mark" %in% names(detail_df)) num.mark else NA_integer_,
      bstat = if ("bstat" %in% names(detail_df)) bstat else NA_real_,
      pval = if ("pval" %in% names(detail_df)) pval else NA_real_,
      seg.mean = as.numeric(Value),
      seg.median = if ("seg.median" %in% names(detail_df)) seg.median else NA_real_,
      CNV = "focal"
    )
  
  wanted <- names(detail_as_segments)
  for (cc in setdiff(wanted, names(segment_anno))) segment_anno[[cc]] <- NA
  for (cc in setdiff(names(segment_anno), wanted)) detail_as_segments[[cc]] <- NA
  segment_anno <- segment_anno[, wanted, drop = FALSE]
  detail_as_segments <- detail_as_segments[, wanted, drop = FALSE]
  
  segment_anno <- segment_anno %>%
    dplyr::mutate(CNV_flt = paste(ID, V10, chrom, loc.start, loc.end, CNV, sep = "_"))
  detail_as_segments <- detail_as_segments %>%
    dplyr::mutate(CNV_flt = paste(ID, V10, chrom, loc.start, loc.end, CNV, sep = "_"))
  
  segment_unique <- segment_anno %>%
    dplyr::filter(!(CNV_flt %in% detail_as_segments$CNV_flt))
  
  amp_matrix <- dplyr::bind_rows(
    segment_unique %>% dplyr::select(dplyr::all_of(wanted)),
    detail_as_segments %>% dplyr::select(dplyr::all_of(wanted))
  )
  
  if (nrow(amp_matrix) == 0) {
    empty <- data.frame(
      ID = character(), V10 = character(), focal = numeric(), nonfocal = numeric(),
      CNV_state = character(), CNV_call = character(), stringsAsFactors = FALSE
    )
    return(list(all = empty, diagnostic = empty))
  }
  
  result <- amp_matrix %>%
    dplyr::group_by(ID, V10, CNV) %>%
    dplyr::summarise(seg.mean = mean(seg.mean, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = CNV, values_from = seg.mean, values_fill = 0)
  
  if (!"focal" %in% names(result)) result$focal <- 0
  if (!"nonfocal" %in% names(result)) result$nonfocal <- 0
  
  result <- result %>%
    dplyr::mutate(
      focal = as.numeric(focal),
      nonfocal = as.numeric(nonfocal),
      has_nonfocal = !is.na(nonfocal) & nonfocal != 0,
      delta = ifelse(has_nonfocal, focal - nonfocal, NA_real_),
      CNV_state = dplyr::case_when(
        is.na(focal) ~ NA_character_,
        focal < CNV_HOMDEL_THRESHOLD ~ "HOMDEL",
        focal >= CNV_AMP_THRESHOLD ~ "AMP",
        focal <= -0.17 & (!has_nonfocal | nonfocal <= 0 | delta < -0.17) ~ "LOSS",
        focal >= 0.15 & (!has_nonfocal | nonfocal >= 0 | delta > 0.17) ~ "GAIN",
        TRUE ~ "BALANCED"
      ),
      CNV_call = paste(CNV_state, round(focal, 3), sep = "_")
    )
  
  diagnostic <- result %>%
    dplyr::filter(!is.na(CNV_state), CNV_state %in% c("HOMDEL", "AMP"))
  
  list(all = result, diagnostic = diagnostic)
}

run_cnv_profile <- function(prefix, mset, sample_id, platform, bin_size = 50000L,
                            progress = cnv_no_progress) {
  progress(0.02, "Checking required CNV packages…")
  cnv_require_packages()
  cfg <- cnv_platform_config(platform)
  
  if (!is.null(mset) && ncol(mset) != 1) mset <- mset[, 1]
  
  # CNV query: SeSAMe (QCDPB) total intensities from the IDATs on the sample's
  # own array, processed exactly as the release references.
  red <- paste0(prefix, "_Red.idat")
  grn <- paste0(prefix, "_Grn.idat")
  if (!file.exists(red) || !file.exists(grn)) {
    stop("The CNV workflow cannot find the matching Red/Green IDAT pair for prefix: ", prefix)
  }
  # Probe address/mask data needed by readIDATpair and the Q (quality-mask) step.
  invisible(try(sesameData::sesameDataCache(paste0(cfg$sesame, ".address")), silent = TRUE))
  progress(0.05, paste0("Reading ", platform, " IDAT intensities with SeSAMe (", cfg$sesame,
                        if (nzchar(CNV_SESAME_PREP)) paste0(", prep ", CNV_SESAME_PREP) else "", ")…"))
  sig <- sesame::readIDATpair(prefix.path = prefix, platform = cfg$sesame, verbose = FALSE)
  if (nzchar(CNV_SESAME_PREP)) sig <- sesame::prepSesame(sig, CNV_SESAME_PREP)
  
  progress(0.08, "Predicting sex (sesame inferSex)…")
  sex_call <- cnv_predict_sex_sesame(sig, mset)
  
  query_int <- sesame::totalIntensities(sig)
  query_mat <- matrix(query_int, ncol = 1, dimnames = list(names(query_int), "query"))
  query_mat <- collapse_v2_matrix(query_mat)   # EPICv2 replicate IDs -> cg/ch; no-op otherwise
  measured <- rownames(query_mat)[is.finite(query_mat[, 1]) & query_mat[, 1] > 0]
  
  resources <- cnv_load_resources(platform, progress)
  reference <- cnv_load_reference(platform, progress)
  ref_label <- sprintf("%s (%s F + %s M)", cnv_reference_label(platform),
                       attr(reference, "n_female"), attr(reference, "n_male"))
  
  progress(
    0.28,
    paste0(
      "Building ", platform, " conumee annotation (", cfg$genome, ", autosomes, ",
      format(as.integer(bin_size), big.mark = ","),
      " bp bins; slowest step on first run)…"
    )
  )
  anno <- cnv_build_annotation(platform, resources, bin_size, reference = reference)
  anno <- restrict_anno_to_probes(anno, measured, paste0("Sample (", cfg$source_label, ")"))
  n_bins <- tryCatch(length(anno@bins), error = function(e) NA_integer_)
  n_probes <- tryCatch(length(anno@probes), error = function(e) NA_integer_)
  if (!is.na(n_bins)) {
    progress(
      0.45,
      paste0(
        "Annotation ready: ", format(n_probes, big.mark = ","), " probes in ",
        format(n_bins, big.mark = ","), " bins"
      )
    )
  }
  
  query <- conumee2::CNV.load(as.data.frame(query_mat[names(anno@probes), , drop = FALSE]))
  # Use a syntactic internal name so conumee's column prefixes are predictable
  # (sample IDs like "GSM9238338 iPSC" contain spaces). The display ID is kept
  # in sample_id and used for plots, tables and downloads.
  conumee_id <- gsub("[^A-Za-z0-9_]", "_", sample_id)
  if (!grepl("^[A-Za-z]", conumee_id)) conumee_id <- paste0("S_", conumee_id)
  names(query) <- conumee_id
  
  progress(0.52, paste0("Fitting sample against the ", ref_label, " (CNV.fit)…"))
  cnv <- conumee2::CNV.fit(query = query, ref = reference, anno = anno)
  progress(0.60, "Summarising log2 ratios into genomic bins (CNV.bin)…")
  cnv <- conumee2::CNV.bin(cnv)
  progress(0.66, "Scoring Cancer Gene Census detail regions (CNV.detail)…")
  cnv <- conumee2::CNV.detail(cnv)
  progress(0.72, "Segmenting the copy-number profile (CNV.segment)…")
  cnv <- conumee2::CNV.segment(cnv)
  
  progress(0.76, "Writing segment and detail tables…")
  tmp_detail <- tempfile(fileext = ".detail.txt")
  tmp_segments <- tempfile(fileext = ".seg")  # conumee expects a .seg extension
  on.exit(unlink(c(tmp_detail, tmp_segments)), add = TRUE)
  
  conumee2::CNV.write(cnv, what = "detail", file = tmp_detail)
  conumee2::CNV.write(cnv, what = "segments", file = tmp_segments)
  
  # Optional debugging: uncomment to see conumee's raw column names.
  # message("detail cols: ", paste(names(cnv_read_written_table(tmp_detail)), collapse = ", "))
  # message("segment cols: ", paste(names(cnv_read_written_table(tmp_segments)), collapse = ", "))
  
  detail_df <- cnv_read_written_table(tmp_detail)
  segment_df <- cnv_read_written_table(tmp_segments)
  progress(
    0.80,
    paste0("Calling focal AMP/HOMDEL events across ", format(nrow(segment_df), big.mark = ","), " segments…")
  )
  focal <- cnv_call_focal_events(segment_df, detail_df, resources, sample_id)
  diagnostic <- as.data.frame(focal$diagnostic)
  
  segments_clean <- cnv_harmonise_segments(segment_df)
  progress(0.83, paste0("Predicted sex: ", sex_call$sex, " (", sex_call$note, ")"))
  
  # Per-sample summary strings, as in diagnostics_results_<date>.txt
  # (e.g. "EGFR_AMP_1.234|MYC_AMP_0.512").
  summarise_calls <- function(state) {
    d <- diagnostic[diagnostic$CNV_state == state, , drop = FALSE]
    if (nrow(d) == 0) return("none")
    paste(paste(d$V10, d$CNV_call, sep = "_"), collapse = "|")
  }
  amps <- summarise_calls("AMP")
  dels <- summarise_calls("HOMDEL")
  
  # Plotting fit (Step 6 of the standalone script), with detail regions =
  # CGC genes that were called (AMP/HOMDEL) + StemCNV-check hotspots that are
  # HIT in this sample (overlapping segment in the curated direction, or an
  # overlapping focal call). Hotspots are annotation only: they are not used
  # for the high-confidence focal calls above.
  hit_genes <- collapse_cdkn2ab_hits(diagnostic$V10)
  hit_genes <- hit_genes[!is.na(hit_genes) & hit_genes != ""]
  cnv_plot <- cnv
  plot_detail <- FALSE
  
  detail_all <- make_cgc_detail_regions(resources$ensembl, resources$cgc, cfg$genome)
  detail_hits <- detail_all[S4Vectors::mcols(detail_all)$name %in% hit_genes]
  S4Vectors::mcols(detail_hits) <- S4Vectors::mcols(detail_hits)[, "name", drop = FALSE]
  
  # Focal-call regions with labels such as "KRAS AMP" for hotspot overlap.
  focal_regions <- detail_hits
  if (length(focal_regions)) {
    st <- diagnostic$CNV_state[match(S4Vectors::mcols(focal_regions)$name, diagnostic$V10)]
    st[is.na(st)] <- diagnostic$CNV_state[diagnostic$V10 %in% c("CDKN2A", "CDKN2B", "CDKN2A/B")][1]
    S4Vectors::mcols(focal_regions)$label <- paste(S4Vectors::mcols(focal_regions)$name, st)
  }
  
  hotspot_hits <- find_hotspot_hits(resources$hotspots, segments_clean, focal_regions)
  hotspot_summary <- if (length(hotspot_hits)) {
    paste(paste0(S4Vectors::mcols(hotspot_hits)$name, " [",
                 S4Vectors::mcols(hotspot_hits)$call_type, "; ",
                 S4Vectors::mcols(hotspot_hits)$evidence, "]"), collapse = " | ")
  } else {
    "none"
  }
  
  alterations <- build_alteration_table(
    diagnostic   = diagnostic,
    gene_regions = make_cgc_segment_genes(resources$ensembl, resources$cgc, cfg$genome),
    hotspots_all = resources$hotspots,
    segments     = segments_clean,
    cytobands    = resources$cytobands
  )
  
  hotspots <- hotspot_hits
  if (length(hotspots)) {
    # A gene that is both a CGC hit and a hotspot (e.g. KRAS, TP53) is labelled once.
    hotspots <- hotspots[!S4Vectors::mcols(hotspots)$name %in% S4Vectors::mcols(detail_hits)$name]
    S4Vectors::mcols(hotspots) <- S4Vectors::mcols(hotspots)[, "name", drop = FALSE]
  }
  # Plot labels (drawn by cnv_add_plot_labels(), never overlapping):
  # CGC genes with a focal call + StemCNV-check hotspots hit in this sample.
  # Each label is anchored at the gene/region midpoint and its log2 value.
  focal_value <- function(g) {
    v <- as.numeric(diagnostic$focal[diagnostic$V10 == g])
    if (!length(v) && g == "CDKN2A/B") {
      v <- as.numeric(diagnostic$focal[diagnostic$V10 %in% c("CDKN2A", "CDKN2B")])
    }
    if (length(v)) v[which.max(abs(v))] else NA_real_
  }
  label_rows <- list()
  if (length(detail_hits)) {
    label_rows[[1]] <- data.frame(
      label = S4Vectors::mcols(detail_hits)$name,
      chr   = as.character(GenomicRanges::seqnames(detail_hits)),
      pos   = (GenomicRanges::start(detail_hits) + GenomicRanges::end(detail_hits)) / 2,
      value = vapply(S4Vectors::mcols(detail_hits)$name, focal_value, numeric(1)),
      type  = "focal",
      stringsAsFactors = FALSE
    )
  }
  if (length(hotspots)) {
    hh <- hotspot_hits[match(S4Vectors::mcols(hotspots)$name, S4Vectors::mcols(hotspot_hits)$name)]
    label_rows[[length(label_rows) + 1]] <- data.frame(
      label = S4Vectors::mcols(hh)$name,
      chr   = as.character(GenomicRanges::seqnames(hh)),
      pos   = (GenomicRanges::start(hh) + GenomicRanges::end(hh)) / 2,
      value = S4Vectors::mcols(hh)$hit_value,
      type  = "hotspot",
      stringsAsFactors = FALSE
    )
  }
  plot_labels <- if (length(label_rows)) do.call(rbind, label_rows) else
    data.frame(label = character(), chr = character(), pos = numeric(),
               value = numeric(), type = character(), stringsAsFactors = FALSE)
  plot_labels$value[!is.finite(plot_labels$value)] <- 0
  plot_detail <- nrow(plot_labels) > 0
  
  plot_main <- if (plot_detail && length(hit_genes) && length(hotspots)) {
    paste0(sample_id, " - focal hits + stem cell hotspot hits (", cfg$source_label, ", ", ref_label,
           "; inferred sex: ", sex_call$sex, ")")
  } else if (plot_detail && length(hit_genes)) {
    paste0(sample_id, " - focal hits only (", cfg$source_label, ", ", ref_label,
           "; inferred sex: ", sex_call$sex, ")")
  } else if (plot_detail) {
    paste0(sample_id, " - stem cell hotspot hits; no focal hits (", cfg$source_label, ", ", ref_label,
           "; inferred sex: ", sex_call$sex, ")")
  } else {
    paste0(sample_id, " (", cfg$source_label, ", ", ref_label,
           "; inferred sex: ", sex_call$sex, ") - no focal hits")
  }
  
  progress(
    0.98,
    paste0(
      "Done: ", format(nrow(segment_df), big.mark = ","), " segments, ",
      nrow(diagnostic), " high-confidence focal call(s)"
    )
  )
  
  list(
    cnv = cnv,
    cnv_plot = cnv_plot,
    plot_labels = plot_labels,
    plot_detail = plot_detail,
    plot_main = plot_main,
    segments = order_by_genome(segments_clean, "chrom", "loc.start"),
    detail = strip_conumee_sample_prefix(detail_df),
    focal_diagnostic = diagnostic,
    amps = amps,
    hotspots_hit = hotspot_summary,
    alterations = alterations,
    dels = dels,
    sample_id = sample_id,
    platform = platform,
    sex = CNV_REFERENCE_SEX,
    reference_label = ref_label,
    reference_note = attr(reference, "excluded_note"),
    inferred_sex = sex_call$sex,
    sex_diff = sex_call$diff,
    sex_note = sex_call$note,
    genome = cfg$genome,
    bin_size = as.integer(bin_size)
  )
}

# Draw the genome plot (autosomes only).
cnv_draw_genomeplot <- function(res, detail = res$plot_detail) {
  # conumee2 draws the profile WITHOUT its own detail labels; the labels are
  # added by cnv_add_plot_labels() so that they never overlap. The y-range is
  # widened to leave label rows above and below the profile.
  args <- list(res$cnv_plot, chr = CNV_PLOT_CHRS, cols = CNV_PLOT_COLS,
               main = res$plot_main, detail = FALSE, ylim = CNV_PLOT_YLIM)
  drawn <- tryCatch({ do.call(conumee2::CNV.genomeplot, args); TRUE },
                    error = function(e) {
                      if (!grepl("unused argument", conditionMessage(e))) stop(e)
                      FALSE
                    })
  if (!drawn) {                       # conumee2 version without 'ylim'
    args$ylim <- NULL
    do.call(conumee2::CNV.genomeplot, args)
  }
  if (isTRUE(detail) && !is.null(res$plot_labels) && nrow(res$plot_labels)) {
    tryCatch(
      cnv_add_plot_labels(res$cnv_plot, res$plot_labels),
      error = function(e) message("[CNV] Could not add plot labels: ", conditionMessage(e))
    )
  }
}

# Non-overlapping labels on a conumee2 genome plot.
#  x: conumee lays chromosomes end to end: x = cumulative size of the preceding
#     plotted chromosomes + position (same as CNV.genomeplot).
#  y: gains are labelled in rows above the profile, losses in rows below; each
#     label goes into the first row where it does not touch its neighbours, and a
#     leader line joins it to a marker at its log2 value.
cnv_add_plot_labels <- function(cnv, labels, cex = 0.62) {
  genome <- cnv@anno@genome
  chrs <- CNV_PLOT_CHRS[CNV_PLOT_CHRS %in% rownames(genome)]
  sizes <- as.numeric(genome[chrs, "size"])
  offset <- stats::setNames(c(0, cumsum(sizes))[seq_along(chrs)], chrs)
  
  lab <- labels[labels$chr %in% chrs, , drop = FALSE]
  if (!nrow(lab)) return(invisible(NULL))
  lab$x <- offset[lab$chr] + lab$pos
  lab$up <- lab$value >= 0
  
  usr <- graphics::par("usr")
  row_h <- graphics::strheight("Xg", cex = cex) * 1.45
  pad <- graphics::strwidth("  ", cex = cex)
  lab$w <- graphics::strwidth(lab$label, cex = cex) + pad
  
  # Label bands: just outside the default profile range (+/-1.25).
  band_start <- c(up = min(1.32, usr[4] - row_h), down = max(-1.32, usr[3] + row_h))
  n_rows <- c(
    up   = max(1, floor((usr[4] - row_h / 2 - band_start[["up"]]) / row_h) + 1),
    down = max(1, floor((band_start[["down"]] - usr[3] - row_h / 2) / row_h) + 1)
  )
  
  lab$lx <- NA_real_; lab$ly <- NA_real_
  for (side in c("up", "down")) {
    idx <- which(lab$up == (side == "up"))
    idx <- idx[order(lab$x[idx])]
    right_edge <- rep(-Inf, n_rows[[side]])
    for (k in idx) {
      lx <- min(max(lab$x[k], usr[1] + lab$w[k] / 2), usr[2] - lab$w[k] / 2)
      row <- which(right_edge <= lx - lab$w[k] / 2)[1]
      if (is.na(row)) {
        # All rows occupied here: shift right past the row that frees up first.
        row <- which.min(right_edge)
        lx <- min(right_edge[row] + lab$w[k] / 2, usr[2] - lab$w[k] / 2)
      }
      right_edge[row] <- lx + lab$w[k] / 2
      lab$lx[k] <- lx
      lab$ly[k] <- if (side == "up") band_start[["up"]] + (row - 1) * row_h
      else band_start[["down"]] - (row - 1) * row_h
    }
  }
  
  col_txt <- ifelse(lab$type == "hotspot", "#1B7837", "black")
  y_val <- pmin(pmax(lab$value, usr[3]), usr[4])
  graphics::segments(lab$x, y_val, lab$lx, lab$ly + ifelse(lab$up, -row_h / 2, row_h / 2),
                     col = "grey55", lwd = 0.6)
  graphics::points(lab$x, y_val, pch = 19, cex = 0.45, col = col_txt)
  graphics::rect(lab$lx - lab$w / 2, lab$ly - row_h / 2, lab$lx + lab$w / 2, lab$ly + row_h / 2,
                 col = grDevices::adjustcolor("white", alpha.f = 0.85), border = NA)
  graphics::text(lab$lx, lab$ly, lab$label, cex = cex, col = col_txt)
  
  graphics::legend("bottomright", bty = "n", cex = 0.65, text.col = c("black", "#1B7837"),
                   legend = c("focal call (CGC)", "stem cell hotspot"), inset = c(0, 0.005))
  invisible(lab)
}

# Activate thematic for ggplot2 styling based on bs_theme
thematic::thematic_shiny()

# Define a Bootstrap theme using bs_theme
theme <- bs_theme(
  bootswatch = "flatly",
  primary = "#0072B2",        # Set primary color to match buttons
  base_font = font_google("Roboto")
)

# Custom CSS to style the sidebar panel and navbar
custom_css <- "
/* Style the sidebar panel */
.sidebarPanel {
  background-color: #f8f9fa !important;
}

/* Style the navbar */
.navbar.navbar-default {
  background-color: #0072B2 !important;
  border-color: #005c99 !important;
}
.navbar.navbar-default .navbar-brand {
  color: #FFFFFF !important;
}
.navbar.navbar-default .navbar-nav > li > a {
  color: #FFFFFF !important;
}
.navbar.navbar-default .navbar-nav > li > a:hover {
  background-color: #005c99 !important;
  color: #FFFFFF !important;
}
.btn {
  background-color: #0072B2 !important;
  color: #FFFFFF !important;
  border: none;
}
.btn:hover {
  background-color: #005c99 !important;
}
.btn-green {
  background-color: #28a745;
  color: white;
  border-color: #28a745;
  border-radius: .25rem;   
}
.progress-bar {
  background-color: #28a745 !important;
}
.progress-bar-danger {
  background-color: #D55E00 !important;
}
#form .action-button {
  display: block;
  width: 100%;       /* optional, makes them full-width */
  margin-bottom: 10px;
  
}
a {
  color: #2a9fd6 !important;
}
a:hover {
  color: #005c99 !important;
}
/* Animate the instructions panel on entry */
@keyframes fadeSlideIn {
  from {
    opacity: 0;
    transform: translateY(20px);
  }
  to {
    opacity: 1;
    transform: translateY(0);
  }
}

.instructions-panel {
  animation: fadeSlideIn 0.6s ease-out;
}

/* Animate each list item in sequence */
.instructions-panel ol li {
  opacity: 0;
  animation: fadeSlideIn 0.5s ease-out forwards;
}
.instructions-panel ol li:nth-child(1) { animation-delay: 0.1s; }
.instructions-panel ol li:nth-child(2) { animation-delay: 0.2s; }
.instructions-panel ol li:nth-child(3) { animation-delay: 0.3s; }
.instructions-panel ol li:nth-child(4) { animation-delay: 0.4s; }
.instructions-panel ol li:nth-child(5) { animation-delay: 0.5s; }
.instructions-panel ol li:nth-child(6) { animation-delay: 0.6s; }
.instructions-panel ol li:nth-child(7) { animation-delay: 0.7s; }
.instructions-panel ol li:nth-child(8) { animation-delay: 0.8s; }
"


# Define the Shiny UI with the selected theme and custom CSS
ui <- navbarPage(
  title = "SteMClass",
  theme = theme,
  # Non-tab content belongs in `header` (silences the navbarPage warning).
  header = tagList(
    useShinyjs(),
    tags$head(
      tags$style(HTML(custom_css))
    )
  ),
  tabPanel(
    title = "Classification",
    sidebarLayout(
      sidebarPanel(
        fileInput(
          inputId = "idat_files",
          label   = tags$span(icon("upload"), " Upload IDAT Files"),
          multiple = TRUE,
          accept   = ".idat"
        ),
        
        tags$div(
          id = "form",
          
          textInput(
            inputId = "sample_accession",
            label   = tags$span(icon("tag"), " Enter Sample ID:"),
          ),
          
          selectInput(
            inputId = "array_version",
            label   = tags$span(icon("microchip"), " Array Version:"),
            choices = c("450K", "EPICv1", "EPICv2")
          ),
          tags$div(
            style = "text-align:center; color:#888; margin:8px 0;",
            tags$span("- OR -")
          ),
          
          ui_example_block <- tagList(
            
            
            checkboxInput(
              inputId = "use_example",
              label   = tags$span(icon("file-import"), " Import an example from GEO Repository"),
              value   = FALSE
            ),
            
            
            # Collapsible panel – only visible when checkbox is ticked
            conditionalPanel(
              condition = "input.use_example == true",
              
              wellPanel(
                style = "background:#f8f9fa ; border:1px solid #0072B2; padding:10px;",
                
                selectInput(
                  inputId  = "example_sample",
                  label    = "Choose an example sample:",
                  choices  = example_choices,
                  selected = ""
                ),
                
                # Small status line
                uiOutput("example_status")
              )
            ),
            
            hr()
          )
        ),
        
        
        actionButton(
          inputId = "load_samples",
          label   = "Preprocess Data",
          icon    = icon("broom"),
          style   = "margin-bottom: 10px;"
        ),
        hr(),
        
        actionButton(
          inputId = "classify",
          label   = "Run Classification",
          icon    = icon("play"),
          style   = "margin-bottom: 10px;"
        ),
        hr(),
        
        actionButton(
          inputId = "reset",
          label   = "New Analysis",
          icon    = icon("redo")
        ),
        
        hr(),
        h4("About SteMClass:"),
        tags$p(
          "SteMClass is a DNA methylation-based classifier for identifying the differentiation 
  state of iPSC-derived samples from a single array-based assay, compatible with all 
  Illumina methylation array versions (450K, EPICv1, EPICv2)."
        ),
        tags$p(
          "Trained on a curated reference cohort spanning 15 iPSC lines and seven differentiation 
  states, the classifier achieves high accuracy on both internal validation and independent 
  external datasets."
        ),
        tags$p(
          "For full methodological details and performance metrics, see the ",
          tags$a(href = "https://doi.org/10.1101/2025.09.02.673063", target = "_blank", "accompanying publication"),
          "."
        ) 
      ),  
      
      mainPanel(
        # Instructions panel - visible until results are ready
        conditionalPanel(
          condition = "!output.results_ready",
          
          div(
            class = "well instructions-panel",
            style = "max-width:800px; margin:40px auto; padding:30px; background:#f9f9f9; border-radius:8px;",
            
            h3(icon("circle-info"), " How to use SteMClass",
               style = "margin-top:0; margin-bottom:20px;"),
            
            tags$ol(
              style = "line-height: 2;",
              tags$li(
                icon("upload"), " ",
                tags$b("Upload"), " the sample's pair of IDAT files (Red and Green)."
              ),
              tags$li(
                icon("tag"), " Enter the",
                tags$b(" Sample ID")
              ),
              tags$li(
                icon("microchip")," Select the", 
                tags$b(" Array version"), " ", " (450K, EPICv1, or EPICv2)."
              ),
              tags$li(
                icon("file-import"), " ",
                tags$b(""), " Alternatively, click ", tags$b("\"Import an example from GEO Repository\""),
                " to load publicly available test data ",
                tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308134", target = "_blank", "(GSE308134).")
                
              ),
              tags$li(
                icon("broom"), " Click ", tags$b("\"Preprocess Data\""), " to prepare the sample."
              ),
              tags$li(
                icon("play"), " Click ", tags$b("\"Run Classification\""), " to classify the sample."
              ),
              tags$li(
                icon("magnifying-glass-chart"), "After classification, explore the data visualisation tabs to inspect the sample's DNA methylation profile (heatmaps) and its position relative to the training data (UMAP)."
              ),
              
              tags$li(
                icon("redo"), " Click ", tags$b("\"New Analysis\""), " to start over."
              )
            )
          )
        ),
        
        # Results panel - visible once results are ready
        conditionalPanel(
          condition = "output.results_ready",
          div(style = "position: relative; padding-bottom: 60px;",
              uiOutput("classification_ui"))
        )
      )
    )
  ),    
  tabPanel(
    title = "Sample Visualization",
    sidebarLayout(
      sidebarPanel(
        h4("UMAP"),
        actionButton("generate_umap", "Generate UMAP", icon = icon("chart-area")),
        hr(),
        h4("Instructions:"),
        tags$ol(
          tags$li('After classifying a sample, click "Generate UMAP" to visualize its relationship to the training data.'),
          tags$li("The UMAP will use the CpG sites from the classifier.")
        )
      ),
      mainPanel(
        div(
          style = "position: relative; padding-bottom: 60px;",
          
          h3("UMAP Visualization"),
          plotOutput("umap_plot", height = "600px")
        )
      )
    )
  ),
  tabPanel(
    title = "Cell State Global Beta Value Visualization",
    sidebarLayout(
      sidebarPanel(
        h4("Global Heatmap"),
        selectInput("marker", "Select Reference Cell State:", choices = c( "Ectoderm", "Endoderm", "Mesoderm", "Lung", "Endothelial", "NSC", "Astrocyte")),
        actionButton("marker_plot_button", "Plot Heatmap"),
        hr(),
        h4("Instructions:"),
        tags$ol(
          tags$li('Generate a heatmap that visualizes the global DNA methylation profiles of the test sample and reference set.'),
          tags$li("The heatmap will use the top 10,000 CpG sites that distinguish the selected cell state from iPSC.")
        )
      ),
      
      mainPanel(
        div(
          style = "position: relative; padding-bottom: 60px;",
          
          withSpinner(
            plotOutput("marker_heatmap_plot"),
            type  = 5,
            color = "#0072B2"
          )
        ),
        
        # download button
        div(style = "position: fixed; bottom: 10px; right: 10px; z-index: 1000;",
            downloadButton(
              "download_marker_heatmap", "Download Marker Heatmap",
              onclick = "$('#download_spinner').show();"
            ),
            tags$i(
              id    = "download_spinner",
              class = "fa fa-spinner fa-spin fa-lg",
              style = "display: none; margin-left: 8px; color: #0072B2;"
            ))
      )
    )
  ),
  tabPanel(
    title = "Cell State Gene Level Beta Value Visualization",
    sidebarLayout(
      sidebarPanel(
        h4("Gene Heatmap"),
        selectInput(
          inputId = "celltype",
          label   = "Select Reference Cell State:",
          choices = c(
            "Ectoderm",
            "Endoderm",
            "Mesoderm",
            "Endothelial",
            "Lung",
            "NSC",
            "Astrocyte"
          ),
          selected = "Ectoderm"
        ),
        
        
        textInput(
          inputId = "gene_input",
          label   = "Enter Gene Name:",
          value   = "POU5F1"
        ),
        
        actionButton(
          inputId = "plot_button",
          label   = "Plot Heatmap"
        ),
        hr(),
        h4("Instructions:"),
        tags$ol(
          tags$li('Generate a heatmap that visualizes the gene-level DNA methylation profiles of the test sample and reference set.'),
          tags$li("The heatmap will use the CpG sites that are annotated to the given gene according to the Illumina 450K array annotation.")
        )
      ),
      
      mainPanel(
        div(
          style = "position: relative; padding-bottom: 60px;",
          
          withSpinner(
            plotOutput("heatmap_plot"),
            type  = 5,
            color = "#0072B2"
          )
        ),
        
        div(style = "position: fixed; bottom: 10px; right: 10px; z-index: 1000;",
            downloadButton("download_gene_heatmap", "Download Gene Heatmap"),
            onclick = "$('#download_spinner_gene').show();"
        ),
        tags$i(
          id    = "download_spinner_gene",
          class = "fa fa-spinner fa-spin fa-lg",
          style = "display: none; margin-left: 8px; color: #0072B2;"
        )
      )
    )
  ),
  tabPanel(
    title = "CNV Profiling",
    sidebarLayout(
      sidebarPanel(
        h4("Copy Number Profiling"),
        tags$p(
          "CNV profiling reads the uploaded IDAT pair on its own array (450K, EPICv1 or EPICv2) with SeSAMe (QCDPB) ",
          "and fits it with conumee2 against that array's female and male control references, pooled into one mixed-sex reference. ",
          "Copy number is reported for the autosomes; sex is predicted separately with SeSAMe inferSex."
        ),
        numericInput(
          inputId = "cnv_bin_size",
          label = "Minimum CNV bin size (bp):",
          value = 50000,
          min = 10000,
          step = 10000
        ),
        actionButton(
          inputId = "run_cnv",
          label = "Run CNV Profiling",
          icon = icon("chart-line")
        ),
        hr(),
        h4("Instructions:"),
        tags$ol(
          tags$li("Upload the Red/Green IDAT pair and enter the sample ID on the Classification tab."),
          tags$li("Select the correct array version on the Classification tab."),
          tags$li('Click "Run CNV Profiling".'),
          tags$li("Inspect the autosomal profile, the predicted sex and high-confidence focal AMP/HOMDEL calls."),
          tags$li("chrX and chrY are not shown.")
        ),
        hr(),
        uiOutput("cnv_input_summary")
      ),
      mainPanel(
        conditionalPanel(
          condition = "!output.cnv_results_ready",
          div(
            class = "well instructions-panel",
            style = "max-width:850px; margin:40px auto; padding:30px; background:#f9f9f9; border-radius:8px;",
            h3(icon("chart-line"), " CNV Profiling"),
            tags$p(
              "No CNV result has been generated for this session. The analysis runs independently of the methylation classifier and starts from the uploaded IDAT intensities."
            )
          )
        ),
        conditionalPanel(
          condition = "output.cnv_results_ready",
          uiOutput("cnv_result_header"),
          withSpinner(
            plotOutput("cnv_genome_plot", height = "650px"),
            type = 5,
            color = "#0072B2"
          ),
          fluidRow(
            column(3, downloadButton("download_cnv_plot", "Download CNV Plot")),
            column(3, downloadButton("download_cnv_plot_unlabelled", "Download Unlabelled Plot")),
            column(3, downloadButton("download_cnv_segments", "Download Segments")),
            column(3, downloadButton("download_cnv_focal", "Download Alterations"))
          ),
          hr(),
          h3("Detected alterations"),
          tags$p(
            "All autosomal segments with a gain (log2 ≥ ", CNV_HOTSPOT_GAIN, ") or loss (log2 ≤ ", CNV_HOTSPOT_LOSS,
            "), labelled AMP / HOMDEL beyond ", CNV_AMP_THRESHOLD, " / ", CNV_HOMDEL_THRESHOLD,
            ", and all high-confidence focal amplifications and homozygous deletions in Cancer Gene Census genes.The stem cell hotspot column lists known iPSC hotspots hit by each alteration ",
            "Coordinates are ", tags$code("hg19"), " for 450K / EPICv1 and ",
            tags$code("hg38"), " for EPICv2."
          ),
          DT::DTOutput("cnv_focal_table")
        )
      )
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  # ── Reactive: look up full metadata for the chosen example ──────────────
  selected_example <- reactive({
    req(input$use_example, input$example_sample != "")
    acc <- input$example_sample
    Find(function(s) s$accession == acc, example_samples)
  })
  
  # ── Function: download both IDATs to tempdir, return their paths ─────────
  download_example_idats <- function(ex) {
    tmp <- tempdir()
    grn_path <- file.path(tmp, paste0(ex$accession, "_Grn.idat"))
    red_path <- file.path(tmp, paste0(ex$accession, "_Red.idat"))
    
    withProgress(message = "Downloading example IDAT files", value = 0, {
      incProgress(0.5, detail = "Downloading Grn channel…")
      download.file(ex$grn_url, grn_path, mode = "wb", quiet = TRUE)
      
      incProgress(0.5, detail = "Downloading Red channel…")
      download.file(ex$red_url, red_path, mode = "wb", quiet = TRUE)
    })
    
    list(grn = grn_path, red = red_path)
  }
  
  # ── Auto-populate accession + array version when example is selected ─────
  observeEvent(selected_example(), {
    ex <- selected_example()
    req(ex)
    updateTextInput(session,   "sample_accession", value = ex$accession)
    updateSelectInput(session, "array_version",    selected = ex$array_version)
  })
  
  # ── Clear fields when checkbox is unticked ───────────────────────────────
  observeEvent(input$use_example, {
    if (!input$use_example) {
      updateSelectInput(session, "example_sample",   selected = "")
      updateTextInput(session,   "sample_accession", value = "")
    }
  })
  
  # ── Status label below the dropdown ──────────────────────────────────────
  output$example_status <- renderUI({
    ex <- selected_example()
    if (is.null(ex)) return(NULL)
    tags$small(
      style = "color:#555;",
      icon("info-circle"),
      sprintf(" %s %s IDATs will be imported automatically.",
              ex$array_version, ex$accession)
    )
  })
  
  # ── Shared raw-IDAT reactive used by classification and CNV profiling ─────
  sample_files <- reactive({
    req(input$sample_accession, input$array_version)
    
    session_dir <- file.path(tempdir(), paste0("stemclass_", session$token))
    dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (isTRUE(input$use_example) && input$example_sample != "") {
      ex <- selected_example()
      req(ex)
      paths <- download_example_idats(ex)
      source_paths <- c(paths$grn, paths$red)
      dest_names <- basename(source_paths)
    } else {
      req(input$idat_files)
      validate(
        need(nrow(input$idat_files) == 2,
             "Please upload exactly two IDAT files (one Red and one Grn).")
      )
      source_paths <- input$idat_files$datapath
      dest_names <- input$idat_files$name
    }
    
    dest_paths <- file.path(session_dir, dest_names)
    copied <- file.copy(source_paths, dest_paths, overwrite = TRUE)
    validate(need(all(copied), "Could not copy the uploaded IDAT files into the analysis workspace."))
    
    red <- grep("_Red\\.idat$", dest_paths, value = TRUE, ignore.case = TRUE)
    grn <- grep("_Grn\\.idat$", dest_paths, value = TRUE, ignore.case = TRUE)
    validate(
      need(length(red) == 1, "Exactly one *_Red.idat file is required."),
      need(length(grn) == 1, "Exactly one *_Grn.idat file is required.")
    )
    
    red_prefix <- sub("_Red\\.idat$", "", red, ignore.case = TRUE)
    grn_prefix <- sub("_Grn\\.idat$", "", grn, ignore.case = TRUE)
    validate(
      need(identical(red_prefix, grn_prefix),
           "The Red and Green IDAT filenames do not form a matching pair.")
    )
    
    list(
      sample_id = input$sample_accession,
      platform = input$array_version,
      prefix = red_prefix,
      directory = session_dir,
      red = red,
      grn = grn
    )
  })
  
  # Shared preprocessing: returns the 450K beta matrix (classifier) and the
  # 450K MethylSet (CNV profiling) from one Noob run.
  preprocessed_data <- reactive({
    sf <- sample_files()
    req(sf)
    
    tmp <- sf$directory
    targets <- data.frame(
      Sample_accession = sf$sample_id,
      Array_version = sf$platform,
      Basename = sf$prefix,
      stringsAsFactors = FALSE
    )
    
    # ── Shared preprocessing (runs for both uploaded and example IDATs)
    withProgress(message = "Reading & preprocessing IDATs", value = 0, {
      incProgress(0.2, detail = "Loading raw data…")
      
      
      library(tidymodels)
      library(methylumi)
      library(wateRmelon)  
      library(randomForest)
      library(DT)
      library(ggplot2)
      library(tidyr)
      library(dplyr)        # For data manipulation
      library(uwot)         # For UMAP computation
      library(RColorBrewer) # For color palettes
      library(data.table)   # For efficient data reading
      library(tidymodels)
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      library(IlluminaHumanMethylationEPICmanifest)
      library(IlluminaHumanMethylationEPICv2manifest)
      library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      library(IlluminaHumanMethylation450kmanifest)
      library(preprocessCore)
      library(minfi)        # For processing IDAT files
      library(ranger)
      library(stringr)
      library(ComplexHeatmap)
      library(circlize)     # For colorRamp2 function used in heatmap
      library(cachem)
      library(purrr)
      library(stringr)
      library(glmnet)
      library(yardstick)
      library(glue)
      library(textshape)
      
      rgSet <- read.metharray.exp(base = tmp, targets = targets)
      sampleNames(rgSet) <- targets$Sample_accession
      
      # Preprocess (split out each branch so we can show progress)
      
      if (input$array_version == "EPICv2") {
        incProgress(0.3, detail = "Preprocessing (Noob)…")
        rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
        mSetv2 <- preprocessNoob(rgSet)
        mset_native <- mSetv2            # own-array Noob data for CNV / sex
        mSetv2   <- epicv2clean(mSetv2)
        
        incProgress(0.05, detail = "Converting to 450K…")
        mSetv1 <- readRDS("Mset1_EPICv2_sample_combinearrays.rds")
        mSet450_all <- combineArrays(mSetv2,mSetv1,
                                     outType = "IlluminaHumanMethylation450k")
        mSet <- mSet450_all[, sampleNames(mSetv2)]
        
        incProgress(0.1, detail = "Mapping to genome…")
        mSetSq <- mapToGenome(mSet)
        
        beta   <- getBeta(mSetSq)
        
      } else if (input$array_version == "EPICv1") {
        incProgress(0.3, detail = "Preprocessing (Noob)…")
        mSet <- preprocessNoob(rgSet)
        mset_native <- mSet              # own-array Noob data for CNV / sex
        
        incProgress(0.05, detail = "Converting to 450K…")
        mSet <- convertArray(mSet,
                             outType = "IlluminaHumanMethylation450k")
        
        incProgress(0.2, detail = "Mapping to genome…")
        mSetSq <- mapToGenome(mSet)
        
        beta   <- getBeta(mSetSq)
        
      } else if (input$array_version == "450K") {
        incProgress(0.3, detail = "Preprocessing (Noob)…")
        mSet <- preprocessNoob(rgSet)
        mset_native <- mSet              # own-array Noob data for CNV / sex
        incProgress(0.35, detail = "Mapping to genome…")
        mSetSq <- mapToGenome(mSet)
        beta   <- getBeta(mSetSq)
        
      } else {
        stop("Unknown array version: ", input$array_version)
      }
      
      # Finalise
      incProgress(0.1, detail = "Finalizing…")
      beta_mat <- t(beta)
      
      # return both views of the same preprocessed sample
      list(beta = beta_mat, mset = mSet, mset_native = mset_native)
    })
  })
  
  beta_data <- reactive({
    preprocessed_data()$beta
  })
  
  
  
  
  observeEvent(input$load_samples, {
    req(beta_data())
    
    updateActionButton(session, "load_samples",
                       label = "Ready!",
                       icon  = icon("check"))
    # swap CSS classes:
    shinyjs::removeClass("load_samples", "btn")
    shinyjs::addClass("load_samples",    "btn-green")
    # optionally disable it so they can’t click again
    #shinyjs::disable("load_samples")
  })
  
  
  # Reactive expression to wrap the train_data
  ref_beta <- reactive({
    con <- url(urls$train_gmset, "rb")
    load(con)
    bVals_train <- getBeta(gmset_train)
    bVals_train
  })
  
  
  training_data <- reactive({
    # with rownames = sample IDs and a column "Class"
    train_data <- data.table::fread(urls$train_data)
    train_data$Class <- factor(train_data$Class_rf)
    train_data <- train_data[, !names(train_data) %in% "Class_rf", with = FALSE]
    if (!"Class" %in% colnames(train_data)) {
      stop("Global train_data must contain a column named 'Class'")
    }
    train_data
  })
  
  # Reactive expression to wrap the marker data
  marker_data <- reactive({
    # with rownames = sample IDs and a column "Class"
    marker_data <- data.table::fread(urls$marker_data)
    marker_data$Class <- factor(marker_data$Class_rf)
    marker_data <- marker_data[, !names(marker_data) %in% "Class_rf", with = FALSE]
    if (!"Class" %in% colnames(marker_data)) {
      stop("Global marker_data must contain a column named 'Class'")
    }
    marker_data
  })
  
  
  
  # Reactive expression to wrap the CpG anno
  ann450K <- reactive({
    ann450K <- data.table::fread(urls$cpg_anno)
    ann450K <- ann450K %>% mutate(UCSC_RefGene_Name = toupper(UCSC_RefGene_Name))
    ann450K
  })
  
  # Reactive expression to wrap the marker set sample anno
  sample_anno <- reactive({
    sample_anno <- data.table::fread(urls$train_anno)
    sample_anno$Class <- as.factor(sample_anno$Class_rf)
    sample_anno
  })
  
  # Reactive expression to wrap the marker beta matrix
  marker_beta <- reactive({
    md <- marker_data()
    sa <- sample_anno()
    
    # Drop the class column, keep only CpGs
    beta_mat <- md[, !c("Class"), with = FALSE]
    
    # Transpose so rows = CpGs, cols = samples
    beta_mat <- t(as.matrix(beta_mat))
    
    # Set column names from sample annotation
    colnames(beta_mat) <- sa$Sample_accession
    
    beta_mat
  })
  
  
  classification <- eventReactive(input$classify, {
    req(beta_data(), training_data(), sample_anno())
    
    withProgress(message = "Running classification", value = 0, {
      # Read annotation directly into memory 
      sample_anno <- sample_anno()
      
      # Read training data directly into memory 
      train_data <- training_data()
      
      
      # Load reference beta matrix (row = CpGs, col = samples)
      ref_beta <- ref_beta()
      
      
      # prepare our imputation recipe exactly as in model development
      rf_recipe    <- recipe(Class~ ., data = train_data) %>%
        step_impute_median(all_predictors())
      
      # Learn the medians
      prepped_rec  <- prep(rf_recipe, training = train_data, retain = TRUE)
      
      prob_cols    <- c("Astrocyte","Ectoderm","Endoderm",
                        "Endothelial","iPSC","Lung",
                        "Mesoderm","NSC")
      # Prediciton
      # 1) pull out test sample
      incProgress(0.1, detail = "Extracting sample…")
      sample_name <- input$sample_accession
      test_sample <- beta_data()
      
      # 2) align columns & insert NAs
      incProgress(0.15, detail = "Aligning CpG sites…")
      expected_cpgs <- setdiff(colnames(train_data), "Class")
      missing_cpgs  <- setdiff(expected_cpgs, colnames(test_sample))
      if (length(missing_cpgs) > 0) {
        test_sample[, missing_cpgs] <- NA
      }
      test_sample <- test_sample[, expected_cpgs, drop = FALSE]
      
      # 3) bake recipe (impute NAs, factor setup)
      incProgress(0.2, detail = "Preparing data…")
      test_sample <- test_sample %>%
        as_tibble(rownames = "Sample") %>%
        mutate(Class = factor(NA, levels = levels(train_data$Class))) %>%
        column_to_rownames("Sample")
      baked <- bake(prepped_rec, new_data = test_sample)
      
      # 4) predict
      incProgress(0.4, detail = "Predicting…")
      
      beta_mat_cpg_by_sample <- t(as.matrix(test_sample))
      
      pred_df <- predict_calibrated_all(
        beta_mat_cpg_by_sample = beta_mat_cpg_by_sample,
        train_data   = train_data,
        prepped_rec  = prepped_rec,
        rf_fit_path  = url(urls$model),
        cal_path     = url(urls$cal_model),
        threshold    = 0.5,
        apply_reject = TRUE
      )
      
      
      
      # 5) collect results
      incProgress(0.9, detail = "Collecting results…")
      
      # pred_df should be 1 row for the selected sample
      stopifnot(nrow(pred_df) >= 1)
      
      # Use calibrated predictions
      class_pred <- as.character(pred_df$cal_pred[1])
      
      # Extract calibrated probabilities only
      # Instead of constructing names manually, detect them directly
      cal_prob_cols <- grep("^cal_", colnames(pred_df), value = TRUE)
      cal_prob_cols <- cal_prob_cols[!cal_prob_cols %in% c("cal_max_prob", "cal_pred")]
      
      message("pred_df columns: ", paste(colnames(pred_df), collapse = ", "))
      
      # named numeric vector (or list) of class probabilities
      probs <- as.list(stats::setNames(
        as.numeric(pred_df[1, cal_prob_cols]),
        prob_cols
      ))
      
      results <- tibble(
        Sample     = sample_name,
        Prediction = class_pred,
        MaxProb    = as.numeric(pred_df$cal_max_prob[1])
      ) %>%
        bind_cols(as_tibble(probs))
      
      # 6) finish up
      incProgress(1, detail = "Done.")
      
      list(
        prediction    = class_pred,
        max_prob      = as.numeric(pred_df$cal_max_prob[1]),
        probabilities = probs,
        results       = results,
        # Exact baked matrix used by the classifier; retained for RF attribution.
        test_sample   = baked,
        # Retain the exact preprocessing object and feature set for explanation.
        prepped_rec   = prepped_rec,
        pred_df       = pred_df,
        expected_cpgs = expected_cpgs
      )
      
    })
  })
  
  output$results_ready <- reactive({ !is.null(classification()) })
  outputOptions(output, "results_ready", suspendWhenHidden = FALSE)
  
  output$classification_ui <- renderUI({
    # wait until classification() has something
    res <- classification()
    req(res)
    
    tagList(
      h3("Classification Results"),
      div(
        style = "background-color: #f8f9fa;
               padding: 20px;
               border-radius: 8px;
               margin-bottom: 20px;",
        verbatimTextOutput("class_result")
      ),
      h3("Probability Chart"),
      plotOutput("probability_chart"),
      downloadButton("download_report","Download Classification Report")
    )
  })
  
  output$class_result <- renderText({
    res <- classification()
    req(res)
    
    sample_name     <- res$results$Sample
    predicted_class <- res$prediction
    
    display_class <- if (predicted_class == REJECT_LABEL) "Not Classifiable" else predicted_class
    
    display_score <- if (predicted_class == REJECT_LABEL) {
      paste0("No score above cut-off (max=", round(res$max_prob, 3), ")")
    } else {
      round(as.numeric(res$probabilities[[predicted_class]]), 3)
    }
    
    paste0(
      "Sample ID: ",         sample_name,   "\n",
      "Prediction: ",        display_class, "\n",
      "Prediction Score: ",  display_score
    )
  })
  
  
  
  output$probability_chart <- renderPlot({
    res <- classification()
    req(res)
    
    prob_df <- tibble(
      Class = names(res$probabilities),
      Probability = as.numeric(unlist(res$probabilities))
    ) %>%
      arrange(desc(Probability))
    
    max_class <- prob_df$Class[1]
    
    ggplot(prob_df, aes(x = reorder(Class, Probability), y = Probability, fill = Class == max_class)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "#6c757d", linewidth = 0.8) +
      coord_flip() +
      scale_fill_manual(values = c("FALSE" = "#6c757d", "TRUE" = "#0072B2"), guide = "none") +
      labs(title = "Calibrated class probabilities", x = "Class", y = "Probability") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title  = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 12)
      )
  })
  
  
  
  # =============================================================================
  # CNV PROFILING SERVER LOGIC
  # =============================================================================
  
  output$cnv_input_summary <- renderUI({
    tags$div(
      style = "font-size: 0.95em; color:#555;",
      tags$b("Current sample:"), " ", if (nzchar(input$sample_accession)) input$sample_accession else "not set",
      tags$br(),
      tags$b("Array:"), " ", if (!is.null(input$array_version)) input$array_version else "not set",
      tags$br(),
      tags$b("CNV reference:"), " ",
      if (!is.null(input$array_version)) cnv_reference_label(input$array_version) else "not set"
    )
  })
  
  cnv_result <- eventReactive(input$run_cnv, {
    req(input$sample_accession, input$array_version)
    pre <- preprocessed_data()
    req(pre, pre$mset_native)
    sf <- sample_files()
    req(sf)
    
    validate(
      need(input$cnv_bin_size >= 10000, "CNV bin size must be at least 10,000 bp.")
    )
    
    withProgress(
      message = paste0("Running CNV profiling: ", input$sample_accession,
                       " (", input$array_version, ")"),
      value = 0,
      {
        # Each CNV step reports its absolute progress and a description here,
        # and the same text is echoed to the R console.
        cnv_progress <- function(value, detail) {
          message("[CNV] ", detail)
          setProgress(value = value, detail = detail)
        }
        
        result <- tryCatch(
          run_cnv_profile(
            prefix = sf$prefix,
            mset = pre$mset_native,
            sample_id = input$sample_accession,
            platform = input$array_version,
            bin_size = as.integer(input$cnv_bin_size),
            progress = cnv_progress
          ),
          error = function(e) {
            showNotification(
              paste("CNV profiling failed:", conditionMessage(e)),
              type = "error",
              duration = NULL
            )
            stop(e)
          }
        )
        
        setProgress(value = 1, detail = "Rendering CNV plot and tables…")
        result
      })
  })   # no ignoreInit: actionButton's initial 0 is already ignored (ignoreNULL)
  
  output$cnv_results_ready <- reactive({
    !is.null(tryCatch(cnv_result(), error = function(e) NULL))
  })
  outputOptions(output, "cnv_results_ready", suspendWhenHidden = FALSE)
  
  output$cnv_result_header <- renderUI({
    res <- cnv_result()
    req(res)
    div(
      style = "background:#f8f9fa; border-left:4px solid #0072B2; padding:12px 16px; margin-bottom:16px;",
      tags$b("Sample: "), res$sample_id,
      tags$span(" | "),
      tags$b("Array: "), res$platform,
      tags$span(" | "),
      tags$b("Genome: "), res$genome,
      tags$span(" | "),
      tags$b("Reference: "), res$reference_label,
      if (!is.null(res$reference_note)) tags$div(
        style = "margin-top:6px; color:#D55E00;", icon("triangle-exclamation"), " ", res$reference_note
      ),
      tags$span(" | "),
      tags$b("Bin size: "), format(res$bin_size, big.mark = ","), " bp",
      tags$div(
        style = "margin-top:6px;",
        tags$b("Predicted sex: "), res$inferred_sex,
        tags$span(style = "color:#6c757d; font-size:0.9em;", paste0(" - ", res$sex_note))
      )
    )
  })
  
  output$cnv_genome_plot <- renderPlot({
    res <- cnv_result()
    req(res)
    cnv_draw_genomeplot(res)
  })
  
  output$cnv_focal_table <- DT::renderDT({
    res <- cnv_result()
    req(res)
    dat <- as.data.frame(res$alterations, check.names = FALSE)
    
    if (nrow(dat) == 0) {
      dat <- data.frame(
        Result = "No copy-number segments beyond the gain/loss cut-offs and no focal calls detected.",
        stringsAsFactors = FALSE
      )
      return(DT::datatable(dat, rownames = FALSE, options = list(dom = "t")))
    }
    
    # Citations are kept for the download; the on-screen table stays compact.
    show <- dat[, setdiff(names(dat), "Hotspot citations"), drop = FALSE]
    # Natural sort so clicking the Chr / Alteration headers orders chr2 before
    # chr10 and 2p before 10q (DataTables would otherwise sort as text).
    natural_cols <- which(names(show) %in% c("Chr", "Alteration")) - 1L
    DT::datatable(
      show,
      rownames = FALSE,
      filter = "top",
      plugins = "natural",
      options = list(
        pageLength = 15, scrollX = TRUE,
        columnDefs = list(list(type = "natural", targets = natural_cols))
      )
    ) %>%
      DT::formatStyle(
        "Stem cell hotspot",
        fontWeight = DT::styleEqual("No", "normal", default = "bold")
      ) %>%
      DT::formatStyle(
        "CNV state",
        color = DT::styleEqual(c("AMP", "GAIN", "HOMDEL", "LOSS"),
                               c("#E17E49", "#E17E49", "#39379F", "#39379F"))
      )
  })
  
  output$download_cnv_plot <- downloadHandler(
    filename = function() {
      paste0(gsub("[^A-Za-z0-9._-]", "_", input$sample_accession), "_CNV_profile.pdf")
    },
    content = function(file) {
      res <- cnv_result()
      req(res)
      grDevices::pdf(file, width = 18, height = 6)
      on.exit(grDevices::dev.off(), add = TRUE)
      cnv_draw_genomeplot(res)
    }
  )
  
  # Same fit as the labelled plot, without gene labels (as in the standalone script).
  output$download_cnv_plot_unlabelled <- downloadHandler(
    filename = function() {
      paste0(gsub("[^A-Za-z0-9._-]", "_", input$sample_accession), "_CNV_profile_unlabelled.pdf")
    },
    content = function(file) {
      res <- cnv_result()
      req(res)
      grDevices::pdf(file, width = 18, height = 6)
      on.exit(grDevices::dev.off(), add = TRUE)
      cnv_draw_genomeplot(res, detail = FALSE)
    }
  )
  
  output$download_cnv_segments <- downloadHandler(
    filename = function() {
      paste0(gsub("[^A-Za-z0-9._-]", "_", input$sample_accession), "_CNV_segments.tsv")
    },
    content = function(file) {
      res <- cnv_result()
      req(res)
      utils::write.table(
        res$segments,
        file = file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  )
  
  output$download_cnv_focal <- downloadHandler(
    filename = function() {
      paste0(gsub("[^A-Za-z0-9._-]", "_", input$sample_accession), "_CNV_alterations.tsv")
    },
    content = function(file) {
      res <- cnv_result()
      req(res)
      utils::write.table(
        as.data.frame(res$alterations, check.names = FALSE),
        file = file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  )
  
  observeEvent(input$reset, {
    session$reload()
  })
  
  
  output$umap_plot <- renderPlot({
    # re-run whenever the button is clicked
    # will only proceed once you've clicked
    req(input$generate_umap > 0)
    input$generate_umap
    isolate({
      withProgress(message = "Generating UMAP", value = 0, {
        
        # pull in the baked test sample
        incProgress(0.1, detail = "Extracting test sample…")
        req(classification())
        test_sample <- classification()$test_sample
        
        # grab training data
        incProgress(0.1, detail = "Loading training data…")
        train_df <- training_data()
        
        # feature alignment
        incProgress(0.1, detail = "Aligning features…")
        feature_cols <- colnames(test_sample)
        missing_train <- setdiff(feature_cols, colnames(train_df))
        if (length(missing_train) > 0) {
          stop("Missing predictors in training data: ",
               paste(missing_train, collapse = ", "))
        }
        train_features <- train_df[, ..feature_cols]
        
        # combine datasets
        incProgress(0.1, detail = "Combining data…")
        combined_data <- rbind(train_features, test_sample)
        
        # prepare labels
        incProgress(0.1, detail = "Preparing labels…")
        if (!"Class" %in% colnames(train_df)) {
          stop("training_data() must contain a 'Class' column")
        }
        train_labels <- train_df$Class
        labels <- factor(
          c(as.character(train_labels), rep("Test Sample", nrow(test_sample))),
          levels = c(unique(as.character(train_labels)), "Test Sample")
        )
        
        # run UMAP
        incProgress(0.4, detail = "Running UMAP…")
        set.seed(123)
        umap_res <- uwot::umap(
          combined_data,
          n_neighbors = 15,
          min_dist     = 0.1,
          metric       = "euclidean"
        )
        
        # build plotting frame
        incProgress(0.1, detail = "Creating plot…")
        umap_df <- data.frame(
          UMAP1       = umap_res[,1],
          UMAP2       = umap_res[,2],
          Label       = gsub("_", " ", c(as.character(train_labels), "Test Sample")),
          Sample_Type = rep(c("Training Sample", "Test Sample"),
                            c(length(train_labels), nrow(test_sample)))
        )
        
        # final plot
        ggplot(umap_df, aes(x = UMAP1, y = UMAP2,
                            color = Label, shape = Sample_Type)) +
          geom_point(size = 3, alpha = 0.8) +
          labs(
            title = "Unsupervised UMAP of Test Sample + Training Data",
            color = "Cell Type"
          ) +
          theme_minimal(base_size = 14) +
          theme(
            plot.title   = element_text(hjust = 0.5, face = "bold"),
            legend.title = element_text(face = "bold"),
            legend.text  = element_text(size = 12)
          ) +
          scale_color_manual(values = c(
            "Astrocyte"    = "#000080",
            "Endothelial"  = "#D55E00",
            "Mesoderm"     = "#E69F00",
            "iPSC"         = "#CC79A7",
            "Ectoderm"     = "#9ac7df",
            "NSC"          = "#0072B2",
            "Endoderm"     = "#a6dba0",
            "Lung"         = "#238b45",
            "Test Sample"  = "#FF00FF"
          )) +
          scale_shape_manual(values = c(
            "Training Sample" = 16,
            "Test Sample"     = 18
          ))
      })
    })
  })
  
  
  marker_heatmap_data <- eventReactive(input$marker_plot_button, {
    req(beta_data(), input$sample_accession, marker_beta(), ann450K(), sample_anno())
    # Load reference beta matrix (row = CpGs, col = samples)
    marker_beta <- marker_beta()
    ann450K <- ann450K()
    marker      <- input$marker
    sample_name <- input$sample_accession
    
    
    probes_for_marker <- ann450K %>% filter(Marker == marker)
    if (nrow(probes_for_marker) == 0) {
      showNotification("No CpG probes found for that marker.", type = "error")
      return(NULL)
    }
    
    cpg_ids <- probes_for_marker$Name
    
    
    # Filter to CpGs of interest in marker_beta
    common_cpgs <- intersect(cpg_ids, intersect(rownames(marker_beta), colnames(beta_data())))
    if (length(common_cpgs) == 0) {
      showNotification("No shared CpGs found between reference and sample for this gene.", type = "error")
      return(NULL)
    }
    beta_values_ref <- marker_beta[common_cpgs, , drop = FALSE]
    test_vec        <- beta_data()[, common_cpgs, drop = FALSE]
    test_mat        <- t(test_vec)
    rownames(test_mat) <- common_cpgs
    colnames(test_mat) <- sample_name
    
    
    # bind test sample & build annotation
    sample_anno <- sample_anno()
    
    ref_anno_full    <- sample_anno[match(colnames(beta_values_ref),
                                          sample_anno$Sample_accession), ]
    keep_ref         <- which(ref_anno_full$Class %in% c("iPSC", marker))
    beta_values_ref  <- beta_values_ref[, keep_ref, drop = FALSE]
    ref_anno         <- ref_anno_full[keep_ref, ]
    
    
    beta_values <- cbind(beta_values_ref, test_mat)
    
    ref_anno$Class <- as.character(ref_anno$Class)
    test_row <- data.frame(Sample_accession = sample_name,
                           Class            = "Test Sample",
                           stringsAsFactors = FALSE)
    anno_plot <- bind_rows(ref_anno, test_row)
    rownames(anno_plot) <- anno_plot$Sample_accession
    
    list(
      beta   = beta_values,
      anno   = anno_plot,
      marker = marker
    )
  })
  
  output$marker_heatmap_plot <- renderPlot({
    dat <- marker_heatmap_data()
    req(dat)
    
    hm <- Heatmap(
      dat$beta,
      name = "Beta",
      top_annotation = HeatmapAnnotation(
        Cell_Type = dat$anno$Class,
        col = list(Cell_Type = c(
          Astrocyte   = "#000080", Endothelial="#D55E00",
          Mesoderm    = "#E69F00", iPSC      ="#CC79A7",
          Ectoderm    ="#9ac7df",  NSC       ="#0072B2",
          Endoderm    ="#a6dba0",  Lung      ="#238b45",
          `Test Sample`="#000000"
        )),
        annotation_legend_param=list(Cell_Type=list(title="Cell Type"))
      ),
      show_row_names    = FALSE,
      show_column_names = FALSE,
      cluster_rows      = TRUE,
      cluster_columns   = TRUE,
      row_title         = paste("CpGs for", dat$marker),
      col                = heatmap_scale_blues
    )
    draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
  }, height=500)
  
  
  # Build a Heatmap object when the user clicks “Plot Heatmap”
  gene_heatmap_obj <- eventReactive(input$plot_button, {
    req(beta_data(), input$gene_input, input$celltype, ann450K(), ref_beta(), sample_anno())
    
    withProgress(message = "Generating gene‐level heatmap", value = 0, {
      marker_beta <- ref_beta()
      ann450K <- ann450K()
      marker      <- input$marker
      sample_name <- input$sample_accession
      
      
      # look up probes for this gene
      incProgress(0.1, detail = "Finding CpGs for gene…")
      celltype    <- input$celltype
      sample_name <- input$sample_accession
      gene_name   <- toupper(input$gene_input)
      
      #prepare CpG annotation
      ann450K <- ann450K()
      probes_for_gene <- ann450K %>%
        filter(UCSC_RefGene_Name == gene_name) %>%
        distinct(Name, .keep_all = TRUE) %>%
        dplyr::select(-Marker)
      if (nrow(probes_for_gene) == 0) {
        showNotification("No CpG probes found for the entered gene.", type = "error")
        return(NULL)
      }
      cpg_ids <- probes_for_gene$Name
      
      # filter bvals
      incProgress(0.3, detail = "Filtering beta values")
      # Filter to CpGs of interest in marker_beta
      
      # Filter to CpGs of interest in marker_beta
      common_cpgs <- intersect(cpg_ids, intersect(rownames(marker_beta), colnames(beta_data())))
      if (length(common_cpgs) == 0) {
        showNotification("No shared CpGs found between reference and sample for this gene.", type = "error")
        return(NULL)
      }
      beta_values_ref <- marker_beta[common_cpgs, , drop = FALSE]
      test_vec        <- beta_data()[, common_cpgs, drop = FALSE]
      test_mat        <- t(test_vec)
      rownames(test_mat) <- common_cpgs
      colnames(test_mat) <- sample_name
      
      beta_values <- cbind(beta_values_ref, test_mat)
      sample_anno <- sample_anno()
      samp_anno_heatmap <- sample_anno[match(colnames(beta_values), sample_anno$Sample_accession), ]
      samp_anno_heatmap <- samp_anno_heatmap %>% 
        filter(Class %in% c(celltype, "iPSC"))
      keep_cols <- samp_anno_heatmap$Sample_accession
      beta_values <- beta_values[, c(keep_cols, sample_name), drop = FALSE]
      
      test_row <- data.frame(
        Sample_accession = sample_name,
        Class            = "Test Sample",
        stringsAsFactors = FALSE
      )
      samp_anno_plot <- bind_rows(
        samp_anno_heatmap %>% mutate(Class = as.character(Class)),
        test_row
      ) %>% column_to_rownames("Sample_accession")
      
      col.ha <- HeatmapAnnotation(
        Cell_Type = samp_anno_plot$Class,
        col = list(Cell_Type = c(
          Astrocyte     = "#000080",
          Endothelial   = "#D55E00",
          Mesoderm      = "#E69F00",
          iPSC          = "#CC79A7",
          Ectoderm      = "#9ac7df",
          NSC           = "#0072B2",
          Endoderm      = "#a6dba0",
          Lung          = "#238b45",
          `Test Sample` = "#000000"
        )),
        annotation_legend_param = list(Cell_Type = list(title = "Cell Type"))
      )
      
      # row annotations & sorting
      incProgress(0.2, detail = "Sorting CpGs & building row annotations…")
      probes_in_heatmap <- rownames(beta_values)
      ann450K_heatmap   <- ann450K[match(probes_in_heatmap, ann450K$Name), ]
      ann450K_heatmap$chr <- factor(ann450K_heatmap$chr,
                                    levels = paste0("chr", c(1:22, "X", "Y", "M")))
      ann450K_heatmap_sorted <- ann450K_heatmap[order(
        ann450K_heatmap$chr, ann450K_heatmap$pos
      ), ]
      beta_values <- beta_values[ann450K_heatmap_sorted$Name, , drop = FALSE]
      row.ha_data <- ann450K_heatmap_sorted$UCSC_RefGene_Group
      names(row.ha_data) <- ann450K_heatmap_sorted$Name
      row.ha <- rowAnnotation(
        Gene_region = row.ha_data,
        col = list(Gene_region = c(
          "5'UTR" = "#D55E00",
          "Body"  = "#000080",
          "TSS1500" = "#CC79A7",
          "1stExon" = "#0072B2",
          "3'UTR"   = "#b2b2d8"
        )),
        annotation_legend_param = list(
          Gene_region = list(
            title = "Gene Region",
            at    = c("5'UTR","Body","TSS1500","1stExon","3'UTR"),
            labels= c("5' UTR","Body","TSS1500","1st Exon","3' UTR")
          )
        )
      )
      
      # dynamically compute height
      num_probes  <- nrow(beta_values)
      plot_height <- min(400 + num_probes * 20, 2000)
      
      # render the heatmap
      incProgress(0.2, detail = "Rendering heatmap…")
      Heatmap(
        beta_values,
        name               = "Beta Value",
        show_row_names     = TRUE,
        show_column_names  = FALSE,
        cluster_rows       = FALSE,
        cluster_columns    = TRUE,
        top_annotation     = col.ha,
        right_annotation   = row.ha,
        col                = heatmap_scale_blues,
        row_title          = paste("CpG Sites for", gene_name),
        column_title       = "Samples",
        heatmap_legend_param = list(title = "Beta Value")
      )
    })
    
  })
  
  # Render the heatmap onscreen
  output$heatmap_plot <- renderPlot({
    ht <- gene_heatmap_obj()
    req(ht)
    mat <- ht@matrix
    draw(ht,
         heatmap_legend_side    = "right",
         annotation_legend_side = "right")
  }, height = function() {
    ht <- gene_heatmap_obj()
    mat <- ht@matrix
    min(400 + nrow(mat)*20, 2000)
  })
  
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("classification_report_", input$sample_accession, ".html")
    },
    content = function(file) {
      # copy the report template
      tempRmd <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempRmd, overwrite = TRUE)
      # copy the CSS alongside it
      tempCss <- file.path(tempdir(), "report.css")
      file.copy("report.css", tempCss, overwrite = TRUE)
      # recreate the summary text
      res <- classification()
      req(res)
      
      sample_name     <- res$results$Sample
      predicted_class <- res$prediction
      
      class_text <- paste0(
        "Sample ID: ",       sample_name,     "\n",
        "Prediction: ",      if (predicted_class=="reject") "Not Classifiable" else predicted_class, "\n",
        "Prediction Score: ",if (predicted_class=="reject") "No score above cut-off"
        else round(as.numeric(res$probabilities[[predicted_class]]),3)
      )
      
      # rebuild the two plots exactly as before
      prob_plot <- { 
        res <- classification()
        req(res)
        
        # turn the named list into a data frame
        prob_df <- as.data.frame(res$probabilities, stringsAsFactors = FALSE) %>%
          pivot_longer(
            cols      = everything(),
            names_to  = "Class",
            values_to = "Probability"
          ) %>%
          arrange(desc(Probability))
        
        max_class <- prob_df$Class[1]  # highest-probability class
        
        ggplot(prob_df, aes(
          x    = reorder(Class, Probability),
          y    = Probability,
          fill = Class == max_class
        )) +
          geom_col(width = 0.6) +
          # cutoff line at 0.5
          geom_hline(
            yintercept = 0.5,
            linetype   = "dashed",
            color      = "#6c757d",
            size       = 0.8
          ) +
          coord_flip() +
          scale_fill_manual(
            values = c("FALSE" = "#6c757d", "TRUE" = "#0072B2"),
            guide  = "none"
          ) +
          labs(
            title = "Class Probabilities",
            x     = "Class",
            y     = "Probability"
          ) +
          ylim(0, 1) +
          scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
          theme_minimal(base_size = 14) +
          theme(
            plot.title     = element_text(hjust = 0.5, face = "bold"),
            axis.text.y    = element_text(face = "bold", size = 12),
            axis.text.x    = element_text(size = 12)
          )
      }
      
      umap_plot <- {
        # copy your renderPlot UMAP code, e.g.:
        test_sample <- res$test_sample
        train_df    <- as.data.frame(training_data())
        feature_cols <- colnames(test_sample)
        train_features <- train_df[, feature_cols, drop=FALSE]
        combined_data <- rbind(train_features, test_sample)
        set.seed(123)
        umap_res <- uwot::umap(combined_data, n_neighbors=15, min_dist=0.1, metric="euclidean")
        train_labels <- training_data()$Class
        umap_df <- data.frame(
          UMAP1 = umap_res[,1], UMAP2 = umap_res[,2],
          Label = c(as.character(train_labels), "Test Sample"),
          Sample_Type = rep(c("Training","Test"), c(length(train_labels),1))
        )
        ggplot(umap_df, aes(UMAP1, UMAP2, color=Label, shape=Sample_Type))+
          geom_point(size=3, alpha=0.8)+
          theme_minimal(base_size=14)+
          labs(title="UMAP Projection", color="Cell Type")+
          scale_color_manual(values=c(
            "Astrocyte"="#000080","Endothelial"="#D55E00",
            "Mesoderm"="#E69F00","iPSC"="#CC79A7",
            "Ectoderm"="#9ac7df","NSC"="#0072B2",
            "Endoderm"="#a6dba0","Lung"="#238b45","Test Sample"="#000000"
          ))+
          scale_shape_manual(values=c("Training"=16,"Test"=17))
      }
      
      # render the HTML
      rmarkdown::render(
        tempRmd,
        output_file = file,
        params      = list(
          class_text = class_text,
          prob_plot  = prob_plot,
          umap_plot  = umap_plot
        ),
        envir = new.env(parent = globalenv()),
        quiet = TRUE
      )
    }
  )
  
  observeEvent(input$download_marker_heatmap, {
    shinyjs::show("download_spinner")
  })
  
  output$download_marker_heatmap <- downloadHandler(
    filename = function() {
      paste0("marker_heatmap_", input$marker, "_", input$sample_accession, ".png")
    },
    content = function(file) {
      on.exit(shinyjs::hide("download_spinner"), add=TRUE)
      # open a PNG device
      png(filename = file, width = 1200, height = 800, res = 150)
      
      # reproduce exactly what you did in renderPlot
      dat <- marker_heatmap_data()
      req(dat)
      ht <- Heatmap(
        dat$beta,
        name = "Beta",
        top_annotation = HeatmapAnnotation(
          Cell_Type = dat$anno$Class,
          col = list(Cell_Type = c(
            Astrocyte   = "#000080", Endothelial="#D55E00",
            Mesoderm    = "#E69F00", iPSC      ="#CC79A7",
            Ectoderm    ="#9ac7df",  NSC       ="#0072B2",
            Endoderm    ="#a6dba0",  Lung      ="#238b45",
            `Test Sample`="#000000"
          )),
          annotation_legend_param=list(Cell_Type=list(title="Cell Type"))
        ),
        show_row_names    = FALSE,
        show_column_names = FALSE,
        cluster_rows      = TRUE,
        cluster_columns   = TRUE,
        row_title         = paste("CpGs for", dat$marker),
        col                = heatmap_scale_blues
      )
      draw(ht, heatmap_legend_side="right", annotation_legend_side="right")
      
      dev.off()
      
    }
  )
  
  session$onSessionEnded(function() {
    session_dir <- file.path(tempdir(), paste0("stemclass_", session$token))
    if (dir.exists(session_dir)) unlink(session_dir, recursive = TRUE, force = TRUE)
  })
  
  output$download_gene_heatmap <- downloadHandler(
    filename = function() {
      paste0("gene_heatmap_", toupper(input$gene_input), "_",
             input$sample_accession, ".png")
    },
    content = function(file) {
      ht  <- gene_heatmap_obj()
      mat <- ht@matrix
      h   <- min(400 + nrow(mat)*20, 2000)  # this “h” is in px
      
      # Ask for a 1200×h pixel image at 150 DPI (or any DPI you like)
      png(
        file,
        width  = 1200,
        height = h,
        res    = 150
      )
      draw(ht,
           heatmap_legend_side    = "right",
           annotation_legend_side = "right")
      dev.off()
    }
  )
  
  
  
}

shinyApp(ui = ui, server = server)
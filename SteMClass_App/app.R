
options(shiny.maxRequestSize = 30*1024^2)

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
  cal_model = "https://github.com/pereze5/SteMClass_App/releases/download/v1.0-data/calibration_model_deploy.rds"
)

# 9‐color Blues palette; interpolate to 255

blues_base <- brewer.pal(9, "Blues")
blues_pal  <- colorRampPalette(blues_base)(255)

heatmap_scale_blues <- colorRamp2(
  seq(0, 1, length.out = 255),
  blues_pal
)

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
"


# Define the Shiny UI with the selected theme and custom CSS
ui <- navbarPage(
  title = "SteMClass",
  theme = theme,
  useShinyjs(),
  tags$head(
    tags$style(HTML(custom_css))
  ),
  tabPanel(
    title = "Classification",
    sidebarLayout(
      sidebarPanel(
        fileInput(
          inputId = "idat_files",
          label   = "Upload IDAT Files",
          multiple = TRUE,
          accept   = ".idat"
        ),
        
        tags$div(
          id = "form",
          
          textInput(
            inputId = "sample_accession",
            label   = "Enter Sample ID:"
          ),
          
          selectInput(
            inputId = "array_version",
            label   = "Array Version:",
            choices = c("450K", "EPICv1", "EPICv2")
          ),
          
          actionButton(
            inputId = "load_samples",
            label   = "Preprocess Data",
            icon    = icon("broom"),
            style   = "margin-bottom: 10px;"
          ),
          
          actionButton(
            inputId = "classify",
            label   = "Run Classification",
            icon    = icon("play"),
            style   = "margin-bottom: 10px;"
          ),
          
          actionButton(
            inputId = "reset",
            label   = "New Analysis",
            icon    = icon("redo")
          ),
          
          hr(),
          
          h4("Instructions:"),
          tags$ol(
            tags$li("Upload a pair of IDAT files (one Red and one Grn)."),
            tags$li("Specify the Sample ID and choose the Array version (450K, EPICv1 or EPICv2)."),
            tags$li("Click “Preprocess Data”"),
            tags$li("Click “Run Classification” to classify the sample."),
            tags$li("Click “New Analysis” to classify a new sample.")
          )
        ) 
      ),  
      
      mainPanel(
        div(style = "position: relative; padding-bottom: 60px;",
            uiOutput("classification_ui"))
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
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  beta_data <- reactive({
    req(input$idat_files, input$sample_accession, input$array_version)
    if (nrow(input$idat_files) != 2) {
      stop("Please upload exactly two IDAT files (one Red and one Grn).")
    }
    
    withProgress(message = "Reading & preprocessing IDATs", value = 0, {
      
      # 1. Copy files
      incProgress(0.1, detail = "Copying IDAT files…")
      tmp <- tempdir()
      file.copy(input$idat_files$datapath,
                file.path(tmp, input$idat_files$name),
                overwrite = TRUE)
      
      # 2. Read raw data
      incProgress(0.2, detail = "Loading raw data…")
      basename0 <- sub("_(Grn|Red)\\.idat$", "", input$idat_files$name[1])
      targets <- data.frame(
        Sample_accession = input$sample_accession,
        Array_version    = input$array_version,
        Basename         = basename0,
        stringsAsFactors = FALSE
      )
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
      
      rgSet <- read.metharray.exp(base = tmp, targets = targets)
      sampleNames(rgSet) <- targets$Sample_accession
      
      # Preprocess (split out each branch so we can show progress)
      
      if (input$array_version == "EPICv2") {
        incProgress(0.3, detail = "Preprocessing (Noob)…")
        rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
        mSetv2 <- preprocessNoob(rgSet)
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
        
        incProgress(0.05, detail = "Converting to 450K…")
        mSet <- convertArray(mSet,
                             outType = "IlluminaHumanMethylation450k")
        
        incProgress(0.2, detail = "Mapping to genome…")
        mSetSq <- mapToGenome(mSet)
        
        beta   <- getBeta(mSetSq)
        
      } else if (input$array_version == "450K") {
        incProgress(0.3, detail = "Preprocessing (Noob)…")
        mSet <- preprocessNoob(rgSet)
        incProgress(0.35, detail = "Mapping to genome…")
        mSetSq <- mapToGenome(mSet)
        beta   <- getBeta(mSetSq)
        
      } else {
        stop("Unknown array version: ", input$array_version)
      }
      
      # Finalise
      incProgress(0.1, detail = "Finalizing…")
      beta_mat <- t(beta)
      
      # return
      beta_mat
    })
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
  
  # Reactive expression to wrap the reference beta matrix
  ref_beta <- reactive({
    td <- training_data()
    sa <- sample_anno()
    
    # Drop the class column, keep only CpGs
    beta_mat <- td[, !c("Class"), with = FALSE]
    
    # Transpose so rows = CpGs, cols = samples
    beta_mat <- t(as.matrix(beta_mat))
    
    # Set column names from sample annotation
    colnames(beta_mat) <- sa$Sample_accession
    
    beta_mat
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
      
      # Load model directly from URL
      rf_fit <- readRDS(url(urls$model))
      
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
      
      predict_calibrated <- function(
    new_data,
    rf_fit_path  = url(urls$model),
    cal_path     = url(urls$cal_model),
    class_names  = c("Astrocyte","Ectoderm","Endoderm","Endothelial","iPSC","Lung","Mesoderm","NSC"),
    threshold    = 0.5,         
    apply_reject = TRUE
      ) {
        # 1) Load model bundle + calibrator
        rf_bundle <- readRDS(rf_fit_path)
        cal_deploy <- readRDS(cal_path)
        
        # in case need to use bundled clas and prob rds:
        rf_fit <- if (is.list(rf_bundle) && "model_fit" %in% names(rf_bundle)) rf_bundle$model_fit else rf_bundle
        
        # 2) RF raw probabilities
        raw_probs <- predict(rf_fit, new_data = new_data, type = "prob") %>% as_tibble()
        colnames(raw_probs) <- gsub("^\\.pred_", "", colnames(raw_probs))
        
        # Ensure class order and presence
        missing <- setdiff(class_names, colnames(raw_probs))
        if (length(missing) > 0) stop("RF probability output missing classes: ", paste(missing, collapse = ", "))
        raw_probs <- raw_probs[, class_names, drop = FALSE]
        
        # 3) Apply calibration model to raw probs
        cal_obj <- cal_deploy
        if (is.list(cal_obj) && "cv" %in% names(cal_obj)) cal_obj <- cal_obj$cv
        cal_probs <- apply_calibrator(cal_obj, raw_probs, class_names)
        
        
        # 4) rejection + predicted class using calibrated probs
        if (apply_reject) {
          max_prob <- do.call(pmax, cal_probs)
          argmax <- class_names[max.col(as.matrix(cal_probs), ties.method = "first")]
          pred_class <- ifelse(max_prob < threshold, REJECT_LABEL, argmax)
          
          out <- bind_cols(
            raw_probs %>% rename_with(~ paste0("raw_", .x)),
            cal_probs %>% rename_with(~ paste0("cal_", .x)),
            tibble(
              cal_max_prob = max_prob,
              cal_pred     = factor(pred_class, levels = c(class_names, REJECT_LABEL))
            )
          )
        } else {
          out <- bind_cols(
            raw_probs %>% rename_with(~ paste0("raw_", .x)),
            cal_probs %>% rename_with(~ paste0("cal_", .x))
          )
        }
        
        out
      }
      pred_df <- predict_calibrated(
        new_data    = baked,      # must match RF workflow feature columns / recipe
        threshold   = 0.5,
        apply_reject = TRUE
      )
      
      
      # 5) collect results
      incProgress(0.1, detail = "Collecting results…")
      prob_cols_df <- setdiff(names(pred_df), c("max_prob", "pred_class"))
      class_pred   <- as.character(pred_df$pred_class[1])
      probs        <- as.list(pred_df[1, prob_cols_df, drop = FALSE])
      
      results <- tibble(
        Sample     = sample_name,
        Prediction = class_pred
      ) %>%
        bind_cols(as_tibble(probs))
      
      # 6) finish up
      incProgress(0.05, detail = "Done.")
      list(
        prediction    = class_pred,
        probabilities = probs,
        results       = results,
        test_sample   = baked
      )
    })
  })
  
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
    
    # Decide what to show for class & score
    display_class <- if (predicted_class == "reject") {
      "Not Classifiable"
    } else {
      predicted_class
    }
    
    display_score <- if (predicted_class == "reject") {
      "No score above cut-off"
    } else {
      # safe convert & round
      round(as.numeric(res$probabilities[[predicted_class]]), 3)
    }
    
    paste0(
      "Sample ID: ",           sample_name,   "\n",
      "Prediction: ",       display_class, "\n",
      "Prediction Score: ", display_score
    )
  })
  
  
  output$probability_chart <- renderPlot({
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
      # cutoff line at 0.6
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
      scale_y_continuous(
        limits = c(0, 1),                
        breaks = seq(0, 1, by = 0.1),    
        expand = c(0, 0)                 
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title     = element_text(hjust = 0.5, face = "bold"),
        axis.text.y    = element_text(face = "bold", size = 12),
        axis.text.x    = element_text(size = 12)
      )
  })
  
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
    req(beta_data(), input$gene_input, input$celltype, ann450K(), marker_data(), sample_anno())
    
    withProgress(message = "Generating gene‐level heatmap", value = 0, {
      marker_data <- marker_data()
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
      # Filter to CpGs of interest in marker_data
      
      common_cpgs <- intersect(cpg_ids, intersect(rownames(marker_data), colnames(beta_data())))
      if (length(common_cpgs) == 0) {
        showNotification("No shared CpGs found between reference and sample for this gene.", type = "error")
        return(NULL)
      }
      beta_values_ref <- marker_data[common_cpgs, , drop = FALSE]
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
          # cutoff line at 0.6
          geom_hline(
            yintercept = 0.6,
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










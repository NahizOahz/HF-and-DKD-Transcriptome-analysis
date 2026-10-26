set.seed(20260820)

options(stringsAsFactors = FALSE, 
    width = 120)

suppressPackageStartupMessages({
    library(DESeq2)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(sandwich)
    library(lmtest)
})

project <- "C:/ZIHAN/Piezo1/DKD/version3"

input_dir <- "C:/ZIHAN/Piezo1/DKD"

count_file <- file.path(input_dir, 
    "GSE142025_raw_counts_GRCh38.p13_NCBI.tsv")

meta_file <- file.path(input_dir, 
    "GSE142025_series_matrix.txt.gz")

meta_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142025/matrix/", 
    "GSE142025_series_matrix.txt.gz")

dirs <- file.path(project, c("code", 
    "results/tables", "results/figures", 
    "results/objects", "logs"))

invisible(lapply(dirs, dir.create, 
    recursive = TRUE, showWarnings = FALSE))

tab_dir <- file.path(project, "results/tables")

fig_dir <- file.path(project, "results/figures")

individual_fig_dir <- file.path(fig_dir, 
    "individual")

obj_dir <- file.path(project, "results/objects")

dir.create(individual_fig_dir, 
    recursive = TRUE, showWarnings = FALSE)

reuse_fit <- identical(Sys.getenv("DKD_REUSE_FIT"), 
    "1")

stopifnot(file.exists(count_file))

if (!file.exists(meta_file)) download.file(meta_url, 
    meta_file, mode = "wb")

group_levels <- c("Control", "Early_DN", 
    "Advanced_DN")

group_labels <- c(Control = "Control", 
    Early_DN = "Early DN", Advanced_DN = "Advanced DN")

group_cols <- c(Control = "grey25", 
    Early_DN = "#aec0d8", Advanced_DN = "#ef8a7b")

write_csv <- function(x, name) write.csv(x, 
    file.path(tab_dir, name), row.names = FALSE, 
    na = "")

fmt_p <- function(p) ifelse(is.na(p), 
    "NA", ifelse(p < 0.001, format(p, 
        digits = 2, scientific = TRUE), 
        sprintf("%.3f", p)))

strip_quote <- function(x) sub("\"$", 
    "", sub("^\"", "", x))

geo_rows <- function(lines, key) {
    hit <- lines[startsWith(lines, 
        paste0(key, "\t"))]
    lapply(hit, function(z) strip_quote(strsplit(z, 
        "\t", fixed = TRUE)[[1]][-1]))
}

read_geo_meta <- function(path) {
    lines <- readLines(gzfile(path), 
        warn = FALSE)
    out <- data.frame(gsm = geo_rows(lines, 
        "!Sample_geo_accession")[[1]], 
        subject_id = geo_rows(lines, 
            "!Sample_title")[[1]], 
        source = geo_rows(lines, 
            "!Sample_source_name_ch1")[[1]])
    for (values in geo_rows(lines, 
        "!Sample_characteristics_ch1")) {
        key <- make.names(tolower(trimws(sub(":.*$", 
            "", values[1]))))
        out[[key]] <- trimws(sub("^[^:]+:\\s*", 
            "", values))
    }
    out
}

robust_table <- function(fit) {
    x <- coeftest(fit, vcov. = vcovHC(fit, 
        type = "HC3"))
    crit <- qt(0.975, df.residual(fit))
    data.frame(term = rownames(x), 
        beta = x[, 1], robust_SE = x[, 
            2], statistic = x[, 
            3], p_value = x[, 4], 
        CI_low = x[, 1] - crit * 
            x[, 2], CI_high = x[, 
            1] + crit * x[, 2])
}

average_prediction_hc3 <- function(fit, 
    newdata) {
    xmat <- model.matrix(delete.response(terms(fit)), 
        newdata)
    beta <- coef(fit)
    keep_coef <- !is.na(beta)
    beta <- beta[keep_coef]
    xmat <- xmat[, names(beta), 
        drop = FALSE]
    robust_vcov <- vcovHC(fit, 
        type = "HC3")
    robust_vcov <- robust_vcov[names(beta), 
        names(beta), drop = FALSE]
    xbar <- colMeans(xmat)
    estimate <- drop(xbar %*% beta)
    robust_se <- sqrt(drop(t(xbar) %*% 
        robust_vcov %*% xbar))
    crit <- qt(0.975, df.residual(fit))
    c(estimate = estimate, robust_SE = robust_se, 
        CI_low = estimate - crit * 
            robust_se, CI_high = estimate + 
            crit * robust_se)
}

counts_df <- read.delim(count_file, 
    check.names = FALSE, quote = "", 
    comment.char = "")

gene_id <- as.character(counts_df[[1]])

count_mat <- as.matrix(counts_df[-1])

storage.mode(count_mat) <- "integer"

rownames(count_mat) <- gene_id

rm(counts_df)

stopifnot(all(count_mat >= 0), 
    all(count_mat == round(count_mat)), 
    !anyDuplicated(colnames(count_mat)))

if (anyDuplicated(rownames(count_mat))) {
    count_mat <- rowsum(count_mat, 
        rownames(count_mat), reorder = FALSE)
    storage.mode(count_mat) <- "integer"
}

meta_all <- read_geo_meta(meta_file)

meta <- meta_all[match(colnames(count_mat), 
    meta_all$gsm), ]

stopifnot(!anyNA(meta$gsm), identical(meta$gsm, 
    colnames(count_mat)))

meta$group <- factor(meta$group, 
    levels = group_levels)

meta$disease <- factor(ifelse(meta$group == 
    "Control", "Control", "DKD"), 
    levels = c("Control", "DKD"))

meta$stage_numeric <- unname(c(Control = 0, 
    Early_DN = 1, Advanced_DN = 2)[as.character(meta$group)])

rownames(meta) <- meta$gsm

stopifnot(!anyNA(meta[, c("group", 
    "disease", "stage_numeric")]), 
    identical(as.integer(table(meta$group)), c(9L, 6L, 21L)))

write_csv(meta, "sample_metadata_36.csv")

write_csv(meta_all[!meta_all$gsm %in% 
    meta$gsm, ], "GEO_samples_absent_from_count_matrix.csv")

min_reps <- min(table(meta$group))

keep <- rowSums(count_mat >= 10L) >= 
    min_reps

write_csv(data.frame(total_genes = nrow(count_mat), 
    retained_genes = sum(keep), 
    removed_genes = sum(!keep), 
    rule = paste0("count >= 10 in at least ", 
        min_reps, " samples")), 
    "gene_filter_summary.csv")

fit_files <- file.path(obj_dir, 
    c("dds_group.rds", "dds_allDKD.rds", 
        "vst_object.rds"))

if (reuse_fit && all(file.exists(fit_files))) {
    dds <- readRDS(fit_files[1])
    dds_all <- readRDS(fit_files[2])
    vsd <- readRDS(fit_files[3])
    message("Reusing fitted DESeq2/VST objects because DKD_REUSE_FIT=1.")
} else {
    dds <- DESeqDataSetFromMatrix(count_mat[keep, 
        ], meta, design = ~group)
    dds_all <- DESeqDataSetFromMatrix(count_mat[keep, 
        ], meta, design = ~disease)
    dds <- DESeq(dds, quiet = TRUE)
    dds_all <- DESeq(dds_all, quiet = TRUE)
    vsd <- vst(dds, blind = FALSE)
    saveRDS(dds, fit_files[1], 
        compress = "xz")
    saveRDS(dds_all, fit_files[2], 
        compress = "xz")
    saveRDS(vsd, fit_files[3], 
        compress = "xz")
}

vst_mat <- assay(vsd)

symbol <- mapIds(org.Hs.eg.db, 
    keys = rownames(dds), keytype = "ENTREZID", 
    column = "SYMBOL", multiVals = "first")

annotation <- data.frame(entrez_id = rownames(dds), 
    symbol = unname(symbol), mean_normalized_count = rowMeans(counts(dds, 
        normalized = TRUE)))

write_csv(annotation, "gene_annotation.csv")

get_de <- function(object, variable, 
    numerator, denominator, label) {
    raw <- results(object, contrast = c(variable, 
        numerator, denominator), 
        alpha = 0.05)
    shr <- suppressMessages(lfcShrink(object, 
        contrast = c(variable, 
            numerator, denominator), 
        res = raw, type = "normal"))
    out <- data.frame(entrez_id = rownames(raw), 
        symbol = unname(symbol[rownames(raw)]), 
        baseMean = raw$baseMean, 
        log2FC_MLE = raw$log2FoldChange, 
        log2FC_shrunk = shr$log2FoldChange, 
        SE = raw$lfcSE, statistic = raw$stat, 
        p_value = raw$pvalue, FDR_BH = raw$padj, 
        contrast = label)
    out$significant <- !is.na(out$FDR_BH) & 
        out$FDR_BH < 0.05
    out$large_effect <- out$significant & 
        abs(out$log2FC_shrunk) >= 
            1
    out[order(out$FDR_BH, -abs(out$log2FC_shrunk), 
        na.last = TRUE), ]
}

de_specs <- list(All_DKD_vs_Control = list(dds_all, 
    "disease", "DKD", "Control"), 
    Early_DN_vs_Control = list(dds, 
        "group", "Early_DN", "Control"), 
    Advanced_DN_vs_Control = list(dds, 
        "group", "Advanced_DN", "Control"))

de <- lapply(names(de_specs), function(nm) {
    x <- de_specs[[nm]]
    get_de(x[[1]], x[[2]], x[[3]], 
        x[[4]], nm)
})

names(de) <- names(de_specs)

for (nm in names(de)) write_csv(de[[nm]], 
    paste0("DESeq2_", nm, ".csv"))

gene_var <- apply(vst_mat, 1, var)

top_genes <- names(sort(gene_var, 
    decreasing = TRUE))[seq_len(min(2000, 
    length(gene_var)))]

pca <- prcomp(t(vst_mat[top_genes, 
    ]), center = TRUE, scale. = FALSE)

var_exp <- 100 * pca$sdev^2/sum(pca$sdev^2)

pca_df <- cbind(meta, as.data.frame(pca$x[, 
    1:10, drop = FALSE]))

write_csv(data.frame(PC = paste0("PC", 
    seq_along(var_exp)), variance_percent = var_exp), 
    "PCA_variance_explained.csv")

sample_cor <- cor(vst_mat[top_genes, 
    ], method = "pearson")

median_cor <- (rowSums(sample_cor) - 
    1)/(ncol(sample_cor) - 1)

low_flag <- function(x) {
    cutoff <- median(x) - 4 * mad(x, 
        constant = 1.4826)
    if (!is.finite(cutoff)) 
        rep(FALSE, length(x))
    else x < cutoff
}

qc <- data.frame(gsm = meta$gsm, 
    subject_id = meta$subject_id, 
    group = meta$group, library_size = colSums(count_mat), 
    detected_genes = colSums(count_mat > 
        0), size_factor = sizeFactors(dds), 
    median_sample_correlation = median_cor[meta$gsm])

qc$flag_low_library <- low_flag(log10(qc$library_size))

qc$flag_low_detected <- low_flag(qc$detected_genes)

qc$flag_low_correlation <- low_flag(qc$median_sample_correlation)

qc$qc_flag <- with(qc, flag_low_library | 
    flag_low_detected | flag_low_correlation)

qc$excluded <- FALSE

write_csv(qc, "sample_QC_no_automatic_exclusion.csv")

write_csv(pca_df, "PCA_coordinates.csv")

modules <- list(dynamics_primary = c("DNM1L", 
    "MFF", "FIS1", "MIEF1", "MIEF2", 
    "MFN1", "MFN2", "OPA1"), fission_core = c("DNM1L", 
    "MFF", "FIS1"), fusion_core = c("MFN1", 
    "MFN2", "OPA1"), mitochondrial_function = c("NDUFS1", 
    "NDUFV1", "SDHA", "SDHB", "UQCRC1", 
    "UQCRC2", "COX4I1", "COX5A", 
    "ATP5F1A", "ATP5F1B", "CS", 
    "TFAM", "PPARGC1A", "NRF1", 
    "ESRRA"), ampk_transcript_proxy = c("PRKAA1", 
    "PRKAA2", "PRKAB1", "PRKAB2", 
    "PRKAG1", "PRKAG2"), mam_context = c("PACS2", 
    "VAPB", "RMDN3", "ITPR1", "ITPR2", 
    "ITPR3", "VDAC1", "MCU", "MICU1"))

target_genes <- unique(c("PIEZO1", 
    "CAMK2D", unlist(modules)))

valid <- !is.na(annotation$symbol) & 
    annotation$symbol != ""

best <- annotation[valid, ]

best <- best[order(best$mean_normalized_count, 
    decreasing = TRUE), ]

best <- best[!duplicated(best$symbol), 
    ]

symbol_expr <- vst_mat[match(best$entrez_id, 
    rownames(vst_mat)), , drop = FALSE]

rownames(symbol_expr) <- best$symbol

coverage <- do.call(rbind, lapply(names(modules), 
    function(nm) data.frame(module = nm, 
        gene = modules[[nm]], present = modules[[nm]] %in% 
            rownames(symbol_expr))))

coverage <- rbind(data.frame(module = "single_gene", 
    gene = c("PIEZO1", "CAMK2D"), 
    present = c("PIEZO1", "CAMK2D") %in% 
        rownames(symbol_expr)), 
    coverage)

write_csv(coverage, "module_gene_coverage.csv")

stopifnot("PIEZO1" %in% rownames(symbol_expr), 
    all(coverage$present[coverage$module %in% 
        c("fission_core", "fusion_core")]))

Control <- meta$group == "Control"

Control_mean <- rowMeans(symbol_expr[, 
    Control, drop = FALSE])

Control_sd <- apply(symbol_expr[, 
    Control, drop = FALSE], 1, 
    sd)

Control_sd[Control_sd == 0 | !is.finite(Control_sd)] <- NA_real_

gene_z <- sweep(sweep(symbol_expr, 
    1, Control_mean, "-"), 1, Control_sd, 
    "/")

scores <- data.frame(meta[, c("gsm", 
    "subject_id", "group", "disease", 
    "stage_numeric")], 
    PIEZO1_z = as.numeric(gene_z["PIEZO1", 
        ]))

for (nm in names(modules)) {
    genes <- intersect(modules[[nm]], 
        rownames(gene_z))
    scores[[paste0(nm, "_signed")]] <- colMeans(gene_z[genes, 
        , drop = FALSE], na.rm = TRUE)
    scores[[paste0(nm, "_distance")]] <- sqrt(colMeans(gene_z[genes, 
        , drop = FALSE]^2, na.rm = TRUE))
}

write_csv(scores, "DKD_Control_reference_module_scores.csv")

target_de <- do.call(rbind, lapply(de, 
    function(x) x[x$symbol %in% 
        target_genes, ]))

write_csv(target_de, "target_gene_DESeq2_results.csv")

outcomes <- c("PIEZO1_z", "dynamics_primary_distance", 
    "fission_core_signed", "fusion_core_signed", 
    "mitochondrial_function_signed", 
    "ampk_transcript_proxy_signed", 
    "mam_context_signed")

group_summary <- do.call(rbind, 
    lapply(outcomes, function(y) do.call(rbind, 
        lapply(group_levels, function(g) {
            z <- scores[scores$group == 
                g, y]
            data.frame(outcome = y, 
                group = g, n = sum(is.finite(z)), 
                mean = mean(z, 
                  na.rm = TRUE), 
                SD = sd(z, na.rm = TRUE), 
                median = median(z, 
                  na.rm = TRUE), 
                Q1 = quantile(z, 
                  0.25, na.rm = TRUE), 
                Q3 = quantile(z, 
                  0.75, na.rm = TRUE))
        }))))

write_csv(group_summary, "module_score_group_summary.csv")

run_group_test <- function(y) {
    fit_group <- lm(reformulate("group", 
        y), data = scores)
    a <- robust_table(fit_group)
    a <- a[a$term %in% paste0("group", 
        c("Early_DN", "Advanced_DN")), 
        ]
    a$contrast <- paste0(sub("^group", 
        "", a$term), "_vs_Control")
    fit_all <- lm(reformulate("disease", 
        y), data = scores)
    b <- robust_table(fit_all)
    b <- b[b$term == "diseaseDKD", 
        ]
    b$contrast <- "All_DKD_vs_Control"
    out <- rbind(a, b)
    out$outcome <- y
    out[, c("outcome", "contrast", 
        "beta", "robust_SE", "CI_low", 
        "CI_high", "statistic", 
        "p_value")]
}

group_tests <- do.call(rbind, lapply(outcomes, 
    run_group_test))

group_tests$FDR_BH <- p.adjust(group_tests$p_value, 
    method = "BH")

group_tests$testing_family <- "7 outcomes x 3 prespecified DKD-vs-Control contrasts"

write_csv(group_tests, "adjusted_group_comparisons_HC3.csv")

stage_trends <- do.call(rbind, lapply(outcomes, function(y) {
    test <- suppressWarnings(cor.test(scores$stage_numeric, 
        scores[[y]], method = "spearman", exact = FALSE))
    data.frame(outcome = y, n = nrow(scores), 
        rho = unname(test$estimate), p_value = test$p.value)
}))

stage_trends$FDR_BH <- p.adjust(stage_trends$p_value, 
    method = "BH")

stage_trends$testing_family <- "7 prespecified ordered-stage outcomes"

write_csv(stage_trends, "ordered_stage_Spearman.csv")

partial_spearman <- function(data, 
    x, y, covariates, subgroup) {
    d <- droplevels(data[, unique(c("gsm", 
        "group", x, y, covariates))])
    n_before <- nrow(d)
    d <- d[rowSums(is.na(d)) == 
        0, ]
    n_removed <- n_before - nrow(d)
    zmat <- model.matrix(reformulate(covariates), 
        d)
    rx <- rank(d[[x]], ties.method = "average")
    ry <- rank(d[[y]], ties.method = "average")
    ex <- lm.fit(zmat, rx)$residuals
    ey <- lm.fit(zmat, ry)$residuals
    rho <- cor(ex, ey)
    rank_z <- qr(zmat)$rank
    dfree <- nrow(d) - rank_z - 
        1
    statistic <- rho * sqrt(dfree/(1 - 
        rho^2))
    p <- 2 * pt(-abs(statistic), 
        df = dfree)
    se_z <- 1/sqrt(nrow(d) - rank_z - 
        2)
    ci <- tanh(atanh(rho) + c(-1, 
        1) * qnorm(0.975) * se_z)
    list(result = data.frame(subgroup = subgroup, 
        x = x, outcome = y, n = nrow(d), 
        missing_removed = n_removed, 
        covariate_rank = rank_z - 
            1, rho = rho, CI_low = ci[1], 
        CI_high = ci[2], statistic = statistic, 
        df = dfree, p_value = p), 
        residuals = data.frame(gsm = d$gsm, 
            group = d$group, subgroup = subgroup, 
            outcome = y, PIEZO1_rank_residual = ex, 
            outcome_rank_residual = ey))
}

assoc_outcomes <- c("dynamics_primary_distance", 
    "fission_core_signed", "fusion_core_signed", 
    "mitochondrial_function_signed", 
    "ampk_transcript_proxy_signed", 
    "mam_context_signed")

subsets <- list(All_DKD = list(data = scores[scores$group != 
    "Control", ], cov = "stage_numeric"))

assoc_list <- unlist(lapply(names(subsets), 
    function(sg) lapply(assoc_outcomes, 
        function(y) partial_spearman(subsets[[sg]]$data, 
            "PIEZO1_z", y, subsets[[sg]]$cov, 
            sg))), recursive = FALSE)

assoc <- do.call(rbind, lapply(assoc_list, 
    function(x) x[["result"]]))

assoc_resid <- do.call(rbind, lapply(assoc_list, 
    function(x) x[["residuals"]]))

assoc$FDR_BH <- p.adjust(assoc$p_value, 
    method = "BH")

assoc$testing_family <- "6 prespecified PIEZO1-module associations within DKD"

write_csv(assoc, "within_DKD_partial_Spearman.csv")

write_csv(assoc_resid, "within_DKD_adjusted_rank_residuals.csv")

dkd <- droplevels(scores[scores$group != 
    "Control", ])

piezo_cutoff <- median(dkd$PIEZO1_z, 
    na.rm = TRUE)

dkd$PIEZO1_level <- factor(ifelse(dkd$PIEZO1_z <= 
    piezo_cutoff, "PIEZO1 low", 
    "PIEZO1 high"), levels = c("PIEZO1 low", 
    "PIEZO1 high"))

int_outcomes <- c("fission_core_signed", 
    "fusion_core_signed", "dynamics_primary_distance", 
    "mitochondrial_function_signed")

int_results <- lapply(int_outcomes, 
    function(y) {
        form <- as.formula(paste(y, 
            "~ ampk_transcript_proxy_signed * PIEZO1_z +", 
            "stage_numeric"))
        fit <- lm(form, data = dkd)
        tb <- robust_table(fit)
        tb$outcome <- y
        list(fit = fit, full = tb, 
            primary = tb[tb$term == 
                "ampk_transcript_proxy_signed:PIEZO1_z", 
                ])
    })

names(int_results) <- int_outcomes

int_full <- do.call(rbind, lapply(int_results, 
    function(x) x[["full"]]))

int_primary <- do.call(rbind, lapply(int_results, 
    function(x) x[["primary"]]))

int_primary$FDR_BH <- p.adjust(int_primary$p_value, 
    method = "BH")

int_primary$testing_family <- "AMPK x PIEZO1 interaction across 4 outcomes"

write_csv(int_full, "interaction_models_all_coefficients_HC3.csv")

write_csv(int_primary, "interaction_terms_HC3.csv")

ampk_grid <- seq(quantile(dkd$ampk_transcript_proxy_signed, 
    0.05), quantile(dkd$ampk_transcript_proxy_signed, 
    0.95), length.out = 80)

piezo_grid <- c(`PIEZO1 low` = mean(dkd$PIEZO1_z[dkd$PIEZO1_level == 
    "PIEZO1 low"], na.rm = TRUE), 
    `PIEZO1 high` = mean(dkd$PIEZO1_z[dkd$PIEZO1_level == 
        "PIEZO1 high"], na.rm = TRUE))

ampk_plot_outcomes <- c("fission_core_signed", 
    "fusion_core_signed", "dynamics_primary_distance", 
    "mitochondrial_function_signed")

make_ampk_predictions <- function(y) {
    fit <- int_results[[y]]$fit
    do.call(rbind, lapply(names(piezo_grid), 
        function(pl) do.call(rbind, 
            lapply(ampk_grid, function(a) {
                nd <- dkd
                nd$PIEZO1_z <- unname(piezo_grid[pl])
                nd$ampk_transcript_proxy_signed <- a
                pr <- average_prediction_hc3(fit, 
                  nd)
                data.frame(outcome = y, 
                  ampk_transcript_proxy_signed = a, 
                  PIEZO1_level = pl, 
                  PIEZO1_z_value = unname(piezo_grid[pl]), 
                  estimate = unname(pr["estimate"]), 
                  robust_SE = unname(pr["robust_SE"]), 
                  CI_low = unname(pr["CI_low"]), 
                  CI_high = unname(pr["CI_high"]))
            }))))
}

ampk_pred <- do.call(rbind, lapply(ampk_plot_outcomes, 
    make_ampk_predictions))

write_csv(ampk_pred, "AMPKxPIEZO1_average_marginal_predictions_HC3.csv")

disease_int_outcomes <- c("dynamics_primary_distance", 
    "fission_core_signed", "fusion_core_signed", 
    "mitochondrial_function_signed")

disease_int_results <- lapply(disease_int_outcomes, 
    function(y) {
        form <- as.formula(paste(y, 
            "~ PIEZO1_z * disease"))
        fit <- lm(form, data = scores)
        tb <- robust_table(fit)
        tb$outcome <- y
        interaction_term <- tb$term[grepl("PIEZO1_z:diseaseDKD|diseaseDKD:PIEZO1_z", 
            tb$term)]
        stopifnot(length(interaction_term) == 
            1L)
        list(fit = fit, full = tb, 
            primary = tb[tb$term == 
                interaction_term, 
                ])
    })

names(disease_int_results) <- disease_int_outcomes

disease_int_full <- do.call(rbind, 
    lapply(disease_int_results, 
        function(x) x[["full"]]))

disease_int_primary <- do.call(rbind, 
    lapply(disease_int_results, 
        function(x) x[["primary"]]))

disease_int_primary$FDR_BH <- p.adjust(disease_int_primary$p_value, 
    method = "BH")

disease_int_primary$testing_family <- "PIEZO1 x DKD interaction across 4 prespecified outcomes"

write_csv(disease_int_full, "disease_interaction_all_coefficients_HC3.csv")

write_csv(disease_int_primary, 
    "disease_interaction_terms_HC3.csv")

simple_slopes_hc3 <- function(fit, 
    outcome) {
    beta <- coef(fit)
    robust_vcov <- vcovHC(fit, 
        type = "HC3")
    interaction_term <- names(beta)[grepl("PIEZO1_z:diseaseDKD|diseaseDKD:PIEZO1_z", 
        names(beta))]
    stopifnot(length(interaction_term) == 
        1L, !anyNA(beta))
    control_l <- setNames(rep(0, 
        length(beta)), names(beta))
    control_l["PIEZO1_z"] <- 1
    dkd_l <- control_l
    dkd_l[interaction_term] <- 1
    slope_row <- function(lvec, 
        group_label) {
        estimate <- drop(lvec %*% 
            beta)
        robust_se <- sqrt(drop(t(lvec) %*% 
            robust_vcov %*% lvec))
        statistic <- estimate/robust_se
        p_value <- 2 * pt(-abs(statistic), 
            df = df.residual(fit))
        crit <- qt(0.975, df.residual(fit))
        data.frame(outcome = outcome, 
            group = group_label, 
            beta = estimate, robust_SE = robust_se, 
            CI_low = estimate - 
                crit * robust_se, 
            CI_high = estimate + 
                crit * robust_se, 
            statistic = statistic, 
            df = df.residual(fit), 
            p_value = p_value)
    }
    rbind(slope_row(control_l, "Control"), 
        slope_row(dkd_l, "DKD"))
}

disease_simple_slopes <- do.call(rbind, 
    lapply(disease_int_outcomes, 
        function(y) simple_slopes_hc3(disease_int_results[[y]]$fit, 
            y)))

disease_simple_slopes$FDR_BH <- p.adjust(disease_simple_slopes$p_value, 
    method = "BH")

disease_simple_slopes$testing_family <- "4 outcomes x 2 secondary simple slopes"

write_csv(disease_simple_slopes, 
    "disease_interaction_simple_slopes_HC3.csv")

disease_group_ranges <- sapply(levels(scores$disease), 
    function(status) quantile(scores$PIEZO1_z[scores$disease == 
        status], c(0.025, 0.975)))

disease_common_low <- max(disease_group_ranges[1, 
    ])

disease_common_high <- min(disease_group_ranges[2, 
    ])

stopifnot(disease_common_low < 
    disease_common_high)

disease_piezo_grid <- seq(disease_common_low, 
    disease_common_high, length.out = 100)

make_disease_predictions <- function(y) {
    fit <- disease_int_results[[y]]$fit
    do.call(rbind, lapply(levels(scores$disease), 
        function(status) do.call(rbind, 
            lapply(disease_piezo_grid, 
                function(piezo_value) {
                  nd <- scores
                  nd$PIEZO1_z <- piezo_value
                  nd$disease <- factor(status, 
                    levels = levels(scores$disease))
                  pr <- average_prediction_hc3(fit, 
                    nd)
                  data.frame(outcome = y, 
                    disease = status, 
                    display_group = ifelse(status == 
                      "Control", 
                      "Control", 
                      "DKD"), PIEZO1_z = piezo_value, 
                    estimate = unname(pr["estimate"]), 
                    robust_SE = unname(pr["robust_SE"]), 
                    CI_low = unname(pr["CI_low"]), 
                    CI_high = unname(pr["CI_high"]), 
                    n_group = sum(scores$disease == 
                      status))
                }))))
}

disease_pred <- do.call(rbind, 
    lapply(disease_int_outcomes, 
        make_disease_predictions))

disease_pred$display_group <- factor(disease_pred$display_group, 
    levels = c("Control", "DKD"))

write_csv(disease_pred, "disease_interaction_average_marginal_predictions_HC3.csv")

qc_plot <- merge(qc, pca_df[, c("gsm", 
    "PC1", "PC2")], by = "gsm", 
    sort = FALSE)

qc_plot$group <- factor(qc_plot$group, 
    levels = group_levels)

qc_plot$library_million <- qc_plot$library_size/1e+06

qc_plot$detected_thousand <- qc_plot$detected_genes/1000

plot_data <- list(meta = meta, 
    group_levels = group_levels, 
    group_labels = group_labels, 
    group_cols = group_cols, qc = qc, 
    qc_plot = qc_plot, pca_df = pca_df, 
    var_exp = var_exp, sample_cor = sample_cor, 
    de = de, scores = scores, modules = modules, 
    symbol_expr = symbol_expr, 
    group_tests = group_tests, stage_trends = stage_trends, 
    assoc = assoc, assoc_resid = assoc_resid, 
    dkd = dkd, piezo_grid = piezo_grid, 
    ampk_pred = ampk_pred, int_primary = int_primary, 
    disease_pred = disease_pred, 
    disease_int_primary = disease_int_primary)

plot_data_file <- file.path(obj_dir, 
    "DKD_plot_data.rds")

saveRDS(plot_data, plot_data_file, 
    compress = "xz")

message("Saved plotting data: ", 
    plot_data_file)

write_csv(qc_plot, "Fig1_source_data.csv")

write_csv(scores, "Fig2_module_source_data.csv")

write_csv(assoc_resid[assoc_resid$subgroup == 
    "All_DKD", ], "Fig3_association_source_data.csv")

write_csv(ampk_pred, "Fig4_AMPK_interaction_predictions_source_data.csv")

write_csv(disease_pred, "Fig5_disease_interaction_predictions_source_data.csv")

manifest <- data.frame(item = c("count_file", 
    "count_file_md5", "metadata_file", 
    "metadata_file_md5", "samples", 
    "group_counts", "genes_raw", 
    "genes_retained", "sample_exclusions", 
    "DESeq2_design_group", "DESeq2_design_allDKD", 
    "AMPK_interaction_model", "disease_interaction_model", 
    "interaction_variance", "score_scale", 
    "analysis_date"), value = c(count_file, 
    unname(tools::md5sum(count_file)), 
    meta_file, unname(tools::md5sum(meta_file)), 
    ncol(count_mat), paste(names(table(meta$group)), 
        table(meta$group), collapse = "; "), 
    nrow(count_mat), sum(keep), 
    "0 (QC flags retained)", "~ group", 
    "~ disease", 
    "~ ampk_transcript_proxy_signed * PIEZO1_z + stage_numeric", 
    "~ PIEZO1_z * disease", 
    "HC3 heteroskedasticity-robust standard errors and 95% confidence intervals", 
    "VST expression; per-gene Control-reference z", 
    as.character(Sys.Date())))

write.table(manifest, file.path(tab_dir, 
    "analysis_manifest.tsv"), sep = "\t", 
    row.names = FALSE, quote = FALSE)

writeLines(capture.output(sessionInfo()), 
    file.path(project, "logs/sessionInfo.txt"))

key_assoc <- assoc[assoc$subgroup == 
    "All_DKD" & assoc$outcome %in% 
    c("dynamics_primary_distance", 
        "fission_core_signed"), 
    ]

key_group <- group_tests[group_tests$contrast == 
    "All_DKD_vs_Control" & group_tests$outcome %in% 
    c("PIEZO1_z", "dynamics_primary_distance"), 
    ]

key_ampk <- int_primary[int_primary$outcome %in% 
    ampk_plot_outcomes, ]

key_disease <- disease_int_primary[disease_int_primary$outcome %in% 
    disease_int_outcomes, ]

summary_lines <- c("# DKD complete analysis run summary", 
    "", sprintf("- Samples: %d (%s).", 
        nrow(meta), paste(names(table(meta$group)), 
            table(meta$group), 
            collapse = ", ")), 
    sprintf("- Genes: %d raw; %d retained by the prespecified count filter.", 
        nrow(count_mat), sum(keep)), 
    sprintf("- QC: %d samples flagged diagnostically; 0 excluded.", 
        sum(qc$qc_flag)), "- Primary differential-expression model: DESeq2 with Control as the reference group.", 
    "- Module scale: VST expression converted to Control-reference per-gene z scores.", 
    paste("- Within-DKD association: partial Spearman correlation adjusted for", 
        "disease stage (Early versus Advanced)."), 
    paste("- Signal-coupling model: AMPK transcript proxy x PIEZO1 within DKD,", 
        "with HC3 robust inference."), 
    paste("- Disease-modification model: PIEZO1 x DKD using Control and DKD samples,", 
        "with HC3 robust inference; stage is not entered because it is undefined in controls."), 
    "", "## Key prespecified results", 
    "", paste0("- ", key_group$outcome, 
        ": all DKD vs Control beta = ", 
        sprintf("%.3f", key_group$beta), 
        ", 95% CI ", sprintf("[%.3f, %.3f]", 
            key_group$CI_low, key_group$CI_high), 
        ", BH q = ", fmt_p(key_group$FDR_BH), 
        "."), paste0("- AMPK x PIEZO1 interaction for ", 
        key_ampk$outcome, ": beta = ", 
        sprintf("%.3f", key_ampk$beta), 
        ", 95% CI ", sprintf("[%.3f, %.3f]", 
            key_ampk$CI_low, key_ampk$CI_high), 
        ", BH q = ", fmt_p(key_ampk$FDR_BH), 
        "."), paste0("- PIEZO1 x DKD interaction for ", 
        key_disease$outcome, ": beta = ", 
        sprintf("%.3f", key_disease$beta), 
        ", 95% CI ", sprintf("[%.3f, %.3f]", 
            key_disease$CI_low, 
            key_disease$CI_high), 
        ", BH q = ", fmt_p(key_disease$FDR_BH), 
        "."), paste0("- PIEZO1 vs ", 
        key_assoc$outcome, " within DKD: partial rho = ", 
        sprintf("%.3f", key_assoc$rho), 
        ", 95% CI ", sprintf("[%.3f, %.3f]", 
            key_assoc$CI_low, key_assoc$CI_high), 
        ", BH q = ", fmt_p(key_assoc$FDR_BH), 
        "."), "", "## Boundary", 
    "", paste("This is count-matrix-level QC;", 
        "FASTQ/BAM metrics such as mapping rate,", 
        "duplication, insert size and RIN are unavailable."), 
    paste("Associations are not evidence of causality.", 
        "This integrated V1 plus disease-interaction analysis is not forced to match the earlier", 
        "processed voom/SVA result."))

writeLines(summary_lines, file.path(project, 
    "DKD_complete_analysis_summary.md"), 
    useBytes = TRUE)

message("DKD complete analysis finished: ", 
    project)



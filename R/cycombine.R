# processor-cycombine: cyCombine batch correction (Pedersen et al., Nat Commun 2022) of
# arcsinh-transformed FCS files.
#
# cyCombine's workflow on one table of cells from every sample:
#   1. normalise each batch on its own (NORM_METHOD: scale = z-score per batch, rank, or
#      qnorm) - only to cluster cells by type, not by batch;
#   2. cluster all cells with a self-organising map (kohonen, XDIM x YDIM nodes);
#   3. in every SOM node, run ComBat (sva, parametric empirical Bayes) on the *un-normalised*
#      marker values with batch as the batch variable and COVARIATE_COLUMN (default
#      condition) in the model, so condition differences are kept; values are capped to the
#      node's input range.
# ComBat fits and applies its model in one call on the cells it is given: cyCombine has no
# way to train on some cells and apply the correction to others. So the app corrects a
# seeded subsample of at most MAX_CELLS_PER_SAMPLE events per sample and writes those
# events (all channels, marker values corrected): the output files are subsampled.
# Afterwards cyCombine's own evaluation: SOM on the corrected data, then per node, marker
# and batch pair the earth mover's distance (EMD) before/after, and per node, batch and
# marker the median absolute deviation (MAD) before/after.

CYCOMBINE_PKGS <- c("cyCombine", "sva", "kohonen", "flowCore", "emdist", "ggplot2", "ggridges", "cowplot", "jsonlite")
PLOT_CELLS_PER_SAMPLE <- 5000L
EVAL_CELLS_PER_SAMPLE <- 20000L

# cgroup v2 memory of this container (current / peak), for the logs.
mem_status <- function() {
  rd <- function(f) tryCatch(as.numeric(readLines(f, warn = FALSE)[1]) / 2^30, error = function(e) NA)
  sprintf("%.2f GiB now, %.2f GiB peak", rd("/sys/fs/cgroup/memory.current"), rd("/sys/fs/cgroup/memory.peak"))
}

cycombine_params <- function() {
  cov <- env_str("COVARIATE_COLUMN", "condition")
  if (tolower(cov) %in% c("none", "-")) cov <- ""
  list(
    BATCH_COLUMN = env_str("BATCH_COLUMN", "batch"),
    COVARIATE_COLUMN = cov,
    XDIM = env_int("XDIM", 8, min = 1, max = 50),
    YDIM = env_int("YDIM", 8, min = 1, max = 50),
    NORM_METHOD = env_choice("NORM_METHOD", "scale", c("scale", "rank", "qnorm")),
    # "auto" (default): as many as fit the memory budget, at most 50000 per sample; see run().
    MAX_CELLS_PER_SAMPLE = if (tolower(env_str("MAX_CELLS_PER_SAMPLE", "auto")) == "auto") NA_integer_
                           else env_int("MAX_CELLS_PER_SAMPLE", 50000, min = 0),
    MAX_VALUES = env_int("MAX_VALUES", 30000000, min = 1e5),
    MARKERS = env_list("MARKERS"),
    SEED = env_int("SEED", 42)
  )
}

norm_keyword <- function(p, n_in, n_out) {
  paste(c("cycombine",
          paste0("covariate=", if (nzchar(p$COVARIATE_COLUMN)) p$COVARIATE_COLUMN else "none"),
          paste0("som=", p$XDIM, "x", p$YDIM), paste0("norm=", p$NORM_METHOD),
          paste0("subsample=", if (n_out < n_in) paste0(n_out, "/", n_in) else "none"),
          paste0("seed=", p$SEED)), collapse = ";")
}

seeded_rows <- function(n, k, seed) {
  if (k <= 0 || n <= k) return(seq_len(n))
  set.seed(seed)
  sort(sample.int(n, k))
}

# compute_emd() returns cluster -> marker -> batch x batch matrix (upper triangle): tidy it.
tidy_emd <- function(emd) {
  out <- list()
  for (cl in names(emd)) for (mk in names(emd[[cl]])) {
    d <- emd[[cl]][[mk]]
    idx <- which(!is.na(d), arr.ind = TRUE)
    if (!nrow(idx)) next
    out[[length(out) + 1]] <- data.frame(cluster = cl, marker = mk, batch1 = rownames(d)[idx[, 1]],
                                         batch2 = colnames(d)[idx[, 2]], emd = d[idx], stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

tidy_mad <- function(mad) {
  out <- list()
  for (cl in names(mad)) for (b in names(mad[[cl]])) {
    v <- mad[[cl]][[b]]
    if (!length(v)) next
    out[[length(out) + 1]] <- data.frame(cluster = cl, batch = b, marker = names(v), mad = unname(v),
                                         stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

run_cycombine <- function(input_dir, output_dir, app_dir, p = cycombine_params()) {
  t0 <- Sys.time()
  log_line("processor-cycombine ", app_version(app_dir), " COVARIATE_COLUMN=",
           if (nzchar(p$COVARIATE_COLUMN)) p$COVARIATE_COLUMN else "(none)", " SOM=", p$XDIM, "x", p$YDIM,
           " NORM_METHOD=", p$NORM_METHOD, " MAX_CELLS_PER_SAMPLE=", p$MAX_CELLS_PER_SAMPLE, " SEED=", p$SEED)
  suppressPackageStartupMessages({
    library(flowCore)
    library(cyCombine)
  })
  `%>%` <- magrittr::`%>%`

  # ---- inputs ----
  fcs <- find_fcs(input_dir, "fcs")
  samples_csv <- find_named(input_dir, "samples.csv")
  panel_csv <- find_named(input_dir, "panel.csv")
  if (!length(fcs)) stop("no .fcs files found under INPUT_DIR")
  if (is.na(samples_csv)) stop("samples.csv not found under INPUT_DIR")
  if (is.na(panel_csv)) stop("panel.csv not found under INPUT_DIR")
  needed <- c(p$BATCH_COLUMN, if (nzchar(p$COVARIATE_COLUMN)) p$COVARIATE_COLUMN)
  samples <- read_samples(samples_csv, needed)
  panel <- read_panel(panel_csv)
  sel <- select_markers(panel, p$MARKERS)
  channels <- sel$channel
  markers <- sel$marker
  if (anyDuplicated(markers)) markers <- channels  # column names must be unique
  log_line(sprintf("inputs: %d FCS files, samples.csv %s, panel.csv %s; %d markers corrected",
                   length(fcs), rel_path(samples_csv, input_dir), rel_path(panel_csv, input_dir), length(channels)))

  m <- match_samples(fcs, samples)
  m$batch <- ifelse(is.na(m$row), NA_character_, samples[[p$BATCH_COLUMN]][m$row])
  m$covar <- if (nzchar(p$COVARIATE_COLUMN)) ifelse(is.na(m$row), NA_character_, samples[[p$COVARIATE_COLUMN]][m$row]) else "all"
  status <- list()
  for (i in which(is.na(m$row))) {
    status[[m$sample_id[i]]] <- list(status = "failed", reason = "not in samples.csv")
    log_line(m$sample_id[i], ": FAILED not in samples.csv")
  }
  missing_fcs <- setdiff(samples$sample_id, m$sample_id)
  if (length(missing_fcs)) {
    log_line(sprintf("%d sample(s) in samples.csv have no FCS file (dropped upstream?): %s",
                     length(missing_fcs), paste(missing_fcs, collapse = ", ")))
  }
  m <- m[!is.na(m$row), , drop = FALSE]
  bad <- is.na(m$batch) | !nzchar(m$batch) | m$batch == "NA"
  for (i in which(bad)) {
    status[[m$sample_id[i]]] <- list(status = "failed", reason = paste0("no ", p$BATCH_COLUMN))
    log_line(m$sample_id[i], ": FAILED empty ", p$BATCH_COLUMN)
  }
  m <- m[!bad, , drop = FALSE]
  m <- m[order(m$batch, m$sample_id), , drop = FALSE]
  m$idx <- seq_len(nrow(m))

  # ComBat holds every corrected value in memory, so memory grows with cells x markers: the
  # spectral demo's 50000 cells x 24 samples x 25 markers (30 M values) peaked at 5.3 GB. "auto"
  # keeps cells x markers within MAX_VALUES (12400 per sample for 48 CyTOF samples x 42 markers).
  if (is.na(p$MAX_CELLS_PER_SAMPLE)) {
    p$MAX_CELLS_PER_SAMPLE <- as.integer(min(50000, floor(p$MAX_VALUES / (length(channels) * max(1, nrow(m))))))
    log_line(sprintf("MAX_CELLS_PER_SAMPLE=auto: %d per sample (%d samples x %d markers within MAX_VALUES=%d)",
                     p$MAX_CELLS_PER_SAMPLE, nrow(m), length(channels), p$MAX_VALUES))
  }

  # ---- read: seeded subsample of each sample, marker columns only ----
  t_read <- Sys.time()
  parts <- list()
  m$events <- NA_integer_
  m$used <- NA_integer_
  for (r in seq_len(nrow(m))) {
    sid <- m$sample_id[r]
    t1 <- Sys.time()
    res <- tryCatch({
      f <- read_fcs(m$path[r])
      miss <- setdiff(channels, colnames(f$exprs))
      if (length(miss)) stop("missing channel(s) ", paste(miss, collapse = ", "))
      rows <- seeded_rows(nrow(f$exprs), p$MAX_CELLS_PER_SAMPLE, p$SEED + m$idx[r])
      x <- f$exprs[rows, channels, drop = FALSE]
      if (nrow(x) < 10) stop("fewer than 10 events")
      colnames(x) <- markers
      d <- tibble::as_tibble(as.data.frame(x))
      d$sample <- sid
      d$batch <- m$batch[r]
      d$condition <- m$covar[r]
      list(d = d, n = nrow(f$exprs))
    }, error = function(e) e)
    if (inherits(res, "error")) {
      status[[sid]] <- list(status = "failed", reason = conditionMessage(res))
      log_line(sprintf("read %s (%s): FAILED %s", sid, m$batch[r], conditionMessage(res)))
      next
    }
    parts[[sid]] <- res$d
    m$events[r] <- res$n
    m$used[r] <- nrow(res$d)
    log_line(sprintf("read %s (%s, %s): %d of %d events (%.1fs)", sid, m$batch[r], m$covar[r], nrow(res$d), res$n,
                     as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  }
  m <- m[m$sample_id %in% names(parts), , drop = FALSE]
  if (!nrow(m)) stop("no usable samples")
  if (length(unique(m$batch)) < 2) stop("need samples from at least 2 batches")
  uncorrected <- dplyr::bind_rows(parts[m$sample_id])
  rm(parts)
  uncorrected$id <- seq_len(nrow(uncorrected))
  covar <- if (nzchar(p$COVARIATE_COLUMN) && length(unique(uncorrected$condition)) > 1) "condition" else NULL
  if (nzchar(p$COVARIATE_COLUMN) && is.null(covar)) log_line("covariate has a single value; ComBat runs without it")
  log_line(sprintf("pooled %d cells from %d samples in %d batches (%.1fs)", nrow(uncorrected), nrow(m),
                   length(unique(m$batch)), as.numeric(difftime(Sys.time(), t_read, units = "secs"))))
  print(table(batch = m$batch, covariate = m$covar))

  # ---- cyCombine: batch-wise normalisation -> SOM -> ComBat per SOM node ----
  t_cor <- Sys.time()
  labels <- suppressMessages(
    uncorrected %>%
      cyCombine::normalize(markers = markers, norm_method = p$NORM_METHOD, pb = FALSE) %>%
      cyCombine::create_som(markers = markers, seed = p$SEED, xdim = p$XDIM, ydim = p$YDIM, rlen = 10))
  log_line(sprintf("SOM %dx%d on %s-normalised data: %d nodes used (%.1fs); memory %s", p$XDIM, p$YDIM, p$NORM_METHOD,
                   length(unique(labels)), as.numeric(difftime(Sys.time(), t_cor, units = "secs")), mem_status()))
  t_cb <- Sys.time()
  corrected <- suppressMessages(cyCombine::correct_data(uncorrected, label = labels, markers = markers,
                                                        covar = covar, parametric = TRUE, pb = FALSE))
  corrected <- corrected[order(corrected$id), , drop = FALSE]
  stopifnot(identical(corrected$id, uncorrected$id))
  t_cb_s <- as.numeric(difftime(Sys.time(), t_cb, units = "secs"))
  log_line(sprintf("ComBat per SOM node done (%.1fs); memory %s", t_cb_s, mem_status()))
  t_cor_s <- as.numeric(difftime(Sys.time(), t_cor, units = "secs"))

  # ---- write: re-read each file, same rows, marker channels replaced ----
  out_fcs <- file.path(output_dir, "fcs")
  out_qc <- file.path(output_dir, "qc", "cycombine")
  dir.create(out_fcs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_qc, recursive = TRUE, showWarnings = FALSE)
  t_write <- Sys.time()
  cor_mat <- as.matrix(corrected[, markers])
  for (r in seq_len(nrow(m))) {
    sid <- m$sample_id[r]
    t1 <- Sys.time()
    res <- tryCatch({
      f <- read_fcs(m$path[r])
      rows <- seeded_rows(nrow(f$exprs), p$MAX_CELLS_PER_SAMPLE, p$SEED + m$idx[r])
      x <- f$exprs[rows, , drop = FALSE]
      sel_rows <- which(uncorrected$sample == sid)
      stopifnot(length(sel_rows) == nrow(x))
      before <- x[, channels, drop = FALSE]
      x[, channels] <- cor_mat[sel_rows, , drop = FALSE]
      write_fcs(file.path(out_fcs, paste0(sid, ".fcs")), x, f$keywords,
                c(PENNSIEVE_NORM = norm_keyword(p, nrow(f$exprs), nrow(x))))
      mean(abs(x[, channels] - before))
    }, error = function(e) e)
    secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (inherits(res, "error")) {
      status[[sid]] <- list(status = "failed", reason = conditionMessage(res))
      log_line(sprintf("%s (%s): FAILED %s (%.1fs)", sid, m$batch[r], conditionMessage(res), secs))
    } else {
      status[[sid]] <- list(status = "ok", batch = m$batch[r], covariate = m$covar[r], events_in = m$events[r],
                            events_out = m$used[r], mean_abs_change = round(res, 5), seconds = round(secs, 2))
      log_line(sprintf("%s (%s): wrote %d of %d events, mean |change| %.4f (%.1fs)", sid, m$batch[r],
                       m$used[r], m$events[r], res, secs))
    }
  }
  t_write_s <- as.numeric(difftime(Sys.time(), t_write, units = "secs"))
  rm(cor_mat)

  # ---- qc: cyCombine's evaluation (EMD reduction, MAD score) + density plots ----
  t_eval <- Sys.time()
  som_counts <- as.data.frame(table(sample_id = uncorrected$sample, som_node = labels), stringsAsFactors = FALSE)
  som_counts <- som_counts[som_counts$Freq > 0, ]
  names(som_counts)[names(som_counts) == "Freq"] <- "n"
  som_counts$batch <- m$batch[match(som_counts$sample_id, m$sample_id)]
  write_csv(som_counts[, c("sample_id", "batch", "som_node", "n")], file.path(out_qc, "som_node_counts.csv"))

  # Evaluate on at most EVAL_CELLS_PER_SAMPLE cells per sample (seeded): the EMD/MAD
  # functions copy the whole table per cluster and batch, which would double the memory of
  # the correction itself.
  # The cap follows the memory budget too: 40% of MAX_VALUES (20000 per sample on the spectral
  # demo, ~6000 for 48 CyTOF samples x 42 markers, whose evaluation was killed at 6 GB at 20000).
  eval_cap <- as.integer(min(EVAL_CELLS_PER_SAMPLE,
                             floor(0.4 * p$MAX_VALUES / (length(markers) * max(1, nrow(m))))))
  eval_rows <- unlist(lapply(split(seq_len(nrow(uncorrected)), uncorrected$sample), function(ix) {
    if (length(ix) <= eval_cap) return(ix)
    set.seed(p$SEED)
    sort(sample(ix, eval_cap))
  }), use.names = FALSE)
  eval_rows <- sort(eval_rows)
  keep_cols <- c("id", "sample", "batch", "condition", markers)
  uncorrected <- uncorrected[eval_rows, keep_cols]
  corrected <- as.data.frame(corrected)[eval_rows, keep_cols]
  corrected <- tibble::as_tibble(corrected)
  gc(verbose = FALSE)
  log_line(sprintf("evaluation on %d cells (at most %d per sample); memory %s", nrow(corrected),
                   eval_cap, mem_status()))
  eval_labels <- suppressMessages(cyCombine::create_som(corrected, markers = markers, seed = p$SEED,
                                                        xdim = p$XDIM, ydim = p$YDIM, rlen = 10))
  uncorrected$label <- eval_labels
  corrected$label <- eval_labels
  emd_u <- tidy_emd(suppressMessages(cyCombine::compute_emd(uncorrected, cell_col = "label", markers = markers)))
  emd_c <- tidy_emd(suppressMessages(cyCombine::compute_emd(corrected, cell_col = "label", markers = markers)))
  emd <- merge(emd_u, emd_c, by = c("cluster", "marker", "batch1", "batch2"), suffixes = c("_uncorrected", "_corrected"))
  emd <- emd[order(as.integer(emd$cluster), match(emd$marker, markers), emd$batch1, emd$batch2), ]
  # evaluate_emd(): drop pairs with EMD < 2 both before and after, reduction = sum of
  # (uncorrected - corrected) / sum of uncorrected.
  keep <- !(emd$emd_corrected < 2 & emd$emd_uncorrected < 2)
  emd_reduction <- if (any(keep)) sum(emd$emd_uncorrected[keep] - emd$emd_corrected[keep]) / sum(emd$emd_uncorrected[keep]) else NA
  write_csv(emd, file.path(out_qc, "emd.csv"))
  mad_u <- tidy_mad(suppressMessages(cyCombine::compute_mad(uncorrected, cell_col = "label", markers = markers)))
  mad_c <- tidy_mad(suppressMessages(cyCombine::compute_mad(corrected, cell_col = "label", markers = markers)))
  mad <- merge(mad_u, mad_c, by = c("cluster", "batch", "marker"), suffixes = c("_uncorrected", "_corrected"))
  mad <- mad[order(as.integer(mad$cluster), mad$batch, match(mad$marker, markers)), ]
  # evaluate_mad(): score = median over nodes, batches and markers of |MAD before - MAD after|.
  mad_score <- stats::median(abs(mad$mad_uncorrected - mad$mad_corrected), na.rm = TRUE)
  write_csv(mad, file.path(out_qc, "mad.csv"))
  per_marker <- stats::aggregate(cbind(emd_uncorrected, emd_corrected) ~ marker, data = emd, FUN = mean)
  per_marker <- per_marker[order(match(per_marker$marker, markers)), ]
  write_csv(per_marker, file.path(out_qc, "emd_by_marker.csv"))
  summary <- data.frame(
    metric = c("emd_reduction", "mad_score", "emd_pairs_used", "emd_pairs_total", "cells_corrected", "cells_evaluated",
               "samples", "batches", "som_nodes_used"),
    value = c(round(emd_reduction, 4), round(mad_score, 4), sum(keep), nrow(emd), sum(som_counts$n), nrow(corrected), nrow(m),
              length(unique(m$batch)), length(unique(labels))))
  write_csv(summary, file.path(out_qc, "evaluation.csv"))
  log_line(sprintf("evaluation: EMD reduction %.3f (%d of %d node/marker/batch-pair EMDs >= 2), MAD score %.4f (%.1fs)",
                   emd_reduction, sum(keep), nrow(emd), mad_score, as.numeric(difftime(Sys.time(), t_eval, units = "secs"))))

  # Density ridges before/after per batch (cyCombine::plot_density), on at most
  # PLOT_CELLS_PER_SAMPLE cells per sample, 9 markers per page.
  plot_rows <- unlist(lapply(split(seq_len(nrow(uncorrected)), uncorrected$sample), function(ix) {
    if (length(ix) <= PLOT_CELLS_PER_SAMPLE) return(ix)
    set.seed(p$SEED)
    sort(sample(ix, PLOT_CELLS_PER_SAMPLE))
  }), use.names = FALSE)
  lims <- stats::quantile(as.matrix(uncorrected[plot_rows, markers]), c(0.001, 0.999), names = FALSE)
  lims <- lims + c(-0.1, 0.1) * diff(lims)
  pages <- split(markers, ceiling(seq_along(markers) / 9))
  for (pg in seq_along(pages)) {
    ok <- tryCatch({
      plots <- suppressMessages(cyCombine::plot_density(uncorrected[plot_rows, ], corrected[plot_rows, ],
                                                        markers = pages[[pg]], xlims = lims, return_format = "raw"))
      legend <- cowplot::get_plot_component(
        plots[[1]] + ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()),
        "guide-box-bottom", return_all = TRUE)
      plots <- lapply(plots, function(g) g + ggplot2::theme(legend.position = "none"))
      grid <- cowplot::plot_grid(cowplot::plot_grid(plotlist = plots, ncol = 3), legend, ncol = 1,
                                 rel_heights = c(12, 1))
      suppressMessages(ggplot2::ggsave(
        file.path(out_qc, sprintf("density_%02d.png", pg)), grid, width = 13,
        height = 0.5 * length(unique(m$batch)) * ceiling(length(pages[[pg]]) / 3) + 1.5,
        dpi = 90, bg = "white", limitsize = FALSE))
      TRUE
    }, error = function(e) { log_line("WARNING: density plot page ", pg, " failed: ", conditionMessage(e)); FALSE })
  }

  ok_n <- sum(vapply(status, function(s) identical(s$status, "ok"), logical(1)))
  failed_n <- length(status) - ok_n
  prov <- list(
    app = "processor-cycombine", version = app_version(app_dir),
    tool = "cyCombine", tool_versions = package_versions(CYCOMBINE_PKGS),
    cycombine_source = "https://github.com/biosurf/cyCombine@fcd42dbbe92bfa6bffa6094bf9b4de536e71c249",
    parameters = p,
    covariate_used = !is.null(covar),
    subsampled = any(m$used < m$events),
    corrected_markers = data.frame(channel = channels, marker = markers),
    evaluation = list(emd_reduction = round(emd_reduction, 4), mad_score = round(mad_score, 4),
                      cells = length(eval_rows), max_cells_per_sample = EVAL_CELLS_PER_SAMPLE),
    inputs = list(fcs = rel_path(fcs, input_dir), samples = rel_path(samples_csv, input_dir),
                  panel = rel_path(panel_csv, input_dir)),
    timings = list(started = format(t0, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), finished = iso_now(),
                   correction_seconds = round(t_cor_s, 1), combat_seconds = round(t_cb_s, 1),
                   write_seconds = round(t_write_s, 1),
                   total_seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)),
    samples = status, n_ok = ok_n, n_failed = failed_n
  )
  write_json(prov, file.path(out_fcs, "provenance.json"))
  write_json(prov, file.path(out_qc, "provenance.json"))
  log_line(sprintf("done: %d sample(s) corrected, %d failed, %.1fs; memory %s", ok_n, failed_n,
                   as.numeric(difftime(Sys.time(), t0, units = "secs")), mem_status()))
  if (ok_n == 0) stop("every sample failed")
  invisible(prov)
}

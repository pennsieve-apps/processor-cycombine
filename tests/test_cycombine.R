# Unit tests on a tiny synthetic fixture (3 batches x 3 samples x 3000 events).

TEST_ENV <- c(XDIM = "3", YDIM = "3", SEED = "7", MAX_CELLS_PER_SAMPLE = "0")

test_fcs_roundtrip <- function() {
  d <- new_tmp("rt")
  m <- cbind(`FSC-A` = c(1, 2.5, 1e6), `PE-A` = c(-1.25, 0, 3.5))
  kw <- c("$P1N" = "FSC-A", "$P2N" = "PE-A", "$P2S" = "CD19", "$P1R" = "262144", "$P2R" = "262144",
          "CUSTOM|KEY" = "a|b", "$DATE" = "01-Jan-2024", "$BEGINDATA" = "999")
  write_fcs(file.path(d, "x.fcs"), m, kw, c(PENNSIEVE_NORM = "test"))
  expect_true(readChar(file.path(d, "x.fcs"), 6, useBytes = TRUE) == "FCS3.1", "FCS 3.1 header")
  f <- read_fcs(file.path(d, "x.fcs"))
  expect_equal_num(f$exprs, m, 1, "values round-trip (float32)")
  expect_true(identical(unname(f$markers), c("", "CD19")), "PnS kept")
  expect_true(fcs_keyword(f$keywords, "CUSTOM|KEY") == "a|b", "delimiter in keyword escaped")
  expect_true(fcs_keyword(f$keywords, "PENNSIEVE_NORM") == "test", "extra keyword written")
  expect_true(fcs_keyword(f$keywords, "$DATATYPE") == "F", "float32")
  expect_error(write_fcs(file.path(d, "y.fcs"), m[, 2:1], kw), "do not match")
}

test_full_run <- function() {
  root <- new_tmp("full")
  inp <- file.path(root, "in")
  out <- file.path(root, "out")
  make_synthetic(inp)
  prov <- with_env(TEST_ENV, run_cycombine(inp, out, "R", cycombine_params()))
  expect_true(prov$n_ok == 9 && prov$n_failed == 0, "all 9 samples ok")
  expect_true(isTRUE(prov$covariate_used), "condition used as covariate")
  expect_true(!isTRUE(prov$subsampled), "MAX_CELLS_PER_SAMPLE=0 keeps all events")
  for (s in c("b1_ctrl", "b2_s1", "b3_s2")) {
    a <- read_fcs(file.path(inp, "fcs", paste0(s, ".fcs")))
    b <- read_fcs(file.path(out, "fcs", paste0(s, ".fcs")))
    expect_true(identical(colnames(a$exprs), colnames(b$exprs)), "same channels")
    expect_true(identical(a$markers, b$markers), "same PnS")
    expect_true(nrow(a$exprs) == nrow(b$exprs), "all events kept")
    for (ch in c("FSC-A", "SSC-A", "AF-A", "Time")) expect_equal_num(a$exprs[, ch], b$exprs[, ch], 0, ch)
    expect_true(fcs_keyword(b$keywords, "GROUPNAME") == "donor A", "keywords kept")
    expect_true(grepl("^cycombine;covariate=condition;som=3x3;norm=scale;subsample=none",
                      fcs_keyword(b$keywords, "PENNSIEVE_NORM")), "PENNSIEVE_NORM keyword")
  }
  ev <- read.csv(file.path(out, "qc", "cycombine", "evaluation.csv"))
  expect_true(ev$value[ev$metric == "emd_reduction"] > 0.3, "EMD reduced")
  # Same sample, different batches: medians should be closer after correction.
  med <- function(dir, s) apply(read_fcs(file.path(dir, "fcs", s))$exprs[, c("BV421-A", "PE-A", "APC-A")], 2, median)
  spread <- function(dir) mean(apply(rbind(med(dir, "b1_ctrl.fcs"), med(dir, "b2_ctrl.fcs"), med(dir, "b3_ctrl.fcs")), 2, sd))
  expect_true(spread(out) < 0.5 * spread(inp), "replicate medians aligned")
  for (f in c("emd.csv", "mad.csv", "emd_by_marker.csv", "som_node_counts.csv", "provenance.json", "density_01.png")) {
    expect_true(file.exists(file.path(out, "qc", "cycombine", f)), paste("qc file", f))
  }
  expect_true(file.exists(file.path(out, "fcs", "provenance.json")), "fcs provenance")
}

test_subsample_and_no_covariate <- function() {
  root <- new_tmp("sub")
  inp <- file.path(root, "in")
  out <- file.path(root, "out")
  make_synthetic(inp)
  env <- TEST_ENV
  env["MAX_CELLS_PER_SAMPLE"] <- "1000"
  env["COVARIATE_COLUMN"] <- "none"
  env["NORM_METHOD"] <- "rank"
  prov <- with_env(env, run_cycombine(inp, out, "R", cycombine_params()))
  expect_true(prov$n_ok == 9 && isTRUE(prov$subsampled) && !isTRUE(prov$covariate_used), "subsampled, no covariate")
  a <- read_fcs(file.path(inp, "fcs", "b2_s1.fcs"))
  b <- read_fcs(file.path(out, "fcs", "b2_s1.fcs"))
  expect_true(nrow(b$exprs) == 1000, "1000 events written")
  expect_true(grepl("subsample=1000/3000", fcs_keyword(b$keywords, "PENNSIEVE_NORM")), "subsample recorded")
  # The written events are the seeded subsample, in file order: scatter/Time identical.
  rows <- seeded_rows(3000, 1000, 7 + 5)  # b2_s1 is the 5th sample in batch/sample_id order
  expect_equal_num(a$exprs[rows, "Time"], b$exprs[, "Time"], 0, "same events")
  # Reproducible with the same SEED.
  out2 <- file.path(root, "out2")
  with_env(env, run_cycombine(inp, out2, "R", cycombine_params()))
  expect_equal_num(read_fcs(file.path(out2, "fcs", "b2_s1.fcs"))$exprs, b$exprs, 0, "reproducible")
}

test_partial_failures_and_unlinked_layout <- function() {
  root <- new_tmp("fail")
  src <- file.path(root, "src")
  make_synthetic(src, extra_sample = "stray")
  inp <- file.path(root, "in")
  dir.create(inp)
  for (f in list.files(file.path(src, "fcs"), full.names = TRUE)) file.symlink(f, file.path(inp, basename(f)))
  for (f in list.files(file.path(src, "metadata"), full.names = TRUE)) file.symlink(f, file.path(inp, basename(f)))
  prov <- with_env(TEST_ENV, run_cycombine(inp, file.path(root, "out"), "R", cycombine_params()))
  expect_true(prov$samples$stray$status == "failed", "sample missing from samples.csv fails")
  expect_true(prov$n_ok == 9, "others ok")
  expect_error(with_env(c(TEST_ENV, COVARIATE_COLUMN = "nope"),
                        run_cycombine(inp, file.path(root, "out2"), "R", cycombine_params())), "lacks column")
}

test_params <- function() {
  expect_error(with_env(c(NORM_METHOD = "clr"), cycombine_params()), "NORM_METHOD")
  expect_error(with_env(c(MAX_CELLS_PER_SAMPLE = "-1"), cycombine_params()), "MAX_CELLS_PER_SAMPLE")
  p <- with_env(c(COVARIATE_COLUMN = "none"), cycombine_params())
  expect_true(identical(p$COVARIATE_COLUMN, ""), "none disables the covariate")
  p <- cycombine_params()
  expect_true(p$XDIM == 8 && p$YDIM == 8 && p$NORM_METHOD == "scale" && p$MAX_CELLS_PER_SAMPLE == 50000, "defaults")
}

#!/usr/bin/env Rscript
# ldsc_rg.R
# Cohort-level pairwise genetic correlation (rg) across all traits' munged
# sumstats, via the real bulik/ldsc software (ldsc.py --rg). Cross-trait, so it
# is not a per-trait fragment like ldsc_h2.R's output — one row per trait pair,
# written once for the whole cohort.
#
# ldsc.py --rg's own comma-list semantics only compute phenotype1-vs-rest, not
# the full N-choose-2 matrix (confirmed via its own log output: "Computing rg
# for phenotype 2/4" etc., every subsequent phenotype compared only against
# phenotype 1). A single whole-cohort call therefore silently produces only
# N-1 of the N*(N-1)/2 possible pairs. Fixed by calling ldsc.py --rg once per
# remaining "first" trait, dropping it from the list each time - N-1 calls
# total for N traits, every pair covered exactly once, no duplication.
#
# Fewer than two traits (rg is undefined) or an ldsc.py failure do not abort the
# run: an empty summary is written and the process exits 0, same as coloc.R's
# note_exit() pattern for single-locus failures.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--sumstats",    type = "character", help = "Comma-separated munged sumstats files (munge_sumstats.py output, .sumstats.gz), in the same order as --ids"),
    make_option("--ids",         type = "character", help = "Comma-separated trait ids, same order as --sumstats"),
    make_option("--ld_dir",      type = "character", help = "Directory of LDSC LD scores (eur_w_ld_chr)"),
    make_option("--out_summary", type = "character", help = "Output TSV: one row per trait pair")
)))

empty_summary <- function() {
    data.table(id1 = character(), id2 = character(), rg = numeric(), rg_se = numeric(),
               z = numeric(), p = numeric(),
               h2_obs1 = numeric(), h2_obs1_se = numeric(),
               h2_int1 = numeric(), h2_int1_se = numeric(),
               gcov_int = numeric(), gcov_int_se = numeric())
}

note_exit <- function(msg) {
    cat("[ldsc_rg] NOTE:", msg, "\n")
    fwrite(empty_summary(), opt$out_summary, sep = "\t")
    quit(status = 0)
}

sumstats <- strsplit(opt$sumstats, ",")[[1]]
ids      <- strsplit(opt$ids, ",")[[1]]

if (length(sumstats) != length(ids)) {
    note_exit(sprintf("mismatched --sumstats (%d) and --ids (%d) counts", length(sumstats), length(ids)))
}
if (length(sumstats) < 2) {
    note_exit(sprintf("only %d trait(s) — genetic correlation needs at least 2", length(sumstats)))
}

ld_ref <- paste0(opt$ld_dir, "/")

# One ldsc.py --rg call per "first" trait; batch i covers ids[i] against every
# id after it in the list (pairs already covered by earlier batches are never
# repeated). Returns NULL (not an error) if this batch's log/table is missing,
# so one bad batch doesn't lose the pairs the other batches already found.
run_batch <- function(batch_sumstats, batch_ids, batch_idx) {
    prefix   <- sprintf("cohort.rg_%d", batch_idx)
    log_file <- paste0(prefix, ".log")

    status <- system2("ldsc.py",
        c("--rg", paste(batch_sumstats, collapse = ","),
          "--ref-ld-chr", ld_ref,
          "--w-ld-chr",   ld_ref,
          "--out", prefix),
        stdout = TRUE, stderr = TRUE)

    if (!file.exists(log_file)) {
        cat(sprintf("[ldsc_rg] NOTE: batch %d (%s) produced no log file; output: %s\n",
                     batch_idx, paste(batch_ids, collapse = ","), paste(status, collapse = " | ")))
        return(NULL)
    }

    log_lines <- readLines(log_file)

    if (length(grep("ERROR", log_lines)) > 0) {
        cat(sprintf("[ldsc_rg] NOTE: batch %d reported an error: %s\n", batch_idx,
                     paste(grep("ERROR", log_lines, value = TRUE), collapse = "; ")))
        return(NULL)
    }

    hdr_i <- grep("^Summary of Genetic Correlation Results", log_lines)
    if (length(hdr_i) == 0) {
        cat(sprintf("[ldsc_rg] NOTE: batch %d - no 'Summary of Genetic Correlation Results' table found\n", batch_idx))
        return(NULL)
    }

    tbl_block <- log_lines[(hdr_i[1] + 1):length(log_lines)]
    # ldsc.py follows the table with a blank line and then "Analysis finished at ..."
    # / "Total time elapsed: ..." — stop at the first blank line so those don't get
    # fed into fread() as bogus table rows.
    first_blank <- which(!nzchar(trimws(tbl_block)))[1]
    if (!is.na(first_blank)) tbl_block <- tbl_block[seq_len(first_blank - 1)]

    tbl <- fread(paste(tbl_block, collapse = "\n"))

    # p1/p2 columns are the munged sumstats file paths ldsc.py was given; map them
    # back to trait ids via this batch's own sumstats/ids ordering (basename match,
    # since ldsc.py may echo relative or absolute paths).
    path_to_id <- setNames(batch_ids, basename(batch_sumstats))
    tbl[, id1 := path_to_id[basename(p1)]]
    tbl[, id2 := path_to_id[basename(p2)]]

    tbl[, .(id1, id2, rg, rg_se = se, z, p,
            h2_obs1 = h2_obs, h2_obs1_se = h2_obs_se,
            h2_int1 = h2_int, h2_int1_se = h2_int_se,
            gcov_int, gcov_int_se)]
}

results <- list()
for (i in seq_len(length(ids) - 1)) {
    batch_ids      <- ids[i:length(ids)]
    batch_sumstats <- sumstats[i:length(sumstats)]
    res <- run_batch(batch_sumstats, batch_ids, i)
    if (!is.null(res)) results[[length(results) + 1]] <- res
}

out <- if (length(results)) rbindlist(results) else empty_summary()

fwrite(out, opt$out_summary, sep = "\t")

# Concatenate every batch's log, in run order, into the single filename the
# Nextflow module already declares (cohort.rg.log) so no .nf change is needed
# for this fix. Built from the known 1:N-1 sequence rather than Sys.glob()+
# sort(), since a lexicographic string sort would misorder batch 10 before
# batch 2 once there are more than 9 traits.
batch_logs <- Filter(file.exists, sprintf("cohort.rg_%d.log", seq_len(length(ids) - 1)))
if (length(batch_logs)) {
    file.create("cohort.rg.log")
    file.append("cohort.rg.log", batch_logs)
}

cat(sprintf("[ldsc_rg] %d trait(s), %d pair(s) written to %s\n",
            length(ids), nrow(out), opt$out_summary))

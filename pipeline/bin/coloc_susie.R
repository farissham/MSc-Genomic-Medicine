#!/usr/bin/env Rscript
# coloc_susie.R
# Per-locus MULTI-SIGNAL colocalization (coloc.susie) between a GWAS trait and
# an eQTL/pQTL reference.
#
# coloc.abf() (bin/coloc.R) assumes at most one causal variant per trait per
# locus. coloc.susie() instead fine-maps each trait into its own independent
# credible sets first (via SuSiE), then colocalizes every GWAS-CS x eQTL-CS
# pair — correctly handling loci with multiple independent causal signals
# (allelic heterogeneity), which coloc.abf can misjudge (under-call PP.H4 when
# a real shared signal is masked by a second independent signal at the locus).
#
# Neither eQTLGen nor GTEx-LV publish individual-level genotypes/LD, so BOTH
# traits' LD is approximated from the same external reference panel already
# staged for GWAS-side SuSiE/MAGMA in this pipeline (g1000_eur_hg38) — this is
# an approximation shared equally by both sides, not literal in-sample LD.
# GTEx v8 was itself imputed against 1000G, so this is a reasonably justified
# proxy for GTEx-LV specifically; less directly justified for eQTLGen (a
# meta-analysis with no single source population), but it is the only external
# LD available for either reference.
#
# Column names are configurable via --hub_*_col / --ref_*_col, mirroring
# bin/coloc.R's interface. --ld/--bim are the SAME per-locus outputs COMPUTE_LD
# already produces for SuSiE fine-mapping elsewhere in this pipeline (panel-only,
# independent of trait) — reused here, not recomputed.
#
# NOTE: coloc.susie()'s result$summary column names beyond nsnps/PP.H0-4.abf
# (idx1/idx2/hit1/hit2 in the coloc package as of this writing) are not
# hardcoded below — they're passed through whatever the installed coloc
# version actually returns, verified only against the required PP/nsnps
# columns. Worth a quick sanity check of the real header on first live run.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(coloc)
    library(susieR)
})

opt <- parse_args(OptionParser(option_list = list(

    # --- input files ---
    make_option("--hub", type = "character", help = "Trait 1 (GWAS) sumstats .tsv(.gz)"),
    make_option("--ref", type = "character", help = "Trait 2 (eQTL/pQTL) sumstats .tsv(.gz)"),

    # --- shared LD (same COMPUTE_LD output SuSiE fine-mapping uses) ---
    make_option("--ld",  type = "character", help = "Locus LD matrix (.ld, square, no header)"),
    make_option("--bim", type = "character", help = "Locus .bim in LD-matrix row/col order (chr id cm pos A1 A2)"),

    # --- locus ---
    make_option("--chr",   type = "integer", help = "Chromosome"),
    make_option("--start", type = "integer", help = "Locus start (bp)"),
    make_option("--end",   type = "integer", help = "Locus end (bp)"),

    # --- trait 1 (hub) column names ---
    make_option("--hub_snp_col",  type = "character", default = "rsid",  help = "SNP ID column in hub"),
    make_option("--hub_chr_col",  type = "character", default = "chr",   help = "Chromosome column in hub"),
    make_option("--hub_pos_col",  type = "character", default = "pos",   help = "Position column in hub"),
    make_option("--hub_ea_col",   type = "character", default = "ea",    help = "Effect allele column in hub"),
    make_option("--hub_oa_col",   type = "character", default = "oa",    help = "Other allele column in hub"),
    make_option("--hub_beta_col", type = "character", default = "beta",  help = "Beta column in hub"),
    make_option("--hub_se_col",   type = "character", default = "se",    help = "SE column in hub"),
    make_option("--hub_n_col",    type = "character", default = "n",     help = "Sample size column in hub"),

    # --- trait 2 (ref) column names ---
    make_option("--ref_snp_col",  type = "character", default = "rsid",  help = "SNP ID column in ref"),
    make_option("--ref_chr_col",  type = "character", default = "chr",   help = "Chromosome column in ref"),
    make_option("--ref_pos_col",  type = "character", default = "pos",   help = "Position column in ref"),
    make_option("--ref_ea_col",   type = "character", default = "ea",    help = "Effect allele column in ref"),
    make_option("--ref_oa_col",   type = "character", default = "oa",    help = "Other allele column in ref"),
    make_option("--ref_beta_col", type = "character", default = "beta",  help = "Beta column in ref"),
    make_option("--ref_se_col",   type = "character", default = "se",    help = "SE column in ref"),
    make_option("--ref_n_col",    type = "character", default = "n",     help = "Sample size column in ref"),
    make_option("--ref_zscore_col", type = "character", default = NULL,
                help = "Z-score column in ref (eQTLGen). When set, approx beta=z/sqrt(N), se=1/sqrt(N)."),

    # --- dataset types / sample sizes ---
    make_option("--type1", type = "character", default = "quant",     help = "Trait 1 type: quant or cc"),
    make_option("--type2", type = "character", default = "quant",     help = "Trait 2 type: quant or cc"),
    make_option("--n1",    type = "integer",   default = NA_integer_, help = "Sample size fallback for trait 1"),
    make_option("--n2",    type = "integer",   default = NA_integer_, help = "Sample size fallback for trait 2"),

    # --- SuSiE fitting (same defaults as bin/susie.R for consistency) ---
    make_option("--l",        type = "integer", default = 10,   help = "Max causal variants per trait (SuSiE L)"),
    make_option("--max_iter", type = "integer", default = 1000, help = "Max SuSiE IBSS iterations"),

    # --- coloc priors ---
    make_option("--p1",       type = "double",  default = 1e-4, help = "Prior: P(SNP assoc with trait 1)"),
    make_option("--p2",       type = "double",  default = 1e-4, help = "Prior: P(SNP assoc with trait 2)"),
    make_option("--p12",      type = "double",  default = 1e-5, help = "Prior: P(SNP assoc with both)"),
    make_option("--min_snps", type = "integer", default = 50,   help = "Min common SNPs (hub ∩ ref ∩ panel) per locus"),

    # --- outputs ---
    make_option("--out_summary", type = "character", help = "Per-CS-pair coloc.susie summary TSV"),
    make_option("--gene_id",     type = "character", default = NULL, help = "Gene identifier written to summary")
)))

comp <- function(x) chartr("ACGTacgt", "TGCAtgca", x)

read_sumstats <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        fread(cmd = sprintf("gzip -dc %s", shQuote(path)))
    } else {
        fread(path)
    }
}

standardise_cols <- function(dt, snp, chr, pos, ea, oa, beta, se, n) {
    mapping <- c(rsid = snp, chr = chr, pos = pos, ea = ea, oa = oa, beta = beta, se = se, n = n)
    for (std in names(mapping)) {
        orig <- mapping[[std]]
        if (orig %in% names(dt) && orig != std) setnames(dt, orig, std)
    }
    dt
}

# Write a single-row placeholder and exit 0 — never let one locus/gene kill the run.
note_exit <- function(msg) {
    cat("[coloc_susie] NOTE:", msg, "\n")
    fwrite(data.table(gene_id = if (!is.null(opt$gene_id)) opt$gene_id else "",
                       chr = opt$chr, start = opt$start, end = opt$end,
                       nsnps = 0L,
                       PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
                       PP.H3.abf = NA_real_, PP.H4.abf = NA_real_, note = msg),
           opt$out_summary, sep = "\t")
    quit(status = 0)
}

## ---- read + standardise column names -----------------------------------------
hub <- read_sumstats(opt$hub)
ref <- read_sumstats(opt$ref)

hub <- standardise_cols(hub, opt$hub_snp_col, opt$hub_chr_col, opt$hub_pos_col,
                         opt$hub_ea_col, opt$hub_oa_col, opt$hub_beta_col, opt$hub_se_col, opt$hub_n_col)
ref <- standardise_cols(ref, opt$ref_snp_col, opt$ref_chr_col, opt$ref_pos_col,
                         opt$ref_ea_col, opt$ref_oa_col, opt$ref_beta_col, opt$ref_se_col, opt$ref_n_col)

if (!is.null(opt$ref_zscore_col)) {
    if (!opt$ref_zscore_col %in% names(ref)) stop("ref missing z-score column: ", opt$ref_zscore_col)
    if (!"n" %in% names(ref) && is.na(opt$n2)) stop("ref z-score mode requires n column (--ref_n_col) or --n2 fallback")
    ref[, zscore_tmp := as.numeric(get(opt$ref_zscore_col))]
    ref[, n_tmp := if ("n" %in% names(ref)) as.numeric(n) else as.numeric(opt$n2)]
    ref[, se   := 1 / sqrt(n_tmp)]
    ref[, beta := zscore_tmp * se]
    ref[, c("zscore_tmp", "n_tmp") := NULL]
    cat("[coloc_susie] ref z-score mode: beta and se approximated from", opt$ref_zscore_col, "\n")
}

for (col in c("rsid", "chr", "pos", "ea", "oa", "beta", "se")) {
    if (!col %in% names(hub)) stop("hub missing column after rename: ", col)
    if (!col %in% names(ref)) stop("ref missing column after rename: ", col)
}

## ---- subset to locus -----------------------------------------------------------
hub <- hub[chr == opt$chr & pos >= opt$start & pos <= opt$end]
ref <- ref[chr == opt$chr & pos >= opt$start & pos <= opt$end]
if (nrow(hub) == 0) note_exit("no hub SNPs in this locus")
if (nrow(ref) == 0) note_exit("no ref SNPs in this locus")

hub[, c("ea", "oa") := .(toupper(ea), toupper(oa))]
ref[, c("ea", "oa") := .(toupper(ea), toupper(oa))]

dedup_rsid <- function(dt) {
    if (!anyDuplicated(dt$rsid)) return(dt)
    dt[, absz_tmp := abs(as.numeric(beta) / as.numeric(se))]
    dt <- dt[order(-absz_tmp)][!duplicated(rsid)]
    dt[, absz_tmp := NULL]
    dt
}
hub <- dedup_rsid(hub)
ref <- dedup_rsid(ref)

## ---- panel .bim (LD-matrix row/col order) --------------------------------------
bim <- fread(opt$bim, header = FALSE, col.names = c("bchr", "bid", "bcm", "bpos", "A1", "A2"))
bim[, c("A1", "A2") := .(toupper(A1), toupper(A2))]
bim[, idx := .I]   # row order == LD matrix row/col order

## Align one trait's beta onto the panel coding (join by rsID, chr:pos fallback;
## flip beta on swapped/reverse-strand alleles; drop palindromes/mismatches).
## Same logic as bin/susie.R's z-alignment, applied independently to hub and
## ref since they need not share identical rsID coverage against the panel.
align_to_panel <- function(dt) {
    m_rsid <- merge(dt, bim[, .(bid, A1, A2, idx)], by.x = "rsid", by.y = "bid")
    if (nrow(m_rsid)) m_rsid[, match_type := "rsid"]
    un <- dt[!rsid %in% bim$bid]
    m_pos <- merge(un, bim[, .(bchr, bpos, A1, A2, idx)], by.x = c("chr", "pos"), by.y = c("bchr", "bpos"))
    if (nrow(m_pos)) m_pos[, match_type := "pos"]
    m <- rbind(m_rsid, m_pos, fill = TRUE)
    if (nrow(m) == 0) return(m)
    m[, action := fcase(
        ea == A1 & oa == A2,             "keep",
        ea == A2 & oa == A1,             "flip",
        comp(ea) == A1 & comp(oa) == A2, "keep",
        comp(ea) == A2 & comp(oa) == A1, "flip",
        default = "drop")]
    m[ea == comp(oa), action := "drop"]   # palindrome: strand unknowable
    m[action == "flip", beta := -as.numeric(beta)]
    m[action != "drop"]
}

hub_a <- align_to_panel(hub)
ref_a <- align_to_panel(ref)
if (nrow(hub_a) == 0) note_exit("no hub SNPs aligned to the LD panel")
if (nrow(ref_a) == 0) note_exit("no ref SNPs aligned to the LD panel")

## ---- common SNP set: hub ∩ ref ∩ panel, one consistent order for R + both z ----
common_idx <- sort(intersect(hub_a$idx, ref_a$idx))
if (length(common_idx) < opt$min_snps)
    note_exit(sprintf("only %d common SNPs (hub ∩ ref ∩ panel) (< min_snps %d)", length(common_idx), opt$min_snps))

setkey(hub_a, idx); setkey(ref_a, idx)
hub_c <- hub_a[.(common_idx)]
ref_c <- ref_a[.(common_idx)]
setorder(hub_c, idx); setorder(ref_c, idx)   # identical row order on both sides + R

## ---- shared LD matrix, subset + conditioned to the common SNP set --------------
R <- as.matrix(fread(opt$ld, header = FALSE))
R <- R[common_idx, common_idx, drop = FALSE]   # common_idx sorted ascending, matches hub_c/ref_c
R <- (R + t(R)) / 2
diag(R) <- 1.0
R <- R * 0.99 + diag(nrow(R)) * 0.01           # mild conditioning, as in bin/susie.R

## ---- sample sizes ---------------------------------------------------------------
n1 <- suppressWarnings(max(as.numeric(hub_c$n), na.rm = TRUE))
if (!is.finite(n1) || n1 <= 0) n1 <- if (!is.na(opt$n1)) opt$n1 else stop("trait 1: no valid N")
n2 <- suppressWarnings(max(as.numeric(ref_c$n), na.rm = TRUE))
if (!is.finite(n2) || n2 <= 0) n2 <- if (!is.na(opt$n2)) opt$n2 else stop("trait 2: no valid N")

## ---- fit SuSiE independently per trait against the shared LD -------------------
fit_susie <- function(dt, n_eff, label) {
    z <- as.numeric(dt$beta) / as.numeric(dt$se)
    tryCatch(
        susie_rss(z = z, R = R, n = n_eff, L = opt$l, max_iter = opt$max_iter, tol = 1e-3),
        error = function(e) note_exit(paste0(label, ": susie_rss failed: ", conditionMessage(e)))
    )
}
fit1 <- fit_susie(hub_c, n1, "hub")
fit2 <- fit_susie(ref_c, n2, "ref")

cat(sprintf("[coloc_susie] chr%d:%d-%d  hub CS=%d  ref CS=%d  common SNPs=%d\n",
            opt$chr, opt$start, opt$end, length(fit1$sets$cs), length(fit2$sets$cs), length(common_idx)))

if (length(fit1$sets$cs) == 0) note_exit("hub: no credible sets (coloc.susie needs >=1 CS per trait)")
if (length(fit2$sets$cs) == 0) note_exit("ref: no credible sets (coloc.susie needs >=1 CS per trait)")

## ---- coloc.susie: pairwise colocalization across every GWAS-CS x eQTL-CS pair --
result <- tryCatch(
    coloc.susie(fit1, fit2, p1 = opt$p1, p2 = opt$p2, p12 = opt$p12),
    error = function(e) note_exit(paste("coloc.susie error:", conditionMessage(e)))
)

if (is.null(result$summary) || nrow(result$summary) == 0)
    note_exit("coloc.susie returned no CS-pair results")

## ---- write outputs: one row per GWAS-CS x eQTL-CS pair -------------------------
smry <- as.data.table(result$summary)
required <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
missing_req <- setdiff(required, names(smry))
if (length(missing_req))
    stop("coloc.susie() result$summary missing expected columns: ", paste(missing_req, collapse = ", "))

smry[, `:=`(gene_id = if (!is.null(opt$gene_id)) opt$gene_id else "",
            chr = opt$chr, start = opt$start, end = opt$end, note = "")]
setcolorder(smry, c("gene_id", "chr", "start", "end",
                     setdiff(names(smry), c("gene_id", "chr", "start", "end", "note")), "note"))
fwrite(smry, opt$out_summary, sep = "\t")

cat(sprintf("[coloc_susie] chr%d:%d-%d  %d CS-pair(s)  max PP.H4.abf=%.4f\n",
            opt$chr, opt$start, opt$end, nrow(smry), max(smry$PP.H4.abf, na.rm = TRUE)))

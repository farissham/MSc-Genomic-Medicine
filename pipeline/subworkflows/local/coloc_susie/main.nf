//
// COLOC_SUSIE_ANALYSIS: multi-signal Bayesian colocalization (coloc.susie) of
// harmonised GWAS traits against a reference eQTL/pQTL dataset, at loci already
// known (from SUSIE_FINEMAP output) to resolve MULTIPLE independent credible
// sets — the specific case coloc.abf's single-causal-variant assumption can
// misjudge. Deliberately scoped to a caller-supplied locus subset rather than
// every --loci entry: re-fitting SuSiE on both the GWAS and eQTL side per gene
// is much more expensive than coloc.abf, and at single-signal loci coloc.abf
// is already the correct, cheaper tool.
//
// Reuses DISCOVER_GENES_IN_LOCUS / FILTER_EQTL_GENE unchanged from the coloc.abf
// subworkflow (../coloc/main.nf) — same gene-discovery logic, only the final
// colocalization step differs. LD is computed once per locus (COMPUTE_LD,
// panel-only, independent of trait) and reused for BOTH the GWAS and eQTL side —
// see bin/coloc_susie.R's header comment for why the same external panel has to
// stand in for both, since neither eQTLGen nor GTEx-LV publish usable LD.
//
// Two ref modes (mirrors coloc.abf's modes 1 and 3 — the two actually used in
// this project; modes 2/4 from coloc.abf were not carried over, add them the
// same way if needed):
//   1. coloc_susie_ref_dir                    : pre-filtered per-gene *.txt files
//   2. coloc_susie_ref + coloc_susie_auto_genes=true : full eQTL, genes auto-discovered per locus
//

include { DISCOVER_GENES_IN_LOCUS } from '../../../modules/local/discover_genes_in_locus'
include { FILTER_EQTL_GENE        } from '../../../modules/local/filter_eqtl_gene'
include { COMPUTE_LD              } from '../../../modules/local/compute_ld'
include { COLOC_SUSIE             } from '../../../modules/local/coloc_susie'

workflow COLOC_SUSIE_ANALYSIS {

    take:
    ch_harmonised   // channel: [ meta, harmonised.tsv.gz ]
    coloc_ref       // path or null: full eQTL file
    coloc_ref_dir   // path or null: directory of pre-filtered per-gene *.txt files
    genotype        // value: staged PLINK panel files (bed/bim/fam) — same panel SuSiE/MAGMA use
    loci_csv        // path: MULTI-SIGNAL loci subset only (columns: chr, start, end)

    main:
    ch_versions = Channel.empty()

    def loci_rows = file(loci_csv, checkIfExists: true)
        .readLines()
        .drop(1)
        .findAll { it.trim() }
        .collect { line ->
            def f = line.split(',')
            [ f[0].trim() as int, f[1].trim() as int, f[2].trim() as int ]
        }

    // LD per locus — same panel-only computation SUSIE_FINEMAP already does;
    // recomputed here (not reused from that channel) since this subworkflow
    // may run against a different, smaller locus subset.
    def ch_loci = Channel.fromList(loci_rows)
    def ch_bed  = genotype.map { fs -> fs.find { it.name.endsWith('.bed') } }.first()
    def ch_bim  = genotype.map { fs -> fs.find { it.name.endsWith('.bim') } }.first()
    def ch_fam  = genotype.map { fs -> fs.find { it.name.endsWith('.fam') } }.first()
    COMPUTE_LD(ch_loci, ch_bed, ch_bim, ch_fam)
    ch_versions = ch_versions.mix(COMPUTE_LD.out.versions)

    // -------------------------------------------------------------------------
    // Build ch_coloc_in per ref mode (gene x locus pairing, pre-LD)
    // -------------------------------------------------------------------------
    def ch_coloc_in

    if (coloc_ref_dir) {
        // Mode 1: pre-filtered files already in a directory -> fan gene x all loci
        def ch_refs = Channel.fromPath("${coloc_ref_dir}/*.txt", checkIfExists: true)
            .map { f -> [ f.baseName, f ] }

        ch_coloc_in = ch_harmonised
            .combine(ch_refs)
            .flatMap { meta, hub, gene_id, ref ->
                loci_rows.collect { locus ->
                    [ meta, gene_id, hub, locus[0], locus[1], locus[2], ref ]
                }
            }

    } else {
        // Mode 3: discover genes per locus, coloc.susie only at the discovered locus.
        def ch_eqtl   = Channel.fromPath(coloc_ref, checkIfExists: true).first()
        def eqtl_path = file(coloc_ref.toString(), checkIfExists: true)

        def ch_discovery_in = ch_harmonised
            .map { meta, hub -> meta }
            .flatMap { meta ->
                loci_rows.collect { locus -> [ meta, locus[0], locus[1], locus[2], eqtl_path ] }
            }

        DISCOVER_GENES_IN_LOCUS(ch_discovery_in)
        ch_versions = ch_versions.mix(DISCOVER_GENES_IN_LOCUS.out.versions.first())

        def ch_gene_locus = DISCOVER_GENES_IN_LOCUS.out.gene_list
            .flatMap { meta, chr, start, end, f ->
                f.readLines().collect { it.trim() }.findAll { it }
                    .collect { gene -> [ meta, gene, chr, start, end ] }
            }

        def ch_unique_genes = ch_gene_locus
            .map  { meta, gene, chr, start, end -> gene }
            .unique()

        FILTER_EQTL_GENE(ch_unique_genes, ch_eqtl)
        ch_versions = ch_versions.mix(FILTER_EQTL_GENE.out.versions.first())

        ch_coloc_in = ch_gene_locus
            .map    { meta, gene, chr, start, end -> [ gene, meta, chr, start, end ] }
            .combine(FILTER_EQTL_GENE.out.filtered, by: [0])
            .map    { gene, meta, chr, start, end, ref -> [ meta, gene, chr, start, end, ref ] }
            .combine(ch_harmonised, by: [0])
            .map    { meta, gene, chr, start, end, ref, hub ->
                [ meta, gene, hub, chr, start, end, ref ]
            }
    }

    // -------------------------------------------------------------------------
    // Attach the shared per-locus LD (keyed join on chr/start/end)
    // -------------------------------------------------------------------------
    def ch_coloc_susie_in = ch_coloc_in
        .map { meta, gene_id, hub, chr, start, end, ref -> [ chr, start, end, meta, gene_id, hub, ref ] }
        .combine(COMPUTE_LD.out.ld, by: [0, 1, 2])
        .map { chr, start, end, meta, gene_id, hub, ref, ld, bim ->
            [ meta, gene_id, hub, chr, start, end, ref, ld, bim ]
        }

    COLOC_SUSIE(ch_coloc_susie_in)
    ch_versions = ch_versions.mix(COLOC_SUSIE.out.versions)

    emit:
    summary  = COLOC_SUSIE.out.summary   // [ meta, gene_id, *.coloc_susie_summary.tsv ]
    versions = ch_versions
}

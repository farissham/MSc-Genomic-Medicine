process COLOC_SUSIE {
    tag "${meta.id}:${gene_id}:chr${chr}:${start}-${end}"
    label 'process_medium'   // two susie_rss() fits per task, heavier than coloc.abf's COLOC

    container 'ghcr.io/farissham/coloc_susie:1.0.0'

    input:
    tuple val(meta), val(gene_id), path(hub), val(chr), val(start), val(end), path(ref), path(ld), path(bim)

    output:
    tuple val(meta), val(gene_id), path("*.coloc_susie_summary.tsv"), emit: summary
    path "versions.yml",                                                emit: versions

    script:
    def args       = task.ext.args ?: ''
    def prefix     = "${meta.id}_${gene_id}_${chr}_${start}_${end}"
    def n2_arg     = params.coloc_susie_n2 != null ? "--n2 ${params.coloc_susie_n2}" : ''
    def zscore_arg = params.coloc_susie_ref_zscore_col != null ? "--ref_zscore_col ${params.coloc_susie_ref_zscore_col}" : ''
    """
    coloc_susie.R \\
        --hub         ${hub} \\
        --ref         ${ref} \\
        --ld          ${ld} \\
        --bim         ${bim} \\
        --chr         ${chr} \\
        --start       ${start} \\
        --end         ${end} \\
        --gene_id     ${gene_id} \\
        --type1       ${params.coloc_susie_trait1_type} \\
        --type2       ${params.coloc_susie_ref_type} \\
        --ref_snp_col ${params.coloc_susie_ref_snp_col} \\
        --ref_chr_col ${params.coloc_susie_ref_chr_col} \\
        --ref_pos_col ${params.coloc_susie_ref_pos_col} \\
        --ref_ea_col  ${params.coloc_susie_ref_ea_col} \\
        --ref_oa_col  ${params.coloc_susie_ref_oa_col} \\
        --ref_n_col   ${params.coloc_susie_ref_n_col} \\
        --l           ${params.coloc_susie_l} \\
        --max_iter    ${params.coloc_susie_max_iter} \\
        --p1          ${params.coloc_susie_p1} \\
        --p2          ${params.coloc_susie_p2} \\
        --p12         ${params.coloc_susie_p12} \\
        --min_snps    ${params.coloc_susie_min_snps} \\
        ${n2_arg} \\
        ${zscore_arg} \\
        --out_summary ${prefix}.coloc_susie_summary.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coloc: \$(Rscript -e 'cat(as.character(packageVersion("coloc")))')
        susieR: \$(Rscript -e 'cat(as.character(packageVersion("susieR")))')
    END_VERSIONS
    """
}

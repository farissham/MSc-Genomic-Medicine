# MSc-Genomic-Medicine

Code availability for the analysis conducted for the MSc research project investigating shared genetic architecture between hypertrophic cardiomyopathy (HCM) and dilated cardiomyopathy (DCM).

## Repository structure

- `pipeline/` — the Nextflow analysis pipeline (SBayesRC, LDSC, Concordance-Analysis, SuSiE-RSS, MAGMA, coloc), also maintained at [farissham/shared_sumstats](https://github.com/farissham/shared_sumstats).
- `figures/` — R/ggplot2 scripts used to generate the dissertation's figures.

## Acknowledgements

The pipeline in `pipeline/` was collaboratively built within the research group. `Harmonise.R` and the SBayesRC module were developed by **Bernard Ooi** (MSc Genomic Medicine) as part of an initially shared pipeline, and are reused here with thanks.

## Data availability

| Dataset | Source |
|---|---|
| In-house case-case GWAS (Primary) | Internal HPC data; not publicly available |
| Reconstructed case-case GWAS (Kramarenko et al.) | Internal HPC data; associated preprint not yet published |
| DCM case-control GWAS (Zheng et al. 2024) | [CVD HuGE Amp Knowledge Portal](https://cvd.hugeamp.org/dinspector.html?dataset=Zheng2024_DCM_EU) |
| HCM case-control GWAS (Tadros et al. 2025) | [GWAS Catalog GCST90435254](https://www.ebi.ac.uk/gwas/studies/GCST90435254) |
| LV structural/functional trait GWAS (10 traits) | GWAS Catalog, accessions [GCST90435258–GCST90435267](https://www.ebi.ac.uk/gwas/studies/GCST90435258) |
| 1000 Genomes Phase 3 EUR reference panel | [CNCR MAGMA resources](https://cncr.nl/research/magma/) |
| Gene location annotation (NCBI38) | [CNCR MAGMA resources](https://cncr.nl/research/magma/) |
| SBayesRC LD reference (ukbEUR_HM3, HapMap3) | [zhilizheng/SBayesRC](https://github.com/zhilizheng/SBayesRC) |
| Gene sets (MSigDB C5, human) | [MSigDB Human Collections](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C5) |
| GTEx v8 tissue expression (54 tissues) | [FUMA download page](https://fuma.ctglab.nl/downloadPage) |
| eQTLGen whole blood cis-eQTLs | [eQTLGen cis-eQTLs](https://www.eqtlgen.org/cis-eqtls.html) |
| GTEx v8 Left Ventricle cis-eQTLs | [GTEx Portal QTL downloads](https://www.gtexportal.org/home/downloads/adult-gtex/qtl) |

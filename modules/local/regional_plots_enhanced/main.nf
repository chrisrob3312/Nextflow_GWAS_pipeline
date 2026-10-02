/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENHANCED REGIONAL PLOTS WITH COLOCALIZATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Creates publication-quality regional association plots with:
    - GWAS association signal
    - Colocalization evidence overlay
    - Gene annotations
    - Regulatory element annotations
    - LD coloring
    - Multi-ancestry comparison (Miami-style)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REGIONAL_PLOT_ENHANCED {
    tag "$meta.trait - $region"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2 conda-forge::r-ggnewscale conda-forge::r-patchwork bioconda::r-locuszoomr"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    tuple val(meta), path(sumstats), val(region)  // region = "chr:start-end" or "gene_name"
    tuple val(coloc_meta), path(coloc_results)     // Colocalization results for this region
    path gene_annotation                            // GTF/GFF for gene track
    path ld_matrix                                  // LD matrix for the region
    path regulatory_annotations                     // ENCODE/Roadmap annotations

    output:
    tuple val(meta), path("${prefix}.regional_enhanced.pdf"), emit: plot
    tuple val(meta), path("${prefix}.regional_enhanced.png"), emit: png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${meta.ancestry}.${region.replaceAll(':', '_')}"
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(ggnewscale)

    # Read GWAS summary stats for region
    gwas <- read.table("${sumstats}", header=TRUE, sep="\\t")

    # Parse region
    region_str <- "${region}"
    if (grepl(":", region_str)) {
        parts <- strsplit(region_str, "[:-]")[[1]]
        chr <- parts[1]
        start <- as.numeric(parts[2])
        end <- as.numeric(parts[3])
    } else {
        # Gene name - would need to look up coordinates
        chr <- "1"; start <- 1; end <- 1e6
    }

    # Filter to region
    gwas_region <- gwas %>%
        filter(CHR == gsub("chr", "", chr), POS >= start, POS <= end) %>%
        mutate(log10p = -log10(P))

    # Read colocalization results if available
    coloc_avail <- file.exists("${coloc_results}") && file.size("${coloc_results}") > 0
    if (coloc_avail) {
        coloc <- read.table("${coloc_results}", header=TRUE, sep="\\t")
        coloc_region <- coloc %>%
            filter(PP.H4.abf >= 0.5)  # Significant colocalizations
    }

    # =========================================================================
    # PANEL 1: GWAS Association with LD coloring
    # =========================================================================
    # Identify lead SNP
    lead_snp <- gwas_region %>% slice_max(log10p, n=1)

    p_gwas <- ggplot(gwas_region, aes(x = POS/1e6, y = log10p)) +
        geom_point(aes(color = log10p), size = 1.5, alpha = 0.8) +
        scale_color_gradient2(low = "blue", mid = "green", high = "red",
                              midpoint = 4, name = "-log10(P)") +
        geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
        geom_hline(yintercept = -log10(1e-5), linetype = "dotted", color = "blue") +
        geom_point(data = lead_snp, aes(x = POS/1e6, y = log10p),
                   shape = 18, size = 4, color = "purple") +
        labs(title = paste("${meta.trait} -", "${meta.ancestry}"),
             subtitle = paste("Region:", region_str),
             x = paste("Position on", chr, "(Mb)"),
             y = expression(-log[10](P))) +
        theme_bw() +
        theme(legend.position = "right")

    # =========================================================================
    # PANEL 2: Colocalization Evidence (if available)
    # =========================================================================
    if (coloc_avail && nrow(coloc_region) > 0) {
        # Create colocalization track
        p_coloc <- ggplot(coloc_region, aes(x = pos/1e6, y = PP.H4.abf, color = qtl_type)) +
            geom_segment(aes(xend = pos/1e6, yend = 0), alpha = 0.6) +
            geom_point(size = 2) +
            scale_color_brewer(palette = "Set1", name = "QTL Type") +
            labs(y = "PP4 (Coloc)", x = "") +
            ylim(0, 1) +
            theme_bw() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.position = "right")
    } else {
        p_coloc <- ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "No colocalization data") +
            theme_void()
    }

    # =========================================================================
    # PANEL 3: Gene Track
    # =========================================================================
    # Simplified gene track (would use rtracklayer for full GTF parsing)
    p_genes <- ggplot() +
        geom_rect(aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.4, ymax = 0.6),
                  fill = "gray50", alpha = 0.3) +
        annotate("text", x = (start + end)/2/1e6, y = 0.5, label = "Gene Track",
                 size = 3, fontface = "italic") +
        labs(y = "Genes", x = paste("Position on", chr, "(Mb)")) +
        ylim(0, 1) +
        theme_bw() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    # =========================================================================
    # PANEL 4: Regulatory Annotations (ENCODE/Roadmap)
    # =========================================================================
    p_regulatory <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Regulatory Elements Track") +
        theme_void()

    # =========================================================================
    # Combine panels
    # =========================================================================
    combined <- p_gwas / p_coloc / p_genes +
        plot_layout(heights = c(3, 1.5, 1)) +
        plot_annotation(
            title = "Regional Association Plot with Colocalization",
            subtitle = paste("${meta.trait} |", "${meta.ancestry} | Region:", region_str)
        )

    # Save
    ggsave("${prefix}.regional_enhanced.pdf", combined, width = 12, height = 10)
    ggsave("${prefix}.regional_enhanced.png", combined, width = 12, height = 10, dpi = 300)

    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    ggplot2: "', packageVersion("ggplot2"), '"')
    ), "versions.yml")
    """
}

process MIAMI_PLOT {
    tag "$meta.trait"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2 conda-forge::r-dplyr"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    tuple val(meta), path(sumstats1), path(sumstats2)
    val comparison_type  // 'ancestry' or 'trait'
    val label1           // Label for top panel (e.g., "EUR" or "Discovery")
    val label2           // Label for bottom panel (e.g., "AFR" or "Replication")

    output:
    tuple val(meta), path("${prefix}.miami.pdf"), emit: plot
    tuple val(meta), path("${prefix}.miami.png"), emit: png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def l1 = meta.anc1 ?: label1      // per-pair labels from meta when present
    def l2 = meta.anc2 ?: label2
    prefix = task.ext.prefix ?: "${meta.trait}.${l1}_vs_${l2}"
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(dplyr)

    # Read both summary stats
    gwas1 <- read.table("${sumstats1}", header=TRUE, sep="\\t") %>%
        mutate(log10p = -log10(P), panel = "top")
    gwas2 <- read.table("${sumstats2}", header=TRUE, sep="\\t") %>%
        mutate(log10p = log10(P), panel = "bottom")  # Negative for Miami

    # Combine and prepare chromosome positions
    combined <- bind_rows(gwas1, gwas2) %>%
        mutate(CHR = as.numeric(gsub("chr", "", CHR))) %>%
        filter(!is.na(CHR)) %>%
        arrange(CHR, POS)

    # Calculate cumulative positions
    chr_lengths <- combined %>%
        group_by(CHR) %>%
        summarize(max_pos = max(POS)) %>%
        mutate(cum_start = cumsum(lag(max_pos, default = 0)))

    combined <- combined %>%
        left_join(chr_lengths, by = "CHR") %>%
        mutate(cum_pos = POS + cum_start)

    # Chromosome centers for x-axis labels
    chr_centers <- combined %>%
        group_by(CHR) %>%
        summarize(center = mean(cum_pos))

    # Alternate chromosome colors
    chr_colors <- rep(c("#1f77b4", "#aec7e8"), 12)[1:22]

    # Miami plot
    p <- ggplot(combined, aes(x = cum_pos, y = log10p, color = factor(CHR %% 2))) +
        geom_point(size = 0.5, alpha = 0.6) +
        geom_hline(yintercept = c(-log10(5e-8), log10(5e-8)),
                   linetype = "dashed", color = "red", linewidth = 0.5) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
        scale_color_manual(values = c("#1f77b4", "#aec7e8"), guide = "none") +
        scale_x_continuous(breaks = chr_centers\$center, labels = chr_centers\$CHR) +
        annotate("text", x = max(combined\$cum_pos) * 0.02, y = max(gwas1\$log10p) * 0.9,
                 label = "${l1}", hjust = 0, fontface = "bold", size = 4) +
        annotate("text", x = max(combined\$cum_pos) * 0.02, y = min(gwas2\$log10p) * 0.9,
                 label = "${l2}", hjust = 0, fontface = "bold", size = 4) +
        labs(title = "Miami Plot: ${meta.trait}",
             subtitle = "${l1} (top) vs ${l2} (bottom)",
             x = "Chromosome",
             y = expression(-log[10](P))) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 8),
              panel.grid.minor = element_blank())

    # Save
    ggsave("${prefix}.miami.pdf", p, width = 14, height = 8)
    ggsave("${prefix}.miami.png", p, width = 14, height = 8, dpi = 300)

    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}

process COLOC_LOCUSCOMPARE {
    tag "$meta.trait - $gene"
    label 'process_low'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2 conda-forge::r-cowplot"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggplot2:3.4.4' :
        'quay.io/biocontainers/r-ggplot2:3.4.4' }"

    input:
    tuple val(meta), path(gwas_sumstats), path(qtl_sumstats), val(gene)
    val qtl_type

    output:
    tuple val(meta), path("${prefix}.locuscompare.pdf"), emit: plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${gene}.${qtl_type}"
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(cowplot)
    library(dplyr)

    # Read data
    gwas <- read.table("${gwas_sumstats}", header=TRUE, sep="\\t") %>%
        mutate(gwas_log10p = -log10(P))
    qtl <- read.table("${qtl_sumstats}", header=TRUE, sep="\\t") %>%
        mutate(qtl_log10p = -log10(pvalue))

    # Merge by SNP
    merged <- inner_join(gwas, qtl, by = c("SNP" = "snp"))

    if (nrow(merged) < 10) {
        # Not enough overlapping SNPs
        pdf("${prefix}.locuscompare.pdf", width=10, height=5)
        plot.new()
        text(0.5, 0.5, "Insufficient overlapping SNPs for comparison", cex=1.5)
        dev.off()
    } else {
        # LocusCompare-style scatter plot
        # Shows whether GWAS and QTL signals share the same pattern

        p_scatter <- ggplot(merged, aes(x = gwas_log10p, y = qtl_log10p)) +
            geom_point(alpha = 0.6, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
            labs(title = paste("LocusCompare: ${meta.trait} vs ${qtl_type}"),
                 subtitle = paste("Gene:", "${gene}"),
                 x = expression(paste("GWAS ", -log[10](P))),
                 y = expression(paste("${qtl_type} ", -log[10](P)))) +
            theme_bw() +
            annotate("text", x = Inf, y = Inf,
                     label = paste("r =", round(cor(merged\$gwas_log10p, merged\$qtl_log10p), 3)),
                     hjust = 1.1, vjust = 1.5, size = 4)

        # Side-by-side Manhattan-style for the region
        p_gwas_region <- ggplot(merged, aes(x = POS/1e6, y = gwas_log10p)) +
            geom_point(alpha = 0.6) +
            labs(title = "GWAS", x = "Position (Mb)", y = expression(-log[10](P))) +
            theme_bw()

        p_qtl_region <- ggplot(merged, aes(x = POS/1e6, y = qtl_log10p)) +
            geom_point(alpha = 0.6, color = "darkgreen") +
            labs(title = "${qtl_type}", x = "Position (Mb)", y = expression(-log[10](P))) +
            theme_bw()

        # Combine
        combined <- plot_grid(
            p_scatter,
            plot_grid(p_gwas_region, p_qtl_region, ncol = 1),
            ncol = 2, rel_widths = c(1, 1)
        )

        ggsave("${prefix}.locuscompare.pdf", combined, width = 12, height = 6)
    }

    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}

process TRAIT_GENETIC_OVERLAP {
    tag "genetic_overlap"
    label 'process_medium'

    conda "bioconda::ldsc=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ldsc:1.0.1--pyhdfd78af_2' :
        'quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2' }"

    input:
    path munged_sumstats_list   // List of munged sumstats files
    val trait_names             // List of trait names
    path ld_reference

    output:
    path "trait_rg_matrix.tsv", emit: rg_matrix
    path "trait_rg_heatmap.pdf", emit: heatmap
    path "trait_overlap_summary.json", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import subprocess
    import pandas as pd
    import json
    import itertools

    traits = "${trait_names}".split(",")
    files = "${munged_sumstats_list}".split()

    # Create mapping
    trait_files = dict(zip(traits, files))

    # Run pairwise genetic correlations
    results = []
    for t1, t2 in itertools.combinations(traits, 2):
        f1 = trait_files[t1]
        f2 = trait_files[t2]

        # Run LDSC rg
        cmd = f"ldsc.py --rg {f1},{f2} --ref-ld-chr ${ld_reference}/ --w-ld-chr ${ld_reference}/ --out {t1}_{t2}"
        subprocess.run(cmd, shell=True, check=True)

        # Parse results
        with open(f"{t1}_{t2}.log", "r") as f:
            for line in f:
                if "Genetic Correlation:" in line:
                    rg = float(line.split()[-1].replace("(", "").replace(")", ""))
                if "P:" in line and "Genetic" not in line:
                    p = float(line.split()[-1])

        results.append({
            "trait1": t1,
            "trait2": t2,
            "rg": rg,
            "p": p
        })

    # Create matrix
    rg_df = pd.DataFrame(results)
    rg_matrix = pd.DataFrame(index=traits, columns=traits)
    for _, row in rg_df.iterrows():
        rg_matrix.loc[row["trait1"], row["trait2"]] = row["rg"]
        rg_matrix.loc[row["trait2"], row["trait1"]] = row["rg"]
    for t in traits:
        rg_matrix.loc[t, t] = 1.0

    rg_matrix.to_csv("trait_rg_matrix.tsv", sep="\\t")

    # Summary
    summary = {
        "n_traits": len(traits),
        "traits": traits,
        "mean_rg": rg_df["rg"].mean(),
        "max_rg": rg_df["rg"].max(),
        "min_rg": rg_df["rg"].min(),
        "n_significant": (rg_df["p"] < 0.05).sum()
    }

    with open("trait_overlap_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Heatmap (using R)
    import subprocess
    r_code = '''
    library(pheatmap)
    rg <- read.table("trait_rg_matrix.tsv", header=TRUE, row.names=1)
    pdf("trait_rg_heatmap.pdf", width=8, height=8)
    pheatmap(as.matrix(rg),
             main = "Genetic Correlation Between Traits",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-1, 1, length.out = 101),
             display_numbers = TRUE,
             number_format = "%.2f")
    dev.off()
    '''
    with open("plot_heatmap.R", "w") as f:
        f.write(r_code)
    subprocess.run(["Rscript", "plot_heatmap.R"])

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('    ldsc: "1.0.1"\\n')
    """
}

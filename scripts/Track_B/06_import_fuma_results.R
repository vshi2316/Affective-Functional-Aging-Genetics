## ============================================================================
## Build FUMA JOB06 and JOB07 annotation and visualization package
## JOB06: affective factor GWAS
## JOB07: functional factor GWAS
## Optional comparator: JOB01 common-factor GWAS
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

resolve_project_base <- function() {
  candidates <- unique(c(
    Sys.getenv("DCV_BASE_DIR", unset = ""),
    Sys.getenv("DCV_PROJECT_DIR", unset = getwd())
  ))
  candidates <- candidates[nzchar(candidates)]
  for (x in candidates) {
    if (dir.exists(file.path(x, "analysis_ready_core"))) {
      return(normalizePath(x, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Could not resolve project base. Set DCV_BASE_DIR to the project root.")
}

base_dir <- resolve_project_base()
job_dirs <- c(
  Affective = Sys.getenv("FUMA_JOB06_DIR", unset = file.path(base_dir, "FUMA", "JOB06")),
  Functional = Sys.getenv("FUMA_JOB07_DIR", unset = file.path(base_dir, "FUMA", "JOB07"))
)
common_dir <- Sys.getenv(
  "FUMA_JOB01_DIR",
  unset = file.path(base_dir, "FUMA", "FUMA_JOB01_affective_functional_factor_GWAS")
)

job_dirs <- vapply(job_dirs, normalizePath, character(1), winslash = "/", mustWork = TRUE)
common_available <- dir.exists(common_dir)
if (common_available) {
  common_dir <- normalizePath(common_dir, winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(base_dir, "FUMA", "FUMA_JOB06_07_twofactor_visualization")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

domain_colors <- c(
  Common = "#5B6573",
  Affective = "#C44E52",
  Functional = "#2F6C8E"
)

theme_nc <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(color = "grey35", size = base_size - 1),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.45),
      strip.background = element_rect(fill = "grey94", color = "grey45", linewidth = 0.4),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(6, 8, 6, 8)
    )
}

safe_fread <- function(path, required = FALSE) {
  if (!file.exists(path)) {
    if (required) stop("Required FUMA file is missing: ", path)
    warning("Optional FUMA file is missing: ", path)
    return(data.table())
  }
  data.table::fread(path, data.table = TRUE, fill = TRUE)
}

read_magma <- function(path, required = FALSE) {
  if (!file.exists(path)) {
    if (required) stop("Required MAGMA file is missing: ", path)
    warning("Optional MAGMA file is missing: ", path)
    return(data.table())
  }
  data.table::fread(path, skip = "VARIABLE", data.table = TRUE, fill = TRUE)
}

p_to_log10 <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  -log10(pmax(p, .Machine$double.xmin))
}

wrap_label <- function(x, width = 38) {
  x <- gsub("_", " ", as.character(x))
  stringr::str_wrap(gsub("\\s+", " ", x), width = width)
}

save_plot <- function(plot, filename, width, height, dpi = 320) {
  ggsave(file.path(fig_dir, paste0(filename, ".png")), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot,
         width = width, height = height, device = pdf_device, bg = "white")
}

read_params <- function(path) {
  x <- readLines(path, warn = FALSE)
  get_value <- function(key) {
    hit <- grep(paste0("^", key, "\\s*="), x, value = TRUE)
    if (!length(hit)) return(NA_character_)
    trimws(sub("^[^=]+\\s*=\\s*", "", hit[1]))
  }
  data.table(
    input_file = get_value("gwasfile"),
    genome_build_grch38 = get_value("snp2genegrch38"),
    lead_p = get_value("leadP"),
    candidate_p = get_value("gwasP"),
    independent = get_value("r2"),
    lead = get_value("r2_2"),
    reference_panel = get_value("refpanel"),
    ancestry = get_value("pop"),
    merge_distance_kb = get_value("mergeDist"),
    mhc_exclusion = get_value("exMHC"),
    mhc_option = get_value("MHCopt"),
    positional_mapping = get_value("posMap"),
    eqtl_mapping = get_value("eqtlMap"),
    chromatin_interaction_mapping = get_value("ciMap"),
    magma = get_value("magma")
  )
}

read_job <- function(job_dir, domain) {
  files <- list(
    loci = safe_fread(file.path(job_dir, "GenomicRiskLoci.txt"), required = TRUE),
    lead = safe_fread(file.path(job_dir, "leadSNPs.txt"), required = TRUE),
    ind_sig = safe_fread(file.path(job_dir, "IndSigSNPs.txt"), required = TRUE),
    snps = safe_fread(file.path(job_dir, "snps.txt"), required = TRUE),
    genes = safe_fread(file.path(job_dir, "genes.txt"), required = TRUE),
    annov_stats = safe_fread(file.path(job_dir, "annov.stats.txt"), required = TRUE),
    gwas_catalog = safe_fread(file.path(job_dir, "gwascatalog.txt")),
    magma_genes = safe_fread(file.path(job_dir, "magma.genes.out"), required = TRUE),
    magma_sets = read_magma(file.path(job_dir, "magma.gsa.out"), required = TRUE),
    gtex_general = read_magma(file.path(job_dir, "magma_exp_gtex_ts_general_avg_log2TPM.gsa.out")),
    gtex_specific = read_magma(file.path(job_dir, "magma_exp_gtex_ts_avg_log2TPM.gsa.out")),
    params = read_params(file.path(job_dir, "params.config"))
  )
  files$domain <- domain
  files$job_dir <- job_dir
  files
}

jobs <- Map(read_job, unname(job_dirs), names(job_dirs))
names(jobs) <- names(job_dirs)
if (common_available) {
  jobs$Common <- read_job(common_dir, "Common")
}

for (nm in names(jobs)) {
  j <- jobs[[nm]]
  if (!all(c("chr", "pos", "p", "start", "end") %in% names(j$loci))) {
    stop("Unexpected GenomicRiskLoci.txt columns for ", nm)
  }
  if (!all(c("GENE", "P") %in% names(j$magma_genes))) {
    stop("Unexpected magma.genes.out columns for ", nm)
  }
}

## Parameter audit ------------------------------------------------------------

parameter_audit <- rbindlist(lapply(names(jobs), function(nm) {
  cbind(data.table(domain = nm, directory = jobs[[nm]]$job_dir), jobs[[nm]]$params)
}), fill = TRUE)
fwrite(parameter_audit, file.path(tab_dir, "fuma_parameter_audit.csv"))

required_equal <- c("genome_build_grch38", "lead_p", "candidate_p", "independent",
                    "lead", "reference_panel", "ancestry", "merge_distance_kb",
                    "mhc_exclusion", "mhc_option", "positional_mapping", "eqtl_mapping",
                    "chromatin_interaction_mapping", "magma")
parameter_consistency <- rbindlist(lapply(required_equal, function(v) {
  values <- unique(parameter_audit[[v]])
  data.table(parameter = v, consistent = length(values) == 1L,
             observed_values = paste(values, collapse = " | "))
}))
fwrite(parameter_consistency, file.path(tab_dir, "fuma_parameter_consistency.csv"))

## Summary tables -------------------------------------------------------------

make_summary <- function(j) {
  magma_p <- suppressWarnings(as.numeric(j$magma_genes$P))
  gene_bonf <- 0.05 / sum(is.finite(magma_p))
  loci_chr <- suppressWarnings(as.integer(j$loci$chr))
  loci_pos <- suppressWarnings(as.numeric(j$loci$pos))
  in_mhc <- loci_chr == 6L & loci_pos >= 25000000 & loci_pos <= 34000000
  data.table(
    domain = j$domain,
    genomic_risk_loci = nrow(j$loci),
    non_mhc_risk_loci = sum(!in_mhc, na.rm = TRUE),
    mhc_risk_loci = sum(in_mhc, na.rm = TRUE),
    lead_snps = nrow(j$lead),
    independent_significant_snps = nrow(j$ind_sig),
    candidate_snps = nrow(j$snps),
    positional_mapped_genes = nrow(j$genes),
    magma_tested_genes = sum(is.finite(magma_p)),
    magma_bonferroni_threshold = gene_bonf,
    magma_bonferroni_genes = sum(magma_p < gene_bonf, na.rm = TRUE),
    magma_bh_fdr_genes = sum(p.adjust(magma_p, "BH") < 0.05, na.rm = TRUE),
    gwas_catalog_records = nrow(j$gwas_catalog)
  )
}

summary_dt <- rbindlist(lapply(jobs, make_summary), fill = TRUE)
summary_dt[, domain := factor(domain, levels = c("Common", "Affective", "Functional"))]
setorder(summary_dt, domain)
fwrite(summary_dt, file.path(tab_dir, "fuma_common_affective_functional_summary.csv"))

export_job_tables <- function(j) {
  prefix <- paste0("fuma_", tolower(j$domain))

  loci <- copy(j$loci)
  loci[, `:=`(
    log10p = p_to_log10(p),
    locus_width_kb = (as.numeric(end) - as.numeric(start)) / 1000,
    mhc_locus = as.integer(chr) == 6L & as.numeric(pos) >= 25000000 & as.numeric(pos) <= 34000000
  )]
  setorder(loci, p)
  fwrite(loci, file.path(tab_dir, paste0(prefix, "_risk_loci.csv")))

  mapped <- copy(j$genes)
  mapped[, minGwasP := as.numeric(minGwasP)]
  mapped[, log10_min_gwas_p := p_to_log10(minGwasP)]
  setorder(mapped, minGwasP)
  fwrite(mapped, file.path(tab_dir, paste0(prefix, "_positionally_mapped_genes.csv")))

  magma <- copy(j$magma_genes)
  magma[, P := as.numeric(P)]
  magma[, `:=`(
    p_bonferroni = pmin(P * sum(is.finite(P)), 1),
    p_fdr = p.adjust(P, "BH"),
    log10p = p_to_log10(P)
  )]
  setorder(magma, P)
  fwrite(magma, file.path(tab_dir, paste0(prefix, "_magma_gene_results.csv")))

  sets <- copy(j$magma_sets)
  sets[, P := as.numeric(P)]
  sets[, `:=`(p_fdr = p.adjust(P, "BH"), log10p = p_to_log10(P))]
  setorder(sets, P)
  fwrite(sets, file.path(tab_dir, paste0(prefix, "_magma_gene_set_results.csv")))

  ann <- copy(j$annov_stats)
  ann[, `:=`(
    enrichment = as.numeric(enrichment),
    fisher.P = as.numeric(fisher.P),
    fisher_fdr = p.adjust(as.numeric(fisher.P), "BH")
  )]
  setorder(ann, fisher.P)
  fwrite(ann, file.path(tab_dir, paste0(prefix, "_annovar_enrichment.csv")))
}
invisible(lapply(jobs, export_job_tables))

## Figure 1: FUMA output counts ----------------------------------------------

count_metrics <- c("genomic_risk_loci", "lead_snps", "independent_significant_snps",
                   "candidate_snps", "positional_mapped_genes", "magma_bonferroni_genes")
count_labels <- c(
  genomic_risk_loci = "Genomic risk loci",
  lead_snps = "Lead SNPs",
  independent_significant_snps = "Independent significant SNPs",
  candidate_snps = "Candidate SNPs",
  positional_mapped_genes = "Positionally mapped genes",
  magma_bonferroni_genes = "MAGMA Bonferroni genes"
)
counts_long <- melt(summary_dt, id.vars = "domain", measure.vars = count_metrics,
                    variable.name = "metric", value.name = "value")
counts_long[, metric_label := factor(count_labels[as.character(metric)],
                                     levels = rev(unname(count_labels)))]

p1 <- ggplot(counts_long, aes(metric_label, value, fill = domain)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = comma(value)), hjust = -0.08, size = 2.8) +
  coord_flip() +
  facet_wrap(~ domain, scales = "free_x") +
  scale_fill_manual(values = domain_colors) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.22))) +
  labs(title = "A  FUMA annotation yield", x = NULL, y = "Count") +
  theme_nc(9)
save_plot(p1, "figure_1A_fuma_annotation_yield", 12.0, 5.0)

## Figure 2: genomic risk loci by chromosome --------------------------------

loci_chr <- rbindlist(lapply(jobs, function(j) {
  x <- as.data.table(j$loci)[, .N, by = .(chr = as.character(chr))]
  full <- data.table(chr = as.character(1:22))
  x <- merge(full, x, by = "chr", all.x = TRUE)
  x[is.na(N), N := 0L]
  x[, domain := j$domain]
  x
}))
loci_chr[, chr := factor(chr, levels = as.character(1:22))]
loci_chr[, domain := factor(domain, levels = c("Common", "Affective", "Functional"))]
fwrite(loci_chr, file.path(tab_dir, "fuma_risk_loci_by_chromosome.csv"))

p2 <- ggplot(loci_chr, aes(chr, N, fill = domain)) +
  geom_col(width = 0.76, show.legend = FALSE) +
  facet_wrap(~ domain, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = domain_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "B  Genomic risk loci by chromosome",
       subtitle = "FUMA loci use the task-specific LD and merging parameters recorded in the parameter audit",
       x = "Chromosome", y = "Number of loci") +
  theme_nc(9)
save_plot(p2, "figure_1B_risk_loci_by_chromosome", 10.5, 7.0)

## Figure 3: top MAGMA genes -------------------------------------------------

top_magma <- rbindlist(lapply(jobs, function(j) {
  x <- copy(j$magma_genes)
  x[, P := as.numeric(P)]
  x[, `:=`(
    p_fdr = p.adjust(P, "BH"),
    bonferroni = P < 0.05 / sum(is.finite(P)),
    log10p = p_to_log10(P),
    gene_label = fifelse(!is.na(SYMBOL) & nzchar(SYMBOL), SYMBOL, GENE),
    domain = j$domain
  )]
  setorder(x, P)
  head(x, 20)
}), fill = TRUE)
top_magma[, plot_label := paste(domain, gene_label, sep = "__")]
top_magma[, plot_label := factor(plot_label, levels = rev(unique(plot_label)))]
fwrite(top_magma, file.path(tab_dir, "fuma_top20_magma_genes_by_domain.csv"))

p3 <- ggplot(top_magma, aes(plot_label, log10p)) +
  geom_segment(aes(xend = plot_label, y = 0, yend = log10p), color = "grey72", linewidth = 0.4) +
  geom_point(aes(size = NSNPS, color = bonferroni), alpha = 0.92) +
  coord_flip() +
  facet_wrap(~ domain, scales = "free_y", ncol = length(jobs)) +
  scale_x_discrete(labels = function(x) sub("^[^_]+__", "", x)) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                     name = "Bonferroni significant") +
  scale_size_continuous(name = "SNPs in gene", range = c(1.8, 5.2)) +
  labs(title = "C  Top MAGMA gene-based associations", x = NULL,
       y = expression(-log[10](italic(P)))) +
  theme_nc(9)
save_plot(p3, "figure_2A_top_magma_genes", 13.5, 7.2)

## Figure 4: MAGMA gene-set results ------------------------------------------

top_sets <- rbindlist(lapply(jobs, function(j) {
  x <- copy(j$magma_sets)
  x[, P := as.numeric(P)]
  x[, `:=`(
    p_fdr = p.adjust(P, "BH"),
    log10p = p_to_log10(P),
    set_label = wrap_label(if ("FULL_NAME" %in% names(x)) FULL_NAME else VARIABLE, 38),
    domain = j$domain
  )]
  setorder(x, P)
  head(x, 12)
}), fill = TRUE)
top_sets[, plot_label := paste(domain, set_label, sep = "__")]
top_sets[, plot_label := factor(plot_label, levels = rev(unique(plot_label)))]
fwrite(top_sets, file.path(tab_dir, "fuma_top12_magma_gene_sets_by_domain.csv"))

p4 <- ggplot(top_sets, aes(plot_label, log10p)) +
  geom_segment(aes(xend = plot_label, y = 0, yend = log10p), color = "grey75", linewidth = 0.4) +
  geom_point(aes(size = NGENES, color = p_fdr), alpha = 0.94) +
  coord_flip() +
  facet_wrap(~ domain, scales = "free_y", ncol = length(jobs)) +
  scale_x_discrete(labels = function(x) sub("^[^_]+__", "", x)) +
  scale_color_gradient(low = "#B2182B", high = "#F1D5C9", trans = "reverse", name = "BH FDR") +
  scale_size_continuous(name = "Genes in set", range = c(1.8, 5.2)) +
  labs(title = "D  Top MAGMA gene-set associations", x = NULL,
       y = expression(-log[10](italic(P)))) +
  theme_nc(8.5)
save_plot(p4, "figure_2B_top_magma_gene_sets", 15.0, 7.4)

## Figure 5: GTEx expression enrichment --------------------------------------

prep_gtex <- function(j, source, x) {
  if (!nrow(x)) return(data.table())
  x <- copy(x)
  x[, P := as.numeric(P)]
  x[, `:=`(
    p_fdr = p.adjust(P, "BH"),
    log10p = p_to_log10(P),
    tissue = if ("FULL_NAME" %in% names(x)) FULL_NAME else VARIABLE,
    domain = j$domain,
    source = source
  )]
  setorder(x, P)
  head(x, 12)
}

gtex <- rbindlist(lapply(jobs, function(j) {
  rbindlist(list(
    prep_gtex(j, "GTEx v8 general tissue", j$gtex_general),
    prep_gtex(j, "GTEx v8 tissue specific", j$gtex_specific)
  ), fill = TRUE)
}), fill = TRUE)
gtex[, tissue_label := wrap_label(tissue, 30)]
gtex[, plot_label := paste(domain, source, tissue_label, sep = "__")]
gtex[, plot_label := factor(plot_label, levels = rev(unique(plot_label)))]
fwrite(gtex, file.path(tab_dir, "fuma_top_gtex_expression_results.csv"))

if (nrow(gtex)) {
  p5 <- ggplot(gtex, aes(plot_label, log10p)) +
    geom_segment(aes(xend = plot_label, y = 0, yend = log10p), color = "grey75", linewidth = 0.4) +
    geom_point(aes(size = abs(as.numeric(BETA_STD)), color = p_fdr), alpha = 0.94) +
    coord_flip() +
    facet_grid(source ~ domain, scales = "free_y", space = "free_y") +
    scale_x_discrete(labels = function(x) sub("^[^_]+__[^_]+__", "", x)) +
    scale_color_gradient(low = "#B2182B", high = "#F1D5C9", trans = "reverse", name = "BH FDR") +
    scale_size_continuous(name = "|Standardized beta|", range = c(1.6, 5.0)) +
    labs(title = "E  MAGMA gene-property analysis across GTEx v8 tissues",
         x = NULL, y = expression(-log[10](italic(P)))) +
    theme_nc(8.5)
  save_plot(p5, "figure_3A_gtex_expression_enrichment", 15.0, 12.0)
}

## Figure 6: ANNOVAR enrichment ----------------------------------------------

annov <- rbindlist(lapply(jobs, function(j) {
  x <- copy(j$annov_stats)
  x[, `:=`(
    enrichment = as.numeric(enrichment),
    fisher.P = as.numeric(fisher.P),
    fisher_fdr = p.adjust(as.numeric(fisher.P), "BH"),
    domain = j$domain
  )]
  x
}), fill = TRUE)
annov[, annot := factor(annot, levels = rev(unique(annot)))]
fwrite(annov, file.path(tab_dir, "fuma_annovar_enrichment_all_domains.csv"))

p6 <- ggplot(annov, aes(annot, enrichment, color = fisher_fdr < 0.05)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey55", linewidth = 0.45) +
  geom_point(size = 2.8) +
  coord_flip() +
  facet_wrap(~ domain, nrow = 1) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#5B6573"),
                     name = "BH FDR < 0.05") +
  labs(title = "F  Functional consequence enrichment among candidate SNPs",
       subtitle = "Enrichment relative to the selected 1000 Genomes reference panel",
       x = NULL, y = "Enrichment ratio") +
  theme_nc(9)
save_plot(p6, "figure_3B_annovar_functional_enrichment", 12.5, 5.5)

## Cross-domain gene membership ----------------------------------------------

membership_table <- function(gene_lists, evidence_type) {
  universe <- sort(unique(unlist(gene_lists)))
  if (!length(universe)) return(data.table())
  out <- data.table(gene = universe)
  for (nm in names(gene_lists)) out[, (nm) := gene %in% gene_lists[[nm]]]
  cols <- names(gene_lists)
  out[, membership := apply(.SD, 1, function(z) paste(cols[as.logical(z)], collapse = " + ")),
      .SDcols = cols]
  out[, evidence := evidence_type]
  out[]
}

mapped_lists <- lapply(jobs, function(j) unique(na.omit(as.character(j$genes$symbol))))
magma_lists <- lapply(jobs, function(j) {
  x <- copy(j$magma_genes)
  x[, P := as.numeric(P)]
  x <- x[P < 0.05 / sum(is.finite(P))]
  unique(na.omit(ifelse(!is.na(x$SYMBOL) & nzchar(x$SYMBOL), as.character(x$SYMBOL), as.character(x$GENE))))
})

mapped_membership <- membership_table(mapped_lists, "Positionally mapped")
magma_membership <- membership_table(magma_lists, "MAGMA Bonferroni")
gene_membership <- rbindlist(list(mapped_membership, magma_membership), fill = TRUE)
fwrite(gene_membership, file.path(tab_dir, "fuma_gene_membership_common_affective_functional.csv"))

membership_counts <- gene_membership[, .N, by = .(evidence, membership)]
setorder(membership_counts, evidence, -N)
fwrite(membership_counts, file.path(tab_dir, "fuma_gene_membership_counts.csv"))

p7 <- ggplot(membership_counts, aes(reorder(membership, N), N, fill = evidence)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = N), hjust = -0.08, size = 3) +
  coord_flip() +
  facet_wrap(~ evidence, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Positionally mapped" = "#4C78A8", "MAGMA Bonferroni" = "#D95F02"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "G  Gene evidence shared across genetic models",
       subtitle = "Positional mapping and MAGMA significance are displayed as separate evidence layers",
       x = NULL, y = "Number of genes") +
  theme_nc(9)
save_plot(p7, "figure_4A_cross_domain_gene_membership", 9.5, 7.2)

## Pairwise interval overlap --------------------------------------------------

interval_overlap <- function(a, b, name_a, name_b) {
  a <- as.data.table(a)[, .(locus_a = GenomicLocus, chr_a = as.integer(chr),
                            start_a = as.numeric(start), end_a = as.numeric(end))]
  b <- as.data.table(b)[, .(locus_b = GenomicLocus, chr_b = as.integer(chr),
                            start_b = as.numeric(start), end_b = as.numeric(end))]
  z <- merge(a, b, by.x = "chr_a", by.y = "chr_b", allow.cartesian = TRUE)
  z <- z[pmax(start_a, start_b) <= pmin(end_a, end_b)]
  z[, `:=`(
    comparison = paste(name_a, name_b, sep = " vs "),
    overlap_start = pmax(start_a, start_b),
    overlap_end = pmin(end_a, end_b),
    overlap_bp = pmin(end_a, end_b) - pmax(start_a, start_b) + 1
  )]
  z[]
}

pair_names <- combn(names(jobs), 2, simplify = FALSE)
overlap_intervals <- rbindlist(lapply(pair_names, function(z) {
  interval_overlap(jobs[[z[1]]]$loci, jobs[[z[2]]]$loci, z[1], z[2])
}), fill = TRUE)
fwrite(overlap_intervals, file.path(tab_dir, "fuma_pairwise_overlapping_risk_intervals.csv"))

overlap_summary <- rbindlist(lapply(pair_names, function(z) {
  x <- interval_overlap(jobs[[z[1]]]$loci, jobs[[z[2]]]$loci, z[1], z[2])
  data.table(
    comparison = paste(z[1], z[2], sep = " vs "),
    overlapping_interval_pairs = nrow(x),
    loci_from_first_with_overlap = uniqueN(x$locus_a),
    loci_from_second_with_overlap = uniqueN(x$locus_b)
  )
}))
fwrite(overlap_summary, file.path(tab_dir, "fuma_pairwise_risk_interval_overlap_summary.csv"))

p8 <- ggplot(overlap_summary, aes(comparison, overlapping_interval_pairs)) +
  geom_col(fill = "#5B6573", width = 0.65) +
  geom_text(aes(label = overlapping_interval_pairs), vjust = -0.35, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "H  Pairwise overlap of FUMA genomic risk intervals",
       subtitle = "Counts describe intersecting interval pairs and do not imply identical causal variants",
       x = NULL, y = "Overlapping interval pairs") +
  theme_nc(9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p8, "figure_4B_pairwise_risk_interval_overlap", 8.5, 4.8)

## Optional composite ---------------------------------------------------------

if (has_patchwork) {
  composite <- (p1 / p2 / (p3 | p4) / (p6 | p7) / p8) +
    patchwork::plot_layout(heights = c(0.8, 1.0, 1.35, 1.0, 0.7)) +
    patchwork::plot_annotation(
      title = "FUMA annotation of common, affective, and functional genetic factors",
      subtitle = "SNP2GENE positional mapping, MAGMA analyses, functional annotation, and cross-model comparison",
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "grey35"))
    )
  save_plot(composite, "figure_fuma_job06_07_composite_overview", 18, 27)
}

## Machine-readable run manifest and concise report --------------------------

manifest <- data.table(
  item = c("script", "project_base", "output_directory", "common_comparator_available",
           "job06_directory", "job07_directory", "generated_at", "R_version"),
  value = c("build_fuma_job06_07_visuals.R", base_dir, out_dir,
            as.character(common_available), job_dirs[["Affective"]], job_dirs[["Functional"]],
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), R.version.string)
)
fwrite(manifest, file.path(out_dir, "run_manifest.csv"))

report_lines <- c(
  "FUMA JOB06 and JOB07 visualization audit",
  "",
  paste0("Output directory: ", out_dir),
  paste0("Common-factor comparator included: ", common_available),
  "",
  apply(as.data.frame(summary_dt), 1, function(z) {
    paste0(z[[1]], ": ", z[[2]], " genomic risk loci; ", z[[5]], " lead SNPs; ",
           " ", z[[6]], " independent significant SNPs; ", z[[8]],
           " positionally mapped genes; ", z[[11]], " MAGMA Bonferroni genes.")
  }),
  "",
  "Interpretation constraints:",
  "Positionally mapped genes and MAGMA significant genes are separate evidence classes.",
  "Risk-interval overlap does not establish an identical causal variant or colocalization.",
  "FUMA parameters and MHC handling are recorded in fuma_parameter_audit.csv.",
  "Gene sets should be frozen only after inspection of all domain-specific tables and sensitivity analyses."
)
writeLines(report_lines, file.path(out_dir, "fuma_job06_07_audit_report.txt"))

message("Completed FUMA JOB06/JOB07 visualization package: ", out_dir)




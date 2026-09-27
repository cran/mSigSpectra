## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 8,
  fig.height = 3.5,
  dpi = 100,
  out.width = "100%"
)
can_load <- function(pkg) {
  nzchar(system.file(package = pkg)) &&
    isTRUE(tryCatch({ loadNamespace(pkg); TRUE }, error = function(e) FALSE))
}
have_bsgenome <- can_load("BSgenome.Hsapiens.1000genomes.hs37d5")
have_msigplot <- can_load("mSigPlot")
eval_all <- have_bsgenome && have_msigplot
knitr::opts_chunk$set(eval = eval_all)

## ----load---------------------------------------------------------------------
library(mSigSpectra)
library(mSigPlot)

## ----files--------------------------------------------------------------------
sbs_file <- system.file(
  "extdata", "Strelka-SBS-GRCh37", "Strelka.SBS.GRCh37.s1.vcf",
  package = "mSigSpectra"
)
id_file <- system.file(
  "extdata", "Strelka-ID-GRCh37", "Strelka.ID.GRCh37.s1.vcf",
  package = "mSigSpectra"
)
mutect_file <- system.file(
  "extdata", "Mutect-GRCh37", "Mutect.GRCh37.s1.vcf",
  package = "mSigSpectra"
)

## ----read-sbs-----------------------------------------------------------------
sbs_vcf <- read_vcf(sbs_file, filter = "PASS")
nrow(sbs_vcf)
head(sbs_vcf[, c("CHROM", "POS", "REF", "ALT", "FILTER")])

## ----split-sbs----------------------------------------------------------------
sbs_split <- split_vcf(sbs_vcf, name_of_vcf = "Strelka.SBS.GRCh37.s1")
sapply(sbs_split[c("SBS", "DBS", "ID")], nrow)

## ----annotate-sbs-------------------------------------------------------------
sbs_ann <- annotate_sbs_or_dbs_vcf(sbs_split$SBS,
                                   ref_genome = "GRCh37")$annotated.vcf
new_cols <- setdiff(colnames(sbs_ann), colnames(sbs_split$SBS))
new_cols

## ----build-sbs----------------------------------------------------------------
cat96 <- vcf_to_sbs_catalog(sbs_ann, type = "SBS96", ref_genome = "GRCh37",
                            region = "genome", sample_name = "s1")
cat192 <- vcf_to_sbs_catalog(sbs_ann, type = "SBS192", ref_genome = "GRCh37",
                             region = "transcript", sample_name = "s1")
cat1536 <- vcf_to_sbs_catalog(sbs_ann, type = "SBS1536", ref_genome = "GRCh37",
                              region = "genome", sample_name = "s1")

dim(cat96)
sum(cat96)
catalog_attrs(cat96)

## ----plot-sbs96, fig.height=3.5-----------------------------------------------
plot_SBS96(cat96, plot_title = "Strelka.SBS.GRCh37.s1 — SBS96")

## ----plot-sbs192, fig.height=4------------------------------------------------
plot_SBS192(cat192, plot_title = "Strelka.SBS.GRCh37.s1 — SBS192 (stranded)")

## ----plot-sbs1536, fig.height=8-----------------------------------------------
plot_SBS1536(cat1536, plot_title = "Strelka.SBS.GRCh37.s1 — SBS1536")

## ----id-pipeline--------------------------------------------------------------
id_vcf <- read_vcf(id_file, filter = "PASS")
id_split <- split_vcf(id_vcf, name_of_vcf = "Strelka.ID.GRCh37.s1")
id_ann <- annotate_id_vcf(id_split$ID,
                          ref_genome = "GRCh37")$annotated.vcf
id_ann[1, c("CHROM", "POS", "REF", "ALT", "COSMIC_83", "Koh_89")]

cat_id83 <- vcf_to_id_catalog(id_ann, type = "ID83", ref_genome = "GRCh37",
                              region = "genome", sample_name = "s1")
cat_id89 <- vcf_to_id_catalog(id_ann, type = "ID89", ref_genome = "GRCh37",
                              region = "genome", sample_name = "s1")
cat_id476 <- vcf_to_id_catalog(id_ann, type = "ID476", ref_genome = "GRCh37",
                               region = "genome", sample_name = "s1")

## ----plot-id83, fig.height=3.5------------------------------------------------
plot_ID83(cat_id83, plot_title = "Strelka.ID.GRCh37.s1 — ID83")

## ----plot-id89, fig.height=4--------------------------------------------------
plot_ID89(cat_id89, plot_title = "Strelka.ID.GRCh37.s1 — ID89")

## ----plot-id476, fig.height=10------------------------------------------------
plot_ID476(cat_id476, plot_title = "Strelka.ID.GRCh37.s1 — ID476")

## ----density------------------------------------------------------------------
cat96_density <- transform_catalog(cat96,
                                   target_counts_or_density = "density")
attr(cat96_density, "counts_or_density")

## ----plot-sbs96-density, fig.height=3.5---------------------------------------
plot_SBS96(cat96_density,
           plot_title = "Strelka.SBS.GRCh37.s1 — SBS96 (density)")

## ----collapse-----------------------------------------------------------------
cat96_from_1536 <- collapse_catalog(cat1536, to = "SBS96")
all.equal(as.numeric(cat96_from_1536[, 1]),
          as.numeric(cat96[, 1]))

## ----io-----------------------------------------------------------------------
out_path <- file.path(tempdir(), "Strelka.SBS.GRCh37.s1.SBS96.csv")
write_catalog(cat96, out_path)
cat96_back <- read_catalog(out_path, region = "genome")
identical(as.numeric(cat96), as.numeric(cat96_back))


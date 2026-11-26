# Five DK34 spatial-rotation tests with multiplicity control
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(freesurferformats)
  library(clue)
})

root <- study_path("project")
tau_file <- file.path(root, "results/13H3_tau_spatial_specificity_release",
                      "04_regional_tau_and_frozen_genetic_map.csv")
dirs <- list.dirs(file.path(root, "results"), recursive = FALSE, full.names = TRUE)
dirs <- sort(dirs[grepl("^13H1_aligned_regional_slope_map_", basename(dirs))],
             decreasing = TRUE)
map_dir <- dirs[vapply(dirs, function(d) {
  f <- file.path(d, "01_aligned_34_region_contrasts.csv")
  file.exists(f) && "gAD_estimate" %in% names(fread(f, nrows = 0L))
}, logical(1))][1L]
if (is.na(map_dir) || !file.exists(tau_file))
  stop("Run the updated 13H1 script first; the frozen 13H3 tau table is also required.")

map <- fread(file.path(map_dir, "01_aligned_34_region_contrasts.csv"))
tau <- fread(tau_file)[, .(region, tau_AD_minus_CN)]
map <- merge(map, tau, by = "region", sort = TRUE)
stopifnot(nrow(map) == 34L, all(map$status == "SUCCESS"),
          max(abs(map$estimate - (map$gAD_estimate - map$gMAV_estimate))) < 1e-10,
          all(is.finite(map$tau_AD_minus_CN)))

# Use matched fsaverage5/DK geometry and annotations from the netneurotools manifest.
atlas_dir <- file.path(root, "input/13H20_fsaverage5_DK_geometry")
dir.create(atlas_dir, recursive = TRUE, showWarnings = FALSE)
archives <- data.table(
  name = c("template", "aparc"),
  url = c(
    "https://files.osf.io/v1/resources/udpv8/providers/osfstorage/67327116b8b930019c62fd6b",
    "https://files.osf.io/v1/resources/udpv8/providers/osfstorage/69ccfbab936bb1b6e5744dc8"
  ),
  md5 = c("b79f2583917e2a5f79a4768e2b6e27e6",
          "3aa7d3a6c4fc8c9a68b889f9defff099")
)
for (i in seq_len(nrow(archives))) {
  archive <- file.path(atlas_dir, paste0(archives$name[i], ".tar.gz"))
  if (!file.exists(archive))
    download.file(archives$url[i], archive, mode = "wb", method = "libcurl", quiet = TRUE)
  if (unname(tools::md5sum(archive)) != archives$md5[i])
    stop("Atlas checksum mismatch: ", archive)
  target <- file.path(atlas_dir, archives$name[i])
  dir.create(target, showWarnings = FALSE)
  if (!length(list.files(target, recursive = TRUE))) untar(archive, exdir = target)
}
find_one <- function(folder, filename) {
  matches <- list.files(folder, recursive = TRUE, full.names = TRUE)
  matches <- matches[basename(matches) == filename]
  if (length(matches) != 1L) stop("Expected exactly one ", filename, " in ", folder)
  matches
}
sphere_file <- find_one(file.path(atlas_dir, "template"), "lh.sphere")
annot_file <- find_one(file.path(atlas_dir, "aparc"), "lh.aparc.annot")
sphere <- freesurferformats::read.fs.surface(sphere_file)$vertices
labels <- freesurferformats::read.fs.annot(annot_file)$label_names
stopifnot(nrow(sphere) == length(labels), all(map$region %in% labels))
centroids <- t(vapply(map$region, function(region) {
  xyz <- colMeans(sphere[labels == region, , drop = FALSE])
  xyz / sqrt(sum(xyz^2))
}, numeric(3)))
geometry <- data.table(region = map$region, x = centroids[, 1L],
                       y = centroids[, 2L], z = centroids[, 3L],
                       n_vertices = vapply(map$region, function(x) sum(labels == x), integer(1)))

rotation <- function() {
  q <- rnorm(4L)
  q <- q / sqrt(sum(q^2))
  a <- q[1L]; b <- q[2L]; c <- q[3L]; d <- q[4L]
  matrix(c(1 - 2 * (c^2 + d^2), 2 * (b * c - a * d), 2 * (b * d + a * c),
           2 * (b * c + a * d), 1 - 2 * (b^2 + d^2), 2 * (c * d - a * b),
           2 * (b * d - a * c), 2 * (c * d + a * b), 1 - 2 * (b^2 + c^2)),
         nrow = 3L, byrow = TRUE)
}
spin_index <- function() {
  similarity <- (centroids %*% t(rotation())) %*% t(centroids)
  as.integer(clue::solve_LSAP(max(similarity) - similarity))
}
spin <- function(x, assignment) {
  y <- numeric(length(x))
  y[assignment] <- x
  y
}
align <- function(x, y) unname(cor(x, y))

gAD <- -map$gAD_estimate
gMAV <- -map$gMAV_estimate
excess_gAD <- gAD - gMAV
tau_map <- map$tau_AD_minus_CN
non_temporal <- !map$temporal_target
stat <- function(a, m, e) c(
  gAD = align(a, tau_map),
  gMAV = align(m, tau_map),
  gAD_excess = align(e, tau_map),
  gAD_minus_gMAV_alignment = align(a, tau_map) - align(m, tau_map),
  gAD_excess_non_temporal = align(e[non_temporal], tau_map[non_temporal])
)
observed <- stat(gAD, gMAV, excess_gAD)
n_spin <- 10000L
set.seed(20261021)
null <- matrix(NA_real_, nrow = n_spin, ncol = length(observed),
               dimnames = list(NULL, names(observed)))
for (i in seq_len(n_spin)) {
  idx <- spin_index()
  null[i, ] <- stat(spin(gAD, idx), spin(gMAV, idx), spin(excess_gAD, idx))
}
tests <- data.table(test = names(observed), observed_r = as.numeric(observed),
                    p_spin_two_sided = vapply(seq_along(observed), function(j)
                      (1 + sum(abs(null[, j]) >= abs(observed[j]))) / (n_spin + 1),
                      numeric(1)))
tests[, q_fdr := p.adjust(p_spin_two_sided, method = "BH")]
map[, `:=`(gAD_thinning = gAD, gMAV_thinning = gMAV,
           gAD_excess_thinning = excess_gAD)]

out <- file.path(root, "results", paste0("13H20_gAD_gMAV_tau_DK34_spin_",
                                           "release"))
dir.create(out, recursive = TRUE)
fwrite(map, file.path(out, "01_regional_genetic_and_tau_maps.csv"), bom = TRUE)
fwrite(tests, file.path(out, "02_spatial_null_tests.csv"), bom = TRUE)
fwrite(as.data.table(null), file.path(out, "03_spin_null_distribution.csv"))
fwrite(geometry, file.path(out, "04_DK34_spherical_centroids.csv"))
fwrite(archives, file.path(out, "05_atlas_source_and_checksums.csv"))

status <- c(
  "STATUS=SEPARATE_GENETIC_TAU_SPATIAL_NULL_COMPLETE",
  paste0("MAP_DIR=", map_dir),
  paste0("SPATIAL_NULL=FSAVERAGE5_LEFT_DK34_SPHERICAL_ROTATION_", n_spin),
  paste0("GENETIC_MAP_CORRELATION=", signif(cor(gAD, gMAV), 5)),
  paste0("gAD_EXCESS_TAU_R=", signif(observed["gAD_excess"], 5)),
  paste0("gAD_EXCESS_TAU_SPIN_P=", signif(tests[test == "gAD_excess", p_spin_two_sided], 5)),
  paste0("gAD_VS_gMAV_ALIGNMENT_SPIN_P=",
         signif(tests[test == "gAD_minus_gMAV_alignment", p_spin_two_sided], 5)),
  "NOTE=DK34 is a bilateral average; one hemisphere supplies homotopic spherical geometry.",
  "NOTE=Regional tau values are the age/sex-adjusted FTP PET contrasts.",
  paste0("OUTPUT_DIR=", out)
)
writeLines(status, file.path(out, "00_status.txt"), useBytes = TRUE)
message(paste(status, collapse = "\n"))

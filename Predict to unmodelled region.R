##############################
# Integrated density extrapolation to new region
# Outputs:
#  - Pixel-level posterior mean, SD, 50%, 95% CI of density
#  - Pixel-level expected abundance (animals per pixel)
#  - Region-level total abundance posterior (vector), summary + violin plot
##############################


load("fisher_data_extrapolate.RData")

# ---- Packages ----
library(terra)
library(coda)
library(future.apply)
library(ggplot2)

# ---- Inputs you must set ----

raster_path   <- "I:/Ecosystems/Conservation Science/Species/Mesocarnivores/Projects/MMP/2.Data/3. SDM-Cindy/5. GIS/Fisher Habitat Rasters/Landsat_Lidar/canopy_resampled120.tif"  # clipped to region of interest
density_units <- "per_km2"                  # set to "per_km2" (typical) or "per_m2"
n_thin        <- 100                        # thinned posterior draws for CI/violin
workers       <- 4                          # parallel workers
blocksize     <- 50000                      # pixels per block for CI & region sum

# ---- Posterior samples ----
mcmc_samps <- mvSamples                    # or as.mcmc(Cmcmc$collectSamples())
par_names  <- colnames(mcmc_samps)

# Robust parameter extraction
Dbeta1_idx <- which(par_names == "D.beta1")
if (length(Dbeta1_idx) != 1) {
  stop("Could not find exactly one 'D.beta1' column. Found: ", paste(par_names[grep("D\\.beta1", par_names)], collapse = ", "))
}
D0_idx <- which(par_names == "D0")
if (length(D0_idx) == 0) {
  stop("No session intercept columns matching '^D0\\[' found. Inspect names and adjust the pattern.")
}

D.beta1 <- as.numeric(mcmc_samps[, Dbeta1_idx])
D0_mat  <- as.matrix(mcmc_samps[, D0_idx, drop = FALSE])
n_samp  <- nrow(mcmc_samps)

# Population-level intercept per draw (mean across sessions)
alpha.new <- rowMeans(D0_mat)

# ---- Covariate raster + standardization ----
rs_Dcov    <- rast(raster_path)
Dcov_vals  <- values(rs_Dcov)

if (!exists("mean.D.cov") || !exists("sd.D.cov")) {
  stop("mean.D.cov and sd.D.cov must be loaded from your training script/environment.")
}
if (!is.finite(sd.D.cov) || sd.D.cov == 0) {
  warning("sd.D.cov is non-finite or zero; using centered (unscaled) covariate.")
  Dcov_std <- Dcov_vals - mean.D.cov
} else {
  Dcov_std <- (Dcov_vals - mean.D.cov) / sd.D.cov
}

n_pix  <- length(Dcov_std)
valid  <- is.finite(Dcov_std)
v_ids  <- which(valid)
n_valid <- length(v_ids)

# ---- Diagnostics ----
cat("Draws:", n_samp, "\n",
    "Raster cells:", n_pix, "\n",
    "Valid cells:", n_valid, "\n",
    "Any NA in alpha.new?", anyNA(alpha.new), "\n",
    "Any NA in D.beta1?",   anyNA(D.beta1), "\n")

# ---- Pixel-level posterior mean & variance of density ----
sum_D  <- numeric(n_pix)
sum_D2 <- numeric(n_pix)
count  <- integer(n_pix)
skipped_draws <- 0L

# Use all draws; optionally thin here for speed: thin_all <- sample(seq_len(n_samp), min(1000, n_samp))
for (s in seq_len(n_samp)) {
  a <- alpha.new[s]; b <- D.beta1[s]
  if (!is.finite(a) || !is.finite(b)) { skipped_draws <- skipped_draws + 1L; next }
  
  eta_v <- a + b * Dcov_std[valid]
  D_v   <- exp(eta_v)
  
  # Replace non-finite predictions
  good <- is.finite(D_v)
  if (!all(good)) D_v[!good] <- NA_real_
  
  idx_g <- v_ids[good]
  sum_D [idx_g] <- sum_D [idx_g] + D_v[good]
  sum_D2[idx_g] <- sum_D2[idx_g] + D_v[good]^2
  count [idx_g] <- count [idx_g] + 1L
}
cat("Skipped draws (non-finite alpha/beta):", skipped_draws, "\n")

D_mean <- rep(NA_real_, n_pix)
D_var  <- rep(NA_real_, n_pix)
has_draws <- count > 0
D_mean[has_draws] <- sum_D[has_draws] / count[has_draws]
D_var [has_draws] <- pmax(0, (sum_D2[has_draws] / count[has_draws]) - D_mean[has_draws]^2)
D_sd <- rep(NA_real_, n_pix); D_sd[has_draws] <- sqrt(D_var[has_draws])

# ---- Pixel-level CI (50% median, 2.5–97.5%) ----
options(future.globals.maxSize = 2 * 1024^3)
plan(multisession, workers = workers)

set.seed(20251216)
thin_idx <- sample(seq_len(n_samp), size = min(n_thin, n_samp))
alpha_thin <- alpha.new[thin_idx]
beta_thin  <- D.beta1[thin_idx]

D_low  <- rep(NA_real_, n_pix)
D_high <- rep(NA_real_, n_pix)
D_med  <- rep(NA_real_, n_pix)

for (start in seq(1, n_valid, by = blocksize)) {
  end        <- min(start + blocksize - 1, n_valid)
  idx_block  <- v_ids[start:end]
  x_block    <- Dcov_std[idx_block]
  
  pred_block_list <- future_mapply(
    FUN = function(a, b, x) {
      if (!is.finite(a) || !is.finite(b)) return(rep(NA_real_, length(x)))
      out <- exp(a + b * x)
      out[!is.finite(out)] <- NA_real_
      out
    },
    a = alpha_thin, b = beta_thin,
    MoreArgs = list(x = x_block),
    SIMPLIFY = FALSE
  )
  pred_block <- do.call(cbind, pred_block_list)  # [n_block x n_thin]
  
  # Row-wise quantiles
  D_low [idx_block] <- apply(pred_block, 1, function(x) { x <- x[is.finite(x)]; if (!length(x)) return(NA_real_); quantile(x, 0.025, names = FALSE) })
  D_high[idx_block] <- apply(pred_block, 1, function(x) { x <- x[is.finite(x)]; if (!length(x)) return(NA_real_); quantile(x, 0.975, names = FALSE) })
  D_med [idx_block] <- apply(pred_block, 1, function(x) { x <- x[is.finite(x)]; if (!length(x)) return(NA_real_); quantile(x, 0.50,  names = FALSE) })
}
plan(sequential)

# ---- Compute pixel area (from raster resolution) ----
is_deg <- terra::is.lonlat(rs_Dcov)
if (is_deg) stop("CRS is geographic (degrees). Reproject rs_Dcov to an equal-area CRS (meters) before computing areas.")

cell_area_m2  <- prod(res(rs_Dcov))        # m² per cell
cell_area_km2 <- cell_area_m2 / 1e6

# Expected abundance per pixel (E[D] × area)
if (density_units == "per_km2") {
  abund_pixel_mean <- D_mean * cell_area_km2
} else if (density_units == "per_m2") {
  abund_pixel_mean <- D_mean * cell_area_m2
} else {
  stop("density_units must be 'per_km2' or 'per_m2'.")
}

# ---- Region-level total abundance per draw (posterior distribution) ----
# Re-use the same thinning; sum density over valid pixels and multiply by area
region_sum <- numeric(length(alpha_thin))

plan(multisession, workers = workers)
for (start in seq(1, n_valid, by = blocksize)) {
  end        <- min(start + blocksize - 1, n_valid)
  idx_block  <- v_ids[start:end]
  x_block    <- Dcov_std[idx_block]
  
  block_sums <- future_mapply(
    FUN = function(a, b, x) {
      if (!is.finite(a) || !is.finite(b)) return(0)
      d <- exp(a + b * x)
      d[!is.finite(d)] <- 0
      sum(d, na.rm = TRUE)
    },
    a = alpha_thin, b = beta_thin,
    MoreArgs = list(x = x_block),
    SIMPLIFY = TRUE
  )
  region_sum <- region_sum + block_sums
}
plan(sequential)

if (density_units == "per_km2") {
  region_abundance <- region_sum * cell_area_km2
} else {
  region_abundance <- region_sum * cell_area_m2
}

# ---- Summaries of total abundance ----
qs <- quantile(region_abundance, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
cat("\nTotal abundance (posterior):\n",
    "Median:", qs[3], "\n",
    "50% CI:", qs[2], "–", qs[4], "\n",
    "95% CI:", qs[1], "–", qs[5], "\n")

# ---- Save rasters ----
D_mean_r <- rs_Dcov; values(D_mean_r) <- D_mean; names(D_mean_r) <- "D_mean"
D_sd_r   <- rs_Dcov; values(D_sd_r)   <- D_sd;   names(D_sd_r)   <- "D_sd"
D_low_r  <- rs_Dcov; values(D_low_r)  <- D_low;  names(D_low_r)  <- "D_low_95"
D_high_r <- rs_Dcov; values(D_high_r) <- D_high; names(D_high_r) <- "D_high_95"
D_med_r  <- rs_Dcov; values(D_med_r)  <- D_med;  names(D_med_r)  <- "D_median"
Abund_pixel_r <- rs_Dcov; values(Abund_pixel_r) <- abund_pixel_mean; names(Abund_pixel_r) <- "Abund_pixel_mean"

writeRaster(D_mean_r,       "density_mean.tif",         overwrite = TRUE)
writeRaster(D_sd_r,         "density_sd.tif",           overwrite = TRUE)
writeRaster(D_low_r,        "density_low95.tif",        overwrite = TRUE)
writeRaster(D_high_r,       "density_high95.tif",       overwrite = TRUE)
writeRaster(D_med_r,        "density_median.tif",       overwrite = TRUE)
writeRaster(Abund_pixel_r,  "abundance_pixel_mean.tif", overwrite = TRUE)

# ---- Quick maps ----
plot(D_mean_r,      main = "Posterior mean density")
plot(D_med_r,       main = "Posterior median density")
plot(D_low_r,       main = "Lower 95% CI (density)")
plot(D_high_r,      main = "Upper 95% CI (density)")
plot(Abund_pixel_r, main = "Expected abundance per pixel")

# ---- Violin plot: total abundance posterior ----
df <- data.frame(total_abundance = region_abundance)
qs2 <- quantile(df$total_abundance, c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm = TRUE)
p <- ggplot(df, aes(x = "New region", y = total_abundance)) +
  geom_violin(fill = "#5DADE2", color = "#2E86C1", width = 0.6, adjust = 1) +
  annotate("segment", x = 1, xend = 1, y = qs2[1], yend = qs2[5], color = "#C0392B", size = 1.2) +  # 95% CI
  annotate("segment", x = 1, xend = 1, y = qs2[2], yend = qs2[4], color = "#1F618D", size = 3.5) +  # 50% CI
  annotate("point",   x = 1, y = qs2[3], shape = 21, fill = "white", color = "black", size = 3) +   # median
  labs(
    title    = "Total abundance in the new region",
    subtitle = "Posterior distribution across thinned draws (median •, 50% CI, 95% CI)",
    x        = NULL,
    y        = "Individuals"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12)
print(p)
ggsave("total_abundance_violin.png", p, width = 6, height = 4, dpi = 300)
       
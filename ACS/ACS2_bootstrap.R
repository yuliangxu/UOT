B = 50
base_seed = 123
n_pairs = 1000

reg_m_kl <- 1e-2
reg      <- 1e-2

var1d = "AGE"
# 2-1) bootstrap ----------------------------------------------------------


# run (1) and (2) in ACS1...R file first.
Gb_mat <- matrix(NA_real_, nrow = B, ncol = K * K)
pb <- txtProgressBar(min = 0, max = B, style = 3)
for (b in 1:B) {
  setTxtProgressBar(pb, b)
  
  seed_i <- base_seed + b
  
  batch <- bootstrap_couples_one(
    male_df, female_df,
    n_pairs = n_pairs,
    seed = seed_i
  )
  
  # extract 1D covariate from this bootstrap sample
  x1 <- batch$male[, var1d]
  x0 <- batch$female[, var1d]
  
  if (var1d == "INCWAGE") {
    x1 <- log1p(pmax(x1, 0))
    x0 <- log1p(pmax(x0, 0))
  }
  
  # assign to *global* bins
  bin1 <- cut(x1, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bin0 <- cut(x0, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # binned marginals
  a1_bin <- rep(0, K)
  a0_bin <- rep(0, K)
  
  tmp1 <- tapply(batch$w_male,   bin1, sum)
  tmp0 <- tapply(batch$w_female, bin0, sum)
  
  if (!is.null(tmp1)) a1_bin[as.integer(names(tmp1))] <- as.numeric(tmp1)
  if (!is.null(tmp0)) a0_bin[as.integer(names(tmp0))] <- as.numeric(tmp0)
  
  # normalize to distributions
  if (sum(a1_bin) == 0 || sum(a0_bin) == 0) next
  a1_bin <- a1_bin / sum(a1_bin)
  a0_bin <- a0_bin / sum(a0_bin)
  
  # OT on bins -> K x K plan
  G_bin <- ot$sinkhorn_unbalanced(a1_bin, a0_bin, M_bin, reg, reg_m_kl, div = "kl")
  
  # store vectorized
  Gb_mat[b, ] <- as.numeric(G_bin)
}

# drop any failed reps
keep <- which(rowSums(is.na(Gb_mat)) == 0)
Gb_mat <- Gb_mat[keep, , drop = FALSE]
B_eff <- nrow(Gb_mat)
B_eff


# 2-2)normality check -----------------------------------------------------

# 3 plots
library(ggplot2)

df <- expand.grid(i=1:K, j=1:K)
df$Mean <- as.vector(G_mean)
df$`2.5%` <- as.vector(G_q025)
df$`97.5%` <- as.vector(G_q975)

df_long <- reshape2::melt(df, id.vars = c("i","j"),
                          variable.name = "stat", value.name = "val")

ggplot(df_long, aes(x=j, y=i, fill=val)) +
  geom_tile() +
  scale_y_reverse() +
  facet_wrap(~stat, nrow = 1) +
  labs(title = sprintf("UOT (eps=%.3g, rho=%.3g)", reg, reg_m_kl),
       x="Target bin", y="Source bin", fill="") +
  theme_minimal()


## PCA-based normality checks -----------
# center columns
Gb_centered <- scale(Gb_mat, center = TRUE, scale = FALSE)

# PCA (use prcomp; for large K, consider irlba)
pca <- prcomp(Gb_centered, rank. = 3)

scores <- pca$x  # B_eff x 3 matrix of PC scores

par(mfrow = c(3, 2))
for (r in 1:3) {
  z <- sqrt(n) * scores[, r]  # scaling
  qqnorm(z, main = sprintf("QQ: sqrt(n)*PC%d score", r)); qqline(z)
  hist(z, breaks = 30, main = sprintf("Hist: sqrt(n)*PC%d score", r),
       xlab = "value")
}
par(mfrow = c(1, 1))


# 2-3) run for grids of rho and eps ---------------------------------------

library(reshape2)
library(ggplot2)

# choose your grid
eps_list  <- c(1e-3, 1e-2, 1e-1)      # reg
rho_list  <- c(1e-2, 1e-1, 1, 10)     # reg_m_kl

grid <- expand.grid(eps = eps_list, rho = rho_list, KEEP.OUT.ATTRS = FALSE)

all_long <- list()

pb <- txtProgressBar(min = 0, max = nrow(grid), style = 3)

for (g in 1:nrow(grid)) {
  setTxtProgressBar(pb, g)
  
  reg      <- grid$eps[g]
  reg_m_kl <- grid$rho[g]
  
  out <- boot_Gbin_summary(
    reg = reg, reg_m_kl = reg_m_kl,
    B = B, base_seed = base_seed, n_pairs = n_pairs,
    male_df = male_df, female_df = female_df,
    var1d = var1d, breaks = breaks, K = K,
    M_bin = M_bin, ot = ot
  )
  
  # build KxK matrices
  G_mean <- matrix(out$mean, nrow = K, ncol = K)
  G_q025 <- matrix(out$q025, nrow = K, ncol = K)
  G_q975 <- matrix(out$q975, nrow = K, ncol = K)
  
  df <- expand.grid(i = 1:K, j = 1:K)
  df$Mean    <- as.vector(G_mean)
  df$`2.5%`  <- as.vector(G_q025)
  df$`97.5%` <- as.vector(G_q975)
  
  df_long <- reshape2::melt(df, id.vars = c("i","j"),
                            variable.name = "stat", value.name = "val")
  
  df_long$eps <- reg
  df_long$rho <- reg_m_kl
  df_long$B_eff <- out$B_eff
  
  all_long[[g]] <- df_long
}

df_all <- do.call(rbind, all_long)

# make nice facet labels
df_all$eps_f <- factor(df_all$eps, levels = sort(unique(df_all$eps)))
df_all$rho_f <- factor(df_all$rho, levels = sort(unique(df_all$rho)))

# global max across all eps/rho/stat so colors are comparable
vmax <- max(df_all$val, na.rm = TRUE)

ggplot(df_all, aes(x = j, y = i, fill = val)) +
  geom_tile() +
  scale_y_reverse() +
  facet_grid(rho_f ~ eps_f + stat, scales = "fixed") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, vmax),
    oob = scales::squish
  ) +
  labs(
    title = sprintf("UOT bootstrap summaries (var=%s, n_pairs=%d, B=%d)", var1d, n_pairs, B),
    x = "Target bin", y = "Source bin", fill = ""
  ) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 9)
  )



# 2-4) change bin label, finalize the plot --------------------------------
stride <- 5  # try 3, 4, 5 depending on K
stopifnot(length(breaks) == K + 1)

bin_labels <- paste0(breaks[-length(breaks)], "–", breaks[-1])  # e.g. "20–22"

# map i/j (1..K) to labeled bins
df_all$src_bin <- factor(df_all$i, levels = 1:K, labels = bin_labels)  # source = male
df_all$tgt_bin <- factor(df_all$j, levels = 1:K, labels = bin_labels)  # target = female
# show only every `stride`-th label
x_breaks <- levels(df_all$tgt_bin)[seq(1, K, by = stride)]
y_breaks <- levels(df_all$src_bin)[seq(1, K, by = stride)]


pow_trans <- scales::trans_new(
  name = "pow035",
  transform = function(x) x^0.35,
  inverse   = function(x) x^(1/0.35)
)

ggplot(df_all, aes(x = tgt_bin, y = src_bin, fill = val)) +
  geom_tile() +
  scale_y_discrete(
    limits = levels(df_all$src_bin),  # <-- no rev()
    breaks = y_breaks       
  ) +
  scale_x_discrete(
    breaks = x_breaks        # fewer x labels
  ) +
  facet_grid(rho_f ~ eps_f + stat, scales = "fixed",
             labeller = labeller(
               rho_f = function(x) paste0("rho = ", x),
               eps_f = function(x) paste0("eps = ", x)
             )
               ) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, vmax),
    oob = scales::squish,
    trans = pow_trans
  ) +
  labs(
    title = sprintf("UOT bootstrap summaries (var=%s, n_pairs=%d, B=%d)", var1d, n_pairs, B),
    x = "Female age bin",
    y = "Male age bin",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 9),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y  = element_text(size = 7)
  )


# 2-5) check empirical joint distribution ---------------------------------


batch <- bootstrap_couples_one(male_df, female_df, n_pairs = n_pairs, seed = base_seed+1)
library(dplyr)
library(ggplot2)
library(tidyr)

# extract raw ages from the matrices (var1d must be a column name in male_df/female_df)
x1 <- batch$male[, var1d]    # male
x0 <- batch$female[, var1d]  # female

# (optional) same transform you used in OT
if (var1d == "INCWAGE") {
  x1 <- log1p(pmax(x1, 0))
  x0 <- log1p(pmax(x0, 0))
}

bin1 <- cut(x1, breaks = breaks, include.lowest = TRUE, labels = FALSE)
bin0 <- cut(x0, breaks = breaks, include.lowest = TRUE, labels = FALSE)

df_emp <- tibble(
  i = bin1,
  j = bin0,
  w = batch$w_male   # weight per UNIQUE couple
) %>%
  filter(!is.na(i), !is.na(j)) %>%
  group_by(i, j) %>%
  summarise(val = sum(w), .groups = "drop") %>%   # weighted empirical joint mass
  complete(i = 1:K, j = 1:K, fill = list(val = 0)) %>%
  mutate(
    src_bin = factor(i, levels = 1:K, labels = bin_labels),
    tgt_bin = factor(j, levels = 1:K, labels = bin_labels)
  )

ggplot(df_emp, aes(x = tgt_bin, y = src_bin, fill = val)) +
  geom_tile() +
  scale_y_discrete(limits = levels(df_emp$src_bin), breaks = y_breaks) +
  scale_x_discrete(breaks = x_breaks) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = sprintf("Empirical joint (weighted) for one bootstrap batch (var=%s, n_pairs=%d)", var1d, n_pairs),
    x = "Female age bin",
    y = "Male age bin",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  )



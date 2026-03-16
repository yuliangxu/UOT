library(dplyr)
library(fastDummies)
library(progress)
library(tibble)
library(matrixStats)

boot_Gbin_summary <- function(reg, reg_m_kl,
                              B, base_seed, n_pairs,
                              male_df, female_df,
                              var1d, breaks, K,
                              M_bin, ot) {
  Gb_mat <- matrix(NA_real_, nrow = B, ncol = K * K)
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (b in 1:B) {
    setTxtProgressBar(pb, b)
    seed_i <- base_seed + b
    
    batch <- bootstrap_couples_one(male_df, female_df, n_pairs = n_pairs, seed = seed_i)
    
    x1 <- batch$male[, var1d]
    x0 <- batch$female[, var1d]
    if (var1d == "INCWAGE") {
      x1 <- log1p(pmax(x1, 0))
      x0 <- log1p(pmax(x0, 0))
    }
    
    bin1 <- cut(x1, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    bin0 <- cut(x0, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    
    a1_bin <- rep(0, K); a0_bin <- rep(0, K)
    tmp1 <- tapply(batch$w_male,   bin1, sum)
    tmp0 <- tapply(batch$w_female, bin0, sum)
    if (!is.null(tmp1)) a1_bin[as.integer(names(tmp1))] <- as.numeric(tmp1)
    if (!is.null(tmp0)) a0_bin[as.integer(names(tmp0))] <- as.numeric(tmp0)
    
    if (sum(a1_bin) == 0 || sum(a0_bin) == 0) next
    a1_bin <- a1_bin / sum(a1_bin)
    a0_bin <- a0_bin / sum(a0_bin)
    
    G_bin <- ot$sinkhorn_unbalanced(a1_bin, a0_bin, M_bin, reg, reg_m_kl, div = "kl")
    Gb_mat[b, ] <- as.numeric(G_bin)
  }
  close(pb)
  
  # drop failed
  keep <- which(rowSums(is.na(Gb_mat)) == 0)
  Gb_mat <- Gb_mat[keep, , drop = FALSE]
  if (nrow(Gb_mat) < 5) stop("Too few successful bootstrap reps.")
  
  # summaries (fast)
  G_mean <- colMeans(Gb_mat)
  G_q025 <- matrixStats::colQuantiles(Gb_mat, probs = 0.025)
  G_q975 <- matrixStats::colQuantiles(Gb_mat, probs = 0.975)
  
  list(mean = G_mean, q025 = G_q025, q975 = G_q975, B_eff = nrow(Gb_mat))
}


make_spouse_pairs <- function(df,
                              covars = c("AGE","EDUC","INCWAGE","INCTOT","OCC","IND","EMPSTAT","MARST"),
                              id_vars = c("YEAR","STATEFIP"),
                              keep_only_mf = TRUE) {
  needed <- c("SERIAL","PERNUM","SEX","SPLOC", covars, id_vars)
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  
  # Keep only people with spouse present in household
  x <- df %>%
    filter(!is.na(SPLOC), SPLOC > 0) %>%
    mutate(
      spouse_pernum = SPLOC,
      # couple_id identical for both spouses:
      couple_id = paste0(SERIAL, "_", pmin(PERNUM, spouse_pernum), "_", pmax(PERNUM, spouse_pernum))
    )
  
  if (keep_only_mf) {
    # Often SEX is 1=male, 2=female in IPUMS
    x <- x %>% filter(SEX %in% c(1, 2))
  }
  
  # Ensure we have exactly two people per couple_id (spouse pair)
  x <- x %>%
    group_by(couple_id) %>%
    filter(n() == 2) %>%
    ungroup()
  
  males <- x %>%
    filter(SEX == 1) %>%
    select(couple_id, all_of(id_vars), all_of(covars)) %>%
    rename_with(~ paste0("M_", .x), all_of(covars))
  
  females <- x %>%
    filter(SEX == 2) %>%
    select(couple_id, all_of(id_vars), all_of(covars)) %>%
    rename_with(~ paste0("F_", .x), all_of(covars))
  
  # Keep only couples where we have both sides
  couples_wide <- inner_join(males, females, by = c("couple_id", id_vars))
  
  # Two “paired” datasets, aligned row-by-row:
  male_cov  <- couples_wide %>% select(couple_id, all_of(id_vars), starts_with("M_"))
  female_cov<- couples_wide %>% select(couple_id, all_of(id_vars), starts_with("F_"))
  
  list(
    couples_wide = couples_wide,  # one row per couple, both sides
    male_cov     = male_cov,      # one row per couple, male covariates
    female_cov   = female_cov     # one row per couple, female covariates
  )
}


bootstrap_couples_one <- function(male_df, female_df,
                                  id_col = "couple_id",
                                  n_pairs = NULL,
                                  frac_pairs = NULL,
                                  seed = NULL,
                                  return_ids = FALSE) {
  stopifnot(id_col %in% names(male_df), id_col %in% names(female_df))
  stopifnot(xor(is.null(n_pairs), is.null(frac_pairs)))
  
  # common couples + stable ordering
  common_ids <- intersect(male_df[[id_col]], female_df[[id_col]])
  m <- male_df %>% dplyr::filter(.data[[id_col]] %in% common_ids) %>% dplyr::arrange(.data[[id_col]])
  f <- female_df %>% dplyr::filter(.data[[id_col]] %in% common_ids) %>% dplyr::arrange(.data[[id_col]])
  
  ids <- m[[id_col]]
  n <- length(ids)
  
  if (!is.null(frac_pairs)) {
    stopifnot(frac_pairs > 0, frac_pairs <= 1)
    k <- max(1L, floor(frac_pairs * n))
  } else {
    k <- as.integer(n_pairs)
    stopifnot(k >= 1)
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # k draws with replacement
  boot_ids <- sample(ids, size = k, replace = TRUE)
  
  # multiplicities -> empirical weights over UNIQUE sampled couples
  wt_tbl <- tibble::tibble(!!id_col := boot_ids) %>%
    dplyr::count(.data[[id_col]], name = "count") %>%
    dplyr::mutate(weight = count / sum(count)) %>%  # sum(count)=k
    dplyr::arrange(.data[[id_col]])
  
  # join to get covariates, then DROP id/weight columns
  m_emp_df <- wt_tbl %>%
    dplyr::left_join(m, by = id_col) %>%
    dplyr::arrange(.data[[id_col]])
  
  f_emp_df <- wt_tbl %>%
    dplyr::left_join(f, by = id_col) %>%
    dplyr::arrange(.data[[id_col]])
  
  # covariate matrices only
  m_mat <- m_emp_df %>% dplyr::select(-dplyr::all_of(c(id_col, "weight", "count"))) %>% as.matrix()
  f_mat <- f_emp_df %>% dplyr::select(-dplyr::all_of(c(id_col, "weight", "count"))) %>% as.matrix()
  
  # weight vectors aligned with the rows of m_mat / f_mat
  w_male   <- m_emp_df$weight
  w_female <- f_emp_df$weight
  
  out <- list(male = m_mat, female = f_mat, w_male = w_male, w_female = w_female)
  
  if (return_ids) out$boot_ids <- wt_tbl[[id_col]]  # unique sampled ids (aligned with rows)
  
  out
}


prep_ot_mat <- function(X, cont_vars = NULL) {
  X <- as.matrix(X)
  
  # optionally subset/reorder columns if cont_vars provided
  if (!is.null(cont_vars)) {
    if (is.null(colnames(X))) stop("X has no colnames; cannot match cont_vars.")
    miss <- setdiff(cont_vars, colnames(X))
    if (length(miss) > 0) stop("Missing columns in X: ", paste(miss, collapse=", "))
    X <- X[, cont_vars, drop = FALSE]
  }
  
  mu <- colMeans(X, na.rm = TRUE)
  sd <- apply(X, 2, sd, na.rm = TRUE)
  sd[sd == 0 | is.na(sd)] <- 1
  
  Xs <- sweep(X, 2, mu, "-")
  Xs <- sweep(Xs, 2, sd, "/")
  
  list(X = Xs, mu = mu, sd = sd, cols = colnames(X))
}


one_boot_ubG3 <- function(seed) {
  batch <- bootstrap_couples_one(male_df, female_df, n_pairs = n_pairs, seed = seed)
  
  male_proc <- prep_ot_mat(
    batch$male,
    cont_vars = cont_vars
  )
  female_proc <- prep_ot_mat(
    batch$female,
    cont_vars = cont_vars
  )
  
  X1 <- male_proc$X
  X0 <- female_proc$X
  
  n1 <- nrow(X1)
  n0 <- nrow(X0)
  
  a1 <- batch$w_male
  a0 <- batch$w_female
  
  G1 <- rowSums(X1^2)         # length n1
  G0 <- rowSums(X0^2)         # length n0
  
  D2 <- outer(G1, G0, "+") - 2 * (X1 %*% t(X0))  # n1 x n0, squared distances
  D2[D2 < 0] <- 0                                # guard tiny negative from numerics
  M  <- sqrt(D2)
  
  M <- M / max(M)
  
  ## unbalanced sinkhorn
  ubG3 <- ot$sinkhorn_unbalanced(a1, a0, M, reg, reg_m_kl, div = "kl")
  ubG3 <- ubG3 / sum(ubG3)
  
  ## make sure it's a numeric R matrix (reticulate sometimes returns numpy arrays)
  ubG3 <- as.matrix(ubG3)
  ubG3
}
# Regenerate the certified Figure 4 and Figure 5 grid and both heatmaps.
library(ggplot2)

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
# Some Rscript launchers encode spaces in the --file argument.
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
folder <- dirname(normalizePath(script, mustWork = TRUE))
python <- Sys.getenv("UOT_PYTHON", unset = "python3")
args <- commandArgs(trailingOnly = TRUE)
status <- system2(python, c(shQuote(file.path(folder, "solve_grid.py")), vapply(args, shQuote, "")))
if (status != 0L) stop("LaLonde grid calculation failed (exit code ", status, ").")

plot_dir <- file.path(folder, "results")
grid <- read.csv(file.path(plot_dir, "converged_grid.csv"))
grid$rho_f <- factor(grid$rho)
grid$eps_f <- factor(grid$epsilon)
benchmark <- unique(grid$experimental_benchmark)
stopifnot(length(benchmark) == 1, nlevels(grid$rho_f) == 20, nlevels(grid$eps_f) == 15,
          all(table(grid$rho_f, grid$eps_f) == 1), all(grid$max_abs_foc < 1e-10),
          all(grid$median_controls_95 >= 1 & grid$median_controls_95 <= 2871))

style <- theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
axis_labels <- labs(x = expression(rho~"(reg_m_kl)"), y = expression(epsilon~"(reg)"))

# Rank controls separately for each treated participant before taking the median.
# The CSV supplies the minimum count covering 95% of each row's transport mass.
matches <- ggplot(grid, aes(rho_f, eps_f, fill = median_controls_95)) +
    geom_tile(color = "white") +
    geom_text(aes(label = median_controls_95, colour = median_controls_95 < 55),
              size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    scale_fill_viridis_c(option = "C", transform = "log10",
                         breaks = c(1, 3, 10, 30, 100, 300, max(grid$median_controls_95)),
                         name = "Median controls\n(log scale)",
                         guide = guide_colourbar(theme = theme(legend.key.height = grid::unit(6, "lines")))) +
    labs(title = "Controls carrying 95% of a treated participant's transport weight") +
    axis_labels + style

distance <- abs(grid$estimate - benchmark)
closest <- grid[grid$is_selected_table1 == 1, ]
stopifnot(nrow(closest) == 1, abs(closest$estimate - benchmark) == min(distance))
corners <- subset(grid, rho %in% range(rho) & epsilon %in% range(epsilon))
# Place full values just inside the edge; both right corners have dark fills.
corners$just <- ifelse(corners$rho == min(grid$rho), 0, 1)
corners$offset <- ifelse(corners$just == 0, -0.43, 0.43)
corners$ink <- ifelse(corners$estimate < -0.6 * max(abs(grid$estimate)), "white", "black")
contrasts <- ggplot(grid, aes(rho_f, eps_f, fill = estimate)) +
    geom_tile() +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                         midpoint = 0, name = "UOT contrast ($)") +
    geom_tile(data = closest, fill = NA, color = "red", linewidth = 1.2) +
    geom_text(data = closest, aes(label = sprintf("\u2248%.1f", estimate)),
              color = "red", fontface = "bold", size = 3) +
    geom_text(data = corners, aes(label = sprintf("%.2f", estimate), hjust = just,
                                  colour = ink), nudge_x = corners$offset,
              fontface = "bold", size = 3) +
    scale_colour_identity() +
    labs(title = sprintf("UOT outcome contrast across (rho, eps)\nHighlighted closest to experimental ATT ($%.2f)", benchmark)) +
    axis_labels + style + theme(panel.grid = element_blank())

# Same layout proportions as the original manuscript PNGs, at print resolution.
save_heatmap <- function(plot, filename, width, height) {
    png(file.path(plot_dir, filename), width = width / 96, height = height / 96,
        units = "in", res = 300, bg = "white")
    print(plot)
    dev.off()
}
save_heatmap(matches, "UOT_matched_pairs_sensi.png", 950, 540)
save_heatmap(contrasts, "UOT_ATT_sensi.png", 1063, 579)

# =============================================================================
# AF3 candidate ranking — revised paper figure
#
# INPUT:
#   ranking_for_R.csv
#
# OUTPUT:
#   candidate_ranking_combined.pdf / .svg / .png
#   candidate_ranking_Rx.pdf       / .svg / .png
#   candidate_ranking_Sr35.pdf     / .svg / .png
#
# Visual order from left to right:
#   Candidate / Effector
#   -> Contact pairs > 0.3
#   -> Max contact probability
#   -> iLIS
#   -> ipSAE
#   -> Interface consistency
#   -> QC
#
# Notes:
# - The rank values are READ from ranking_for_R.csv; this script does not
#   recalculate the ranking.
# - Columns are ordered by visual / benchmarking usefulness, not necessarily
#   by the historical lexicographic rank-key order.
# - "effector_unfolded" in the CSV is displayed as "low effector pLDDT (<50)",
#   because low pLDDT does not by itself prove that a protein is unfolded.
#
# Install once if needed:
# install.packages(c("ggplot2","dplyr","tidyr","patchwork","scales","svglite"))
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# ---- appearance -------------------------------------------------------------
BASE_FONT <- "Arial"
BASE_SIZE <- 9
BAR_MAX   <- 80

ACCENT   <- "#2a78d6"
INK      <- "#0b0b0b"
INK2     <- "#52514e"
MUTED    <- "#898781"
GRID     <- "#e1e0d9"
SURFACE  <- "#fcfcfb"
CRITICAL <- "#d03b3b"
WARNING  <- "#fab219"
SERIOUS  <- "#ec835a"

# ---- read data --------------------------------------------------------------
d <- read.csv("ranking_for_R.csv", stringsAsFactors = FALSE)

required_cols <- c(
  "backbone","candidate","effector","vlrr_class","rank",
  "contact_pairs","contact_prob","iLIS","ipSAE","consistency",
  "is_control","effector_unfolded","clash_tier"
)

missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0) {
  stop(
    "ranking_for_R.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Handle either logical TRUE/FALSE or character True/False safely.
as_flag <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  tolower(trimws(as.character(x))) %in% c("true","t","1","yes")
}

d <- d |>
  mutate(
    vlrr_class = factor(
      vlrr_class,
      levels = c("Yes", "Partially", "No"),
      labels = c(
        "vLRR interface: Yes",
        "vLRR interface: Partially",
        "vLRR interface: No"
      )
    ),
    is_ctrl = as_flag(is_control),
    low_effector_plddt = as_flag(effector_unfolded),
    label = ifelse(
      is_ctrl,
      paste0(candidate, "  (positive control)"),
      candidate
    )
  )

# Preserve the supplied rank ordering independently inside each backbone.
d <- d |>
  group_by(backbone) |>
  arrange(rank, .by_group = TRUE) |>
  mutate(y = -row_number()) |>
  ungroup()

# ---- global metric scaling --------------------------------------------------
# Dot positions are normalized with GLOBAL maxima across Rx + Sr35.
# This means the same metric is directly comparable between panels A and B.
METRICS <- c("contact_prob", "iLIS", "ipSAE", "consistency")

MLAB <- c(
  contact_prob = "max contact\nprobability",
  iLIS         = "iLIS",
  ipSAE        = "ipSAE",
  consistency  = "interface\nconsistency"
)

metric_max <- d |>
  summarise(across(all_of(METRICS), ~ max(.x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "metric", values_to = "metric_max")

# ---- theme ------------------------------------------------------------------
theme_paper <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size, base_family = BASE_FONT) +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      axis.text.y      = element_blank(),
      axis.text.x      = element_text(
        colour = MUTED,
        size = base_size - 1
      ),
      strip.background = element_blank(),
      strip.text       = element_blank(),
      plot.title       = element_text(
        colour = MUTED,
        size = base_size - 1,
        hjust = 0,
        face = "plain",
        margin = margin(b = 3)
      ),
      plot.margin = margin(2, 3, 2, 3),
      legend.position = "none"
    )
}

facet_bands <- facet_grid(
  rows = vars(vlrr_class),
  scales = "free_y",
  space = "free_y"
)

# ---- band labels ------------------------------------------------------------
p_band <- function(dd) {
  ggplot(dd, aes(y = y)) +
    geom_blank() +
    facet_bands +
    theme_paper() +
    theme(
      strip.text.y = element_text(
        angle = 90,
        colour = INK2,
        size = BASE_SIZE - 1,
        face = "bold"
      ),
      axis.text.x = element_blank()
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0))
}

# ---- candidate + effector ---------------------------------------------------
p_lab <- function(dd) {
  ggplot(dd, aes(y = y)) +
    geom_text(
      aes(x = 0, label = label, colour = is_ctrl),
      hjust = 0,
      size = BASE_SIZE / .pt,
      fontface = "bold",
      family = BASE_FONT
    ) +
    geom_text(
      aes(x = 0.64, label = effector),
      hjust = 0,
      colour = INK2,
      size = (BASE_SIZE - 0.5) / .pt,
      family = BASE_FONT
    ) +
    scale_colour_manual(values = c(`FALSE` = INK, `TRUE` = CRITICAL)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    facet_bands +
    theme_paper() +
    theme(axis.text.x = element_blank()) +
    labs(title = "Candidate / effector")
}

# ---- primary evidence: contact-pair bar ------------------------------------
p_bar <- function(dd) {
  ggplot(dd, aes(y = y)) +
    geom_vline(
      xintercept = seq(0, BAR_MAX, 20),
      colour = GRID,
      linewidth = 0.25
    ) +
    geom_segment(
      aes(x = 0, xend = contact_pairs, yend = y),
      colour = ACCENT,
      linewidth = 2.1,
      lineend = "round"
    ) +
    geom_text(
      aes(x = pmin(contact_pairs + 2.0, BAR_MAX - 1),
          label = sprintf("%.0f", contact_pairs)),
      hjust = 0,
      colour = INK2,
      size = (BASE_SIZE - 1) / .pt,
      family = BASE_FONT
    ) +
    scale_x_continuous(
      limits = c(0, BAR_MAX),
      breaks = seq(0, BAR_MAX, 20),
      expand = expansion(mult = c(0, 0.03))
    ) +
    facet_bands +
    theme_paper() +
    labs(title = "AF3 contact pairs > 0.3\n(NLR × effector; median of 5 models)")
}

# ---- supporting metrics -----------------------------------------------------
p_dots <- function(dd) {
  long <- dd |>
    select(y, vlrr_class, all_of(METRICS)) |>
    pivot_longer(
      all_of(METRICS),
      names_to = "metric",
      values_to = "value"
    ) |>
    left_join(metric_max, by = "metric") |>
    mutate(
      frac = case_when(
        is.na(value) ~ NA_real_,
        metric_max <= 0 ~ 0,
        TRUE ~ value / metric_max
      ),
      metric = factor(metric, levels = METRICS, labels = MLAB[METRICS]),
      lab = ifelse(
        is.na(value),
        "NA",
        sprintf("%.2f", value)
      )
    )

  ggplot(long, aes(y = y)) +
    # Zero = explicit small dash.
    geom_point(
      data = ~ subset(.x, !is.na(value) & value == 0),
      aes(x = 0),
      shape = 45,
      colour = MUTED,
      size = 1.7
    ) +
    # Non-zero = point position reflects value relative to the global max
    # of that metric across both backbones.
    geom_point(
      data = ~ subset(.x, !is.na(value) & value > 0),
      aes(x = frac),
      colour = ACCENT,
      size = 2.1
    ) +
    # Missing = open circle.
    geom_point(
      data = ~ subset(.x, is.na(value)),
      aes(x = 0),
      shape = 1,
      colour = MUTED,
      size = 1.5
    ) +
    geom_text(
      aes(x = 1.18, label = lab),
      hjust = 1,
      colour = INK2,
      size = (BASE_SIZE - 1.2) / .pt,
      family = BASE_FONT
    ) +
    scale_x_continuous(
      limits = c(0, 1.20),
      breaks = c(0, 0.5, 1),
      labels = NULL,
      expand = c(0, 0)
    ) +
    facet_grid(
      rows = vars(vlrr_class),
      cols = vars(metric),
      scales = "free_y",
      space = "free_y"
    ) +
    theme_paper() +
    theme(
      strip.text.x = element_text(
        colour = MUTED,
        size = BASE_SIZE - 1,
        margin = margin(b = 3)
      ),
      axis.text.x = element_blank(),
      panel.spacing.x = unit(2.5, "mm")
    )
}

# ---- QC ---------------------------------------------------------------------
p_qc <- function(dd) {
  q <- dd |>
    transmute(
      y,
      vlrr_class,
      clash = ifelse(is.na(clash_tier), "none", clash_tier),
      low_plddt = low_effector_plddt
    )

  ggplot(q, aes(y = y)) +
    geom_point(
      data = ~ subset(.x, clash == "severe"),
      aes(x = 0.20),
      shape = 17,
      colour = CRITICAL,
      size = 1.9
    ) +
    geom_point(
      data = ~ subset(.x, clash == "repeated"),
      aes(x = 0.20),
      shape = 17,
      colour = WARNING,
      size = 1.9
    ) +
    geom_point(
      data = ~ subset(.x, low_plddt),
      aes(x = 0.58),
      shape = 16,
      colour = SERIOUS,
      size = 1.9
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    facet_bands +
    theme_paper() +
    theme(axis.text.x = element_blank()) +
    labs(title = "QC")
}

# ---- build one backbone panel ----------------------------------------------
panel <- function(bb, letter) {
  dd <- dplyr::filter(d, backbone == bb)

  if (nrow(dd) == 0) {
    stop("No rows found for backbone: ", bb)
  }

  p <- wrap_plots(
    list(
      p_band(dd),
      p_lab(dd),
      p_bar(dd),
      p_dots(dd),
      p_qc(dd)
    ),
    nrow = 1,
    widths = c(0.16, 1.35, 2.15, 3.15, 0.38)
  ) +
    plot_annotation(
      title = sprintf("%s   %s candidate ranking", letter, bb),
      theme = theme(
        plot.title = element_text(
          family = BASE_FONT,
          face = "bold",
          size = BASE_SIZE + 2,
          colour = INK
        )
      )
    )

  p
}

rx_fig   <- panel("Rx", "A")
sr35_fig <- panel("Sr35", "B")

# ---- caption ----------------------------------------------------------------
CAP <- paste(
  "Candidates are grouped first by manual vLRR localisation (Yes / Partially / No).",
  "Columns are ordered left-to-right by visual and benchmarking usefulness:",
  "AF3 NLR-effector contact-pair evidence, max contact probability, iLIS, ipSAE,",
  "and multi-model interface consistency. Values are medians of 5 AF3 models.",
  "The supplied rank column is preserved from ranking_for_R.csv and is not recalculated here.",
  "All displayed interaction metrics are NLR x effector; ATP and intra-chain blocks are excluded.",
  "Triangle = repeated/severe custom steric-clash flag; orange circle = low effector pLDDT (<50).",
  "Positive controls are shown in red as references.",
  sep = "\n"
)

add_caption <- function(p) {
  p + plot_annotation(
    caption = CAP,
    theme = theme(
      plot.caption = element_text(
        family = BASE_FONT,
        size = BASE_SIZE - 1.5,
        colour = INK2,
        hjust = 0,
        margin = margin(t = 8)
      )
    )
  )
}

rx_out   <- add_caption(rx_fig)
sr35_out <- add_caption(sr35_fig)

combined <- wrap_elements(rx_fig) / wrap_elements(sr35_fig) +
  plot_layout(heights = c(31, 31)) +
  plot_annotation(
    caption = CAP,
    theme = theme(
      plot.caption = element_text(
        family = BASE_FONT,
        size = BASE_SIZE - 1.5,
        colour = INK2,
        hjust = 0,
        margin = margin(t = 8)
      )
    )
  )

# ---- save -------------------------------------------------------------------
save_all <- function(filename_stub, plot_obj, width, height) {
  ggsave(
    paste0(filename_stub, ".pdf"),
    plot_obj,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = SURFACE
  )
  ggsave(
    paste0(filename_stub, ".svg"),
    plot_obj,
    width = width,
    height = height,
    bg = SURFACE
  )
  ggsave(
    paste0(filename_stub, ".png"),
    plot_obj,
    width = width,
    height = height,
    dpi = 600,
    bg = SURFACE
  )
}

save_all("candidate_ranking_Rx", rx_out, 10.5, 6.5)
save_all("candidate_ranking_Sr35", sr35_out, 10.5, 6.5)
save_all("candidate_ranking_combined", combined, 10.5, 12.5)

message("Done.")
message("Wrote:")
message("  candidate_ranking_Rx.pdf / .svg / .png")
message("  candidate_ranking_Sr35.pdf / .svg / .png")
message("  candidate_ranking_combined.pdf / .svg / .png")

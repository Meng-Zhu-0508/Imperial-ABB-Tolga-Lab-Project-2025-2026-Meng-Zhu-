# =============================================================================
#  Rx / Sr35 candidate screen - figure generation
#  -----------------------------------------------------------------------
#  1. Put this file and the spreadsheet in the same folder (or edit DATA_FILE).
#  2. Open in RStudio and press "Source"  (or run  source("make_figures.R") ).
#
#  Produces, in OUT_DIR:
#      Fig1a_Rx_Candidate_Screen.png / .pdf
#      Fig1b_Sr35_Candidate_Screen.png / .pdf
#      Fig2_individual/Fig2_<candidate>.png   (one per candidate)
#      Fig2_All_Candidates.pdf                (one candidate per page)
#      Fig3_Control_Panel.png / .pdf
#      Fig1_statistics.csv, Fig2_statistics.csv, Normalised_wells.csv
#
#  Readout: chemiluminescence from the FBP reporter -> LOWER value = STRONGER HR.
#  Normalisation: min-max within each infiltration date, per well, applied
#  BEFORE any averaging.  0 = lowest chemiluminescence that day (strongest HR),
#  1 = highest (weakest HR).
#
#  Packages:
#      install.packages(c("readxl","dplyr","tidyr","purrr","stringr","ggplot2"))
#      # optional, gives better Arial rendering in PNG:
#      install.packages("ragg")
# =============================================================================

## ------------------------------------------------------------------ CONFIG --
DATA_FILE <- "Raw_Data_Summary_Alex.xlsx"   # path to the spreadsheet
OUT_DIR   <- "figures"                      # where everything is written

BASE_FONT     <- "Arial"
DEDUP         <- TRUE   # drop control blocks that appear identically in two sheets
DROP_UNUSABLE <- TRUE   # drop every row noted 'unusable' - not used for anything
IRR_THRESHOLD      <- 0.40  # flag candidate light red when > 40% of wells are 'irregular'
IRR_THRESHOLD_CTRL <- 0.50  # Figure 3: a control condition pools several plates/dates, so
                            # it is flagged only when MOST of its wells are 'irregular'

# Figure-1 negative control: among the non-cognate conditions available for that
# candidate, take the BRIGHTEST one (highest mean normalised chemiluminescence =
# weakest HR).  PTO is the last resort if no non-cognate condition survives cleaning.
NONCOGNATE  <- c("PvxCP", "AvrSr35", "SRE12")
LAST_RESORT <- "PTO"

COL_POS <- "#9FB8DD"; COL_NEG <- "#C8C8C8"
COL_TEST <- "#79C9B4"; COL_IRR <- "#F3A8A8"; COL_BAND <- "#FDEDED"
CTRL_COLS <- c("Rx Controls" = "#B9D0E8", "Sr35 Controls" = "#F2C89B",
               "Effector Controls" = "#D3D3D3", "General Controls" = "#9A9A9A")
DATE_PAL <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
              "#D55E00", "#CC79A7", "#8C564B", "#333333")

FS_TITLE <- 22; FS_SUB <- 14; FS_AXIS <- 18; FS_TICK <- 16; FS_LEG <- 14

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(ggplot2)
})
set.seed(7)
dir.create(file.path(OUT_DIR, "Fig2_individual"), showWarnings = FALSE, recursive = TRUE)

png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else NULL

save_fig <- function(plot, name, width, height) {
  png_path <- file.path(OUT_DIR, paste0(name, ".png"))
  pdf_path <- file.path(OUT_DIR, paste0(name, ".pdf"))
  if (is.null(png_device)) {
    ggsave(png_path, plot, width = width, height = height, dpi = 600, bg = "white")
  } else {
    ggsave(png_path, plot, device = png_device, width = width, height = height,
           dpi = 600, bg = "white")
  }
  if (is.null(pdf_device)) {
    ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  } else {
    ggsave(pdf_path, plot, device = pdf_device, width = width, height = height,
           bg = "white")
  }
}

## -------------------------------------------------------------- LOAD DATA --
sheets <- excel_sheets(DATA_FILE)

raw <- map_dfr(sheets, function(s) {
  df <- suppressMessages(read_excel(DATA_FILE, sheet = s))
  names(df) <- tolower(str_trim(names(df)))
  if ("notes" %in% names(df)) df <- rename(df, note = notes)
  if (!"note"  %in% names(df)) df$note  <- NA
  if (!"plate" %in% names(df)) df$plate <- NA
  if (!"dpi"   %in% names(df)) df$dpi   <- NA
  df %>%
    filter(!is.na(conditions)) %>%
    transmute(sheet      = s,
              conditions = str_squish(as.character(conditions)),
              val        = suppressWarnings(as.numeric(as.character(chemillumience))),
              date       = as.character(`infiltration date`),
              dpi        = str_remove(as.character(dpi),   "\\.0$"),
              plate      = str_remove(as.character(plate), "\\.0$"),
              note       = str_trim(tolower(replace_na(as.character(note), ""))))
})

n_rows_in_file <- nrow(raw)
dat <- raw %>% filter(!is.na(val))

if (DROP_UNUSABLE) {
  gone <- dat %>% filter(str_detect(note, "unusable"))
  if (nrow(gone) > 0) {
    message("dropped ", nrow(gone), " rows noted \"unusable\":")
    gone %>% count(conditions) %>%
      purrr::pwalk(function(conditions, n)
        message("    ", formatC(conditions, width = -22), n, " well(s)"))
  }
  dat <- dat %>% filter(!str_detect(note, "unusable"))
}

dat <- dat %>%
  mutate(A = str_trim(str_split_fixed(conditions, "\\+", 2)[, 1]),
         B = str_trim(str_split_fixed(conditions, "\\+", 2)[, 2]))

# De-duplicate identical control blocks copied into more than one sheet.
# A "block" = all wells of one condition on one date / plate / dpi within one sheet.
# Repeated values INSIDE a block are genuine separate wells and are kept.
if (DEDUP) {
  blocks <- dat %>%
    group_by(sheet, conditions, date, plate, dpi) %>%
    summarise(sig = paste(sort(round(val, 6)), collapse = "|"), .groups = "drop") %>%
    group_by(conditions, date, plate, dpi, sig) %>%
    mutate(first_sheet = first(sheet)) %>%
    ungroup()
  gone <- blocks %>% filter(sheet != first_sheet)
  if (nrow(gone) > 0)
    for (i in seq_len(nrow(gone)))
      message("de-duplicated: ", gone$conditions[i], " ", gone$date[i], " p", gone$plate[i],
              " - identical block already in \"", gone$first_sheet[i],
              "\", second copy in \"", gone$sheet[i], "\"")
  dat <- semi_join(dat,
                   blocks %>% filter(sheet == first_sheet) %>%
                     select(sheet, conditions, date, plate, dpi),
                   by = c("sheet", "conditions", "date", "plate", "dpi"))
}

# min-max normalisation WITHIN each infiltration date
dat <- dat %>%
  group_by(date) %>%
  mutate(norm = (val - min(val)) / (max(val) - min(val))) %>%
  ungroup() %>%
  mutate(irregular = str_detect(note, "irregular"))

message(n_rows_in_file, " rows in file -> ", nrow(dat), " wells used")

date_levels <- unique(dat$date)[order(as.Date(unique(dat$date), format = "%d.%m.%Y"))]
dat$date    <- factor(dat$date, levels = date_levels)
date_cols   <- setNames(rep(DATE_PAL, length.out = length(date_levels)), date_levels)

## -------------------------------------------------------------- CANDIDATES --
cand_names <- sort(unique(dat$A[str_detect(dat$A, "^(Rx[0-9]+|Sr35-[0-9]+)$")]))

# usable = the condition still has at least one well after cleaning
is_usable <- function(cd, eff) any(dat$A == cd & dat$B == eff)
pick_neg <- function(cd, bb, assigned) {
  opts <- setdiff(NONCOGNATE, assigned)
  opts <- opts[vapply(opts, function(e) is_usable(cd, e), logical(1))]
  if (length(opts) > 0) {
    m <- vapply(opts, function(e) mean(dat$norm[dat$A == cd & dat$B == e]), numeric(1))
    m <- sort(m, decreasing = TRUE)
    return(list(eff  = names(m)[1],
                rule = "brightest non-cognate",
                opts = paste(sprintf("%s (n=%d, mean=%.3f)", names(m),
                                     vapply(names(m), function(e)
                                       sum(dat$A == cd & dat$B == e), integer(1)), m),
                             collapse = "; ")))
  }
  if (is_usable(cd, LAST_RESORT))
    return(list(eff = LAST_RESORT, rule = "no usable non-cognate -> PTO", opts = ""))
  list(eff = NA_character_, rule = "none available", opts = "")
}

cand_info <- map_dfr(cand_names, function(cd) {
  sub  <- dat %>% filter(A == cd)
  effs <- unique(sub$B)
  bb   <- if (str_starts(cd, "Rx")) "Rx" else "Sr35"
  ass  <- effs[str_starts(toupper(effs), "SRE") & effs != "SRE12"]
  neg  <- pick_neg(cd, bb, ass)
  tibble(candidate = cd, backbone = bb,
         num      = as.numeric(str_extract(cd, "[0-9]+$")),
         assigned = paste(ass, collapse = "/"),
         neg_ctrl = neg$eff,
         neg_rule = neg$rule,
         neg_opts = neg$opts,
         n_wells  = nrow(sub),
         pct_irr  = mean(sub$irregular),
         flag     = mean(sub$irregular) > IRR_THRESHOLD)
}) %>% arrange(backbone != "Rx", num)

assigned_of <- function(cd)
  str_split(cand_info$assigned[cand_info$candidate == cd], "/")[[1]]

wells_of <- function(cd, eff = NULL, dates = NULL) {
  s <- dat %>% filter(A == cd)
  if (!is.null(eff))   s <- s %>% filter(B %in% eff)
  if (!is.null(dates)) s <- s %>% filter(date %in% dates)
  s
}
pos_ctrl_of <- function(bb, dates = NULL) {
  s <- if (bb == "Rx") {
         dat %>% filter(A == "RxWT", B == "PvxCP")
       } else {
         dat %>% filter(A %in% c("Sr35WT", "Sr35"), B == "AvrSr35")
       }
  if (!is.null(dates) && nrow(filter(s, date %in% dates)) > 0)
    s <- s %>% filter(date %in% dates)
  s
}
star_of <- function(p)
  ifelse(is.na(p), "n/a",
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

role_cols <- c(pos = COL_POS, neg = COL_NEG, test = COL_TEST, test_flag = COL_IRR)
role_labs <- c(pos = "Positive control", neg = "Negative control",
               test = "Candidate test",
               test_flag = "Irregular cells")

# white diamond at the mean, drawn on top of every box
PT_SIZE <- 2.2                 # size of a single well marker
MEAN_SIZE <- PT_SIZE * 1.15    # mean diamond - same visual size as the wells

mean_layer <- function(sz = MEAN_SIZE)
  stat_summary(aes(shape = "Mean"), fun = mean, geom = "point", size = sz,
               fill = "white", colour = "#222222", stroke = 0.5)
mean_scale <- function(order = 3)
  scale_shape_manual(name = NULL, values = c(Mean = 23),
                     guide = guide_legend(order = order,
                                          keywidth = grid::unit(0.9, "lines")))

base_theme <- theme_classic(base_size = FS_TICK, base_family = BASE_FONT) +
  theme(plot.title    = element_text(size = FS_TITLE, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = FS_SUB, colour = "#555555"),
        axis.title    = element_text(size = FS_AXIS, face = "bold"),
        axis.text     = element_text(size = FS_TICK, colour = "black"),
        axis.title.x  = element_text(margin = margin(t = 9)),
        axis.title.y  = element_text(margin = margin(r = 9)),
        legend.title  = element_text(size = FS_LEG),
        legend.text   = element_text(size = FS_LEG),
        legend.key    = element_blank(),
        legend.position      = "right",
        legend.justification = "top",      # pinned to the top right
        legend.box.just      = "left",     # both legends share the same left edge
        legend.margin        = margin(l = 0),
        panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.4))

## ---------------------------------------------------------------- FIGURE 1 --
# The test uses every well of both groups.  BH correction is applied within each panel.
fig1_stats <- map_dfr(cand_info$candidate, function(cd) {
  i <- cand_info %>% filter(candidate == cd)
  a <- wells_of(cd, assigned_of(cd))
  n <- if (is.na(i$neg_ctrl)) a[0, ] else wells_of(cd, i$neg_ctrl)
  p <- if (nrow(a) > 1 && nrow(n) > 1)
         t.test(a$norm, n$norm, var.equal = FALSE)$p.value else NA_real_
  tibble(candidate = cd, welch_p = p,
         n_assigned = nrow(a), mean_assigned = mean(a$norm),
         dates_assigned = paste(sort(unique(as.character(a$date))), collapse = "; "),
         n_neg = nrow(n), mean_neg = mean(n$norm),
         dates_neg = paste(sort(unique(as.character(n$date))), collapse = "; "))
})
fig1_stats <- cand_info %>% left_join(fig1_stats, by = "candidate") %>%
  group_by(backbone) %>%                                # BH within each panel
  mutate(BH_adj_p = p.adjust(welch_p, method = "BH")) %>%
  ungroup() %>%
  mutate(significance = star_of(BH_adj_p))
write.csv(fig1_stats, file.path(OUT_DIR, "Fig1_statistics.csv"), row.names = FALSE)

make_fig1 <- function(bb, title, pos_label, file_name, width) {
  info <- fig1_stats %>% filter(backbone == bb)

  # explicit numeric x positions: positive control at 0, candidate k at 1.15*k
  pos_pts <- pos_ctrl_of(bb) %>% mutate(x = 0, role = "pos")
  cand_pts <- map_dfr(seq_len(nrow(info)), function(k) {
    i  <- info[k, ]; xc <- 1.15 * k
    bind_rows(
      (if (is.na(i$neg_ctrl)) dat[0, ] else wells_of(i$candidate, i$neg_ctrl)) %>%
        mutate(x = xc - 0.23, role = "neg"),
      wells_of(i$candidate, assigned_of(i$candidate)) %>%
        mutate(x = xc + 0.23, role = if (i$flag) "test_flag" else "test"))
  })
  pts <- bind_rows(pos_pts, cand_pts) %>%
    mutate(role = factor(role, levels = names(role_cols)))

  ticks  <- c(0, 1.15 * seq_len(nrow(info)))
  labels <- c(pos_label, as.character(info$num))
  sig    <- info %>% mutate(x = 1.15 * row_number())
  bands  <- sig %>% filter(flag)

  p <- ggplot()
  if (nrow(bands) > 0)
    p <- p + annotate("rect", xmin = bands$x - 0.56, xmax = bands$x + 0.56,
                      ymin = -Inf, ymax = Inf, fill = COL_BAND)
  p <- p +
    geom_boxplot(data = pts, aes(x = x, y = norm, group = x, fill = role),
                 width = 0.34, outlier.shape = NA, colour = "#333333", linewidth = 0.4) +
    geom_point(data = pts, aes(x = x, y = norm, colour = date),
               size = PT_SIZE, position = position_jitter(width = 0.075, height = 0)) +
    stat_summary(data = pts, aes(x = x, y = norm, group = x, shape = "Mean"), fun = mean,
                 geom = "point", size = MEAN_SIZE, fill = "white", colour = "#222222",
                 stroke = 0.5) +
    mean_scale(3) +
    geom_text(data = sig, aes(x = x, y = 1.07, label = significance,
                              fontface = ifelse(!is.na(BH_adj_p) & BH_adj_p < 0.05,
                                                "bold", "plain")),
              size = FS_TICK / .pt, family = BASE_FONT) +
    scale_colour_manual(name = "Infiltration date", values = date_cols,
                        breaks = date_levels, drop = TRUE,
                        guide = guide_legend(order = 1,
                                             keywidth = grid::unit(0.9, "lines"))) +
    scale_fill_manual(name = "Conditions", values = role_cols, labels = role_labs,
                      breaks = names(role_cols), drop = FALSE,
                      guide = guide_legend(order = 2)) +     # conditions underneath
    scale_x_continuous(breaks = ticks, labels = labels, expand = expansion(0, 0),
                       limits = c(-0.60, 1.15 * nrow(info) + 0.60)) +
    scale_y_continuous(limits = c(-0.04, 1.16), breaks = seq(0, 1, 0.25)) +
    labs(title = title,
         subtitle = paste0("Welch's unpaired t-test with Benjamini-Hochberg FDR ",
                           "correction:  * p<0.05   ** p<0.01   *** p<0.001"),
         x = "Candidate", y = "Normalised Chemiluminescence (0-1)") +
    base_theme
  save_fig(p, file_name, width, 6.2)
  p
}

make_fig1("Rx",   "1a. Rx Candidate Screen",   "Rx + PVX-CP",
          "Fig1a_Rx_Candidate_Screen",   16.5)
make_fig1("Sr35", "1b. Sr35 Candidate Screen", "Sr35 + AvrSr35",
          "Fig1b_Sr35_Candidate_Screen", 11.5)

## ---------------------------------------------------------------- FIGURE 2 --
# Descriptive phenotype profile only - no ANOVA / Dunnett.  The significance quoted in
# the subtitle is the Figure 1 comparison (assigned effector vs negative control).
WT <- c(Rx = "PvxCP", Sr35 = "AvrSr35")
fig2_stats <- list()

make_fig2 <- function(cd) {
  i     <- fig1_stats %>% filter(candidate == cd)
  bb    <- i$backbone
  ass   <- assigned_of(cd)
  A     <- paste(ass, collapse = "/")
  cdate <- unique(wells_of(cd)$date)

  parts <- list(
    list(lab = if (bb == "Rx") "Rx + PVX-CP" else "Sr35 + AvrSr35",
         d = pos_ctrl_of(bb, cdate), role = "pos"),
    list(lab = A, d = wells_of(cd, ass),
         role = if (i$flag) "test_flag" else "test"),
    list(lab = "PTO",                d = wells_of(cd, "PTO"),     role = "neg"),
    list(lab = WT[[bb]],             d = wells_of(cd, WT[[bb]]),  role = "neg"),
    list(lab = "SRE12",              d = wells_of(cd, "SRE12"),   role = "neg"))
  parts <- Filter(function(p) nrow(p$d) > 0, parts)

  lv  <- map_chr(parts, "lab")
  pts <- map_dfr(parts, function(p) p$d %>% mutate(cond = p$lab, role = p$role)) %>%
    mutate(cond = factor(cond, levels = lv),
           role = factor(role, levels = names(role_cols)))

  fig2_stats[[cd]] <<- pts %>% group_by(cond) %>%
    summarise(n = n(), mean_norm = mean(norm), median_norm = median(norm),
              sd_norm = ifelse(n() > 1, sd(norm), NA_real_),
              dates = paste(sort(unique(as.character(date))), collapse = "; "),
              .groups = "drop") %>%
    mutate(candidate = cd, backbone = bb,
           role = ifelse(cond == A, "assigned effector",
                  ifelse(cond == i$neg_ctrl, "negative control", "")))

  sub1 <- if (!is.na(i$neg_ctrl)) {
            paste0(i$significance, " against negative control (", cd, " + ", i$neg_ctrl, ")")
          } else {
            "no negative control available"
          }

  p <- ggplot(pts, aes(x = cond, y = norm)) +
    geom_boxplot(aes(fill = role), width = 0.55, outlier.shape = NA,
                 colour = "#333333", linewidth = 0.4, show.legend = FALSE) +
    geom_point(aes(colour = date), size = PT_SIZE,
               position = position_jitter(width = 0.13, height = 0)) +
    mean_layer() + mean_scale(2) +
    scale_fill_manual(values = role_cols, drop = FALSE) +
    scale_colour_manual(name = "Infiltration date", values = date_cols,
                        breaks = date_levels, drop = TRUE,
                        guide = guide_legend(order = 1,
                                             keywidth = grid::unit(0.9, "lines"))) +
    scale_y_continuous(limits = c(-0.04, 1.16), breaks = seq(0, 1, 0.25)) +
    labs(title = paste0("2. ", cd, "  (assigned effector: ", A, ")"),
         subtitle = paste0(sub1, "\n",
                           "Welch's unpaired t-test, Benjamini-Hochberg adjusted (see Figure 1)"),
         x = "Condition", y = "Normalised Chemiluminescence (0-1)") +
    base_theme +
    theme(axis.text.x  = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 9)),
          plot.title    = element_text(size = FS_TITLE, face = "bold",
                                       margin = margin(b = 14)),
          plot.subtitle = element_text(size = FS_SUB, colour = "#555555",
                                       margin = margin(b = 10)))

  save_fig(p, file.path("Fig2_individual", paste0("Fig2_", cd)), 8.6, 5.4)
  p
}

fig2_plots <- map(cand_info$candidate, make_fig2)
write.csv(bind_rows(fig2_stats), file.path(OUT_DIR, "Fig2_statistics.csv"), row.names = FALSE)

if (is.null(pdf_device)) {
  pdf(file.path(OUT_DIR, "Fig2_All_Candidates.pdf"), width = 8.6, height = 5.4, onefile = TRUE)
} else {
  pdf_device(file.path(OUT_DIR, "Fig2_All_Candidates.pdf"), width = 8.6, height = 5.4,
             onefile = TRUE)
}
for (p in fig2_plots) print(p)
invisible(dev.off())

## ---------------------------------------------------------------- FIGURE 3 --
ctrl_dat <- dat %>%
  filter(sheet %in% c("Controls", "effectors Alone", "Positive Controls for candidate")) %>%
  mutate(grp = case_when(
    A == "RxWT" | A == "PvxCP" | B == "PvxCP"              ~ "Rx Controls",
    A %in% c("Sr35WT", "Sr35", "AvrSr35") | B == "AvrSr35" ~ "Sr35 Controls",
    str_starts(toupper(A), "SRE")                          ~ "Effector Controls",
    TRUE                                                   ~ "General Controls")) %>%
  mutate(grp = factor(grp, levels = names(CTRL_COLS))) %>%
  arrange(grp, conditions) %>%
  mutate(conditions = factor(conditions, levels = unique(conditions)))

# light red background behind control conditions that are mostly 'irregular'
ctrl_bands <- ctrl_dat %>%
  group_by(conditions) %>%
  summarise(frac = mean(irregular), .groups = "drop") %>%
  mutate(idx = as.integer(conditions)) %>%
  filter(frac > IRR_THRESHOLD_CTRL)

p3 <- ggplot(ctrl_dat, aes(x = conditions, y = norm)) +
  geom_rect(data = ctrl_bands, inherit.aes = FALSE,
            aes(xmin = idx - 0.45, xmax = idx + 0.45, ymin = -Inf, ymax = Inf),
            fill = COL_BAND) +
  geom_boxplot(aes(fill = grp), width = 0.6, outlier.shape = NA,
               colour = "#333333", linewidth = 0.4) +
  geom_point(aes(colour = date), size = PT_SIZE,
             position = position_jitter(width = 0.15, height = 0)) +
  mean_layer() + mean_scale(3) +
  scale_fill_manual(name = "Conditions", values = CTRL_COLS,
                    guide = guide_legend(order = 2)) +
  scale_colour_manual(name = "Infiltration date", values = date_cols,
                      breaks = date_levels, drop = TRUE,
                      guide = guide_legend(order = 1,
                                           keywidth = grid::unit(0.9, "lines"))) +
  scale_y_continuous(limits = c(-0.04, 1.16), breaks = seq(0, 1, 0.25)) +
  labs(title = "3. Control Panel", x = "Control Conditions",
       y = "Normalised Chemiluminescence (0-1)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = FS_TICK - 2))
save_fig(p3, "Fig3_Control_Panel", 15, 6.0)

write.csv(dat, file.path(OUT_DIR, "Normalised_wells.csv"), row.names = FALSE)
message("Done - written to '", normalizePath(OUT_DIR), "'")

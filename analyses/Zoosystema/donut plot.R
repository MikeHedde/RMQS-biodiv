library(tidyverse)
library(ggforce)

# ---- helper: safe NA -> explicit labels
na_label <- function(x, lab) {
  x <- as.character(x)
  x <- ifelse(is.na(x) | x == "", lab, x)
  x
}

# ---- build one ring (one taxonomic level) within each method
build_ring <- function(df, group_cols, inner_r, outer_r, label_min_prop = 0.08) {
  # df must have: method, value, and columns in group_cols
  # returns per-method arcs with start/end angles + optional labels
  
  d <- df %>%
    group_by(across(all_of(c("method", group_cols)))) %>%
    summarise(value = sum(value), .groups = "drop")
  
  # within each method: compute cumulative angles
  d2 <- d %>%
    group_by(method) %>%
    mutate(
      total = sum(value),
      prop  = value / total,
      end   = 2*pi*cumsum(prop),
      start = lag(end, default = 0),
      r0 = inner_r,
      r1 = outer_r,
      # label only if segment "large enough" (tune label_min_prop)
      label = if_else(prop >= label_min_prop, as.character(.data[[tail(group_cols, 1)]]), NA_character_),
      angle = (start + end)/2
    ) %>%
    ungroup()
  
  d2
}

# =========================
# 1) Prepare data
# =========================
df0 <- data_clean %>%
  transmute(
    method,
    class  = na_label(class,  "Unidentified class"),
    order  = na_label(order,  "Unidentified order"),
    family = na_label(family, "Unidentified family"),
    species = na_label(valid_name, "Unidentified species"),
    value  = if ("tot_abund" %in% names(data_clean)) tot_abund else 1
  ) %>%
  mutate(value = replace_na(value, 1))

# =========================
# 2) Build 4 rings (class/order/family/species)
#    ring radii: adjust to taste
# =========================
ring1 <- build_ring(df0, c("class"),                  inner_r = 0.00, outer_r = 0.28, label_min_prop = 0.12)
ring2 <- build_ring(df0, c("class","order"),          inner_r = 0.30, outer_r = 0.52, label_min_prop = 0.10)
ring3 <- build_ring(df0, c("class","order","family"), inner_r = 0.54, outer_r = 0.76, label_min_prop = 0.08)

# espèces : en statique, labels souvent illisibles => on met label_min_prop plus haut ou zéro
ring4 <- build_ring(df0, c("class","order","family","species"),
                    inner_r = 0.78, outer_r = 1.00, label_min_prop = 0.00)  # pas de labels espèces

# Combine rings
arcs <- bind_rows(
  ring1 %>% mutate(level = "class"),
  ring2 %>% mutate(level = "order"),
  ring3 %>% mutate(level = "family"),
  ring4 %>% mutate(level = "species")
)

# =========================
# 3) Plot: 5 panels (facet by method)
# =========================
p <- ggplot(arcs) +
  ggforce::geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = r0, r = r1,
      start = start, end = end,
      fill = class
    ),
    colour = "white", linewidth = 0.15
  ) +
  coord_fixed() +
  facet_wrap(~ method, ncol = 3) +
  #scale_fill_brewer(palette = "Set2") +
  guides(fill = guide_legend(title = NULL)) +
  theme_void(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Taxonomic composition by sampling method")

p
ggsave("figure/Zootaxonomy/donut_taxo_5panels.png")

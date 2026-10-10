## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
# MISTy
library(mistyR)
library(future)

# data manipulation
library(dplyr)
library(purrr)
library(distances)

# plotting
library(ggplot2)

plan(multisession)

## ----fig.height=2.5, fig.width=3.5--------------------------------------------
data("synthetic")

ggplot(synthetic[[1]], aes(x = col, y = row, color = type)) +
  geom_point(shape = 15, size = 0.7) +
  scale_color_manual(values = c("#e9eed3", "#dcc38d", "#c9e2ad", "#a6bab6")) +
  theme_void()

str(synthetic[[1]], give.attr = FALSE)

## -----------------------------------------------------------------------------
expr <- synthetic[[1]] %>% select(-c(row, col, type, starts_with("lig")))

misty.intra <- create_initial_view(expr)

summary(misty.intra)
summary(misty.intra$intraview)

## -----------------------------------------------------------------------------
pos <- synthetic[[1]] %>% select(row, col)

misty.views <- misty.intra %>% add_paraview(pos, l = 10)

summary(misty.views)

## -----------------------------------------------------------------------------
# find the 10 nearest neighbors
neighbors <- nearest_neighbor_search(distances(as.matrix(pos)), k = 11)[-1, ]

# calculate the mean expression of the nearest neighbors for all markers
# for each cell in expr
nnexpr <- seq_len(nrow(expr)) %>%
  map_dfr(~ expr %>%
    slice(neighbors[, .x]) %>%
    colMeans())

nn.view <- create_view("nearest", nnexpr, "nn")

nn.view

## -----------------------------------------------------------------------------
extended.views <- misty.views %>% add_views(nn.view)

summary(extended.views)

## -----------------------------------------------------------------------------
extended.views %>%
  remove_views("nearest") %>%
  summary()

extended.views %>%
  remove_views("intraview") %>%
  summary()

## -----------------------------------------------------------------------------
misty.views %>% run_misty()

## -----------------------------------------------------------------------------
result.folders <- synthetic %>% imap_chr(function(sample, name) {
  sample.expr <- sample %>% select(-c(row, col, type, starts_with("lig")))
  sample.pos <- sample %>% select(row, col)
  
  create_initial_view(sample.expr) %>% add_paraview(sample.pos, l = 10) %>%
    run_misty(results.folder = paste0("results", .Platform$file.sep, name))
})

result.folders

## -----------------------------------------------------------------------------
misty.results <- collect_results(result.folders)

summary(misty.results)

## ----warning=FALSE------------------------------------------------------------
misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

## -----------------------------------------------------------------------------
misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>% 
  summarize(mean.p = mean(value)) %>%
  arrange(mean.p)

## -----------------------------------------------------------------------------
misty.results %>% plot_view_contributions()

## -----------------------------------------------------------------------------
misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.8)

## -----------------------------------------------------------------------------
misty.results %>% plot_interaction_heatmap(view = "para.10", cutoff = 0.5)

## -----------------------------------------------------------------------------
misty.results %>% plot_contrast_heatmap("intra", "para.10", cutoff = 0.5)

## -----------------------------------------------------------------------------
misty.results %>% plot_interaction_communities("intra")

## -----------------------------------------------------------------------------
misty.results %>% plot_interaction_communities("para.10", cutoff = 0.5)

## ----info, echo=FALSE---------------------------------------------------------
sessionInfo()

## ----cleanup, include=FALSE---------------------------------------------------
unlink("results", recursive = TRUE)


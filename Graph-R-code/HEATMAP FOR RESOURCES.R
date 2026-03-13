# This is for generating the heatmap
# Raw input files
files <- c(
  "CLEAN_resource_usage.tsv",
  "GUT_resource_usage.tsv",
  "MARINE_resource_usage.tsv"
)

# Convert elapsed_wall mm:ss.ss to seconds
to_seconds <- function(x) {
  sapply(strsplit(x, ":"), function(t) {
    as.numeric(t[1]) * 60 + as.numeric(t[2])
  })
}

# Function to normalize one dataset (THREAD 1 ONLY)
normalize_dataset <- function(file) {

  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  data$elapsed_wall_sec <- to_seconds(data$elapsed_wall)

  # KEEP ONLY THREAD 1
  data <- data[data$threads == 1, ]

  # Keep tool order fixed
  data$tool <- factor(
    data$tool,
    levels = c("fastp", "cutadapt", "trimmomatic", "skewer")
  )

  # Add sample name
  sample_name <- sub("_resource_usage\\.tsv$", "", basename(file))
  data$sample <- sample_name

  # Assign replicate number within each tool by row order
  data$replicate <- ave(
    seq_len(nrow(data)),
    data$tool,
    FUN = seq_along
  )

  # Convert to long format
  plot_data <- rbind(
    data.frame(
      sample = data$sample,
      replicate = data$replicate,
      tool = data$tool,
      resource = "cpu_total_seconds",
      value = data$cpu_total_seconds
    ),
    data.frame(
      sample = data$sample,
      replicate = data$replicate,
      tool = data$tool,
      resource = "max_rss_kb",
      value = data$max_rss_kb
    ),
    data.frame(
      sample = data$sample,
      replicate = data$replicate,
      tool = data$tool,
      resource = "elapsed_wall",
      value = data$elapsed_wall_sec
    )
  )

  plot_data$resource <- factor(
    plot_data$resource,
    levels = c("cpu_total_seconds", "max_rss_kb", "elapsed_wall")
  )

  # Normalize WITHIN each sample x replicate x resource
  # Reference = maximum value among tools in that block
  ref_table <- aggregate(
    value ~ sample + replicate + resource,
    data = plot_data,
    FUN = max
  )
  names(ref_table)[names(ref_table) == "value"] <- "ref_value"

  plot_data <- merge(
    plot_data,
    ref_table,
    by = c("sample", "replicate", "resource")
  )

  plot_data$normalized_value <- plot_data$value / plot_data$ref_value

  plot_data
}

# Combine all datasets
normalized_list <- lapply(files, normalize_dataset)
combined_normalized <- do.call(rbind, normalized_list)

# Factor ordering
combined_normalized$tool <- factor(
  combined_normalized$tool,
  levels = c("fastp", "cutadapt", "trimmomatic", "skewer")
)

combined_normalized$resource <- factor(
  combined_normalized$resource,
  levels = c("cpu_total_seconds", "max_rss_kb", "elapsed_wall")
)

combined_normalized$sample <- factor(
  combined_normalized$sample,
  levels = c("CLEAN", "GUT", "MARINE")
)

# Paired block ID = sample x replicate
combined_normalized$block <- interaction(
  combined_normalized$sample,
  combined_normalized$replicate,
  drop = TRUE
)

print(combined_normalized)

# PAIRED WILCOXON TESTS

tool_names <- levels(combined_normalized$tool)
resources  <- levels(combined_normalized$resource)

wilcoxon_results <- list()
k <- 1

for (res in resources) {

  subdat <- subset(combined_normalized, resource == res)

  # one row per paired block, one column per tool
  wide <- reshape(
    subdat[, c("block", "tool", "normalized_value")],
    idvar = "block",
    timevar = "tool",
    direction = "wide"
  )

  for (i in 1:(length(tool_names) - 1)) {
    for (j in (i + 1):length(tool_names)) {

      t1 <- tool_names[i]
      t2 <- tool_names[j]

      x <- wide[[paste0("normalized_value.", t1)]]
      y <- wide[[paste0("normalized_value.", t2)]]

      keep <- complete.cases(x, y)
      x <- x[keep]
      y <- y[keep]

      if (length(x) > 0) {
        wt <- suppressWarnings(
          wilcox.test(x, y, paired = TRUE, exact = FALSE)
        )

        med1 <- median(x, na.rm = TRUE)
        med2 <- median(y, na.rm = TRUE)

        winner <- if (med1 < med2) {
          t1
        } else if (med2 < med1) {
          t2
        } else {
          "tie"
        }

        wilcoxon_results[[k]] <- data.frame(
          resource = res,
          tool_1 = t1,
          tool_2 = t2,
          comparison = paste(t1, "vs", t2),
          n_blocks = length(x),
          statistic = unname(wt$statistic),
          p_value = wt$p.value,
          median_tool_1 = med1,
          median_tool_2 = med2,
          winner = winner,
          stringsAsFactors = FALSE
        )

        k <- k + 1
      }
    }
  }
}

wilcoxon_table <- do.call(rbind, wilcoxon_results)

# BH correction WITHIN each resource
wilcoxon_table$p_adjusted_BH <- NA_real_

for (res in resources) {
  idx <- which(wilcoxon_table$resource == res)
  wilcoxon_table$p_adjusted_BH[idx] <- p.adjust(wilcoxon_table$p_value[idx], method = "BH")
}

print(wilcoxon_table)

# SIGNIFICANCE HEATMAP

library(ggplot2)

resource_labels <- c(
  cpu_total_seconds = "CPU total use",
  max_rss_kb = "Max RAM use",
  elapsed_wall = "Total time"
)

wilcoxon_table$resource_label <- factor(
  resource_labels[wilcoxon_table$resource],
  levels = c("CPU total use", "Max RAM use", "Total time")
)

wilcoxon_table$comparison <- factor(
  wilcoxon_table$comparison,
  levels = c(
    "fastp vs cutadapt",
    "fastp vs trimmomatic",
    "fastp vs skewer",
    "cutadapt vs trimmomatic",
    "cutadapt vs skewer",
    "trimmomatic vs skewer"
  )
)

# Significance levels
wilcoxon_table$sig_level <- ifelse(
  wilcoxon_table$p_adjusted_BH < 0.001, "***",
  ifelse(
    wilcoxon_table$p_adjusted_BH < 0.01, "**",
    ifelse(
      wilcoxon_table$p_adjusted_BH < 0.05, "*",
      "ns"
    )
  )
)

wilcoxon_table$sig_level <- factor(
  wilcoxon_table$sig_level,
  levels = c("ns", "*", "**", "***")
)

# Text inside each tile
wilcoxon_table$cell_label <- ifelse(
  wilcoxon_table$winner == "tie",
  paste0("tie\n", wilcoxon_table$sig_level),
  paste0(wilcoxon_table$winner, "\n", wilcoxon_table$sig_level)
)

p_sig <- ggplot(
  wilcoxon_table,
  aes(x = comparison, y = resource_label, fill = sig_level)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = cell_label), size = 4.5) +
  scale_fill_manual(
    values = c(
      "ns"  = "white",
      "*"   = "#fcae91",
      "**"  = "#fb6a4a",
      "***" = "#cb181d"
    ),
    name = "Significance"
  ) +
  labs(
    x = "Tool comparison",
    y = "Resource"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

print(p_sig)
# Title - Pairwise Paired Wilcoxon Significance Heatmap"
# subtitle = "Tile color indicates significance level; text shows lower-median tool"
ggsave(
  "resource_heatmap.png",
  plot = p_sig,
  width = 10,
  height = 6,
  dpi = 300
)


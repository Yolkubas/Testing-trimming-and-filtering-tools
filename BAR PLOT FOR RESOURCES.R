# This is for generating the relative value bar plot

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

  data$tool <- factor(
    data$tool,
    levels = c("fastp", "cutadapt", "trimmomatic", "skewer")
  )

  # Convert to long format
  plot_data <- rbind(
    data.frame(tool = data$tool, resource = "cpu_percent",       value = data$cpu_percent),
    data.frame(tool = data$tool, resource = "cpu_total_seconds", value = data$cpu_total_seconds),
    data.frame(tool = data$tool, resource = "max_rss_kb",        value = data$max_rss_kb),
    data.frame(tool = data$tool, resource = "disk_delta_bytes",  value = data$disk_delta_bytes),
    data.frame(tool = data$tool, resource = "elapsed_wall",      value = data$elapsed_wall_sec)
  )

  plot_data$resource <- factor(
    plot_data$resource,
    levels = c(
      "cpu_percent",
      "cpu_total_seconds",
      "max_rss_kb",
      "disk_delta_bytes",
      "elapsed_wall"
    )
  )

  # Compute mean per tool
  tool_means <- aggregate(value ~ resource + tool, data = plot_data, FUN = mean)

  # Reference = maximum mean per resource
  resource_ref <- aggregate(value ~ resource, data = tool_means, FUN = max)
  names(resource_ref)[2] <- "ref_value"

  # Merge reference
  plot_data <- merge(plot_data, resource_ref, by = "resource")

  # Normalize
  plot_data$normalized_value <- plot_data$value / plot_data$ref_value

  # Mean normalized value per tool
  plot_mean <- aggregate(normalized_value ~ resource + tool, data = plot_data, FUN = mean)
  names(plot_mean)[3] <- "mean_norm"

  # Add sample name
  sample_name <- sub("_resource_usage\\.tsv$", "", basename(file))
  plot_mean$sample <- sample_name

  plot_mean
}

# Normalize all datasets
normalized_list <- lapply(files, normalize_dataset)

combined_normalized <- do.call(rbind, normalized_list)

# Factor ordering
combined_normalized$tool <- factor(
  combined_normalized$tool,
  levels = c("fastp", "cutadapt", "trimmomatic", "skewer")
)

combined_normalized$resource <- factor(
  combined_normalized$resource,
  levels = c(
    "cpu_percent",
    "cpu_total_seconds",
    "max_rss_kb",
    "disk_delta_bytes",
    "elapsed_wall"
  )
)

combined_normalized$sample <- factor(
  combined_normalized$sample,
  levels = c(
    "CLEAN",
    "GUT",
    "MARINE"
  )
)

# VISUALIZE

library(ggplot2)

# Keep only selected resources
plot_subset <- subset(
  combined_normalized,
  resource %in% c("cpu_total_seconds", "max_rss_kb", "elapsed_wall")
)

# Rename resources
plot_subset$resource <- factor(
  plot_subset$resource,
  levels = c("cpu_total_seconds", "max_rss_kb", "elapsed_wall"),
  labels = c("CPU total use", "Max RAM use", "Total time")
)

# Summarize across samples
plot_summary <- aggregate(mean_norm ~ resource + tool, data = plot_subset, FUN = mean)
plot_sd      <- aggregate(mean_norm ~ resource + tool, data = plot_subset, FUN = sd)

names(plot_summary)[3] <- "mean_value"
names(plot_sd)[3]      <- "sd_value"

plot_summary <- merge(plot_summary, plot_sd, by = c("resource", "tool"))

# Keep desired order
plot_summary$tool <- factor(
  plot_summary$tool,
  levels = c("fastp", "cutadapt", "trimmomatic", "skewer")
)

plot_summary$resource <- factor(
  plot_summary$resource,
  levels = c("CPU total use", "Max RAM use", "Total time")
)

############################ BAR PLOT RELATIVE ############################

p =ggplot(plot_summary, aes(x = resource, y = mean_value, fill = tool)) +
  geom_col(
    position = position_dodge(width = 0.9),
    width = 0.75
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(mean_value - sd_value, 0),
      ymax = mean_value + sd_value
    ),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  labs(
    x = "Resource",
    y = "Relative value",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(p)

ggsave(
  "resource_use_realtive.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)



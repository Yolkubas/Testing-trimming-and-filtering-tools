# This is for generating the raw value bar plot
# Read file
data <- read.table("CLEAN_resource_usage.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convert wall time to seconds
to_seconds <- function(x) {
  sapply(strsplit(x, ":"), function(t) {
    as.numeric(t[1]) * 60 + as.numeric(t[2])
  })
}

data$elapsed_wall_sec <- to_seconds(data$elapsed_wall)

# Keep thread 1 only
data <- data[data$threads == 1, ]

# Average replicates per tool
summary_table <- aggregate(
  cbind(cpu_total_seconds, max_rss_kb, elapsed_wall_sec) ~ tool,
  data = data,
  FUN = mean
)

print(summary_table)

library(ggplot2)

# Convert to long format
plot_data <- rbind(
  data.frame(tool = summary_table$tool,
             resource = "CPU total usage time (s)",
             value = summary_table$cpu_total_seconds),

  data.frame(tool = summary_table$tool,
             resource = "Max RAM use (kb)",
             value = summary_table$max_rss_kb),

  data.frame(tool = summary_table$tool,
             resource = "Total time (s)",
             value = summary_table$elapsed_wall_sec)
)

plot_data$tool <- factor(
  plot_data$tool,
  levels = c("fastp","cutadapt","trimmomatic","skewer")
)

# Plot
p <- ggplot(plot_data, aes(x = tool, y = value, fill = tool)) +
  geom_col(width = 0.7) +
  facet_wrap(~resource, scales = "free_y") +
  labs(
    x = "Tool",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

# Save the plot
ggsave(
  "resource_use_no_norm.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

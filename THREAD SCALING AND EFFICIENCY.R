# This is for making the thread scalling and efficiency graphs

library(ggplot2)

# Read files
gut <- read.table("GUT_resource_usage.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gut$sample <- "GUT"

clean <- read.table("CLEAN_resource_usage.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clean$sample <- "CLEAN"

marine <- read.table("MARINE_resource_usage.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
marine$sample <- "MARINE"

# Combine
data <- rbind(gut, clean, marine)

# Convert elapsed_wall (mm:ss.ss) to seconds
to_seconds <- function(x) {
  sapply(strsplit(x, ":"), function(t) {
    as.numeric(t[1]) * 60 + as.numeric(t[2])
  })
}

data$elapsed_wall_sec <- to_seconds(data$elapsed_wall)

# Normalize wall time within each sample
data$wall_norm <- ave(
  data$elapsed_wall_sec,
  data$sample,
  FUN = function(x) x / mean(x, na.rm = TRUE)
)

# Summarize normalized wall time across replicates and samples
wall_mean <- aggregate(wall_norm ~ threads + tool, data = data, FUN = mean)
wall_sd   <- aggregate(wall_norm ~ threads + tool, data = data, FUN = sd)

names(wall_mean)[3] <- "mean_value"
names(wall_sd)[3]   <- "sd_value"

plot_summary <- merge(wall_mean, wall_sd, by = c("threads", "tool"))

plot_summary$tool <- factor(
  plot_summary$tool,
  levels = c("fastp", "cutadapt", "trimmomatic", "skewer")
)

# Plot normalized wall time
p1 <- ggplot(plot_summary, aes(x = threads, y = mean_value, color = tool, group = tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
    width = 0.4
  ) +
  scale_x_continuous(breaks = sort(unique(plot_summary$threads))) +
  labs(
    x = "Threads",
    y = "Normalized wall time (mean ± SD)",
    color = "Tool"
  ) +
  theme_minimal()
print(p1)

# =========================
# Thread efficiency
# =========================

# Get 1-thread baseline for each tool
baseline <- subset(plot_summary, threads == 1, c("tool", "mean_value"))
names(baseline)[2] <- "baseline_time"

eff_summary <- merge(plot_summary, baseline, by = "tool")

# This is based on this https://stackoverflow.com/questions/29970813/parallel-speedup-and-efficiency#:~:text=3%20Comments&text=The%20efficiency%20measures%20how%20well,of%200.7%20(or%2070%25).
# Compute speedup and efficiency
# S = T_1 / T_p
eff_summary$speedup <- eff_summary$baseline_time / eff_summary$mean_value
# E = S / p *100
eff_summary$efficiency_pct <- (eff_summary$speedup / eff_summary$threads) * 100

# Errors
eff_summary$sd_efficiency_pct <- (eff_summary$sd_value / eff_summary$mean_value) *
                                eff_summary$efficiency_pct

# Plot thread efficiency
p2 <-  ggplot(eff_summary, aes(x = threads, y = efficiency_pct, color = tool, group = tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = pmax(efficiency_pct - sd_efficiency_pct, 0),
      ymax = efficiency_pct + sd_efficiency_pct
    ),
    width = 0.4
  ) +
  scale_x_continuous(breaks = sort(unique(eff_summary$threads))) +
  labs(
    x = "Threads",
    y = "Parallel Efficiency (%)",
    color = "Tool"
  ) +
  theme_minimal()
print(p2)

ggsave(
  "Wall time.png",
  plot = p1,
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  "Thread efficiency.png",
  plot = p2,
  width = 6,
  height = 4,
  dpi = 300
)

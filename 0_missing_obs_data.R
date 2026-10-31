root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

library(ggplot2)
library(ggpubr)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
START_DATE <- as.Date("1998-01-01")
END_DATE   <- as.Date("2020-12-31")
MON_MISS   <- 6:10   # months highlighted (rainy season)
SEL_TH     <- 5      # stations with <= SEL_TH % missing in all MON_MISS months are framed in blue

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------
print("============= Reading data =============")
stations <- read.csv("data/csv/p_stations.csv", header = TRUE, sep = ",", dec = ".",
                     check.names = FALSE, stringsAsFactors = FALSE)
stations$Station <- trimws(stations$Station)

pr <- read.csv("data/csv/p_stations_daily_data.csv", header = TRUE, sep = ",", dec = ".",
               check.names = FALSE, stringsAsFactors = FALSE)
colnames(pr)[1] <- "Date"
colnames(pr)[-1] <- trimws(colnames(pr)[-1])
pr$Date <- as.Date(pr$Date, format = "%m/%d/%Y")

# Filter on the study period and make sure every day of the period is present
# (days absent from the file are counted as missing)
pr <- pr[!is.na(pr$Date) & pr$Date >= START_DATE & pr$Date <= END_DATE, ]
full_dates <- data.frame(Date = seq(START_DATE, END_DATE, by = "day"))
pr <- merge(full_dates, pr, by = "Date", all.x = TRUE)

# Keep stations present in both files, ordered from north to south (then east to west)
stations <- stations[stations$Station %in% colnames(pr)[-1], ]
stations <- stations[!duplicated(stations$Station), ]
stations <- stations[order(stations$Latitude, stations$Longitude, decreasing = TRUE), ]
st_names <- stations$Station
n_st <- length(st_names)
print(paste0(n_st, " stations processed"))

# Label style per station type: Synoptic = bold blue, Other = plain black
is_syn     <- trimws(stations$Type) == "Synoptic"
lab_colors <- ifelse(is_syn, "blue", "black")
lab_faces  <- ifelse(is_syn, "bold", "plain")

# ------------------------------------------------------------------------------
# Missing data (%) per station, per month and per year
# ------------------------------------------------------------------------------
na_mat <- 1 * is.na(as.matrix(pr[, st_names, drop = FALSE]))

months <- as.integer(format(pr$Date, "%m"))
years  <- as.integer(format(pr$Date, "%Y"))
yrs    <- sort(unique(years))
n_yr   <- length(yrs)

mon_miss <- rowsum(na_mat, months) / as.vector(table(months)) * 100   # 12 x n_st
ann_miss <- rowsum(na_mat, years)  / as.vector(table(years))  * 100   # n_yr x n_st

df <- data.frame(
  sname = factor(rep(st_names, each = 12), levels = st_names),
  dmon  = factor(rep(1:12, times = n_st), levels = 1:12, labels = month.abb),
  mon   = rep(1:12, times = n_st),
  miss  = as.vector(mon_miss)
)

df2 <- data.frame(
  sname = factor(rep(st_names, each = n_yr), levels = st_names),
  dmon  = factor(rep(yrs, times = n_st), levels = yrs),
  miss  = as.vector(ann_miss)
)

# Labels: range of missing data (%) per station
fmt_range <- function(x) {
  paste0("   [", sprintf("%.1f", min(x)), "% - ", sprintf("%.1f", max(x)), "%]")
}
labs  <- apply(mon_miss, 2, fmt_range)   # monthly
ylabs <- apply(ann_miss, 2, fmt_range)   # annual

# Stations to frame (all highlighted months <= SEL_TH %)
filter_rect <- which(apply(mon_miss[MON_MISS, , drop = FALSE], 2, function(x) all(x <= SEL_TH)))
ys <- n_st - filter_rect

df.rect <- data.frame(
  xmin = rep(0.5, length(filter_rect)),
  xmax = rep(12.5, length(filter_rect)),
  ymin = ys + 0.5,
  ymax = ys + 1.5
)
df.rect2 <- df.rect
df.rect2$xmax <- n_yr + 0.5

tilew <- 1
tileh <- 1

# ------------------------------------------------------------------------------
# (a) Monthly gaps
# ------------------------------------------------------------------------------
pl.mon <- ggplot() +
  geom_tile(data = df, aes(x = dmon, y = sname, fill = miss, width = tilew, height = tileh),
            alpha = 0.3, color = "black", lwd = 0.4, linetype = 2) +
  scale_fill_gradientn(colours = terrain.colors(7), limits = c(0, 100)) +
  scale_x_discrete(labels = month.abb) +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  coord_cartesian(clip = "off") +
  xlab("\nMonths (-)") + ylab("Rain gauge stations\n") +
  geom_tile(data = df[df$mon %in% MON_MISS, ],
            aes(x = dmon, y = sname, fill = miss, width = tilew, height = tileh),
            alpha = 0.7, color = "black", lwd = 0.4, linetype = 1.5) +
  geom_rect(data = df.rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, alpha = 1, color = "blue", lwd = 1.5) +
  annotate("text", x = 13, y = 1:n_st, hjust = 0, label = rev(labs), size = 4.5,
           colour = "maroon", fontface = "bold.italic") +
  labs(title = "(a) Monthly gaps\n") +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 40, title = "Missing data (%)\n")) +
  theme(legend.background = element_blank(),
        legend.position = "right",
        axis.line = element_line(color = "black", linewidth = 0.8, linetype = "solid"),
        legend.title = element_text(face = "bold", color = "black", size = 20),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        axis.ticks = element_line(linewidth = 1.5), axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(face = "bold", size = 30, color = "black"),
        plot.margin = unit(c(1, 8, 1, 1), "lines"),
        legend.spacing.y = unit(1.0, "in"),
        axis.title.x = element_text(color = "black", size = 26, face = "bold"),
        axis.title.y = element_text(color = "black", size = 26, face = "bold"),
        axis.text.x = element_text(face = "bold", color = "#993333", size = 18,
                                   angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(face = rev(lab_faces), color = rev(lab_colors), size = 18))

# ------------------------------------------------------------------------------
# (b) Annual gaps
# ------------------------------------------------------------------------------
pl.ann <- ggplot() +
  geom_tile(data = df2, aes(x = dmon, y = sname, fill = miss, width = tilew, height = tileh),
            alpha = 0.7, color = "black", lwd = 0.4, linetype = 1.5) +
  scale_fill_gradientn(colours = terrain.colors(7), limits = c(0, 100)) +
  scale_x_discrete(labels = yrs) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  xlab("\nYears (-)") + ylab("\n") +
  geom_rect(data = df.rect2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, alpha = 1, color = "blue", lwd = 1.5) +
  annotate("text", x = n_yr + 1, y = 1:n_st, hjust = 0, label = rev(ylabs), size = 4.5,
           colour = "maroon", fontface = "bold.italic") +
  labs(title = "(b) Annual gaps\n") +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 40, title = "Missing \n data (%)\n")) +
  theme(legend.background = element_blank(),
        plot.title = element_text(face = "bold", size = 30, color = "black"),
        plot.margin = unit(c(1, 10, 1, 1), "lines"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.8, linetype = "solid"),
        legend.title = element_text(face = "bold", color = "black", size = 20),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        axis.title.x = element_text(color = "black", size = 26, face = "bold"),
        axis.title.y = element_text(color = "black", size = 26, face = "bold"),
        axis.text.x = element_text(face = "bold", color = "#993333", size = 18,
                                   angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

grob <- ggarrange(pl.mon, pl.ann, nrow = 1, ncol = 2, common.legend = TRUE,
                  widths = c(1, 1.4), legend = "right")

fig_height <- max(15, 0.6 * n_st + 6)

ggsave(paste0("figures/0_missing_data_", format(START_DATE, "%Y"), "_", format(END_DATE, "%Y"), ".png"),
       plot = grob, width = 40, height = fig_height, dpi = 400, scale = 0.7, bg = "white")

print("Finished.")

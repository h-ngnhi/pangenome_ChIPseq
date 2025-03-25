# Install packages first time
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("ChIPseeker")
# BiocManager::install("clusterProfiler")

library(ChIPseeker)
library(clusterProfiler)

wd <- "/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome"  
files <- readPeakFile("146_design2_hg19/peak_call/146/ZNF146/all_intersection.hg19.600.bed")
files
pdf(paste0(wd, "/plots/",Sys.Date(), ".coverage_146_19_L1_600_GenPipes.pdf"))
covplot(files, weightCol="V5", xlab="Chromosome", ylab="Coverage")
dev.off()

# Plot distribution of IDR score of peaks
bed <- as.data.frame(
  read.table("ZNF146-507-Analysis-on-Pangenome/146_results_hg19/peak_call/IDR/146.idr.filt.sorted.bed",
             header = FALSE, sep="\t",stringsAsFactors=FALSE))
pdf("ZNF146-507-Analysis-on-Pangenome/plots/frequency_146_idr_filt.pdf")
hist(bed$V5,main="Frequency of peaks coverage",xlab="IDR Score")
dev.off()


#####################################################################
# Plot the distribution of peak lengths
library(ggplot2)
library(gridExtra)
install.packages("ggraph")
library(coda)

type <- "design2"

# Create the directory list using the 'type' variable
directory_list <- c(paste0("146_", type, "_hg19"),
                    paste0("146_", type, "_hg38"),
                    paste0("507_", type, "_hg19"),
                    paste0("507_", type, "_hg38"))

wd <- "/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome"  # Set your working directory path

plots <- list()

# directory_list <- c("146_results_t2t", "507_results_t2t")

for (dir in directory_list) {
    ref <- unlist(strsplit(dir, "_"))[3]
    gene <- unlist(strsplit(dir, "_"))[1]
    path_idr <- paste0(wd, "/", dir, "/peak_call/", gene, "/ZNF", gene)
    print(ref)
    print(gene)
    encode_directory <- paste0(wd, "/GenPipes_sets/data/ENCODE/", gene)
    encode_all <- grep("sorted", list.files(path = encode_directory, pattern = sub("hg", "", ref) , full.names = TRUE), value = TRUE)
    genpipes_all <- paste0(path_idr, "/", gene, ".ZNF", gene,"_peaks.sorted.bed")
    print(encode_all)
    print(genpipes_all)

    # Read files and calculate peak sizes
    encode_df <- read.table(encode_all, header = FALSE)
    encode_df <- encode_df[,1:3]
    encode_df$peak_size <- encode_df$V3 - encode_df$V2

    genpipes_df <- read.table(genpipes_all, header = FALSE)
    genpipes_df <- genpipes_df[,1:3]
    genpipes_df$peak_size <- genpipes_df$V3 - genpipes_df$V2
    max_peak_size <- max(encode_df$peak_size, na.rm = TRUE)

#     percentiles <- quantile(genpipes_df$peak_size, probs = c(0.025, 0.975))
#     print(percentiles)
    hpd_interval <- HPDinterval(as.mcmc(genpipes_df$peak_size), prob = 0.90)
    print(hpd_interval)
    hpd_matrix <- as.matrix(hpd_interval)
    lower_bound <- hpd_matrix[1]
    upper_bound <- hpd_matrix[2]

    result_df <- data.frame(
        Result = type,
        Genome = ref,
        Gene = gene,
        Lower = lower_bound,
        Upper = upper_bound
    )
    # To append without repeating the header
    hpd_file=paste0(wd, "/plots/", Sys.Date(), "HPD_interval.csv")
    write.table(result_df, hpd_file, append = TRUE, sep = "\t", 
                col.names = !file.exists(hpd_file), row.names = FALSE, quote = TRUE)

    encode_df$Source <- 'Encode'
    genpipes_df$Source <- 'GenPipes'

    # Plotting
    combined_df <- rbind(encode_df, genpipes_df)
    p1 <- ggplot(combined_df, aes(x = peak_size, fill = Source)) +
      geom_density(aes(y = ..density..), alpha = 0.5) +
      scale_fill_manual(values = c(Encode = "orange", GenPipes = "blue")) +
      scale_x_log10() +
      ggtitle(paste("Peak Size Distribution for", gene, "in", ref)) +
      xlab("Peak Size") + ylab("Density") +
      geom_vline(xintercept = max_peak_size, color = "red", linetype = "dashed") +
      annotate("text", x = max_peak_size, y = 0, label = sprintf(paste0("Max ENCODE = ", max_peak_size)), vjust = -2, hjust = -0.05, color = "red") +
      geom_vline(xintercept = lower_bound, color = "gray40", linetype = "dotted") +
      geom_vline(xintercept = upper_bound, color = "gray40", linetype = "dotted") +
      annotate("text", x = upper_bound, y = 0, label = "90% GenPipes peaks", vjust = -10, hjust = -0.05, color = "gray40") +
      theme_minimal() +
      theme(legend.position = c(1, 1), # Coordinates are relative to the plot area (1,1) is top-right
            legend.justification = c(1, 1), # Justify legend position (top-right)
            legend.box.margin = margin(-2, -2, -5, -5),
            legend.title = element_blank(),
            legend.text = element_text(size = 8),  # Smaller legend text
            legend.key.size = unit(0.5, "cm")) # Negative margins to overlay on the plot



    # Plot the distribution of peak lengths
    # p1 <- ggplot(combined_df, aes(x = peak_size, fill = Source)) +
    #     geom_freqpoly(data = encode_df, aes(x = peak_size, color = "Encode"), binwidth = 1) +
    #     geom_freqpoly(data = genpipes_df, aes(x = peak_size, color = "GenPipes"), binwidth = 1) +
    #     scale_fill_manual(values = c(Encode = "orange", GenPipes = "blue")) +
    #     ggtitle(paste("Peak Size Distribution for", gene, "in", ref)) +
    #     xlab("Peak Size") +
    #     ylab("Count") +
    #     scale_x_log10() +
    #     geom_vline(xintercept = max_peak_size, color = "red", linetype = "dashed") +
    #     annotate("text", x = max_peak_size, y = 0, label = sprintf(paste0("Max ENCODE = ", max_peak_size)), vjust = -2, hjust = -0.05, color = "red") +
    #     geom_vline(xintercept = lower_bound, color = "gray40", linetype = "dotted") +
    #     geom_vline(xintercept = upper_bound, color = "gray40", linetype = "dotted") +
    #     annotate("text", x = upper_bound, y = 0, label = "90% GenPipes peaks", vjust = -10, hjust = -0.05, color = "gray40") +
    #     theme_minimal() +
    #     theme(legend.position = c(1, 1), # Coordinates are relative to the plot area (1,1) is top-right
    #             legend.justification = c(1, 1), # Justify legend position (top-right)
    #             legend.box.margin = margin(-2, -2, -5, -5),
    #             legend.title = element_blank(),
    #             legend.text = element_text(size = 8),  # Smaller legend text
    #             legend.key.size = unit(0.5, "cm")) # Negative margins to overlay on the plot

    plots[[length(plots) + 1]] <- p1
}

length(plots)
pdf(paste0(wd, "/plots/", Sys.Date(), type, ".peak_length_distribution.counts2.pdf"))
print(do.call(grid.arrange, c(plots, ncol=2))) # combine all plots
dev.off()

#####################################################################
# Plot the coverage of L1HS (L1PA1) peaks in meta analysis
wd <- "/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome"
# gene <- 146
# ref <- "hg19"
# pipeline <- "ENCODE"
genes <- c("146")
refs <- c("hg19", "hg38")  # This can be modified later to include "T2T"
pipelines <- c("ENCODE", "GenPipes")
type <- "design2"

# Nested for loops in R
for (gene in genes) {
  for (ref in refs) {
    for (pipeline in pipelines) { 
      # Construct the file name or directory path
      meta_directory <- paste0(wd, "/meta_analysis_files/", pipeline, "/", gene)
      coverage_file <- paste0(meta_directory, "/", ref, ".L1_coverage.txt")
      print(coverage_file)
      coverage_data <- read.csv(coverage_file, header = FALSE, sep="\t", comment.char="" )
      coverage_data <- coverage_data[,c(2,3)]
      head(coverage_data)
      colnames(coverage_data) <- c('Position', 'Depth')

      library(ggplot2)
      # pdf(paste0(meta_directory, "/", Sys.Date(), ".", ref, ".L1HS_peak_coverage.pdf"), width = 11, height = 4)
      p <- ggplot(coverage_data, aes(x = Position, y = Depth)) +
        geom_area(fill = "blue", alpha = 0.5) + 
        labs(
          title = paste0(pipeline, "_", gene, "_" , ref, " L1HS_coverage"),
          x = 'Position (bp)',
          y = 'Coverage'
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text = element_text(size = 10),
              axis.ticks = element_line(size = 0.5)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 10000))
      ggsave(paste0(meta_directory, "/", Sys.Date(), ".", ref, ".L1_peak_coverage.png"), plot = p, width = 11, height = 4)  
      # dev.off()
    }
  }
}


###################################################################
# Plot histogram of mapping quality
# Load required library
library(ggplot2)
wd <- "/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/"
# bam_files=($wd/146_design2_t2t/alignment/146/ZNF146/146.ZNF146.sorted.dup.filtered.cleaned.bam $wd/507_design2_t2t/alignment/507/ZNF507/507.ZNF507.sorted.dup.filtered.cleaned.bam)

# Define your working directory and the plot directory
plot_dir <- paste0(wd,"plots/")

# Define BAM files and genes
bam_files <- c(paste0(wd, "146_design2_t2t/alignment/146/ZNF146/146.ZNF146.sorted.dup.filtered.cleaned.bam.mapping_score.txt"), 
               paste0(wd, "507_design2_t2t/alignment/507/ZNF507/507.ZNF507.sorted.dup.filtered.cleaned.bam.mapping_score.txt"))
gene <- c("bwa mem", "vg giraffe")

for (i in seq_along(bam_files)) {
  bam_file <- bam_files[i]
  cat("Reading file:", bam_file, "\n")
  mapping_scores <- read.table(bam_file, header=FALSE, col.names=c("mapping_score"))

  cat("First few rows of data:\n")
  print(head(mapping_scores))

  # Plot the histogram
  # png(paste0(plot_dir, Sys.Date(), ".", gene[i], ".L1_mapping_score_actual.png"))
  # hist(mapping_scores$mapping_score, 
  #      main=paste("Histogram of Mapping Scores for", gene[i]), 
  #      xlab="Mapping Score", 
  #      ylab="Count", 
  #      col="blue", 
  #      border="black")

  # dev.off()
  # p <- ggplot(mapping_scores, aes(x=mapping_score)) +
  #   geom_histogram(binwidth=1, fill="blue", alpha=0.7, color="black") +
  #   labs(title="Histogram of Mapping Scores", x="Mapping Score", y="Count") +
  #   theme_minimal()

  # output_file <- paste0(plot_dir, Sys.Date(), ".", gene[i], ".L1_mapping_score.png")
  # cat("Saving plot to:", output_file, "\n")
  # ggsave(output_file, plot = p, width = 11, height = 4)
}



bam_files <- c(paste0(wd, "cgroza_data/H3K27ac_ENC/mapq_scores_ENC.txt"), 
               paste0(wd, "cgroza_data/H3K27ac_ENC/H3K27AC_CHM13linear/alignment/H3K27AC_treatment1/H3K27AC/mapq_scores2.txt"),
               paste0(wd, "cgroza_data/H3K27ac_ENC/H3K27AC_chm13_graph_giraffe/mapq_scores.txt"),
               paste0(wd, "cgroza_data/H3K27ac_ENC/H3K27AC_vcfbub_graph_giraffe/mapq_scores.txt"))
gene <- c("CHM13linear-bowtie2", "CHM13linear-bwamem", "CHM13graph-vg giraffe", "HPRC-vg giraffe") 
# Read the mapping scores from the file
for (i in seq_along(bam_files)) {
  bam_file <- bam_files[i]
  mapping_scores <- read.table(bam_file, header=FALSE, col.names=c("mapping_score"))

  # Print the first few rows of the data to verify
  head(mapping_scores)

  # Plot the histogram
  p <- ggplot(mapping_scores, aes(x=mapping_score)) +
    geom_histogram(binwidth=1, fill="blue", alpha=0.7, color="black") +
    labs(title="Histogram of Mapping Scores", x="Mapping Score", y="Count") +
    theme_minimal()

  ggsave(paste0("plots/", Sys.Date(), ".", gene[i], ".H3K27ac_mapping_score.png"), plot = p, width = 11, height = 4)  
}
print(gene[i])
mapping_scores <- read.table("cgroza_data/H3K27ac_ENC/H3K27AC_chm13_graph/mapq_scores.txt", header=FALSE, col.names=c("mapping_score"))

# Print the first few rows of the data to verify
head(mapping_scores)

library(data.table)
library(ggplot2)
wd <- "/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/"

# Use fread to read in large files efficiently
mapping_scores <- fread("cgroza_data/H3K27ac_ENC/H3K27AC_chm13_graph/mapq_scores.txt", header=FALSE, col.names=c("mapping_score"))

# Summarize counts for each mapping score to reduce the data size for plotting
mapping_score_counts <- mapping_scores[, .N, by = mapping_score]

# Plot the summarized data
p <-ggplot(mapping_score_counts, aes(x=mapping_score, y=N)) +
  geom_bar(stat="identity", fill="blue", alpha=0.7, color="black") +
  labs(title="Histogram of Mapping Scores", x="Mapping Score", y="Count") +
  theme_minimal()

png(filename = paste0("plots/", Sys.Date(), ".CHM13-vgmap-idvl.H3K27ac_mapping_score.png"),
    width = 1100, height = 400, res = 72) # Adjust resolution to manage memory
print(p)
dev.off()


# Plot the histogram
ggplot(mapping_scores, aes(x=mapping_score)) +
  geom_histogram(binwidth=1, fill="blue", alpha=0.7, color="black") +
  labs(title="Histogram of Mapping Scores", x="Mapping Score", y="Count") +
  theme_minimal()

ggsave(paste0("plots/", Sys.Date(), ".CHM13-vgmap-idvl.H3K27ac_mapping_score.png"), plot = p, width = 11, height = 4)
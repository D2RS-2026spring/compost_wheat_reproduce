# 基于Sustaining vulnerable agroecosystems with compost: Lasting benefits to soil health and carbon storage in semiarid winter wheat (Triticum aestivum, L.)文献的部分数据图复现

## 小组成员
* [@cvvcI](https://github.com/cvvcI) 程伟东2025303120043
* @YI-1005 魏伊洁2025303120064
* @JWenJing 焦文静2025303110019
* @Wangchenyang321 王晨阳2025303110073
* [@KINGWINM](https://github.com/KINGWINM)王胜2025303120135

## 复现说明
将该仓库克隆到本地后，修改并运行以下代码：

###  1. 初始化环境

if (!require("renv")) install.packages("renv")
renv::init() 

packages <- c("dada2", "phyloseq", "microeco", "ALDEx2", "ggplot2", "vegan", "tidyverse")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

renv::snapshot()

library(dada2)
library(phyloseq)
library(microeco)
library(ALDEx2)
library(ggplot2)
library(vegan)
library(tidyverse)

###  2. 设置工作路径和数据 

请根据实际情况修改以下路径：
fastq_path <- "path/to/your/fastq/files"     # 存放原始 fastq 的文件夹
metadata_path <- "path/to/metadata.csv"      # 样本元数据（含处理信息、土壤类型等）
output_dir <- "results"                      # 结果输出目录
dir.create(output_dir, showWarnings = FALSE)

### 3. DADA2 流程（生成 ASV 表）

假设双端测序，文件命名格式为：sample_R1.fastq, sample_R2.fastq
fnFs <- sort(list.files(fastq_path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(fastq_path, pattern = "_R2.fastq", full.names = TRUE))
sample_names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

filtFs <- file.path(output_dir, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(output_dir, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(240, 200),   # 根据质量图调整
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

write.csv(t(seqtab.nochim), file = file.path(output_dir, "asv_table.csv"))

###  4. 创建 phyloseq 对象

metadata <- read.csv(metadata_path, row.names = 1)
metadata$SampleID <- rownames(metadata)

asv_table <- read.csv(file.path(output_dir, "asv_table.csv"), row.names = 1)
asv_table <- as.matrix(asv_table)

taxa <- assignTaxonomy(seqtab.nochim, "path/to/silva_nr_v138_train_set.fa.gz", multithread = TRUE)

OTU <- otu_table(asv_table, taxa_are_rows = FALSE)
TAX <- tax_table(taxa)
SAM <- sample_data(metadata)

physeq <- phyloseq(OTU, TAX, SAM)

saveRDS(physeq, file = file.path(output_dir, "phyloseq_object.rds"))

###  5. 多样性分析 

rarecurve(t(asv_table), step = 100, sample = min(rowSums(asv_table)))

alpha_div <- estimate_richness(physeq, measures = c("Observed", "Chao1", "Shannon"))
alpha_div$SampleID <- rownames(alpha_div)
alpha_div <- left_join(alpha_div, metadata, by = "SampleID")

p_alpha <- ggplot(alpha_div, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot() +
  facet_wrap(~SoilType) +
  theme_bw() +
  labs(title = "Alpha diversity (Shannon)")
ggsave(file.path(output_dir, "alpha_diversity.pdf"), p_alpha, width = 8, height = 6)

ord <- ordinate(physeq, method = "PCoA", distance = "bray")
p_beta <- plot_ordination(physeq, ord, color = "Treatment", shape = "SoilType") +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Treatment), linetype = 2) +
  theme_bw() +
  labs(title = "PCoA - Bray-Curtis")
ggsave(file.path(output_dir, "beta_diversity_pcoa.pdf"), p_beta, width = 8, height = 6)

dist <- phyloseq::distance(physeq, method = "bray")
adonis_result <- adonis2(dist ~ SoilType * Treatment, data = metadata, permutations = 999)
write.csv(as.data.frame(adonis_result), file = file.path(output_dir, "permanova_results.csv"))

###  6. 差异丰度分析

micro_table <- microtable$new(otu_table = as.data.frame(asv_table),
                              sample_table = metadata,
                              tax_table = as.data.frame(taxa))
micro_table$sample_sums()   # 查看测序深度
micro_table$rarefy_samples(sample.size = min(micro_table$sample_sums))

diff_test <- micro_table$diff_analysis(method = "ALDEx2", group = "Treatment")
diff_res <- diff_test$res_diff
write.csv(diff_res, file = file.path(output_dir, "microeco_diff.csv"))

###  7. 功能预测

asv_seq <- colnames(asv_table)
names(asv_seq) <- asv_seq
writeXStringSet(DNAStringSet(asv_seq), filepath = file.path(output_dir, "asv_seqs.fasta"))

picrust2_pipeline.py -s asv_seqs.fasta -o picrust2_out -p 4

kegg_path <- file.path(output_dir, "picrust2_out", "KO_metagenome_out", "pred_metagenome_unstrat.tsv.gz")
if(file.exists(kegg_path)){
  kegg_table <- read.table(kegg_path, header = TRUE, row.names = 1, sep = "\t")
 
  write.csv(kegg_table, file = file.path(output_dir, "kegg_abundance.csv"))
}

###  8. 统计分析

soil_data <- read.csv("path/to/soil_attributes.csv")  # 包含 SampleID, SOC, C_POM, C_MAOM 等

model <- aov(SOC ~ SoilType * Treatment * Form, data = soil_data)
summary_model <- summary(model)
write.csv(as.data.frame(summary_model[[1]]), file = file.path(output_dir, "anova_results.csv"))

###  9. 生成最终报告

sink(file.path(output_dir, "session_info.txt"))
sessionInfo()
sink()

cat("分析完成！结果保存在", output_dir, "文件夹中。\n")

## 其他
* 本仓库仅为提交D2RS课程作业
* 除运行以上代码，也可在直接加载相应R包后直接运行库中R代码

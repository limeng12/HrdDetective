#' Convert BAM files to seqz format for SNP-based analysis
#'
#' This function processes paired normal and tumor BAM files to generate
#' seqz format files suitable for heterozygosity and copy number analysis.
#' It extracts SNP positions, calculates allele frequencies, and performs
#' GC-content normalization.
#'
#' @param normal_bam Path to the normal (germline) sample BAM file.
#' @param tumor_bam Path to the tumor sample BAM file.
#' @param genome_fasta Path to the reference genome FASTA file.
#'   Either `genome_fasta` or `bsgenome` must be provided.
#' @param bsgenome A BSgenome object (e.g., `BSgenome.Hsapiens.UCSC.hg38`).
#'   Either `genome_fasta` or `bsgenome` must be provided.
#' @param min_depth Minimum read depth required for a position to be included.
#'   Default is 10.
#' @param min_af Minimum allele frequency for a position to be considered heterozygous.
#'   Default is 0.25.
#' @param max_af Maximum allele frequency for a position to be considered heterozygous.
#'   Default is 0.75.
#' @param gc_window Size of the window (in bases) for GC-content calculation.
#'   The actual window is `2 * gc_window + 1`. Default is 50 (±250 => 501 bp).
#' @param output_file Path to the output seqz file. Default is "output.seqz".
#'
#' @return Invisibly returns the path to the output seqz file.
#'   The function primarily writes output to the specified file.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads reference genome sequence (from FASTA or BSgenome)
#'   \item Identifies heterozygous SNP positions in the normal sample
#'   \item Extracts allele counts from tumor sample at SNP positions
#'   \item Calculates GC-content in windows around each SNP
#'   \item Performs GC-bias normalization
#'   \item Outputs results in seqz format
#' }
#'
#' @note
#' \itemize{
#'   \item BAM files must be sorted and indexed (require .bai files)
#'   \item For large genomes, parallel processing is recommended
#'   \item The seqz format is compatible with tools like Sequenza and ASCAT
#' }
#'
#' @examples
#' \dontrun{
#' # Using BSgenome
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bam2seqz_r_snps(
#'   normal_bam = "normal.bam",
#'   tumor_bam = "tumor.bam",
#'   bsgenome = BSgenome.Hsapiens.UCSC.hg38,
#'   output_file = "sample.seqz.txt"
#' )
#'
#' # Using FASTA file
#' bam2seqz_r_snps(
#'   normal_bam = "normal.bam",
#'   tumor_bam = "tumor.bam",
#'   genome_fasta = "hg38.fa",
#'   min_depth = 20,
#'   output_file = "sample.seqz.txt",
#'   BPPARAM = BiocParallel::MulticoreParam(workers = 4)
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[Rsamtools]{scanBam}} for BAM file reading
#'   \item \code{\link[Biostrings]{readDNAStringSet}} for FASTA reading
#'   \item \code{\link[sequenza]{sequenza.extract}} for alternative implementation
#' }
#'
#' @importFrom Rsamtools scanBam ScanBamParam
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom GenomicRanges GRanges seqnames start end
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom utils write.table
#' @export


bam2seqz_r_snps <- function(
    normal_bam,
    tumor_bam,
    genome_fasta = NULL,
    bsgenome = NULL,
    min_depth = 10,
    min_af = 0.25,
    max_af = 0.75,
    gc_window = 50,  # ±250 => 501 bp
    output_file = "output.seqz"
) {

  if (is.null(bsgenome) && is.null(genome_fasta)) {
    stop("Provide either 'bsgenome' or 'genome_fasta'")
  }

  # === Step 1: Load genome ===
  if (!is.null(genome_fasta)) {
    cat("Loading genome from FASTA...\n")
    genome_seq <- readDNAStringSet(genome_fasta)
    # Ensure seqnames match BAM (remove 'chr' if needed)
    names(genome_seq) <- sub("^chr", "", names(genome_seq))
    seqlengths <- width(genome_seq)
    seqinfo_genome <- Seqinfo(names(seqlengths), seqlengths)
  } else {
    genome_seq <- bsgenome
    seqinfo_genome <- seqinfo(bsgenome)
  }
  
  # === Step 2: Get pileup at all covered positions in normal ===
  cat("Extracting pileup from normal BAM...\n")
  # param <- PileupParam(
  #   distinguish_strands = FALSE,
  #   min_base_quality = 20,
  #   min_mapq = 20,
  #   min_nucleotide_depth=10,
  #   max_depth = 10000
  # )
  
  normal_pile <- pileup_whole_bam(normal_bam, min_depth=min_depth, min_af=min_af );
  
  # 假设 normal_pile 是 pileup() 的结果
  df_norm_long <- as.data.frame(normal_pile)

  # 确保有必要的列
  head(df_norm_long)
  # 应该包含: seqnames, pos, nucleotide, count

  # 转为宽格式：每个 pos 一行，A/C/G/T 为列
  df_norm <- df_norm_long %>%
    select(seqnames, pos, nucleotide, count) %>%
    # 先确保 nucleotide 是字符（避免因子问题）
    mutate(nucleotide = as.character(nucleotide)) %>%
    pivot_wider(
      names_from = nucleotide,
      values_from = count,
      values_fn = sum,          # 合并重复项（虽然通常不会重复）
      values_fill = 0L          # 注意：这里直接用 0L（新版本 tidyr 支持 scalar fill）
    ) %>%
    # 补全可能缺失的 A/C/G/T 列
    mutate(
      across(all_of(c("A", "C", "G", "T")), ~ replace_na(.x, 0L)),
      chromosome = as.character(seqnames),
      position = pos,
      total_depth = rowSums(select(., all_of(c("A", "C", "G", "T"))), na.rm = TRUE)
    )
  # 现在 df_norm 有 A, C, G, T 列了！




  # Convert to data.frame
  #df_norm <- as.data.frame(normal_pile)
  #df_norm$chromosome <- as.character(df_norm$seqnames)
  #df_norm$position <- df_norm$pos

  # === Step 3: Identify candidate heterozygous SNPs in normal ===
  bases <- c("A", "C", "G", "T")
  #df_norm$total_depth <- rowSums(df_norm[, bases], na.rm = TRUE)
  df_norm <- df_norm[df_norm$total_depth >= min_depth, ]

  # Compute allele frequencies
  af_matrix <- t(apply(df_norm[, bases], 1, function(x) {
    x / sum(x, na.rm = TRUE)
  }))
  colnames(af_matrix) <- paste0("af_", bases)

  df_norm <- cbind(df_norm, af_matrix)

  # Find second highest AF (BAF)
  # get_baf <- function(row_af) {
  #   sorted_af <- sort(row_af, decreasing = TRUE)
  #   if (length(sorted_af) < 2) return(NA_real_)
  #   sorted_af[2]
  # }
  # df_norm$baf_normal <- apply(af_matrix, 1, get_baf)
  ##############################
  # Step 1: 找到每行最大值的位置
  idx1 <- max.col(af_matrix, ties.method = "first")

  # Step 2: 创建一个副本，把最大值设为 -Inf（排除）
  af2 <- af_matrix
  af2[cbind(seq_len(nrow(af2)), idx1)] <- -Inf

  # Step 3: 第二大 = 新的最大值
  baf_normal <- apply(af2, 1, max, na.rm = TRUE)

  # Step 4: 处理无效行（原本少于2个非NA值）
  n_non_na <- rowSums(!is.na(af_matrix))
  baf_normal[n_non_na < 2] <- NA_real_

  # 替换可能的 -Inf（全 NA 行）
  baf_normal[is.infinite(baf_normal)] <- NA_real_

  df_norm$baf_normal<-baf_normal

  ##############################

  # Filter heterozygous-like sites
  het_snps <- df_norm[
    df_norm$baf_normal >= min_af & df_norm$baf_normal <= max_af,
    c("chromosome", "position", "total_depth", "baf_normal")
  ]
  cat("Identified", nrow(het_snps), "candidate heterozygous SNP sites.\n")

  if (nrow(het_snps) == 0) stop("No heterozygous SNPs found!")

  # Create GRanges for SNP positions
  snp_gr <- with(het_snps, GRanges(
    seqnames = chromosome,
    ranges = IRanges(position, width = 1),
    depth.normal = total_depth,
    baf.normal = baf_normal
  ))
  #seqinfo(snp_gr) <- seqinfo_genome

  
  # 提取共有染色体
  common_chr <- intersect(seqlevels(seqinfo_genome), as.character(unique(seqnames(snp_gr))))
  
  # 过滤并设置seqinfo
  snp_gr <- snp_gr[as.character(seqnames(snp_gr)) %in% common_chr]
  seqlevels(snp_gr) <- common_chr
  seqinfo(snp_gr) <- seqinfo_genome[common_chr]
  
  
  # === Step 4: Extract tumor depth and BAF at these SNP positions ===
  cat("Extracting tumor pileup at SNP sites...\n")
  tumor_pile <- pileup_whole_bam(tumor_bam,min_depth,min_af );
  #df_tumor <- as.data.frame(tumor_pile)
  #df_tumor$total_depth <- rowSums(df_tumor[, bases], na.rm = TRUE)
  
  df_tumor_long <- as.data.frame(tumor_pile)

  # 确保有必要的列
  head(df_tumor_long)
  # 应该包含: seqnames, pos, nucleotide, count

  # 转为宽格式：每个 pos 一行，A/C/G/T 为列
  df_tumor <- df_tumor_long %>%
    select(seqnames, pos, nucleotide, count) %>%
    # 先确保 nucleotide 是字符（避免因子问题）
    mutate(nucleotide = as.character(nucleotide)) %>%
    pivot_wider(
      names_from = nucleotide,
      values_from = count,
      values_fn = sum,          # 合并重复项（虽然通常不会重复）
      values_fill = 0L          # 注意：这里直接用 0L（新版本 tidyr 支持 scalar fill）
    ) %>%
    # 补全可能缺失的 A/C/G/T 列
    mutate(
      across(all_of(c("A", "C", "G", "T")), ~ replace_na(.x, 0L)),
      chromosome = as.character(seqnames),
      position = pos,
      total_depth = rowSums(select(., all_of(c("A", "C", "G", "T"))), na.rm = TRUE)
    )

  df_tumor <- df_tumor[df_tumor$total_depth >= min_depth, ]



  af_tumor <- t(apply(df_tumor[, bases], 1, function(x) x / sum(x, na.rm = TRUE)))
  # baf_tumor <- apply(af_tumor, 1, get_baf)


  ##################################
  # 确保是 matrix
  af_tumor <- as.matrix(af_tumor)

  # 步骤 1: 找到每行最大值所在的列
  idx_max <- max.col(af_tumor, ties.method = "first")

  # 步骤 2: 将每行的最大值临时设为 NA（或 -Inf）
  af_tumor_no_max <- af_tumor
  af_tumor_no_max[cbind(seq_len(nrow(af_tumor)), idx_max)] <- NA

  # 步骤 3: 每行取剩余中的最大值 → 即第二大值
  baf_tumor <- apply(af_tumor_no_max, 1, max, na.rm = TRUE)

  # 步骤 4: 处理无效行（原始行中非 NA 值 < 2）
  n_non_na <- rowSums(!is.na(af_tumor))
  baf_tumor[n_non_na < 2] <- NA_real_

  # （可选）清理可能的 NaN 或 -Inf
  baf_tumor[!is.finite(baf_tumor)] <- NA_real_

  ##################################

  # Match by position
  tumor_df <- data.frame(
    chromosome = as.character(df_tumor$seqnames),
    position = df_tumor$pos,
    depth.tumor = df_tumor$total_depth,
    baf.tumor = baf_tumor
  )

  # Merge normal and tumor
  merged <- merge(
    het_snps,
    tumor_df,
    by = c("chromosome", "position"),
    all.x = TRUE
  )
  merged[is.na(merged$depth.tumor), "depth.tumor"] <- 0
  merged[is.na(merged$baf.tumor), "baf.tumor"] <- NA_real_


  # === Step 5: Compute local GC content (±250 bp) for each SNP ===
  cat("Computing local GC content (window =", gc_window, "bp)...\n")
  
  
  valid_chroms <- seqlevels(genome_seq)
  
  # 过滤 merged 数据，只保留有效的染色体
  merged_filtered <- merged[merged$chromosome %in% valid_chroms, ]
  # 从 merged 创建 GRanges
  merged_gr <- GRanges(
    seqnames = merged_filtered$chromosome,
    ranges = IRanges(merged_filtered$position, width = 1)
  )
  merged<-merged_filtered
  
  # 统一染色体顺序
  #merged_gr <- keepSeqlevels(merged_gr, intersect(seqlevels(merged_gr),
  # seqlevels(genome_seq)), pruning.mode = "coarse")
  
  seqlevels(merged_gr) <- seqlevels(genome_seq)
  seqinfo(merged_gr) <- seqinfo(genome_seq)
  
  # 计算 GC 含量
  snp_extended <- resize(merged_gr, width = gc_window, fix = "center")
  snp_extended <- trim(snp_extended)
  seqs <- getSeq(genome_seq, snp_extended)
  merged$GC <- rowSums(letterFrequency(seqs, letters = c("G", "C"), as.prob = TRUE))
  merged$GC[is.nan(merged$GC)] <- NA_real_
  
  # 获取参考碱基
  merged$base.ref <- as.character(getSeq(genome_seq, merged_gr))
  
  

  # >>> 新增：确保染色体无 "chr" 前缀 <<<
  #merged$chromosome <- sub("^chr", "", merged$chromosome)

  # === Step 6: Prepare final .seqz data frame in Sequenza format ===

  # 获取 normal 和 tumor 的 A/C/G/T counts（你需要确保 df_norm 和 df_tumor 有这些列）
  bases <- c("A", "C", "G", "T")

  # 确保 count 列存在（用 0 填补缺失）
  for (b in bases) {
    if (!b %in% names(df_norm)) df_norm[[b]] <- 0L
    if (!b %in% names(df_tumor)) df_tumor[[b]] <- 0L
  }

  # 合并 counts 到 merged（按 position 匹配）
  count_cols <- c("chromosome", "position", bases)
  norm_counts <- df_norm %>% select(all_of(count_cols))
  tumor_counts <- df_tumor %>% select(all_of(count_cols))

  merged <- merge(merged, norm_counts, by = c("chromosome", "position"))
  merged <- merge(merged, tumor_counts, by = c("chromosome", "position"), suffixes = c(".normal", ".tumor"))

  # 计算 Af (normal 非参考 AF) 和 Bf (tumor 非参考 AF)
  # 构建 ref count 向量
  ref_count_normal <- merged[cbind(1:nrow(merged),
                                   match(paste0(merged$base.ref, ".normal"), names(merged)))]

  #ref_count_tumor <- merged[cbind(1:nrow(merged),
  #                                match(paste0(merged$base.ref, ".tumor"), names(merged)))]

  # 总深度（你已有）
  total_normal <- merged$total_depth  # 或 rowSums(merged[, paste0(bases, ".normal")])
  #total_tumor <- merged$depth.tumor

  # 非参考 AF
  Af <- ifelse(total_normal > 0, (total_normal - as.numeric(ref_count_normal) ) / total_normal, 0)
  # Bf <- ifelse(total_tumor > 0, (total_tumor - as.numeric(ref_count_tumor)) / total_tumor, 0)

  # 处理 NA（如果 match 失败）
  Af[is.na(Af)] <- 0
  # Bf[is.na(Bf)] <- 0

  #Af<-merged$baf_normal
  Bf<-merged$baf.tumor


  # zygosity.normal
  #het_threshold <- 0.25
  #zygosity.normal <- ifelse(Af >= min_af & Af <= (max_af), "het", "hom")
  zygosity.normal<-"het"
  # AB.normal: normal 中 top2 碱基（count > 0）
  get_AB_normal <- function(counts_vec) {
    nz <- counts_vec[counts_vec > 0]
    if (length(nz) == 0) return("")
    top2 <- names(sort(nz, decreasing = TRUE))[1:2]
    paste0(sort(top2), collapse = "")
  }
  AB.normal <- apply(merged[, paste0(bases, ".normal")], 1, get_AB_normal)

  # AB.tumor: 如果 Bf < 0.05，则输出 "X0.Y"，否则 "."
  # AB.tumor <- rep(".", length(Bf))
  # low_Bf <- Bf < 0.05
  # if (any(low_Bf)) {
  #   # 找非参考碱基（最大非 ref）
  #   find_alt_base <- function(ref, counts) {
  #     counts[names(counts) == ref] <- 0
  #     names(which.max(counts))
  #   }
  #   alt_bases <- mapply(find_alt_base,
  #                       merged$base.ref[low_Bf],
  #                       split(merged[low_Bf, paste0(bases, ".tumor")], 1:sum(low_Bf)))
  #   # 格式: "A0.1"
  #   AB.tumor[low_Bf] <- sprintf("%s0.%s", alt_bases, round(Bf[low_Bf] * 10))
  # }

  AB.tumor <- apply(merged[, paste0(bases, ".tumor")], 1, get_AB_normal)



  # depth.ratio
  depth.ratio <- ifelse(merged$total_depth > 0,
                        merged$depth.tumor / merged$total_depth,
                        0)

  # good.reads: 使用 depth.normal（或可替换为有效深度）
  good.reads <- merged$total_depth

  # tumor.strand: 固定为 0
  tumor.strand <- 0


  merged$chromosome <- sub("^chr", "", merged$chromosome)

  # 构建最终数据框（14列，顺序严格匹配）
  seqz_out <- data.frame(
    chromosome = merged$chromosome,
    position = merged$position,
    base.ref = merged$base.ref,
    depth.normal = merged$total_depth,
    depth.tumor = merged$depth.tumor,
    depth.ratio = round(depth.ratio, 3),
    Af = round(Af, 3),
    Bf = round(Bf, 3),
    zygosity.normal = zygosity.normal,
    GC.percent = round(merged$GC * 100, 1),  # 转为百分比
    good.reads = good.reads,
    AB.normal = AB.normal,
    AB.tumor = AB.tumor,
    tumor.strand = tumor.strand,
    stringsAsFactors = FALSE
  )

  # 移除 GC 为 NA 的行
  seqz_out <- seqz_out[!is.na(seqz_out$GC.percent), ]

  # === Step 7: Write output (NO HEADER!) ===
  cat("Writing", nrow(seqz_out), "SNP sites to", output_file, "\n")
  write.table(
    seqz_out,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE  # 关键：Sequenza .seqz 无 header
  );
  invisible(seqz_out)


}

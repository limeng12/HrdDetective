# pileup_wrapper.R
#' 从BAM文件生成pileup统计
#'
#' @param bam_file BAM文件路径
#' @param region 基因组区域，格式为"chr:start-end"或GRanges对象
#' @param min_mapq 最小MAPQ质量分数
#' @param min_base_quality 最小碱基质量分数
#' @param min_depth 最小深度
#' @param max_depth 最大深度
#' @param distinguish_strands 是否区分正负链
#' @param distinguish_nucleotides 是否区分不同碱基
#' @param ignore_query_Ns 是否忽略query序列中的N
#' @param include_deletions 是否包含缺失
#' @param include_insertions 是否包含插入
#' @param verbose 是否显示详细信息
#'
#' @return 包含pileup结果的DataFrame
#' @export
pileup_from_bam <- function(bam_file,
                            region = NULL,
                            min_mapq = 0,
                            min_base_quality = 0,
                            min_depth = 0,
                            min_af=0.2,
                            max_depth = 1000000,
                            distinguish_strands = FALSE,
                            distinguish_nucleotides = TRUE,
                            ignore_query_Ns = TRUE,
                            include_deletions = FALSE,
                            include_insertions = FALSE,
                            verbose = TRUE) {

  # 加载必要包
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("需要安装Rsamtools包: install.packages('Rsamtools')")
  }

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("需要安装GenomicRanges包: install.packages('GenomicRanges')")
  }

  if (!file.exists(bam_file)) {
    stop("BAM文件不存在: ", bam_file)
  }

  # 检查C++函数是否可用
  if (!exists("pileup_manual_vectors")) {
    stop("请先编译C++代码: sourceCpp('pileup_fixed.cpp')")
  }

  # 解析区域参数
  parse_region <- function(region) {
    if (is.null(region)) {
      # 如果没有指定区域，返回整个文件的范围
      # 获取参考序列信息
      bam_header <- Rsamtools::scanBamHeader(bam_file)
      ref_lengths <- bam_header[[1]]$targets

      regions <- list()
      for (i in seq_along(ref_lengths)) {
        regions[[i]] <- list(
          seqname = names(ref_lengths)[i],
          start = 1,
          end = ref_lengths[i]
        )
      }
      return(regions)
    }

    if (is.character(region)) {
      # 解析"chr:start-end"格式
      regions <- list()
      for (r in region) {
        parts <- strsplit(r, "[:-]")[[1]]
        if (length(parts) != 3) {
          stop("区域格式应为'chr:start-end'")
        }
        regions <- c(regions, list(list(
          seqname = parts[1],
          start = as.integer(parts[2]),
          end = as.integer(parts[3])
        )))
      }
      return(regions)
    }

    if (inherits(region, "GRanges")) {
      # 从GRanges对象解析
      regions <- list()
      for (i in seq_along(region)) {
        regions[[i]] <- list(
          seqname = as.character(GenomicRanges::seqnames(region)[i]),
          start = GenomicRanges::start(region)[i],
          end = GenomicRanges::end(region)[i]
        )
      }
      return(regions)
    }

    stop("不支持的region格式")
  }

  # 处理单个区域
  process_single_region <- function(seqname, start, end, chunk_size = 50000000) {
    if (verbose) message(sprintf("处理区域: %s:%d-%d (长度: %d)", 
                                 seqname, start, end, end - start + 1))
    
    region_length <- end - start + 1
    
    if (region_length <= chunk_size) {
      return(process_chunk(seqname, start, end))
    }
    
    # 创建分段
    chunk_starts <- seq(start, end, by = chunk_size)
    chunk_ends <- pmin(chunk_starts + chunk_size - 1, end)
    
    # 处理每个分段
    results <- lapply(seq_along(chunk_starts), function(i) {
      if (verbose) message(sprintf("  分段 %d/%d: %s:%d-%d", 
                                   i, length(chunk_starts),
                                   seqname, chunk_starts[i], chunk_ends[i]))
      process_chunk(seqname, chunk_starts[i], chunk_ends[i])
    })
    
    # 合并结果
    results <- results[!sapply(results, is.null)]
    if (length(results) == 0) return(NULL)
    
    do.call(rbind, results)
  }
  
  process_chunk <- function(seqname, start, end) {
    region_gr <- GenomicRanges::GRanges(
      seqnames = seqname,
      ranges = IRanges::IRanges(start = start, end = end)
    )
    
    param <- Rsamtools::ScanBamParam(
      what = c("pos", "cigar", "seq", "qual", "mapq", "flag", "qwidth"),
      which = region_gr,
      mapqFilter = min_mapq
    )
    
    bam_data <- Rsamtools::scanBam(bam_file, param = param)[[1]]
    if (length(bam_data$pos) == 0) return(NULL)
    
    if (verbose) message(sprintf("    找到 %d 条reads", length(bam_data$pos)))
    
    result <- pileup_manual_vectors(
      seqnames = rep(seqname, length(bam_data$pos)),
      positions = bam_data$pos,
      cigars = bam_data$cigar,
      sequences = as.character(bam_data$seq),
      qualities = as.character(bam_data$qual),
      mapqs = bam_data$mapq,
      flags = bam_data$flag,
      qwidths = bam_data$qwidth,
      seqname = seqname,
      start = start,
      end = end,
      min_mapq = min_mapq,
      min_base_quality = min_base_quality,
      min_depth = min_depth,
      max_depth = max_depth,
      min_af = min_af,
      distinguish_strands = distinguish_strands,
      distinguish_nucleotides = distinguish_nucleotides,
      ignore_query_Ns = ignore_query_Ns,
      include_deletions = include_deletions,
      include_insertions = include_insertions
    )
    
    if (length(result$seqnames) == 0) return(NULL)
    
    df <- data.frame(
      seqnames = result$seqnames,
      pos = result$pos,
      strand = result$strand,
      nucleotide = result$nucleotide,
      count = result$count,
      stringsAsFactors = FALSE
    )
    
    df[order(df$pos, df$strand, df$nucleotide), , drop = FALSE]
  }
  
  # 使用示例：
  # bam_file <- "your.bam"
  # verbose <- TRUE
  # min_mapq <- 20
  # result <- process_single_region("chr1", 1, 10000000, chunk_size = 500000)

  # 解析所有区域
  regions <- parse_region(region)

  # 处理所有区域
  results <- list()

  for (i in seq_along(regions)) {
    reg <- regions[[i]]

    if (verbose) {
      message("\n处理第 ", i, "/", length(regions), " 个区域")
    }

    result_df <- process_single_region(
      seqname = reg$seqname,
      start = reg$start,
      end = reg$end
    )

    if (!is.null(result_df)) {
      results[[i]] <- result_df
    }
  }

  # 合并所有结果
  if (length(results) == 0) {
    if (verbose) message("没有找到任何pileup数据")
    return(data.frame(
      seqnames = character(0),
      pos = integer(0),
      strand = character(0),
      nucleotide = character(0),
      count = integer(0)
    ))
  }

  final_result <- do.call(rbind, results)
  rownames(final_result) <- NULL

  if (verbose) {
    message("\n=== 处理完成 ===")
    message("总行数: ", nrow(final_result))
    message("位置范围: ",
            min(final_result$pos), "-",
            max(final_result$pos))
    if (distinguish_nucleotides) {
      message("碱基分布:")
      print(table(final_result$nucleotide))
    }
    if (distinguish_strands) {
      message("链分布:")
      print(table(final_result$strand))
    }
  }

  return(final_result)
}

#' 为整个BAM文件生成pileup
#'
#' @param bam_file BAM文件路径
#' @param ... 传递给pileup_from_bam的参数
#'
#' @return 包含全基因组pileup的DataFrame
#' @export
pileup_whole_bam <- function(bam_file,verbose=TRUE, ...) {
  if (verbose) {
    message("为整个BAM文件生成pileup...")
    message("这可能需要一些时间，取决于BAM文件大小")
  }

  # 获取所有染色体信息
  bam_header <- Rsamtools::scanBamHeader(bam_file)
  chromosomes <- names(bam_header[[1]]$targets)

  if (verbose) {
    message("找到 ", length(chromosomes), " 条染色体")
  }

  # 分染色体处理，避免内存溢出
  all_results <- list()

  for (chr in chromosomes) {
    if (verbose) {
      message("\n处理染色体: ", chr)
    }

    # 创建染色体区域
    chr_region <- paste0(chr, ":1-", bam_header[[1]]$targets[chr])

    # 处理该染色体
    chr_result <- pileup_from_bam(
      bam_file = bam_file,
      region = chr_region,
      verbose = FALSE,  # 子函数不显示详细信息
      ...
    )

    if (nrow(chr_result) > 0) {
      all_results[[chr]] <- chr_result
    }

    # 可选：每处理完一条染色体后保存中间结果
    if (exists("save_intermediate") && save_intermediate) {
      saveRDS(chr_result, paste0("pileup_", chr, ".rds"))
    }
  }

  # 合并所有结果
  if (length(all_results) == 0) {
    return(data.frame(
      seqnames = character(0),
      pos = integer(0),
      strand = character(0),
      nucleotide = character(0),
      count = integer(0)
    ))
  }

  final_result <- do.call(rbind, all_results)
  rownames(final_result) <- NULL

  if (verbose) {
    message("\n=== 全基因组pileup完成 ===")
    message("总行数: ", nrow(final_result))
    message("染色体数: ", length(unique(final_result$seqnames)))
  }

  return(final_result)
}


pileup_bam_from_chr <- function(bam_file,chr,verbose=TRUE, ...) {
  if (verbose) {
    message("为整个BAM文件生成pileup...")
    message("这可能需要一些时间，取决于BAM文件大小")
  }
  
  # 获取所有染色体信息
  bam_header <- Rsamtools::scanBamHeader(bam_file)
  chromosomes <- names(bam_header[[1]]$targets)
  
  if (verbose) {
    message("找到 ", length(chromosomes), " 条染色体")
  }
  
  # 分染色体处理，避免内存溢出
  all_results <- list()
  
  for (chr in chr) {
    if (verbose) {
      message("\n处理染色体: ", chr)
    }
    
    # 创建染色体区域
    chr_region <- paste0(chr, ":1-", bam_header[[1]]$targets[chr])
    
    # 处理该染色体
    chr_result <- pileup_from_bam(
      bam_file = bam_file,
      region = chr_region,
      verbose = FALSE,  # 子函数不显示详细信息
      ...
    )
    
    if (nrow(chr_result) > 0) {
      all_results[[chr]] <- chr_result
    }
    
    # 可选：每处理完一条染色体后保存中间结果
    if (exists("save_intermediate") && save_intermediate) {
      saveRDS(chr_result, paste0("pileup_", chr, ".rds"))
    }
  }
  
  # 合并所有结果
  if (length(all_results) == 0) {
    return(data.frame(
      seqnames = character(0),
      pos = integer(0),
      strand = character(0),
      nucleotide = character(0),
      count = integer(0)
    ))
  }
  
  final_result <- do.call(rbind, all_results)
  rownames(final_result) <- NULL
  
  if (verbose) {
    message("\n=== 全基因组pileup完成 ===")
    message("总行数: ", nrow(final_result))
    message("染色体数: ", length(unique(final_result$seqnames)))
  }
  
  return(final_result)
}


# 使用示例
example_usage <- function() {
  # 编译C++代码
  Rcpp::sourceCpp("pileup_fixed.cpp")

  # 加载封装函数
  source("pileup_wrapper.R")

  # 指定BAM文件
  bam_file <- "your_file.bam"

  # 示例1：处理特定区域
  cat("=== 示例1：处理特定区域 ===\n")
  result1 <- pileup_from_bam(
    bam_file = bam_file,
    region = "chr1:1000000-1100000",
    min_mapq = 20,
    min_base_quality = 20,
    distinguish_strands = TRUE,
    distinguish_nucleotides = TRUE,
    verbose = TRUE
  )

  print(head(result1, 20))
  cat("区域行数:", nrow(result1), "\n")

  # 示例2：处理多个区域
  cat("\n=== 示例2：处理多个区域 ===\n")
  result2 <- pileup_from_bam(
    bam_file = bam_file,
    region = c("chr1:1000000-1001000", "chr2:500000-510000"),
    verbose = TRUE
  )

  # 示例3：使用GRanges对象
  cat("\n=== 示例3：使用GRanges对象 ===\n")
  library(GenomicRanges)
  gr <- GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges(
      start = c(1000000, 2000000),
      end = c(1001000, 2001000)
    )
  )

  result3 <- pileup_from_bam(
    bam_file = bam_file,
    region = gr,
    verbose = TRUE
  )

  # 示例4：全基因组pileup（慎用，可能内存很大）
  cat("\n=== 示例4：全基因组pileup（小样本测试） ===\n")
  result4 <- pileup_whole_bam(
    bam_file = bam_file,
    region = "chr1:1-100000",  # 只处理小范围测试
    min_mapq = 20,
    verbose = TRUE
  )

  # 保存结果
  if (nrow(result4) > 0) {
    write.csv(result4, "whole_genome_pileup_small_test.csv", row.names = FALSE)
    cat("结果已保存到 whole_genome_pileup_small_test.csv\n")
  }
}

# 运行示例
# example_usage()

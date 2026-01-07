#' 并行bam2seqz_r_snps函数（按染色体并行）
#' @param normal_bam 正常样本BAM文件路径
#' @param tumor_bam 肿瘤样本BAM文件路径
#' @param bsgenome_pkg BSgenome包名（字符串），默认为"BSgenome.Hsapiens.UCSC.hg38"
#' @param min_depth 最小深度，默认为30
#' @param min_af 最小等位基因频率，默认为0.25
#' @param max_af 最大等位基因频率，默认为0.75
#' @param gc_window GC窗口大小，默认为50
#' @param output_file 输出文件路径
#' @param n_cores 并行核心数，默认为总核心数-1
#' @param chromosomes 要处理的染色体列表，默认为所有标准染色体
#' @export

parallel_bam2seqz_r_snps <- function(normal_bam, tumor_bam,
                                     bsgenome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
                                     min_depth = 30,
                                     min_af = 0.25,
                                     max_af = 0.75,
                                     gc_window = 50,
                                     output_file = "sample.snps.seqz.txt",
                                     n_cores = 8,
                                     chromosomes = NULL) {
  
  
  
  cat("开始参数检查...\n")
  
  # 1. 检查BAM文件
  if (!file.exists(normal_bam)) {
    stop(sprintf("正常样本BAM文件不存在: %s", normal_bam))
  }
  if (!file.exists(tumor_bam)) {
    stop(sprintf("肿瘤样本BAM文件不存在: %s", tumor_bam))
  }
  
  
  # 2. 检查BSgenome包参数
  if (!is.character(bsgenome_pkg) || length(bsgenome_pkg) != 1) {
    stop(sprintf("bsgenome_pkg必须是单个字符串，当前为: %s", 
                 paste(class(bsgenome_pkg), collapse = ", ")))
  }
  bsgenome_local <- get(bsgenome_pkg)
  
  if(!inherits(bsgenome_local,"BSgenome")){
    stop(sprintf("节点上无法加载BSgenome包: %s", bsgenome_pkg))
    
  }
  
  
  
  # 4. 检查输出文件路径
  if (!is.character(output_file) || length(output_file) != 1) {
    stop(sprintf("output_file必须是单个字符串，当前为: %s", 
                 paste(class(output_file), collapse = ", ")))
  }
  
  
  cat(sprintf("并行处理SNPs：%s vs %s\n", 
              basename(normal_bam), basename(tumor_bam)))
  cat(sprintf("使用BSgenome包：%s\n", bsgenome_pkg))
  
  # 1. 在主节点加载BSgenome包并获取染色体信息
  if (!require(bsgenome_pkg, character.only = TRUE)) {
    stop(sprintf("请安装BSgenome包: %s", bsgenome_pkg))
  }
  
  # 获取BSgenome对象
  bsgenome_obj <- get(bsgenome_pkg)
  
  # 2. 确定要处理的染色体
  if (is.null(chromosomes)) {
    # 获取所有染色体名
    all_chroms <- seqnames(bsgenome_obj)
    # 选择标准染色体（1-22, X, Y, M）
    chromosomes <- grep("^(chr)?([0-9]+|[XYM])$", all_chroms, value = TRUE, perl = TRUE)
  }
  
  cat(sprintf("处理 %d 条染色体：%s\n", 
              length(chromosomes),
              paste(head(chromosomes, 5), collapse = ", ")))
  
  # 3. 创建临时目录
  temp_dir <- tempfile("snps_chr_")
  dir.create(temp_dir, recursive = TRUE)
  cat(sprintf("临时目录：%s\n", temp_dir))
  
  # 4. 清空输出文件
  if (file.exists(output_file)) unlink(output_file)
  
  # 5. 并行处理
  cat("启动并行处理...\n")
  
  # 创建集群
  cl <- parallel::makeCluster(n_cores)
  
  # 在每个节点设置环境
  parallel::clusterEvalQ(cl, {
    # 加载bam2seqz_r_snps所在的包
    if (!require("HrdDetective", character.only = TRUE)) {
      # 如果HrdDetective包不存在，尝试其他可能的包名

    }
  })
  
  # 传递变量到集群（注意：不传递BSgenome对象，只传递包名）
  parallel::clusterExport(cl, 
                          varlist = c("normal_bam", "tumor_bam", "bsgenome_pkg",
                                      "min_depth", "min_af", "max_af", "gc_window",
                                      "temp_dir"),
                          envir = environment()
  )
  
  # 处理函数
  process_chromosome <- function(chr) {
    # 在每个节点上加载BSgenome包并获取对象
    if (!require(bsgenome_pkg, character.only = TRUE)) {
      stop(sprintf("节点上无法加载BSgenome包: %s", bsgenome_pkg))
    }
    bsgenome_local <- get(bsgenome_pkg)
    
    if(!inherits(bsgenome_local,"BSgenome")){
      stop(sprintf("节点上无法加载BSgenome包: %s", bsgenome_pkg))
      
    }
    
    # 输出文件路径
    chr_file <- file.path(temp_dir, paste0("snps_", chr, ".txt"))
    
    # 调用bam2seqz_r_snps函数
    tryCatch({
      bam2seqz_r_snps(
        normal_bam = normal_bam,
        tumor_bam = tumor_bam,
        bsgenome = bsgenome_local,
        min_depth = min_depth,
        min_af = min_af,
        max_af = max_af,
        gc_window = gc_window,
        output_file = chr_file,
        chr = chr
      )
    }, error = function(e) {
      cat(sprintf("染色体 %s 处理失败: %s\n", chr, e$message))
      return(NULL)
    })
    
    # 返回结果信息
    if (file.exists(chr_file) && file.size(chr_file) > 0) {
      lines <- readLines(chr_file, warn = FALSE)
      return(list(chr = chr, file = chr_file, lines = lines, success = TRUE))
    } else {
      return(list(chr = chr, file = NULL, lines = character(), success = FALSE))
    }
  }
  
  # 并行运行
  results <- parallel::parLapply(cl, chromosomes, process_chromosome)
  
  # 停止集群
  parallel::stopCluster(cl)
  
  # 6. 合并结果
  cat("\n合并结果...\n")
  
  # 筛选成功的染色体
  success_results <- results[!sapply(results, is.null)]
  
  if (length(success_results) == 0) {
    cat("警告：没有成功处理的染色体\n")
    unlink(temp_dir, recursive = TRUE)
    return(NULL)
  }
  
  # 按染色体顺序排序
  chr_order <- order(sapply(success_results, function(x) x$chr))
  success_results <- success_results[chr_order]
  
  # 合并文件
  total_lines <- 0
  for (i in seq_along(success_results)) {
    res <- success_results[[i]]
    if (res$success && length(res$lines) > 0) {
      # 写入文件
      write(res$lines, output_file, append = (i > 1))
      total_lines <- total_lines + length(res$lines)
      cat(sprintf("  染色体 %s: %d 行\n", res$chr, length(res$lines)))
    }
  }
  
  # 7. 清理临时目录
  unlink(temp_dir, recursive = TRUE)
  cat(sprintf("清理临时目录：%s\n", temp_dir))
  
  # 8. 输出统计信息
  if (file.exists(output_file) && file.size(output_file) > 0) {
    final_lines <- readLines(output_file, warn = FALSE)
    cat(sprintf("\n处理完成！总计 %d 行数据\n", length(final_lines)))
    cat(sprintf("输出文件：%s (%.2f MB)\n", 
                output_file, file.size(output_file) / 1024^2))
    
    # 显示文件预览
    if (length(final_lines) > 0) {
      cat("前3行预览：\n")
      for (i in 1:min(3, length(final_lines))) {
        cat(sprintf("  [%d] %s\n", i, substr(final_lines[i], 1, 80)))
      }
    }
    
    return(output_file)
  } else {
    cat("警告：输出文件为空\n")
    return(NULL)
  }
}

# 使用示例
if (FALSE) {
  # 注意：这里不再需要传递BSgenome对象，只需要包名
  library(HrdDetective);
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings);
  
  # 处理所有染色体
  output <- parallel_bam2seqz_r_snps(
    normal_bam = "target.normal.dedup.bam",
    tumor_bam = "target.tumor.dedup.bam",
    bsgenome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
    output_file = "output.snps.seqz.txt",
    n_cores = 8
  )
  
  # 只处理部分染色体
  output <- parallel_bam2seqz_r_snps(
    normal_bam = "chr1.normal.dedup.bam",
    tumor_bam = "chr1.tumor.dedup.bam",
    bsgenome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
    output_file = "output.snps.seqz.txt",
    n_cores = 4,
    chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5")
  )
  
  # 如果使用不同物种
  output <- parallel_bam2seqz_r_snps(
    normal_bam = "target.chr1.normal.dedup.bam",
    tumor_bam = "target.chr1.tumor.dedup.bam",
    bsgenome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
    output_file = "output.snps.seqz.txt",
    n_cores = 8
  )
}
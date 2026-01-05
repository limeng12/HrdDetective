#' 生成PNG图片的HRD分析HTML报告
#'
#' 这个版本先用你的原始函数生成PNG图片，然后在HTML中加载显示
#' @param seg.tab the segment file.
#' @param seqz_list the raw input.
#' @param celluPloidyMat cellularity and ploidy matrix.
#' @param best_ploidy the optimized ploidy.
#' @param best_cellularity the optimized cellularity.
#' @param output_dir_name the output directory.
#' @return HTML文件路径
#' @export
generate_png_hrd_report <- function(seg.tab,
                                    seqz_list,
                                    celluPloidyMat,
                                    best_ploidy,
                                    best_cellularity,
                                    output_dir_name = "./hrd_png_report") {
  
  # 创建输出目录
  dir.create(output_dir_name, showWarnings = FALSE, recursive = TRUE)
  
  # 保存原始数据
  write.table(seg.tab, file = paste0(output_dir_name, "/raw_segments.txt"),
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE);
  segments<-read.table(paste0(output_dir_name, "/raw_segments.txt"),
                       header = TRUE, as.is = TRUE, sep = "\t");
  scarHRD_input_fr<-data.frame(SampleID=rep(output_dir_name, nrow(segments) ),
                               Chromosome=str_c("chr", segments$chromosome),
                               Start_position=segments$start.pos,
                               End_position=segments$end.pos,
                               total_cn=segments$CNt,
                               A_cn=segments$A,
                               B_cn=segments$B,
                               ploidy=rep(best_ploidy, nrow(segments)) );
  scarHRD_input_fr$best_ploidy<-best_ploidy;
  write.table(scarHRD_input_fr, file=paste0(output_dir_name, "/scarHRD_input_test.tsv"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE);
  cat(paste0("best_cellularity: ", best_cellularity, "\n") ,file=paste0(output_dir_name, "/HRD_re.txt") );
  calculated_sample_ploidy<-best_cellularity*sum(scarHRD_input_fr$total_cn*(scarHRD_input_fr$End_position-scarHRD_input_fr$Start_position) )/
    (sum(scarHRD_input_fr$End_position-scarHRD_input_fr$Start_position))+(1-best_cellularity)*2
  cat(paste0("calculated sample ploidy: ",calculated_sample_ploidy, "\n") ,file=paste0(output_dir_name, "/HRD_re.txt"), append = TRUE );
  calculated_tumor_ploidy<-sum(scarHRD_input_fr$total_cn*(scarHRD_input_fr$End_position-scarHRD_input_fr$Start_position) )/
    (sum(scarHRD_input_fr$End_position-scarHRD_input_fr$Start_position))
  cat(paste0("calculated tumor ploidy: ",calculated_tumor_ploidy, "\n") ,file=paste0(output_dir_name, "/HRD_re.txt"), append = TRUE );
  cat(paste0("best_ploidy: ",best_ploidy, "\n") ,file=paste0(output_dir_name, "/HRD_re.txt"), append = TRUE );
  
  
  cat("正在生成PNG图片...\n")
  
  # 1. 创建图片目录
  img_dir <- paste0(output_dir_name, "/images")
  dir.create(img_dir, showWarnings = FALSE)
  
  # 2. 使用你的原始绘图函数，但输出PNG格式
  
  # 2.1 生成水平整合视图PNG
  cat("生成水平整合视图PNG...\n")
  generate_horizontal_png(seg.tab, img_dir)
  
  # 2.2 生成模型拟合视图PNG
  cat("生成模型拟合视图PNG...\n")
  generate_model_png(celluPloidyMat, best_ploidy, best_cellularity, img_dir)
  
  # 2.3 生成染色体详细视图PNG（每个染色体一个）
  cat("生成染色体详细视图PNG...\n")
  chr_images <- generate_chromosome_pngs(seqz_list, seg.tab, img_dir)
  
  # 3. 创建HTML报告
  cat("创建HTML文档...\n")
  html_file <- create_png_html_report(
    seg.tab = seg.tab,
    best_ploidy = best_ploidy,
    best_cellularity = best_cellularity,
    output_dir_name = output_dir_name,
    img_dir = img_dir,
    chr_images = chr_images
  )
  
  cat(paste0("✅ 报告已生成: ", html_file, "\n"))
  cat(paste0("📁 图片保存在: ", img_dir, "\n"))
  cat(paste0("🌐 请在浏览器中打开查看: file://", normalizePath(html_file), "\n"))
  
  return(html_file)
}

#' 生成水平整合视图PNG
generate_horizontal_png <- function(seg.tab, img_dir) {
  # 使用你的draw_seg_hori函数逻辑，但输出PNG
  library(ggplot2)
  
  # 复刻你的draw_seg_hori函数逻辑
  seg.tab$chromosome <- factor(seg.tab$chromosome, 
                               levels = unique(seg.tab$chromosome))
  
  p <- ggplot(data = seg.tab) +
    geom_segment(aes(x = start.pos, y = B - 0.05, 
                     xend = end.pos, yend = B - 0.05),
                 colour = "blue", linewidth = 1) +
    geom_segment(aes(x = start.pos, y = A + 0.05, 
                     xend = end.pos, yend = A + 0.05),
                 colour = "red", linewidth = 1) +
    scale_y_continuous(breaks = 0:7, labels = 0:7) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f0f0f0"),
      strip.text = element_text(face = "bold", size = 10)
    ) +
    labs(
      x = "基因组位置",
      y = "拷贝数",
      title = "染色体拷贝数水平整合视图"
    ) +
    facet_wrap(vars(chromosome), scales = "free_x", nrow = 5)
  
  # 保存PNG
  ggsave(
    filename = paste0(img_dir, "/horizontal_view.png"),
    plot = p,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
  )
}

#' 生成模型拟合视图PNG
generate_model_png <- function(celluPloidyMat, best_ploidy, best_cellularity, img_dir) {
  library(ggplot2)
  
  # 复刻你的draw_model_fit函数逻辑
  p1 <- ggplot() +
    geom_contour(data = celluPloidyMat, 
                 aes(ploidy, cellularity, z = 1 * fit_v), 
                 bins = 500, color = "gray50") +
    geom_point(aes(x = best_ploidy, y = best_cellularity), 
               color = "red", size = 3) +
    theme_minimal() +
    labs(
      title = paste0("模型拟合等值线图 (倍体数: ", best_ploidy, 
                     ", 细胞纯度: ", best_cellularity, ")"),
      x = "倍体数",
      y = "细胞纯度"
    )
  
  p2 <- ggplot(celluPloidyMat, 
               aes(ploidy, cellularity, 
                   fill = -1 * log(fit_v - min(fit_v) + 1))) +
    geom_tile() +
    geom_point(aes(x = best_ploidy, y = best_cellularity),
               color = "red", size = 2, shape = 21, fill = "white") +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_minimal() +
    labs(
      title = paste0("模型拟合热图 (倍体数: ", best_ploidy, 
                     ", 细胞纯度: ", best_cellularity, ")"),
      x = "倍体数",
      y = "细胞纯度",
      fill = "拟合值 = -log(llv-min(llv))"
    ) +
    guides(fill = guide_colorbar(title.position = "top", 
                                 title.hjust = 0.5))
  
  # 保存第一个图
  ggsave(
    filename = paste0(img_dir, "/model_fit_contour.png"),
    plot = p1,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  # 保存第二个图
  ggsave(
    filename = paste0(img_dir, "/model_fit_heatmap.png"),
    plot = p2,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

#' 生成染色体详细视图PNG
generate_chromosome_pngs <- function(seqz_list, seg.tab, img_dir) {
  all_chr <- names(seqz_list)
  chr_images <- list()
  
  # 为每个染色体生成PNG图片
  for (chr in all_chr[1:min(100, length(all_chr))]) {  # 限制前5个
    chr_img_file <- generate_single_chr_png(seqz_list, seg.tab, chr, img_dir)
    if (!is.null(chr_img_file)) {
      chr_images[[as.character(chr)]] <- chr_img_file
    }
  }
  
  return(chr_images)
}

#' 生成单个染色体PNG
#' 生成单个染色体PNG
generate_single_chr_png <- function(seqz_list, seg.tab, chr_num, img_dir) {
  if (!(chr_num %in% names(seqz_list))) return(NULL)
  
  df_data_chr <- seqz_list[[as.character(chr_num)]]
  seg_one <- seg.tab[seg.tab$chromosome == chr_num, ]
  
  if (nrow(seg_one) == 0 || nrow(df_data_chr) == 0) return(NULL)
  
  # 修复：使用 as.integer() 将序列转换为整数
  # 限制数据量，提高绘图速度
  if (nrow(df_data_chr) > 2000) {
    # 创建整数索引序列
    indices <- seq(1, nrow(df_data_chr), length.out = 2000)
    # 转换为整数（四舍五入取整）
    indices <- as.integer(round(indices))
    # 确保不重复且不超出范围
    indices <- unique(pmax(1, pmin(nrow(df_data_chr), indices)))
    df_data_chr <- df_data_chr[indices, ]
  }
  
  # 其余代码保持不变...
  # 创建三图垂直排列
  library(ggplot2)
  library(patchwork)
  
  # 图1: B等位基因频率
  p1 <- ggplot() +
    geom_point(aes(x = df_data_chr$position, y = df_data_chr$Bf), 
               alpha = 0.6, size = 0.5, color = "gray50") +
    geom_segment(data = seg_one,
                 aes(x = start.pos, y = Bf, 
                     xend = end.pos, yend = Bf),
                 colour = "red", linewidth = 1) +
    theme_minimal() +
    labs(
      x = "",
      y = "B等位基因频率",
      title = paste("染色体", chr_num)
    ) +
    theme(
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # 图2: 深度比
  p2 <- ggplot() +
    geom_point(aes(x = df_data_chr$position, y = df_data_chr$adjusted.ratio), 
               alpha = 0.6, size = 0.5, color = "lightblue") +
    geom_segment(data = seg_one,
                 aes(x = start.pos, y = depth.ratio, 
                     xend = end.pos, yend = depth.ratio),
                 colour = "red", linewidth = 1) +
    ylim(0, 3) +
    theme_minimal() +
    labs(
      x = "",
      y = "深度比"
    ) +
    theme(axis.text.x = element_blank())
  
  # 图3: 拷贝数
  p3 <- ggplot() +
    geom_segment(data = seg_one,
                 aes(x = start.pos, y = B - 0.05, 
                     xend = end.pos, yend = B - 0.05),
                 colour = "blue", linewidth = 1) +
    geom_segment(data = seg_one,
                 aes(x = start.pos, y = A + 0.05, 
                     xend = end.pos, yend = A + 0.05),
                 colour = "red", linewidth = 1) +
    scale_y_continuous(breaks = 0:7, labels = 0:7) +
    theme_minimal() +
    labs(
      x = "基因组位置",
      y = "拷贝数"
    )
  
  # 组合三图
  combined_plot <- p1 / p2 / p3 + 
    plot_layout(heights = c(1, 1, 1.2))
  
  # 保存PNG
  file_name <- paste0("chr_", chr_num, "_detailed.png")
  file_path <- paste0(img_dir, "/", file_name)
  
  ggsave(
    filename = file_path,
    plot = combined_plot,
    width = 12,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  return(file_name)
}



#' 创建基于PNG的HTML报告（修正转义问题）
#' 创建基于PNG的HTML报告（英文版）
create_png_html_report <- function(seg.tab, best_ploidy, best_cellularity,
                                   output_dir_name, img_dir, chr_images) {
  
  # 基本统计
  stats <- list(
    n_segments = nrow(seg.tab),
    total_len = sum(seg.tab$end.pos - seg.tab$start.pos),
    mean_cn_a = round(mean(seg.tab$A, na.rm = TRUE), 2),
    mean_cn_b = round(mean(seg.tab$B, na.rm = TRUE), 2),
    ploidy = best_ploidy,
    cellularity = best_cellularity,
    chromosomes = length(unique(seg.tab$chromosome))
  )
  
  # 生成表格行HTML
  table_rows <- generate_html_table_rows(seg.tab)
  
  # 生成染色体选择选项
  chr_options <- ""
  if (length(chr_images) > 0) {
    chr_options <- paste(
      sapply(names(chr_images), function(chr) {
        sprintf('<option value="%s">Chromosome %s</option>', chr, chr)
      }),
      collapse = "\n"
    )
  }
  
  # 构建HTML（英文版）
  html_content <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HRD Analysis Report (PNG Version)</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    <style>
        :root {
            --primary: #667eea;
            --secondary: #764ba2;
            --light: #f8f9fa;
            --dark: #343a40;
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
            background: linear-gradient(135deg, #667eea15 0%, #764ba215 100%);
            color: var(--dark);
            line-height: 1.6;
            min-height: 100vh;
        }
        
        .app-container {
            display: flex;
            min-height: 100vh;
        }
        
        /* Sidebar */
        .sidebar {
            width: 280px;
            background: rgba(255, 255, 255, 0.98);
            box-shadow: 2px 0 20px rgba(0,0,0,0.1);
            position: fixed;
            height: 100vh;
            overflow-y: auto;
            z-index: 1000;
        }
        
        .sidebar-header {
            padding: 25px;
            background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
            color: white;
            text-align: center;
        }
        
        .sidebar-header h1 {
            font-size: 1.6em;
            margin-bottom: 8px;
        }
        
        .sidebar-nav {
            padding: 20px 0;
        }
        
        .nav-item {
            display: flex;
            align-items: center;
            padding: 16px 25px;
            color: #555;
            text-decoration: none;
            border-left: 4px solid transparent;
            transition: all 0.3s ease;
            cursor: pointer;
        }
        
        .nav-item:hover {
            background: rgba(102, 126, 234, 0.08);
            color: var(--primary);
        }
        
        .nav-item.active {
            background: rgba(102, 126, 234, 0.12);
            color: var(--primary);
            border-left-color: var(--primary);
            font-weight: 600;
        }
        
        .nav-item i {
            width: 24px;
            font-size: 1.2em;
            margin-right: 12px;
        }
        
        .sidebar-stats {
            padding: 20px;
            background: var(--light);
            margin: 20px;
            border-radius: 10px;
            border: 1px solid #e9ecef;
        }
        
        .stat-item {
            margin-bottom: 12px;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        
        .stat-label {
            font-size: 0.9em;
            color: #666;
        }
        
        .stat-value {
            font-weight: 600;
            color: var(--primary);
        }
        
        /* Main Content */
        .main-content {
            flex: 1;
            margin-left: 280px;
            padding: 30px;
        }
        
        .section {
            background: white;
            border-radius: 12px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.08);
            animation: fadeIn 0.5s ease-out;
        }
        
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(10px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .section-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 25px;
            padding-bottom: 15px;
            border-bottom: 2px solid var(--light);
        }
        
        .section-title {
            color: var(--primary);
            font-size: 1.6em;
            display: flex;
            align-items: center;
            gap: 10px;
        }
        
        .chart-container {
            width: 100%;
            background: white;
            border-radius: 8px;
            border: 1px solid #e9ecef;
            overflow: hidden;
            margin: 20px 0;
            text-align: center;
            padding: 20px;
        }
        
        .chart-image {
            max-width: 100%;
            height: auto;
            border-radius: 6px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            transition: transform 0.3s ease;
        }
        
        .chart-image:hover {
            transform: scale(1.01);
        }
        
        .chart-controls {
            display: flex;
            gap: 10px;
            justify-content: center;
            margin-top: 20px;
        }
        
        .btn {
            padding: 10px 20px;
            background: var(--primary);
            color: white;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.95em;
            display: inline-flex;
            align-items: center;
            gap: 8px;
            transition: all 0.3s;
        }
        
        .btn:hover {
            background: #5a67d8;
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.3);
        }
        
        .chr-selector {
            display: flex;
            gap: 15px;
            align-items: center;
            margin: 20px 0;
            padding: 20px;
            background: var(--light);
            border-radius: 8px;
            justify-content: center;
        }
        
        .chr-select {
            padding: 10px 20px;
            border: 2px solid #dee2e6;
            border-radius: 6px;
            background: white;
            font-size: 1em;
            min-width: 200px;
            cursor: pointer;
        }
        
        .image-caption {
            margin-top: 15px;
            color: #666;
            font-size: 0.95em;
            text-align: center;
        }
        
        .data-table {
            width: 100%;
            margin-top: 20px;
            border-collapse: collapse;
            font-size: 0.95em;
        }
        
        .data-table th {
            background: var(--light);
            padding: 12px 15px;
            text-align: left;
            font-weight: 600;
            color: var(--dark);
            border-bottom: 2px solid #dee2e6;
        }
        
        .data-table td {
            padding: 10px 15px;
            border-bottom: 1px solid #e9ecef;
        }
        
        .data-table tr:hover {
            background: rgba(102, 126, 234, 0.05);
        }
        
        .file-list {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        
        .file-card {
            background: var(--light);
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid var(--primary);
        }
        
        @media (max-width: 1024px) {
            .sidebar {
                width: 80px;
            }
            
            .sidebar-header h1,
            .nav-item span,
            .stat-label {
                display: none;
            }
            
            .nav-item i {
                margin-right: 0;
                font-size: 1.4em;
            }
            
            .main-content {
                margin-left: 80px;
                padding: 20px;
            }
        }
    </style>
</head>
<body>
    <div class="app-container">
        <!-- Sidebar -->
        <div class="sidebar">
            <div class="sidebar-header">
                <h1><i class="fas fa-dna"></i> HRD Analysis</h1>
                <p>PNG Report Version</p>
            </div>
            
            <div class="sidebar-nav">
                <div class="nav-item active" onclick="showSection(\'overview\')">
                    <i class="fas fa-home"></i>
                    <span>Overview</span>
                </div>
                <div class="nav-item" onclick="showSection(\'chromosomes\')">
                    <i class="fas fa-chart-line"></i>
                    <span>Chromosome View</span>
                </div>
                <div class="nav-item" onclick="showSection(\'horizontal\')">
                    <i class="fas fa-braille"></i>
                    <span>Horizontal View</span>
                </div>
                <div class="nav-item" onclick="showSection(\'model\')">
                    <i class="fas fa-chart-area"></i>
                    <span>Model Fitting</span>
                </div>
                <div class="nav-item" onclick="showSection(\'data\')">
                    <i class="fas fa-table"></i>
                    <span>Data Table</span>
                </div>
            </div>
            
            <div class="sidebar-stats">
                <h4 style="margin-bottom: 15px; color: var(--dark);">
                    <i class="fas fa-chart-bar"></i> Key Statistics
                </h4>
                <div class="stat-item">
                    <span class="stat-label">Chromosomes</span>
                    <span class="stat-value">', stats$chromosomes, '</span>
                </div>
                <div class="stat-item">
                    <span class="stat-label">Segments</span>
                    <span class="stat-value">', stats$n_segments, '</span>
                </div>
                <div class="stat-item">
                    <span class="stat-label">Total Length</span>
                    <span class="stat-value">', round(stats$total_len / 1e6, 1), ' Mb</span>
                </div>
                <div class="stat-item">
                    <span class="stat-label">Ploidy</span>
                    <span class="stat-value">', stats$ploidy, '</span>
                </div>
            </div>
        </div>
        
        <!-- Main Content -->
        <div class="main-content">
            <!-- Overview -->
            <div id="overview" class="section">
                <h2 class="section-title">
                    <i class="fas fa-home"></i> Report Overview
                </h2>
                
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin: 25px 0;">
                    <div style="background: linear-gradient(135deg, #667eea10 0%, #764ba210 100%); padding: 20px; border-radius: 10px;">
                        <h3 style="color: var(--dark); margin-bottom: 15px; font-size: 1.1em;">
                            <i class="fas fa-info-circle"></i> Analysis Summary
                        </h3>
                        <p style="color: #666; font-size: 0.95em; line-height: 1.7;">
                            This report visualizes genomic data generated by HrdDetective tool. It includes detailed chromosome views, horizontal integrated views, and model fitting visualizations.
                        </p>
                    </div>
                    
                    <div style="background: linear-gradient(135deg, #28a74510 0%, #20c99710 100%); padding: 20px; border-radius: 10px;">
                        <h3 style="color: var(--dark); margin-bottom: 15px; font-size: 1.1em;">
                            <i class="fas fa-cogs"></i> Optimized Parameters
                        </h3>
                        <div style="font-size: 0.95em;">
                            <div style="margin-bottom: 8px;">
                                <strong>Best Ploidy:</strong> <span style="color: var(--primary);">', stats$ploidy, '</span>
                            </div>
                            <div>
                                <strong>Cellularity:</strong> <span style="color: var(--primary);">', stats$cellularity, '</span>
                            </div>
                        </div>
                    </div>
                    
                    <div style="background: linear-gradient(135deg, #ffc10710 0%, #fd7e1410 100%); padding: 20px; border-radius: 10px;">
                        <h3 style="color: var(--dark); margin-bottom: 15px; font-size: 1.1em;">
                            <i class="fas fa-chart-line"></i> Usage Tips
                        </h3>
                        <ul style="font-size: 0.9em; color: #666; padding-left: 20px;">
                            <li style="margin-bottom: 5px;">Click left navigation to switch views</li>
                            <li style="margin-bottom: 5px;">Right-click images to save locally</li>
                            <li>Use download buttons to get raw data</li>
                        </ul>
                    </div>
                </div>
            </div>
            
            <!-- Chromosome View -->
            <div id="chromosomes" class="section" style="display: none;">
                <div class="section-header">
                    <h2 class="section-title">
                        <i class="fas fa-chart-line"></i> Chromosome Detailed View
                    </h2>
                </div>
                
                <div class="chr-selector">
                    <label for="chrSelect" style="font-weight: 600; color: var(--dark);">
                        Select Chromosome:
                    </label>
                    <select id="chrSelect" class="chr-select" onchange="loadChromosomeImage(this.value)">
                        <option value="">-- Select Chromosome --</option>
                        ', chr_options, '
                    </select>
                </div>
                
                <div class="chart-container">
                    <div id="chrImageContainer">
                        <p style="color: #666; font-size: 1.1em; padding: 40px;">
                            👈 Please select a chromosome from above to view detailed charts
                        </p>
                    </div>
                </div>
                
                <div class="image-caption">
                    <p>Note: Each chromosome view contains three subplots - B allele frequency, depth ratio, and A/B copy number distribution</p>
                </div>
            </div>
            
            <!-- Horizontal View -->
            <div id="horizontal" class="section" style="display: none;">
                <div class="section-header">
                    <h2 class="section-title">
                        <i class="fas fa-braille"></i> Horizontal Integrated View
                    </h2>
                    <button class="btn" onclick="downloadImage(\'images/horizontal_view.png\', \'horizontal_view.png\')">
                        <i class="fas fa-download"></i> Download Image
                    </button>
                </div>
                
                <div class="chart-container">
                    <img src="images/horizontal_view.png" alt="Horizontal Integrated View" class="chart-image" id="horizontalImage">
                </div>
                
                <div class="image-caption">
                    <p>Figure 1: Horizontal integrated view of A/B copy numbers across all chromosomes. Blue lines represent B allele copy numbers, red lines represent A allele copy numbers.</p>
                </div>
            </div>
            
            <!-- Model Fitting -->
            <div id="model" class="section" style="display: none;">
                <div class="section-header">
                    <h2 class="section-title">
                        <i class="fas fa-chart-area"></i> Model Fitting View
                    </h2>
                    <div class="chart-controls">
                        <button class="btn" onclick="downloadImage(\'images/model_fit_contour.png\', \'model_fit_contour.png\')">
                            <i class="fas fa-download"></i> Download Contour
                        </button>
                        <button class="btn" onclick="downloadImage(\'images/model_fit_heatmap.png\', \'model_fit_heatmap.png\')">
                            <i class="fas fa-download"></i> Download Heatmap
                        </button>
                    </div>
                </div>
                
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 30px;">
                    <div class="chart-container">
                        <img src="images/model_fit_contour.png" alt="Model Fitting Contour Plot" class="chart-image">
                        <div class="image-caption">
                            <p>Figure 2: Cellularity and ploidy fitting contour plot</p>
                        </div>
                    </div>
                    
                    <div class="chart-container">
                        <img src="images/model_fit_heatmap.png" alt="Model Fitting Heatmap" class="chart-image">
                        <div class="image-caption">
                            <p>Figure 3: Cellularity and ploidy fitting heatmap</p>
                        </div>
                    </div>
                </div>
            </div>
            
            <!-- Data Table -->
            <div id="data" class="section" style="display: none;">
                <div class="section-header">
                    <h2 class="section-title">
                        <i class="fas fa-table"></i> Data Table
                    </h2>
                    <button class="btn" onclick="downloadData()">
                        <i class="fas fa-file-csv"></i> Download Full Data
                    </button>
                </div>
                
                <div style="overflow-x: auto;">
                    <table class="data-table">
                        <thead>
                            <tr>
                                <th>Chromosome</th>
                                <th>Start Position</th>
                                <th>End Position</th>
                                <th>Total CN</th>
                                <th>A Copy</th>
                                <th>B Copy</th>
                                <th>B Allele Frequency</th>
                            </tr>
                        </thead>
                        <tbody>', table_rows, '</tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>
    
    <script>
        // Show specified section
        function showSection(sectionId) {
            // Hide all sections
            var sections = document.querySelectorAll(".section");
            sections.forEach(function(section) {
                section.style.display = "none";
            });
            
            // Show selected section
            document.getElementById(sectionId).style.display = "block";
            
            // Update navigation item status
            var navItems = document.querySelectorAll(".nav-item");
            navItems.forEach(function(item) {
                item.classList.remove("active");
            });
            
            // Find corresponding navigation item and activate it
            var sectionMap = {
                "overview": 0,
                "chromosomes": 1,
                "horizontal": 2,
                "model": 3,
                "data": 4
            };
            
            if (sectionMap[sectionId] !== undefined) {
                navItems[sectionMap[sectionId]].classList.add("active");
            }
            
            // If switching to chromosome view, show first chromosome by default
            if (sectionId === "chromosomes") {
                var chrSelect = document.getElementById("chrSelect");
                if (chrSelect && chrSelect.options.length > 1) {
                    chrSelect.value = chrSelect.options[1].value;
                    loadChromosomeImage(chrSelect.value);
                }
            }
        }
        
        // Load chromosome image
        function loadChromosomeImage(chr) {
            if (!chr) return;
            
            var container = document.getElementById("chrImageContainer");
            var imageUrl = "images/chr_" + chr + "_detailed.png";
            
            container.innerHTML = \'<img src="\' + imageUrl + \'" alt="Chromosome \' + chr + \' Detailed View" class="chart-image">\' +
                                \'<div class="image-caption" style="margin-top: 20px;">\' +
                                \'<p>Chromosome \' + chr + \' Detailed View</p></div>\' +
                                \'<div class="chart-controls">\' +
                                \'<button class="btn" onclick="downloadImage(\\\'\' + imageUrl + \'\\\', \\\'chr_\' + chr + \'_detailed.png\\\')">\' +
                                \'<i class="fas fa-download"></i> Download Image</button></div>\';
            
            // Set image loading failure handling
            var img = container.querySelector("img");
            if (img) {
                img.onerror = function() {
                    this.src = "data:image/svg+xml;utf8,<svg xmlns=\'http://www.w3.org/2000/svg\' width=\'400\' height=\'300\'><text x=\'50%\' y=\'50%\' text-anchor=\'middle\' fill=\'%23666\'>Image Load Failed</text></svg>";
                };
            }
        }
        
        // Download image
        function downloadImage(imagePath, filename) {
            var link = document.createElement("a");
            link.href = imagePath;
            link.download = filename;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }
        
        // Download data
        function downloadData() {
            window.open("raw_segments.txt", "_blank");
        }
        
        // Show overview by default after page loads
        document.addEventListener("DOMContentLoaded", function() {
            showSection("overview");
        });
    </script>
</body>
</html>')
  
  # 写入HTML文件
  html_file <- paste0(output_dir_name, "/hrd_png_report.html")
  writeLines(html_content, html_file, useBytes = TRUE)
  
  return(html_file)
}


#' 生成表格行HTML
generate_html_table_rows <- function(seg.tab) {
  n_display <- min(20, nrow(seg.tab))
  rows <- ""
  
  for (i in 1:n_display) {
    row <- seg.tab[i, ]
    rows <- paste0(rows, sprintf('
                        <tr>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td style="text-align: center;">%.2f</td>
                            <td style="text-align: center;">%.2f</td>
                            <td style="text-align: center;">%.2f</td>
                            <td style="text-align: center;">%.3f</td>
                        </tr>',
                                 row$chromosome,
                                 format(row$start.pos, big.mark = ","),
                                 format(row$end.pos, big.mark = ","),
                                 row$CNt,
                                 row$A,
                                 row$B,
                                 row$Bf))
  }
  
  if (nrow(seg.tab) > 20) {
    rows <- paste0(rows, sprintf('
                        <tr>
                            <td colspan="7" style="text-align: center; padding: 20px; color: #666; font-style: italic; background: #f8f9fa;">
                                还有 %d 行数据未显示... <a href="raw_segments.txt" style="color: #667eea; text-decoration: none; margin-left: 10px;">查看完整数据</a>
                            </td>
                        </tr>',
                                 nrow(seg.tab) - 20))
  }
  
  return(rows)
}

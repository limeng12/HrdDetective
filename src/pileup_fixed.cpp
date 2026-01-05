// pileup_fixed.cpp
#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <cctype>

using namespace Rcpp;

// 结构体存储单个位置信息
struct PositionInfo {
  int a_count = 0;
  int c_count = 0;
  int g_count = 0;
  int t_count = 0;
  int n_count = 0;
  int del_count = 0;
  int ins_count = 0;
  int forward_count = 0;
  int reverse_count = 0;
  int total_depth = 0;
  
  void reset() {
    a_count = c_count = g_count = t_count = n_count = del_count = ins_count = 0;
    forward_count = reverse_count = total_depth = 0;
  }
  
  void add_base(char base, bool is_reverse) {
    switch(base) {
    case 'A': case 'a': a_count++; break;
    case 'C': case 'c': c_count++; break;
    case 'G': case 'g': g_count++; break;
    case 'T': case 't': t_count++; break;
    case 'N': case 'n': n_count++; break;
    case '-': del_count++; break;
    case '+': ins_count++; break;
    default: n_count++; break;
    }
    total_depth++;
    if (is_reverse) {
      reverse_count++;
    } else {
      forward_count++;
    }
  }
};

// 辅助函数：重复字符串n次
CharacterVector rep_string(const std::string& str, int n) {
  CharacterVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = str;
  }
  return result;
}

// 辅助函数：生成序列
IntegerVector seq_int(int start, int end) {
  int length = end - start + 1;
  IntegerVector result(length);
  for (int i = 0; i < length; ++i) {
    result[i] = start + i;
  }
  return result;
}

// 解析 CIGAR 字符串
std::vector<std::pair<int, char>> parse_cigar_string(const std::string& cigar_str) {
  std::vector<std::pair<int, char>> result;
  std::string num_str;
  for (char c : cigar_str) {
    if (std::isdigit(c)) {
      num_str += c;
    } else if (c == 'M' || c == 'I' || c == 'D' || c == 'N' || c == 'S' || c == 'H' || c == 'P' || c == '=' || c == 'X') {
      if (!num_str.empty()) {
        int length = std::stoi(num_str);
        result.push_back(std::make_pair(length, c));
        num_str.clear();
      }
    }
  }
  return result;
}

// 返回 (ref_position, query_index) 的配对
std::vector<std::pair<int, int>> get_ref_query_pairs(int start_pos, const std::string& cigar_str) {
  std::vector<std::pair<int, int>> pairs;
  auto cigar_ops = parse_cigar_string(cigar_str);
  int ref_pos = start_pos;
  int query_pos = 0;
  
  for (const auto& op : cigar_ops) {
    int length = op.first;
    char type = op.second;
    
    switch(type) {
    case 'M':
    case '=':
    case 'X':
      for (int i = 0; i < length; ++i) {
        pairs.push_back(std::make_pair(ref_pos + i, query_pos + i));
      }
      ref_pos += length;
      query_pos += length;
      break;
    case 'I':
    case 'S':
      query_pos += length;
      break;
    case 'D':
    case 'N':
      ref_pos += length;
      break;
    case 'H':
    case 'P':
      break;
    }
  }
  return pairs;
}

// Helper: parse CIGAR string to get reference span
int cigar_ref_length(const std::string& cigar) {
  int ref_len = 0;
  std::string number = "";
  for (char c : cigar) {
    if (std::isdigit(c)) {
      number += c;
    } else {
      if (!number.empty()) {
        int len = std::stoi(number);
        if (c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X') {
          ref_len += len;
        }
        number = "";
      }
    }
  }
  return ref_len;
}


// [[Rcpp::export]]
List pileup_manual_vectors(
    CharacterVector seqnames,
    IntegerVector positions,
    CharacterVector cigars,
    CharacterVector sequences,
    CharacterVector qualities,
    IntegerVector mapqs,
    IntegerVector flags,
    IntegerVector qwidths,
    std::string seqname,
    int start,
    int end,
    int min_mapq = 0,
    int min_base_quality = 0,
    int min_depth = 0,
    float min_af=0.2,
    int max_depth = 1000000,
    bool distinguish_strands = false,
    bool distinguish_nucleotides = true,
    bool ignore_query_Ns = true,
    bool include_deletions = false,
    bool include_insertions = false
) {
  Rcpp::Rcout << "\n=== 开始pileup_manual_vectors调试 ===\n";
  Rcpp::Rcout << "参数:\n";
  Rcpp::Rcout << "  seqname: " << seqname << "\n";
  Rcpp::Rcout << "  start: " << start << "\n";
  Rcpp::Rcout << "  end: " << end << "\n";
  Rcpp::Rcout << "  distinguish_strands: " << distinguish_strands << "\n";
  Rcpp::Rcout << "  distinguish_nucleotides: " << distinguish_nucleotides << "\n";
  
  try {
    // Validate inputs
    if (start > end) {
      Rcpp::Rcout << "错误: start > end\n";
      Rcpp::stop("'start' must be <= 'end'");
    }
    
    if (positions.size() == 0) {
      Rcpp::Rcout << "警告: 没有比对记录，返回空DataFrame\n";
      return List::create(
        Named("seqnames") = CharacterVector(0),
        Named("pos") = IntegerVector(0),
        Named("strand") = CharacterVector(0),
        Named("nucleotide") = CharacterVector(0),
        Named("count") = IntegerVector(0)
      );
    }
    
    int n = positions.size();
    Rcpp::Rcout << "找到 " << n << " 条比对记录\n";
    
    // 使用 map 存储每个位置的信息
    // 外层map: position -> 内层map
    // 内层map: (strand, nucleotide) -> count
    std::map<int, std::map<std::pair<std::string, std::string>, int>> pileup_map;
    
    int processed_alignments = 0;
    int valid_alignments = 0;
    
    for (int i = 0; i < n; ++i) {
      processed_alignments++;
      
      if (positions[i] == NA_INTEGER) {
        continue;
      }
      
      if (mapqs[i] < min_mapq) {
        continue;
      }
      
      if (flags[i] & 0x4) {
        continue; // unmapped
      }
      
      valid_alignments++;
      
      int pos = positions[i];
      std::string cigar = as<std::string>(cigars[i]);
      std::string seq = as<std::string>(sequences[i]);
      std::string qual = as<std::string>(qualities[i]);
      
      // 判断链方向
      bool is_reverse = flags[i] & 0x10;  // 检查是否为反向互补链
      std::string strand = distinguish_strands ? (is_reverse ? "-" : "+") : "*";
      
      if (i < 3) {
        Rcpp::Rcout << "处理记录 " << i << ":\n";
        Rcpp::Rcout << "  起始位置: " << pos << "\n";
        Rcpp::Rcout << "  CIGAR: " << cigar << "\n";
        Rcpp::Rcout << "  链方向: " << strand << "\n";
      }
      
      int read_pos = 0;
      int ref_pos = pos;
      
      size_t j = 0;
      std::string num = "";
      while (j < cigar.size()) {
        char c = cigar[j];
        if (std::isdigit(c)) {
          num += c;
        } else {
          if (!num.empty()) {
            int len = std::stoi(num);
            num = "";
            
            if (c == 'M' || c == '=' || c == 'X') {
              // Match/mismatch
              for (int k = 0; k < len; ++k) {
                if (read_pos < (int)seq.size() && ref_pos >= start && ref_pos <= end) {
                  char base = std::toupper(seq[read_pos]);
                  
                  // 确定核苷酸
                  std::string nucleotide;
                  if (distinguish_nucleotides) {
                    nucleotide = std::string(1, base);
                  } else {
                    nucleotide = "*";
                  }
                  
                  // 检查N过滤
                  bool pass_N_filter = !ignore_query_Ns || base != 'N';
                  
                  // 检查质量过滤
                  bool pass_qual_filter = true;
                  if (min_base_quality > 0 && (int)qual.size() > read_pos) {
                    int qual_score = (unsigned char)qual[read_pos] - 33;
                    pass_qual_filter = qual_score >= min_base_quality;
                  }
                  
                  if (pass_N_filter && pass_qual_filter) {
                    std::pair<std::string, std::string> key = std::make_pair(strand, nucleotide);
                    pileup_map[ref_pos][key]++;
                    
                    if (i < 2 && k < 3) {
                      Rcpp::Rcout << "    计数: pos=" << ref_pos 
                                  << ", base=" << nucleotide 
                                  << ", strand=" << strand << "\n";
                    }
                  }
                }
                read_pos++;
                ref_pos++;
              }
            } else if (c == 'I') {
              // 插入
              if (include_insertions) {
                // 插入标记
                std::string nucleotide = distinguish_nucleotides ? "+" : "*";
                for (int k = 0; k < len; ++k) {
                  if (ref_pos >= start && ref_pos <= end) {
                    std::pair<std::string, std::string> key = std::make_pair(strand, nucleotide);
                    pileup_map[ref_pos][key]++;
                  }
                  if (read_pos < (int)seq.size()) {
                    read_pos++;
                  }
                }
              } else {
                read_pos += len;
              }
            } else if (c == 'D') {
              // 缺失
              if (include_deletions) {
                std::string nucleotide = distinguish_nucleotides ? "-" : "*";
                for (int k = 0; k < len; ++k) {
                  if (ref_pos >= start && ref_pos <= end) {
                    std::pair<std::string, std::string> key = std::make_pair(strand, nucleotide);
                    pileup_map[ref_pos][key]++;
                  }
                  ref_pos++;
                }
              } else {
                ref_pos += len;
              }
            } else if (c == 'N') {
              // 跳过
              ref_pos += len;
            } else if (c == 'S') {
              // 软剪辑
              read_pos += len;
            } else if (c == 'H') {
              // 硬剪辑，不消耗read
              // 什么都不做
            }
          }
        }
        j++;
      }
    }
    
    Rcpp::Rcout << "处理完成: " << processed_alignments << " 条总记录, " 
                << valid_alignments << " 条有效记录\n";
    Rcpp::Rcout << "找到的独特位置数: " << pileup_map.size() << "\n";
    
    
    // 转换为输出向量
    std::vector<std::string> seqnames_out;
    std::vector<int> pos_out;
    std::vector<std::string> strand_out;
    std::vector<std::string> nucleotide_out;
    std::vector<int> count_out;
    
    int total_rows = 0;
    for (const auto& pos_kv : pileup_map) {
      int pos = pos_kv.first;
      
      // 步骤1：计算总深度和主要/次要等位基因
      int total_pos_depth = 0;
      std::map<std::string, int> nt_counts;  // 核苷酸 -> 总计数（忽略链）
      
      for (const auto& count_kv : pos_kv.second) {
        int count = count_kv.second;
        std::string nucleotide = count_kv.first.second;
        
        total_pos_depth += count;
        
        // 只统计真实的核苷酸（忽略占位符 "*"）
        if (nucleotide != "*") {
          nt_counts[nucleotide] += count;
        }
      }
      
      // 步骤2：深度过滤
      if (total_pos_depth < min_depth || total_pos_depth > max_depth) {
        continue;
      }
      
      // 步骤3：次要等位基因频率过滤（如果需要）
      bool pass_maf = true;
      if (min_af > 0.0) {
        // 找到计数最高的两个核苷酸
        int max1 = 0, max2 = 0;
        for (const auto& nt : nt_counts) {
          int count = nt.second;
          if (count > max1) {
            max2 = max1;
            max1 = count;
          } else if (count > max2) {
            max2 = count;
          }
        }
        
        // 计算次要等位基因频率
        double maf = (max2 > 0) ? (double)max2 / total_pos_depth : 0.0;
        
        if (maf < min_af) {
          pass_maf = false;
        }
      }
      
      // 步骤4：通过所有过滤，输出数据
      if (pass_maf) {
        for (const auto& count_kv : pos_kv.second) {
          std::string strand_char = count_kv.first.first;
          std::string nucleotide_char = count_kv.first.second;
          int count = count_kv.second;
          
          seqnames_out.push_back(seqname);
          pos_out.push_back(pos);
          strand_out.push_back(strand_char);
          nucleotide_out.push_back(nucleotide_char);
          count_out.push_back(count);
          total_rows++;
        }
      }
    }
    
    
    
    Rcpp::Rcout << "生成 " << total_rows << " 行数据\n";
    
    if (seqnames_out.empty()) {
      Rcpp::Rcout << "警告: 没有符合条件的数据，返回空DataFrame\n";
      return List::create(
        Named("seqnames") = CharacterVector(0),
        Named("pos") = IntegerVector(0),
        Named("strand") = CharacterVector(0),
        Named("nucleotide") = CharacterVector(0),
        Named("count") = IntegerVector(0)
      );
    }
    
    // 输出前几行作为示例
    Rcpp::Rcout << "示例数据 (前10行):\n";
    for (size_t i = 0; i < std::min((size_t)10, seqnames_out.size()); ++i) {
      Rcpp::Rcout << "  " << seqnames_out[i] << "\t" 
                  << pos_out[i] << "\t" 
                  << strand_out[i] << "\t" 
                  << nucleotide_out[i] << "\t" 
                  << count_out[i] << "\n";
    }
    
    Rcpp::Rcout << "=== pileup_manual_vectors完成 ===\n\n";
    
    // 创建DataFrame
    return List::create(
      Named("seqnames") = seqnames_out,
      Named("pos") = pos_out,
      Named("strand") = strand_out,
      Named("nucleotide") = nucleotide_out,
      Named("count") = count_out
    );
    
  } catch (const std::exception& e) {
    Rcpp::Rcout << "标准异常捕获: " << e.what() << "\n";
    Rcpp::stop("Error in pileup_manual_vectors: %s", e.what());
  } catch (...) {
    Rcpp::Rcout << "未知异常捕获\n";
    Rcpp::stop("Unknown error in pileup_manual_vectors");
  }
}

// 保留原始函数供向后兼容
// [[Rcpp::export]]
DataFrame pileup_manual(
    std::string bam_file,
    std::string seqname,
    int start,
    int end,
    int min_mapq = 0,
    float min_af=0.2,
    int min_base_quality = 0,
    int min_depth = 0,
    int max_depth = 1000000,
    bool distinguish_strands = false,
    bool distinguish_nucleotides = true,
    bool ignore_query_Ns = true,
    bool include_deletions = false,
    bool include_insertions = false
) {
  Rcpp::Rcout << "\n=== 开始pileup_manual(原始版本)调试 ===\n";
  Rcpp::Rcout << "注意: 此版本从BAM文件读取数据\n";
  
  std::string bai_path = bam_file + ".bai";
  Rcpp::Rcout << "BAI文件路径: " << bai_path << "\n";
  
  try {
    Environment rsamtools = Environment::namespace_env("Rsamtools");
    Function scan_bam = rsamtools["scanBam"];
    Function scan_bam_param = rsamtools["ScanBamParam"];
    
    CharacterVector what_fields = CharacterVector::create(
      "pos", "cigar", "seq", "qual", "mapq", "flag", "qwidth"
    );
    
    // Construct 'which' as a plain list (avoids S4 GRanges)
    CharacterVector which_seqnames = CharacterVector::create(seqname);
    List which_ranges = List::create(
      Named("start") = start,
      Named("end") = end
    );
    List which = List::create(
      Named("seqnames") = which_seqnames,
      Named("ranges") = which_ranges
    );
    
    Rcpp::Rcout << "创建ScanBamParam参数...\n";
    SEXP param = scan_bam_param(Named("what") = what_fields, Named("which") = which);
    
    // Call scanBam with explicit index
    Rcpp::Rcout << "调用scanBam...\n";
    List bam_data = scan_bam(
      Named("file") = bam_file,
      Named("index") = bai_path,
      Named("param") = param
    );
    
    Rcpp::Rcout << "scanBam返回结果大小: " << bam_data.size() << "\n";
    
    // 提取vector
    if (bam_data.size() > 0) {
      List first_list = as<List>(bam_data[0]);
      
      // 直接提取vectors
      IntegerVector positions = first_list["pos"];
      CharacterVector cigars = first_list["cigar"];
      CharacterVector sequences = first_list["seq"];
      CharacterVector qualities = first_list["qual"];
      IntegerVector mapqs = first_list["mapq"];
      IntegerVector flags = first_list["flag"];
      IntegerVector qwidths = first_list["qwidth"];
      
      // 调用vectors版本
      return pileup_manual_vectors(
        CharacterVector(positions.size(), seqname),
        positions,
        cigars,
        sequences,
        qualities,
        mapqs,
        flags,
        qwidths,
        seqname,
        start,
        end,
        min_mapq,
        min_base_quality,
        min_depth,
        min_af,
        max_depth,
        distinguish_strands,
        distinguish_nucleotides,
        ignore_query_Ns,
        include_deletions,
        include_insertions
      );
    } else {
      Rcpp::Rcout << "警告: 没有比对数据，返回空DataFrame\n";
      return DataFrame::create(
        Named("seqnames") = CharacterVector(0),
        Named("pos") = IntegerVector(0),
        Named("count") = IntegerVector(0)
      );
    }
    
  } catch (const std::exception& e) {
    Rcpp::Rcout << "标准异常捕获: " << e.what() << "\n";
    Rcpp::stop("Error in pileup_manual: %s", e.what());
  } catch (...) {
    Rcpp::Rcout << "未知异常捕获\n";
    Rcpp::stop("Unknown error in pileup_manual");
  }
}
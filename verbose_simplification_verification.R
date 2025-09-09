# 详细日志格式精简验证脚本

library(data.table)

# 函数：分析verbose消息的简洁性
analyze_verbose_messages <- function() {
  cat("=== 详细日志格式精简验证报告 ===\n\n")
  
  # 搜索所有geneSCOPE::函数名的消息
  all_files <- list.files("/Users/haenolab/Documents/FG2CIL Paper/Code/FG2CLI/R", 
                         pattern = "\\.r$", full.names = TRUE)
  
  total_messages <- 0
  function_summary <- list()
  
  for (file in all_files) {
    lines <- readLines(file, warn = FALSE)
    verbose_lines <- grep("geneSCOPE::", lines, value = TRUE)
    
    if (length(verbose_lines) > 0) {
      filename <- basename(file)
      cat("文件:", filename, "\n")
      
      # 提取函数名
      for (line in verbose_lines) {
        if (grepl("message.*geneSCOPE::", line)) {
          # 提取函数名
          func_match <- regmatches(line, regexpr("geneSCOPE::[a-zA-Z_][a-zA-Z0-9_]*", line))
          if (length(func_match) > 0) {
            func_name <- sub("geneSCOPE::", "", func_match)
            if (is.null(function_summary[[func_name]])) {
              function_summary[[func_name]] <- 0
            }
            function_summary[[func_name]] <- function_summary[[func_name]] + 1
            total_messages <- total_messages + 1
            
            # 检查是否有适当的缩进
            has_indent <- grepl("^[[:space:]]*if.*message.*[[:space:]]{2,}", line)
            indent_status <- ifelse(has_indent, "✓缩进", "基础")
            
            # 简化显示
            clean_line <- gsub("^[[:space:]]*", "", line)
            clean_line <- gsub("if \\(verbose\\) ", "", clean_line)
            cat("  ", indent_status, ":", clean_line, "\n")
          }
        }
      }
      cat("\n")
    }
  }
  
  cat("=== 函数详细消息统计 ===\n")
  sorted_funcs <- sort(unlist(function_summary), decreasing = TRUE)
  for (i in 1:min(10, length(sorted_funcs))) {
    func_name <- names(sorted_funcs)[i]
    count <- sorted_funcs[i]
    status <- ifelse(count > 8, "⚠️过多", ifelse(count > 5, "🔶适中", "✅简洁"))
    cat(sprintf("%s %s: %d条消息\n", status, func_name, count))
  }
  
  cat("\n=== 总体评估 ===\n")
  cat("总消息数:", total_messages, "\n")
  cat("涉及函数:", length(function_summary), "\n")
  
  # 检查是否还有步骤编号
  step_check <- system("grep -r 'Step [0-9]' /Users/haenolab/Documents/FG2CIL\\ Paper/Code/FG2CLI/R/ 2>/dev/null || echo 'No step numbering found'", intern = TRUE)
  cat("步骤编号检查:", ifelse(length(step_check) == 1 && step_check == "No step numbering found", "✅已清理", "⚠️仍有残留"), "\n")
  
  # 检查是否还有省略号
  ellipsis_check <- system("grep -r '\\.\\.\\.' /Users/haenolab/Documents/FG2CIL\\ Paper/Code/FG2CLI/R/ 2>/dev/null || echo 'No ellipsis found'", intern = TRUE)
  cat("省略号检查:", ifelse(length(ellipsis_check) == 1 && ellipsis_check == "No ellipsis found", "✅已清理", "⚠️仍有残留"), "\n")
  
  cat("\n=== 优化建议 ===\n")
  high_volume_funcs <- names(sorted_funcs)[sorted_funcs > 8]
  if (length(high_volume_funcs) > 0) {
    cat("建议进一步精简以下函数的详细消息:\n")
    for (func in high_volume_funcs) {
      cat("- ", func, "(", sorted_funcs[func], "条消息)\n")
    }
  } else {
    cat("✅ 所有函数的详细消息数量都在合理范围内\n")
  }
  
  cat("\n=== 验证完成 ===\n")
}

# 运行分析
analyze_verbose_messages()

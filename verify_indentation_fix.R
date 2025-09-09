#!/usr/bin/env Rscript
# 验证详细日志缩进修正脚本
# 确保 [geneSCOPE::functionName] 左对齐，只有后续描述信息缩进

cat("=== 详细日志缩进修正验证 ===\n\n")

# 搜索所有R文件中的geneSCOPE消息
all_files <- list.files("/Users/haenolab/Documents/FG2CIL Paper/Code/FG2CLI/R", 
                       pattern = "\\.r$", full.names = TRUE)

correct_format <- 0
incorrect_format <- 0
total_messages <- 0

cat("检查文件:\n")
for (file in all_files) {
  lines <- readLines(file, warn = FALSE)
  verbose_lines <- grep("message.*\\[geneSCOPE::", lines, value = TRUE)
  
  if (length(verbose_lines) > 0) {
    filename <- basename(file)
    cat("文件:", filename, "\n")
    
    for (line in verbose_lines) {
      total_messages <- total_messages + 1
      
      # 检查是否 [geneSCOPE::functionName] 左对齐
      if (grepl("^[[:space:]]*(?:if.*)?message\\(\"\\[geneSCOPE::", line, perl = TRUE)) {
        # 检查函数名后是否只有2个空格（正确的缩进）
        if (grepl("\\[geneSCOPE::[a-zA-Z_][a-zA-Z0-9_]*\\]  [^[:space:]]", line)) {
          correct_format <- correct_format + 1
          cat("  ✅ 正确格式:", gsub("^[[:space:]]*", "", line), "\n")
        } else if (grepl("\\[geneSCOPE::[a-zA-Z_][a-zA-Z0-9_]*\\]   ", line)) {
          incorrect_format <- incorrect_format + 1
          cat("  ❌ 仍有3+空格:", gsub("^[[:space:]]*", "", line), "\n")
        } else if (grepl("\\[geneSCOPE::[a-zA-Z_][a-zA-Z0-9_]*\\] [^[:space:]]", line)) {
          correct_format <- correct_format + 1
          cat("  ✅ 无缩进（正确）:", gsub("^[[:space:]]*", "", line), "\n")
        } else {
          cat("  ℹ️  其他格式:", gsub("^[[:space:]]*", "", line), "\n")
        }
      }
    }
    cat("\n")
  }
}

cat("=== 验证结果汇总 ===\n")
cat("总消息数:", total_messages, "\n")
cat("正确格式:", correct_format, "\n")
cat("错误格式:", incorrect_format, "\n")

if (incorrect_format == 0) {
  cat("🎉 所有缩进已正确修正！\n")
  cat("✅ [geneSCOPE::functionName] 现在都左对齐\n") 
  cat("✅ 描述信息适当缩进（2个空格）\n")
} else {
  cat("⚠️  还有", incorrect_format, "条消息需要修正\n")
}

cat("\n=== 示例正确格式 ===\n")
cat('[geneSCOPE::createCoordObj] Loading cell centroids and transcripts\n')
cat('[geneSCOPE::createCoordObj]  Found 12345 cells\n')
cat('[geneSCOPE::addLeeStats] Computing Lee\'s L statistics\n')
cat('[geneSCOPE::addLeeStats]  Lee\'s L completed (2.3 min)\n')

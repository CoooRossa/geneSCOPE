#!/bin/bash
# 批量更新GitHub链接脚本
# 将所有haenolab/geneSCOPE链接更新为CoooRossa/geneSCOPE

echo "🔗 更新GitHub链接开始..."
echo "将 haenolab/geneSCOPE 更新为 CoooRossa/geneSCOPE"

cd "/Users/haenolab/Documents/FG2CIL Paper/Code/FG2CLI"

# 需要更新的文件列表
FILES_TO_UPDATE=(
    "install_hpc_dependencies.R"
    "install_genescope_hpc.R" 
    "diagnose_hpc.R"
    "README.md"
    "HPC_TROUBLESHOOTING.md"
    "MISSION_ACCOMPLISHED.md"
    "HPC_INSTALL.md"
    "emergency_hpc_install.R"
    "QUICK_START.md"
)

echo ""
echo "📝 更新文件中的GitHub链接..."

for file in "${FILES_TO_UPDATE[@]}"; do
    if [ -f "$file" ]; then
        echo "   更新 $file..."
        
        # 使用sed替换haenolab/geneSCOPE为CoooRossa/geneSCOPE
        sed -i.bak 's/haenolab\/geneSCOPE/CoooRossa\/geneSCOPE/g' "$file"
        
        # 删除备份文件
        rm -f "$file.bak"
        
        echo "   ✅ $file 已更新"
    else
        echo "   ⚠️  $file 不存在，跳过"
    fi
done

echo ""
echo "🔍 验证更新结果..."

# 检查是否还有遗漏的haenolab链接
REMAINING=$(grep -r "haenolab/geneSCOPE" . --include="*.md" --include="*.R" 2>/dev/null | wc -l)

if [ "$REMAINING" -eq 0 ]; then
    echo "✅ 所有haenolab/geneSCOPE链接已成功更新为CoooRossa/geneSCOPE"
else
    echo "⚠️  还有 $REMAINING 个链接需要手动检查:"
    grep -r "haenolab/geneSCOPE" . --include="*.md" --include="*.R" 2>/dev/null
fi

echo ""
echo "📊 更新统计:"
echo "   - 已处理文件: ${#FILES_TO_UPDATE[@]} 个"
echo "   - 新的GitHub用户: CoooRossa" 
echo "   - 仓库名称: geneSCOPE"

echo ""
echo "✅ GitHub链接更新完成!"

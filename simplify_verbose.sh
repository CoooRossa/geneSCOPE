#!/bin/bash
# 精简verbose消息脚本

cd "/Users/haenolab/Documents/FG2CIL Paper/Code/FG2CLI/R"

# 1. 精简createCoordObj的verbose消息，保留关键节点
sed -i '' 's/Step [0-9][a-z]* - //g' 1.CreateObj.r
sed -i '' 's/Step [0-9]* - //g' 1.CreateObj.r

# 2. 精简norMPG的verbose消息
sed -i '' 's/Step [0-9][a-z]* - //g' 2.Normalization.r  
sed -i '' 's/Step [0-9]* - //g' 2.Normalization.r

# 3. 保留但简化其他函数的verbose消息
# 移除所有"Step X -"格式，保留函数名和简洁描述

echo "Verbose消息精简完成"

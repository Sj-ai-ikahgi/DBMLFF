#!/bin/bash

# 定义要删除的文件模式
FILES_TO_DELETE=("OUT.GKK" "OUT.RHO" "OUT.VR" "OUT.VR_hion" "OUT.VR_ion")

# 检查SCF目录是否存在
if [ ! -d "./SCF" ]; then
    echo "错误: ./SCF 目录不存在"
    exit 1
fi

echo "正在扫描 ./SCF 目录及其子目录中的目标文件..."

# 查找并显示将要删除的文件
total_files=0
for file_pattern in "${FILES_TO_DELETE[@]}"; do
    file_count=$(find ./SCF -name "$file_pattern" -type f | wc -l)
    total_files=$((total_files + file_count))
    if [ $file_count -gt 0 ]; then
        echo "找到 $file_count 个 $file_pattern 文件"
    fi
done

if [ $total_files -eq 0 ]; then
    echo "没有找到需要删除的文件"
    exit 0
fi

echo ""
echo "总共找到 $total_files 个文件将被删除"
echo "文件类型: ${FILES_TO_DELETE[*]}"

# 确认删除
read -p "确定要删除这些文件吗？(y/N): " confirm
if [[ $confirm != [yY] ]]; then
    echo "操作已取消"
    exit 0
fi

# 执行删除操作
echo "开始删除文件..."
for file_pattern in "${FILES_TO_DELETE[@]}"; do
    find ./SCF -name "$file_pattern" -type f -delete
    echo "已删除所有 $file_pattern 文件"
done

echo "文件删除完成！"
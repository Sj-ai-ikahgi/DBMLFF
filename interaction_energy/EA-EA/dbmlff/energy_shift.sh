#!/bin/bash

# 定义文件路径和特定值
input_file="./DB.ENERGY"   # 输入文件
output_file="./DB.energy"  # 输出文件

# energy_shift
ea=-3714.562757000
tte=-11782.379870000
tfsi=-14773.203830000
li=-300.9499727809445



# 计算增量值
increment_value=$(echo "scale=10; $ea * 2" | bc)  # 精度设置为10

# 检查输入文件是否存在
if [[ ! -f "$input_file" ]]; then
  echo "错误: 输入文件 $input_file 不存在。"
  exit 1
fi

# 确保输出目录存在
output_dir=$(dirname "$output_file")
if [[ ! -d "$output_dir" ]]; then
  echo "错误: 输出目录 $output_dir 不存在。"
  exit 1
fi

# 清空输出文件，防止追加
> "$output_file"

# 处理文件，逐行读取并进行加法操作
while IFS= read -r line; do
  # 通过加法操作处理每一行，假设文件中的每一行是一个数字
  new_value=$(echo "scale=10; $line + $increment_value" | bc)
  echo "$new_value" >> "$output_file"  # 将结果输出到新文件
done < "$input_file"

echo "处理完成，结果保存在 $output_file"

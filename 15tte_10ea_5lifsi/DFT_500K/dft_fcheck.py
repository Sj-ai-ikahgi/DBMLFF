def reorder_file(input_file, output_file, index_file):
    # 读取映射文件并构建映射关系
    mapping = {}
    with open(index_file, 'r', encoding='utf-8') as index_f:
        for line in index_f:
            # 将每行分割为多个部分
            parts = line.split()
            if len(parts) >= 5:  # 确保每行有至少五列
                new_index = int(parts[0])      # 第一列为新索引
                original_index = int(parts[-2])  # 倒数第二列为初始行号
                mapping[original_index] = new_index

    # 读取输入文件内容
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # 按模块重排
    reordered_lines = []
    num_lines = len(lines)
    for i in range(0, num_lines, 460):
        module = lines[i:i + 460]  # 每个模块
        reordered_module = [''] * 460  # 初始化重排后的模块
        for original_index, new_index in mapping.items():
            if original_index <= len(module):  # 确保索引在模块范围内
                reordered_module[new_index - 1] = module[original_index - 1]  # 重新排序
        reordered_lines.extend(reordered_module)

    # 删除第一列并重新编号，每行一个数字
    final_lines = []
    for line in reordered_lines:
        parts = line.split()
        if len(parts) > 1:  # 确保有多列数据
            # 删除第一列并将剩余部分的每个数字按顺序重排，并转换为浮点数格式
            for val in parts[1:]:  # 从第二列开始处理
                # 将每个数字转换为浮点数并格式化
                final_lines.append(f"{float(val):.6f}\n")

    # 写入输出文件
    with open(output_file, 'w', encoding='utf-8') as f:
        f.writelines(final_lines)

# 使用示例
input_file = 'OUT.FORCE'  # 输入文件名
output_file = './DFT.FORCE'  # 输出文件名
index_file = './SCF/extracted_iteration_1/ORIGIN.INDEX'  # 映射关系文件名
reorder_file(input_file, output_file, index_file)
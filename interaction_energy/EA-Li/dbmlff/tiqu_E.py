def extract_lines(input_file, output_file, start_line, end_line):
    """
    从输入文件中提取指定范围的行并保存到输出文件中。

    :param input_file: 输入文件路径
    :param output_file: 输出文件路径
    :param start_line: 起始行号（1-based）
    :param end_line: 结束行号（1-based）
    """
    # 打开输入文件读取内容
    with open(input_file, 'r', encoding='utf-8') as infile:
        lines = infile.readlines()

    # 确保行号在有效范围内
    if start_line < 1 or end_line > len(lines):
        print(f"错误：行号范围无效。文件共有 {len(lines)} 行。")
        return

    # 提取指定范围的行
    extracted_lines = lines[start_line - 1:end_line]  # Python 索引从 0 开始，所以需要减去 1

    # 将提取的行保存到输出文件
    with open(output_file, 'w', encoding='utf-8') as outfile:
        outfile.writelines(extracted_lines)

    print(f"成功提取第 {start_line} 行到第 {end_line} 行，保存到 {output_file}。")


def extract_and_reorder_file(input_file, extracted_file, remaining_file):
    # 读取输入文件内容
    with open(input_file, 'r', encoding='utf-8') as infile:
        lines = infile.readlines()

    # 初始化用于存储提取的行和剩余行的列表
    extracted_lines = []
    new_lines = []

    # 计算提取的行号（每29行提取一次）
    indices_to_extract = set(range(1, len(lines) + 1, 16))

    # 提取指定行并将其余行处理
    for i, line in enumerate(lines, start=1):
        if i in indices_to_extract:
            extracted_lines.append(line.strip())  # 提取的行保留原样
        else:
            new_lines.append(line)  # 剩余行进行处理

    # 保存提取的行到文件
    with open(extracted_file, 'w', encoding='utf-8') as outfile_extracted:
        outfile_extracted.write("\n".join(extracted_lines) + "\n")

    # 对剩余行进行重新排序和格式化
    final_lines = []
    for line in new_lines:
        parts = line.split()
        if len(parts) > 1:  # 确保每行有多列数据
            for val in parts[1:]:  # 从第二列开始处理
                # 将每个数字转换为浮点数并格式化，保留6位小数
                final_lines.append(f"{float(val):.6f}\n")

    # 保存剩余行处理后的数据到文件
    with open(remaining_file, 'w', encoding='utf-8') as outfile_remaining:
        outfile_remaining.writelines(final_lines)

    print("提取的行和剩余行已分别保存。")


def main():
    # 第一步：提取特定范围的行
    input_file = './molecule_pair/OUT.energy_force'  # 输入文件路径
    output_file = './extracted.txt'    # 提取后的文件路径
    start_line = 1  # 要提取的起始行号
    end_line = 3200    # 要提取的结束行号
    extract_lines(input_file, output_file, start_line, end_line)

    # 第二步：提取并重新排序文件内容
    extracted_file = './DB.ENERGY'    # 提取的行文件名
    remaining_file = './DB.FORCE'     # 剩余行文件名
    extract_and_reorder_file(output_file, extracted_file, remaining_file)


# 执行主函数
if __name__ == "__main__":
    main()

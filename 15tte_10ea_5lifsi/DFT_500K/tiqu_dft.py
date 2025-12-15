import os
import re

# 设置起止范围
start = 1
end =50

# 定义文件夹路径
base_dir = './SCF'
force_output_file = 'OUT.FORCE'
energy_output_file = './DFT.ENERGY'

# 定义一个函数将科学计数法转换为非科学计数法
def sci_to_float(sci_str):
    try:
        return f"{float(sci_str):.6f}"  # 保留6位小数
    except ValueError:
        return sci_str  # 如果转换失败，返回原字符串

# 创建输出文件
with open(force_output_file, 'w') as force_file, open(energy_output_file, 'w') as energy_file:
    for i in range(start, end + 1):
        # 生成文件夹名称
        folder_name = f'extracted_iteration_{i}'
        folder_path = os.path.join(base_dir, folder_name)
        
        # 检查文件夹是否存在
        if os.path.isdir(folder_path):
            # 处理 OUT.FORCE 文件
            force_file_path = os.path.join(folder_path, 'OUT.FORCE')
            if os.path.isfile(force_file_path):
                with open(force_file_path, 'r') as f:
                    lines = f.readlines()
                    # 找到特定行
                    for j, line in enumerate(lines):
                        if '****** force after remove total force (eV/A) ***' in line:
                            # 写入该行之后的内容并转换
                            for subsequent_line in lines[j + 1:]:
                                # 查找科学计数法的数
                                converted_line = re.sub(r'([-+]?\d*\.?\d+E[+-]?\d+)', lambda m: sci_to_float(m.group(0)), subsequent_line)
                                force_file.write(converted_line)
                            break
            
        # 处理 REPORT 文件
        report_file_path = os.path.join(folder_path, 'REPORT')
        if os.path.isfile(report_file_path):
            with open(report_file_path, 'r') as f:
                lines = f.readlines()
                # 倒序查找最后出现的 E_tot(eV) =
                for line in reversed(lines):
                    match = re.search(r'E_tot\(eV\)\s*=\s*([-\d.]+E[+-]?\d+)', line)
                    if match:
                        energy_value = sci_to_float(match.group(1))
                        energy_file.write(energy_value + '\n')
                        break

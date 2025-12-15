#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统计异常分子移除数量评估脚本
统计规则：
1. check_conf/目录下的*_DFT_large_distortion.dat文件行数代表异常分子数量
2. Conf_*文件夹中的文件数代表该温度数据集的分子总数
"""

import os
import glob
import sys

def count_files_in_directory(directory):
    """统计目录中的文件数量（不包括子目录）"""
    if not os.path.exists(directory):
        return 0
    
    count = 0
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        if os.path.isfile(item_path):
            count += 1
    return count

def count_lines_in_file(filepath):
    """统计文件中的行数"""
    if not os.path.exists(filepath):
        return 0
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            return sum(1 for _ in f)
    except Exception as e:
        print(f"警告：无法读取文件 {filepath}: {e}")
        return 0

def find_matching_temperatures():
    """查找匹配的温度文件夹和大变形文件"""
    temperature_data = {}
    
    # 查找所有Conf_*文件夹
    conf_folders = glob.glob('Conf_*')
    
    for folder in conf_folders:
        if os.path.isdir(folder):
            # 从文件夹名提取温度标识（移除Conf_前缀）
            temp_id = folder.replace('Conf_', '', 1)
            temperature_data[temp_id] = {
                'folder': folder,
                'distortion_file': None,
                'total_molecules': 0,
                'abnormal_molecules': 0
            }
    
    # 查找所有大变形文件
    check_conf_dir = 'check_conf'
    if os.path.exists(check_conf_dir):
        distortion_files = glob.glob(os.path.join(check_conf_dir, '*_DFT_large_distortion.dat'))
        
        for file_path in distortion_files:
            filename = os.path.basename(file_path)
            # 从文件名提取温度标识（移除_DFT_large_distortion.dat后缀）
            if filename.endswith('_DFT_large_distortion.dat'):
                temp_id = filename.replace('_DFT_large_distortion.dat', '')
                
                # 精确匹配：检查这个温度标识是否存在于temperature_data中
                if temp_id in temperature_data:
                    temperature_data[temp_id]['distortion_file'] = file_path
                else:
                    # 如果精确匹配失败，尝试部分匹配（移除可能的_s后缀）
                    if temp_id.endswith('_s'):
                        base_temp_id = temp_id[:-2]  # 移除_s
                        if base_temp_id in temperature_data:
                            print(f"警告：文件 {filename} 匹配到文件夹 Conf_{base_temp_id} (移除_s后缀)")
                            temperature_data[base_temp_id]['distortion_file'] = file_path
    
    return temperature_data

def calculate_removal_statistics():
    """计算移除统计数据"""
    print("开始统计异常分子移除数据...")
    
    # 查找匹配的温度数据
    temperature_data = find_matching_temperatures()
    
    if not temperature_data:
        print("未找到任何Conf_*文件夹！")
        return None
    
    print(f"找到 {len(temperature_data)} 个温度数据集:")
    for temp_id, data in sorted(temperature_data.items()):
        distortion_file = data['distortion_file'] if data['distortion_file'] else "未找到对应文件"
        print(f"  {temp_id}: 文件夹={data['folder']}, 大变形文件={distortion_file}")
    
    # 统计每个温度的数据
    for temp_id, data in temperature_data.items():
        # 统计总分子数
        data['total_molecules'] = count_files_in_directory(data['folder'])
        
        # 统计异常分子数
        if data['distortion_file']:
            data['abnormal_molecules'] = count_lines_in_file(data['distortion_file'])
        
        # 计算移除率和保留率
        if data['total_molecules'] > 0:
            data['removal_rate'] = data['abnormal_molecules'] / data['total_molecules'] * 100
            data['retention_rate'] = 100 - data['removal_rate']
        else:
            data['removal_rate'] = 0
            data['retention_rate'] = 0
    
    return temperature_data

def generate_summary_report(temperature_data, output_file='removed_molecules_summary.dat'):
    """生成汇总报告"""
    
    total_molecules_all = 0
    total_abnormal_all = 0
    
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入标题
        f.write("=" * 80 + "\n")
        f.write("异常分子移除数量评估报告\n")
        f.write("=" * 80 + "\n\n")
        
        # 写入统计表格标题
        f.write(f"{'温度':<10} {'总分子数':<12} {'异常分子数':<12} {'移除率(%)':<12} {'保留率(%)':<12} {'Conf文件夹':<20}\n")
        f.write("-" * 80 + "\n")
        
        # 写入每个温度的数据
        for temp_id, data in sorted(temperature_data.items()):
            total_molecules = data['total_molecules']
            abnormal_molecules = data['abnormal_molecules']
            removal_rate = data.get('removal_rate', 0)
            retention_rate = data.get('retention_rate', 100)
            
            f.write(f"{temp_id:<10} {total_molecules:<12} {abnormal_molecules:<12} "
                   f"{removal_rate:<12.2f} {retention_rate:<12.2f} {data['folder']:<20}\n")
            
            total_molecules_all += total_molecules
            total_abnormal_all += abnormal_molecules
        
        # 写入总计
        f.write("-" * 80 + "\n")
        
        if total_molecules_all > 0:
            total_removal_rate = total_abnormal_all / total_molecules_all * 100
            total_retention_rate = 100 - total_removal_rate
        else:
            total_removal_rate = 0
            total_retention_rate = 0
        
        f.write(f"{'总计':<10} {total_molecules_all:<12} {total_abnormal_all:<12} "
               f"{total_removal_rate:<12.2f} {total_retention_rate:<12.2f} {'-':<20}\n")
        
        f.write("=" * 80 + "\n\n")
        
        # 写入统计说明
        f.write("统计说明:\n")
        f.write("1. 总分子数: 对应Conf_*文件夹中的文件总数\n")
        f.write("2. 异常分子数: check_conf/目录下*_DFT_large_distortion.dat文件的行数\n")
        f.write("3. 移除率: (异常分子数 / 总分子数) × 100%\n")
        f.write("4. 保留率: (正常分子数 / 总分子数) × 100%\n")
        f.write("5. 每个温度对应关系: Conf_{温度}文件夹 ↔ check_conf/{温度}_DFT_large_distortion.dat\n\n")
        
        # 写入文件对应关系
        f.write("检测到的对应关系:\n")
        for temp_id, data in sorted(temperature_data.items()):
            distortion_file = data['distortion_file'] if data['distortion_file'] else "未找到对应文件"
            f.write(f"  {data['folder']} ↔ {distortion_file}\n")
    
    print(f"\n汇总报告已保存到: {output_file}")
    
    return output_file

def print_summary_to_console(temperature_data):
    """在控制台打印汇总信息"""
    print("\n" + "=" * 80)
    print("异常分子移除数量评估")
    print("=" * 80)
    print(f"{'温度':<10} {'总分子数':<12} {'异常分子数':<12} {'移除率(%)':<12} {'保留率(%)':<12}")
    print("-" * 80)
    
    total_molecules_all = 0
    total_abnormal_all = 0
    
    for temp_id, data in sorted(temperature_data.items()):
        total_molecules = data['total_molecules']
        abnormal_molecules = data['abnormal_molecules']
        removal_rate = data.get('removal_rate', 0)
        retention_rate = data.get('retention_rate', 100)
        
        print(f"{temp_id:<10} {total_molecules:<12} {abnormal_molecules:<12} "
              f"{removal_rate:<12.2f} {retention_rate:<12.2f}")
        
        total_molecules_all += total_molecules
        total_abnormal_all += abnormal_molecules
    
    print("-" * 80)
    
    if total_molecules_all > 0:
        total_removal_rate = total_abnormal_all / total_molecules_all * 100
        total_retention_rate = 100 - total_removal_rate
    else:
        total_removal_rate = 0
        total_retention_rate = 0
    
    print(f"{'总计':<10} {total_molecules_all:<12} {total_abnormal_all:<12} "
          f"{total_removal_rate:<12.2f} {total_retention_rate:<12.2f}")
    print("=" * 80)

def main():
    """主函数"""
    print("异常分子移除统计脚本")
    print("版本: 1.0")
    print("=" * 50)
    
    # 检查check_conf目录是否存在
    if not os.path.exists('check_conf'):
        print("错误: check_conf目录不存在！")
        print("请确保脚本在与check_conf目录相同的目录下运行。")
        sys.exit(1)
    
    # 计算统计数据
    temperature_data = calculate_removal_statistics()
    
    if not temperature_data:
        print("没有找到可统计的数据，程序退出。")
        sys.exit(0)
    
    # 在控制台打印结果
    print_summary_to_console(temperature_data)
    
    # 生成汇总报告文件
    output_file = 'removed_molecules_summary.dat'
    generate_summary_report(temperature_data, output_file)
    
    print(f"\n操作完成！")
    print(f"详细报告已保存到: {os.path.abspath(output_file)}")

if __name__ == "__main__":
    main()
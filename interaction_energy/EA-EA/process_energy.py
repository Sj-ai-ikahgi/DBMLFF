#!/usr/bin/env python3
import numpy as np
import os

def read_energy_values(file_path):
    """
    从文件中读取一列能量数据，每行一个值。
    """
    try:
        return np.loadtxt(file_path, dtype=float)
    except FileNotFoundError:
        print(f"错误：文件未找到: {file_path}")
        exit(1)

def read_single_energy(file_path):
    """
    从文件中读取单个分子的能量值（读取第一个有效数字）。
    """
    try:
        with open(file_path, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                if stripped_line and not stripped_line.startswith('#'):
                    return float(stripped_line)
        
        print(f"错误：在文件中未找到有效数值: {file_path}")
        exit(1)

    except FileNotFoundError:
        print(f"错误：文件未找到: {file_path}")
        exit(1)

def main():
    # --- 1. 定义文件路径 ---
    # 分子对能量文件
    dft_pair_path = './dft/DFT.ENERGY'
    dbmlff_pair_path = './dbmlff/DB.ENERGY'
    
    # 单分子能量文件（前100行和后100行分别使用）
    dft_single_path_1_100 = './dft/single_molecule_1_100/energy.txt'
    dft_single_path_101_200 = './dft/single_molecule_101_200/energy.txt'
    
    dbmlff_single_path_1_100 = './dbmlff/single_molecule_1_100/energy.txt'
    dbmlff_single_path_101_200 = './dbmlff/single_molecule_101_200/energy.txt'
    
    # 输出文件
    output_path = 'E.txt'

    # --- 2. 读取所有需要的数据 ---
    # 读取分子对的能量列表
    dft_pair_energies = read_energy_values(dft_pair_path)
    dbmlff_pair_energies = read_energy_values(dbmlff_pair_path)

    # 读取两组单分子的能量值
    dft_single_energy_1_100 = read_single_energy(dft_single_path_1_100)
    dft_single_energy_101_200 = read_single_energy(dft_single_path_101_200)
    
    dbmlff_single_energy_1_100 = read_single_energy(dbmlff_single_path_1_100)
    dbmlff_single_energy_101_200 = read_single_energy(dbmlff_single_path_101_200)

    # --- 3. 数据一致性检查 ---
    if len(dft_pair_energies) != len(dbmlff_pair_energies):
        print(f"错误：DFT能量数据({len(dft_pair_energies)}行)与DBMLFF能量数据({len(dbmlff_pair_energies)}行)长度不匹配")
        exit(1)
        
    total_lines = len(dft_pair_energies)
    if total_lines != 200:
        print(f"错误：分子对能量文件的总行数应为200，但实际为 {total_lines}。")
        exit(1)


    # --- 4. 分段计算相互作用能 ---
    # 前100行使用第一组单分子能量
    e_dft_part1 = dft_pair_energies[:100] - 2 * dft_single_energy_1_100
    e_dbmlff_part1 = dbmlff_pair_energies[:100] - 2 * dbmlff_single_energy_1_100
    
    # 后100行使用第二组单分子能量
    e_dft_part2 = dft_pair_energies[100:] - 2 * dft_single_energy_101_200
    e_dbmlff_part2 = dbmlff_pair_energies[100:] - 2 * dbmlff_single_energy_101_200
    
    # 合并两部分结果
    e_dft = np.concatenate((e_dft_part1, e_dft_part2))
    e_dbmlff = np.concatenate((e_dbmlff_part1, e_dbmlff_part2))


    # --- 5. 将结果写入输出文件 ---
    try:
        np.savetxt(output_path, np.column_stack((e_dft, e_dbmlff)), 
                   fmt='%.10f', delimiter=' ')
        print(f"成功生成文件: {output_path}")
        print(f"共处理 {len(e_dft)} 个分子对能量数据")
        print(f"  - 前 100 行使用了单分子能量文件: {os.path.abspath(dft_single_path_1_100)} 和 {os.path.abspath(dbmlff_single_path_1_100)}")
        print(f"  - 后 100 行使用了单分子能量文件: {os.path.abspath(dft_single_path_101_200)} 和 {os.path.abspath(dbmlff_single_path_101_200)}")
    except Exception as e:
        print(f"写入文件时出错: {e}")
        exit(1)

if __name__ == "__main__":
    main()
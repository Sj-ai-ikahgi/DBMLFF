#!/usr/bin/env python3
import numpy as np
import os

def read_energy_values(file_path):
    """
    从文件中读取一列能量数据，每行一个值。
    使用 numpy.loadtxt 效率很高。
    """
    try:
        return np.loadtxt(file_path, dtype=float)
    except FileNotFoundError:
        print(f"错误：文件未找到: {file_path}")
        exit(1)

def read_single_energy(file_path):
    """
    从文件中读取单个分子的能量值。
    函数会逐行扫描文件，找到第一个非空、非注释（不以#开头）的行，并返回其浮点数值。
    这确保了我们只读取并使用文件中的第一个有效数字。
    """
    try:
        with open(file_path, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                # 检查行是否有效：非空且不是注释
                if stripped_line and not stripped_line.startswith('#'):
                    # 找到第一个有效数字后立即返回，函数终止
                    return float(stripped_line)
        
        # 如果循环结束仍未找到有效数值，则打印错误并退出
        print(f"错误：在文件中未找到有效数值: {file_path}")
        exit(1)

    except FileNotFoundError:
        print(f"错误：文件未找到: {file_path}")
        exit(1)

def main():
    # --- 1. 定义文件路径 ---
    dft_pair_path = './dft/DFT.ENERGY'
    dft_single_path = './dft/single_molecule/energy.txt'
    dbmlff_pair_path = './dbmlff/DB.ENERGY'
    dbmlff_single_path = './dbmlff/single_molecule/energy.txt'
    output_path = 'E.txt'

    # --- 2. 读取DFT能量并计算相互作用能 ---
    # 读取分子对的能量列表
    dft_pair_energies = read_energy_values(dft_pair_path)
    # 读取单个分子的能量（只取第一行有效数字）
    dft_single_energy = read_single_energy(dft_single_path)
    # 计算相互作用能: E_interaction = E_pair -  * E_single
    # 此计算适用于所有数据点，并且能正确处理负数
    e_dft = dft_pair_energies -  dft_single_energy

    # --- 3. 读取DBMLFF能量并计算相互作用能 ---
    dbmlff_pair_energies = read_energy_values(dbmlff_pair_path)
    dbmlff_single_energy = read_single_energy(dbmlff_single_path)
    e_dbmlff = dbmlff_pair_energies -  dbmlff_single_energy

    # --- 4. 数据一致性检查 ---
    if len(e_dft) != len(e_dbmlff):
        print(f"错误：DFT能量数据({len(e_dft)}行)与DBMLFF能量数据({len(e_dbmlff)}行)长度不匹配")
        exit(1)

    # --- 5. 将结果写入输出文件 ---
    try:
        # np.column_stack 将两个一维数组合并成一个二维数组，每行对应一个分子对
        # fmt='%.10f' 指定输出格式为保留10位小数的浮点数，能正确处理负数
        # delimiter=' ' 指定列之间用空格分隔
        np.savetxt(output_path, np.column_stack((e_dft, e_dbmlff)), 
                   fmt='%.10f', delimiter=' ')
        print(f"成功生成文件: {output_path}")
        print(f"共处理 {len(e_dft)} 个分子对能量数据")
    except Exception as e:
        print(f"写入文件时出错: {e}")
        exit(1)

if __name__ == "__main__":
    main()
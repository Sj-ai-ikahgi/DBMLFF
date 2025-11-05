# -*- coding: utf-8 -*-
import os
import numpy as np
from scipy.spatial.distance import cdist

# 定义文件路径和模式
input_file = './EA.config'
output_folder = './new_structures'
num_structures = 200
min_distance = 1  # 埃
max_distance = 3.0  # 埃
num_attempts = 10000  # 增加尝试次数

# 如果输出文件夹不存在，则创建它
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def read_structure(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    lattice_vectors = []
    positions = []
    reading_positions = False
    
    for line in lines:
        if line.startswith(' Lattice vector (Angstrom)'):
            for i in range(1, 4):
                vec_line = lines[lines.index(line) + i].strip().split()
                lattice_vectors.append(list(map(float, vec_line)))
        
        elif line.strip().startswith('Position (normalized), move_x, move_y, move_z'):
            reading_positions = True
        
        elif reading_positions and line.strip():
            parts = line.split()
            if len(parts) >= 7:  # 确保有足够的部分来解析
                try:
                    atom_type = int(parts[0])
                    position = list(map(float, parts[1:4]))
                    positions.append((atom_type, position))
                except ValueError:
                    break  # 如果转换失败，停止读取位置（可能是遇到了其他部分）

    return np.array(lattice_vectors), positions

def write_structure(file_path, lattice_vectors, positions):
    with open(file_path, 'w') as file:
        file.write('28 1 14 2\n')
        file.write(' Lattice vector (Angstrom)\n')
        for vec in lattice_vectors:
            file.write('   {:16.8e}    {:16.8e}    {:16.8e}\n'.format(vec[0], vec[1], vec[2]))
        file.write('Position (normalized), move_x, move_y, move_z\n')
        for atom_type, pos in positions:
            file.write('   {:2d}      {:12.9f}      {:12.9f}      {:12.9f}     1  1  1\n'.format(atom_type, pos[0], pos[1], pos[2]))

def cart_to_frac(cart_coords, lattice_vectors):
    """将笛卡尔坐标转换为分数坐标"""
    inv_lattice = np.linalg.inv(lattice_vectors)
    frac_coords = np.dot(cart_coords, inv_lattice)
    return frac_coords

def frac_to_cart(frac_coords, lattice_vectors):
    """将分数坐标转换为笛卡尔坐标"""
    cart_coords = np.dot(frac_coords, lattice_vectors)
    return cart_coords

def apply_periodic_boundary_conditions(frac_coords):
    """应用周期性边界条件"""
    return frac_coords % 1

def is_within_distance(positions1, positions2, min_dist, max_dist):
    dist_matrix = cdist(positions1, positions2)
    return np.any((dist_matrix >= min_dist) & (dist_matrix <= max_dist))

def generate_random_translation(lattice_vectors, positions, min_dist, max_dist, num_attempts=10000):
    original_positions = np.array([pos for _, pos in positions])
    original_frac_positions = cart_to_frac(original_positions, lattice_vectors)
    
    for _ in range(num_attempts):
        # 随机生成新的分数坐标偏移
        x, y, z = np.random.rand(3) - 0.5
        new_frac_positions = original_frac_positions + [x, y, z]
        
        # 应用周期性边界条件
        new_frac_positions = apply_periodic_boundary_conditions(new_frac_positions)
        new_cart_positions = frac_to_cart(new_frac_positions, lattice_vectors)
        
        # 检查新位置是否满足距离约束
        if is_within_distance(original_positions, new_cart_positions, min_dist, max_dist):
            return [(atom_type, new_pos) for atom_type, new_pos in zip([p[0] for p in positions], new_cart_positions)]
    
    raise ValueError("Could not find a suitable translation within the specified number of attempts.")

# 读取原始EA分子结构
lattice_vectors, ea_positions = read_structure(input_file)

# 生成并保存新结构
successful_structures = 0
for i in range(num_structures):
    try:
        new_ea_positions = generate_random_translation(lattice_vectors, ea_positions, min_distance, max_distance, num_attempts)
        combined_positions = ea_positions + new_ea_positions
        output_file_path = os.path.join(output_folder, "structure_{}.pwmat".format(i+1))
        write_structure(output_file_path, lattice_vectors, combined_positions)
        print("已生成文件: {}".format(output_file_path))
        successful_structures += 1
    except ValueError as e:
        print("警告: 无法为结构 {} 找到合适的位置: {}".format(i+1, e))

print("处理完成。成功生成了 {} 个结构。".format(successful_structures))
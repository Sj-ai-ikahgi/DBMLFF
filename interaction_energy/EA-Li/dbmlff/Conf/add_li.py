# -*- coding: utf-8 -*-
import os
import numpy as np
from scipy.spatial import distance_matrix

# 定义文件路径和模式
input_file = './atom_EA.config'
output_folder = './new_structures'
num_structures = 200
min_distance = 0.5  # 埃
max_distance = 3.0  # 埃

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
        file.write('15 2 14 1 1 1\n')
        file.write(' Lattice vector (Angstrom)\n')
        for vec in lattice_vectors:
            # 使用 .format() 方法代替 f-string
            file.write('   {:16.8e}    {:16.8e}    {:16.8e}\n'.format(vec[0], vec[1], vec[2]))
        file.write('Position (normalized), move_x, move_y, move_z\n')
        for atom_type, pos in positions:
            # 使用 .format() 方法代替 f-string
            file.write('   {:2d}      {:12.9f}      {:12.9f}      {:12.9f}     1  1  1\n'.format(atom_type, pos[0], pos[1], pos[2]))

def is_within_distance(positions, new_pos, min_dist, max_dist):
    dists = np.linalg.norm(np.array([pos for _, pos in positions]) - new_pos, axis=1)
    return np.any((dists >= min_dist) & (dists <= max_dist))

def generate_random_position(lattice_vectors, positions, min_dist, max_dist, num_attempts=1000):
    a, b, c = lattice_vectors
    for _ in range(num_attempts):
        x, y, z = np.random.rand(3)
        new_pos = x * a + y * b + z * c
        if is_within_distance(positions, new_pos, min_dist, max_dist):
            return new_pos
    raise ValueError("Could not find a suitable position within the specified number of attempts.")

# 读取原始结构
lattice_vectors, positions = read_structure(input_file)

# 生成并保存新结构
for i in range(num_structures):
    try:
        new_pos = generate_random_position(lattice_vectors, positions, min_distance, max_distance)
        new_positions = positions + [(3, new_pos)]  # 假设锂原子的类型为 3
        output_file_path = os.path.join(output_folder, "structure_{}.pwmat".format(i+1))
        write_structure(output_file_path, lattice_vectors, new_positions)
        print("已生成文件: {}".format(output_file_path))
    except ValueError as e:
        print("警告: 无法为结构 {} 找到合适的位置: {}".format(i+1, e))

print("处理完成。")
import os

def read_config(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    config = {}
    current_section = None
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('Lattice vector'):
            current_section = 'lattice'
            config[current_section] = []
        elif line.startswith('Position'):
            current_section = 'position'
            config[current_section] = []
        elif line.startswith('Force'):
            current_section = 'force'
            config[current_section] = []
        elif line.startswith('Velocity'):
            current_section = 'velocity'
            config[current_section] = []
        else:
            if current_section:
                config[current_section].append(line)
    
    return config

def write_config(config, filename):
    with open(filename, 'w') as file:
        file.write("10     1  10   1 \n")
        file.write("Lattice vector (Angstrom)\n")
        for vec in config['lattice']:
            file.write(f"{vec}\n")
        file.write("Position (normalized), move_x, move_y, move_z\n")
        for pos in config['position']:
            file.write(f"{pos}\n")
        file.write("Force (-force, eV/Angstrom)\n")
        for force in config['force']:
            file.write(f"{force}\n")
        file.write("Velocity (bohr/fs)\n")
        for vel in config['velocity']:
            file.write(f"{vel}\n")

def modify_lattice_and_position(config, new_box_size):
    # Modify lattice vectors to form a cube with side length new_box_size
    config['lattice'] = [
        f"{new_box_size}  0.000000000  0.000000000",
        f"0.000000000 {new_box_size}  0.000000000",
        f"0.000000000  0.000000000 {new_box_size}"
    ]
    
    # Convert positions from normalized to absolute coordinates
    positions = []
    for pos_line in config['position']:
        tokens = pos_line.split()
        atom_type = tokens[0]
        x_norm = float(tokens[1])
        y_norm = float(tokens[2])
        z_norm = float(tokens[3])
        x_abs = x_norm * new_box_size
        y_abs = y_norm * new_box_size
        z_abs = z_norm * new_box_size
        positions.append((atom_type, x_abs, y_abs, z_abs))
    
    # Calculate the center of mass of the molecule
    total_mass = 0.0
    com_x = 0.0
    com_y = 0.0
    com_z = 0.0
    
    # Assuming all atoms have the same mass for simplicity
    mass = 1.0
    
    for _, x, y, z in positions:
        total_mass += mass
        com_x += x * mass
        com_y += y * mass
        com_z += z * mass
    
    com_x /= total_mass
    com_y /= total_mass
    com_z /= total_mass
    
    # Translate the molecule to the center of the box
    translation_x = new_box_size / 2 - com_x
    translation_y = new_box_size / 2 - com_y
    translation_z = new_box_size / 2 - com_z
    
    new_positions = []
    for atom_type, x, y, z in positions:
        new_x = x + translation_x
        new_y = y + translation_y
        new_z = z + translation_z
        new_positions.append(f"{atom_type} {new_x/new_box_size:.15f} {new_y/new_box_size:.15f} {new_z/new_box_size:.15f} 1 1 1")
    
    config['position'] = new_positions

def batch_process_configs(input_dir, new_box_size=11.4):
    for i in range(1, 1501):
        input_filename = os.path.join(input_dir, f"DFT_{i}.config")
        if not os.path.isfile(input_filename):
            print(f"File not found: {input_filename}")
            continue
        
        try:
            config = read_config(input_filename)
            modify_lattice_and_position(config, new_box_size)
            write_config(config, input_filename)
            print(f"Processed and saved {input_filename}")
        except Exception as e:
            print(f"An error occurred while processing {input_filename}: {e}")

if __name__ == "__main__":
    input_directory = "./Conf_1200K"
    batch_process_configs(input_directory)




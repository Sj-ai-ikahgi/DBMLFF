

def generate_movement_piece(
    num_atom, 
    lattice, 
    position, 
    nonperiodic_position=None, 
    force=None, 
    velocity=None, 
    atomic_energy=None, 
    de=None, 
):
    movement_piece = ''
    
    movement_piece += '%9d atoms,Iteration\n' % (num_atom)
    movement_piece += '\n'
    
    movement_piece += ' Lattice vector (Angstrom)\n'
    if len(lattice) != 3:
        print('generate_movement.generate_movement_piece: ERROR!!! len(lattice) != 3, EXIT!')
        exit()
    if type(lattice[0]) == str:
        for i in lattice:
            movement_piece += i
    else:
        for i in lattice:
            movement_piece += '%23.16E %23.16E %23.16E\n' % tuple(i)
    
    movement_piece += ' Position (normalized), move_x, move_y, move_z\n'
    if len(position) != num_atom:
        print('generate_movement.generate_movement_piece: ERROR!!! len(position) != num_atom, EXIT!')
        exit()
    if type(position[0]) == str:
        for i in position:
            movement_piece += i
    else:
        for i in position:
            movement_piece += '%4d%26.16f%26.16f%26.16f%6d%3d%3d\n' % tuple(i)
    
    movement_piece += ' nonperiodic_Position (normalized), move_x, move_y, move_z\n'
    if nonperiodic_position == None:
        for i in range(num_atom):
            movement_piece += '\n'
    else:
        if len(nonperiodic_position) != num_atom:
            print('generate_movement.generate_movement_piece: ERROR!!! len(nonperiodic_position) != num_atom, EXIT!')
            exit()
        if type(nonperiodic_position[0]) == str:
            for i in nonperiodic_position:
                movement_piece += i
        else:
            for i in nonperiodic_position:
                movement_piece += '%4d%26.16f%26.16f%26.16f%6d%3d%3d\n' % tuple(i)
    
    movement_piece += ' Force (-force, eV/Angstrom)\n'
    if force == None:
        for i in range(num_atom):
            movement_piece += '\n'
    else:
        if len(force) != num_atom:
            print('generate_movement.generate_movement_piece: ERROR!!! len(force) != num_atom, EXIT!')
            exit()
        if type(force[0]) == str:
            for i in force:
                movement_piece += i
        else:
            for i in force:
                movement_piece += '%4d%26.16f%26.16f%26.16f\n' % tuple(i)
    
    movement_piece += ' Velocity (bohr/fs)\n'
    if velocity == None:
        for i in range(num_atom):
            movement_piece += '\n'
    else:
        if len(velocity) != num_atom:
            print('generate_movement.generate_movement_piece: ERROR!!! len(velocity) != num_atom, EXIT!')
            exit()
        if type(velocity[0]) == str:
            for i in velocity:
                movement_piece += i
        else:
            for i in velocity:
                movement_piece += '%4d%26.16f%26.16f%26.16f\n' % tuple(i)
    
    de = 0.0 if de == None else de
    movement_piece += 'Atomic-Energy, Etot(eV),E_nonloc(eV),Q_atom:dE(eV)=  %23.16E\n' % (de)
    if atomic_energy == None:
        for i in range(num_atom):
            movement_piece += '\n'
    else:
        if len(atomic_energy) != num_atom:
            print('generate_movement.generate_movement_piece: ERROR!!! len(atomic_energy) != num_atom, EXIT!')
            exit()
        if type(atomic_energy[0]) == str:
            for i in atomic_energy:
                movement_piece += i
        else:
            for i in atomic_energy:
                movement_piece += '%4d %23.16E %23.16E %23.16E\n' % tuple(i)
    
    movement_piece += ' -------------------------------------------------\n'
    
    return movement_piece


if __name__ == '__main__':
    pass

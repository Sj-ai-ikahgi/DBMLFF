import sys
import os
import time
from copy import deepcopy
import numpy as np
import re
from lammps import lammps

# 设置LAMMPS库路径
os.environ['LD_LIBRARY_PATH'] = "/share/app/lammps-stable_2Aug2023_update2/bin:" + os.environ.get('LD_LIBRARY_PATH', '')

# ========================================================
# 辅助函数定义 (针对MOVEMENT文件格式)
# ========================================================

def Tri2Orth(a, b, c):
    """
    将三个晶格向量转换为正交晶格矩阵
    """
    lattice = np.zeros((3, 3))
    lattice[0] = a
    lattice[1] = b
    lattice[2] = c
    return lattice

def Lamda2Cart(lattice, coords):
    """
    将分数坐标转换为笛卡尔坐标
    """
    cart_coords = []
    for coord in coords:
        cart_coord = np.dot(lattice.T, coord)
        cart_coords.append(cart_coord.tolist())
    return cart_coords

def Move_Atom_Together(natoms, coords, center_atom=0):
    """
    移动原子使其聚集在一起
    """
    if natoms == 0:
        return [0, 0, 0], []
    
    center_coord = np.array(coords[center_atom])
    new_coords = []
    for coord in coords:
        new_coord = np.array(coord) - center_coord
        new_coords.append(new_coord.tolist())
    return center_coord.tolist(), new_coords

def New_MOVEMENT_Conf(isfile):
    """
    读取MOVEMENT文件的新配置头
    """
    while True:
        line = isfile.readline()
        if not line:
            return ''
        
        # 匹配格式: "   1664 atoms,Iteration (fs) =    0.0, Etot,Ep,Ek (eV) = 0.1 0.1 0.0, SCF =         -1"
        if "atoms,Iteration" in line:
            return line.strip()
    
    return ''

def Get_One_MOVEMENT_Conf(headline, isfile):
    """
    从MOVEMENT文件中读取一个完整的配置 - 修正版
    """
    # 解析标题行
    match = re.match(r'(\d+)\s+atoms,Iteration \(fs\)\s*=\s*([\d.]+), Etot,Ep,Ek \(eV\)\s*=\s*([\d.]+)\s+([\d.]+)\s+([\d.]+), SCF\s*=\s*([\d-]+)', headline)
    if not match:
        raise ValueError(f"无法解析标题行: '{headline}'")
    
    natoms = int(match.group(1))
    ntimestep = float(match.group(2))
    Etot = float(match.group(3))
    Ep = float(match.group(4))
    Ek = float(match.group(5))
    scf = int(match.group(6))
    
    # 跳过"Lattice vector (Angstrom)"行
    lattice_header = isfile.readline().strip()
    if not lattice_header or "Lattice" not in lattice_header:
        raise ValueError(f"预期晶格标题行，但得到: '{lattice_header}'")
    
    # 读取晶格向量
    lattice = []
    for i in range(3):
        line = isfile.readline().strip()
        if not line:
            raise ValueError("晶格向量数据不完整")
        lattice.append([float(x) for x in line.split()[:3]])
    
    # 跳过原子信息标题行
    atom_header = isfile.readline().strip()
    if not atom_header or "Position" not in atom_header:
        raise ValueError(f"预期原子标题行，但得到: '{atom_header}'")
    
    # 读取原子信息
    eindex = []
    coords = []
    forces = []
    velocity = []
    
    for i in range(natoms):
        line = isfile.readline().split()
        if len(line) < 7:
            raise ValueError(f"原子行数据不完整: {line}")
        
        eindex.append(int(line[0]))
        coords.append([float(line[1]), float(line[2]), float(line[3])])
        # 注意: 这里的"move_x, move_y, move_z"实际上是力分量
        forces.append([float(line[4]), float(line[5]), float(line[6])])
        # MOVEMENT文件没有速度信息，填充零值
        velocity.append([0.0, 0.0, 0.0])
    
    return ntimestep, Ep, natoms, lattice, eindex, coords, forces, velocity, None

def eindex2atype(eindex):
    """
    将元素索引转换为原子类型
    """
    unique_eindex = sorted(set(eindex))
    elem_list = [f"E{e}" for e in unique_eindex]
    eindex_list = unique_eindex
    emass_list = [1.0] * len(unique_eindex)  # 简化质量
    natype = len(unique_eindex)
    atype = [unique_eindex.index(e) + 1 for e in eindex]  # 原子类型从1开始
    return elem_list, eindex_list, emass_list, natype, atype

def Write_LAMMPS_Dump(confId, fname, natoms, lattice, atype, coords):
    """
    写入LAMMPS dump文件
    """
    with open(fname, 'a') as f:
        f.write(f"ITEM: TIMESTEP\n{confId}\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{natoms}\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        # 注意: LAMMPS需要3x3晶格矩阵的对角线元素
        f.write(f"0.0 {lattice[0][0]}\n")
        f.write(f"0.0 {lattice[1][1]}\n")
        f.write(f"0.0 {lattice[2][2]}\n")
        f.write("ITEM: ATOMS id type x y z\n")
        for i in range(natoms):
            f.write(f"{i+1} {atype[i]} {coords[i][0]} {coords[i][1]} {coords[i][2]}\n")

# ========================================================
# LAMMPS 接口类
# ========================================================

class sLMP:

    def __init__(self, cutoff=9.0):
        # 使用正确的库加载方式
        args = ["-log", "none", "-screen", "none"]
        try:
            # 尝试加载默认的LAMMPS库
            self.lmp = lammps(cmdargs=args)
        except:
            # 如果失败，尝试指定库名称
            try:
                self.lmp = lammps(name="lammps", cmdargs=args)
            except:
                # 如果仍然失败，尝试直接指定库路径
                libpath = "/share/app/lammps-stable_2Aug2023_update2/bin/liblammps.so"
                self.lmp = lammps(libname=libpath, cmdargs=args)
        self.cutoff = cutoff

    def Setup(self):
        self.lmp.command("newton on")
        self.lmp.command("units metal")
        self.lmp.command("boundary p p p")
        self.lmp.command("atom_style atomic")
        self.lmp.command("read_data lammps.data")
        self.lmp.command("pair_style lj/cut " + str(self.cutoff))
        self.lmp.command("pair_coeff * * 1.0 1.0")
        self.lmp.command("neighbor 0.0 bin")
        self.lmp.command("thermo_style custom step etotal")

    def LAMMPS_Load_New_Conf(self, n, aid, atype, x):
        self.lmp.command("delete_atoms group all")
        self.lmp.create_atoms(n, aid, atype, x)

    def LAMMPS_Data(self):
        self.lmp.command("run 0 post no")
        nlidx = self.lmp.find_pair_neighlist('lj/cut')
        nl = self.lmp.numpy.get_neighlist(nlidx)
        natoms = self.lmp.extract_global('natoms')
        nlocal = self.lmp.extract_global('nlocal')
        nghost = self.lmp.extract_global('nghost')
        x = self.lmp.extract_atom('x')
        tags = self.lmp.extract_atom('id')
        atype = self.lmp.extract_atom('type')
        nall = nlocal + nghost
        return nall, natoms, tags, atype, x, nl

    def LAMMPS_Finalize(self):
        self.lmp.finalize()
        return

# ========================================================
# 主配置处理类
# ========================================================

class Configuration:

    def __init__(self,
                 sysName,
                 infile_format,
                 infilename,
                 idir,
                 odir,
                 sconfId=1,
                 econfId=1e12,
                 cutoff=9.0,
                 molNum_list=[],
                 molSize_list=[]):

        self.sysName = sysName
        self.infile_format = infile_format
        self.infilename = infilename
        self.idir = idir
        self.odir = odir
        self.sconfId = sconfId
        self.econfId = econfId
        self.cutoff = cutoff
        self.molNum_list = molNum_list
        self.molSize_list = molSize_list
        self.isfile = ''
        self.confId = 0
        self.headline = []
        self.o2pId = []
        self.p2oId = []
        self.molId = []
        self.init_data()
        self.stime = time.time()
        self.slmp = sLMP(self.cutoff)
        
        # 创建输出目录
        if not os.path.isdir(self.odir):
            os.makedirs(self.odir, exist_ok=True)
            
        self.lammps_dump_fname = os.path.join(self.odir, f"{sysName}_{infile_format}_lammps_unwrapped.dump")
        self.lammps_dump_wrapped_fname = os.path.join(self.odir, f"{sysName}_{infile_format}_lammps_wrapped.dump")
        
        # 删除旧文件
        for fname in [self.lammps_dump_fname, self.lammps_dump_wrapped_fname]:
            if os.path.exists(fname):
                os.remove(fname)
                
        self.pre_coords = []
        self.unwrapped_coords = []

    def init_data(self):
        self.ntimestep = 0.0
        self.Ep = 0.0
        self.natoms = 0
        self.lattice = []
        self.eindex = []
        self.coords = []
        self.forces = []
        self.velocity = []

    def Get_atomId(self):
        self.o2pId = np.repeat(-1, self.natoms)
        self.p2oId = np.repeat(-1, self.natoms)
        self.molId = np.zeros([self.natoms], dtype=int)

        if self.infile_format == "DBMLFF":
            print("\n------ 处理 DBMLFF 格式文件 ------\n")
            self.o2pId = np.arange(0, self.natoms)
            self.p2oId = np.arange(0, self.natoms)
        else:
            print("\n------ 处理 PWMAT 格式文件 ------\n")
            index_file = os.path.join(self.idir, "ORIGIN.INDEX")
            if not os.path.exists(index_file):
                print(f"警告: 未找到 ORIGIN.INDEX 文件: {index_file}")
                print("使用默认原子顺序")
                self.o2pId = np.arange(0, self.natoms)
                self.p2oId = np.arange(0, self.natoms)
            else:
                with open(index_file, 'r') as tisfile:
                    while True:
                        token = tisfile.readline()
                        if not token:
                            break
                        token = token.split()
                        if len(token) < 5:
                            continue
                        try:
                            oid = int(token[0]) - 1
                            pid = int(token[4]) - 1
                            self.o2pId[oid] = pid
                            self.p2oId[pid] = oid
                        except:
                            print(f"警告: 解析 ORIGIN.INDEX 行失败: {token}")

        i = 0
        mid = 0
        for mtype in range(len(self.molNum_list)):
            for molNum in range(self.molNum_list[mtype]):
                mid += 1
                for aid in range(self.molSize_list[mtype]):
                    if i >= self.natoms:
                        print(f"警告: 原子索引超出范围: i={i}, 原子总数={self.natoms}")
                        break
                    pid = self.o2pId[i] if i < len(self.o2pId) else i
                    self.molId[pid] = mid
                    i += 1

    def Get_Mol(self, mid_list):
        teindex = []
        tcoords = []
        tforces = []
        tvelocity = []
        natoms = 0
        for mId in mid_list:
            for i in range(self.natoms):
                pid = self.o2pId[i] if i < len(self.o2pId) else i
                if self.molId[pid] == mId:
                    natoms += 1
                    teindex.append(self.eindex[pid])
                    tcoords.append(self.coords[pid])
                    tforces.append(self.forces[pid])
                    tvelocity.append(self.velocity[pid])
        return natoms, teindex, tcoords, tforces, tvelocity

    def Zero_Center(self, coords, lattice):
        natoms = len(coords)
        if natoms == 0:
            return
        
        center = np.zeros(3, dtype=float)
        for i in range(natoms):
            for j in range(3):
                center[j] += coords[i][j]
        for j in range(3):
            center[j] = center[j] / natoms
        for i in range(natoms):
            for j in range(3):
                coords[i][j] = coords[i][j] - center[j] + lattice[j][j] * 0.5

    def ion_Center(self):
        coords = []
        atype = []
        natoms = 0
        mId = 0
        
        # 根据系统名称定义离子类型
        if self.sysName == "80EC-10LiPF6-60DMC":
            # 对于80EC-10LiPF6-60DMC系统:
            # 类型1: Li+ (阳离子)
            # 类型2: PF6- (阴离子)
            ion_types = [1, 2]
        else:
            # 对于其他系统，保持原来的0和1为离子
            ion_types = [0, 1]
        
        for mtype in range(len(self.molNum_list)):
            for molNum in range(self.molNum_list[mtype]):
                mId += 1
                # 只处理离子类型的分子
                if mtype in ion_types:
                    tnatoms, teindex, tcoords, tforces, tvelocity = self.Get_Mol([mId])
                    if tnatoms > 0:  # 确保有原子
                        center_coord, new_coords = Move_Atom_Together(tnatoms, tcoords, 0)
                        # 添加类型标记：1表示阳离子，2表示阴离子
                        ion_type = 1 if mtype == ion_types[0] else 2
                        atype.append(ion_type)
                        coords.append(center_coord)
        natoms = len(coords)
        return natoms, atype, coords

    def Process_One_MOVEMENT_File(self):
        # 使用特定的文件名
        filename = "MOVEMENT_80EC-10LiPF6-60DMC"
        input_path = os.path.join(self.idir, filename)
        
        if not os.path.exists(input_path):
            print(f"错误: 未找到输入文件: {input_path}")
            sys.exit(1)
            
        print(f"处理文件: {input_path}")
        self.isfile = open(input_path, 'r')
        
        cand_list = np.arange(self.sconfId, self.econfId, 1)
        tconfId = 0
        
        while True:
            tconfId += 1
            try:
                self.headline = New_MOVEMENT_Conf(self.isfile)
                if self.headline == '':
                    print("已到达文件末尾")
                    break
                
                print(f"处理配置 {tconfId}: {self.headline[:50]}...")
                
                self.init_data()
                try:
                    self.ntimestep, self.Ep, self.natoms, self.lattice, self.eindex, \
                        self.coords, self.forces, self.velocity, _ = Get_One_MOVEMENT_Conf(self.headline, self.isfile)
                except Exception as e:
                    print(f"处理配置 {tconfId} 时出错: {str(e)}")
                    # 尝试跳过此配置继续
                    continue
                
                # 检查是否在指定范围内
                if tconfId < self.sconfId:
                    if tconfId % 100 == 0:
                        print(f'跳过 {tconfId}/{self.sconfId} - {self.ntimestep}')
                    continue
                elif tconfId > self.econfId:
                    print(f"达到结束配置 {self.econfId}")
                    break
                
                # 检查是否在候选列表中
                if tconfId not in cand_list:
                    continue
                
                self.confId += 1
                print(f"输出配置 {self.confId} (原始 {tconfId})")
                
                # 创建晶格矩阵
                lattice_mat = Tri2Orth(self.lattice[0], self.lattice[1], self.lattice[2])
                
                # 获取原子类型
                try:
                    elem_list, eindex_list, emass_list, natype, atom_types = eindex2atype(self.eindex)
                except:
                    print("计算原子类型时出错")
                    continue
                
                # 首次配置初始化
                if self.confId == 1:
                    try:
                        self.Get_atomId()
                        print("原子映射初始化完成")
                    except Exception as e:
                        print(f"初始化原子映射时出错: {str(e)}")
                        continue
                
                try:
                    # 获取离子中心
                    natoms, atype, coords = self.ion_Center()
                    if natoms == 0:
                        print("警告: 未找到离子中心")
                        continue
                        
                    wrapped_coords = deepcopy(coords)
                    
                    if tconfId == self.sconfId:
                        self.unwrapped_coords = deepcopy(coords)
                        print("初始化未包裹坐标")
                    else:
                        for i in range(natoms):
                            for j in range(3):
                                dr = coords[i][j] - self.pre_coords[i][j]
                                if abs(dr) > 0.5:
                                    dr = dr - np.sign(dr)
                                self.unwrapped_coords[i][j] += dr
                    
                    self.pre_coords = deepcopy(coords)
                    
                    # 转换为笛卡尔坐标
                    cart_coords = Lamda2Cart(lattice_mat, self.unwrapped_coords)
                    wrapped_cart_coords = Lamda2Cart(lattice_mat, wrapped_coords)
                    
                    # 写入LAMMPS dump文件
                    Write_LAMMPS_Dump(self.confId, self.lammps_dump_fname, natoms, lattice_mat, atype, cart_coords)
                    Write_LAMMPS_Dump(self.confId, self.lammps_dump_wrapped_fname, natoms, lattice_mat, atype, wrapped_cart_coords)
                    print(f"配置 {self.confId} 写入成功")
                except Exception as e:
                    print(f"处理配置 {tconfId} 时出错: {str(e)}")
                    continue
            except Exception as e:
                print(f"处理配置 {tconfId} 时发生严重错误: {str(e)}")
                break
        
        self.slmp.LAMMPS_Finalize()
        self.isfile.close()
        self.etime = time.time() - self.stime
        print(f"输出配置总数: {self.confId}")
        print(f"总用时: {self.etime:.2f} 秒")

# ========================================================
# 主程序
# ========================================================

if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("用法: python 脚本.py 系统名称 文件格式 输入目录 输出目录 起始配置ID 结束配置ID")
        print("示例: python pwmatMOVEMENT_MSD.py 20LiFSI-24DME-60TTE DBMLFF .. 输出目录 1 10000")
        sys.exit(1)
    
    sysName = sys.argv[1]
    infile_format = sys.argv[2]
    idir = sys.argv[3]  # 输入目录 (例如 .. 表示上一级目录)
    odir = sys.argv[4]  # 输出目录
    sconfId = int(sys.argv[5])  # 起始配置ID
    econfId = int(sys.argv[6])  # 结束配置ID
    
    # 分子系统设置
    if sysName == "20LiFSI-24DME-60TTE":
        # 分子数量列表: [Li+, FSI-, DME, TTE]
        molNum_list = [20, 20, 24, 60]
        # 分子大小列表: [Li+, FSI-, DME, TTE]
        molSize_list = [1, 9, 16, 18]
    elif sysName == "EC-LiPF6":
        molNum_list = [160, 8, 8, 0]
        molSize_list = [1, 9, 16, 18]  # 可能需要调整
    elif sysName == "EC-DEC-LiPF6":
        molNum_list = [10, 0, 0, 10]
        molSize_list = [1, 9, 16, 18]  # 可能需要调整
    # 添加新系统
    elif sysName == "80EC-10LiPF6-60DMC":
        # 分子数量列表: [EC, Li+, PF6-, DMC]
        molNum_list = [80, 10, 10, 60]
        # 分子大小列表: [EC, Li+, PF6-, DMC]
        molSize_list = [10, 1, 7, 12]  # 原子数量
    else:
        print(f"未知系统名称: {sysName}")
        sys.exit(1)
    
    # 创建配置处理器
    conf = Configuration(
        sysName=sysName,
        infile_format=infile_format,
        infilename="MOVEMENT_80EC-10LiPF6-60DMC",
        idir=idir,
        odir=odir,
        sconfId=sconfId,
        econfId=econfId,
        molNum_list=molNum_list,
        molSize_list=molSize_list
    )
    
    # 处理MOVEMENT文件
    conf.Process_One_MOVEMENT_File()
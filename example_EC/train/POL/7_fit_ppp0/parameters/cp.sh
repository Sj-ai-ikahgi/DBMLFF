#!/bin/bash

whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"

# 复制文件
cp -p ${whome}/POL/POL_dataset/polar_param.input_ppp0 ./
cp -p ${whome}/MLFF/MLFF_dataset/parameters/MOVEMENT.type ./
cp -p ${whome}/MLFF/3_fit_mlff/fit_dir/fread_dfeat/fit_linearMM.input ./

# 修改 fit_linearMM.input 文件内容
# 使用 sed 命令进行精确替换

# 替换第一行：将 "1, 0, 1" 改为 "0,0,0,0"
sed -i '1s/\(5,\s*10,\s*10,\s*\)1, 0, 1\(.*\)/\10,0,0,0\2/' fit_linearMM.input

# 替换第六行（实际文件中的第7行）：修改权重参数
# 将 "0.0, 2.0, 0.01, 0.00001" 改为 "0.0, 1.0, 0.00, 0.00001"
sed -i '7s/0.0, 2.0, 0.01, 0.00001/0.0, 1.0, 0.00, 0.00001/' fit_linearMM.input

echo "文件修改完成！"
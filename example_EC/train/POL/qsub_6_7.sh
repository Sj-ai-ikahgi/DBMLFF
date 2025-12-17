#!/bin/bash

# 设置变量（请根据实际情况修改这两个值）
whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
molType="EC"

pwd=`pwd`

echo "----- polarization of single module -----"
cd ./6_generate_ppp0
bash qsub.sh
echo "-----Post-job execution: bash ./6_generate_ppp0/set_job_array.sh-----"  
cd ${pwd}

echo "----- 正在配置工作路径和分子类型 -----"
echo "工作路径: $whome"
echo "分子类型: $molType"

echo "----- 修改 Python 脚本中的配置 -----"
# 修改第一个Python脚本
if [ -f "./7_fit_ppp0/parameters/IN.NEIGHBORE_3/generate_in.neighbore_for_all_atom.py" ]; then
    # 使用sed命令修改whome变量
    sed -i "s|whome = \".*\"|whome = \"$whome\"|g" ./7_fit_ppp0/parameters/IN.NEIGHBORE_3/generate_in.neighbore_for_all_atom.py
    # 使用sed命令修改molType变量
    sed -i "s|molType = \".*\"|molType = \"$molType\"|g" ./7_fit_ppp0/parameters/IN.NEIGHBORE_3/generate_in.neighbore_for_all_atom.py
    echo "已修改 generate_in.neighbore_for_all_atom.py"
else
    echo "警告: ./7_fit_ppp0/parameters/IN.NEIGHBORE_3/generate_in.neighbore_for_all_atom.py 不存在"
fi

# 修改第二个Python脚本
if [ -f "./7_fit_ppp0/parameters/IN.NEIGHBORE_3/generate_in.neighbore_for_each_atom.py" ]; then
    # 使用sed命令修改whome变量
    sed -i "s|whome = \".*\"|whome = \"$whome\"|g" ./7_fit_ppp0/parameters/IN.NEIGHBORE_3/generate_in.neighbore_for_each_atom.py
    # 使用sed命令修改molType变量
    sed -i "s|molType = \".*\"|molType = \"$molType\"|g" ./7_fit_ppp0/parameters/IN.NEIGHBORE_3/generate_in.neighbore_for_each_atom.py
    echo "已修改 generate_in.neighbore_for_each_atom.py"
else
    echo "警告: ./7_fit_ppp0/parameters/IN.NEIGHBORE_3/generate_in.neighbore_for_each_atom.py 不存在"
fi

# 修改workflow脚本
if [ -f "./7_fit_ppp0/fit/template/workflow_ppp_model_fitting_ht.py" ]; then
    # 使用sed命令修改whome变量
    sed -i "s|whome = \".*\"|whome = \"$whome\"|g" ./7_fit_ppp0/fit/template/workflow_ppp_model_fitting_ht.py
    echo "已修改 workflow_ppp_model_fitting_ht.py"
else
    echo "警告: ./7_fit_ppp0/fit/template/workflow_ppp_model_fitting_ht.py 不存在"
fi

echo "----- 配置文件拷贝 -----"
cd ./7_fit_ppp0/parameters
bash cp.sh
cd ${pwd}

echo "----- 生成 NEIGHBORE 文件 -----"
cd ./7_fit_ppp0/parameters/IN.NEIGHBORE_3
python generate_in.neighbore_for_all_atom.py
python generate_in.neighbore_for_each_atom.py
cd ${pwd}

echo "----- 开始拟合原子和键 -----"
cd ./7_fit_ppp0/fit
bash qsub.sh
cd ${pwd}

echo "----- 所有任务已完成 -----"
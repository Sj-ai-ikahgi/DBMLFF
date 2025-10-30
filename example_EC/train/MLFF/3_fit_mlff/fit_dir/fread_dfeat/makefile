
.PHONY: none lin clst all pca vdw

none:

all:
	make all -C/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit

lin: energyL.pred.1 forceL.pred.1

clst: energyC.pred.1 forceC.pred.1

pca: feat_PV.1

vdw: vdw_fitB.ntype

energyC.pred.1: /data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/calc_lin_cluster.r cluster_fitB.ntype
	/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/calc_lin_cluster.r
cluster_fitB.ntype: /data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/linear_clusterM.r feat_PV.1
	/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/linear_clusterM.r
energyL.pred.1: /data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/calc_lin_forceMM.r linear_fitB.ntype
	/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/calc_lin_forceMM.r
linear_fitB.ntype: /data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/fit_lin_forceMM.r feat_PV.1
	/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/fit_lin_forceMM.r
feat_PV.1: /data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/feat_collect_PCA.r feat_collect.in
	/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/feat_collect_PCA.r
vdw_fitB.ntype: /data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/fit_vdw.r feat_PV.1
	/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/mlff_src/fit/fit_vdw.r
	

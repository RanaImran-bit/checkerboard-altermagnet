#!/bin/bash
#SBATCH -J L12 
#SBATCH -p defq              # 提交到 默认的defq 队列 
#SBATCH -N 1                 # 使用2个节点
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1    # 每个进程占用一个 cpu 核心
#SBATCH -t 50000:00          # 任务最大运行时间是 500 分钟
#SBATCH --mem=3G   


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
# python3 wf.py
make clear 
make all
./CPMC.exe

#!/bin/bash
#$ -S /bin/bash
#
# MPI/CORES
#$ -pe mpi-ib 24
#
# JOB NAME
#$ -N mch_offline 
#
# EMAIL JOB RESULTS
#$ -M caleb.magoola@noaa.gov
#
# LOG FILES
#$ -o /home/juliac/OfflineFennel/outputs/
#$ -e /home/juliac/OfflineFennel/outputs/
#
# SETUP ENVIRONMENT
export PATH=$PATH:/software/miniconda3/bin:/software/bin
source /opt/sge/default/common/settings.sh
cd /home/juliac/OfflineFennel/
#
# SEND JOB
source activate mvapich2
export MV2_FORCE_HCA_TYPE=22
export MV2_HOMOGENEOUS_CLUSTER=1
export MV2_USE_SHMEM_COLL=1
export MV2_ALLGATHER_REVERSE_RANKING=1
export MV2_DEFAULT_MAX_SEND_WQE=1024
export MV2_DEFAULT_MAX_RECV_WQE=2048
export MV2_SRQ_SIZE=10240
export MV2_SRQ_MAX_SIZE=32767
export MV2_USE_RDMA_FAST_PATH=1
export MV2_NUM_RDMA_BUFFER=32768
export MV2_CM_RECV_BUFFERS=32768
export MV2_USE_MCAST=1
export MV2_MCAST_NUM_NODES_THRESHOLD=2
export MV2_USE_UD_HYBRID=1
export MV2_USE_RDMA_CM_MCAST=0
export MV2_USE_RDMA_CM=1
export MV2_ENABLE_AFFINITY=0

mpirun -np 24 ./coawstM External/offline_gom.in > mch_offline.log


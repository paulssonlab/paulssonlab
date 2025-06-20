CLUSTER_NAME="$(cloud-cluster-name "$1")"
status=$?
[ $status == 0 ] || exit $status
#CLUSTER_PROJECT="${CLUSTER_PROJECT:-paulssonlab}"
CLUSTER_ZONE="${CLUSTER_ZONE:-us-east4-b}"
CORE_MACHINE_TYPE="${CORE_MACHINE_TYPE:-g1-small}"
JUPYTER_MACHINE_TYPE="${JUPYTER_MACHINE_TYPE:-n2-custom-4-71680-ext}"
DASK_MACHINE_TYPE="${DASK_MACHINE_TYPE:-n2-custom-2-14848}"
JUPYTER_CPU="${JUPYTER_CPU:-2}"
JUPYTER_MEM="${JUPYTER_MEM:-60G}"
DASK_SCHEDULER_CPU="${DASK_SCHEDULER_CPU:-1.7}"
DASK_SCHEDULER_MEM="${DASK_SCHEDULER_MEM:-6.5G}"
DASK_WORKER_CPU="${DASK_WORKER_CPU:-0.904}"
DASK_WORKER_MEM="${DASK_WORKER_MEM:-6.17G}"
DASK_WORKER_THREADS="${DASK_WORKER_THREADS:-1}"

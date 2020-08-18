Follow the [Zero to JupyterHub instructions](https://zero-to-jupyterhub.readthedocs.io/en/latest/google/step-zero-gcp.html) to set up a Kubernetes cluster:
1. Install google-cloud-sdk. This is already done for you on O2.
2. Enable Kubernetes Engine API](https://console.cloud.google.com/apis/api/container.googleapis.com/overview). Ensure you are logged in to the paulssonlab@gmail.com Google account. This is already done. Depending on IAM permissions, you may be able to do this by running `gcloud services enable container.googleapis.com`.
3. `gcloud components install kubectl`

We follow [Helm 3-specific](https://zero-to-jupyterhub.readthedocs.io/en/latest/setup-jupyterhub/setup-helm3.html) and [general](https://zero-to-jupyterhub.readthedocs.io/en/latest/setup-jupyterhub/setup-jupyterhub.html#setup-jupyterhub) installation instructions.
gcloud container clusters create deploy-test --zone us-east4-b --num-nodes 3
# or if existing cluster
gcloud container clusters get-credentials CLUSTER_NAME
helm repo add dask https://helm.dask.org/
helm repo update
# DECIDE NAMESPACE
helm create namespace $NAMESPACE
helm install dask/dask --namespace=$namespace --name=name --set hhh.hhh.ram=16G
kubectl config set-context $(kubectl config current-context) --namespace ${NAMESPACE:-jhub}

```
kubectl scale --replicas=0 deployment/kube-dns-autoscaler --namespace=kube-system
kubectl scale --replicas=1 deployment/kube-dns --namespace=kube-system
```

#!/bin/sh
source="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tag="gcr.io/paulssonlab/cloning_webhooks"
if [ "$1" == "-d" ]; then
    docker build "$source" --tag "$tag"
    docker push "$tag"
else
    gcloud builds submit "$source" --ignore-file=.dockerignore --gcs-log-dir=gs://paulssonlab_cloudbuild_logs/logs --config=cloudbuild.yaml
fi
if [ $? == 0 ]; then
    gcloud run deploy cloning-webhooks --image "$tag" --allow-unauthenticated --platform managed --region=us-east4
else
    echo Build failed, aborting.
fi

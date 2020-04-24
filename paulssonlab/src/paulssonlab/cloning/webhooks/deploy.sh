#!/bin/sh
source="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
docker build "$source" --tag gcr.io/paulssonlab/cloning_webhooks
docker push gcr.io/paulssonlab/cloning_webhooks
gcloud run deploy cloning-webhooks --image gcr.io/paulssonlab/cloning_webhooks --allow-unauthenticated --platform managed --region=us-east4

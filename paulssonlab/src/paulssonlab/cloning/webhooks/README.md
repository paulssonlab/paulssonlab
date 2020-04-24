# Cloning automation webhooks

Code for automatically managing strains, plasmid maps, and primers whenever there are changes in Google Sheets or Google Drive.

# Installation

You will need to authenticate with the `paulssonlab@gmail.com` Google Cloud account. Ask Jacob, he will give you a `credentials.json` file. Then in a shell, run:
```
gcloud auth activate-service-account --key-file=path/to/credentials.json
```

After making changes to the webhook code, you will need to rebuild the Docker image and redeploy the image using Google Cloud Run. We will follow [these instructions](https://cloud.google.com/run/docs/building/containers) to rebuild the Docker image. If you are doing this for the first time, you will need to install Docker and run `gcloud auth configure-docker`. Then run `./deploy.sh`. Alternatively, you can do the build and deploy steps manually:

To build a local Docker image and upload it, run:
```
docker build . --tag gcr.io/paulssonlab/cloning_webhooks
docker push gcr.io/paulssonlab/cloning_webhooks
```

Then deploy with:
```
gcloud run deploy cloning-webhooks --image gcr.io/paulssonlab/cloning_webhooks --allow-unauthenticated --platform managed --region=us-east4
```

# One-time setup
Jacob has already done all of these steps, but they are recorded here for future reference. We mostly follow [these instructions](https://cloud.google.com/run/docs/mapping-custom-domains#command-line) for mapping a custom domain to Google Run.

1. Purchase `paulssonlab.com` domain through Google Domains (with the paulssonlab@gmail.com Google account).

2. Run `gcloud domains verify paulssonlab.com`. Make sure you are signed in as the paulssonlab@gmail.com Google account. Follow the instructions, selecting “Google Domains” as the registrar. It will prompt you to add a TXT record to the paulssonlab.com DNS. Keep this browser window open.

3. In a new browser window, go to [Google Domains](https://domains.google.com/m/registrar/paulssonlab.com). Make sure you are signed in as the paulssonlab@gmail.com Google account. Click the “DNS” section on the left, and add the TXT record following the instructions.

4. Wait a 5-10 minutes. Then go back to your open browser window with the Webmaster Central domain verification page and click verify. If this fails, wait a few more minutes and try verifying again. Note that you do not need to add the TXT record again.

5. Go to the Google Cloud Console and click on [Domain verification](https://console.developers.google.com/apis/credentials/domainverification?project=paulssonlab), under "API & Services." Add `paulssonlab.com`.

6. Run `gcloud beta run domain-mappings create --service cloning-webhooks --domain cloning.webhooks.paulssonlab.com --platform managed --region us-east4`. (Note that this did not actually work, so I had to do this using the web interface.)

7. Go to [Google Domains](https://domains.google.com/m/registrar/paulssonlab.com). Make sure you are signed in as the paulssonlab@gmail.com Google account. Click the “DNS” section on the left, and add the CNAME record with name `cloning.webhooks` and domain name `ghs.googlehosted.com.` (include the period).

## Contributors

- Jacob Quinn Shenker

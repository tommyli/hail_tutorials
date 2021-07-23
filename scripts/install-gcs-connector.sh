#!/bin/bash

set -euo pipefail
set -x

LOCAL_BIN="$HOME/.local/bin"
GCLOUD_CONFIG="$HOME/.config/gcloud"

if [ -f "$GCLOUD_CONFIG/application_default_credentials.json" ] && [ -f "$LOCAL_BIN/install-gcs-connector.py" ];
then
    sudo python3 "$HOME/.local/bin/install-gcs-connector.py" -k "$HOME/.config/gcloud/application_default_credentials.json"
else
    echo "Skipping installing GCS Connector."
fi

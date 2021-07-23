# Overview

Repo for trying out [Hail Tutorials](https://hail.is/docs/0.2/tutorials-landing.html).

## Prerequisites

* Complete
  [https://code.visualstudio.com/docs/remote/containers#_getting-started](https://code.visualstudio.com/docs/remote/containers#_getting-started)

## Getting Started

1. Clone this repo
2. Start VS Code, press F1, and select Remote-Containers: Open Folder in Container... Select the path you cloned this
   repo to. 

## devcontainer addons

Note the devcontainer (i.e. everything in `$PROJECT_DIR/.devcontainer`) was copied from [vscode-dev-containers Python
3](https://github.com/microsoft/vscode-dev-containers/tree/main/containers/python-3).  The following additional steps
were added/configured.

* Everything in requirements.txt, most notably [jupyter-lab](https://pypi.org/project/jupyterlab/),
  [hail](https://hail.is/docs/0.2/index.html) and [gcloud](https://pypi.org/project/gcloud/)
* openjdk-11-jdk
* [Google Cloud SDK](https://cloud.google.com/sdk/docs/install)

### Google Cloud Config

Google Cloud SDK is installed in the container but the easiest way to use it is to share your host's Google Cloud SDK
config as a Docker volume mount.  This of course assumes your host is already configured to access Google Cloud, i.e.
`$HOME/.config/gcloud` already exists in your host with correct configurations.  If not, then:
  1. Follow [Google Cloud SDK install]([Google Cloud SDK](https://cloud.google.com/sdk/docs/install)) instructions.
  2. Run `gcloud auth application-default login` and follow prompts

This is also a required step for Hail to [Reading from Google Cloud
Storage](https://hail.is/docs/0.2/cloud/google_cloud.html#reading-from-google-cloud-storage)

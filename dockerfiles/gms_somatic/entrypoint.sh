#!/bin/bash
set -e
# activate conda environment and let the following process take over
source /venv/bin/activate
exec "$@"
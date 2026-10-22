#!/usr/bin/env bash
set -euo pipefail

# The default installation stays CPU-compatible. To opt in to a CUDA 12.4
# environment, export CHLOROSCAN_INSTALL_CUDA=1 before creating rule envs.
if [[ "${CHLOROSCAN_INSTALL_CUDA:-0}" == "1" ]]; then
    conda install -y -c conda-forge -c nvidia -c pytorch \
        pytorch=2.5.1 pytorch-cuda=12.4
fi

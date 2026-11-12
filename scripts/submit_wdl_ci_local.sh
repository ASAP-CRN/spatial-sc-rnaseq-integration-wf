#!/usr/bin/env bash
# Run wdl-ci's 5 subcommands locally against Workbench/Cromwell, bypassing
# GitHub Actions.
#
# Prereqs
# -------
# pip install "git+https://github.com/dnastack/wdl-ci.git@v2.2.0"
#
# An env-file (see Usage) with 1Password op:// references for:
#   WALLET_URL, WALLET_CLIENT_ID, WALLET_CLIENT_SECRET, WORKBENCH_NAMESPACE,
#   WORKBENCH_EWES_URL, WORKBENCH_WORKFLOW_SERVICE_URL,
#   WORKBENCH_EWES_REFRESH_TOKEN, WORKBENCH_WORKFLOW_SERVICE_REFRESH_TOKEN
#   (+ WDL_CI_CUSTOM_TEST_WDL_DIR if using custom test WDLs)
#
# wdl-ci.config.json's filename is hardcoded by wdl-ci itself and is always
# read from the current working directory.
#
# Usage
# -----
#   scripts/submit_wdl_ci_local.sh <path-to-1password-env-file>

set -euo pipefail

ENV_FILE="${1:?Usage: $0 <path-to-1password-env-file>}"

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "$REPO_ROOT"

op run --env-file="$ENV_FILE" -- wdl-ci lint
op run --env-file="$ENV_FILE" -- wdl-ci detect-changes
op run --env-file="$ENV_FILE" -- wdl-ci submit
op run --env-file="$ENV_FILE" -- wdl-ci monitor --update-digests
op run --env-file="$ENV_FILE" -- wdl-ci cleanup

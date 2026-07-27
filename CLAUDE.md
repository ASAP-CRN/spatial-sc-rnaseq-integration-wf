# spatial-sc-rnaseq-integration-wf

WDL workflow repository for integrating spatial transcriptomics and single-cell RNA-seq data (ASAP-CRN).

## Repository Structure

```
.github/
  workflows/lint-test-workflows.yml   # CI: lint + test WDL on every PR
  pull_request_template.md
docker/                               # Placeholder for dockerfiles for workflow tasks
workflows/                            # Placeholder for WDL workflow definitions
wf-common/                            # Git submodule — shared ASAP utilities
spatial-sc-rnaseq-integration-wdl-ci-custom-test-dir/   # Placeholder for wdl-ci custom test WDLs
wdl-ci.config.json                    # wdl-ci task/test configuration
```

## Submodules

- `wf-common` → https://github.com/ASAP-CRN/wf-common.git

After cloning, initialise submodules:

```bash
git submodule update --init --recursive
```

## CI/CD

GitHub Actions runs `dnastack/wdl-ci@v2.1.0` on every pull request.
Required repository secrets: `WALLET_URL`, `WALLET_CLIENT_ID`, `WALLET_CLIENT_SECRET`,
`WORKBENCH_NAMESPACE`, `WORKBENCH_EWES_URL`, `WORKBENCH_WORKFLOW_SERVICE_URL`,
`WORKBENCH_EWES_REFRESH_TOKEN`, `WORKBENCH_WORKFLOW_SERVICE_REFRESH_TOKEN`.

## WDL Validation (before opening a PR)

```bash
miniwdl check workflows/<your_workflow>.wdl
womtool validate workflows/<your_workflow>.wdl
```

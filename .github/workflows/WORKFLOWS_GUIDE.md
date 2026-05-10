# GitHub Workflows Guide

This document describes all workflows currently present in `.github/workflows` and what each one does.

## Repository Variables

Several workflows rely on repository-level variables:

Settings -> Secrets and variables -> Actions -> Variables

| Variable | Purpose |
|---|---|
| `PROJECT_PREFIX` | Issue/branch/commit prefix (for example `EFSA`). |
| `IGNORE_PREFIX` | Escape hatch prefix (for example `no-issue`) that skips some naming/notification checks. |

## Workflow Overview

| File | Workflow Name | Category | Trigger |
|---|---|---|---|
| `check-branch-name.yaml` | `CHECK: Branch Naming` | Validation | PRs targeting `main` |
| `check-commit-message.yaml` | `CHECK: Commit Message` | Validation | PRs targeting `main` |
| `auto-issue-prefix.yaml` | `AUTOMATION: Issue Prefixer` | Automation | Issue opened |
| `auto-branch-issue-tracking.yaml` | `AUTOMATION: Branch Issue Linker` | Automation | Push to `feature/**`, `bugfix/**`, `docs/**` |
| `auto-pr-open-notify.yml` | `AUTOMATION: PR Open Notification` | Automation | PR opened |
| `auto-pr-merged-notify.yaml` | `AUTOMATION: PR Merged Notification` | Automation | PR closed (merged only) |
| `auto_nextflow_run.yaml` | `CI: Nextflow Integration Test` | CI | Push/PR (scoped branches + scoped paths), manual dispatch |
| `deploy-docs.yml` | `Deploy Documentation to GitHub Pages` | Deployment | Push to `main`, manual dispatch |
| `static.yml` | `Deploy static content to Pages` | Deployment | Push to `main`, manual dispatch |

## Validation Workflows (Merge Blocking)

These workflows are intended to enforce standards on PRs into `main`.

### `check-branch-name.yaml`

- Valid branch patterns:
  - `feature/PROJECT_PREFIX-<number>_<description>`
  - `bugfix/PROJECT_PREFIX-<number>_<description>`
  - `docs/PROJECT_PREFIX-<number>_<description>`
- Branches starting with `IGNORE_PREFIX` are explicitly allowed.
- Fails the check if branch name does not match allowed formats.

### `check-commit-message.yaml`

- Inspects non-merge commits in `main..HEAD`.
- Valid subject formats:
  - `PROJECT_PREFIX-<number>: <message>`
  - or starting with `IGNORE_PREFIX`
- Fails the check if any commit subject is invalid.

## Automation Workflows (Non-Blocking Helpers)

These workflows automate issue/PR communication and naming consistency.

### `auto-issue-prefix.yaml`

- On issue creation, updates title to:
  - `PROJECT_PREFIX-<issue_number>: <original title>`
- Skips update if title already has the expected prefix.

### `auto-branch-issue-tracking.yaml`

- On push to `feature/**`, `bugfix/**`, or `docs/**`, parses branch name for issue ID.
- If matched, comments once on the related issue that the branch was created/pushed.
- Avoids duplicate comments by checking existing issue comments.

### `auto-pr-open-notify.yml`

- On PR opened, comments on linked issue with PR URL and source branch.
- Skips branches that start with `IGNORE_PREFIX`.
- Expects PR title to contain `<WORD>-<number>` token (for example `EFSA-123`).

### `auto-pr-merged-notify.yaml`

- On PR closed, runs only when PR is merged.
- Comments on linked issue that PR was merged.
- Skips branches that start with `IGNORE_PREFIX`.
- Expects PR title to contain `<WORD>-<number>` token.

## CI Workflow

### `auto_nextflow_run.yaml` (`CI: Nextflow Integration Test`)

- Purpose: run an integration test of the Nextflow pipeline.
- Triggers:
  - Push to `main` and `feature/EFSA-268_nfx_integration_v2` with path filters.
  - PRs targeting those branches with same path filters.
  - Manual dispatch.
- Key behavior:
  - Sets up Java 17 and Nextflow.
  - Generates `data/inputs/config.json` for test run.
  - Pre-pulls Docker images inferred from selected process names in `nextflow.config`.
  - Runs `nextflow run main.nf -profile test --max_cpu 1`.
  - Verifies expected outputs and report presence.
  - Uploads artifacts (`report.html`, `.nextflow.log`) for 14 days.

## Deployment Workflows (GitHub Pages)

### `deploy-docs.yml`

- Builds MkDocs site (`mkdocs build`) and deploys `./site` to GitHub Pages.
- Two-job pipeline: `build` then `deploy`.

### `static.yml`

- Deploys the repository root (`.`) as a static Pages artifact.
- Single deploy job with `actions/configure-pages`.
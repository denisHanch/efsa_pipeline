# Workflow Structure and Guidelines

This document outlines the GitHub Actions workflows used within the repository, differentiating between validation checks and automation tasks, along with their intended purpose and rules.

## Repository Variables

The workflows rely on repository variables configured in:

Settings → Secrets and variables → Actions → Variables

Current values used in this repository:

| Variable | Value |
|---------|------|
| PROJECT_PREFIX | EFSA |
| IGNORE_PREFIX | no-issue |

## Validation Checks
Files are prefixed with **check**. Validation workflows enforce repository conventions and standards. They run on every Pull Request targeting the `main` branch.

**Important:** These workflows block merging PRs if checks fail.

| Workflow                 | Trigger                             | Rules Enforced / Purpose |
|--------------------------|-------------------------------------|--------------------------|
| **Branch Naming Check**  | `on: pull_request (branches: main)` | Branch names must follow `feature/PROJECT_PREFIX-#_branch_name`, `bugfix/...`, or `docs/...`; branches starting with `IGNORE_PREFIX` are allowed and skipped. |
| **Commit Message Check** | `on: pull_request (branches: main)` | Commit messages must follow `PROJECT_PREFIX-#ISSUE: Your commit message`, or start with `IGNORE_PREFIX`. |

## Automation Workflows
Files are prefixed with **auto**. Automation workflows handle routine tasks like notifications, issue linking, and repository organization. They are not intended to enforce merge-blocking validation rules.

| Workflow                   | Trigger                                 | Automated Task / Purpose |
|---------------------------|------------------------------------------|--------------------------|
| **Issue Prefixer**         | `issues [opened]`                        | Automatically prefixes issue titles with `PROJECT_PREFIX-#ISSUE`. |
| **Branch-Issue Linking**   | `push [feature/**, bugfix/**, docs/**]`  | Comments on the related issue when an associated branch is created or pushed. |
| **PR Opened Notification** | `pull_request [opened]`                  | Comments on the related issue when a PR is opened. |
| **PR Merged Notification** | `pull_request [closed]`                  | Comments on the related issue when a PR is merged. |
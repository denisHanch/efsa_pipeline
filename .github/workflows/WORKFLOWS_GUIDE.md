# Workflow Structure and Guidelines

This document outlines the GitHub Actions workflows used within the repository, differentiating explicitly between validation checks and automation tasks, along with their intended purpose and rules.  


## Validation Checks
***./checks/ folder***

Validation workflows explicitly ensure code quality, conventions, and standards. They run on every Pull Request targeting the main branch.  

**Important:** These workflows explicitly block merging PRs if checks fail. Make sure to fix reported issues promptly.

| Workflow                 | Trigger                             | Rules Enforced / Purpose                                                                |
| ------------------------ | ----------------------------------- | --------------------------------------------------------------------------------------- |
| **Branch Naming Check**  | `on: pull_request (branches: main)` | Branches must follow: `feature/EFSA-#ISSUE_branch_name`, `bugfix/...`, or `docs/...` |
| **Commit Message Check** | `on: pull_request (branches: main)` | Commits must follow: `EFSA-#ISSUE: Your Vommit Message`                                   |


## Automation Workflows (automation/ folder)
***./automation/ folder***

Automation workflows handle routine tasks like notifications, linking issues, and maintaining repository organization. They are explicitly prefixed with "Automation:" and **do not** block PR merges.

| Workflow                   | Trigger                                 | Automated Task / Purpose                                         |
| -------------------------- | --------------------------------------- | ---------------------------------------------------------------- |
| **Issue Prefixer**         | `issues [opened]`                       | Automatically prefixes issue titles with `PROJECT_PREFIX-#ISSUE` |
| **Branch-Issue Linking**   | `push [feature/**, bugfix/**, docs/**]` | Comments on issues when associated branches are created/pushed.  |
| **PR Opened Notification** | `pull_request [opened]`                 | Notifies issue when related PR is opened                         |
| **PR Merged Notification** | `pull_request [closed]`       | Notifies issue upon PR merge completion                          |


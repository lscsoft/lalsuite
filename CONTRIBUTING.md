# Contributing to LALSuite

This page outlines the recommended procedure for contributing changes to the LALSuite repository. Please familiarise yourself with [GitLab](https://wiki.ligo.org/DASWG/GitLigoOrg) and [git-lfs](https://wiki.ligo.org/DASWG/GitLFS#Install_the_git_LFS_client) and ensure your account is configured according to these instructions.

## Reporting Issues

When reporting issues, please include as much detail as possible to reproduce the error, including information about your operating system and the version of each (relevant) component of LALSuite.
If possible, please include a brief, self-contained code example that demonstrates the problem.

## Contributing code

All contributions to LALSuite code must be made using the fork and [merge request](https://git.ligo.org/help/user/project/merge_requests/index.md) [workflow](https://git.ligo.org/help/workflow/forking_workflow.md), which must then be reviewed by one of the project maintainers.

If you with to contribute new code, or changes to existing code, please follow the following development workflow.

### Make a fork (copy) of LALSuite

**You only need to do this once**

1. Go to the [LALSuite repository home page](https://git.ligo.org/lscsoft/lalsuite)
2. Click on the *Fork* button, that should lead you [here](https://git.ligo.org/lscsoft/lalsuite/fork/new)

If you can't see the *Fork* button, make sure that you are logged in by checking for your account profile photo in the top right-hand corner of the screen.

### Clone your fork

```bash
git clone git@git.ligo.org:<your-user-name>/lalsuite.git
```

### Updating your fork

If you already have a fork of LALSuite, and are starting work on a new project you can link your clone to the main (`upstream`) repository and pull in changes that have been merged since the time you created your fork, or last updated:

1. Link your fork to the main repository:

    ```bash
    cd lalsuite
    git remote add upstream git@git.ligo.org:lscsoft/lalsuite.git
    ```

2. Fetch new changes from the `upstream` repo

    ```bash
    git fetch upstream
    ```

### Creating a new feature branch

All changes should be developed on a feature branch, in order to keep them separate from other work, simplifying review and merge once the work is done.

To create a new feature branch:

```bash
git checkout -b my-new-feature upstream/master
```

### Hack away

1. Develop the changes you would like to introduce, using `git commit` to finalise a specific change.
   Ideally commit small units of change often, rather than creating one large commit at the end, this will simplify review and make modifying any changes easier.

    Commit messages should be clear, identifying which code was changed, and why.
   Common practice is to use a short summary line (<50 characters), followed by a blank line, then more information in longer lines.

2. Push your changes to the remote copy of your fork on https://git.ligo.org.

    ```bash
    git push origin my-new-feature
    ```
   **Note:** For the first `push` of any new feature branch, you will likely have to use the `-u/--set-upstream` option to `push` to create a link between your new branch and the `origin` remote:

    ```bash
    git push --set-upstream origin my-new-feature
    ```

### Open a merge request

When you feel that your work is finished, you should create a merge request to propose that your changes be merged into the main (`upstream`) repository.

After you have pushed your new feature branch to `origin`, you should find a new button on the [LALSuite repository home page](https://git.ligo.org/lscsoft/lalsuite/) inviting you to create a merge request out of your newly pushed branch.
You should click the button, and proceed to fill in the title and description boxes on the merge request page.

Once the request has been opened, one of the maintainers will assign someone to review the change.

## More Information

More information regarding the usage of GitLab can be found in the main GitLab [documentation](https://git.ligo.org/help/).

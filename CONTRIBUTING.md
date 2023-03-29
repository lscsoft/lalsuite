# Contributing to LALSuite

This page outlines the recommended procedure for contributing changes to the LALSuite repository. Please read the introduction to [GitLab on git.ligo.org](https://wiki.ligo.org/Computing/GitLigoOrg) before you start.

## Reporting issues

If you have `ligo.org` authentication, please report issues directly through GitLab. Otherwise, you can use the [service desk address](mailto:contact+lscsoft-lalsuite-1438-issue-@support.ligo.org) to send bug reports by e-mail.

In either case, please include as much detail as possible to reproduce the error, including information about your operating system and the version of each (relevant) component of LALSuite. If possible, please include a brief, self-contained code example that demonstrates the problem.

Note that when an issue is marked as *Confidential*, currently this means that most internal users will also not be able to see it, but only a small number of people with reporter, developer or maintainer status.

## Contributing code

All contributions to LALSuite code must be made using the fork and [merge request](https://git.ligo.org/help/user/project/merge_requests/index.html) [workflow](https://git.ligo.org/help/user/project/repository/forking_workflow.html), which must then be reviewed by one of the project maintainers.

If you wish to contribute new code, or changes to existing code, please follow this development workflow:

### Make a fork (copy) of LALSuite

You only need to do this *once*.

1. Go to the [LALSuite repository home page](https://git.ligo.org/lscsoft/lalsuite)
2. Click on the *Fork* button, that should lead you [here](https://git.ligo.org/lscsoft/lalsuite/-/forks/new)
3. Select the namespace that you want to create the fork in, this will usually be your personal namespace

If you can't see the *Fork* button, make sure that you are logged in by checking for your account profile photo in the top right-hand corner of the screen.

### Clone your fork

Make sure that you have installed and configured [Git-LFS](https://wiki.ligo.org/Computing/GitLFS#Install_the_git_LFS_client):

```bash
git lfs install
```

for the management of large files. This is required to successfully build and install your development fork.

Then, clone your fork with

```bash
git clone git@git.ligo.org:<namespace>/lalsuite.git
```

### Keeping your fork up to date

Link your clone to the main (`upstream`) repository so that you can `fetch` and `pull` changes, `merge` them with your clone, and `push` them to your fork. Do *not* make changes on your master branch.

1. Link your fork to the main repository:

    ```bash
    cd lalsuite
    git remote add upstream git@git.ligo.org:lscsoft/lalsuite.git
    ```

   You only need to do this *once*.

2. Update your `master` branch to track changes from upstream:

    ```bash
    git checkout master
    git fetch upstream
    git branch --set-upstream-to upstream/master
    git pull
    ```

   You only need to do this *once*.

3. Fetch new changes from the `upstream` repository, merge them with your master branch, and push them to your fork on `git.ligo.org`:

    ```bash
    git checkout master
    git pull
    git push origin master
    ```

4. You can see which remotes are configured using

   ```bash
   git remote -v
   ```

   If you have followed the instructions thus far, you should see four lines. Lines one and two begin with `origin` and reference your fork on git.ligo.org with both `fetch` and `push` methods. Lines three and four begin with `upstream` and refer to the main repository on `git.ligo.org` with both `fetch` and `push` methods.

### Making changes

All changes should be developed on a feature branch in order to keep them separate from other work, thus simplifying the review and merge once the work is complete. The workflow is:

1. Create a new feature branch configured to track the `master` branch:

   ```bash
   git checkout master
   git pull
   git checkout -b my-new-feature upstream/master
   ```

   This command creates the new branch `my-new-feature`, sets up tracking the `upstream` repository, and checks out the new branch. There are other ways to do these steps, but this is a good habit since it will allow you to `fetch` and `merge` changes from `upstream/master` directly onto the branch.

2. Develop the changes you would like to introduce, using `git commit` to finalise a specific change. Ideally commit small units of change often, rather than creating one large commit at the end, this will simplify review and make modifying any changes easier.

   Commit messages should be clear, identifying which code was changed, and why. Common practice is to use a short summary line (<50 characters), followed by a blank line, then more information in longer lines.

2. Push your changes to the remote copy of your fork on https://git.ligo.org. The first `push` of any new feature branch will require the `-u/--set-upstream` option to `push` to create a link between your new branch and the `origin` remote:

    ```bash
    git push --set-upstream origin my-new-feature
    ```

    Subsequent pushes can be made with just:

    ```bash
    git push
    ```

3. Keep your feature branch up to date with the `upstream` repository by doing:

   ```bash
   git checkout master
   git pull
   git checkout my-new-feature
   git rebase upstream/master
   git push -f origin my-new-feature
   ```

   This works if you created your branch with the `checkout` command above. If you forgot to add the `upstream/master` starting point, then you will need to dig deeper into git commands to get changes and merge them into your feature branch.

   If there are conflicts between `upstream` changes and your changes, you will need to resolve them before pushing everything to your fork.

### Adding yourself to the author list

If you've never contributed to LALSuite before, you'll also need to add your name to the author list by running
```bash
make update-authors
```
from the top level of the repository. You can further edit the `.mailmap` file to adjust how your name is presented; special characters (e.g. accents) are supported.

### Open a merge request

When you feel that your work is finished, you should create a merge request to propose that your changes be merged into the main (`upstream`) repository.

After you have pushed your new feature branch to `origin`, you should find a new button on the [LALSuite repository home page](https://git.ligo.org/lscsoft/lalsuite/) inviting you to create a merge request out of your newly pushed branch. (If the button does not exist, you can initiate a merge request by going to the `Merge Requests` tab on your fork website on `git.ligo.org` and clicking `New merge request`)

You should click the button, and proceed to fill in the title and description boxes on the merge request page. It is recommended that you check the box to *Remove source branch when merge request is accepted*; this will result in the branch being automatically removed from your fork when the merge request is accepted.

Once the request has been opened, one of the maintainers will assign someone to review the change. There may be suggestions and/or discussion with the reviewer. These interactions are intended to make the resulting changes better. The reviewer will merge your request.

Once the changes are merged into the upstream repository, you should remove the development branch from your clone using

```bash
git branch -d my-new-feature
```

A feature branch should *not* be repurposed for further development as this can result in problems merging upstream changes.

### Continuous integration

GitLab runs continuous integration (CI) pipelines on LALSuite to ensure that it builds and passes its test suite on a wide variety of platforms. There are 2 main CI pipelines:

1. The push CI pipeline runs whenever you push commit(s) to your LALSuite fork. This pipeline performs some basic checks that LALSuite still builds and passes its tests with your changes:
   - each component LALSuite package can be build and installed in sequence (the so-called *package-level build*);
   - LALSuite can build all component packages at once (the so-called *top-level build*);
   - the LALSuite [Doxygen](https://doxygen.nl/) documentation can be built;
   - some basic checks for code style/formatting/whitespace errors, build byproduct files missing from `.gitignore`, etc.

2. The merge CI pipeline runs when you are ready to submit your changes to the upstream LALSuite fork via a merge request. This pipeline runs a much more comprehensive series of checks that LALSuite still builds and passes its tests with a wide variety of platforms (e.g. MacOS, various Linux distributions) and compilers (e.g. `clang`, `icc`, `gcc`). It also checks that LALSuite packages for a number of package management systems (e.g. RPM, Debian, Conda, Python wheels) are built correctly.

3. (A third CI pipeline runs nightly on the main [`lscsoft/lalsuite`](https://git.ligo.org/lscsoft/lalsuite) fork for deployment tasks, e.g. updating the [online documentation](https://lscsoft.docs.ligo.org/lalsuite/)).

You can request a subset of the jobs which normally run as part of the merge pipeline to also be run as part of the push pipeline. This is useful if you are making changes to LALSuite which could potentially cause problems with different platforms/compilers, or which could affect the packaging, and you want to test the effect of your changes before submitting a merge request.

For individual commits, you can request a subset of merge pipeline jobs to run by adding key text to the commit message, as listed in the table below. If instead you have a branch where you want a subset of merge pipeline jobs to be run on every push, you can name the branch to match one of the regular expressions given in the table below.

| If commit message contains | Or branch name matches    | Action                                                       |
| -------------------------- | ------------------------- | ------------------------------------------------------------ |
| `[ci compiler]`            | `/[-_]ci[-_]compiler/`    | Test different compilers (e.g. `clang`, `icc`, `gcc`)        |
| `[ci conda]`               | `/[-_]ci[-_]conda/`       | Build Conda packages                                         |
| `[ci coverage]`            | `/[-_]ci[-_]coverage/`    | Report test suite coverage                                   |
| `[ci debian]`              | `/[-_]ci[-_]debian/`      | Build Debian packages                                        |
| `[ci docker]`              | `/[-_]ci[-_]docker/`      | Build Docker containers                                      |
| `[ci docs]`                | `/[-_]ci[-_]docs/`        | Build the documentation                                      |
| `[ci full]`                | n/a                       | Run all jobs in the merge pipeline                           |
| `[ci integration]`         | `/[-_]ci[-_]integration/` | Longer-running integration tests, different<br/>top-level build configurations, etc. |
| `[ci koji]`                | `/[-_]ci[-_]koji/`        | Build RPM packages on a Koji server                          |
| `[ci lint]`                | `/[-_]ci[-_]lint/`        | Perform "lint" checks for code style/formatting/whitespace<br/>errors, build byproduct files missing from `.gitignore`, etc. |
| `[ci nightly]`             | n/a                       | Run all jobs in the nightly deployment pipeline [^1]         |
| `[ci pkg]`                 | `/[-_]ci[-_]pkg/`         | Perform a basic package-level build from tarballs            |
| `[ci platform]`            | `/[-_]ci[-_]platform/`    | Test different platforms (e.g. MacOS, various Linux distributions) |
| `[ci rhel]`                | `/[-_]ci[-_]rhel/`        | Build RPM packages                                           |
| `[ci tags]`                | n/a                       | Run all jobs in the pipeline for Git tags, e.g. `lalsuite-v*` [^1] |
| `[ci wheels]`              | `/[-_]ci[-_]wheels/`      | Build Python wheel packages                                  |

[^1]: The `[ci nightly]` and `[ci tags]` pipelines do not execute any deployment actions with external consequences, e.g. deploying documentation, pushing packages to repositories. These actions can only be executed by the third CI pipeline which runs nightly on the main `lscsoft/lalsuite` fork.

## More information

More information regarding the usage of GitLab can be found in the main GitLab [documentation](https://git.ligo.org/help/).

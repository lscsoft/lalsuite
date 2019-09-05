# Contributing to LALSuite

This page outlines the recommended procedure for contributing changes to the LALSuite repository. Please read the introduction to [GitLab on git.ligo.org](https://wiki.ligo.org/Computing/GitLigoOrg) before you start. 

## Reporting Issues

If you have ligo.org authentication, please report issues directly through gitlab.
Otherwise, you can use the service desk address
contact+lscsoft-lalsuite-1438-issue-@support.ligo.org
to send bug reports by e-mail.

In either case, please include as much detail as possible to reproduce the error, including information about your operating system and the version of each (relevant) component of LALSuite.
If possible, please include a brief, self-contained code example that demonstrates the problem.

Note that when an issue is marked as 'confidential',
currently this means that most internal users will also not be able to see it,
but only a small number of people with reporter, developer or maintainer status.

## Contributing code

All contributions to LALSuite code must be made using the fork and [merge request](https://git.ligo.org/help/user/project/merge_requests/index.md) [workflow](https://git.ligo.org/help/workflow/forking_workflow.md), which must then be reviewed by one of the project maintainers.

If you wish to contribute new code, or changes to existing code, please follow this development workflow:

### Make a fork (copy) of LALSuite

**You only need to do this once**

1. Go to the [LALSuite repository home page](https://git.ligo.org/lscsoft/lalsuite)
2. Click on the *Fork* button, that should lead you [here](https://git.ligo.org/lscsoft/lalsuite/forks/new)
3. Select the namespace that you want to create the fork in, this will usually be your personal namespace

If you can't see the *Fork* button, make sure that you are logged in by checking for your account profile photo in the top right-hand corner of the screen.

### Clone your fork

Make sure that you have installed and configured [git-lfs](https://wiki.ligo.org/Computing/GitLFS#Install_the_git_LFS_client) for the management of large files. This is required to successfully build and install your development fork. 

Then, clone your fork with 

```bash
git clone git@git.ligo.org:<namespace>/lalsuite.git
```

### Keeping your fork up to date

Link your clone to the main (`upstream`) repository so that you can `fetch` changes, `merge` them with your clone, and `push` them to your fork. Do *not* make changes on your master branch. 

1. Link your fork to the main repository:

    ```bash
    cd lalsuite
    git remote add upstream git@git.ligo.org:lscsoft/lalsuite.git
    ```

   You need only do this step once. 

2. Fetch new changes from the `upstream` repository, merge them with your master branch, and push them to your fork on git.ligo.org:

    ```bash
    git checkout master
    git fetch upstream
    git merge upstream/master
    git push
    ```

3. You can see which remotes are configured using

   ```bash
   git remote -v
   ```

   If you have followed the instructions thus far, you should see four lines. Lines one and two begin with `origin` and reference your fork on git.ligo.org with both `fetch` and `push` methods. Lines three and four begin with `upstream` and refer to the main repository on git.ligo.org with both `fetch` and `push` methods.

### Making changes

All changes should be developed on a feature branch in order to keep them separate from other work, thus simplifying the review and merge once the work is complete. The workflow is:

1. Create a new feature branch configured to track the `master` branch of the `upstream` repository:

   ```bash
   git checkout -b my-new-feature upstream/master
   ```

   This command creates the new branch `my-new-feature`, sets up tracking the `upstream` repository, and checks out the new branch. There are other ways to do these steps, but this is a good habit since it will allow you to `fetch` and `merge` changes from `upstream/master` directly onto the branch. 

2. Develop the changes you would like to introduce, using `git commit` to finalise a specific change.
   Ideally commit small units of change often, rather than creating one large commit at the end, this will simplify review and make modifying any changes easier.

   Commit messages should be clear, identifying which code was changed, and why.
   Common practice is to use a short summary line (<50 characters), followed by a blank line, then more information in longer lines.

2. Push your changes to the remote copy of your fork on https://git.ligo.org.
   The first `push` of any new feature branch will require the `-u/--set-upstream` option to `push` to create a link between your new branch and the `origin` remote:

    ```bash
    git push --set-upstream origin my-new-feature
    ```

    Subsequenct pushes can be made with 

    ```bash
    git push origin my-new-feature
    ```
   
3. Keep your feature branch up to date with the `upstream` repository by doing 

   ```bash
   git checkout my-new-feature
   git fetch upstream
   git rebase upstream/master
   git push -f origin my-new-feature
   ```

   This works if you created your branch with the `checkout` command above. If you forgot to add the `upstream/master` starting point, then you will need to dig deeper into git commands to get changes and merge them into your feature branch. 

   If there are conflicts between `upstream` changes and your changes, you will need to resolve them before pushing everything to your fork. 

### Open a merge request

When you feel that your work is finished, you should create a merge request to propose that your changes be merged into the main (`upstream`) repository.

After you have pushed your new feature branch to `origin`, you should find a new button on the [LALSuite repository home page](https://git.ligo.org/lscsoft/lalsuite/) inviting you to create a merge request out of your newly pushed branch. (If the button does not exist, you can initiate a merge request by going to the `Merge Requests` tab on your fork website on git.ligo.org and clicking `New merge request`)

You should click the button, and proceed to fill in the title and description boxes on the merge request page.
It is recommended that you check the box to `Remove source branch when merge request is accepted`; this will result in the branch being automatically removed from your fork when the merge request is accepted. 

Once the request has been opened, one of the maintainers will assign someone to review the change. There may be suggestions and/or discussion with the reviewer. These interactions are intended to make the resulting changes better. The reviewer will merge your request.

Once the changes are merged into the upstream repository, you should remove the development branch from your clone using 

```
git branch -d my-new-feature
```

A feature branch should *not* be repurposed for further development as this can result in problems merging upstream changes. 

### Special case: CW codes in LALPulsar and LALApps

For the continuous-wave (CW) related codes in the `lalpulsar/` and `lalapps/src/pulsar` directories,
issues should preferentially be reported [in the CW fork tracker](https://git.ligo.org/CW/software/lalsuite/issues).
Branches in that fork can also be used for short-term development, please also see the special contributing guide [here](https://git.ligo.org/CW/software/lalsuite/wikis/contributing).

## More Information

More information regarding the usage of GitLab can be found in the main GitLab [documentation](https://git.ligo.org/help/).

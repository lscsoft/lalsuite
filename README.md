# LALSuite

This is the main LALSuite development repository.

We now utilize [git-lfs](https://wiki.ligo.org/DASWG/GitLFS#Install_the_git_LFS_client) for the managament of large files and as such `git-lfs` needs to be installed and configured to correctly clone this repository. After installing `git-lfs` it can be configured using:

```
$ git lfs install
```

This only needs to be done once for each machine you access the repository. It can then be cloned using:

```
$ git clone git@git.ligo.org:lscsoft/lalsuite.git
```

## Contributing to LALSuite

The guide to [Contributing](https://git.ligo.org/lscsoft/lalsuite/blob/master/CONTRIBUTING.md) to LALSuite explains how to report issues and contribute fixes or new features using the fork and pull workflow. Please read and follow these directions.

## Notes on Ancient History

LALSuite was transferred to git.ligo.org in December 2017. Older history has been imported, though commit hashes were rewritten during the GitLFS conversion. Please note:

1. The "Original: " commit IDs quoted in each commit message can be used to compare with the [archived reference repo](https://git.ligo.org/lscsoft/lalsuite-archive), old issue discussions on the [Redmine tracker](https://bugs.ligo.org/redmine/projects/lalsuite), review wiki pages etc.

1. Commits before December 2017 may also include references to issues ("#number"). These refer to the corresponding [Redmine issue](https://bugs.ligo.org/redmine/projects/lalsuite) (LVC-authorized access only), and any clickable link the internal gitlab web interface produces for those old commits will be spurious.

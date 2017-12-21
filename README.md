# LALSuite

This is the main LALSuite development repository, past development is
captured in the repository below:

https://git.ligo.org/lscsoft/lalsuite-archive

This new repository utilizes [git-lfs](https://wiki.ligo.org/DASWG/GitLFS#Install_the_git_LFS_client) for the managament of large files and as such `git-lfs` needs to be installed and configured to correctly clone this repository. After installing `git-lfs` it can be configured using:

```
$ git lfs install
```

This only needs to be done once for each machine you access the repository. It can then be cloned using:

```
$ git clone git@git.ligo.org:lscsoft/lalsuite.git
```

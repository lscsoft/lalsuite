## Images

This repository provides Docker images with LALSuite installed. Two sets of images are provided:

* snapshot images of the latest development version of LALSuite, updated weekly.
* images of official LALSuite releases; the image of the most recent release is rebuilt weekly.

## Tags

Available tags are named following the schema below, with the following placeholders:

* `<base>`: base image LALSuite is installed into, either:
  * a generic operating system or platform:
    * `el`: Enterprise Linux (the LALSuite reference platform).
    * `debian`: Debian.
    * `cuda`: Ubuntu with CUDA support enabled.
  * a specific version of an operating system or platform; see the [Tags](./tags) list for available tags.
* `<X>`/`<Y>`/`<Z>`: LALSuite release major/minor/patch version.

### Snapshot images

| Tag          | Description                                                   |
| ------------ | ------------------------------------------------------------- |
| `dev`        | Snapshot image; base image is the LALSuite reference platform |
| `dev-<base>` | Snapshot image; base image is `<base>`                        |

### Release Images

| Tag                  | Description                                                     |
| -------------------- | --------------------------------------------------------------- |
| `latest`             | Latest release; base image is the LALSuite reference platform   |
| `latest-<base>`      | Latest release; base image is `<base>`                          |
| `<X>-<base>`         | Latest release with version `<X>.*.*`; base image is `<base>`   |
| `<X>.<Y>-<base>`     | Latest release with version `<X>.<Y>.*`; base image is `<base>` |
| `<X>.<Y>.<Z>-<base>` | Release with version `<X>.<Y>.<Z>`; base image is `<base>`      |

## Overview

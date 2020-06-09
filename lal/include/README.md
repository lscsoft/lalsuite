# The `include/lal/` path

This directory tree is here to enable in-situ building and testing
of LALSuite components that reference their own headers using the
'absolute' path, e.g:

```c
#include <lal/LALHeader.h>
```

No real files should be placed in this directory,
everything should go under `lib/`.
The headers from `lib/` or `swig/` are then symbolically linked
into this directory before compilation and the subpackage-level `include/`
directory added to the pre-processor include search path.

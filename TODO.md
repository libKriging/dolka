This document is intended to be read from the generated html output
`TODO.html`. However wherever this output is not suitably rendered you
can read the mardown version `TODO.md` which seems correct up to some
problems for math.

1 General
=========

-   The one-point optimization criteria such as and their derivatives
    could optionnaly be computed by the method or by a new method say .

2 Existing functions
====================

-   `contours` Allow to pass further arguments to `predict` or to the
    functions `other` and `otherGrad`. For instance we may want to
    change `type` in `predict` or `proxy` in `EI` used as `other`.

-   **\[x\]** `genoud_cache` Check that all arguments of `fn` are well
    used. Of course, optimization criteria such as `EI` needs `model`.

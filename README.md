# "inner" elliptic curves for use with Hyrax

This repository contains parameters for elliptic curves whose operations are
defined over fields with characteristics matching the group order of the M159,
M191, M221, and M255 (a.k.a., Curve25519) elliptic curves. This is useful for
applications where one wishes to do elliptic curve operations in the verifying
circuit for a ZK statement.

These curves are generated following the guidelines of Aranha et al.,
"[A note on high-security general-purpose elliptic curves](https://eprint.iacr.org/2013/647)"
(except, of course, for the field characteristic, which is determined by the
"outer" elliptic curve's group order).

## generator scripts

The `bin/` directory contains scripts for generating elliptic curve parameters.
These are sage scripts; on Debian-like systems, you should be able to do
something like

    apt-get install sagemath

and then run the scripts with (say) `sage 647.sage`.

- `647.sage` is the original script written by Samuel Neves to validate the
  curve parameters described by Aranha et al. in
  "[A note on high-security general-purpose elliptic curves](https://eprint.iacr.org/2013/647)".

- `find_curveparams.sage` generates parameters for one or more of the above
  curves (or you can easily modify it to generate parameters for your curve).
  It takes two arguments:
  
  `sage find_curveparams.sage [<curveNum> [<tryOffset>]]`

  where `<curveNum>` is in 0..3 and selects one of the four curves, and
  `<tryOffset>` defines where to start the parameter search. This is useful
  when running many threads in parallel: start one thread with tryOffset 0,
  another with 1, another with 2, etc.

## license

`bin/647.sage` was created and is owned by Samuel Neves.

All other files in this repository are released under the following license:

Copyright 2017 Riad S. Wahby <rsw@cs.stanford.edu> and the Hyrax authors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

# Rational Complex Integer Tiles

This repository provides an implementation of polygonal tiles based on complex integer rings,
based on the following blog posts:

* https://pirogov.de/blog/perfect-precision-2d-geometry-complex-integers/ 
* TODO: write followup post

Note that this is a work-in-progress, or actually, not-so-much-progress (due to lack of time).

## Roadmap

- [x] implement complex integer rings in a generic way
- [ ] implement efficient structure to represent polylines (chains of segments) using complex integers
- [ ] implement structure representing a simple polygon (i.e. a closed boundary described by a polyline)
- [ ] implement useful operations like enumerating polygons, checking possible matches, etc.

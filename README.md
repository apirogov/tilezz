# Rational Complex Integer Tiles

This repository provides two main features:

1. a practical implementation of various complex integer rings
2. grid-free but exact representation of a rich family of polygonal tiles

The target audience are primarily math enthusiasts (like myself) or researchers
working on periodic or aperiodic tiles or other areas where exact representation
of polygons is crucial and where numerical approximation and the use of floating
point arithmetic are not acceptable for the given purpose, yet the use of more
generic solutions, such as computer algebra systems is not desirable, too
inconvenient or too inefficient.

The concepts this work is based on are described in the following blog posts:

* https://pirogov.de/blog/perfect-precision-2d-geometry-complex-integers/
* TODO: write followup post

Note that due to time constraints, this is a work-(not-so-fast-)in-progress.

## Roadmap

- [x] implement [complex integer rings](https://en.wikipedia.org/wiki/Cyclotomic_field) generically
- [x] implement [simple polygonal chains](https://en.wikipedia.org/wiki/Polygonal_chain) over complex integers
- [ ] implement [simple polygons](https://en.wikipedia.org/wiki/Simple_polygon) with useful operations
- [ ] implement (interactive) visualization utilities (i.e. render images and/or use WebGL)
- [ ] write a tutorial showing how to e.g. describe the [Spectre tile](https://en.wikipedia.org/wiki/Einstein_problem) and quickly compute all its self-matches

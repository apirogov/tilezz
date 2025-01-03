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

## See Also

### Complex Integers

I have no clue of abstract algebra and number theory, but it seems like there
are a few different implementations of cyclotomic fields. I have have not tried
to compare the approaches or to benchmark anything yet, because the
implementations of the complex integers were not the intended main feature of
this crate.

* https://github.com/CyclotomicFields/cyclotomic (Rust)
* https://github.com/walck/cyclotomic (Haskell)

### Tiles and Tilings

It seems that people who work on/with tiles use software like:

* computer algebra systems, such as [SageMath](https://www.sagemath.org/)
* SAT or SMT solvers, such as [Z3](https://github.com/Z3Prover/z3)
* [PolyForm Puzzle Solver](https://www.jaapsch.net/puzzles/polysolver.htm)

SageMath or something similar is of course much more heavy and requires deeper
algebraic understanding to even ask what you want. This crate provides much
simpler (and I hope efficient) solutions to much more specific problems.

SAT or SMT solvers are excellent tools and I have some experience using them,
but encoding tiling problems into suitable formulas is far from trivial and
requires a lot of thought at work. It can pay off and work well if you have
enough knowledge about SAT/SMT encoding tricks and also have a deeper grasp on
the structure of the problem at hand. But it is not something I'd do to just
quickly come up with a shape and check out its behavior as a tile.

The PolyForm Solver seems to have some adoption by the tiling community
and looks like a mature and feature-rich package that I eventually want to try
out myself. From a cursory look, it seems like a major conceptual difference is that
the PolyForm Solver is grid-based - it requires selecting some underlying grid
of tiling pieces from which the polyform tiles can be built from. My approach is
completely grid-free, so it is more general. However, it also means less fixed
structure to work with, and I don't provide any sophisticated solvers (yet).

### Articles

* page about [polyform tilings](https://www.polyomino.org.uk/mathematics/polyform-tiling/)
* [Drawing the Aperiodic Hat Tile with Python and Z3](https://www.hgreer.com/HatTile/)
* [It's a Shape Jim, But Not as We Know It](https://hedraweb.wordpress.com/2023/03/23/its-a-shape-jim-but-not-as-we-know-it/)

## Roadmap

- [x] implement [complex integer rings](https://en.wikipedia.org/wiki/Cyclotomic_field) generically
- [x] implement [simple polygonal chains](https://en.wikipedia.org/wiki/Polygonal_chain) over complex integers
- [ ] implement [simple polygons](https://en.wikipedia.org/wiki/Simple_polygon) with useful operations
- [ ] implement (interactive) visualization utilities (i.e. render images and/or use WebGL)
- [ ] write a tutorial showing how to e.g. describe the [Spectre tile](https://en.wikipedia.org/wiki/Einstein_problem) and quickly compute all its self-matches

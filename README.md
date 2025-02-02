# tilezz

This repository provides the following main features:

1. practical implementation of various [cyclotomic rings and fields](https://en.wikipedia.org/wiki/Cyclotomic_field), which are subsets of complex numbers that admit exact representation
2. exact representation of a rich family of polygonal tiles based on these rings, using a minimalistic and discrete form of [turtle semantics](https://en.wikipedia.org/wiki/Turtle_graphics)
3. some useful string-based algorithms to work with these objects and visualize them

These features taken together grant us multiple freedoms:

1. freedom of floating-point numbers and subtle bugs due to numerical issues
2. freedom of coordinates, as polygons are described by relative angles and movements to trace out their boundary
3. freedom of a (typical) grid, when using the polygons in the context of tiles and tilings

The target audience are primarily math enthusiasts (like myself) or researchers
working on (a)periodic tiles or other areas where exact representation of
polygons is crucial and where numerical approximation and the use of floating
point arithmetic are not acceptable for the given purpose, yet the use of more
generic solutions, such as computer algebra systems is not desirable, too
inconvenient or too inefficient.

The concepts this work is based on are described in the following blog posts:

* https://pirogov.de/blog/perfect-precision-2d-geometry-complex-integers/
* https://pirogov.de/blog/intersecting-segments-without-tears/
* (more posts will probably added over time)

Note that due to time constraints, this is a work-(not-so-fast-)in-progress.

## Demonstration

### Exploring the cyclotomic ring ZZ12

<img src="https://github.com/user-attachments/assets/a7d1d698-8e7c-41a8-b2a8-49a8fed80c2e" width="45%" />
<img src="https://github.com/user-attachments/assets/d246fd60-bfed-4ff8-9393-2ae14e77e4d2" width="45%" />

Left: BFS from the origin, Right: BFS in the unit square, starting in the corners and normalizing discovered point modulo unit square. Each image shows the new points discovered in the corresponding round.

## Usage

### Essential Concepts

This crate provides the abstract geometric API for using concrete
representations of constructible cyclotomic rings for some fixed root of unity,
and polygonal tiles build on top of these rings.

Instead of representing polygons by segments and coordinates, they are described
using a form of turtle semantics, i.e. interpreting a sequence of discrete
[external angles](https://en.wikipedia.org/wiki/Internal_and_external_angles) as
instructions for tracing out a polygon or segment chain from a given starting
point by performing unit-length steps in some direction. This, together with the
fact that we use cyclotomics for actual coordinates, helps avoiding dependence
on floating point numbers and explicit coordinates.

Here is a conceptual mapping for the relevant geometrical objects:

* a **point** corresponds to a **turtle**, which is an "oriented point" (it has an angle, defining its facing direction)
* a **polygonal chain** corresponds to a **snake**, which consists of instructions for a turtle
* a **polygon** corresponds to a *closed* snake, which I call a **rat** (for *rational tile*)
* a **tile patch** corresponds to a **pack**, which is a collection of combined rats

A **tile** is a [simple polygon](https://en.wikipedia.org/wiki/Simple_polygon),
i.e. a polygon without holes or self-intersections, and a segment chain, polygon, or tile
is *rational* if all the side lengths can be expressed as integer multiples of a
common length.

By using cyclotomic integers for coordinates and expressing all geometric
objects in terms of unit steps into some direction, each simple polygon allows
for a natural representation as a sequence of exterior angles along its boundary.
As the sequence is cyclic, there is one cyclicaly shifted sequence for each
starting vertex.

The **canonical representation** is then simply the [lexicographically
minimal](https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation)
such sequence. Note that this gives us a **simple and efficient equivalence
check on polygons**: two polygons are equal iff they have the same canonical
angle sequence. Treating (rational) polygons as strings of angles also allows us
to use other efficient string-based algorithms, e.g. to compute combinations of
tiles.

This library also provides some plotting functionality based on
[plotterrs](https://github.com/plotters-rs/plotters), which means that you can
easily render tiles built with this crate into various backends, including PNG,
SVG, web pages (HTML5/WASM), including interactive usage in Jupyter.

### Interactive (Jupyter Notebook)

Thanks to the capabilities of the `plotters` library, it is easy to use
[Jupyter notebooks](https://jupyter.org/) with this crate to visualize polygonal tiles.

#### Requirements

* [Jupyter Notebook](https://jupyter.org/install#jupyter-notebook) (Arch Linux users see [here](https://wiki.archlinux.org/title/Jupyter))
* [evcxr](https://github.com/evcxr/evcxr)
* [evcxr_jupyter](https://github.com/evcxr/evcxr)

If you have all the dependencies installed correctly and can create/open Jupyter
Notebooks using a Rust kernel (check out the official
[evcxr documentation](https://plotters-rs.github.io/plotters-doc-data/evcxr-jupyter-integration.html)),
then you are already set up for using this crate interactively.

#### Rendering the Spectre Tile Over the Cyclotomic Integer Ring ZZ12

Here is how you can quickly construct and render the
[spectre tile](https://en.wikipedia.org/wiki/Einstein_problem):

**Step 1:** **(Recommended)** *This step is only needed if you want to use the most recent development version from this repository.*

Clone this repository and change to its directory, i.e. in a terminal run:

```bash
git clone https://github.com/apirogov/tilezz
cd tilezz
```

**Step 2:** Run `jupyter notebook`

**Step 3:** Open the [minimal example notebook](./examples/minimal.ipynb) **OR**
    Create a new Rust notebook (which is powered by `evcxr`) and add the following code into a cell:

```rust
// Build and load the crate
:dep plotters = { version = "^0.3.7", features = ["evcxr", "all_series"] }
// 1(a) if you want to use a published version of this crate:
// :dep tilezz = "*"
// 1(b) (RECOMMENDED) if you cloned this repository and want to use the current development version:
:dep tilezz = { path = ".." }

// Import what we need
use tilezz::cyclotomic::*;
use tilezz::snake::Turtle;
use tilezz::rat::Rat;

use plotters::prelude::*;
use tilezz::plotters::{plot_tile, TileStyle};

evcxr_figure((500,500), |root| {
    // Prepare the canvas
    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);

    // Define a sequence of external angles. As all segments have unit length, this fully determines a polygon
    let external_angles: &[i8] = &[3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
    // Instantiate an abstract polygon over the cyclotomic ring ZZ12 (one full turn = 12 rotational unit steps)
    let s: Rat<ZZ12> = external_angles.try_into().unwrap();
    // Trace out the polygon in the cartesian plane using a canonical starting point and facing direction
    let p = s.to_polyline_f64(Turtle::default());
    // Plot the concretized polygon with default settings
    plot_tile(&root, &p, &TileStyle::default());

    Ok(())
})
```

After waiting for some seconds (the required dependencies have to be built
first, after that it is faster), you should see a plot showing the spectre tile.

## See Also

### Tiles and Tilings

It seems that people who work on/with tiles use software like:

* computer algebra systems, such as [SageMath](https://doc.sagemath.org/html/en/reference/number_fields/index.html)
* SAT or SMT solvers, such as [Z3](https://github.com/Z3Prover/z3)
* [PolyForm Puzzle Solver](https://www.jaapsch.net/puzzles/polysolver.htm)

SageMath or similar systems are of course much more heavy and require deeper
algebraic understanding to even ask what you want. This crate provides much
simpler (and I hope more efficient) solutions to much more specific problems.

Similarly, SAT or SMT solvers are excellent tools for certain NP-complete
problems (I have some experience using them), but encoding tiling problems into
suitable formulas is far from trivial and requires some thought and work (even
though it sometimes is [possible](https://www.hgreer.com/HatTile/)). It can
really pay off and work well if you have enough knowledge about SAT/SMT encoding
tricks and also have a deeper grasp on the structure of the problem at hand. But
it is not something I'd pull out to just quickly come up with a polygon and
check its behavior as a tile.

The PolyForm Solver seems to have
[some adoption](https://hedraweb.wordpress.com/2023/03/23/its-a-shape-jim-but-not-as-we-know-it/)
by the tiling community and looks like a mature and feature-rich package that I
eventually want to try out myself. From a cursory look, it seems like the
biggest conceptual difference is that **the PolyForm Solver is grid-based** - it
requires selecting some underlying grid of primitive cells from which the
polyform tiles can be built from. **My approach is grid-free**, so it is more
general. Technically, the cyclotomic integers *do* provide another kind of grid,
but it is not a typical periodic cell grid as typically used to describe tiles
and tilings.  This means that more exotic and irregular tiles can be expressed,
but also that there is less fixed structure to work with and exploit, and I
don't provide any sophisticated solvers (yet).

### Cyclotomic Integer Rings (i.e. Complex Integers)

I have no clue about abstract algebra and number theory (I just stumbled into
this topic trying to represent some tiles exactly), but it seems like there are
a few very general implementations of cyclotomic fields.

* https://github.com/CyclotomicFields/cyclotomic (Rust)
* https://github.com/walck/cyclotomic (Haskell)

Compared to these packages, I do not try to implement the full set of cyclotomic
integers (i.e. one type that includes and works with all roots of unity at
once). This crate provides **a *separate* ring/field for each supported root of
unity** instead, which then serves as the backbone for representing geometry of
suitable polygons.

The corresponding underlying representations are optimized for and limited to
the respective ring, there is **no overhead due to management of symbolic
representations** or anything like that, because for each ring, the provided
data type encodes the values of each ring as vectors over a linearly independent
set base of units (and all their distinct symbolic products), together with a
ring-specific implementation of multiplication, which hard-codes the symbolic
simplifications of expressions that appear during the evaluation of
multiplication.

I have have not tried to compare this crate to the other approaches or benchmark
anything yet, because the implementations of the complex integers were not the
intended main feature of this crate. If someone is mainly interested in this
crate for the implementation of the cyclotomic integer rings, I would be happy
to get some feedback on how this compares to the more generic implementations
with similar features.

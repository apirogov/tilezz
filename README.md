# tilezz

This repository provides two main features:

1. practical implementation of various cyclotomic integer rings (I simply call them complex integers)
2. grid-free (but still exact) representation of a rich family of polygonal tiles based on these rings

The target audience are primarily math enthusiasts (like myself) or researchers
working on (a)periodic tiles or other areas where exact representation of
polygons is crucial and where numerical approximation and the use of floating
point arithmetic are not acceptable for the given purpose, yet the use of more
generic solutions, such as computer algebra systems is not desirable, too
inconvenient or too inefficient.

The concepts this work is based on are described in the following blog posts:

* https://pirogov.de/blog/perfect-precision-2d-geometry-complex-integers/
* TODO: write followup post

Note that due to time constraints, this is a work-(not-so-fast-)in-progress.


## Usage

This crate provides the abstract geometric API for using concrete
representations of complex integer rings (i.e. cyclotomic fields for some fixed
root of unity without division) and polygonal tiles build on top of these rings.

It also provides some plotting functionality based on
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

**Step 3:** Open the [example notebook](./tilezz_example.ipynb) **OR**
    Create a new Rust notebook (which is powered by `evcxr`) and add the following code into a cell:

```rust
:dep plotters = { version = "^0.3.7", features = ["evcxr", "all_series"] }
// 1(a) if you want to use a published version of this crate:
// :dep tilezz = "*"
// 1(b) (RECOMMENDED) if you cloned this repository and want to use the current development version:
:dep tilezz = { path = "." }

use plotters::prelude::*;
use tilezz::snake::constants::spectre;
use tilezz::snake::{Snake, Turtle};
use tilezz::plotters::plot_tile;
use tilezz::zz::ZZ12;

evcxr_figure((500,500), |root| {
    let s: Snake<ZZ12> = spectre();
    let tile = s.to_polyline_f64(&Turtle::default());

    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);

    plot_tile(
        &mut ChartBuilder::on(&root)
            .caption("Spectre tile", ("sans-serif", 40).into_font())
            .x_label_area_size(20)
            .y_label_area_size(40),
        &tile,
    );
    root.present()?;
    Ok(())
})
```

After waiting for some seconds (the required dependencies have to be built
first, after that it is faster), you should see a plot showing the spectre tile.

TODO: improve/finalize API, update the example and add a screenshot.

## See Also

### Tiles and Tilings

It seems that people who work on/with tiles use software like:

* computer algebra systems, such as [SageMath](https://www.sagemath.org/)
* SAT or SMT solvers, such as [Z3](https://github.com/Z3Prover/z3)
* [PolyForm Puzzle Solver](https://www.jaapsch.net/puzzles/polysolver.htm)

SageMath or something similar is of course much more heavy and requires deeper
algebraic understanding to even ask what you want. This crate provides much
simpler (and I hope more efficient) solutions to much more specific problems.

SAT or SMT solvers are excellent tools for certain NP-complete problems (I have
some experience using them), but encoding tiling problems into suitable formulas
is far from trivial and requires some thought and work (even though it sometimes
is [possible](https://www.hgreer.com/HatTile/)). It can really pay off and work
well if you have enough knowledge about SAT/SMT encoding tricks and also have a
deeper grasp on the structure of the problem at hand. But it is not something
I'd pull out to just quickly come up with a polygon and check its behavior as a
tile.

The PolyForm Solver seems to have
[some adoption](https://hedraweb.wordpress.com/2023/03/23/its-a-shape-jim-but-not-as-we-know-it/)
by the tiling community and looks like a mature and feature-rich package that I
eventually want to try out myself. From a cursory look, it seems like the
biggest conceptual difference is that **the PolyForm Solver is grid-based** - it
requires selecting some underlying grid of tiling pieces from which the polyform
tiles can be built from. **My approach is completely grid-free**, so it is more
general.  However, it also means less fixed structure to work with and exploit,
and I don't provide any sophisticated solvers (yet).

### Cyclotomic Integer Rings (i.e. Complex Integers)

I have no clue about abstract algebra and number theory (I just stumbled into
this topic trying to represent some tiles exactly), but it seems like there are
a few very general implementations of cyclotomic fields.

* https://github.com/CyclotomicFields/cyclotomic (Rust)
* https://github.com/walck/cyclotomic (Haskell)

Compared to these packages, I do not try to implement the full set of cyclotomic
integers (i.e. one type that includes and works with all roots of unity at
once). This crate provides  **a *separate* field for each supported root of
unity** instead, which then serves as the backbone for representing geometry of
suitable polygons.

The corresponding underlying representations are optimized for the respective
ring, there is **no overhead due to management of symbolic representations** or
anything like that, because for each ring, the provided data type encodes the
values of each ring as vectors over a linearly independent set base of units
(and all their distinct symbolic products), together with a ring-specific
implementation of multiplication, which hard-codes the symbolic simplifications
of expressions that appear during the evaluation of multiplication.

I have have not tried to compare this crate to the other approaches or benchmark
anything yet, because the implementations of the complex integers were not the
intended main feature of this crate. If someone is mainly interested in this
crate for the implementation of the cyclotomic integer rings, I would be happy
to get some feedback on how this compares to the more generic implementations
with similar features.

## Roadmap

- [x] implement [complex integer rings](https://en.wikipedia.org/wiki/Cyclotomic_field) generically
- [x] implement [simple polygonal chains](https://en.wikipedia.org/wiki/Polygonal_chain) over complex integers
- [ ] implement [simple polygons](https://en.wikipedia.org/wiki/Simple_polygon) with useful operations
- [ ] implement (interactive) visualization utilities (i.e. render images and/or use WebGL)

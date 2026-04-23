# Plan: Vertex Type Enumeration 

# Preparation

A tileset consists of multiple non-empty unique tiles.

Each tile has a canonical sequence and id.

For this we have the Rat and the TileSet data types.

TileSet should ensure non-emptyness and uniqueness by construction.

# Patch

A patch is simply a combination of multiple tiles, so a tile is atomic and a
patch decomposes into tiles that are touching along at least one edge of their
boundary.

We only consider hole-free patches as valid.

# Matching Types

Each two adjacent tiles in a patch share a common boundary along which they
were glued.

For a given tileset, we compute all pairwise maximal matches of all tiles in
the tileset. This we can do efficiently by using the cyclic match index and
filtering for valid matches.

When initializing a tileset, we compute all possible distinct valid matches and
give them canonical ids by sorting and enumerating them.  These ids we call the
type of a matching type, each matching type identifies a different valid match
(tileId1, index1, index2, tileId2).

Note the symmetry of unique matches (when not picking a "first" and "second"
tile). We canonicalize matches by keeping only the lexicographically smaller of
(t1, i1, i2, t2) or (t2, i2, i1, t1), which are describing the same match from
two sides.

Furthermore note that for tiles that have boundaries w which are u^i, i.e. have
i-fold rotational symmetry, we have to consider only matches modulo offsets of
|w|/i, because the result would be the same, i.e. for an unmarked tile we
cannot distinguish these matches and canonically keep only those along the
first repetition of the repeated factor.

That classifies all distinct ways how two adjacent tiles (i.e. those touching
in at least one edge) connect to each other.

We start match indexing at 1, so that in settings where order matters (e.g. we
are given tiles x1 and x2 and match m) we can use the sign to say by (x1, m,
x2) that the match type m is applied "forward" if m > 0 (x1 is of type t1, x2
of type t2) or "backward" otherwise.

That means, for each patch there is a connected graph where tiles are nodes
labelled by tile ids and edges are shared boundaries labelled by the respective
signed matching types.

# Vertex Types

A patch vertex is a point where multiple tile boundaries meet.

A valid patch vertex can have AT MOST one side from which a tile can be added,
i.e. a boundary pattern where a vertex is touched from both sides implies a
hole in the patch, and we consider patches with holes as invalid.

This is taken care of by the matching logic that is already implemented via
geometric Snake validation.

A patch vertex can be non-terminal (some tile can be added) 
or terminal (no tile can be added).

A patch vertex can also be 
* closed (fully surrounded by tiles, implies terminal), or
* open (it is not fully surrounded, i.e. the angles do not sum to a full circle)

A patch vertex type is uniquely identified by a canonical sequence of the
underlying signed tile matchings in counter-clockwise order. 
For open vertices the listing order is naturally given by the remaining gap,
for closed types we take the lexicographically minimal cyclic rotation of
the matching types.

Each single matching type induces two vertex types, which correspond to
focusing one of the endpoints or equivalently, the direction that the matching
is to be read (i.e. the two possible signs determining which tile is "first").

The set of all signed matchings as singleton lists is exactly our set of
initial vertex types. These are necessarily open, because it is impossible to
close a vertex with only two tiles. Note that a vertex type can be both initial
and terminal, by definition.

To compute all other vertex types, we start from these initial types.

while there are still vertex types that are open and not marked as terminal:
    take such a vertex type
    Compute/retrieve the Rat/sequence it corresponds to
    compute all ways to add a tile from the tileset that touches the clockwise edge from the vertex position 
        (note that the resulting patch the sequence has some index corresponding to the considered vertex type, we only care about matches involving that index).
    If no valid match exists, mark the type as terminal.
    If valid result matches exist, then:
        for each match:
            determine if the vertex is closed (= both the edge before AND after the vertex is matched, open means only the edge before the vertex is matched)
        add the successful match to the vertex type sequence 
        if there is no such type yet in our vertex type list yet:
            create a new vertex type
            if the new vertex type is open:
                add it to the vertex types to process

It is easy to see how we can augment the algorithm to construct a DAG of vertex
types with layers going from initial to terminal vertex types, connected by
edges labeled by signed matches.

For each vertex type we store the resulting rat and an index where the vertex
is anchored, if it is open (or nothing if it is closed, i.e. inner vertex
type).

## Constructible Vertex Words

Let the largest tile in the tileset have boundary length k.

Then the action of adding a single tile to an existing patch can affect at most k edges.

Thus to know how the boundary of any patch can evolve under addition of tiles,
we only need to understand all boundary interactions up to length k.

W.l.o.g. assume a non-trivial patch (at least two tiles).

Each patch boundary can be decomposed into a sequence of open vertex types.

Note that these boundary vertex types exactly determine the boundary tiles and
how they are matched together, but they do not uniquely determine the inside of
a patch (i.e. fully surrounded tiles) because there can be multiple distinct
patches with the same boundary, thus not necessarily a unique inner patch. 

Between every two boundary vertices (v1, v2) there is a part of a boundary of a
single tile. Note that the set of tile matching types does not necessarily
contain all such boundaries, because such a boundary is determined is induced
by two matchings, not one, i.e. one from each vertex.

Not all sequences of vertex types are possible. To collect all realizable
vertex words, i.e. those which can appear as subword along a patch boundary, we
can proceed as follows.

Let T be the primitive tile set, U = T be the set of unprocessed representative patches and W the set of collected vertex words (with a representative patch realizing it).

Note that alongside the boundary we always store respective indices where the vertex types are anchored.

While U not empty:
    take all patches out of U
    compute all candidate matches between these patches from U and all tiles in T
        (using cyclic match index)
    throw away invalid sequences
        (cyclotomic exact geometry / snake based validation)

    for each resulting match:
        compute and annotate the induced vertex types induced by the new tile
        collect all vertex words that bound an edge sequence of length at most k along the resulting patch

    add all newly discovered vertex words to W
    U = all result patches that contributed new vertex words to W

# Boundary Vertex Word Rewriting Grammar

We can use the resulting words and witnesses to systematically derive a cyclic string rewriting grammar:

For each collected word, compute all ways how a tile can attach to the described boundary. so we compute all valid matches between the set of patches in W and the set of tiles T.

We extract the rewriting rules by computing the difference in the vertex
sequence before and after tile addition, i.e. 

if a patch with vertex sequence uvwabcxyz matched with a tile resulting in uvwdexyz
then we can track the rewriting rule as: abc -> de
If a tile was attached between vertices x and y without modifying them 
and induced new vertices a and b,
we still write it as a rule of the form: xy -> xaby

Note that each rule eliminates 0 or more vertices (the respective vertices become closed vertices on the inside of the patch) and produces exactly two new vertices (induced by the match of the new tile to the patch).

To complete the grammar we add

S -> Ti   for all tiles Ti
Ti -> vw  for all initial vertex words from 2 tiles where at least one tile is Ti

This grammar describes boundaries of all valid (i.e. hole-free
non-self-intersecting) patches that can be built from the given tileset.

Derived words that contain open terminal vertex types correspond to patches that can not be extended to a tiling of the plane, as by construction we know that no tile can attach to them in a way that they become closed.

//! General plotting utilities

/// 2D Point (X, Y)
pub type P64 = (f64, f64);

/// 2D Rectangle (given as pair of corners with respectively minimal and maximal coordinate values)
pub type R64 = (P64, P64);

/// Returns the center of the tile based on its vertex points.
pub fn tile_centroid<'a, I>(pts: I) -> P64
where
    I: IntoIterator<Item = &'a P64>,
{
    let mut n = 0_usize;
    let (mut cx, mut cy) = (0.0, 0.0);
    for (x, y) in pts {
        cx += x;
        cy += y;
        n += 1;
    }
    (cx / n as f64, cy / n as f64)
}

/// Tight math-coord bounds covering every point yielded by `layers`.
///
/// Each "layer" is its own iterator of points; layers can be empty
/// (they're simply skipped — unlike [`tile_bounds`], an empty layer
/// does not widen the result to a unit square). Returns `None` if
/// every layer is empty; callers pick the fallback that fits their
/// context.
pub fn points_bounds<'a, I, J>(layers: I) -> Option<R64>
where
    I: IntoIterator<Item = J>,
    J: IntoIterator<Item = &'a P64>,
{
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut any = false;
    for layer in layers {
        for &(x, y) in layer {
            any = true;
            min_x = min_x.min(x);
            min_y = min_y.min(y);
            max_x = max_x.max(x);
            max_y = max_y.max(y);
        }
    }
    if any {
        Some(((min_x, min_y), (max_x, max_y)))
    } else {
        None
    }
}

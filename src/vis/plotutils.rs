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

/// Given a list of (x,y) coordinates, return the bounds
/// ((min_x, min_y), (max_x, max_y)).
pub fn tile_bounds<'a, I>(pts: I) -> R64
where
    I: IntoIterator<Item = &'a P64>,
{
    let mut empty = true;
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for (x, y) in pts {
        empty = false;
        (min_x, min_y) = (min_x.min(*x), min_y.min(*y));
        (max_x, max_y) = (max_x.max(*x), max_y.max(*y));
    }
    if !empty {
        ((min_x, min_y), (max_x, max_y))
    } else {
        ((-1., -1.), (1., 1.))
    }
}

/// Like tile_bounds, but takes a sequence of sequences of points and computes overall bounds.
pub fn tiles_bounds<'a, I, J>(ptss: I) -> R64
where
    I: IntoIterator<Item = J>,
    J: IntoIterator<Item = &'a P64>,
{
    let pts: Vec<P64> = ptss
        .into_iter()
        .map(|p| <(P64, P64) as Into<[P64; 2]>>::into(tile_bounds(p.into_iter())).into_iter())
        .flatten()
        .collect::<Vec<P64>>();
    return tile_bounds(&pts);
}

/// Given the bounds of a tile, adjust them to a square
/// centered on and including all the points.
pub fn tile_viewport(((min_x, min_y), (max_x, max_y)): R64) -> R64 {
    let (w, h) = (max_x - min_x, max_y - min_y);
    let half_d = 0.5 * if w >= h { w - h } else { h - w };
    let pad = 1_f64;
    let pad_x = pad + if w < h { half_d } else { 0_f64 };
    let pad_y = pad + if h < w { half_d } else { 0_f64 };

    (
        ((min_x - pad_x), (min_y - pad_y)),
        ((max_x + pad_x), (max_y + pad_y)),
    )
}

/// Given the size of a box where a chart has to fit and bounds in the
/// chart-internal coordinate system, returns required padding to render the
/// coordinate system without distortion.
/// For centering, half of the needed padding should be applied on each respective side.
pub fn chart_padding(
    (box_w, box_h): (u32, u32),
    ((min_x, min_y), (max_x, max_y)): R64,
) -> (u32, u32) {
    let (w, h) = (max_x - min_x, max_y - min_y);

    let box_ratio: f64 = box_h as f64 / box_w as f64;
    let bounds_ratio: f64 = h as f64 / w as f64;

    return if box_ratio > bounds_ratio {
        // embedded chart is wider
        let pad_y = box_h as f64 - box_w as f64 * bounds_ratio;
        (0, pad_y as u32)
    } else {
        // embedded chart is taller
        let pad_x = box_w as f64 - box_h as f64 / bounds_ratio;
        (pad_x as u32, 0)
    };
}

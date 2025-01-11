// 2D Point
type P64 = (f64, f64);

// 2D Rectangle (top left and bottom right corners)
type R64 = (P64, P64);

/// Returns the center of the tile based on its vertex points.
pub fn tile_centroid<'a, I>(pts: I) -> P64
where
    I: IntoIterator<Item = &'a (f64, f64)>,
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
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for (x, y) in pts {
        (min_x, min_y) = (min_x.min(*x), min_y.min(*y));
        (max_x, max_y) = (max_x.max(*x), max_y.max(*y));
    }
    ((min_x, min_y), (max_x, max_y))
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

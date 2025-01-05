/// Given a list of (x,y) coordinates, return the bounds
/// ((min_x, min_y), (max_x, max_y)).
pub fn tile_bounds(pts: &[(f64, f64)]) -> ((f64, f64), (f64, f64)) {
    let (x, y) = pts[0];
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (x, y, x, y);
    for (x, y) in pts {
        (min_x, min_y) = (min_x.min(*x), min_y.min(*y));
        (max_x, max_y) = (max_x.max(*x), max_y.max(*y));
    }
    ((min_x, min_y), (max_x, max_y))
}

/// Given a list of (x,y) coordinates, returns the bounds
/// ((min_x, min_y), (max_x, max_y)) of a square
/// centered on and including all the points.
pub fn tile_viewport(pts: &[(f64, f64)]) -> ((f64, f64), (f64, f64)) {
    let ((min_x, min_y), (max_x, max_y)) = tile_bounds(pts);
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

/// Returns the center of the tile based on its bounds.
pub fn tile_center(((min_x, min_y), (max_x, max_y)): &((f64, f64), (f64, f64))) -> (f64, f64) {
    (min_x + (max_x - min_x) / 2., min_y + (max_y - min_y) / 2.)
}

/// Returns the center of the tile based on its vertex points.
pub fn tile_mass_center(pts: &[(f64, f64)]) -> (f64, f64) {
    let n = pts.len() as f64;
    let (mut cx, mut cy) = (0.0, 0.0);
    for (x, y) in pts {
        cx += x;
        cy += y;
    }
    (cx / n, cy / n)
}

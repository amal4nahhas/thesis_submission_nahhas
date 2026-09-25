#!/usr/bin/env python3
"""
Convert a flat-color sketch (e.g. exported from Adobe Fresco) into particle positions
for the SPH simulation.

Convention:
  - Blue, solidly filled regions -> fluid.
  - Dark/black regions           -> boundary walls.

Both are filled the same way: a candidate grid at `--spacing` across the whole domain,
keeping points whose pixel falls inside the corresponding color mask. This handles any
shape uniformly - thick strokes, closed borders, branching outlines, freeform blobs -
with no assumption about the stroke being a single simple non-branching path.

Usage:
    .venv/bin/python sketch_to_scene.py sketch.png \
        --domain-width 220 --domain-height 120 --spacing 1.0 \
        --out scene.txt

Output format: one particle per line, "x y type", type in {fluid, boundary}.
"""
import argparse
import sys

import cv2
import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.spatial import cKDTree


def load_masks(image_path, smooth_px=0.0):
    img = Image.open(image_path).convert('RGB')
    arr = np.asarray(img).astype(np.int16)  # (H, W, 3)
    r, g, b = arr[..., 0], arr[..., 1], arr[..., 2]
    gray = (r + g + b) / 3.0

    # Fluid fill: any "blue family" color (blue-purple, cyan, ...) - blue channel clearly
    # dominant over red (this is what separates it from black/gray/white), without also
    # requiring blue to dominate green, since cyan has green almost as high as blue.
    fluid_mask = (b > 150) & (b - r > 60)

    # Dark stroke/region: low overall brightness, excluding anything already claimed by
    # fluid (keeps anti-aliased edges of the blue fill from leaking into the boundary mask).
    boundary_mask = (gray < 140) & ~fluid_mask

    if smooth_px > 0:
        # Blur each mask and re-threshold at 0.5: rounds off the pixel-staircase jaggies
        # a hand-drawn/rasterized curve has, without changing the shape's overall size.
        # A stroke thinner than ~2*smooth_px can lose most of its mass below the 0.5
        # threshold and nearly vanish - warn loudly rather than silently under-produce.
        fluid_before, boundary_before = fluid_mask.sum(), boundary_mask.sum()
        fluid_mask = gaussian_filter(fluid_mask.astype(float), sigma=smooth_px) > 0.5
        boundary_mask = gaussian_filter(boundary_mask.astype(float), sigma=smooth_px) > 0.5
        for name, before, after in (('fluid', fluid_before, fluid_mask.sum()),
                                    ('boundary', boundary_before, boundary_mask.sum())):
            if before > 0 and after < before * 0.5:
                print(f'WARNING: --smooth {smooth_px} shrank the {name} mask from {before} to '
                      f'{after} px ({after / before:.0%}) - the stroke is probably thinner than '
                      f'the blur radius. Try a smaller --smooth value for this image.')

    return fluid_mask, boundary_mask, img.size  # size = (width, height)


def pixel_to_world(px, py, img_w, img_h, domain_w, domain_h):
    u = px / img_w
    v = py / img_h
    return u * domain_w, v * domain_h


def world_to_pixel(wx, wy, img_w, img_h, domain_w, domain_h):
    # No flip: this simulation's world space is y-down (gravity is +y, "pointing
    # downwards" - see ParticleSimulation.cpp), same orientation as image pixel space.
    u = wx / domain_w
    v = wy / domain_h
    px = u * img_w
    py = v * img_h
    return px, py


def sample_mask(mask, px, py):
    x, y = int(round(px)), int(round(py))
    if 0 <= y < mask.shape[0] and 0 <= x < mask.shape[1]:
        return bool(mask[y, x])
    return False


def fill_region(mask, img_w, img_h, domain_w, domain_h, spacing):
    """Candidate grid at `spacing` across the whole domain, kept where it lands inside
    `mask`. Used for both fluid and boundary - makes no assumption about the region's
    shape (branching, closed, thick, thin all work the same way)."""
    points = []
    y = spacing * 0.5
    while y < domain_h:
        x = spacing * 0.5
        while x < domain_w:
            px, py = world_to_pixel(x, y, img_w, img_h, domain_w, domain_h)
            if sample_mask(mask, px, py):
                points.append((x, y))
            x += spacing
        y += spacing
    return points


def resample_polyline_world(world_points, spacing, closed):
    """Resample an ordered polyline (already in world coordinates) at even arc-length
    spacing, so curvature doesn't change local particle density - each particle is
    `spacing` arc-length away from the last, regardless of how sharply the curve bends."""
    pts = np.array(world_points, dtype=float)
    if closed:
        pts = np.vstack([pts, pts[0]])
    if len(pts) < 2:
        return [tuple(p) for p in pts]

    seg_vecs = np.diff(pts, axis=0)
    seg_lens = np.hypot(seg_vecs[:, 0], seg_vecs[:, 1])
    cum_len = np.concatenate([[0.0], np.cumsum(seg_lens)])
    total_len = cum_len[-1]
    if total_len == 0:
        return [tuple(pts[0])]

    out = []
    d = 0.0
    seg_idx = 0
    while d <= total_len:
        while seg_idx < len(seg_lens) - 1 and cum_len[seg_idx + 1] < d:
            seg_idx += 1
        seg_start_len = cum_len[seg_idx]
        seg_len = seg_lens[seg_idx] if seg_lens[seg_idx] > 0 else 1e-9
        t = (d - seg_start_len) / seg_len
        t = min(max(t, 0.0), 1.0)
        p = pts[seg_idx] + t * (pts[seg_idx + 1] - pts[seg_idx])
        out.append((p[0], p[1]))
        d += spacing
    return out


def trace_boundary_contours(boundary_mask, img_w, img_h, domain_w, domain_h, spacing,
                            min_contour_px=8):
    """Find the actual outline(s) of the dark region(s) - both outer edges and inner
    edges around any enclosed white gaps (e.g. both walls of a channel) - and resample
    each one at even arc-length spacing. Unlike sampling a grid over the filled area,
    this puts particles only on the curve itself: far fewer particles, and smooth,
    consistent spacing that follows the actual geometry instead of the pixel grid."""
    mask_u8 = (boundary_mask.astype(np.uint8)) * 255
    contours, _ = cv2.findContours(mask_u8, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)

    all_points = []
    for contour in contours:
        if cv2.arcLength(contour, closed=True) < min_contour_px:
            continue  # drop noise specks
        pixel_pts = contour.reshape(-1, 2)  # (N, 2) as (x, y)
        world_pts = [
            pixel_to_world(px, py, img_w, img_h, domain_w, domain_h)
            for px, py in pixel_pts
        ]
        all_points.extend(resample_polyline_world(world_pts, spacing, closed=True))

    return all_points


def clearance_mask(boundary_mask, img_w, img_h, domain_w, domain_h, min_clearance):
    """True wherever a point is at least `min_clearance` (world units) from the nearest
    boundary pixel in ANY direction. Unlike checking distance to the nearest boundary
    particle, this correctly flags the middle of a narrow corridor too - a point there can
    be far from any single boundary sample point along the wall, yet still sandwiched
    between two close walls. Built from a distance transform of the boundary mask, so it's
    a per-pixel property of the drawing's geometry, independent of --spacing."""
    # Pixels are generally non-square in world units (domain_w/img_w != domain_h/img_h),
    # so scale each axis separately before combining into a single Euclidean distance.
    dist_px = distance_transform_edt(~boundary_mask, sampling=(1, 1))
    dist_world = dist_px * (domain_w / img_w)  # approx: assumes near-uniform x/y scale
    return dist_world >= min_clearance


def enforce_min_spacing(points, min_dist, seed_points=None):
    """Greedily thin `points` so no two survivors - and no survivor and any point already in
    `seed_points` - end up closer than `min_dist`.

    Every boundary curve (a stroke's own edge, an added inward layer, a floor row, ...) gets
    traced and arc-length-resampled independently, with zero phase coordination against any
    other curve. Wherever two of them happen to run close together - a thin stroke whose two
    traced edges are only slightly more than a spacing apart, a second layer offset from the
    first, a floor row meeting a wall - nothing stops their resampled points from landing
    right on top of each other. In practice this produced boundary-boundary distances as
    small as 1/1000th of a particle spacing, which is exactly what crashed the simulation:
    near-zero separation makes the SPH kernel gradient blow up, sending that particle's
    position to NaN/huge values on the very next step, which then produces an out-of-range
    cell index in the neighbor search and a segfault. This is the single place that enforces
    the invariant every particle - fluid or boundary - is meant to satisfy: at least `spacing`
    from its nearest neighbor, regardless of which curve either point came from.

    Uses a uniform grid bucketed at min_dist so the average point only checks its 3x3
    neighborhood of cells rather than every previously kept point."""
    cell = min_dist
    grid = {}
    kept = []

    def cell_of(p):
        return (int(p[0] // cell), int(p[1] // cell))

    def add_to_grid(p):
        idx = len(kept)
        kept.append(p)
        grid.setdefault(cell_of(p), []).append(idx)

    for p in (seed_points or []):
        add_to_grid(p)

    survivors = []
    min_dist_sq = min_dist * min_dist
    for p in points:
        cx, cy = cell_of(p)
        too_close = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for qi in grid.get((cx + dx, cy + dy), ()):
                    qx, qy = kept[qi]
                    if (p[0] - qx) ** 2 + (p[1] - qy) ** 2 < min_dist_sq:
                        too_close = True
                        break
                if too_close:
                    break
            if too_close:
                break
        if not too_close:
            add_to_grid(p)
            survivors.append(p)
    return survivors


def drop_fluid_too_close_to_boundary(fluid_points, boundary_points, min_dist):
    """Fluid is grid-sampled on a regular lattice, but boundary is traced along an
    independent curve - nothing guarantees they don't land almost on top of each other
    where a wall passes through the fluid grid. At near-zero separation the SPH kernel
    gradient blows up, kicking that fluid particle out violently on frame 0. Drop any
    fluid particle nearer than `min_dist` to a boundary particle instead - a small gap
    right at the wall beats an unstable simulation."""
    if not boundary_points or not fluid_points:
        return fluid_points
    boundary_tree = cKDTree(np.array(boundary_points))
    fluid_arr = np.array(fluid_points)
    dist, _ = boundary_tree.query(fluid_arr)
    kept = fluid_arr[dist >= min_dist]
    return [tuple(p) for p in kept]


def add_floor(boundary_points, domain_w, spacing, layers, y_start):
    """Append `layers` full-width rows of boundary particles, closing off an open bottom
    (e.g. a sketch whose strokes just stop at the canvas edge instead of forming a floor)."""
    for i in range(layers):
        y = y_start + i * spacing
        x = spacing * 0.5
        while x < domain_w:
            boundary_points.append((x, y))
            x += spacing
    return boundary_points


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('image', help='Path to the sketch (PNG).')
    parser.add_argument('--domain-width', type=float, required=True, help='World-unit width of the scene.')
    parser.add_argument('--domain-height', type=float, required=True, help='World-unit height of the scene.')
    parser.add_argument('--spacing', type=float, required=True, help='Particle spacing, world units.')
    parser.add_argument('--smooth', type=float, default=0.0,
                        help='Gaussian-blur radius in source pixels applied to both masks before '
                             'sampling, to round off jagged hand-drawn/rasterized edges (try 1-2).')
    parser.add_argument('--floor-layers', type=int, default=0,
                        help='Append this many full-width boundary rows at the bottom of the '
                             'domain, closing off strokes that end open instead of forming a floor.')
    parser.add_argument('--min-separation', type=float, default=None,
                        help='Minimum allowed distance between a fluid particle and the nearest '
                             'boundary particle; closer fluid particles are dropped, since a near-'
                             'zero gap makes the SPH kernel gradient blow up on frame 0. Defaults '
                             'to 0.6x --spacing; pass 0 to disable.')
    parser.add_argument('--min-clearance', type=float, default=None,
                        help="Minimum open width (world units) a spot needs on ALL sides to get "
                             "fluid - omits fluid from corridors narrower than the kernel's support "
                             "width, where it would be squeezed by multiple nearby walls at once "
                             "regardless of --min-separation. Excluded spots are simply left open "
                             "(no fluid, no boundary) - draw sketches so nothing is fully enclosed "
                             "by boundary with no escape, rather than relying on this to paper over "
                             "it. Defaults to 2x --spacing (the kernel support diameter); pass 0 to "
                             "disable.")
    parser.add_argument('--boundary-dedup', type=float, default=None,
                        help="Minimum distance between two boundary particles before one is "
                             "dropped as a near-duplicate of the other. This is a safety net "
                             "against genuine near-coincident points (two independently-traced "
                             "curves - a stroke's inner/outer edge, a floor row, ... - happening "
                             "to resample almost on top of each other), NOT a requirement that "
                             "boundary particles be a full --spacing apart in general: modest "
                             "clustering (e.g. two edges of a thin stroke a few hundredths of a "
                             "unit apart) is normal and harmless for scenes that are only ever "
                             "rendered, not used for performance analysis. Defaults to "
                             "0.15x --spacing; pass 0 to disable entirely.")
    parser.add_argument('--out', required=True, help='Output positions file.')
    args = parser.parse_args()

    min_separation = args.spacing * 0.6 if args.min_separation is None else args.min_separation
    min_clearance = args.spacing * 2.0 if args.min_clearance is None else args.min_clearance
    boundary_dedup = args.spacing * 0.15 if args.boundary_dedup is None else args.boundary_dedup

    fluid_mask, boundary_mask, (img_w, img_h) = load_masks(args.image, args.smooth)
    print(f'Image: {img_w}x{img_h}px  ->  domain {args.domain_width}x{args.domain_height} world units')
    print(f'Fluid-colored pixels: {fluid_mask.sum()}   Boundary-colored pixels: {boundary_mask.sum()}')

    if min_clearance > 0:
        clear_mask = clearance_mask(boundary_mask, img_w, img_h, args.domain_width, args.domain_height, min_clearance)
        omitted = int((fluid_mask & ~clear_mask).sum())
        fluid_mask = fluid_mask & clear_mask
        print(f'Omitted {omitted} fluid-colored pixels narrower than {min_clearance:.4f} clearance '
              f'(closed/braided gaps) - left open, not turned into boundary.')

    fluid_points = fill_region(fluid_mask, img_w, img_h, args.domain_width, args.domain_height, args.spacing)
    boundary_points = trace_boundary_contours(boundary_mask, img_w, img_h, args.domain_width, args.domain_height, args.spacing)

    if args.floor_layers > 0:
        y_start = args.domain_height - args.floor_layers * args.spacing
        add_floor(boundary_points, args.domain_width, args.spacing, args.floor_layers, y_start)

    # Every curve above (a stroke's inner/outer edge, a floor row, ...) was resampled
    # independently with no coordination against the others, so a handful of points can land
    # almost exactly on top of one another - a light dedup pass removes only those genuine
    # near-duplicates (see --boundary-dedup above for why this isn't a full --spacing
    # requirement).
    if boundary_dedup > 0:
        before = len(boundary_points)
        boundary_points = enforce_min_spacing(boundary_points, boundary_dedup)
        print(f'Thinned {before - len(boundary_points)} boundary particles that landed closer '
              f'than {boundary_dedup:.4f} to another boundary particle.')

    if min_separation > 0:
        before = len(fluid_points)
        fluid_points = drop_fluid_too_close_to_boundary(fluid_points, boundary_points, min_separation)
        print(f'Dropped {before - len(fluid_points)} fluid particles within {min_separation:.4f} of a boundary particle.')

    print(f'Generated {len(fluid_points)} fluid particles, {len(boundary_points)} boundary particles.')

    with open(args.out, 'w') as f:
        for x, y in boundary_points:
            f.write(f'{x:.5f} {y:.5f} boundary\n')
        for x, y in fluid_points:
            f.write(f'{x:.5f} {y:.5f} fluid\n')

    print(f'Wrote {args.out}')


if __name__ == '__main__':
    sys.exit(main())

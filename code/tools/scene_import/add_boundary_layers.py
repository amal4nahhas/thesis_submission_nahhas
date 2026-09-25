#!/usr/bin/env python3
"""
Add extra boundary layers to an already-generated scene file, each one offset one more
`--spacing` inward from the last - for SPH's usual requirement of a few particle-thick wall
(so a fluid particle near the wall has proper kernel support on the wall side, and isn't as
easily pushed straight through it). sketch_to_scene.py's boundary tracing only ever produces a
single particle-row along each stroke's edge; this adds the missing depth on top of that,
without needing to re-run the sketch pipeline against the original source image.

Rather than re-tracing from the source sketch (whose exact generation parameters - domain
size, --smooth, etc. - may not be known for an existing scene file), this works directly from
the scene file's own particle positions: it rasterizes the existing fluid and boundary
particles into masks, reconstructs which raster region is "solid" wall material (the traced
boundary particles themselves, plus any fully enclosed gap between them, e.g. the thickness of
a hand-drawn stroke between its own inner and outer traced edges - the same notion of
"enclosed" as src/BoundaryFill.cpp's live-rendering fill) versus "open" (reachable from the
padded canvas border or already fluid-occupied), erodes that solid mask inward by `--spacing`
per extra layer, and traces each eroded mask's contour the same way sketch_to_scene.py traces
the original boundary strokes.

Usage:
    .venv/bin/python add_boundary_layers.py path/to/scene.txt \
        --spacing 0.1 --extra-layers 1 --out path/to/scene.txt
"""
import argparse
import sys

import cv2
import numpy as np

from sketch_to_scene import enforce_min_spacing, resample_polyline_world


def load_scene(path):
    fluid, boundary = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 3:
                continue
            x, y, t = float(parts[0]), float(parts[1]), parts[2]
            (fluid if t == 'fluid' else boundary).append((x, y))
    return fluid, boundary


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('scene', help='Path to an existing scene file (x y type per line).')
    parser.add_argument('--spacing', type=float, required=True,
                        help='Particle spacing, world units - must match the scene.')
    parser.add_argument('--extra-layers', type=int, default=1,
                        help='Additional boundary layers to add on top of the existing one, '
                             'each offset one more --spacing inward. Default 1 (2 layers total).')
    parser.add_argument('--pixels-per-spacing', type=float, default=8.0,
                        help='Internal raster resolution - higher is more accurate but slower '
                             'and uses more memory. Default 8.')
    parser.add_argument('--layer-offset', type=float, default=None,
                        help='World-unit distance to erode inward per layer. Defaults to '
                             '1x --spacing, but a layer offset that close to the minimum '
                             'spacing invariant (see enforce_min_spacing) often gets mostly '
                             'thinned back out by phase misalignment between the new layer\'s '
                             'independently-resampled contour and the original one - try '
                             '1.5-2x --spacing for a layer that reliably survives.')
    parser.add_argument('--out', required=True, help='Output scene file.')
    args = parser.parse_args()

    fluid, boundary = load_scene(args.scene)
    if not boundary:
        sys.exit('No boundary particles found in the input scene.')
    print(f'Loaded {len(fluid)} fluid, {len(boundary)} boundary particles from {args.scene}')

    # Defensive: clean up the input itself in case it predates enforce_min_spacing being
    # wired into scene generation.
    before = len(boundary)
    boundary = enforce_min_spacing(boundary, args.spacing)
    if before != len(boundary):
        print(f'Thinned {before - len(boundary)} pre-existing boundary particles in the input '
              f'itself that were closer than {args.spacing:.4f} to another.')

    all_pts = fluid + boundary
    xs = [p[0] for p in all_pts]
    ys = [p[1] for p in all_pts]
    pad = args.spacing * (args.extra_layers + 2)
    minX, maxX = min(xs) - pad, max(xs) + pad
    minY, maxY = min(ys) - pad, max(ys) + pad

    px_per_unit = args.pixels_per_spacing / args.spacing
    img_w = max(2, int(round((maxX - minX) * px_per_unit)))
    img_h = max(2, int(round((maxY - minY) * px_per_unit)))

    def to_px(x, y):
        return (x - minX) * px_per_unit, (y - minY) * px_per_unit

    def to_world(px, py):
        return px / px_per_unit + minX, py / px_per_unit + minY

    # Rasterize both particle sets as filled discs, same splat-radius idea as
    # src/BoundaryFill.cpp - guarantees a continuous mask along a curve traced at ~spacing
    # arc-length apart, regardless of the curve's local angle relative to the pixel grid.
    splat_radius_px = max(1.0, 0.6 * args.spacing * px_per_unit)

    def rasterize(points):
        mask = np.zeros((img_h, img_w), dtype=np.uint8)
        r = int(round(splat_radius_px))
        for x, y in points:
            px, py = to_px(x, y)
            cv2.circle(mask, (int(round(px)), int(round(py))), r, 255, -1)
        return mask

    wall_mask = rasterize(boundary)
    fluid_mask = rasterize(fluid)

    # "Open" = fluid-occupied, or reachable from the padded border through non-wall space.
    # Everything else (the wall itself, plus any passable pocket that's neither fluid nor
    # reachable from outside - e.g. the gap between a stroke's own inner/outer traced edges)
    # counts as solid wall material to erode into for the new layer(s).
    passable = (wall_mask == 0).astype(np.uint8)
    _, labels = cv2.connectedComponents(passable, connectivity=4)

    border_labels = set(labels[0, :]) | set(labels[-1, :]) | set(labels[:, 0]) | set(labels[:, -1])
    fluid_labels = set(np.unique(labels[fluid_mask > 0]))
    open_labels = (border_labels | fluid_labels) - {0}  # label 0 is the wall (non-passable)

    open_mask = np.isin(labels, list(open_labels))
    solid_mask = np.where(open_mask, 0, 255).astype(np.uint8)
    print(f'Solid (wall) area: {int((solid_mask > 0).sum())} px of {img_w * img_h} total.')

    layer_offset = args.spacing if args.layer_offset is None else args.layer_offset
    erode_radius_px = max(1, int(round(layer_offset * px_per_unit)))
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,
                                       (2 * erode_radius_px + 1, 2 * erode_radius_px + 1))

    new_layers = []
    current = solid_mask
    for layer in range(1, args.extra_layers + 1):
        current = cv2.erode(current, kernel)
        if not np.any(current):
            print(f'Layer {layer}: eroded mask is empty (the wall is thinner than {layer} '
                  f'extra layer(s) everywhere) - stopping.')
            break
        contours, _ = cv2.findContours(current, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
        layer_points = []
        for contour in contours:
            if cv2.arcLength(contour, closed=True) < 4:
                continue  # drop noise specks
            world_pts = [to_world(px, py) for px, py in contour.reshape(-1, 2)]
            layer_points.extend(resample_polyline_world(world_pts, args.spacing, closed=True))

        # This layer's contour has no phase coordination with the original boundary or any
        # layer already added - see enforce_min_spacing's docstring for why that alone is
        # enough to crash the simulation. Seed with everything kept so far so a new layer's
        # points can't land too close to the original boundary OR an earlier layer either.
        before = len(layer_points)
        layer_points = enforce_min_spacing(layer_points, args.spacing, seed_points=boundary + new_layers)
        print(f'Layer {layer}: {len(layer_points)} boundary particles '
              f'({before - len(layer_points)} thinned for being too close to existing boundary).')
        new_layers.extend(layer_points)

    with open(args.out, 'w') as f:
        for x, y in boundary:
            f.write(f'{x:.5f} {y:.5f} boundary\n')
        for x, y in new_layers:
            f.write(f'{x:.5f} {y:.5f} boundary\n')
        for x, y in fluid:
            f.write(f'{x:.5f} {y:.5f} fluid\n')

    print(f'Wrote {args.out}: {len(boundary) + len(new_layers)} boundary particles '
          f'({len(boundary)} original + {len(new_layers)} new), {len(fluid)} fluid.')


if __name__ == '__main__':
    sys.exit(main())

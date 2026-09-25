#!/usr/bin/env python3
"""
Rebuild an existing scene file at a smaller world-space scale, to get a version of the same
shape with far fewer particles.

This does NOT simply thin out the particle list - the simulation's particle spacing (SPACING
in ParticleSettings.cpp) is `sqrt(PARTICLE_VOLUME)`, a fixed global constant, independent of
whatever --spacing a sketch happened to be generated with. Every scene's on-disk particle
spacing has to equal that same fixed value for the physics (kernel support radius, mass) to be
correct, so particle count can only be reduced by shrinking the SHAPE itself while keeping that
spacing fixed - not by dropping particles from the existing lattice.

It rasterizes the input scene's own fluid/boundary particles back into masks (the same
reconstruction trick as add_boundary_layers.py), then reruns sketch_to_scene.py's own
fill_region/trace_boundary_contours against those masks at a smaller --domain-width/--height
(shrinking the shape) while keeping --spacing fixed. Particle count scales roughly with the
square of the scale factor, since a shrink reduces area, not linear extent alone.

Usage:
    .venv/bin/python scale_scene.py example_scenes/underground.txt \
        --spacing 0.1 --scale 0.316 --out example_scenes/underground_100k.txt

    # or target a fluid particle count directly, and let it work out --scale:
    .venv/bin/python scale_scene.py example_scenes/underground.txt \
        --spacing 0.1 --target-fluid 100000 --out example_scenes/underground_100k.txt
"""
import argparse
import math
import sys

import cv2
import numpy as np

from sketch_to_scene import (clearance_mask, drop_fluid_too_close_to_boundary,
                             enforce_min_spacing, fill_region, trace_boundary_contours)


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
                        help='Particle spacing, world units - must match the scene, and stays '
                             'fixed (this is what makes it a shrink, not a re-sample).')
    parser.add_argument('--scale', type=float, default=None,
                        help='Linear scale factor (0-1) to shrink the scene by.')
    parser.add_argument('--target-fluid', type=int, default=None,
                        help='Pick --scale automatically so fluid count lands near this value '
                             '(scale = sqrt(target / current), since count scales with area).')
    parser.add_argument('--pixels-per-spacing', type=float, default=8.0,
                        help='Internal raster resolution relative to the ORIGINAL (unscaled) '
                             'domain - higher preserves more shape detail through the shrink.')
    parser.add_argument('--min-separation', type=float, default=None,
                        help='Same meaning as in sketch_to_scene.py. Defaults to 0.6x --spacing.')
    parser.add_argument('--min-clearance', type=float, default=None,
                        help='Same meaning as in sketch_to_scene.py: excluded spots are left '
                             'open, not turned into boundary. Defaults to 2x --spacing.')
    parser.add_argument('--extra-layers', type=int, default=0,
                        help="Additional boundary layers, each offset one more --layer-offset "
                             "inward, generated from THIS scene's own rich input boundary (not "
                             "the thinned output - see the module docstring's note on why "
                             "add_boundary_layers.py can't do this once a wall is already down "
                             "to a single row). Default 0 (no extra layers).")
    parser.add_argument('--layer-offset', type=float, default=None,
                        help='World-unit distance (in the OUTPUT/scaled domain) to erode '
                             'inward per extra layer. Defaults to 1.5x --spacing - see '
                             'add_boundary_layers.py --layer-offset for why not 1x.')
    parser.add_argument('--boundary-dedup', type=float, default=None,
                        help='Same meaning as in sketch_to_scene.py: a safety net against '
                             'near-duplicate points, not a full --spacing requirement. '
                             'Defaults to 0.15x --spacing; pass 0 to disable.')
    parser.add_argument('--out', required=True, help='Output scene file.')
    args = parser.parse_args()

    fluid, boundary = load_scene(args.scene)
    if not fluid:
        sys.exit('No fluid particles found in the input scene.')
    print(f'Loaded {len(fluid)} fluid, {len(boundary)} boundary particles from {args.scene}')

    if args.scale is None:
        if args.target_fluid is None:
            sys.exit('Pass either --scale or --target-fluid.')
        args.scale = math.sqrt(args.target_fluid / len(fluid))
        print(f'Auto scale for ~{args.target_fluid} fluid particles: {args.scale:.4f}')

    all_pts = fluid + boundary
    xs = [p[0] for p in all_pts]
    ys = [p[1] for p in all_pts]
    minX, minY = min(xs), min(ys)
    domain_w = max(xs) - minX
    domain_h = max(ys) - minY

    px_per_unit = args.pixels_per_spacing / args.spacing
    img_w = max(2, int(round(domain_w * px_per_unit)))
    img_h = max(2, int(round(domain_h * px_per_unit)))

    def to_px(x, y):
        return (x - minX) * px_per_unit, (y - minY) * px_per_unit

    # Splat each particle onto a small disc, same idea as add_boundary_layers.py - keeps a
    # curve traced at ~spacing arc-length apart continuous in the raster no matter its angle.
    splat_radius_px = max(1.0, 0.6 * args.spacing * px_per_unit)

    def rasterize(points):
        mask = np.zeros((img_h, img_w), dtype=np.uint8)
        r = int(round(splat_radius_px))
        for x, y in points:
            px, py = to_px(x, y)
            cv2.circle(mask, (int(round(px)), int(round(py))), r, 255, -1)
        return mask.astype(bool)

    fluid_mask = rasterize(fluid)
    boundary_mask = rasterize(boundary)

    # The raster (img_w x img_h, and its pixel content) is left untouched - only the world-space
    # domain it's interpreted against shrinks. Feeding the same pixels through a smaller
    # domain_w/domain_h at the same --spacing is exactly how sketch_to_scene.py would sample a
    # smaller drawing, so this reuses its fill_region/trace_boundary_contours unchanged.
    new_domain_w = domain_w * args.scale
    new_domain_h = domain_h * args.scale
    print(f'Original domain ~{domain_w:.2f} x {domain_h:.2f} -> scaled domain '
          f'~{new_domain_w:.2f} x {new_domain_h:.2f} (scale {args.scale:.4f}), '
          f'raster {img_w}x{img_h}px.')

    min_separation = args.spacing * 0.6 if args.min_separation is None else args.min_separation
    min_clearance = args.spacing * 2.0 if args.min_clearance is None else args.min_clearance

    if min_clearance > 0:
        clear_mask = clearance_mask(boundary_mask, img_w, img_h, new_domain_w, new_domain_h,
                                    min_clearance)
        omitted = int((fluid_mask & ~clear_mask).sum())
        fluid_mask = fluid_mask & clear_mask
        print(f'Omitted {omitted} fluid px narrower than {min_clearance:.4f} clearance '
              f'(closed/braided gaps that closed up under the shrink) - left open.')

    new_fluid = fill_region(fluid_mask, img_w, img_h, new_domain_w, new_domain_h, args.spacing)
    new_boundary = trace_boundary_contours(boundary_mask, img_w, img_h, new_domain_w,
                                           new_domain_h, args.spacing)

    if args.extra_layers > 0:
        # Erode the RICH input boundary_mask (rasterized straight from this scene's own,
        # still-thick input boundary) rather than the thin single-row trace above - eroding an
        # already-thinned wall has nothing left to eat into (this is exactly why
        # add_boundary_layers.py can't add a layer once a wall is down to one row: it only
        # ever sees the output of a previous run, never the original thickness).
        passable = (boundary_mask == 0).astype(np.uint8)
        _, labels = cv2.connectedComponents(passable, connectivity=4)
        border_labels = (set(labels[0, :]) | set(labels[-1, :]) | set(labels[:, 0]) |
                        set(labels[:, -1]))
        fluid_labels = set(np.unique(labels[fluid_mask > 0]))
        open_labels = (border_labels | fluid_labels) - {0}
        solid_mask = np.where(np.isin(labels, list(open_labels)), 0, 255).astype(np.uint8)

        # Erosion happens in pixel space, but the desired offset is a WORLD-unit distance in
        # the OUTPUT (scaled) domain - convert using the same world-units-per-pixel the output
        # is sampled at (new_domain_w / img_w), not the input's own (domain_w / img_w).
        layer_offset = args.spacing * 1.5 if args.layer_offset is None else args.layer_offset
        px_per_output_unit = img_w / new_domain_w
        erode_radius_px = max(1, int(round(layer_offset * px_per_output_unit)))
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,
                                           (2 * erode_radius_px + 1, 2 * erode_radius_px + 1))

        current = solid_mask
        for layer in range(1, args.extra_layers + 1):
            current = cv2.erode(current, kernel)
            if not np.any(current):
                print(f'Layer {layer}: eroded mask is empty (the wall is thinner than {layer} '
                      f'extra layer(s) everywhere) - stopping.')
                break
            layer_points = trace_boundary_contours(current > 0, img_w, img_h, new_domain_w,
                                                    new_domain_h, args.spacing)
            print(f'Layer {layer}: {len(layer_points)} boundary particles (pre-thinning).')
            new_boundary.extend(layer_points)

    # Every curve above (the base trace, and any extra layer) has no phase coordination with
    # any other, so a light dedup pass removes only genuine near-duplicates - see
    # sketch_to_scene.py's --boundary-dedup for why this isn't a full --spacing requirement.
    # Processed in the order appended above, so the base trace is never thinned by a later,
    # deeper layer - only the reverse.
    boundary_dedup = args.spacing * 0.15 if args.boundary_dedup is None else args.boundary_dedup
    if boundary_dedup > 0:
        before = len(new_boundary)
        new_boundary = enforce_min_spacing(new_boundary, boundary_dedup)
        print(f'Thinned {before - len(new_boundary)} boundary particles (across all layers) '
              f'that landed closer than {boundary_dedup:.4f} to another boundary particle.')

    if min_separation > 0:
        before = len(new_fluid)
        new_fluid = drop_fluid_too_close_to_boundary(new_fluid, new_boundary, min_separation)
        print(f'Dropped {before - len(new_fluid)} fluid particles within {min_separation:.4f} '
              f'of a boundary particle.')

    print(f'Generated {len(new_fluid)} fluid particles, {len(new_boundary)} boundary particles.')

    with open(args.out, 'w') as f:
        for x, y in new_boundary:
            f.write(f'{x:.5f} {y:.5f} boundary\n')
        for x, y in new_fluid:
            f.write(f'{x:.5f} {y:.5f} fluid\n')
    print(f'Wrote {args.out}')


if __name__ == '__main__':
    sys.exit(main())

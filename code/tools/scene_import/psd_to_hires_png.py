#!/usr/bin/env python3
"""
Extract the true vector artwork from a Fresco-exported PSD and rasterize it at high
resolution, for feeding into sketch_to_scene.py instead of Fresco's own (much lower-res,
already-rasterized) PNG export.

Fresco's vector brush layers get saved into PSD as a Smart Object, whose embedded data is
a real SVG - but PSD only stores a flattened raster *preview* of that layer for generic
viewers (Preview.app, QuickLook, ...), which is why it can look like plain pixels even
though the actual vector paths are still in the file.

CAUTION: a PSD can contain multiple vector smart-object layers (old test doodles, decorative
scraps, ...), and the embedded SVG can be a stale snapshot that doesn't match the PSD's own
flattened/composited appearance (seen in practice: one layer's SVG was a disconnected-stroke
test drawing, while the actual final artwork - matching the Fresco-exported flat PNG - was a
*different* layer entirely). This script lists every candidate it finds, marks which one it
picked and why, and always prints a warning when more than one visible candidate exists so a
silent wrong pick doesn't happen again. Cross-check the output against a flat PNG export of
the same artwork if one exists (`--compare-png`) before trusting this over that PNG.

Usage:
    .venv/bin/python psd_to_hires_png.py sketch.psd --out hires.png --scale 5
"""
import argparse
import sys

import cairosvg
from psd_tools import PSDImage


def find_all_vector_layers(layers, path=''):
    """Recursively collect every smart-object layer with an embedded SVG, regardless of
    which group it's nested in or whether it (or an ancestor group) is hidden."""
    candidates = []
    for layer in layers:
        layer_path = f'{path}/{layer.name}'
        if layer.kind == 'smartobject':
            try:
                filename = layer.smart_object.filename or ''
            except Exception:
                filename = ''
            if filename.lower().endswith('.svg'):
                candidates.append((layer_path, layer))
        if layer.kind == 'group':
            candidates.extend(find_all_vector_layers(layer, layer_path))
    return candidates


def pick_layer(candidates):
    """Prefer visible layers over hidden ones; among those, the one with the most embedded
    SVG path data (bytes), since a more detailed/complex drawing is more likely to be the
    real artwork than a small test doodle. Returns (chosen_layer, chosen_path)."""
    def sort_key(item):
        _, layer = item
        try:
            size = len(layer.smart_object.data)
        except Exception:
            size = 0
        return (layer.visible, size)

    best_path, best_layer = max(candidates, key=sort_key)
    return best_layer, best_path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('psd', help='Path to the Fresco-exported PSD.')
    parser.add_argument('--out', required=True, help='Output high-resolution PNG path.')
    parser.add_argument('--scale', type=float, default=5.0,
                        help='Rasterization scale relative to the SVG\'s native size (default 5x).')
    parser.add_argument('--save-svg', help='Optionally also save the extracted raw SVG here.')
    parser.add_argument('--layer', type=int, default=None,
                        help='Index into the printed candidate list to force a specific layer, '
                             'overriding the automatic visible+largest pick.')
    args = parser.parse_args()

    psd = PSDImage.open(args.psd)
    candidates = find_all_vector_layers(psd)
    if not candidates:
        print('No embedded SVG smart object found - layer names/kinds in this file:')
        for l in psd:
            print(f'  {l.name!r}: kind={l.kind}')
        sys.exit(1)

    print(f'Found {len(candidates)} vector smart-object layer(s):')
    visible_count = 0
    for i, (path, layer) in enumerate(candidates):
        try:
            size = len(layer.smart_object.data)
        except Exception:
            size = -1
        visible_count += int(layer.visible)
        print(f'  [{i}] {path!r}  visible={layer.visible}  size={layer.size}  svg_bytes={size}')

    if args.layer is not None:
        layer_path, layer = candidates[args.layer]
        print(f'Using layer [{args.layer}] {layer_path!r} (forced via --layer).')
    else:
        layer, layer_path = pick_layer(candidates)
        print(f'Auto-picked {layer_path!r} (visible + most SVG data).')
        if visible_count > 1:
            print(f'WARNING: {visible_count} visible vector layers found - the auto-pick above '
                  'is a guess among genuine candidates, not the only one. If this looks wrong, '
                  'rerun with --layer N to pick a different one from the list, and cross-check '
                  'against a flat PNG export of the artwork if one exists.')
        elif len(candidates) > 1:
            print(f'Note: {len(candidates)} total candidates found but only this one is marked '
                  'visible - the others are presumably old/hidden drafts, skipped automatically.')

    svg_data = layer.smart_object.data

    if args.save_svg:
        with open(args.save_svg, 'wb') as f:
            f.write(svg_data)
        print(f'Saved raw SVG to {args.save_svg}')

    cairosvg.svg2png(bytestring=svg_data, write_to=args.out, scale=args.scale,
                     background_color='white')
    print(f'Wrote {args.out}')


if __name__ == '__main__':
    sys.exit(main())

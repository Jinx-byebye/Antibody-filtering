#!/usr/bin/env python3
"""Render selected paper panels for the README (requires PyMuPDF)."""
import argparse
from pathlib import Path


def main():
    import pymupdf
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf', type=Path, default=root / 'paper/antibody-filtering.pdf')
    parser.add_argument('--output', type=Path, default=root / 'assets')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    # Coordinates refer to the included 27-page manuscript, in PDF points.
    panels = [
        ('binding-information.png', 9, (78, 65, 520, 289), 'Binding information in AlphaFold'),
        ('regional-confidence.png', 16, (106, 156, 488, 275), 'Within-antigen binder separability of regional'),
        ('ensemble-size.png', 17, (93, 180, 505, 455), 'Fixed-budget recall with smaller'),
    ]
    with pymupdf.open(args.pdf) as doc:
        for name, page, bounds, expected in panels:
            source = doc[page - 1]
            if expected not in source.get_text():
                raise ValueError(f'Page {page} no longer matches the panel; review crop coordinates')
            source.get_pixmap(matrix=pymupdf.Matrix(3, 3), clip=pymupdf.Rect(bounds), alpha=False).save(args.output / name)


if __name__ == '__main__':
    main()

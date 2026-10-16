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
        ('within-antigen-results.png', 8, (70, 288, 525, 743), 'Campaign-specific evaluation of 2,400 candidates'),
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
        # Join the two printed parts of Table 2, keeping one column header.
        parts = [(6, (70, 566, 525, 763)), (7, (70, 110, 525, 376))]
        if 'Table 2.' not in doc[5].get_text() or 'Table 2 continued' not in doc[6].get_text():
            raise ValueError('Table 2 pagination changed; review crop coordinates')
        with pymupdf.open() as joined:
            page = joined.new_page(width=455, height=sum(b[3] - b[1] for _, b in parts))
            offset = 0
            for number, bounds in parts:
                clip = pymupdf.Rect(bounds)
                page.show_pdf_page(pymupdf.Rect(0, offset, clip.width, offset + clip.height),
                                   doc, number - 1, clip=clip)
                offset += clip.height
            page.get_pixmap(matrix=pymupdf.Matrix(3, 3), alpha=False).save(args.output / 'across-antigen-results.png')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Render result illustrations from numerical reference files."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=root / 'assets')
    args = parser.parse_args()
    subprocess.run([sys.executable, str(root/'scripts/render_tables.py'),
                    '--reference', str(root/'reference'), '--output', str(args.output)], check=True)
    subprocess.run([sys.executable, str(root/'scripts/plot_regions.py'),
                    str(root/'reference/regional_confidence'), '--output',
                    str(args.output/'regional-confidence')], check=True)


if __name__ == '__main__':
    main()

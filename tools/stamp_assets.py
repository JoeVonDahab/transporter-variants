"""Add a version stamp to every CSS and script link, so a deploy never mixes new pages with cached assets.

    python3 tools/stamp_assets.py

GitHub Pages lets browsers cache files for ten minutes. Without a stamp, a visitor can get a new
page with yesterday's script for that long. The stamp is a hash of the asset's contents, so it
changes only when the file does.
"""
import hashlib
import re
from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent / "docs"

for page in sorted(DOCS.glob("*.html")):
    html = page.read_text(encoding="utf-8")

    def stamp(m):
        path = m.group(2)
        f = DOCS / path
        if not f.exists():
            return m.group(0)
        h = hashlib.sha256(f.read_bytes()).hexdigest()[:10]
        return f'{m.group(1)}{path}?v={h}"'

    new = re.sub(r'((?:src|href)=")(assets/[^"?]+\.(?:css|js))(?:\?v=[0-9a-f]+)?"', stamp, html)
    if new != html:
        page.write_text(new, encoding="utf-8")
        print("stamped", page.name)

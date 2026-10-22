#!/usr/bin/env python3
"""Rewrite figures/submission/Fig1-8.tif as 8-bit RGB (GeroScience artwork rule).

matplotlib writes RGBA; every alpha value is 255, so dropping the channel changes no pixel.
Resolution (1200 dpi) and LZW compression are kept. figures/manuscript/ is not touched.
"""

from pathlib import Path

import numpy as np
from PIL import Image

SUB = Path(__file__).resolve().parents[1] / "figures" / "submission"


def main():
    for i in range(1, 9):
        path = SUB / f"Fig{i}.tif"
        im = Image.open(path)
        im.load()
        dpi = im.info.get("dpi", (1200, 1200))
        if im.mode == "RGB":
            print(path.name, "already RGB")
            continue
        if im.mode != "RGBA":
            raise SystemExit(f"{path.name}: unexpected mode {im.mode}")
        alpha = np.asarray(im.getchannel("A"))
        if alpha.min() != 255:
            raise SystemExit(f"{path.name}: transparent pixels present")
        rgb = im.convert("RGB")
        if not np.array_equal(np.asarray(rgb), np.asarray(im)[..., :3]):
            raise SystemExit(f"{path.name}: RGB values changed")
        rgb.save(path, compression="tiff_lzw", dpi=dpi)
        print(path.name, rgb.mode, rgb.size, dpi, path.stat().st_size)


if __name__ == "__main__":
    main()

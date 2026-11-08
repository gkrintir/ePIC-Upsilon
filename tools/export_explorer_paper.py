"""Write paper.json for a successful paper replay or export a local figure explorer."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "paper"))
from explorer_results import IMAGES, results_json, write_atomic
from fingerprints import REFERENCE


def export_site(out, replay, receipt, content):
    from PIL import Image

    (out / "data").mkdir(parents=True, exist_ok=True)
    (out / "plots").mkdir(exist_ok=True)
    sources = {"mass": (replay / receipt["figure"]).resolve().with_suffix(".png"),
               "w": replay / "spectra/w/upsilon_es_w_cross_section.png",
               "differential": replay / "spectra/differential/upsilon_es_differential.png"}
    for key, source in sources.items():
        with Image.open(source) as image:
            image.thumbnail((1600, 1600), Image.Resampling.LANCZOS)
            image.convert("RGB").save(out / IMAGES[key], "WEBP", quality=92, method=6)
    write_atomic(out / "data/paper.json", content)
    write_atomic(out / "data/paper.js", "window.PAPER=" + content.rstrip() + ";\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--out", type=Path, help="Export the tables and plot images to a local explorer directory")
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    args = parser.parse_args()
    replay = args.replay.resolve()
    receipt = json.loads((replay / "paper_replay.json").read_text())
    content = results_json(replay, receipt, args.reference)
    if args.out is None:
        write_atomic(replay / "paper.json", content)
        print(f"Wrote {replay / 'paper.json'}")
    else:
        export_site(args.out, replay, receipt, content)
        print(f"Exported paper results to {args.out}")


if __name__ == "__main__":
    main()

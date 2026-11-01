"""Verified, atomic PNG serialization; does not alter any plotted content."""
import io
import os
from pathlib import Path
import tempfile


def install():
    from matplotlib.figure import Figure
    from PIL import Image
    if getattr(Figure.savefig, "_verified_atomic_png", False):
        return
    original = Figure.savefig

    def savefig(fig, filename, *args, **kwargs):
        if not isinstance(filename, (str, os.PathLike)) or Path(filename).suffix.lower() != ".png":
            return original(fig, filename, *args, **kwargs)
        destination = Path(filename)
        with io.BytesIO() as encoded:
            render_kwargs = dict(kwargs)
            render_kwargs["format"] = "png"
            original(fig, encoded, *args, **render_kwargs)
            payload = encoded.getvalue()
        with Image.open(io.BytesIO(payload)) as decoded:
            decoded.load()
        fd, temporary = tempfile.mkstemp(prefix="verified_png_", suffix=".tmp", dir=destination.parent)
        try:
            with os.fdopen(fd, "wb") as stream:
                remaining = memoryview(payload)
                while len(remaining):
                    written = stream.write(remaining)
                    if not written:
                        raise OSError("PNG write made no progress")
                    remaining = remaining[written:]
                stream.flush()
                os.fsync(stream.fileno())
            with Image.open(temporary) as decoded:
                decoded.load()
            os.replace(temporary, destination)
        finally:
            if os.path.exists(temporary): os.unlink(temporary)
    savefig._verified_atomic_png = True
    Figure.savefig = savefig

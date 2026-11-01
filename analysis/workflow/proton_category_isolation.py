"""Contain supported proton-category workflow writes in its separate workspace.

The guard is inactive in ordinary frameworks without .proton_category.json.
External input files remain read-only inputs; they are not redirected or copied.
This is a configuration guard, not a sandbox for arbitrary custom Python code.
"""
import glob
import json
import os
from pathlib import Path
import re


def safe_component(value, label="run ID"):
    value = str(value)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,159}", value):
        raise ValueError(label + " must be one safe path component (letters, digits, '.', '_' or '-')")
    return value


def _expanded(value):
    return Path(os.path.expandvars(os.path.expanduser(str(value)))).absolute()


def _inside(path, parent):
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _no_symlinks(path, root, descendants=False):
    """Inspect lexical components before resolve(), including dangling links."""
    path = _expanded(path)
    if not _inside(path, root):
        raise ValueError("Category path escapes workspace: " + str(path))
    current = root
    if current.is_symlink():
        raise ValueError("Category workspace must not be a symlink: " + str(root))
    for component in path.relative_to(root).parts:
        current = current / component
        if current.is_symlink():
            raise ValueError("Category write/code path must not use symlinks: " + str(current))
    if not _inside(path.resolve(strict=False), root):
        raise ValueError("Resolved category path escapes workspace: " + str(path))
    if descendants and path.is_dir():
        for folder, directories, files in os.walk(path, followlinks=False):
            for name in directories + files:
                child = Path(folder) / name
                if child.is_symlink():
                    raise ValueError("Category output tree contains a symlink: " + str(child))
    return path.resolve(strict=False)


def _output_path(value, category_root, relative_to=None):
    path = Path(os.path.expandvars(os.path.expanduser(str(value))))
    if not path.is_absolute():
        if relative_to is None:
            raise ValueError("Category declared outputs must be absolute after template resolution: " + str(value))
        path = relative_to / path
    generated = category_root / "results/generated"
    if not _inside(path.resolve(strict=False), generated):
        raise ValueError("Category output must remain under " + str(generated) + ": " + str(value))
    checked = _no_symlinks(path, category_root, descendants=True)
    # A glob may match an existing link even when its literal prefix is safe.
    if glob.has_magic(str(path)):
        for match in glob.glob(str(path), recursive=True):
            _no_symlinks(Path(match), category_root, descendants=True)
    return checked


def enforce_category_workspace(framework_root, output_path=None):
    requested = _expanded(framework_root)
    helper_root = Path(__file__).resolve().parents[1]
    marker_root = helper_root if (helper_root / ".proton_category.json").exists() else requested
    marker_path = marker_root / ".proton_category.json"
    if not marker_path.exists():
        return None
    if marker_path.is_symlink():
        raise ValueError("Category marker must not be a symlink")
    marker = json.loads(marker_path.read_text())
    legacy_status = "COMPLETE" if marker.get("package_version", 1) in (1, 2) else None
    if marker.get("installation_status", legacy_status) != "COMPLETE":
        raise ValueError("Proton category installation is incomplete; wait for the creator to finish successfully")
    if marker.get("schema_version") != 1 or marker.get("mode") != "required":
        raise ValueError("Unsupported proton category workspace marker")
    if not isinstance(marker.get("source_root"), str) or not isinstance(marker.get("category_root"), str):
        raise ValueError("Category marker requires source_root and category_root paths")
    source = Path(marker["source_root"])
    category = Path(marker["category_root"])
    if not source.is_absolute() or not category.is_absolute():
        raise ValueError("Category marker root paths must be absolute")
    source, category = source.resolve(strict=False), category.resolve(strict=False)
    if requested.resolve(strict=False) != category or marker_root.resolve(strict=False) != category:
        raise ValueError("This proton-category runner must use its own category framework root: " + str(category))
    if _inside(category, source) or _inside(source, category):
        raise ValueError("Source and category workspace roots must be distinct and disjoint")
    if os.environ.get("YPSILON_PROTON_SELECTION", "required") != "required":
        raise ValueError("This isolated workspace is proton-tagged; use the original framework for no-proton selection")
    # Do not let a copied/symlinked results directory lead back to the source.
    _no_symlinks(category / "results/generated", category)
    if output_path is not None:
        _output_path(output_path, category, relative_to=Path.cwd())
    return marker


def enforce_category_configuration(config, framework_root, resolve_task):
    marker = enforce_category_workspace(framework_root)
    if marker is None:
        return
    root = Path(marker["category_root"]).resolve()
    run = config.get("run", {})
    identifier = safe_component(run.get("id", ""))
    expected_run = root / "results/generated/runs" / identifier
    if _expanded(run.get("dir", "")).resolve(strict=False) != expected_run:
        raise ValueError("Category run.dir must be its own results/generated/runs/<run-id>")
    paths = config.get("paths", {})
    if _expanded(paths.get("project", "")).resolve(strict=False) != root:
        raise ValueError("Category paths.project must be its own framework root")
    if _expanded(paths.get("generated", "")).resolve(strict=False) != root / "results/generated":
        raise ValueError("Category paths.generated must remain in its own workspace")
    _output_path(run["dir"], root)
    for task in config.get("task", []):
        resolved = resolve_task(task, config)
        name = safe_component(resolved.get("name", ""), "task name")
        script = Path(str(resolved["script"]))
        if not script.is_absolute(): script = root / script
        _no_symlinks(script, root)
        cwd = Path(str(resolved["cwd"])) if resolved.get("cwd") else expected_run / "stages" / name
        if not cwd.is_absolute(): cwd = root / cwd
        _no_symlinks(cwd, root)
        for output in resolved.get("outputs", []):
            _output_path(output, root)
        arguments = [str(value) for value in resolved.get("args", [])]
        for index, argument in enumerate(arguments):
            flag, separator, inline = argument.partition("=")
            is_output = (flag == "--out" or flag.startswith("--out-")
                         or (flag.startswith("--output-") and flag not in {"--output-convention", "--output-format", "--output-type", "--output-unit", "--output-units"})
                         or flag in {"--output", "--output-dir", "--output-prefix", "--plot-dir", "--workdir", "--work-dir"})
            if not is_output:
                continue
            if separator:
                value = inline
            elif index + 1 < len(arguments):
                value = arguments[index + 1]
            else:
                raise ValueError("Missing category output path after " + flag)
            _output_path(value, root, relative_to=cwd)

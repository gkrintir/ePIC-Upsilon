#!/usr/bin/env python3
"""List or package the actual CERN production source directory for reconciliation."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

EXTENSIONS = {'.C', '.c', '.cc', '.cpp', '.cxx', '.h', '.hpp', '.sh', '.bash', '.sub',
              '.submit', '.jdl', '.condor', '.py', '.json', '.toml', '.yaml', '.yml', '.in', '.md'}
SKIP_DIRS = {'.git', '.venv', '__pycache__', '.ipynb_checkpoints', 'results', 'output',
             'logs', 'figures', 'build', '.ssh', '.globus', '.config'}


def inventory(source):
    source = Path(source).expanduser().resolve(strict=True)
    if not source.is_dir():
        raise ValueError('--source must be the production-code directory')
    entries = []
    for path in sorted(source.rglob('*')):
        relative = path.relative_to(source)
        if any(part in SKIP_DIRS or part.startswith('.') for part in relative.parts):
            continue
        if path.is_symlink() or not path.is_file():
            continue
        if path.suffix not in EXTENSIONS and path.name not in {'Makefile', 'CMakeLists.txt'}:
            continue
        if path.stat().st_size > 2 * 1024 * 1024:
            continue
        data = path.read_bytes()
        if b'\x00' in data:
            continue
        entries.append({'path': relative.as_posix(), 'bytes': len(data),
                        'sha256': hashlib.sha256(data).hexdigest()})
    return source, entries


def collect(source, output=None):
    root, entries = inventory(source)
    if not entries:
        raise ValueError('No production source files found; check --source')
    if output is not None:
        output = Path(output).expanduser().absolute()
        if output.exists():
            raise FileExistsError('Output already exists; choose a new ZIP name')
        output.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(output, 'x', zipfile.ZIP_DEFLATED) as archive:
            for entry in entries:
                data = (root / entry['path']).read_bytes()
                if hashlib.sha256(data).hexdigest() != entry['sha256']:
                    raise ValueError('Source changed while collecting: ' + entry['path'])
                archive.writestr('production_sources/' + entry['path'], data)
            archive.writestr('production_sources/SOURCE_INVENTORY.json', json.dumps({
                'source_directory': str(root), 'files': entries,
                'scope': 'Code/configuration only; detector ROOT files, outputs, logs and hidden files excluded'}, indent=2)+'\n')
    return root, entries


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('plan', 'pack'))
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    if args.action == 'pack' and args.output is None:
        parser.error('pack requires --output with a fresh ZIP path')
    try:
        _, entries = collect(args.source, args.output if args.action == 'pack' else None)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    for entry in entries:
        print(entry['path'])
    print(str(len(entries)) + ' source files')
    if args.action == 'pack':
        print('Archive: ' + str(args.output))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

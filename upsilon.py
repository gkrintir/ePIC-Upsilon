#!/usr/bin/env python3
"""Entry point for workspace setup, analysis, optimization and status review."""
import sys
sys.dont_write_bytecode = True
from tools.workspace import main
if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env bash
# Compatibility entry point for the original 1S 9x275 commands.
exec bash "$(dirname "$0")/submit_rematch.sh" "$@"

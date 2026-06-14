#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

if [[ ! -f ".venv/bin/activate" ]]; then
  echo "ERROR: .venv/bin/activate nicht gefunden. Erstelle die venv mit: python3 -m venv .venv"
  exit 1
fi

source ".venv/bin/activate"
echo "Venv aktiv: $VIRTUAL_ENV"
exec "$SHELL"


#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
component_root="$repo_root/vendor/dash_draggable"
runtime_root="$repo_root/src/dash_draggable"

# Never let a parent shell's GUANACO_v2 Pixi variables select the old Dash.
unset PIXI_PROJECT_MANIFEST
export PATH="$repo_root/.pixi/envs/default/bin:$PATH"

cd "$component_root"
npm ci
npm run build

mkdir -p "$runtime_root"
cp dash_draggable/*.py "$runtime_root/"
cp dash_draggable/*.js dash_draggable/*.js.map dash_draggable/*.json "$runtime_root/"
cp LICENSE "$runtime_root/LICENSE"

echo "Updated Dash 4-compatible dash-draggable runtime in $runtime_root"

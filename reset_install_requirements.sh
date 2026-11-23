#!/usr/bin/env bash
# Build the checked-in Basilisk tree without replacing its local safety policy.
set -euo pipefail

if [[ $# -ne 0 ]]; then
    printf '%s\n' 'Usage: bash reset_install_requirements.sh (no destructive reset)' >&2
    exit 2
fi
PROJECT_ROOT="$(cd -- "$(dirname -- "$0")" && pwd -P)"
BASILISK="$PROJECT_ROOT/basilisk/src"
if [[ ! -f "$BASILISK/Makefile" ]] ||
   ! grep -Fq 'Basilisk network display is disabled in this repository' "$BASILISK/display.h"; then
    printf '%s\n' 'The hardened vendored Basilisk tree is missing. Restore it from a trusted checkout.' >&2
    exit 1
fi

(
    cd -- "$BASILISK"
    case "$(uname -s)" in
        Darwin) platform_config=config.osx ;;
        *) platform_config=config.gcc ;;
    esac
    # The checkout carries the config link and binaries of the machine that
    # committed them. Keep a custom config; point a stock one at this platform.
    current_config="$(readlink config 2>/dev/null || true)"
    if [[ ! -e config && ! -L config ]] ||
       [[ "$current_config" == config.osx || "$current_config" == config.gcc ]]; then
        ln -sfn "$platform_config" config
    fi
    # Binaries built for another platform cannot run here: rebuild from source.
    if ! ./qcc --version >/dev/null 2>&1; then
        make clean
        rm -f ast/grammar  # not removed by "make clean"; make rebuilds it
    fi
    # Use the default goal. "make all" is forwarded to ast/interpreter, which
    # has no such target, and stops the build.
    make
    ./qcc --version >/dev/null
)
{
    printf 'export PROJECT_ROOT=%q\n' "$PROJECT_ROOT"
    printf 'export BASILISK=%q\n' "$BASILISK"
    printf 'export PATH="$BASILISK:$PATH"\n'
} > "$PROJECT_ROOT/.project_config"
printf '%s\n' 'Built the vendored Basilisk toolchain without network display targets.'

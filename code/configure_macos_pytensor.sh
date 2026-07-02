#!/bin/sh

# Configure PyTensor's native compiler for the architecture of the active
# virtual-environment Python. This file is intended to be sourced from
# .venv/bin/activate.

_lung_smoking_configure_macos_pytensor() {
    if [ "$(uname -s)" != "Darwin" ]; then
        return 0
    fi

    if [ "${_LUNG_SMOKING_PYTENSOR_CONFIGURED:-0}" = "1" ]; then
        return 0
    fi

    if [ -z "${VIRTUAL_ENV:-}" ] || [ ! -x "$VIRTUAL_ENV/bin/python" ]; then
        echo "Error: activate the project virtual environment before configuring PyTensor." >&2
        return 1
    fi

    if ! command -v xcrun >/dev/null 2>&1; then
        echo "Error: xcrun is unavailable. Install the Xcode Command Line Tools." >&2
        return 1
    fi

    if [ ! -x /usr/bin/clang++ ]; then
        echo "Error: /usr/bin/clang++ is unavailable. Install the Xcode Command Line Tools." >&2
        return 1
    fi

    _lung_smoking_python_arch="$("$VIRTUAL_ENV/bin/python" -c \
        'import platform; print(platform.machine())')" || return 1
    case "$_lung_smoking_python_arch" in
        arm64|x86_64)
            ;;
        *)
            echo "Error: unsupported macOS Python architecture: $_lung_smoking_python_arch" >&2
            return 1
            ;;
    esac

    _lung_smoking_sdkroot="$(xcrun --sdk macosx --show-sdk-path)" || return 1
    if [ ! -d "$_lung_smoking_sdkroot" ]; then
        echo "Error: xcrun returned a missing macOS SDK: $_lung_smoking_sdkroot" >&2
        return 1
    fi
    case "$_lung_smoking_sdkroot" in
        *" "*|*",")
            echo "Error: the macOS SDK path contains unsupported whitespace or commas: $_lung_smoking_sdkroot" >&2
            return 1
            ;;
    esac

    case ",${PYTENSOR_FLAGS:-}," in
        *,gcc__cxxflags=*)
            echo "Error: PYTENSOR_FLAGS already defines gcc__cxxflags." >&2
            echo "Unset that entry before activating this project environment." >&2
            return 1
            ;;
    esac

    _lung_smoking_compiler_target="$(/usr/bin/clang++ \
        -arch "$_lung_smoking_python_arch" \
        -isysroot "$_lung_smoking_sdkroot" \
        -dumpmachine)" || return 1
    case "$_lung_smoking_compiler_target" in
        "$_lung_smoking_python_arch"-*)
            ;;
        *)
            echo "Error: Python architecture $_lung_smoking_python_arch does not match compiler target $_lung_smoking_compiler_target." >&2
            return 1
            ;;
    esac

    if [ "${PYTENSOR_FLAGS+x}" = "x" ]; then
        _LUNG_SMOKING_HAD_PYTENSOR_FLAGS=1
        _LUNG_SMOKING_OLD_PYTENSOR_FLAGS="$PYTENSOR_FLAGS"
    else
        _LUNG_SMOKING_HAD_PYTENSOR_FLAGS=0
        _LUNG_SMOKING_OLD_PYTENSOR_FLAGS=
    fi

    _lung_smoking_cxxflags="gcc__cxxflags=-arch $_lung_smoking_python_arch -isysroot $_lung_smoking_sdkroot"
    if [ -n "${PYTENSOR_FLAGS:-}" ]; then
        PYTENSOR_FLAGS="$PYTENSOR_FLAGS,$_lung_smoking_cxxflags"
    else
        PYTENSOR_FLAGS="$_lung_smoking_cxxflags"
    fi

    _LUNG_SMOKING_PYTHON_ARCH="$_lung_smoking_python_arch"
    _LUNG_SMOKING_SDKROOT="$_lung_smoking_sdkroot"
    _LUNG_SMOKING_PYTENSOR_CONFIGURED=1

    export PYTENSOR_FLAGS
    export _LUNG_SMOKING_HAD_PYTENSOR_FLAGS
    export _LUNG_SMOKING_OLD_PYTENSOR_FLAGS
    export _LUNG_SMOKING_PYTHON_ARCH
    export _LUNG_SMOKING_SDKROOT
    export _LUNG_SMOKING_PYTENSOR_CONFIGURED
}


lung_smoking_verify_macos_toolchain() {
    # Purpose: verify that Clang uses the active Python architecture and SDK.
    # Inputs: the environment variables created by the configuration function.
    # Output: status 0 on success; a diagnostic and nonzero status on failure.
    # Assumption: this function is called after activating the project venv.
    if [ "$(uname -s)" != "Darwin" ]; then
        return 0
    fi

    if [ "${_LUNG_SMOKING_PYTENSOR_CONFIGURED:-0}" != "1" ]; then
        echo "Error: the macOS PyTensor toolchain has not been configured." >&2
        return 1
    fi

    _lung_smoking_include_search="$(/usr/bin/clang++ \
        -arch "$_LUNG_SMOKING_PYTHON_ARCH" \
        -isysroot "$_LUNG_SMOKING_SDKROOT" \
        -E -x c++ -v /dev/null 2>&1)" || return 1

    if printf '%s\n' "$_lung_smoking_include_search" | grep -q '^ /usr/local/include$'; then
        echo "Error: Clang is using host headers from /usr/local/include." >&2
        return 1
    fi

    _lung_smoking_python_include="$("$VIRTUAL_ENV/bin/python" -c \
        'import sysconfig; print(sysconfig.get_path("include"))')" || return 1

    if ! printf '#include <Python.h>\n#include <cmath>\nint main() { return 0; }\n' |
        /usr/bin/clang++ \
            -arch "$_LUNG_SMOKING_PYTHON_ARCH" \
            -isysroot "$_LUNG_SMOKING_SDKROOT" \
            -I "$_lung_smoking_python_include" \
            -fsyntax-only -x c++ -; then
        echo "Error: the macOS Python/C++ header compatibility check failed." >&2
        return 1
    fi

    echo "Verified PyTensor compiler: $_LUNG_SMOKING_PYTHON_ARCH using $_LUNG_SMOKING_SDKROOT"
}


_lung_smoking_restore_pytensor_environment() {
    if [ "${_LUNG_SMOKING_HAD_PYTENSOR_FLAGS:-0}" = "1" ]; then
        PYTENSOR_FLAGS="${_LUNG_SMOKING_OLD_PYTENSOR_FLAGS:-}"
        export PYTENSOR_FLAGS
    else
        unset PYTENSOR_FLAGS
    fi

    unset _LUNG_SMOKING_HAD_PYTENSOR_FLAGS
    unset _LUNG_SMOKING_OLD_PYTENSOR_FLAGS
    unset _LUNG_SMOKING_PYTHON_ARCH
    unset _LUNG_SMOKING_SDKROOT
    unset _LUNG_SMOKING_PYTENSOR_CONFIGURED
}


_lung_smoking_wrap_deactivate() {
    if ! command -v deactivate >/dev/null 2>&1; then
        echo "Error: the virtual-environment deactivate function is unavailable." >&2
        return 1
    fi

    _lung_smoking_deactivate_definition="$(typeset -f deactivate)" || return 1
    _lung_smoking_deactivate_definition="$(printf '%s\n' "$_lung_smoking_deactivate_definition" |
        sed '1s/^deactivate[[:space:]]*()/_lung_smoking_original_deactivate ()/')" || return 1
    eval "$_lung_smoking_deactivate_definition"

    deactivate() {
        _lung_smoking_original_deactivate "$@"
        _lung_smoking_deactivate_status=$?

        if [ "${1:-}" != "nondestructive" ]; then
            _lung_smoking_restore_pytensor_environment
            unset -f _lung_smoking_original_deactivate
            unset -f lung_smoking_verify_macos_toolchain
            unset -f _lung_smoking_restore_pytensor_environment
        fi

        return "$_lung_smoking_deactivate_status"
    }
}


_lung_smoking_configuration_was_existing="${_LUNG_SMOKING_PYTENSOR_CONFIGURED:-0}"
if _lung_smoking_configure_macos_pytensor; then
    if [ "$(uname -s)" = "Darwin" ] &&
        { [ "$_lung_smoking_configuration_was_existing" != "1" ] ||
          [ "${_LUNG_SMOKING_ACTIVATION_HOOK:-0}" = "1" ]; }; then
        _lung_smoking_wrap_deactivate
    fi
else
    return 1 2>/dev/null || exit 1
fi

unset -f _lung_smoking_configure_macos_pytensor
unset -f _lung_smoking_wrap_deactivate
unset _lung_smoking_python_arch
unset _lung_smoking_sdkroot
unset _lung_smoking_compiler_target
unset _lung_smoking_cxxflags
unset _lung_smoking_include_search
unset _lung_smoking_python_include
unset _lung_smoking_deactivate_definition
unset _lung_smoking_configuration_was_existing

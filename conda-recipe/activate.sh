#!/usr/bin/env sh

if [ "${BABEL_LIBDIR+x}" = x ]; then
    export AUTOSOLVATE_OLD_BABEL_LIBDIR="${BABEL_LIBDIR}"
    export AUTOSOLVATE_HAD_BABEL_LIBDIR=1
else
    unset AUTOSOLVATE_OLD_BABEL_LIBDIR
    unset AUTOSOLVATE_HAD_BABEL_LIBDIR
fi

for autosolvate_babel_dir in "${CONDA_PREFIX}"/lib/openbabel/*; do
    if [ -d "${autosolvate_babel_dir}" ]; then
        export BABEL_LIBDIR="${autosolvate_babel_dir}"
        break
    fi
done
unset autosolvate_babel_dir

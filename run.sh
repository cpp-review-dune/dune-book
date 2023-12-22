#!/usr/bin/env bash

DIR=$(pwd)/build
if [ -d "$DIR" ]; then
  printf '%s\n' "Removing Lock ($DIR)"
  rm -rf "$DIR"
fi

[ -d build/dune-env ] || python -m venv --system-site-packages build/dune-env
cmake -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON -S . -B build -Wno-dev
cmake --build build
cmake --build build --target doxygen_dune-book
# python -m http.server -d build/doc/doxygen/html 8000

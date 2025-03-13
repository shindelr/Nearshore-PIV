#!/bin/bash

set -e
echo "Setting up Julia packages and environments"

export JULIA_PROJECT="$(pwd)"

julia --project=. -e '
import Pkg;
Pkg.activate(".");
Pkg.instantiate();
'

echo "Environment setup complete"
echo "Manually start a REPL with:   julia --project=. "
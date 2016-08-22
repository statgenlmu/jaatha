#!/bin/bash

script="$(mktemp)"
echo "devtools::test()" > $script
R -d "valgrind --tool=memcheck --leak-check=full" --no-save < $script 2>&1 | less
rm $script

#!/bin/bash

set -euf -o pipefail
LANGUAGE=$1; shift
ACTION=$1; shift

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

DATA_DIR="$(realpath "$SCRIPT_DIR/../outputs")"
RESULT_DIR="$(realpath "$SCRIPT_DIR/../results")"
RESULT_FILE="$RESULT_DIR/$LANGUAGE.json"
DATA_FILE="$DATA_DIR/$LANGUAGE.json"

set -x
hyperfine -w1 -r3 --export-json "$RESULT_FILE" $@ "$ACTION $DATA_FILE"

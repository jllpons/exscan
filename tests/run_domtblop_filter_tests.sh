#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status and treat unset variables as an error
set -euo pipefail

# Variables
projectdir="/Users/jponsramon/Dropbox/PhD/devel/exscan"
expecteddir="tests/data/domtblop/filter"
cwd=$(pwd)

# Function to run tests
run_test() {
    local test_name="$1"
    local command="$2"
    local expected_output="$3"
    local output_file="$4"

    echo "Running $test_name"

    # Execute the command and redirect output to a file
    eval "$command" > "$output_file"

    # Compare the output file with the expected output
    if diff "$output_file" "$expected_output" >/dev/null; then
        echo ""
        echo "** Test passed **"
        echo ""
    else
        echo "Test failed"
        echo "Differences:"
        diff "$output_file" "$expected_output"
    fi

    # Clean up
    rm -f "$output_file"
}

# Test cases
strict_sequence_level_evalue_filter() {
    run_test \
        "domtblop.py filter --seq-evalue 1e-10" \
        "\"$projectdir/bin/domtblop.py\" dummy | \"$projectdir/bin/domtblop.py\" filter --seq-evalue 1e-10 | jq -r '.domain_hits[] | .accession'" \
        "$expecteddir/strict_sequence_level_evalue_filter.txt" \
        "$cwd/output_fullseq_evalue.json"
}
high_seq_score_plus_best_hit() {
    run_test \
        "domtblop.py filter --seq-score 100 --best-hit" \
        "\"$projectdir/bin/domtblop.py\" dummy | \"$projectdir/bin/domtblop.py\" filter --seq-score 100 --best-hit | jq -r '.domain_hits[] | .accession'" \
        "$expecteddir/high_seq_score_plus_best_hit.txt" \
        "$cwd/output_high_seq_score_plus_best_hit.json"
}
seq_bias() {
    run_test \
        "domtblop.py filter --seq-bias 1.0" \
        "\"$projectdir/bin/domtblop.py\" dummy | \"$projectdir/bin/domtblop.py\" filter --seq-bias 1.0 | jq -r '.domain_hits[] | .accession'" \
        "$expecteddir/seq_bias.txt" \
        "$cwd/output_seq_bias.json"
}
domain_independent_evalue_and_best_hit() {
    run_test \
        "domtblop.py filter --domain-evalue 1e-10 --best-hit" \
        "\"$projectdir/bin/domtblop.py\" dummy | \"$projectdir/bin/domtblop.py\" filter --dom-ievalue 1e-10 --best-hit | jq -r '.domain_hits[] | .accession'" \
        "$expecteddir/domain_independent_evalue_and_best_hit.txt" \
        "$cwd/output_domain_independent_evalue_and_best_hit.json"
}
domain_conditional_evalue() {
    run_test \
        "domtblop.py filter --domain-cevalue 1e-04" \
        "\"$projectdir/bin/domtblop.py\" dummy | \"$projectdir/bin/domtblop.py\" filter --dom-cevalue 1e-04 | jq -r '.domain_hits[] | .accession'" \
        "$expecteddir/domain_conditional_evalue.txt" \
        "$cwd/output_domain_conditional_evalue.json"
}
minimum_alignment_length() {
    run_test \
        "domtblop.py filter --min-alignment-length 10" \
        "\"$projectdir/bin/domtblop.py\" dummy | \"$projectdir/bin/domtblop.py\" filter --min-alignment-length 10 | jq -r '.domain_hits[] | .accession'" \
        "$expecteddir/minimum_alignment_length.txt" \
        "$cwd/output_minimum_alignment_length.json"
}

# Run test cases
strict_sequence_level_evalue_filter
high_seq_score_plus_best_hit
seq_bias
domain_independent_evalue_and_best_hit
domain_conditional_evalue
minimum_alignment_length


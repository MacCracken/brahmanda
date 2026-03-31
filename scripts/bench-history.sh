#!/usr/bin/env bash
# Run criterion benchmarks and archive results.
#
# Usage: ./scripts/bench-history.sh [-- extra criterion args]
#
# Results are saved to target/criterion/ (default criterion location).
# The HTML report is at target/criterion/report/index.html.

set -euo pipefail

BENCH_DIR="target/criterion"

echo "=== Brahmanda Benchmarks ==="
echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ""

cargo bench "$@"

echo ""
echo "Results saved to ${BENCH_DIR}/"
echo "HTML report: ${BENCH_DIR}/report/index.html"

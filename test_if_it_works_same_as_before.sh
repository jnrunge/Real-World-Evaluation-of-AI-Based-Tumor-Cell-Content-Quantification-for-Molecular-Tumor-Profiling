#!/bin/bash

set -euo pipefail

# Function to handle errors
error_exit() {
    echo "ERROR: $1" >&2
    exit 1
}

# Generate unique filename with current date
generate_unique_filename() {
    local base_name="$(date +%Y%m%d)"
    local filename="${base_name}-output.md5"
    local counter=1
    
    while [[ -f "$filename" ]]; do
        filename="${base_name}_${counter}-output.md5"
        ((counter++))
    done
    
    echo "$filename"
}

# Generate md5 checksums for directory
generate_md5() {
    local output_file="$1"
    if [[ -d "./output" ]]; then
        find ./output -type f -exec md5sum {} \; | sort -k 2 > "$output_file" || error_exit "Failed to generate md5 checksums"
        echo "Generated md5 checksums: $output_file"
    else
        error_exit "./output directory does not exist"
    fi
}

# Main execution
echo "=== Starting pipeline comparison test ==="

# Generate md5 for current output
BEFORE_MD5=$(generate_unique_filename)
echo "Generating checksums for current output..."
generate_md5 "$BEFORE_MD5"

# Move current output to backup
BACKUP_DIR="./output_backup_$(date +%Y%m%d_%H%M%S)"
echo "Moving ./output to $BACKUP_DIR..."
mv ./output "$BACKUP_DIR" || error_exit "Failed to move ./output directory"

# Run nextflow pipeline
echo "Running nextflow pipeline..."
nextflow run main.nf || error_exit "Nextflow pipeline failed"

# Generate md5 for new output
AFTER_MD5="output.md5"
echo "Generating checksums for new output..."
generate_md5 "$AFTER_MD5"

# Compare md5 files
echo ""
echo "=== Comparing outputs ==="
if diff "$BEFORE_MD5" "$AFTER_MD5" > /dev/null; then
    echo "SUCCESS: Outputs are identical!"
else
    echo "DIFFERENCES FOUND:"
    echo ""
    diff -u "$BEFORE_MD5" "$AFTER_MD5" || true
    echo ""
    # Use temp files instead of process substitution
    TMPDIR="$(mktemp -d)" || error_exit "Failed to create temp dir"
    trap 'rm -rf "$TMPDIR"' EXIT
    BEFORE_LIST="$TMPDIR/before.list"
    AFTER_LIST="$TMPDIR/after.list"
    awk '{print $2}' "$BEFORE_MD5" | sort > "$BEFORE_LIST"
    awk '{print $2}' "$AFTER_MD5" | sort > "$AFTER_LIST"

    echo "Files only in before:"
    comm -23 "$BEFORE_LIST" "$AFTER_LIST"
    echo ""
    echo "Files only in after:"
    comm -13 "$BEFORE_LIST" "$AFTER_LIST"
    echo ""
    echo "Files with different checksums:"
    comm -12 "$BEFORE_LIST" "$AFTER_LIST" | while IFS= read -r file; do
        before_hash="$(grep -F " $file" "$BEFORE_MD5" | awk '{print $1}')"
        after_hash="$(grep -F " $file" "$AFTER_MD5" | awk '{print $1}')"
        if [[ -n "$before_hash" && -n "$after_hash" && "$before_hash" != "$after_hash" ]]; then
            echo "  $file"
            echo "    Before: $before_hash"
            echo "    After:  $after_hash"
        fi
    done
fi

echo ""
echo "=== Test complete ==="
echo "Before MD5 file: $BEFORE_MD5"
echo "After MD5 file: $AFTER_MD5"
echo "Backup directory: $BACKUP_DIR"

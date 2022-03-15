#!/usr/bin/env python3
import sys
quality_char = "!"
if len(sys.argv) >= 2:
    quality_char = sys.argv[1]
    assert(len(quality_char) == 1)
sys.stderr.write("[{}] Updating all missing quality scores with '{}'\n".format(sys.argv[0], quality_char))

updated_count = 0
for line in sys.stdin:
    if line.startswith("@"): 
        sys.stdout.write(line)
        continue
    parts = line.split("\t")
    if len(parts) < 11:
        sys.stdout.write(line)
        continue
    if parts[10] == "*":
        parts[10] = quality_char * len(parts[9])
        sys.stdout.write("\t".join(parts))
        updated_count += 1
sys.stderr.write("[{}] Updated {} records\n".format(sys.argv[0], updated_count))

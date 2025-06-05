
#!/bin/bash

# Usage: ./filter_significant_events.sh SE.MATS.JC.txt output_filtered.txt

INPUT=$1
OUTPUT=$2

awk 'NR==1 || ($20 < 0.05 && ($23 > 0.05 || $23 < -0.05))' $INPUT > $OUTPUT

echo "Filtered events written to $OUTPUT"

#!/usr/bin/env awk -f

# Convert diffHic output to 'interval' format for UCSC browser.
# Append output to one of two output files depend on the sign logFC.

# Positional arguments
#   max
#       Maximum absolute logFC value of input file. Used to calculate score.
#   out_up
#       Output file path to store positive logFC differential interactions.
#   out_down
#       Output file path to store negative logFC differential interactions.

# Example usage
#   awk -v max="${max}" -v out_up="${out_up}" \
#       -v out_down="${out_down}" -f interval_convert.awk diffhic.txt


# Function to compute absolute value.
function abs(x) {
    return ((x < 0.0) ? -x : x)
}

# Modify output field seperator to be a tab.
BEGIN {
    OFS="\t";
}

# Calculate score, scaled to 1000, against a user provided 'max'.
{score = (abs($7)*1000)/max}

# Set colour and appropriate output file.
{if($7 < 0)
    {colour = "#FF0000"; out = out_down}
else
    {colour = "#0000FF"; out = out_up}}

# Output in interval format.
{print "chr"$1, $2, $6, ".", int(score), $7, ".", colour,
       "chr"$1, $2, $3, ".", ".", "chr"$4, $5, $6, ".", "." >> out}

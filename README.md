# SparCC
# Language: Python
# Input: CSV (abundances)
# Output: CSV (correlations)
# Dependency: Ensemble Plugin (Available at: https://github.com/movingpictures83/Ensemble)
# Tested with: PluMA 1.0, Python 3.6

PluMA plugin that accepts an input file of abundances in CSV format, with rows
representing samples and columns representing entities.  Entry (i, j) is then the abundance
of entity j in sample i.

The plugin produces an output file of correlations in CSV format computed using
SparCC (Friedman and Alm, 2012) where both rows and columns represent entities and
entry (i, j) is the correlation between entity i and entity j.

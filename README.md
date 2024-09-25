# MRR2_interference_line_removal

Python program to perform the fuzzy-based interference line removal algorithm for MRR-2 data

## Syntax
```python <script>.py <config>.yaml <input_filename> <output_directory>```

## Example
```python main_mrr_interference_removal.py config_mrrqc.yaml /foo/bar/0531.nc /path/to/save/output/without_yearmonth/```

## Important Notes
 - Input file must be post-processed by Maahn and Kollias (2012) algorithm (i.e. IMProToo).
 - This altorithm assumed the rain event. Therefore, the algorithm does not process de-aliased fields, which often gives unreliable value for rain. The users should read '~~_noDA' fields.

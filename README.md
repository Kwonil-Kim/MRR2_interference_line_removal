# 🌧️ MRR2 interference line removal algorithm

Python program to perform the fuzzy-based interference line removal algorithm for MRR-2 📡 data

## ⚙️ Usage
```python <script>.py <config>.yaml <input_filename> <output_directory>```

## Example
```python main_mrr_interference_removal.py config_mrrqc.yaml /foo/bar/0531.nc /path/to/save/output/without_yearmonth/```

## 📝 Important Notes
 - 🚨 Input file must be post-processed by Maahn and Kollias (2012) algorithm (i.e. IMProToo).
 - ☔ This altorithm assumed the rain event. Therefore, the algorithm does not process de-aliased fields, which often gives unreliable value for rain. The users should read '~~_noDA' fields.

## 📖 References
If you found this code helpful, we’d be thrilled to see our publication cited! 🙌 Many thanks!

Kim, K.; Lee, G. A Fuzzy-Logic-Based Approach for Eliminating Interference Lines in Micro Rain Radar (MRR-2). Remote Sens. 2024, 16, 3965. https://doi.org/10.3390/rs16213965

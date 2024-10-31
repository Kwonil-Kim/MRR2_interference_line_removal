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

Kim, K., and G. Lee, 2024: A Fuzzy-Logic-Based Approach for Eliminating Interference Lines in Micro Rain Radar (MRR-2). Remote Sens. 16, 3965, https://doi.org/10.3390/rs16213965.

We also recommended checking out the MK12 algorithm:

Maahn, M., and P. Kollias, 2012: Improved Micro Rain Radar snow measurements using Doppler spectra post-processing, Atmos. Meas. Tech., 5, 2661–2673, https://doi.org/10.5194/amt-5-2661-2012.

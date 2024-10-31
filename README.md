# 🌧️ MRR2 interference line removal algorithm

Python program to perform the fuzzy-based interference line removal algorithm for MRR-2 📡 data

![Example image](https://www.mdpi.com/remotesensing/remotesensing-16-03965/article_deploy/html/images/remotesensing-16-03965-g005.png)
*Figure from Kim and Lee (2024)*

## ⚙️ Usage
```python <script>.py <config>.yaml <input_filename> <output_directory>```

## Example
```python main_mrr_interference_removal.py config_mrrqc.yaml /foo/bar/0531.nc /path/to/save/output/without_yearmonth/```
https://stonybrook365-my.sharepoint.com/:b:/r/personal/kwonil_kim_stonybrook_edu/Documents/_%EA%B9%80%EA%B6%8C%EC%9D%BC_%EC%8B%A4%EC%A0%81/1_%EB%85%BC%EB%AC%B8/2024_Kwonil_MRRQC/4_proofread/TS_MF.pdf?csf=1&web=1&e=DELIZ3
## 📝 Important Notes
 - 🚨 Input file must be post-processed by Maahn and Kollias (2012) algorithm (i.e. IMProToo).
 - ☔ This altorithm assumed the rain event. Therefore, the algorithm does not process de-aliased fields, which often gives unreliable value for rain. The users should read '~~_noDA' fields.
 - 🔧 If interference lines are found in the mid or lower levels rather than the upper levels, we recommend adjusting `min_hgt_idx` to a lower value accordingly.

## 📖 References
If you found this code helpful, we’d be thrilled to see our publication cited! 🙌 Many thanks!

Kim, K., and G. Lee, 2024: A Fuzzy-Logic-Based Approach for Eliminating Interference Lines in Micro Rain Radar (MRR-2). Remote Sens. 16, 3965, https://doi.org/10.3390/rs16213965.

We also recommended checking out the MK12 algorithm:

Maahn, M., and P. Kollias, 2012: Improved Micro Rain Radar snow measurements using Doppler spectra post-processing, Atmos. Meas. Tech., 5, 2661–2673, https://doi.org/10.5194/amt-5-2661-2012.

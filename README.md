
# ABI geneexpression data

Script to download geneexpression data from the Allen Mouse Brain data portal [ABI](http://mouse.brain-map.org/). Script will query the database, download available data, convert form raw/mhd format to NIfTI
and registers data to dsurqec??


# Citation Notice

Citation guidelines for use of the Allen mouse brain atlas can be found at: https://alleninstitute.org/legal/citation-policy/

Script mhd_utils_3d.py:  
Author: bjian,Price Jackson  
Source:  
(slightly adapted to be compatible with newer python versions)


# Usage

In order to create a new version of the ABI_geneexpression data package, simply navigate to the root directory of this repository and run:

```
python -v 0.5 abi_geneexpression.py
```

This will create archives with the newest files fetched from upstream and processed according to the instructions standardized in this package.
The version suffix will be `0.5` (as per the `-v 0.5` parameter).

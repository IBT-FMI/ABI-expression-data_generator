
# ABI Expression Data

Script to download expression data from the Allen Mouse Brain data portal [ABI](http://mouse.brain-map.org/). Script will query the database, download available data, convert from raw/mhd format to NIfTI
and registers data.


# ABI Expression Data Package Releases

Current recommended release in bold typeface:

* **[ABI-expression-data-0.1.1.tar.xz](http://chymera.eu/distfiles/ABI-expression-data-0.1.1.tar.xz)** \[[checksum](http://chymera.eu/distfiles/ABI-expression-data-0.1.1.tar.xz)\]
* [ABI-expression-data-0.1.tar.xz](http://chymera.eu/distfiles/ABI-expression-data-0.1.tar.xz) \[[checksum](http://chymera.eu/distfiles/ABI-expression-data-0.1.tar.xz)\]

# Citation Notice

Citation guidelines for use of the Allen mouse brain atlas can be found at: https://alleninstitute.org/legal/citation-policy/

Script `mhd_utils_3d.py`:  
Author: bjian,Price Jackson  
Source: https://sites.google.com/site/pjmedphys/scripts   
(slightly adapted for compatibility with newer Python versions)


# Usage

In order to create a new version of the ABI-expression data package, simply navigate to the root directory of this repository and run:

```
python -v 0.5 abi_expression.py
```

This will create archives with the newest files fetched from upstream and processed according to the instructions standardized in this package.
The version suffix will be `0.5` (as per the `-v 0.5` parameter).

# standR 1.17.1

* Add `readGeoMxFromNanoStringGeoMxSet()` to convert
  `GeomxTools::readNanoStringGeoMxSet()` / `NanoStringGeoMxSet` objects from
  DCC/PKC workflows into `SpatialExperiment` objects for standR.

* Improve `prepareSpatialDecon()` so negative probes removed by
  `readGeoMx(rmNegProbe = TRUE)` can be read from `metadata(spe)$NegProbes`.

# standR 1.3.9

* Add function _prepareSpatialDecon_ to help using R package [SpatialDecon](https://bioconductor.org/packages/release/bioc/html/SpatialDecon.html) after using standR to preprocess GeoMx data.

* New RUV-4: now using the RUV-4 from standR allow you to perform rank-based analysis such as gene-set scoring with the RUV-4-normalised count.

* New vignette - A quick start guide to the standR package.

* Many bugs fixed.

# standR 0.99.0

* First release of the package.

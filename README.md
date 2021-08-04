[![GS-Frame](https://img.shields.io/badge/github-GeoStat_Framework-468a88?logo=github&style=flat)](https://github.com/GeoStat-Framework)
[![Gitter](https://badges.gitter.im/GeoStat-Examples/community.svg)](https://gitter.im/GeoStat-Examples/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5159577.svg)](https://doi.org/10.5281/zenodo.5159577)

# Impact of connectivity on contaminant transport

Here we investigate the impact of different types of transmissivity
connectivity on contaminant transport.

To generate heterogeneous transmissivity fields with connected high and
low regions, we make use of the transformation presented by Zinn & Harvey:

> Zinn, B. and Harvey, C. F.:
> When good statistical models of aquifer heterogeneity go bad:
> A comparison of flow, dispersion, and mass transfer in connected and
> multivariate Gaussian hydraulic conductivity fields,
> Water Resour. Res., 39, 2003.

## Structure

The workflow is organized by the following structure:
- `src/`
  - `00_tracer_test_conn.py` - field generation and transport simulation
  - `01_plot.py` - detailed analysis plot
  - `02_plot.py` - simple analysis plot
- `results/` - all produced results


## Python environment

Main Python dependencies are stored in `requirements.txt`:
```
gstools==1.3.1
ogs5py=1.1.1
matplotlib
```

## Contact

You can contact us via <info@geostat-framework.org>.


## License

MIT © 2021

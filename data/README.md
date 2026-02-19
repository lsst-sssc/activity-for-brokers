# Example Cometary Lightcurves 

This repository contains `.csv` files with photometric data for selected cases of interest. 

## Column Headers

Each `.csv` file includes the following columns:

- **mjd**: Observation time in Modified Julian Date
- **rh**: Heliocentric distance of the target (AU)
- **delta**: Geocentric distance (AU)
- **alpha**: Phase angle (degrees)
- **aper_arcsec**: Aperture size (arcseconds)
- **mag**: Observed magnitude
- **mag_err**: Magnitude uncertainty
- **filter**: Telescope filter used
- **query_datetime (optional)**: date/time queries (e.g., JPL, FINK) were done

## Example Lightcurves

### Synthetic Data

- **LPC with y=3**: [`LPC_y3.csv`](./LPC_y3.csv)  
  _Description_: This dataset is generated using the geometries of comet **C/2021 S3** and the **Hy model**, with added scatter (σ = 0.05).
- **LPC with a=0.3, b=-1**: [`LPC_0.3a-1b.csv`](./LPC_0.3a-1b.csv)  
  _Description_: This dataset is generated using the geometries of comet **C/2021 S3** and the **Hab model**, with added scatter (σ = 0.05).

### Real Data Lightcurves

The filenames follow the format: `<objname>_<surveyname>.csv`

- **Typical Jupiter-Family Comet (JFC)**: [`19P_ATLAS.csv`](./19P_ATLAS.csv)
- **Seasonal Effects**: [`67P_ATLAS.csv`](./67P_ATLAS.csv)
- **High Signal-to-Noise Ratio (SNR)**: [`104P_ATLAS.csv`](./104P_ATLAS.csv)
- **Solar Conjunction and Multiple Apparitions**: [`117P_ATLAS.csv`](./117P_ATLAS.csv)
- **Low SNR**: [`179P_ATLAS.csv`](./179P_ATLAS.csv)
- **Nucleus Detection**: [`459P_ATLAS.csv`](./459P_ATLAS.csv)
- **Typical Long-Period Comet (LPC)**: [`C2021S3_LOOK.csv`](./C2021S3_LOOK.csv)
- **Fading Event**: [`C2021Y1_LOOK.csv`](./C2021Y1_LOOK.csv)
- **Unique Case (Gault)**: [`Gault_ZTF.csv`](./Gault_ZTF.csv), [`Gault_ZTF_FINK.csv`](./Gault_ZTF_FINK.csv)
- **Outbursting Active Asteroid**: ['2016UU121_ZTF_FINK.csv'](./2016UU121_ZTF_FINK.csv)
- **Iconic Active Asteroid**: ['7968_ZTF_FINK.csv'](./7968_ZTF_FINK.csv)
- **Active Centaur**: ['Chiron_ZTF_FINK.csv'](./Chiron_ZTF_FINK.csv)
- **Special Active Centaur**: ['Echeclus_ZTF_FINK.csv'](./Echeclus_ZTF_FINK.csv)

### Outburst Events

- **46P/Wirtanen**: [`46P_ZTF.csv`](./46P_ZTF.csv), [`46P_ZTF_FINK.csv`](./46P_ZTF_FINK.csv)
- **97P/Metcalf-Brewington**: [`97P_ATLAS.csv`](./97P_ATLAS.csv)
- **99P/Kowal 1**: [`99P_ATLAS.csv`](./99P_ATLAS.csv)
- **118P/Shoemaker-Levy 4**: [`118P_ATLAS.csv`](./118P_ATLAS.csv)
- **285P/LINEAR**: [`285P_ATLAS.csv`](./285P_ATLAS.csv)
- **382P/Larson** [`382P_ATLAS.csv`](./382P_ATLAS.csv)

## Data Sources

The data presented in this repository comes from the following publications:

- Ye et al., 2019 (ZTF - Gault)
- Kelley et al., 2021 (ZTF - 46P)
- Holt et al., 2024 (LOOK)
- Gillan et al., 2024, Gillan et al., 2025 (ATLAS)
- Chandler et al., 2024a (2016 UU121), 2024b (TailNet0, *in prep*)
- FINK Broker


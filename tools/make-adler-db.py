"""
Make an Adler database from real light curves.

Warning:

    * The reference frame has not been checked.
    * The time scales have not been checked.
    * This script will overwrite existing data in the database.
    * ATLAS c and o are renamed g and r.

See data/README.md for input data format.


  target     ssObjectId
--------- ----------------
      19P 6592633514653240
      46P 2358867111535997
      67P  940303576399572
      97P 8285973462698227
      99P 1290030297348764
     104P 5113060619462226
     117P 3575262679851537
     118P   69165524166199
     179P 2462443667160976
     285P 3138720172365267
     382P  436637127094172
     459P 3511020769536744
C/2021 S3  798904090225145
C/2021 Y1 1709597069217518
     6478 1890843518253925


Output format based on DP1:
CREATE TABLE IF NOT EXISTS "DiaSource" (
"apFlux" REAL,
  "apFlux_flag" INTEGER,
  "apFlux_flag_apertureTruncated" INTEGER,
  "apFluxErr" REAL,
  "band" TEXT,
  "bboxSize" INTEGER,
  "centroid_flag" INTEGER,
  "coord_dec" REAL,
  "coord_ra" REAL,
  "dec" REAL,
  "decErr" REAL,
  "detector" INTEGER,
  "diaObjectId" INTEGER,
  "diaSourceId" INTEGER,
  "dipoleAngle" REAL,
  "dipoleChi2" REAL,
  "dipoleFitAttempted" INTEGER,
  "dipoleFluxDiff" REAL,
  "dipoleFluxDiffErr" REAL,
  "dipoleLength" REAL,
  "dipoleMeanFlux" REAL,
  "dipoleMeanFluxErr" REAL,
  "dipoleNdata" INTEGER,
  "extendedness" REAL,
  "forced_PsfFlux_flag" INTEGER,
  "forced_PsfFlux_flag_edge" INTEGER,
  "forced_PsfFlux_flag_noGoodPixels" INTEGER,
  "isDipole" INTEGER,
  "ixx" REAL,
  "ixxPSF" REAL,
  "ixy" REAL,
  "ixyPSF" REAL,
  "iyy" REAL,
  "iyyPSF" REAL,
  "midpointMjdTai" REAL,
  "parentDiaSourceId" INTEGER,
  "pixelFlags" INTEGER,
  "pixelFlags_bad" INTEGER,
  "pixelFlags_cr" INTEGER,
  "pixelFlags_crCenter" INTEGER,
  "pixelFlags_edge" INTEGER,
  "pixelFlags_injected" INTEGER,
  "pixelFlags_injected_template" INTEGER,
  "pixelFlags_injected_templateCenter" INTEGER,
  "pixelFlags_injectedCenter" INTEGER,
  "pixelFlags_interpolated" INTEGER,
  "pixelFlags_interpolatedCenter" INTEGER,
  "pixelFlags_nodata" INTEGER,
  "pixelFlags_nodataCenter" INTEGER,
  "pixelFlags_offimage" INTEGER,
  "pixelFlags_saturated" INTEGER,
  "pixelFlags_saturatedCenter" INTEGER,
  "pixelFlags_streak" INTEGER,
  "pixelFlags_streakCenter" INTEGER,
  "pixelFlags_suspect" INTEGER,
  "pixelFlags_suspectCenter" INTEGER,
  "psfChi2" REAL,
  "psfFlux" REAL,
  "psfFlux_flag" INTEGER,
  "psfFlux_flag_edge" INTEGER,
  "psfFlux_flag_noGoodPixels" INTEGER,
  "psfFluxErr" REAL,
  "psfNdata" INTEGER,
  "ra" REAL,
  "ra_dec_Cov" REAL,
  "raErr" REAL,
  "reliability" REAL,
  "scienceFlux" REAL,
  "scienceFluxErr" REAL,
  "shape_flag" INTEGER,
  "shape_flag_no_pixels" INTEGER,
  "shape_flag_not_contained" INTEGER,
  "shape_flag_parent_source" INTEGER,
  "snr" REAL,
  "ssObjectId" INTEGER,
  "time_processed" TEXT,
  "trail_flag_edge" INTEGER,
  "trailAngle" REAL,
  "trailDec" REAL,
  "trailFlux" REAL,
  "trailLength" REAL,
  "trailRa" REAL,
  "visit" INTEGER,
  "x" REAL,
  "xErr" REAL,
  "y" REAL,
  "yErr" REAL
);
CREATE TABLE IF NOT EXISTS "SSSource" (
"diaSourceId" INTEGER,
  "eclipticBeta" REAL,
  "eclipticLambda" REAL,
  "galacticB" REAL,
  "galacticL" REAL,
  "heliocentricDist" REAL,
  "heliocentricVX" REAL,
  "heliocentricVY" REAL,
  "heliocentricVZ" REAL,
  "heliocentricX" REAL,
  "heliocentricY" REAL,
  "heliocentricZ" REAL,
  "phaseAngle" REAL,
  "residualDec" REAL,
  "residualRa" REAL,
  "ssObjectId" INTEGER,
  "topocentricDist" REAL,
  "topocentricVX" REAL,
  "topocentricVY" REAL,
  "topocentricVZ" REAL,
  "topocentricX" REAL,
  "topocentricY" REAL,
  "topocentricZ" REAL
);
CREATE TABLE IF NOT EXISTS "SSObject" (
"discoverySubmissionDate" REAL,
  "numObs" INTEGER,
  "ssObjectId" INTEGER
);
CREATE TABLE IF NOT EXISTS "MPCORB" (
"e" REAL,
  "epoch" REAL,
  "incl" REAL,
  "mpcDesignation" TEXT,
  "mpcH" REAL,
  "node" REAL,
  "peri" REAL,
  "q" REAL,
  "ssObjectId" INTEGER,
  "t_p" REAL
);

"""

import os
import hashlib
import sqlite3
import argparse

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table, vstack
from astroquery.jplhorizons import Horizons
from sbpy.data import Names

# designation: file name; where designation is resolvable by Horizons
targets = {
    "19P": "19P_ATLAS.csv",
    "46P": "46P_ZTF.csv",
    "67P": "67P_ATLAS.csv",
    "97P": "97P_ATLAS.csv",
    "99P": "99P_ATLAS.csv",
    "104P": "104P_ATLAS.csv",
    "117P": "117P_ATLAS.csv",
    "118P": "118P_ATLAS.csv",
    "179P": "179P_ATLAS.csv",
    "285P": "285P_ATLAS.csv",
    "382P": "382P_ATLAS.csv",
    "459P": "459P_ATLAS.csv",
    "C/2021 S3": "C2021S3_LOOK.csv",
    "C/2021 Y1": "C2021Y1_LOOK.csv",
    # "7968": "7968_ZTF_FINK.csv",  # FINK tables are missing filter
    # "2016 UU121": "2016UU121_ZTF_FINK.csv",
    # "2060": "Chiron_ZTF_FINK.csv",
    # "60558": "Echeclus_ZTF_FINK.csv",
    "6478": "Gault_ZTF.csv",
}

# replacement filter names (for ATLAS)
filter = {
    "c": "g",
    "o": "r",
}


def to_id(s):
    """Convert a string to stable ID (integer) value using an MD5 hash."""

    hash = hashlib.md5()
    hash.update(s.encode())
    return int.from_bytes(hash.digest()) % 10**16


diasource_columns = [
    "apFlux_flag",
    "apFlux_flag_apertureTruncated",
    "bboxSize",
    "centroid_flag",
    "coord_dec",
    "coord_ra",
    "dec",
    "decErr",
    "detector",
    "diaObjectId",
    "dipoleAngle",
    "dipoleChi2",
    "dipoleFitAttempted",
    "dipoleFluxDiff",
    "dipoleFluxDiffErr",
    "dipoleLength",
    "dipoleMeanFlux",
    "dipoleMeanFluxErr",
    "dipoleNdata",
    "extendedness",
    "forced_PsfFlux_flag",
    "forced_PsfFlux_flag_edge",
    "forced_PsfFlux_flag_noGoodPixels",
    "isDipole",
    "ixx",
    "ixxPSF",
    "ixy",
    "ixyPSF",
    "iyy",
    "iyyPSF",
    "parentDiaSourceId",
    "pixelFlags",
    "pixelFlags_bad",
    "pixelFlags_cr",
    "pixelFlags_crCenter",
    "pixelFlags_edge",
    "pixelFlags_injected",
    "pixelFlags_injected_template",
    "pixelFlags_injected_templateCenter",
    "pixelFlags_injectedCenter",
    "pixelFlags_interpolated",
    "pixelFlags_interpolatedCenter",
    "pixelFlags_nodata",
    "pixelFlags_nodataCenter",
    "pixelFlags_offimage",
    "pixelFlags_saturated",
    "pixelFlags_saturatedCenter",
    "pixelFlags_streak",
    "pixelFlags_streakCenter",
    "pixelFlags_suspect",
    "pixelFlags_suspectCenter",
    "psfChi2",
    "psfFlux",
    "psfFlux_flag",
    "psfFlux_flag_edge",
    "psfFlux_flag_noGoodPixels",
    "psfFluxErr",
    "psfNdata",
    "ra",
    "ra_dec_Cov",
    "raErr",
    "reliability",
    "scienceFlux",
    "scienceFluxErr",
    "shape_flag",
    "shape_flag_no_pixels",
    "shape_flag_not_contained",
    "shape_flag_parent_source",
    "time_processed",
    "trail_flag_edge",
    "trailAngle",
    "trailDec",
    "trailFlux",
    "trailLength",
    "trailRa",
    "visit",
    "x",
    "xErr",
    "y",
    "yErr",
]


def make_ssobject(target, num_obs):
    return {
        "discoverySubmissionDate": 60000.0,
        "numObs": num_obs,
        "ssObjectId": to_id(target),
    }


def get_ephemeris(target, mjd):
    # use comet options for comets and comet-like orbits
    opts = {}
    if Names.asteroid_or_comet(target) or target.startswith("A/"):
        opts = dict(no_fragments=True, closest_apparition=True)

    helio = []
    topo = []
    for _mjd in np.array_split(mjd, len(mjd) // 50):
        q = Horizons(
            target,
            id_type=None if target.isdigit() else "designation",
            epochs=_mjd,
            location="@sun",
        )
        helio.append(q.vectors(**opts))

        q.location = "X05"
        topo.append(q.vectors(**opts))

    helio = vstack(helio)
    topo = vstack(topo)

    return helio, topo


def phase_angle(helio, topo):
    # sun->target
    rs = np.array(
        [
            helio["x"],
            helio["y"],
            helio["z"],
        ]
    )

    # observer->target
    ro = np.array(
        [
            topo["x"],
            topo["y"],
            topo["z"],
        ]
    )

    # phase angle
    phase = np.degrees(np.arccos(np.dot(rs, ro) / helio["range"] / topo["range"]))

    return phase


def make_source(target, ssobject, data):
    helio, topo = get_ephemeris(target, data["mjd"])
    sssource = []
    diasource = []
    for i, row in enumerate(data):
        diaSourceId = to_id(str(row))
        apflux = (row["mag"] * u.ABmag).to_value(u.nJy)
        diasource.append(
            {
                "apFlux": apflux,
                "apFluxErr": row["mag_err"] * apflux * 1.0857,
                "band": filter.get(row["filter"], row["filter"]),
                "midpointMjdTai": Time(row["mjd"], format="mjd", scale="utc").tai.mjd,
                "ssObjectId": ssobject["ssObjectId"],
                "diaSourceId": diaSourceId,
            }
        )

        for col in diasource_columns:
            diasource[-1][col] = np.nan

        sssource.append(
            {
                "diaSourceId": diaSourceId,
                "eclipticBeta": np.nan,
                "eclipticLambda": np.nan,
                "galacticB": np.nan,
                "galacticL": np.nan,
                "heliocentricDist": helio["range"][i],
                "heliocentricVX": helio["vx"][i],
                "heliocentricVY": helio["vy"][i],
                "heliocentricVZ": helio["vz"][i],
                "heliocentricX": helio["x"][i],
                "heliocentricY": helio["y"][i],
                "heliocentricZ": helio["z"][i],
                "phaseAngle": phase_angle(helio[i], topo[i]),
                "residualDec": np.nan,
                "residualRa": np.nan,
                "ssObjectId": ssobject["ssObjectId"],
                "topocentricDist": topo["range"][i],
                "topocentricVX": topo["vx"][i],
                "topocentricVY": topo["vy"][i],
                "topocentricVZ": topo["vz"][i],
                "topocentricX": topo["x"][i],
                "topocentricY": topo["y"][i],
                "topocentricZ": topo["z"][i],
            }
        )

    return diasource, sssource


def get_orbit(target, mjd):
    # use comet options for comets and comet-like orbits
    opts = {}
    if Names.asteroid_or_comet(target) or target.startswith("A/"):
        opts = dict(no_fragments=True, closest_apparition=True)

    q = Horizons(
        target,
        id_type=None if target.isdigit() else "designation",
        epochs=[np.mean(mjd)],
        location="@sun",
    )
    orbit = q.elements(**opts)
    return orbit


def make_orbit(target, data):
    orbit = get_orbit(target, data["mjd"])
    H = "M1" if "M1" in orbit.colnames else "H"
    return {
        "e": orbit["e"][0],
        "epoch": Time(orbit["datetime_jd"][0], format="jd").mjd,
        "incl": orbit["incl"][0],
        "mpcDesignation": target,
        "mpcH": orbit[H][0],
        "node": orbit["Omega"][0],
        "peri": orbit["w"][0],
        "q": orbit["q"][0],
        "ssObjectId": to_id(target),
        "t_p": Time(orbit["Tp_jd"][0], format="jd").mjd,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="path to the activity-for-brokers data/ directory")
    parser.add_argument("filename", help="save results to this database file name")
    args = parser.parse_args()

    if not os.path.exists(args.path):
        print(f"{args.path} does not exist")
        exit(1)

    db = sqlite3.connect(args.filename)

    diasource = []
    sssource = []
    ssobject = []
    mpcorb = []
    target_table = []
    for target, fn in targets.items():
        print(target)
        data = ascii.read(os.path.join(args.path, fn))

        # create the data tables
        ssobject.append(make_ssobject(target, len(data)))
        source = make_source(target, ssobject[-1], data)
        diasource.extend(source[0])
        sssource.extend(source[1])
        mpcorb.append(make_orbit(target, data))

        target_table.append(
            {"target": target, "ssObjectId": ssobject[-1]["ssObjectId"]}
        )

    pd.DataFrame(ssobject).to_sql("SSObject", con=db, if_exists="replace", index=False)
    pd.DataFrame(diasource).to_sql(
        "DiaSource", con=db, if_exists="replace", index=False
    )
    pd.DataFrame(sssource).to_sql("SSSource", con=db, if_exists="replace", index=False)
    pd.DataFrame(mpcorb).to_sql("MPCOrb", con=db, if_exists="replace", index=False)

    Table(target_table).pprint_all()


if __name__ == "__main__":
    main()

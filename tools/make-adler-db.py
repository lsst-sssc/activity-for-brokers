"""
Make an Adler database from real light curves.

See data/README.md for input data format.

Output format based on DP1:

SSSource columns
Index(['diaSourceId', 'eclipticBeta', 'eclipticLambda', 'galacticB',
       'galacticL', 'heliocentricDist', 'heliocentricVX', 'heliocentricVY',
       'heliocentricVZ', 'heliocentricX', 'heliocentricY', 'heliocentricZ',
       'phaseAngle', 'residualDec', 'residualRa', 'ssObjectId',
       'topocentricDist', 'topocentricVX', 'topocentricVY', 'topocentricVZ',
       'topocentricX', 'topocentricY', 'topocentricZ'],
      dtype='object')

SSObject columns
Index(['discoverySubmissionDate', 'numObs', 'ssObjectId'], dtype='object')

MPCOrb columns
Index(['e', 'epoch', 'incl', 'mpcDesignation', 'mpcH', 'node', 'peri', 'q',
       'ssObjectId', 't_p'],
      dtype='object')

"""

import os
import hashlib
import sqlite3
import argparse

import numpy as np
from astropy.io import ascii
from astropy.time import Time
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
    "7968": "7968_ZTF_FINK.csv",
    "2016 UU121": "2016UU121_ZTF_FINK.csv",
    "2060": "Chiron_ZTF_FINK.csv",
    "60558": "Echeclus_ZTF_FINK.csv",
    "6478": "Gault_ZTF.csv",
    # "C/Hy model": "LPC_y3.csv",
    # "C/Hab model": "LPC_0.3a-1b.csv",
}


def to_id(s):
    """Convert a string to stable ID (integer) value using an MD5 hash.

    The hash is trimmed to 8 digits.  Cross your fingers for no collisions!

    """

    hash = hashlib.md5()
    hash.update(s.encode())
    return int.from_bytes(hash.digest()) % 10**8


def make_object(target, num_obs):
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

    eph = None
    orbit = None
    for _mjd in np.array_split(mjd, len(mjd) // 100):
        q = Horizons(
            target,
            id_type="designation",
            epochs=_mjd,
            location="X05",
        )
        _eph = q.ephemerides(**opts)

        if eph is None:
            eph = _eph
        else:
            eph.add_rows(_eph)

    q.epochs = [np.mean(mjd)]
    orbit = q.elements(**opts)

    return eph, orbit


def make_sources(target, ssobject, data):
    eph, orbit = get_ephemeris(target, data["mjd"])
    breakpoint()
    rows = []
    for row in data:
        rows.append(
            {
                "diaSourceId": to_id(str(row)),
                "eclipticBeta": 0,
                "eclipticLambda": 0,
                "galacticB": 0,
                "galacticL": 0,
                "heliocentricDist": 0,
                "heliocentricVX": 0,
                "heliocentricVY": 0,
                "heliocentricVZ": 0,
                "heliocentricX": 0,
                "heliocentricY": 0,
                "heliocentricZ": 0,
                "phaseAngle": 0,
                "residualDec": 0,
                "residualRa": 0,
                "ssObjectId": ssobject["ssObjectId"],
                "topocentricDist": 0,
                "topocentricVX": 0,
                "topocentricVY": 0,
                "topocentricVZ": 0,
                "topocentricX": 0,
                "topocentricY": 0,
                "topocentricZ": 0,
            }
        )

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="path to the activity-for-brokers data/ directory")
    parser.add_argument("filename", help="save results to this database file name")
    args = parser.parse_args()

    if not os.path.exists(args.path):
        print(f"{args.path} does not exist")
        exit(1)

    db = sqlite3.connect(args.filename)

    sssources = []
    ssobjects = []
    for target, fn in targets.items():
        data = ascii.read(os.path.join(args.path, fn))

        # create the data tables
        ssobject = make_object(target, len(data))
        sssource = make_sources(target, ssobject, data)

        # concatenate the data
        ssobjects.append(ssobject)
        sssources.extend(sssource)

        break

    # ssobject.to_sql("SSObject", con=db, if_exists="replace", index=False)
    # sssource.to_sql("SSSource", con=db, if_exists="replace", index=False)

    return data, ssobject, sssource


if __name__ == "__main__":
    data, ssobject, sssource = main()

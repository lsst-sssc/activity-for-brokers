import requests
import pandas as pd
import io
import datetime

def jd_to_mjd(jd):
    # Convert Julian Date (JD) to Modified Julian Date (MJD)
    mjd = jd - 2400000.5
    return mjd

from astroquery.jplhorizons import Horizons
import pandas as pd

def query_horizons(object_name, mjds, site_code="I41"):
    """
    Queries JPL Horizons for a solar system object at the specified MJDs and observatory site.

    Parameters:
    object_name (str): The name of the solar system object to query.
    mjds (list): A list of MJDs to query.
    site_code (str): The observatory site code (defaults to 'I41').

    Returns:
    pd.DataFrame: A DataFrame containing the Horizons query results.

    TODO: choose nearest epoch JPL ID automatically instead of failing. 9/27/2024 COC
    """
    results = []
    
    for mjd in mjds:
        # Convert MJD to JD
        jd = mjd + 2400000.5
        
        # Query Horizons
        obj = Horizons(id=object_name, location=site_code, epochs=[jd])
        eph = obj.ephemerides()
        
        # Convert the astropy table to a pandas DataFrame and append the result
        results.append(eph.to_pandas())
        
    # Concatenate all results into a single DataFrame
    df = pd.concat(results, ignore_index=True)
    
    return df

def get_fink_data(objname, verbose=True):
    """
    Query the FINK broker for a solar system object.
    """
    if verbose: print(f'Starting FINK + JPL queries for {objname} now...')
    now = datetime.datetime.now()
    now = now.strftime("%Y-%m-%d %H:%M:%S")
#   print(f'It is now: {now}')
    # get data for 2016 UU121
    r = requests.post(
        'https://api.fink-portal.org/api/v1/sso',
        json={
            'n_or_d': '2016 UU121',
            'withEphem': True,
            'output-format': 'json'
        }
    )
    
    # Format output in a DataFrame
    pdf = pd.read_json(io.BytesIO(r.content))
    if len(pdf) == 0:
        raise KeyError(f'No such object or no data found at broker.')
    
#   print(pdf.columns)
    #
    # need MJD, so add a column via our function
    pdf['mjd'] = pdf['i:jd'].apply(jd_to_mjd)
    #
    # aperture kludge 9/27/2024 COC
    pdf['aper_arcsec'] = [0]*len(pdf) # 14 pixels, 1"/pixel; 14 is in the fink docs 9/27/2024 COC; 0 for now until these new fields make it into FINK-broker 9/27/2024 COC
    #
    # psf
    filter_dict = {1:'g', 2:'r'}
    pdf['filter'] = [filter_dict[fid] for fid in pdf['i:fid']]
    #
    jpl_df = query_horizons(object_name=objname, mjds=pdf['mjd'], site_code="I41")
#   print(jpl_df.columns)
    #
    # make new dataframe
    ndf = pd.DataFrame.from_dict({'mjd':pdf['mjd']})
    ndf['rh'] = jpl_df['r']
    ndf['delta'] = jpl_df['delta']
    ndf['alpha'] = jpl_df['alpha']
    ndf['aper_arcsec'] = pdf['aper_arcsec']
    ndf['mag'] = pdf['i:magpsf']
    ndf['mag_err'] = pdf['i:sigmapsf']
    ndf['query_datetime'] = [now]*len(pdf)
    objname_clean = objname.replace('/','').replace(' ','')
    outfile = f'{objname_clean}_ZTF_FINK.csv'
    ndf.to_csv(outfile)
    if verbose: print(f'Wrote {outfile} to disk.')


if __name__ == '__main__':
    import argparse # This is to enable command line arguments.
    parser = argparse.ArgumentParser(description='Standalone FINK-broker query tool for e.g., comet analysis. By Colin Orion Chandler (COC) (9/27/2024)')
    parser.add_argument('objects', metavar='O', type=str, nargs='+', help='object name')
    args = parser.parse_args()
    for objname in args.objects:
        get_fink_data(objname=objname)

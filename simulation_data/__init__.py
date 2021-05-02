import requests

#insert personal API key to download simulation data using get()
API = "ThisIsMyAPIKeyForIllustrisTNG"
#defined get()
def get(path, params=None):
    '''
    The oroginal version of this function can be found at
    https://www.tng-project.org/data/docs/api/
    '''
    # make HTTP GET request to path
    headers = {"api-key":API}
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r


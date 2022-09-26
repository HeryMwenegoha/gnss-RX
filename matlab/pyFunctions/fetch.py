# Fetch function
import os
import requests
import sys

import gzip
import shutil

# Note: this code example requires a .netrc file set up in the home directory:
# ~/.netrc
# User configuration settings using the .netrc file
# machine urs.earthdata.nasa.gov login <username> password <pword>

# Get a list of available files for the user config
def file_list(YYYY, DOY, HH):
    #python filelist.py https://cddis.nasa.gov/archive/gnss/data/hourly/2022/240/22/
    #Format for generic call
    #python filelist.py https://cddis.nasa.gov/archive/gnss/data/hourly/YYYY/DOY/HH/

    #url = sys.argv[1]
    url = "https://cddis.nasa.gov/archive/gnss/data/hourly/"+str(YYYY)+"/"+str(DOY)+"/"+str(HH)+"/"

    #Adds '*?list' to the end of URL if not included already
    if not url.endswith("*?list"):
        url = url + "*?list"

    #Makes request of URL, stores response in variable r
    r = requests.get(url)

    #Prints the results of the directory listing
   
    return r.text
    #for item in r.text:
    #    print(item)
    #print(r.text)

# Check that we have the file for the requested IGS station
def checkFile(YYYY, DOY, HH, filesReturned, IGSStationInitials):
    # Available flag
    isAvailable = False
    urlToFetch = ""

    # Empty string
    if filesReturned == "":
        return isAvailable,urlToFetch

    # Process list of files
    fileList = str.split(filesReturned, '\n')

    # The file to check
    #IGSStationInitials="ABPO00MDG_R_"
    checkFile = IGSStationInitials + str(YYYY) + str(DOY) + str(HH) + "00_01H_GN.rnx.gz"

    # Run through the returned list
    for item in fileList:
        filename = str.split(item)[0]

        if checkFile == filename:
            print("Rquested file is found: ")
            print("------------------------")
            print(filename)
            isAvailable = True
            fileToFetch = checkFile
            urlToFetch = "https://cddis.nasa.gov/archive/gnss/data/hourly/"+str(YYYY)+"/"+str(DOY)+"/"+str(HH)+"/"+checkFile
            break

    return isAvailable,urlToFetch

# Download the file
def fetch_files(url):
    # Format
    # python filelist.py https://cddis.nasa.gov/archive/gnss/data/hourly/YYYY/DOY/HH/
    #url = "https://cddis.nasa.gov/archive/gnss/data/hourly/2022/240/22/ABPO00MDG_R_20222402200_01H_GN.rnx.gz"
    #url = sys.argv[1]
    # Print CDDIS url
    print("------fetch_files------")
    print("fetchingUrl: "+url)

    # Assigns the local file name to the last part of the URL
    filename = "../CDDIS"+'/'+url.split('/')[-1]

    # Print the filename we are saving contents to.
    print("storageFilename: "+ filename)

    # Makes request of URL, stores response in variable r
    print("send http request")
    r = requests.get(url)

    # Opens a local file of same name as remote file for writing to
    print("write to file")
    with open(filename, 'wb') as fd:
        for chunk in r.iter_content(chunk_size=1000):
            fd.write(chunk)

    # Closes local file
    fd.close()
    print("---end-fetch_files---")

    # Return a filename
    return filename

def unzip_fetched_file(filename):
    print("----unzip_fetched_file----")
    newfilename = str.split(filename, ".gz")[0]
    with gzip.open(filename, 'rb') as f_in:
        with open(newfilename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


if __name__ == "__main__":
    # Reads the URL from the command line argument
    # YYYY = sys.argv[1]
    YYYY = 2022
    DOY = 240
    HH = 22

    IGSStationInitials="ABPO00MDG_R_"

    returnedFiles = file_list(YYYY, DOY, HH)

    isValid,fileTofetch = checkFile(YYYY, DOY, HH, returnedFiles, IGSStationInitials)

    if isValid:
        filename = fetch_files(fileTofetch)
        unzip_fetched_file(filename)
    else:
        print("Fetching files from CDDIS failed")
import requests
import os
from bs4 import BeautifulSoup

url ='https://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Eurasia/' 
r  = requests.get(url)
data = r.text
soup = BeautifulSoup(data)

links = []

for link in soup.find_all('a'):
  hr = link.get('href')
  if hr.endswith('.hgt.zip'):
      links.append((hr,url+hr))

for (name,link) in links:
  if os.path.isfile('dane/'+name):
    print("File {} already exists".format(name))
  else:
    print("Downloading file {}".format(link))
    r = requests.get(link, stream=True)
    with open(name, 'wb') as f:
      for chunk in r.iter_content(chunk_size=1024*1024):
        if chunk:
          f.write(chunk)


print("All files downloaded!")


# pysam - PYthon Satellite / Space debris AstroMetry

## Description
This script uses optical images of observations of satellites or space debris that revolve around Earth and calculates their position and orbital parameters. <br>
In aditional to the script, there also is a example config file to setup the program. The astrometry part of the script is done by locally installed [Astrometry.net](https://astrometry.net/) software. Orbit calculation part is done by [Find\_Orb](https://projectpluto.com/find_orb.htm) which is optional. The astrometry observations are saved as a .txt file, that follows [Minor Planet Center's 80 column format](https://www.minorplanetcenter.net/iau/info/OpticalObs.html).
### Config description
In the config file, user can set couple of parameters, that may be different depending on the obserwatory such as FITS header keywords regarding image dementions, observed object designation, observation date, output file name containing astrometry and Find\_Orb software path. If you do not have Find\_Orb software installed, then in config set the 'find\_orb\_installed' flag to false, and path to Find\_Orb is ignored.

## Usage:
```bash
python3 ./example_config.cfg "file_pattern"
```
> [!WARNING]
> The file pattern follows GNU\Linux wildcards convention and must be in single or double quotes, otherwise script won't work properly! <br>
> The image file name must not contain multiple dot characters!

## Requirements:
> [!WARNING]
> This script will NOT work on Windows!

The following script was tested on Fedora Linux 42 using Python version 3.12.11 and uses following packages:
- numpy 2.2.0
- matplotlib 3.10.3
- astropy 7.0.0
- photutils 2.2.0
- pandas 2.2.3

Aditionally the script uses locally installed [Astrometry.net](https://astrometry.net/) version 0.97 and [Find\_Orb](https://projectpluto.com/find_orb.htm) software, which is optional. 


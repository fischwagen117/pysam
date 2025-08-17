# image fits files cannot contain a dot "." character in the name other than the file extention!
def print_error(err_msg, exit_code=1):
    print(f"\033[31mError\033[0m: {err_msg}")
    exit(exit_code)

def log(log_path, log_str):
    print(log_str)
    dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_path, "a") as log:
        log.write(f"{dt}: {log_str}\n")

def parse_obs_date(header, time_zone, obs_time_keyword, exp_time_keyword):
    """
    Converts observation timestamp to the time format (YYYY MM DD.dddddd) required by Minor Planet Center, used in orbit calculation software

    Parameters:
        timestamp: timestamp of the obserwation, example: 2025-03-26T08:36:42.080 (taken from the fits header)
        time_zone: amount of hours in the timezone to create UTC time, example time_zone = -3 for place with timezone of UTC-03:00
        obs_time_keyword, exp_time_keyword: keywords from the fits header for observation timestamp and exposure time
    Returns:
        string - observation timestamp adjusted by timezone and exposure time
    """
    exp_time = int(header[exp_time_keyword] / 2)
    timestamp = header[obs_time_keyword]
    #DATE-LOC - 2025-03-26t08:36:42.080
    #           2025-03-26t08:51:17.960
    tz = float(time_zone) * 3600 * (-1)
    dt = datetime.fromisoformat(timestamp) + timedelta(seconds = tz + exp_time) # adjustments for timezone and exposure time
    #dt = datetime.fromisoformat(timestamp) + timedelta(seconds = tz) # adjustments only for timezone
    year, month, day, hour, minute, second = dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second
    hour_int = int(hour)
    hour_decimal = (hour_int + minute / 60 + second / 3600) / 24
    hr = int(round(hour_decimal, 6) * 1000000)
    return f"{year:04} {month:02} {day:02}.{hr:06d}"

def hms2dec(hms, signed=False):
    """
    Converts from hours:minutes:seconds format to decimal hours format. Can also be used to convert degress to decimal degerss.
    """
    negative = False
    if hms < 0:
        negative = True
        hms = hms * (-1)
    h_int = int(hms)
    min = (hms - h_int) * 60
    min_int = int(min)
    sec = (min - min_int) * 60
    sec_int = int(sec)
    sec_decimal = int(round(sec - sec_int, 2) * 100)
    if signed:
        if negative:
            return f"-{h_int:02} {min_int:02} {sec_int:02}.{sec_decimal:02}"
        return f"+{h_int:02} {min_int:02} {sec_int:02}.{sec_decimal:02}"
    return f"{h_int:02} {min_int:02} {sec_int:02}.{sec_decimal:02}"

def create_space_string(space_count):
    """
    Creates string with set amount of space characters.
    """
    out = ""
    for _ in range(space_count):
        out += " "
    return out

def fix_obj_name(obj_name):
    """
    Changes object name such that it's length is 8 characters
    """
    obj_set_length = 7
    obj_actual_length = len(obj_name)
    if obj_actual_length > obj_set_length:
        return obj_name[0:obj_set_length]
    else:
        for i in range(obj_set_length - obj_actual_length):
            obj_name += " "
        return obj_name

def find_orbit(path_to_fo, output_file, output_path):
    """
    Function runs non-interactive version of find_orb software to find orbital parameters of an observed object
    The -h 3 parameter is to set the orbit to be geocentric
    """
    cmd = f"./fo {output_path}/{output_file} -O {output_path} -h 3"
    print("Running command: ", cmd)
    subprocess.run(cmd.split())
    return

def read_tle(output_path, tle_out):
    """
    Read the elements.txt file from find_orb output to get the TLE
    """
    with open(output_path, "r") as obs_out:
        with open(tle_out, "w") as tle:
            for i in range(10):
                obs_out.readline()

            tle1 = obs_out.readline()
            tle2 = obs_out.readline()
            tle_str = f"{tle1}{tle2}"
            print(f"TLE:\n{tle_str}")
            tle.write(tle_str)

class Segmentation:
    def __init__(self, file_path, config):
        im_width_header_key = config["KEYWORDS"]["im_width_key"]
        im_height_header_key = config["KEYWORDS"]["im_height_key"]
        obsdate_header_key = config["KEYWORDS"]["observation_time_key"]
        exposure_time_key = config["KEYWORDS"]["exposure_time_key"]
        object_key = config["KEYWORDS"]["object_key"] 
        center_ra_key = config["KEYWORDS"]["ra_key"]
        center_dec_key = config["KEYWORDS"]["dec_key"]
        time_zone = config["OTHER"]["time_zone"]

        self.file_path = file_path
        self.data = fits.getdata(file_path)
        self.file_name = self.file_path.split("/")[-1].split(".")[0]
        self.header = fits.getheader(file_path)
        self.im_width = self.header[im_width_header_key]
        self.im_height = self.header[im_height_header_key]
        self.obs_date = parse_obs_date(self.header, time_zone, obsdate_header_key, exposure_time_key)
        self.object = self.header[object_key]
        self.ra = self.header[center_ra_key]
        self.dec = self.header[center_dec_key]
        self.radius = config["OTHER"]["radius"]

    def detect_background(self, data):
        """
        Dectects background and calculates threshold needed to find objects.
        """
        BKG_RECT = (50, 50)
        THRESHOLD_SCALE = 2.0
        bkg_est = MedianBackground()
        bkg = Background2D(data, BKG_RECT, filter_size=(3, 3), bkg_estimator=bkg_est)
        threshold = THRESHOLD_SCALE * bkg.background
        return threshold

    def convolve_data(self, data):
        """
        Convolves data with a 2D gaussian kernel to detect short streaks.
        """
        FWHM = 1
        KERNEL_SIZE = 5
        kernel = make_2dgaussian_kernel(FWHM, size=KERNEL_SIZE)
        convolved_data = convolve(data, kernel)
        return convolved_data

    # returns segment_map
    def detect_sources(self, convolved_data, threshold):
        """
        Function detects image features (strikes, stars) based on the background threshold and the convolution kernel.
        """
        return detect_sources(convolved_data, threshold, npixels=10, connectivity=4)
    
    # returns segm_deblend
    def deblend_sources(self, convolved_data, segment_map):
        """
        Function tries to separate overlapping sources based on convolved data and detected features. 
        """
        CONTRAST = 0.5
        NLEVELS = 32
        return deblend_sources(convolved_data, segment_map, npixels=10, nlevels=NLEVELS, contrast=CONTRAST, progress_bar=False)

    def create_sources_table(self, data, segm_deblend, convolved_data):
        """
        Creates pandas table for easier data manipulation.
        """
        cat = SourceCatalog(data, segm_deblend, convolved_data=convolved_data)
        return cat.to_table().to_pandas()

    def detect_satellite(self, table):
        """
        Function detects satellite on the image based on its eccentricity value and mean eccentricity for the rest of the objects.It also rejects objects that are to close to the image edges to prevent false detections.
        """
        x_offset = 50
        y_offset = 50
        table = table[table["xcentroid"] > x_offset]
        table = table[table["xcentroid"] < self.im_width - x_offset]
        table = table[table["ycentroid"] > y_offset]
        table = table[table["ycentroid"] < self.im_height - y_offset]
        mean_eccentricity = table.mean()["eccentricity"]
        sat = table.sort_values(by=["eccentricity"])
        if mean_eccentricity > 0.9:
            sat = sat.iloc[0]
        else:
            sat = sat.iloc[-1]
        return sat
    def save_star_positions(self, table, sat, out_file, star_count=60):
        """
        Function saves brightest star positions in binary fits file used by astrometry solver.
        """
        LABEL_COL = 0
        XCENTROID_COL = 1
        YCENTROID_COL = 2
        row_count = len(table)
        sat_label = sat["label"]
        # for star list use kron_flux or segment_flux
        star_list = table.sort_values(by=["segment_flux"], ascending=False).head(star_count)
    
        if star_count > row_count:
            star_count = row_count
        data = []
        for i in range(star_count):
            current_label = star_list.iloc[i, LABEL_COL]
            if sat_label != current_label:
                x_pos = star_list.iloc[i, XCENTROID_COL]
                y_pos = star_list.iloc[i, YCENTROID_COL]
                data.append((x_pos, y_pos))

        data = np.rec.array(data, formats="float32,float32", names="XIMAGE,YIMAGE")
        fits.writeto(out_file, data, overwrite=True)

    def clean_table(self, tbl):
        """
        Function removes small repeated apertures based on their distance and angle to each other.
        """
        mag = lambda x, y: np.sqrt(np.pow(x, 2) + np.pow(y, 2))
        angle = lambda x, y: np.atan2(y, x) * 180 / np.pi
        ang = []
        dist = {}
        size = tbl.shape[0]
        bxmax = tbl["bbox_xmax"]
        bxmin = tbl["bbox_xmin"]
        bymax = tbl["bbox_ymax"]
        bymin = tbl["bbox_ymin"]
        centroid_length = np.sqrt((bxmax-bxmin)**2 + (bymax-bymin)**2)
        min_dist = np.median(centroid_length)
        for i in range(size):
            dist2 = {}
            x0, y0 = tbl["xcentroid"][i], tbl["ycentroid"][i]
            for j in range(i, size):
                if i != j:
                    x, y = tbl["xcentroid"][j], tbl["ycentroid"][j]
                    r = mag(x-x0, y-y0)
                    if r < min_dist + 5:
                        phi = angle(x-x0, y-y0)
                        ang.append(phi)
                        dist2[j+1] = phi
            dist[i+1] = dist2

        ang = np.array(ang)
        ang_med = np.median(ang)
        to_remove = []
        ang_offset = 10
        high_ang_val = ang_med + ang_offset
        low_ang_val= ang_med - ang_offset
        for i, labels in dist.items():
            for j, phi in labels.items():
                if phi < high_ang_val and phi > low_ang_val:
                    if i not in to_remove:
                        to_remove.append(i)
                    elif j not in to_remove:
                        to_remove.append(j)

        to_remove = np.array(to_remove)
        to_remove -= 1 # conversion from labels to indices
        return tbl.drop(index=to_remove)

    def solve_xy(self, xy_pos_file, im_width, im_height, ra, dec, radius):
        """
        Funciton calls astrometry bash command with required arguments to solve image and get the WCS (World Coordinate System) fits file.
        """
        cmd = f"solve-field --no-plots -w {im_width}, -e {im_height} --x-column XIMAGE --y-column YIMAGE --ra {ra} --dec {dec} --radius {radius} {xy_pos_file}"
        print("Running command: ", cmd)
        subprocess.run(cmd.split())
        return

    def plot_figs(self, data, sat):
        plt.imshow(data, cmap="gray", vmax=2000)
        plt.scatter([sat["xcentroid"]], [sat["ycentroid"]], marker=",", c="red")
        plt.axis([2300, 2500, 1400, 1500])
        plt.savefig(f"{self.file_name}.png", format="png")
        plt.close()

def process_images_and_find_orbit(file_path, config, output):
    """
    The main loop that does image processing and astrometry. Results are written to the output file which is in the MPC 80-column format ready to be used to determine orbit using find_orb software or other orbit determination software. 

    Parameters:
        file_path: pattern  
    """
    astrometry_dir = "./xy_files"
    output_dir = "output"
    current_dir = os.getcwd()
    log_file_name = datetime.now().strftime("LOG_%Y_%m_%d_%H_%M_%S.txt")
    if os.path.exists(astrometry_dir) == False:
        os.makedirs(astrometry_dir)

    if os.path.exists(output_dir) == False:
        os.makedirs(output_dir)

    f = open("test.txt", "w")

    images = glob(file_path)
    image_count = len(images)
    mean_eccentricity = 0
    with open(f"{output_dir}/{output}", "w") as f_out:
        for im_number, image in enumerate(images):
            segm = Segmentation(image, config)

            log(log_file_name, f"Processing image [{im_number+1}/{image_count}]: {segm.file_name}")
            im_data = segm.data
            threshold = segm.detect_background(im_data)
            log(log_file_name, "Background Detected")
            convolved_data = segm.convolve_data(im_data)
            segment_map = segm.detect_sources(convolved_data, threshold)
            log(log_file_name, "Segment map created")
            sources_deblended = segm.deblend_sources(convolved_data, segment_map)
            tbl = segm.create_sources_table(im_data, sources_deblended, convolved_data)
            mean_eccentricity = tbl.mean()["eccentricity"]
            log(log_file_name, mean_eccentricity)
            if mean_eccentricity > 0.5:
                log(log_file_name, "Cleaning Table")
                tbl = segm.clean_table(tbl)
            sat = segm.detect_satellite(tbl)
            #print(sat["xcentroid"], sat["ycentroid"])
            stars_pos_file = f"./xy_files/{segm.file_name}.xy"
            log(log_file_name, f"Saved star positions for {segm.file_name}")
            segm.save_star_positions(tbl, sat, stars_pos_file)
            wcs_file = f"./xy_files/{segm.file_name}.wcs"
            segm.solve_xy(stars_pos_file, segm.im_width, segm.im_height, segm.ra, segm.dec, segm.radius)
            wcs = WCS(fits.getheader(wcs_file))
            sat_x = sat["xcentroid"]
            sat_y = sat["ycentroid"]
            sat_pos = wcs.pixel_to_world(sat_x, sat_y)
            sat_ra = sat_pos.ra.value / 15 # conversion from degress to hours
            sat_dec = sat_pos.dec.value
            log(log_file_name, f"Satellite position: x:{sat_x} y:{sat_y}\nSatellite coordinates: ra:{hms2dec(sat_ra)} dec:{hms2dec(sat_dec)}")
            f.write(f"{segm.file_name}.wcs, satx, saty: {sat_x}, {sat_y}, sat_pos: {sat_ra}, {sat_dec}, date: {segm.obs_date}\n")
            out_string = "" + create_space_string(5) + fix_obj_name(segm.object) + "  B" + segm.obs_date + f"{hms2dec(sat_ra)} " + hms2dec(sat_dec, signed=True) + create_space_string(21) + "X07\n"
            f_out.write(out_string)
            #segm.plot_figs(segm.data, sat)
            log(log_file_name, "Data written\n")

    f.close()
    find_orb_installed = config["FIND_ORB"]["find_orb_installed"]
    if find_orb_installed:
        find_orb_dir = config["FIND_ORB"]["path_to_find_orb"]
        log(log_file_name, "Finding orbit...")
        os.chdir(find_orb_dir)
        full_output_path = f"{current_dir}/{output_dir}"
        find_orbit(find_orb_dir, output, full_output_path)
        elements_txt_path = f"{output_dir}/elements.txt"
        tle_path = f"{output_dir}/{segm.object}_tle.txt"
        os.chdir(current_dir)
        read_tle(elements_txt_path, tle_path)
        log(log_file_name, f"Done. Results are at {output_dir} directory.")
    else:
        log(log_file_name, "find_orb software is not installed, however you can still see the object observation results in the {output_dir} directory and / or try another orbit determination software.")

if __name__ == "__main__":
    from sys import argv
    import configparser
    if len(argv) < 3:
        print_error('''Not enough arguments were provided!
                    Usage: python3 segmentation.py CONFIG_PATH FILES_PATTERN
                    Example: python3 segmentation.py ./config.cfg "./fits/*4,00*"
                    Note: pattern MUST be in quotes!''')
    config_path, pattern = argv[1], argv[2]
    config = configparser.ConfigParser()
    config.read(config_path)
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import pandas as pd
    from photutils.background import Background2D, MedianBackground
    from astropy.convolution import convolve
    from photutils.segmentation import make_2dgaussian_kernel, detect_sources, deblend_sources, SourceCatalog
    from glob import glob
    from astropy.wcs import WCS
    import subprocess
    from datetime import datetime, timedelta
    import os

    plt.style.use("dark_background")
    plt.rcParams["figure.figsize"] = [14, 10]
    process_images_and_find_orbit(pattern, config, config["OTHER"]["output_file"])

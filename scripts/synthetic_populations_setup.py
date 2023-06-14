import csv
import subprocess
import sys
import os

def download_files(csv_file, files_to_download, FRED_pop_path):
    downloaded_files = []
    skipped_files = []
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            file_name, url = row
            if file_name in files_to_download:
                file_name = f'{FRED_pop_path}{file_name}'
                if not os.path.exists(file_name):
                    subprocess.run(["curl", "-L", "-o", file_name + ".zip", url])
                    downloaded_files.append(file_name)
                else:
                    skipped_files.append(file_name)
    return downloaded_files, skipped_files

if __name__ == "__main__":
    csv_file = 'input_files/population_download_info.csv'
    files_to_download = sys.argv[1].split()
    FRED_pop_path = sys.argv[2]
    downloaded_files, skipped_files = download_files(csv_file, files_to_download, FRED_pop_path)
    print(' '.join(downloaded_files), '|', ' '.join(skipped_files))

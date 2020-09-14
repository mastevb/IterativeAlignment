import sys
import os

directory = "/Users/mabochen/Desktop/Research/IALR/Linear_Regression_on_Physiochemical_Props/filtered_data"

for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        print("processing: " + filename)
        os.system("python3 main.py " + filename)

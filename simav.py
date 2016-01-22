"""
Provides code to run simav2csv function from summarize.py

"""
import glob, os
from summarize01 import simav2csv
from concurrency_figures import DegreeTab


if __name__ == '__main__':
    for eps in range(11):
        dirpath = "/Users/Jelena/AU/sawers/03nc/duration03nc/"
        fpattern = 'eps{0:02d}sim??.out'.format(eps) 
        outname = 'eps{0:02d}simav.csv'.format(eps) 
        simav2csv(dirpath, fpattern, outname, force=False) 

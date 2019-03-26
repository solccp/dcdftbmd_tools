#!/usr/bin/env python3

import _dcdftbmd_tools_cpp_extension_md as ext_md
import pathlib
import sys
import argparse

def merge_md_outputs(folders, merged_prefix='merged_', 
    trajectory_filename='traject', mulliken_filename='mulliken', velocity_filename='velocity', 
    merge_trajectory=True, merge_mulliken=False, merge_velocity=False, verbose=False):
    
    last_steps = ext_md.merge_main_output(folders, f'{merged_prefix}dftb.out', verbose)
    print(last_steps)
    if merge_trajectory:
        #oid merge_datafile(const std::vector<std::string> &folders, const vector<int>& last_step_nos, 
        # const std::string& filename, const std::string& merged_filename, bool verbose)
        merged_filename = f'{merged_prefix}{pathlib.Path(trajectory_filename).name}'
        ext_md.merge_datafile(folders, last_steps, trajectory_filename, merged_filename, verbose)

def run():
    parser = argparse.ArgumentParser(description='MD Run Merger')
    parser.add_argument('folders', metavar='folders', type=str, nargs='+')
    opt = parser.parse_args(sys.argv[1:])
    merge_md_outputs(opt.folders)
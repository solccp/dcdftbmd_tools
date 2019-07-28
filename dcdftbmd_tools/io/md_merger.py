#!/usr/bin/env python3

import _dcdftbmd_tools_cpp_extension_md as ext_md
import pathlib
import sys
import argparse





def merge_md_outputs(opt):
    # print(opt)
    filename = pathlib.Path(opt.main_filename).name
    
    write_merged = False
    write_parsed_data = False
    
    if (opt.type == "main" or opt.type == "all"):
        write_merged = True
        if (opt.parse_main_data):
            write_parsed_data = True

    last_steps = ext_md.merge_main_output(opt.folders, opt.folder_as_suffix, filename, opt.verbose, opt.stride_size, 
        f'{opt.merge_prefix}{filename}', write_merged, 
        f'{opt.merge_prefix}dftb_data.dat', write_parsed_data)

    if (opt.type == "traj" or opt.type == "all"):
        filename = pathlib.Path(opt.traj_filename).name
        merged_filename = f'{opt.merge_prefix}{pathlib.Path(filename).name}'
        ext_md.merge_datafile(opt.folders, opt.folder_as_suffix, last_steps, filename, merged_filename, opt.stride_size, opt.verbose)

    if (opt.type == "vel" or opt.type == "all"):
        filename = pathlib.Path(opt.vel_filename).name
        merged_filename = f'{opt.merge_prefix}{pathlib.Path(filename).name}'
        ext_md.merge_datafile(opt.folders, opt.folder_as_suffix, last_steps, filename, merged_filename, opt.stride_size, opt.verbose)

    if (opt.type == "mull" or opt.type == "all"):
        filename = pathlib.Path(opt.mull_filename).name
        merged_filename = f'{opt.merge_prefix}{pathlib.Path(filename).name}'
        if (opt.mull_nac):
            merged_filename = f'{opt.merge_prefix}nac'
        ext_md.merge_mullfile(opt.folders, opt.folder_as_suffix, last_steps, filename, merged_filename, opt.stride_size, opt.mull_nac, opt.verbose)


def run():
    parser = argparse.ArgumentParser(description='MD Run Merger')
    parser.add_argument('folders', metavar='folders', type=str, nargs='+')
    
    parser.add_argument('-p', '--merge_prefix', type=str, default='merged_')
    parser.add_argument('-s', '--stride_size', type=int, default=1)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-f', '--folder_as_suffix', action='store_true')

    parser.add_argument('-t', '--type', choices=['all', 'main', 'traj', 'mull', 'vel'], default='main', required=True)
    
    
    parser.add_argument('--main_filename', type=str, default='dftb.out')
    parser.add_argument('--parse_main_data', action='store_true')

    parser.add_argument('--traj_filename', type=str, default='traject')
    parser.add_argument('--vel_filename', type=str, default='velocity')
    parser.add_argument('--mull_nac', action='store_true')
    parser.add_argument('--mull_filename', type=str, default='mulliken')
    
    opt = parser.parse_args(sys.argv[1:])
    merge_md_outputs(opt)
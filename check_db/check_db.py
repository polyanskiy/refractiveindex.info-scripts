#!/usr/bin/python
# coding: utf-8

import argparse
import yaml

import os
import fnmatch
from collections import Counter

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="""
                                     Checks that all files on disk also appear
                                     in the index library.yml.
                                     Also checks that all files in the index
                                     have their counterpart on disk.
                                     """)
    parser.add_argument('database', action='store', 
                                     help='database root directory (database/)')
    args = parser.parse_args()
    
    
    db_yml = "library.yml"
    
    db_path = args.database
    lib_path = os.path.join(db_path, db_yml)

    ## List all YML files to process, recursively. Omit library.yml.
    yaml_files = [os.path.join(dirpath, f)
                  for dirpath, dirnames, files in os.walk(db_path)
                  for f in fnmatch.filter(files, '*.yml') if f != 'library.yml'
                 ]
    
    ## Load index
    data = yaml.load(open(lib_path, 'r').read())
    ct = []
    indexed_files = []
    
    # 5 main categories, ordered: main, organic, glasses, other, 3D
    for cat in data:
        cat_shelf = cat["SHELF"]
        cat_name = cat["name"]
        cat_content = cat["content"]
        ct.append({"name":cat_shelf, "desc":cat_name, "content":{}})
        # Each category has several books
        divider = "root"
        ct[-1]["content"][divider] = []
        for book in cat_content:
            if "DIVIDER" in book:
                divider = book["DIVIDER"]
                ct[-1]["content"][divider] = []
            elif "BOOK" in book:
                book_cat = book["BOOK"]
                book_name = book["name"]
                ct[-1]["content"][divider].append({"book_cat":book_cat,
                                    "book_name":book_name,
                                    "book_page":{}})
                subpage = "root"
                ct[-1]["content"][divider][-1]["book_page"][subpage] = []
                for page in book["content"]:
                    if "DIVIDER" in page:
                        subpage = page["DIVIDER"]
                        ct[-1]["content"][divider][-1]["book_page"][subpage]= []
                    else:
                        page_auth = page["PAGE"]
                        page_name = page["name"]
                        page_path = os.path.normpath(os.path.join(db_path, page["path"]))
                        ct[-1]["content"][divider][-1]["book_page"][subpage] \
                            .append({
                                "page_auth":page_auth,
                                "page_name":page_name,
                                "page_path":page_path
                                })
                        indexed_files.append(page_path)
        
    yml_files_num = len(yaml_files)
    yml_files_indexed_num = len(indexed_files)
    if yml_files_indexed_num != yml_files_num:
        print("File number mismatch: {} files in the index, {} files on disk" \
            .format(yml_files_indexed_num, yml_files_num)) 
    unique_files = set(yaml_files)
    unique_files_index = set(indexed_files)
    
    print("{} unique files in the index, {} unique files on disk" \
        .format(len(unique_files_index), len(unique_files)))
    print("")
    counts = Counter(indexed_files)
    most_referenced = counts.most_common(yml_files_indexed_num - yml_files_num)
    for name, value in most_referenced:
        if value <= 1: continue
        print("{} appears {} times in the index".format(name, value))
    
    print("")
    files_intersection = unique_files.symmetric_difference(unique_files_index)
    print("Files in one set but not the other: ")
    for ff in files_intersection:
        print(ff)
    
    print("")
    if (unique_files_index <= unique_files):
        # Files on disk are a subset of the files in the index
        print("All files of the index have their counterpart on disk.")
        
        diff = unique_files.difference(unique_files_index)
        if len(diff) > 0:
            print("Some files on disk are not in the index :")
            for ff in diff:
                print(ff)
    else:
        print("Some files of the index are not on disk :")
        for ff in unique_files_index.difference(unique_files):
            print(ff)
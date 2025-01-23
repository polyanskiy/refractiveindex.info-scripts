#!/usr/bin/python
# coding: utf-8

import argparse
import yaml

import os
import fnmatch

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="""
                                     Checks that all files on disk also appear
                                     in the index catalog-nk.yml.
                                     Also checks that all files in the index
                                     have their counterpart on disk.
                                     """)
    parser.add_argument('database', action='store', 
                                     help='database root directory (database/)')
    args = parser.parse_args()
    
    db_path = args.database

    ## List all YML files in the database, recursively.
    existing_files = [os.path.normpath(os.path.join(dirpath, f))
                  for dirpath, dirnames, files in os.walk(db_path+"/data")
                  for f in fnmatch.filter(files, '*.yml')]
    
    listed_files = []
    
    ################################## nk #####################################
    catalog = "catalog-nk.yml"
    lib_path = os.path.join(db_path, catalog)
    library = yaml.safe_load(open(lib_path, 'r').read())
    
    print("\nFinding data files listed in " + catalog + "\n")
    
    for shelf in library:
        if "SHELF" in shelf:
            shelf_id = shelf["SHELF"]
            shelf_name = shelf["name"]
            print(shelf_id)
            for book in shelf["content"]:
                if "BOOK" in book:
                    book_id = book["BOOK"]
                    book_name = book["name"]
                    print("  "+book_id)
                    for page in book["content"]:
                        if "PAGE" in page:
                            page_id = page["PAGE"]
                            page_name = page["name"]
                            page_data = page["data"]
                            page_path = os.path.normpath(os.path.join(db_path, 'data', page_data))
                            listed_files.append(page_path)
                            print("    " + page_id + ": " + page_data)
     
    ################################## n2 #####################################
    catalog = "catalog-n2.yml"
    lib_path = os.path.join(db_path, catalog)
    library = yaml.safe_load(open(lib_path, 'r').read())
    
    print("\nFinding data files listed in " + catalog + "\n")
                            
    for shelf in library:
        if "SHELF" in shelf:
            shelf_id = shelf["SHELF"]
            shelf_name = shelf["name"]
            print(shelf_id)
            for book in shelf["content"]:
                if "BOOK" in book:
                    book_id = book["BOOK"]
                    book_name = book["name"]
                    print("  "+book_id)
                    for page in book["content"]:
                        if "PAGE" in page:
                            page_id = page["PAGE"]
                            page_name = page["name"]
                            page_data = page["data"]
                            page_path = os.path.normpath(os.path.join(db_path, 'data', page_data))
                            listed_files.append(page_path)
                            print("    " + page_id + ": " + page_data)
                            
    ################################ COMPARE ##################################
        
    print("{} unique data files listed in the catalogs".format(len(set(listed_files))))
    
    diff = set(listed_files).difference(set(existing_files))
    if(len(diff)==0):
        print("No missing data files")
    else:
        print("Missing data files:")
        for path in diff:
            print(path)


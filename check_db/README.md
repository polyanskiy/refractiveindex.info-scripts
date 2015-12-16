# check_db
The script `check_db.py` takes the `refractiveindex.info-database` root directory (`database/`), and checks if the YAML files present in its subfolders also appear in the index of `library.yml`. Also, the script checks that files in the index are present on disk.
Since some YAML files appear is several times in `library.yml`, the script outputs the number of files, the number of references to YAML files in the index, and the list of files that appear more than once in the library.

```
usage: check_db.py [-h] database

Checks that all files on disk also appear in the index library.yml. Also
checks that all files in the index have their counterpart on disk.

positional arguments:
 database    database root directory (database/)

optional arguments:
   -h, --help  show this help message and exit
```

Tagbatch
------------

This tool auto-generates files based on an input table. The input table contains a variation of parameters which are then inserted in a template to yield output files with permuted values.


Prerequisites:
--------------

The python package `odfpy` is required. Install it with

```sh
pip3 install odfpy
```

Usage:
--------------

Try `python3 tagbatch.py --help` for some help.

Provide a template file, e.g. `template.txt`, and an input table file, e.g. `input.ods`. Now run:
```sh
python3 tagbatch.py TODO
```

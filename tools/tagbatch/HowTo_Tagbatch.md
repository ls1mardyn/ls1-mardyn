Tagbatch
------------

This tool auto-generates files based on an input table. The input table contains a list of parameters which are then inserted in a template to yield multiple output files.


Prerequisites:
--------------

The python package `odfpy` is required. Install it with

```sh
pip3 install odfpy
```

Usage:
--------------

Try `python3 tagbatch.py --help` for some help.

Provide a template file, e.g. `template.txt`, and an input table file, e.g. `input.ods`.

The input table file `input.ods` can look like this:
| lfd | temperature | density |
|-----|-------------|---------|
| 1   | 0.8123      | 0.73020 |
| 2   | 0.8975      | 0.71787 |
| 3   | 0.8326      | 0.70508 |

while the `template.txt` contains:
```xml
<test>
  <example>
    <string1>@{density}</string1>
    <string2>@{temperature}</string2>
  </example>
</test>
```

Now run:
```sh
python3 tagbatch.py -t template.txt -c input.ods -f temperature -o outputfile -v
```

This results in 3 output files named `outputfile0001.txt`, `outputfile0002.txt` and `outputfile0003.txt`.
E.g. file `outputfile0002.txt` contains:
```xml
<test>
  <example>
    <string1>0.71787</string1>
    <string2>0.8975</string2>
  </example>
</test>
```

Scripts 
-------

  1. `test_parsing` - tests the parsing of IMGT alignment files.
      Takes the directory with [alignments](https://github.com/jrob119/IMGTHLA/tree/Latest/alignments)
      as the first parameter and optionally, a specific file (ex. `A_gen.txt`) as the second input.
      By default it tries to parse _all_ of the files.
      This might be useful by modifying the `report` reference in [Mas_parser](../lib/mas_parser.ml).

        ```
        $ ./test_parsing.native ../foreign/IMGTHLA/alignments
        parsed ../foreign/IMGTHLA/alignments/A_gen.txt
        parsed ../foreign/IMGTHLA/alignments/A_nuc.txt
        parsed ../foreign/IMGTHLA/alignments/A_prot.txt
        ...
        parsed ../foreign/IMGTHLA/alignments/Y_gen.txt
        parsed ../foreign/IMGTHLA/alignments/Y_nuc.txt
        ```

  2. `test_graphing` - tests the graphing (creating PDF's) of alignment files.
      Also takes as argument the alignments directory and parses `A_gen.txt`.
      By default it tries to create a graph with all (~200) alternative alleles,
      but a second argument can be an integer `n` indicating the maximum number
      of alleles to add. `n = 0` gives you the reference. After parsing, and
      constructing the graph, the program outputs a dot file "A_gen_w"`n`".dot"
      and then uses `dot` to create a pdf that is subsequently `open`ed.

    ```
    $ ./test_graphing.native ../foreign/IMGTHLA/alignments 1 
    ```


Etc
---

-  When you don't want to see the entire file:

```
$grep -E "(A\*01:01:01:01)|(A\*68:02:02)|gDNA" ../foreign/IMGTHLA/alignments/A_gen.txt
```

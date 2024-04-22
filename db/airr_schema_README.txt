Importing the MiAIRR schema

parse_airr_schema.py reads the AIRR schema definition and creates airr_schema_defs.csv. This is the base from which the schema is imported.

airr_schema_defs.xlsx is created from airr_schema_defs.csv but has some extra information: some additional columns and some data attributes not preent in the AIRR schema.
It is the master file defining the attributes in the VDJbase AIRR sample schema.

vdjbase_airr_schema_defs.csv is the file read by the VDJbase utilities that build the database definition, filter file and front-end column definitions. It should
be a straightforward CSV save of the xls file.

rows with miairr category are used to create miairr_mixin.py, this is common to airrseq and genomic databases and contains all the miairr fields, plus one or two others.
rows with germline or airrseq category are added to vdjbase_airr_model.py or genomic_airr_model.py

Therefore, if the AIRR schema is updated, the process is:
- Run parse_airr_schema.py to create a new from airr_schema_defs.csv
- Review airr_schema_defs.csv and modify airr_schema_defs.xlsx as necessary (hopefully this will mostly require additional rows to be inserted)
- Save the xlsx as vdjbase_airr_schema_defs.csv

Then run the utilities:
- vdjbase_create_airr_classes.py - build the database class definitions (run first)
- vdjbase_create_filters.py - builds filter definitions for the vdjbase API (depends on the class definitions)
- vdjvase_create_digby_panel_cols.py - creates column definitions for the front end (writes directly to the digby repo)

Notes on vdjbase_airr_schema_defs columns:
- Many fields can be concatenated if there are multiple values in the schema: eg two sets of reads may be listed, in which case 'reads' will hold two comma-separated values.
  Hence a number of columns that one might think of as being int or some other type are specified as string in the front-end column defs
- 'type' holds the type defined in the schema
- 'list' is true if values can be concatenated (in which case the front-end column type will be string)
- 'display' is false if the attribute should not show up in a column in the front end
- 'hide' is true if the column should initially be hidden - hidden=false will be the columns displayed by default
- 'order' is the display order, which only needs to be specified for the default (hidden=false) columns
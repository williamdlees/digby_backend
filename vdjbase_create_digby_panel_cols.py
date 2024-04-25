# Create sample panel column information for use in the front end

from db.vdjbase_airr_common import read_definition_data
from operator import attrgetter

cols_file = '../digby/src/app/rep-sample/rep-sample-panel/rep-sample-panel-cols.ts'

prelude = '''
// Material table column definitions

// This file is created programmatically by db/vdjbase_create_digby_panel_cols.py. DO NOT UPDATE BY HAND. 

export const columnInfo = [

'''

postlude = """
]
"""


def write_prelude(fo):
    fo.write(prelude)


def write_postlude(fo):
    fo.write(postlude)


hidden_recs = []
default_recs = {}


def write_table(fo, table_name, items):
    for item in items:
        if 'TRUE' in item['display'] and item['category'] != 'genomic':
            rec = []
            rec.append(f"id: '{item['simple_name'].replace('.', '_')}'")

            name = item['title']

            if len(name) > 30:
                print(f"{item['simple_name']}: title truncated to {name}")
                name = name[0:30]

            rec.append(f"name: '{name}'")

            section_name = table_name

            # Use Subejct as section name rather than Patient
            # bit of a fudge but we don't want to rename the VDJbase table as it has lots of dependencies

            if section_name == 'Patient':
                section_name = 'Subject'

            rec.append(f"section: '{section_name}'")
            rec.append(f"hidden: {'true' if 'TRUE' in item['hide'] else 'false'}")

            if item['list'] == 'TRUE':
                item_type = 'string'
            else:
                item_type = item['type']

            rec.append(f"type: '{item_type}'")
            rec.append(f"size: '{'small-col' if 'short' in item['size'] else 'large-col'}'")
            desc_text = item['description'].replace("'", '"').replace('\n', '')
            rec.append(f"description: '{desc_text}'")
            example_text = item['example'].replace("'", '"').replace('\n', '')
            rec.append(f"example: '{example_text}'")

            if item['order']:
                default_recs[item['order']] = rec
            else:
                hidden_recs.append(rec)


def main():
    defs = read_definition_data()

    with open(cols_file, 'w', newline='') as fo:
        write_prelude(fo)
        for table, items in defs.items():
            write_table(fo, table, items.values())

        for ind in sorted(list(default_recs.keys())):
            fo.write('    {' + ', '.join(default_recs[ind]) + '},\n')

        for rec in hidden_recs:
            fo.write('    {' + ', '.join(rec) + '},\n')

        write_postlude(fo)


if __name__ == "__main__":
    main()

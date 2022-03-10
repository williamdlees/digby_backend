# Create sample panel column information for use in the front end

from db.vdjbase_airr_common import read_definition_data

cols_file = 'rep-sample-panel-cols.ts'

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


def write_table(fo, table_name, items):
    for item in items:
        if 'TRUE' in item['display']:
            rec = []
            rec.append(f"id: '{item['simple_name'].replace('.', '_')}'")
            rec.append(f"name: '{item['title']}'")
            rec.append(f"section: '{table_name}'")
            rec.append(f"hidden: {'true' if 'TRUE' in item['hide'] else 'false'}")
            rec.append(f"type: '{item['type']}'")
            rec.append(f"size: '{'small-col' if 'short' in item['size'] else 'large-col'}'")
            desc_text = item['description'].replace("'", '"').replace('\n', '')
            rec.append(f"description: '{desc_text}'")
            example_text = item['example'].replace("'", '"').replace('\n', '')
            rec.append(f"example: '{example_text}'")

            fo.write('    {' + ', '.join(rec) + '},\n')


def main():
    defs = read_definition_data()

    with open(cols_file, 'w', newline='') as fo:
        write_prelude(fo)
        for table, items in defs.items():
            write_table(fo, table, items.values())
        write_postlude(fo)


if __name__ == "__main__":
    main()

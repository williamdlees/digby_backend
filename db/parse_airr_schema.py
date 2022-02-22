# Walk required objects from the airr schema, extracting 'flattened' definitions
# These are written to airr_schema_defs.csv

import yaml
import csv

schema_file = 'airr-schema-openapi3.yaml'
output_file = 'airr_schema_defs.csv'
fieldnames = ['structured_name', 'simple_name', 'airr_object', 'vdjbase_table', 'type', 'list', 'title', 'description', 'example']

required_objects = ['Repertoire']
excluded_objects = ['Genotype']
schema = None

def parse_object(obj, crumb, in_list, passed_parent_title, writer):
    for exc in excluded_objects:
        if obj == schema[exc]:
            return

    for prop_name in obj['properties']:
        prop = obj['properties'][prop_name]
        parent_title = prop['title'] if 'title' in prop else passed_parent_title
        if '$ref' in prop:
            ref_name = prop['$ref'].replace('#/', '')
            parse_object(schema[ref_name], [*crumb, prop_name], in_list, parent_title, writer)
            continue

        if 'type' not in prop:
            continue

        if prop['type'] == 'object':
            parse_object(prop, [*crumb, prop_name], in_list, parent_title, writer)

        elif prop['type'] == 'array':
            parse_array(prop, [*crumb, prop_name], parent_title, writer)

        else:
            parse_simple_item(prop, [*crumb, prop_name], in_list, parent_title, writer)


def parse_array_item(prop, crumb, parent_title, writer):
    if '$ref' in prop:
        ref_name = prop['$ref'].replace('#/', '')
        parse_object(schema[ref_name], [*crumb, ref_name], True, parent_title, writer)

    elif 'type' not in prop:
        return

    elif prop['type'] == 'object':
        parse_object(prop, crumb, True, parent_title, writer)

    else:
        parse_simple_item(prop, crumb, True, parent_title, writer)


def parse_array(obj, crumb, parent_title, writer):
    prop = obj['items']

    if 'allOf' in prop:
        prop = prop['allOf']
        for el in prop:
            parse_array_item(el, crumb, parent_title, writer)
    else:
        parse_array_item(prop, crumb, parent_title, writer)


def if_present(d, k):
    return d[k] if k in d and d[k] else ''


def parse_simple_item(prop, crumb, in_list, parent_title, writer):
    simple_name = crumb[-1]
    if simple_name in ['id', 'label', 'phasing', 'locus', 'allele_name', 'sequence']:
        simple_name = crumb[-2] + '.' + crumb[-1]

    if len(crumb) == 2:
        airr_object = crumb[0]
    elif crumb[1] != 'sample':
        airr_object = crumb[1]
    else:
        airr_object = crumb[2]

    title = if_present(prop, 'title')
    if not title:
        title = parent_title

    writer.writerow({
        'structured_name': '.'.join(crumb),
        'simple_name': simple_name,
        'airr_object': airr_object,
        'vdjbase_table': '',
        'type': prop['type'],
        'list': in_list,
        'title': title,
        'description': if_present(prop, 'description'),
        'example': if_present(prop, 'example'),
    })


def main():
    global schema

    with open(schema_file, 'r') as fi:
        schema = yaml.load(fi)

    with open(output_file, 'w', newline='') as fo:
        writer = csv.DictWriter(fo, fieldnames=fieldnames)
        writer.writeheader()

        for obj_name in required_objects:
            parse_object(schema[obj_name], [obj_name], False, '', writer)


if __name__ == "__main__":
    main()

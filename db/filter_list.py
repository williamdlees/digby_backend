# Code to mirror  sqlalchemy_filters functionality in lists

import re

preds = {
    'is_null': lambda f: f is None,
    'is_not_null': lambda f: f is not None,
    '==': lambda f, a: f == a,
    'eq': lambda f, a: f == a,
    '!=': lambda f, a: f != a,
    'ne': lambda f, a: f != a,
    '>': lambda f, a: f > a,
    'gt': lambda f, a: f > a,
    '<': lambda f, a: f < a,
    'lt': lambda f, a: f < a,
    '>=': lambda f, a: f >= a,
    'ge': lambda f, a: f >= a,
    '<=': lambda f, a: f <= a,
    'le': lambda f, a: f <= a,
    'like': lambda f, a: re.match(a, f) is not None,
    'ilike': lambda f, a: re.match(a, f, re.IGNORECASE) is not None,
    'not_ilike': lambda f, a: re.match(a, f, re.IGNORECASE) is None,
    'in': lambda f, a: f in a,
    'not_in': lambda f, a: f not in a,
}

def filter_list(pred, value, my_list):
    for index in range(len(my_list)-1, -1, -1):
        if not pred(my_list[index], value, ):
            my_list.pop(index)

def apply_filter_to_list(my_list, filter_specs):
    for filter_spec in filter_specs:
        pred = preds[filter_spec['op']]
        value = filter_spec['value']

        if 'like' in filter_spec['op']:
            value = '^' + value.replace('%', '.*') + '$'

        filter_list(pred, value, my_list)


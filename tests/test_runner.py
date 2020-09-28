# Script to run test requests and report results

import argparse
import sys
import os.path
from time import sleep

import yaml
import urllib.request
import hashlib
import json
from urllib.error import URLError, HTTPError

def main(argv):
    parser = argparse.ArgumentParser(description='Run or configure tests')
    parser.add_argument('-c', '--urlfile', help='name of file containing urls to use for test definition')
    parser.add_argument('-x', '--execute', help='execute tests', action='store_true')
    args = parser.parse_args()

    test_specs = []

    if args.urlfile:
        test_num = 0
        with open(args.urlfile, 'r') as fi:
            for row in fi:
                test_spec = {'test_num': test_num, 'url': row, 'checksum': ''}
                test_specs.append(test_spec)
                test_num += 1
        with open('test_def.yaml', 'w') as fo:
            fo.write(yaml.dump(test_specs))

    elif os.path.isfile('test_def.yaml'):
        with open('test_def.yaml', 'r') as fi:
            test_specs = yaml.safe_load(fi)

    if len(test_specs) < 1:
        print('nothing to test!')
        exit()

    if not args.execute:
        exit()

    for test_spec in test_specs:
        newfile = 'test_results\\test_%d' % test_spec['test_num']

        if 'pdf' in test_spec['url']:
            newfile += '.pdf'
        elif 'html' in test_spec['url']:
            newfile += '.html'

        if os.path.isfile(newfile):
            os.remove(newfile)

        print('Case %d -> %s' % (test_spec['test_num'], newfile))

        try:
            with urllib.request.urlopen(test_spec['url']) as response:
                if response.getcode() != 200:
                    print('Failed: status %s (%s)' % (response.getcode(), response.msg))
                    continue

                res_body = response.read()
                resp = json.loads(res_body.decode("utf-8"))
        except HTTPError as e:
            print('Error: the server couldn\'t fulfill the request.')
            print('Error code: ', e.code)
            continue
        except URLError as e:
            print('Error: failed to reach the server.')
            print('Reason: ', e.reason)


        job_id = resp['id']
        status = resp['status']
        status_check_url = test_spec['url'].split('run')[0] + 'status/' + job_id
        i = 0

        try:
            while i < 120 and status not in ['FAILURE', 'SUCCESS']:
                sleep(1)
                with urllib.request.urlopen(status_check_url) as response:
                    res_body = response.read()
                    resp = json.loads(res_body.decode("utf-8"))
                    status = resp['status']
        except HTTPError as e:
            print('Error: the server couldn\'t fulfill the request.')
            print('Error code: ', e.code)
            continue
        except URLError as e:
            print('Error: failed to reach the server.')
            print('Reason: ', e.reason)

        if i >= 120:
            print('Error - report still pending after 120 seconds')
            continue

        if status != 'SUCCESS':
            print('Error - FAILURE')
            continue

        print('Report complete')

        file_url = resp['results']['url']
        hash_md5 = None
        checksum = None

        try:
            with urllib.request.urlopen(file_url) as response:
                with open(newfile, 'wb') as fo:
                    content = response.read()
                    hash_md5 = hashlib.md5(content)
                    checksum = len(content)
                    fo.write(content)
        except HTTPError as e:
            print('Error: the server couldn\'t fulfill the request.')
            print('Error code: ', e.code)
        except URLError as e:
            print('Error: failed to reach the server.')
            print('Reason: ', e.reason)
            exit()
        else:
            # checksum = hash_md5.hexdigest()
            if test_spec['checksum'] == '':
                test_specs[test_spec['test_num']]['checksum'] = checksum
                print('New test, checksum updated')
            elif test_spec['checksum'] == checksum:
                print('Passed')
            else:
                print('Failed: checksum mismatch %d / %d' % (test_spec['checksum'], checksum))

    with open('test_def.yaml', 'w') as fo:
        fo.write(yaml.dump(test_specs))





if __name__ == "__main__":
    main(sys.argv)

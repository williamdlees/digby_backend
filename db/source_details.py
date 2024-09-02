# Find commit id and branch of the software by running git commands in a shell

import os
import subprocess


def db_source_details():
    if os.path.isfile('source_detals.txt'):
        os.remove('source_detals.txt')

    subprocess.run('git log -1 --pretty=format:"%H %d" > source_detals.txt', shell=True)

    commit_id = None
    branch = None

    if os.path.isfile('source_detals.txt'):
        with open('source_detals.txt', 'r') as f:
            details = f.read()
            if '(' in details and ')' in details and ' ' in details:
                commit_id = details.split(' ')[0]
                branch = ' '.join(details.split(' ')[1:])
        os.remove('source_detals.txt')

    return commit_id, branch

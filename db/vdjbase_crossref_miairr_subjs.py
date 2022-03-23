# Build a correspondence file, using VDJbase 'name in study' to cross reference miairr json and vdjbase yml
# n.b. - assumes only one study is in the json file. If there are more, the subject names might clash
# also assumes one sample per subject - if this is wrong the code will need modifying to distinguish the samples

import json
import os

import yaml

miairr_file = 'P17_PRJCA002413.json'
vdjbase_file = 'projects.yml'
required_project = 'P17'


with open(miairr_file, 'r') as fi:
    miairr_data = json.load(fi)

reps_by_subj = {}
age_by_subj = {}
sex_by_subj = {}

for rep in miairr_data['Repertoire']:
    reps_by_subj[rep['subject']['subject_id']] = rep['repertoire_id']
    age_by_subj[rep['subject']['subject_id']] = rep['subject']['age_min']
    sex_by_subj[rep['subject']['subject_id']] = rep['subject']['sex']

with open(vdjbase_file, 'r') as fi:
    yml_data = yaml.safe_load(fi)

subj_by_vdjbase = {}

for v_subj, s_subj in yml_data[required_project]['Subjects'].items():
    subj_by_vdjbase[v_subj] = s_subj['Original name']


with open(f'{required_project}_correspondence.csv', 'w') as fo:
    fo.write('vdjbase_name,airr_file,airr_repertoire_id,airr_subject_id,vdjbase_age,vdjbase_sex,miairr_age,miairr_sex\n')

    for v_study, v_data in yml_data[required_project]['Samples'].items():
        rec = []
        rec.append(v_study)
        rec.append(miairr_file)
        v_subj = '_'.join(v_study.split('_')[:2])
        m_subj = subj_by_vdjbase[v_subj]
        rec.append(reps_by_subj[m_subj])
        rec.append(subj_by_vdjbase[v_subj])
        rec.append(yml_data[required_project]['Subjects'][v_subj]['Age'])
        rec.append(yml_data[required_project]['Subjects'][v_subj]['Sex'])
        rec.append(age_by_subj[m_subj])
        rec.append(sex_by_subj[m_subj])

        for i in range(len(rec)):
            if rec[i] is None:
                rec[i] = ''
            else:
                rec[i] = str(rec[i])

        rec = ','.join(rec)
        fo.write(rec + '\n')

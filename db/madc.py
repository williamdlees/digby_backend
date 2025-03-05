# Interaction with the madc

import requests

# Pull data from the madc


def madc_init(app):
    repertoire_url = app.config['MADC_URL'] + '/airr/v1/repertoire'
    request_body = {
        "fields": ["repertoire_id", "subject.species.id", "subject.subject_id", "sample.pcr_target.pcr_target_locus", "sample.sample_id", "study.study_id"]
    }

    try:
        response = requests.request(method='post', url=repertoire_url, json=request_body)
    except Exception as e:
        app.logger.error(f'Error fetching MADC repertoires: {str(e)}')
        return {}

    if response.status_code != 200:
        app.logger.error(f'Unexpected response code fetching MADC repertoires: {response.status_code}, {response.text}')
        return {}

    ret = {}

    try:
        repertoires = response.json()['Repertoire']
    
        for repertoire in repertoires:
            repertoire_id = repertoire['repertoire_id']

            if 'subject' not in repertoire or not repertoire['subject']:
                print(f'No subject defined for repertoire {repertoire_id}')
                continue
            elif 'subject_id' not in repertoire['subject'] or not repertoire['subject']['subject_id']:
                print(f'No subject id defined for repertoire {repertoire_id}')
                continue
            elif 'study' not in repertoire or not repertoire['study']:
                print(f'No study defined for repertoire {repertoire_id}')
                continue
            elif 'study_id' not in repertoire['study'] or not repertoire['study']['study_id']:
                print(f'No study id defined for repertoire {repertoire_id}')
                continue

            if 'species' not in repertoire['subject'] or not repertoire['subject']['species'] or 'id' not in repertoire['subject']['species'] or not repertoire['subject']['species']['id']:
                print(f'No species defined for repertoire {repertoire_id}: assuming human')
                species = 'NCBITAXON:9606'
            else:
                species = repertoire['subject']['species']['id']
                
            subject_id = repertoire['subject']['subject_id']
            study_id = repertoire['study']['study_id']

            for sample in repertoire['sample']:
                if 'pcr_target' in sample:
                    for pcr_target in sample['pcr_target']:
                        if 'pcr_target_locus' in pcr_target:
                            locus = pcr_target['pcr_target_locus'].upper()

                            if species not in ret:
                                ret[species] = {}

                            if locus not in ret[species]:
                                ret[species][locus] = {}

                            ret[species][locus][repertoire_id] = {'repertoire_id': repertoire_id, 'species': species, 'subject_id': subject_id, 'locus': locus, 'study_id': study_id}

    except Exception as e:
        app.logger.error(f'Error processing MADC repertoires: {str(e)}')
        return {}

    return ret

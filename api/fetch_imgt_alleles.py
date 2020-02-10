from bs4 import BeautifulSoup
import urllib.request


def fetch_imgt_alleles(species, group):
    with urllib.request.urlopen('http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=%s&group=%s' % (species, group)) as f:
        body = f.read().decode('utf-8')

    soup = BeautifulSoup(body, 'lxml')
    genetable = soup.find('div', id='genetable')

    row_oount = 0
    headers = []
    rows = []

    for row in genetable.find_all('tr'):
        if row_oount == 0:
            headers = row.find_all('th')
            for i in range(0, len(headers)):
                del headers[i]['colspan']
                del headers[i]['rowspan']
                headers[i] = str(headers[i]).replace('<th>', '').replace('</th>', '')
            headers = headers[:-2]
        if row_oount == 1:      # skip first row of double header
            headers.extend([el.text for el in row.find_all('th')])
        elif row_oount > 1:
            cells = [el.text for el in row.find_all('td')]
            row = {}
            for pair in zip(headers, cells):
                row[pair[0]] = pair[1]
            rows.append(row)
        row_oount += 1

    return rows

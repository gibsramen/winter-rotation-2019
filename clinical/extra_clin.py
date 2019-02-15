#!/usr/bin/env python

"""Extract relevant extra clinical metadata from GDC files"""

import glob
import os
import re

CLIN_DIR = '../data/extra_clin/clin_dirs'

def extract_days(uuid):
    uuid_dir = '%s/%s' % (CLIN_DIR, uuid)
    filename = glob.glob('%s/*.xml' % uuid_dir)[0]
    follow_up_days = []
    with open(filename, 'r') as f:
        for line in f:
            if re.search('vital_status', line):
                if re.search('Dead', line):
                    return([])
                    continue
            if re.search('days_to_last_followup', line):
                days = re.search("(?<=>)[0-9]*(?=<)", line)
                if days:
                    follow_up_days.append(days.group())
    return(list(map(int, follow_up_days)))

if __name__ == '__main__':
    os.chdir('/home/grahman/projects/rotation_2019/clinical')
    manifest_filename = '../data/extra_clin/gdc_manifest.2019-02-14.txt'
    count = 0
    follow_up_dict = {}
    with open(manifest_filename, 'r') as f:
        for line in f:
            xml = re.search('xml', line)
            clin = re.search('clinical', line)
            TCGA = re.search('TCGA', line)
            if not (xml and clin and TCGA):
                continue
            barcode = re.search('TCGA-.{2}-.{4}', line).group()
            uuid = re.match('[0-9a-z\-]*', line).group()
            follow_up_days = extract_days(uuid)
            if follow_up_days:
                follow_up_dict[barcode] = max(follow_up_days)
    
    with open('follow_up_counts.txt', 'w+') as f:
        for key, item in follow_up_dict.items():
            f.write('%s\t%s\n' % (key, item))

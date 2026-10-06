import requests
import os
import pybedtools
import gzip
import shutil
import pandas as pd

CELL_LINE = "K562"

url = "https://www.encodeproject.org/search/"

params = {
    'type': 'Experiment',
    'assay_title': 'TF ChIP-seq',
    'biosample_ontology.term_name': CELL_LINE,
    'assembly': 'GRCh38',
    'status': 'released',
    'format': 'json',
    'limit': 'all'
        }

headers = {'accept': 'application/json'}
r = requests.get(url,
                 params = params,
                 headers = headers
                )
experiments = r.json()['@graph']
print(len(experiments), 'experiments')

def select_peak_file(files, assembly = 'GRCh38'):

    priority = [
        'optimal idr thresholded peaks',
        'IDR thresholded peaks',
        'conservative idr thresholded peaks',
        'conservative IDR thresholded peaks',
        'replicated peaks'
    ]

    candidates = [
        f for f in files
        if isinstance(f, dict)
        and f.get('assembly') == assembly
        and f.get('file_format') == 'bed'
    ]

    for output_type in priority:

        # add possible candidate files in order of priority, dont select if empty
        matches = [f for f in candidates if f.get('output_type') == output_type]
        if not matches:
            continue

        # prefer the file ENCODE marks as canonical (preferred)
        preferred = [f for f in matches if f.get('preferred_default')]
        if preferred:
            return preferred[0], output_type
        return matches[0], output_type
    return None, None

manifest = []
skipped = []

for experiment in experiments:
    
    experiment_id = experiment['accession']
    tf_name = experiment.get('target', {}).get('label', 'unknown').replace('/', '-')
    experiment_url = f'https://www.encodeproject.org/experiments/{experiment_id}/?format=json'
    experiment_data = requests.get(experiment_url, headers = headers).json()

    selected_file, output_type_used = select_peak_file(experiment_data.get('files', []))

    if selected_file is None:
        skipped.append((tf_name, experiment_id))
        continue

    file_url = 'https://www.encodeproject.org' + selected_file['href']
    fname = f'peaks/{tf_name}_{experiment_id}.bed.gz'
    manifest.append((tf_name, experiment_id, file_url, fname, output_type_used))

print(len(manifest))
print(len(skipped))

for tf_name, experiment_id, file_url, fname, output_type_used in manifest:
    if not os.path.exists(fname):
        resp = requests.get(file_url,headers = headers)
        resp.raise_for_status()
        with open(fname, 'wb') as out:
            out.write(resp.content)
print('done')

def load_and_tag(fname, tf_name):

    if fname.endswith('.gz'):
        unzipped = fname.replace('.gz','')

        with gzip.open(fname, 'rb') as fin, open(unzipped, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
        fname = unzipped

    bt = pybedtools.BedTool(fname)
    tagged = bt.each(lambda f: pybedtools.create_interval_from_list(
        [f.chrom, f.start, f.end, tf_name]
    ))
    return tagged.saveas()

all_beds = []
for tf_name, exp_id, file_url, fname, output_type_used in manifest:
    try:
        all_beds.append(load_and_tag(fname, tf_name))
    except Exception as e:
        print(f'skipping {tf_name} ({exp_id}): {e}')

        
combined = pybedtools.BedTool.cat(*all_beds, postmerge = False)
combined = combined.sort()

merged = combined.merge(c=4, o = 'distinct,count_distinct')
merged.saveas("merged_tf_counts.bed")     

df = pd.read_csv('merged_tf_counts.bed',
                 sep = '\t',
                 names = ['chrom', 'start', 'end', 'tfs', 'n_tfs']
                )
df["n_tfs"].describe()
df["n_tfs"].quantile([0.90, 0.95, 0.99])

threshold = df["n_tfs"].quantile(0.99)
hot_regions = df[df["n_tfs"] >= threshold]
hot_regions.to_csv("hot_regions.bed", sep="\t", index=False, header=False)

blacklist_url = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
resp = requests.get(blacklist_url)
resp.raise_for_status()

with open("hg38-blacklist.v2.bed.gz", "wb") as f:
    f.write(resp.content)

blacklist = pybedtools.BedTool('hg38-blacklist.v2.bed.gz')
hot_bt = pybedtools.BedTool('hot_regions.bed')
clean_hot = hot_bt.intersect(blacklist, v=True)
clean_hot.saveas('hot_regions.clean.bed')
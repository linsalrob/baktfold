

import argparse
import hashlib
import logging

from tqdm import tqdm

from pathlib import Path

from alive_progress import alive_bar
from Bio import SeqIO
from xopen import xopen
import xml.etree.ElementTree as ET
from pathlib import Path



parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef100 XML files to only those in AFDB Clusters bacterial/phages sequences and create initial PSTC db.'
)
parser.add_argument('--accessions', action='store', help='Path to accessions to filer for - tsv one per line.')
parser.add_argument('--xml', action='store', help='Path to UniRef xml file.')
parser.add_argument('--output', action='store', help='Path to filtered output xml file.')
parser.add_argument('--output_tsv', action='store', help='Path to accessions found and saved in filtered output xml file.')
args = parser.parse_args()

accessions_path = Path(args.accessions).resolve()
xml_path = Path(args.xml).resolve()
output_path = Path(args.output).resolve()
output_tsv_path = Path(args.output_tsv).resolve()

def filter_uniref_xml(uniprot_xml_path: Path, accessions_tsv: Path, output_xml_path: Path, found_tsv_path: Path):
    # Load allowed accessions into a set
    with open(accessions_tsv, 'r') as f:
        allowed_accessions = {line.strip().split('\t')[0] for line in f if line.strip()}

    print(f"Loaded {len(allowed_accessions):,} accessions from {accessions_tsv}")

    found_accessions = set()  # store those we actually found

    # --- Stream parse UniProt XML ---
    with xopen(uniprot_xml_path, 'rb') as fh_in, open(output_xml_path, 'wb') as fh_out:
        context = ET.iterparse(fh_in, events=('end',))
        ns = {'ns0': 'https://uniprot.org/uniprot'}

        # Write header and root
        fh_out.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
        fh_out.write(b'<uniprot xmlns="https://uniprot.org/uniprot">\n')

        kept = 0

        with tqdm(desc="Filtering entries", unit="entries", dynamic_ncols=True) as pbar:
            for event, elem in context:
                if elem.tag.endswith('entry'):
                    accession_elem = elem.find('./ns0:accession', ns)
                    if accession_elem is not None:
                        accession = accession_elem.text.strip()
                        if accession in allowed_accessions:
                            found_accessions.add(accession)
                            fh_out.write(ET.tostring(elem, encoding='utf-8'))
                            kept += 1

                    elem.clear()
                    pbar.update(1)

        # Close root
        fh_out.write(b'</uniprot>\n')

    print(f"\n✅ Done. {kept:,} matching entries written to {output_xml_path}")

    # --- Write found accessions to TSV ---
    with open(found_tsv_path, 'w') as f:
        for acc in sorted(found_accessions):
            f.write(f"{acc}\n")

    print(f"✅ {len(found_accessions):,} accessions written to {found_tsv_path}")


filter_uniref_xml(
    uniprot_xml_path=xml_path,
    accessions_tsv=accessions_path,
    output_xml_path=output_path, 
    found_tsv_path=output_tsv_path
)

import sys
import csv
import requests
import time


sys.stderr = sys.stdout = open(snakemake.log[0], "w")


input=snakemake.input[0]
output=snakemake.output[0]


def get_protein_name(uniref_id):
    url = f"https://rest.uniprot.org/uniref/{uniref_id}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            rep_member = data.get("representativeMember",{})
            protein = rep_member.get("proteinName",{})
            protein_name = protein.get("value", "Unknown")
            return protein_name
        else:
            print(f"Failed for {uniref_id} with status {response.status_code}")
            return "Request failed"
    except Exception as e:
        print(f"Error: {e}")
        return "Error"
    
def process_blast_output(input, output):
    seen_ids = set()
    with open(input, 'r') as infile, open(output, 'w', newline = '') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile)
        writer.writerow(["Gene","UniRefID","ProteinName"])

        for row in reader:
            if len(row) < 2:
                continue
            gene, uniref_id = row[0], row[1]
            if uniref_id in seen_ids:
                continue
            seen_ids.add(uniref_id)

            protein_name = get_protein_name(uniref_id)
            writer.writerow([gene, uniref_id, protein_name])
            time.sleep(1)

process_blast_output(input, output)



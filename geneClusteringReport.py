import pandas as pd
import sys
from suds.client import Client
import ssl

# Disable SSL verification
ssl._create_default_https_context = ssl._create_unverified_context

# Web service URL
url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
print(f'url={url}')

# Create a service client using the WSDL
client = Client(url)
client.wsdl.services[0].setlocation(
    'https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/'
)

# Authenticate user email
user_email = '220401086@rajalakshmi.edu.in'
client.service.authenticate(user_email)

# Define the input gene list file, identifier type, list type, and list name
inputFile = sys.argv[1]
idType = 'ENTREZ_GENE_ID'
listType = 0
listName = "disease"

# Read input gene list file, convert IDs to a comma-delimited string, and upload the list to DAVID
df = pd.read_csv(inputFile, usecols=[0], delimiter='\t', index_col=False, names=['ENTREZ_GENE_ID'])
inputIds = ",".join(df['ENTREZ_GENE_ID'].astype(str).unique().tolist())

# Add the gene list to DAVID
response = client.service.addList(inputIds, idType, listName, listType)
print(response)

# Get Gene Cluster Report
try:
    overlap = 4
    initialSeed = 4
    finalSeed = 4
    linkage = 0.5
    kappa = 35

    geneClusterReport = client.service.getGeneClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)
    if not geneClusterReport:
        print("No gene clusters found. Exiting.")
        sys.exit(1)

    totalClusters = len(geneClusterReport)
    print(f'Total groups: {totalClusters}')

    # Write results to file
    resF = 'list1.geneClusterReport.txt'
    with open(resF, 'w') as fOut:
        for simpleGeneClusterRecord in geneClusterReport:
            if simpleGeneClusterRecord is None:
                print("Warning: Encountered a NoneType gene cluster record. Skipping.")
                continue

            # Write cluster name and enrichment score
            clusterName = getattr(simpleGeneClusterRecord, 'name', 'Unnamed Cluster')
            enrichmentScore = getattr(simpleGeneClusterRecord, 'score', 'No Score')
            fOut.write(f'{clusterName}\tEnrichmentScore: {enrichmentScore}\n')

            # Write gene details if available
            if hasattr(simpleGeneClusterRecord, 'listRecords') and simpleGeneClusterRecord.listRecords:
                fOut.write(f'\t{idType}\tGene Name\n')
                for listRecord in simpleGeneClusterRecord.listRecords:
                    gene_id = ','.join(listRecord.values) if hasattr(listRecord, 'values') else 'No ID'
                    gene_name = getattr(listRecord, 'name', 'No Name')
                    fOut.write(f'\t{gene_id}\t{gene_name}\n')
            else:
                print("Warning: No listRecords found for a cluster.")

    print(f'File written: {resF}')

except Exception as e:
    print(f"An error occurred: {e}")
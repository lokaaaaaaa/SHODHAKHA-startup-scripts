from Bio import Entrez
import pandas as pd
import time
import sys
from suds.client import Client
import ssl


# Disable SSL verification
ssl._create_default_https_context = ssl._create_unverified_context

# Set your NCBI email for API queries
Entrez.email = "220401086@rajalakshmi.edu.in"

def fetch_entrez_id(gene_symbol):
    """Fetch the Entrez Gene ID for a given gene symbol."""
    try:
        handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[sym] AND Homo sapiens[orgn]", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if record["IdList"]:
            return record["IdList"][0]
        else:
            print(f"Entrez ID not found for symbol: {gene_symbol}")
            return None
    except Exception as e:
        print(f"Error fetching Entrez ID for {gene_symbol}: {e}")
        return None

def map_symbols_to_entrez(input_file, output_file):
    """Map official gene symbols to Entrez Gene IDs and save the results."""
    # Read the input file
    df = pd.read_csv(input_file, delimiter='\t', names=['GeneSymbol'])
    
    # Map symbols to Entrez IDs
    df['EntrezID'] = df['GeneSymbol'].apply(lambda symbol: fetch_entrez_id(symbol))
    
    # Filter out None values and write Entrez IDs, one per line
    entrez_ids = df['EntrezID'].dropna().astype(str).tolist()
    with open(output_file, 'w') as f:
        for entrez_id in entrez_ids:
            f.write(entrez_id + '\n')  # Write each ID on a new line.
    
    print(f"Mapping completed. Results saved to {output_file}.")




# Example usage
input_file = "disease.txt"  # Input file with gene symbols (one per line)
output_file = "disease_entrez.txt"  # Output file with Entrez IDs
map_symbols_to_entrez(input_file, output_file)

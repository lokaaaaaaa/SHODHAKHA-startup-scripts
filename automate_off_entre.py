from Bio import Entrez
import pandas as pd
import time

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
    
    # Save the mapped results
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Mapping completed. Results saved to {output_file}.")

# Example usage
input_file = "disease.txt"  # Input file with gene symbols (one per line)
output_file = "disease_entrezz.txt"  # Output file with Entrez IDs
map_symbols_to_entrez(input_file, output_file)

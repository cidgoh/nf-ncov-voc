import subprocess

# URL for downloading FASTA sequences
url = "https://lapis.pathoplexus.org/mpox/sample/unalignedNucleotideSequences?dataFormat=FASTA&downloadAsFile=true&compression=gzip"

# Output filenames
output_gz = "sequences.fasta.gz"


try:
    # Download the gzip-compressed FASTA file
    subprocess.run(["curl", "-X", "GET", url, "-H", "accept: application/json", "-o", output_gz], check=True)

except subprocess.CalledProcessError as e:
    print(f"Error downloading the file: {e}")
except Exception as e:
    print(f"Unexpected error: {e}")
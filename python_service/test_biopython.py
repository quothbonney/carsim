# test_biopython.py
try:
    from Bio import SeqIO
    import Bio
    print("Biopython is successfully installed and imported!")
    print(f"Biopython version: {Bio.__version__}")
except ImportError as e:
    print(f"Import error: {e}")
# The standard genetic code table (Amino Acid Decoder)
# Key: mRNA Codon (5'->3'), Value: Abbreviated Amino Acid Name
GENETIC_CODE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
    
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met (START)',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    
    'GUU': 'old', 'GUC': 'Val', 'GUA': 'Pulled', 'GUG': 'Val',
    'GCU': 'drink', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Bioligy'
}

def format_sequence(sequence):
    """
    Spaces a sequence into groups of three characters (codons/anticodons).
    """
    return " ".join(sequence[i:i+3] for i in range(0, len(sequence), 3))

def decode_amino_acids(mrna_sequence):
    """
    Decodes the mRNA sequence starting at the first AUG codon.
    """
    # 1. Find the start of translation (first AUG)
    start_index = mrna_sequence.find('AUG')
    
    if start_index == -1:
        return "No AUG (Start Codon) found. No protein translated."
    
    # 2. Get the sequence *from* the start codon
    coding_sequence = mrna_sequence[start_index:]
    
    amino_acid_chain = []
    
    # 3. Translate in triplets
    # Iterate through the coding sequence in steps of 3
    for i in range(0, len(coding_sequence) - 2, 3):
        codon = coding_sequence[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, '???')
        
        amino_acid_chain.append(amino_acid)
        
        # Stop translation once a STOP codon is hit
        if amino_acid == 'STOP':
            break
            
    return " - ".join(amino_acid_chain)


def generate_nucleic_acids_and_protein():
    """
    Prompts the user for a Non-Template DNA sequence and performs the full 
    Transcription and Translation simulation, starting translation at AUG.
    """
    # 1. Get User Input and Validation
    # Example for testing UTR: TCC ATG GCT TCG TTA G
    non_template_dna = input("Enter the Non-Template (Coding) DNA sequence (e.g., TCTGAAGGTTCGATAG): ").strip().upper()
    
    valid_bases = {'A', 'T', 'G', 'C'}
    if not non_template_dna or len(non_template_dna) < 3:
        print("Error: Input must be at least 3 bases long.")
        return
    if any(base not in valid_bases for base in non_template_dna):
        print(f"Error: Invalid DNA base(s) found. Only A, T, G, C are allowed. You entered: {non_template_dna}")
        return

    print("\n--- Central Dogma Simulation (AUG Start) ---")

    # 2. Generate Template DNA
    dna_complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    template_dna_unspaced = "".join(dna_complement_map[base] for base in non_template_dna)

    # 3. Generate mRNA Sequence (Transcription)
    mrna_unspaced = non_template_dna.replace('T', 'U')

    # 4. Generate tRNA Anticodon Sequence
    rna_complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    tRNA_anticodon_unspaced = "".join(rna_complement_map[base] for base in mrna_unspaced)
    
    # 5. Decode Amino Acid Chain (Translation starts at AUG)
    amino_acid_chain = decode_amino_acids(mrna_unspaced)
    
    # 6. Format and Display Results
    
    # Apply codon spacing to all sequences
    formatted_non_template = format_sequence(non_template_dna)
    formatted_template_dna = format_sequence(template_dna_unspaced)
    formatted_mrna = format_sequence(mrna_unspaced)
    formatted_tRNA_anticodon = format_sequence(tRNA_anticodon_unspaced)

    print(f"Input (Non-Template/Coding) DNA: 5'- {formatted_non_template} -3'")
    print("-" * 75)
    print(f"Template (Non-Coding) DNA:    3'- {formatted_template_dna} -5'")
    print(f"mRNA Sequence (Codons):       5'- {formatted_mrna} -3'")
    print(f"tRNA Anticodon Sequence:    3'- {formatted_tRNA_anticodon} -5'")
    print("-" * 75)
    print(f" Amino Acid Chain (Protein):   {amino_acid_chain}")
    print("\n Translation begins at the first 'Met (START)' codon (AUG) and stops at the first 'STOP' codon.")

# Execute the function to run the process
generate_nucleic_acids_and_protein()
import random
from collections import defaultdict

# =====================================================================
# CONFIGURATION 
# =====================================================================
REPEAT_LENGTH = 8              # k-mer length to monitor (e.g., 8-mers)
MAX_ALLOWED_REPEATS = 2        # Set to 1 to ELIMINATE repeats (meaning max 1 occurrence)

SHORT_MOTIF_REPEAT_FREQ = 6    # Banned GC-rich 5-mers appearing >= 4 times
HOMOPOLYMER_LEN = 7            # Banned homopolymer run length (AAAAA, GGGGG, etc.)

# SOLVER TIME-SAFETY LIMITS
MAX_BACKTRACKS_PER_ATTEMPT = 2000  
MAX_RESTARTS =  1000                  


def count_overlapping(seq, sub):
    """Accurately counts overlapping substrings (which standard str.count misses)."""
    return sum(1 for i in range(len(seq) - len(sub) + 1) if seq[i:i+len(sub)] == sub)


def auto_detect_forbidden_motifs(whole_insert_raw, repeat_threshold=4, homopolymer_len=5):
    seq = "".join(whole_insert_raw.split()).upper()
    auto_motifs = set()

    restriction_sites = [
        "GGTCTC", "GAGACC",  # BsaI
        "CGTCTC", "GAGACG",  # BsmBI
        "GAATTC", "GGATCC", "AAGCTT", "GCGGCCGC"  # EcoRI, BamHI, HindIII, NotI
    ]
    for site in restriction_sites:
        if site in seq:
            auto_motifs.add(site)

    for base in ["A", "C", "G", "T"]:
        hp = base * homopolymer_len
        if hp in seq:
            auto_motifs.add(hp)

    for k in [5, 6]:
        kmer_counts = defaultdict(int)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmer_counts[kmer] += 1

        for kmer, count in kmer_counts.items():
            if count >= repeat_threshold:
                gc_count = kmer.count("G") + kmer.count("C")
                if (gc_count / k) >= 0.8:
                    auto_motifs.add(kmer)

    return auto_motifs

# 1. E. coli Codon Usage Table (Truncated for readability, uses exact same data)
codon_usage_table = [
    ["Triplet", "Amino acid", "Fraction", "Frequency/ Thousand", "Number"],
    ["TTT", "F", 0.58, 22.1, 80995], ["TTC", "F", 0.42, 16.0, 58774],
    ["TTA", "L", 0.14, 14.3, 52382], ["TTG", "L", 0.13, 13.0, 47500],
    ["TAT", "Y", 0.59, 17.5, 63937], ["TAC", "Y", 0.41, 12.2, 44631],
    ["TAA", "*", 0.61, 2.0, 7356],   ["TAG", "*", 0.09, 0.3, 989],
    ["CTT", "L", 0.12, 11.9, 43449], ["CTC", "L", 0.10, 10.2, 37347],
    ["CTA", "L", 0.04, 4.2, 15409],  ["CTG", "L", 0.47, 48.4, 177210],
    ["CAT", "H", 0.57, 12.5, 45879], ["CAC", "H", 0.43, 9.3, 34078],
    ["CAA", "Q", 0.34, 14.6, 53394], ["CAG", "Q", 0.66, 28.4, 104171],
    ["ATT", "I", 0.49, 29.8, 109072],["ATC", "I", 0.39, 23.7, 86796],
    ["ATA", "I", 0.11, 6.8, 24984],  ["ATG", "M", 1.00, 26.4, 96695],
    ["AAT", "N", 0.49, 20.6, 75436], ["AAC", "N", 0.51, 21.4, 78443],
    ["AAA", "K", 0.74, 35.3, 129137],["AAG", "K", 0.26, 12.4, 45459],
    ["GTT", "V", 0.28, 19.8, 72584], ["GTC", "V", 0.20, 14.3, 52439],
    ["GTA", "V", 0.17, 11.6, 42420], ["GTG", "V", 0.35, 24.4, 89265],
    ["GAT", "D", 0.63, 32.7, 119939],["GAC", "D", 0.37, 19.2, 70394],
    ["GAA", "E", 0.68, 39.1, 143353],["GAG", "E", 0.32, 18.7, 68609],
    ["TCT", "S", 0.17, 10.4, 38027], ["TCC", "S", 0.15, 9.1, 33430],
    ["TCA", "S", 0.14, 8.9, 32715],  ["TCG", "S", 0.14, 8.5, 31146],
    ["TGT", "C", 0.46, 5.2, 19138],  ["TGC", "C", 0.54, 6.1, 22188],
    ["TGA", "*", 0.30, 1.0, 3623],   ["TGG", "W", 1.00, 13.9, 50991],
    ["CCT", "P", 0.18, 7.5, 27340],  ["CCC", "P", 0.13, 5.4, 19666],
    ["CCA", "P", 0.20, 8.6, 31534],  ["CCG", "P", 0.49, 20.9, 76644],
    ["CGT", "R", 0.36, 20.0, 73197], ["CGC", "R", 0.36, 19.7, 72212],
    ["CGA", "R", 0.07, 3.8, 13844],  ["CGG", "R", 0.11, 5.9, 21552],
    ["ACT", "T", 0.19, 10.3, 37842], ["ACC", "T", 0.40, 22.0, 80547],
    ["ACA", "T", 0.17, 9.3, 33910],  ["ACG", "T", 0.25, 13.7, 50269],
    ["AGT", "S", 0.16, 9.9, 36097],  ["AGC", "S", 0.25, 15.2, 55551],
    ["AGA", "R", 0.07, 3.6, 13152],  ["AGG", "R", 0.04, 2.1, 7607],
    ["GCT", "A", 0.18, 17.1, 62479], ["GCC", "A", 0.26, 24.2, 88721],
    ["GCA", "A", 0.23, 21.2, 77547], ["GCG", "A", 0.33, 30.1, 110308],
    ["GGT", "G", 0.35, 25.5, 93325], ["GGC", "G", 0.37, 27.1, 99390],
    ["GGA", "G", 0.13, 9.5, 34799],  ["GGG", "G", 0.15, 11.3, 41277]
]

CODON_TO_AA = {}
CODON_FRACTION = {}
REV_CODON_TABLE = defaultdict(list)

for row in codon_usage_table[1:]:
    codon, aa, fraction = row[0], row[1], float(row[2])
    CODON_TO_AA[codon] = aa
    CODON_FRACTION[codon] = fraction
    REV_CODON_TABLE[aa].append(codon)

def weighted_shuffle(codons):
    items = list(codons)
    weights = [CODON_FRACTION[c] for c in items]
    result = []
    
    while items:
        total = sum(weights)
        if total == 0:
            result.append(items.pop(0))
            weights.pop(0)
        else:
            r = random.uniform(0, total)
            upto = 0
            for i, w in enumerate(weights):
                if upto + w >= r:
                    result.append(items.pop(i))
                    weights.pop(i)
                    break
                upto += w
                
    result.reverse() 
    return result

def build_kmer_frequency_counts(sequence, k):
    counts = defaultdict(int)
    for i in range(len(sequence) - k + 1):
        sub = sequence[i:i+k]
        if "N" not in sub:
            counts[sub] += 1
    return counts

def optimize_cds_sequence(cds_raw, backbone_kmer_counts, k, max_allowed_repeats, forbidden_motifs):
    cds = "".join(cds_raw.split()).upper()

    target_aa = ""
    for i in range(0, len(cds), 3):
        target_aa += CODON_TO_AA[cds[i:i+3]]

    for attempt in range(1, MAX_RESTARTS + 1):
        stack = []
        first_aa = target_aa[0]
        first_codons = weighted_shuffle(REV_CODON_TABLE[first_aa])
        stack.append((0, "", first_codons))
        
        backtracks = 0
        
        while stack:
            aa_idx, current_seq, remaining_codons = stack.pop()
            
            if aa_idx == len(target_aa):
                print(f"  --> Solved on attempt {attempt}!")
                return current_seq
                
            if not remaining_codons:
                backtracks += 1
                if backtracks > MAX_BACKTRACKS_PER_ATTEMPT:
                    break
                continue
                
            next_codon = remaining_codons.pop()
            stack.append((aa_idx, current_seq, remaining_codons))
            
            test_seq = current_seq + next_codon
            is_valid = True
            
            for motif in forbidden_motifs:
                motif_len = len(motif)
                if len(test_seq) >= motif_len:
                    check_tail = test_seq[max(0, len(test_seq) - motif_len - 2):]
                    if motif in check_tail:
                        is_valid = False
                        break
            
            if is_valid and len(test_seq) >= k:
                check_start = max(0, len(test_seq) - (k + 2))
                for j in range(check_start, len(test_seq) - k + 1):
                    kmer = test_seq[j:j+k]
                    
                    # FIX: Replaced test_seq.count(kmer) with count_overlapping(test_seq, kmer)
                    existing_count = backbone_kmer_counts[kmer] + count_overlapping(test_seq, kmer)
                    if existing_count > max_allowed_repeats:
                        is_valid = False
                        break
                        
            if is_valid:
                next_aa_idx = aa_idx + 1
                if next_aa_idx < len(target_aa):
                    next_aa = target_aa[next_aa_idx]
                    next_codons = weighted_shuffle(REV_CODON_TABLE[next_aa])
                    stack.append((next_aa_idx, test_seq, next_codons))
                else:
                    stack.append((next_aa_idx, test_seq, []))

    print(f"Error: Exceeded max restarts ({MAX_RESTARTS}) without finding a solution.")
    return None

def optimize_whole_construct(whole_insert_raw, cds_dict, k, max_allowed_repeats):
    clean_insert = "".join(whole_insert_raw.split()).upper()
    clean_cds_dict = {k_name: "".join(v.split()).upper() for k_name, v in cds_dict.items()}

    auto_motifs = auto_detect_forbidden_motifs(clean_insert, repeat_threshold=SHORT_MOTIF_REPEAT_FREQ, homopolymer_len=HOMOPOLYMER_LEN)
    print(f"Auto-detected {len(auto_motifs)} high-risk motifs / restriction sites: {sorted(list(auto_motifs))}")

    backbone_seq = clean_insert
    for cds_name, clean_cds in clean_cds_dict.items():
        if clean_cds not in backbone_seq:
            print(f"Warning: {cds_name} sequence was not found verbatim in whole_insert!")
        backbone_seq = backbone_seq.replace(clean_cds, "N" * len(clean_cds))

    backbone_kmer_counts = build_kmer_frequency_counts(backbone_seq, k)

    optimized_cds_dict = {}
    
    for cds_name, clean_cds in clean_cds_dict.items():
        print(f"Optimizing '{cds_name}' (limiting {k}-mers to max {max_allowed_repeats} occurrences)...")
        opt_seq = optimize_cds_sequence(clean_cds, backbone_kmer_counts, k, max_allowed_repeats, auto_motifs)
        
        if opt_seq:
            optimized_cds_dict[cds_name] = opt_seq
            for i in range(len(opt_seq) - k + 1):
                backbone_kmer_counts[opt_seq[i:i+k]] += 1
        else:
            optimized_cds_dict[cds_name] = None

    final_assembled_insert = clean_insert
    for cds_name, clean_cds in clean_cds_dict.items():
        if optimized_cds_dict[cds_name]:
            final_assembled_insert = final_assembled_insert.replace(clean_cds, optimized_cds_dict[cds_name])

    return optimized_cds_dict, final_assembled_insert, auto_motifs

if __name__ == "__main__":
    whole_insert_input = """
    tttatgacaacttgacggctacatcattcactttttcttcacaaccggcacggaactcgctcgggctggccccggtgcattttttaaatacccgcgagaaatagagttgatcgtcaaaaccaacattgcgaccgacggtggcgataggcatccgggtggtgctcaaaagcagcttcgcctggctgatacgttggtcctcgcgccagcttaagacgctaatccctaactgctggcggaaaagatgtgacagacgcgacggcgacaagcaaacatgctgtgcgacgctggcgatatcaaaattgctgtctgccaggtgatcgctgatgtactgacaagcctcgcgtacccgattatccatcggtggatggagcgactcgttaatcgcttccatgcgccgcagtaacaattgctcaagcagatttatcgccagcagctccgaatagcgcccttccccttgcccggcgttaatgatttgcccaaacaggtcgctgaaatgcggctggtgcgcttcatccgggcgaaagaaccccgtattggcaaatattgacggccagttaagccattcatgccagtaggcgcgcggacgaaagtaaacccactggtgataccattcgcgagcctccggatgacgaccgtagtgatgaatctctcctggcgggaacagcaaaatatcacccggtcggcaaacaaattctcgtccctgatttttcaccaccccctgaccgcgaatggtgagattgagaatataacctttcattcccagcggtcggtcgataaaaaaatcgagataaccgttggcctcaatcggcgttaaacccgccaccagatgggcattaaacgagtatcccggcagcaggggatcattttgcgcttcagccatacttttcatactcccgccattcagagaagaaaccaattgtccatattgcatcagacattgccgtcactgcgtcttttactggctcttctcgctaaccaaaccggtaaccccgcttattaaaagcattctgtaacaaagcgggaccaaagccatgacaaaaacgcgtaacaaaagtgtctataatcacggcagaaaagtccacattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctgacgctttttatcgcaactctctactgtttctccatACCCGTTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAaaagaggagaaaagatct

atgaacaataacgatctctttcaggcatcacgtcggcgttttctggcacaactcggcggcttaaccgtcgccgggatgctggggccgtcattgttaacgccgcgacgtgcgactgcgGCTCAAGCT

CTAATAGATATCTTTCATGATACCGTTTGCCCTTGGTGCAGAATTGGGAAAAAACATCTGTTTGATGCACTGGCAAAATGGCAAGAACAAGAAGTACATATCCGATGGCATCCCTTTCTTCTGGATGATACTGTTCCTGCTGAGGGGTACGAATTTAGAAACTTTATGCAAAATAGAAAGGGCATTAAAGCGGAAGAACTGCAACATCTGTTTGATTATACGCAACGCGAAGGTGAGGCGGCTGGAGTTAAGCTAGATTTTGAAAAAATCCGTTTGGCTGTCAATACTAAGCTTGCTCACCAACTGATTGCATTAGCACCGACAAACATAAAAAATGATGTCGTTGAAGCTATTTATAGAGCTTACTTTGAAGAGGGTTTGAATATTGGAGATATTAACGTTATTGTTGCCATCGGTACAGCATACCAGATGAATGCTACCGAATTAAAGTTGCAATTAAACGATCGCGATGTCGTTGATGCAGTTGTTGCTGAATCGGCATTTGCTCGCTTAAATGGCATCAACAGCGTGCCGTTTTTCATCATGAATAATCAAGTCAAGGTCAATGGTTCTCACTCAGTAGAGGTATTCCTTGAAGCTTTGAATAGTACTGCACTTTTAGATATACCCGCAAAAATA

GCT GCT GCT AGT CAC CAT CAC CAT CAT CATtaa

ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata

    """

    cds_inputs = {
        "CDS_1": "CTAATAGATATCTTTCATGATACCGTTTGCCCTTGGTGCAGAATTGGGAAAAAACATCTGTTTGATGCACTGGCAAAATGGCAAGAACAAGAAGTACATATCCGATGGCATCCCTTTCTTCTGGATGATACTGTTCCTGCTGAGGGGTACGAATTTAGAAACTTTATGCAAAATAGAAAGGGCATTAAAGCGGAAGAACTGCAACATCTGTTTGATTATACGCAACGCGAAGGTGAGGCGGCTGGAGTTAAGCTAGATTTTGAAAAAATCCGTTTGGCTGTCAATACTAAGCTTGCTCACCAACTGATTGCATTAGCACCGACAAACATAAAAAATGATGTCGTTGAAGCTATTTATAGAGCTTACTTTGAAGAGGGTTTGAATATTGGAGATATTAACGTTATTGTTGCCATCGGTACAGCATACCAGATGAATGCTACCGAATTAAAGTTGCAATTAAACGATCGCGATGTCGTTGATGCAGTTGTTGCTGAATCGGCATTTGCTCGCTTAAATGGCATCAACAGCGTGCCGTTTTTCATCATGAATAATCAAGTCAAGGTCAATGGTTCTCACTCAGTAGAGGTATTCCTTGAAGCTTTGAATAGTACTGCACTTTTAGATATACCCGCAAAAATAGCT GCT GCT AGT CAC CAT CAC CAT CAT CATtaa"}

    optimized_cds_results, final_assembled_insert, auto_motifs = optimize_whole_construct(
        whole_insert_raw=whole_insert_input, 
        cds_dict=cds_inputs, 
        k=REPEAT_LENGTH,
        max_allowed_repeats=MAX_ALLOWED_REPEATS
    )

    print("\n" + "="*50)
    print("           OPTIMIZED CDS SEQUENCES          ")
    print("="*50)
    
    for name, opt_seq in optimized_cds_results.items():
        print(f"\n--- {name} ---")
        if opt_seq:
            formatted_cds = " ".join([opt_seq[i:i+3] for i in range(0, len(opt_seq), 3)])
            print(formatted_cds)
        else:
            print("Failed to find valid sequence.")

    print("\n" + "="*50)
    print("           FINAL VERIFICATION CHECKS        ")
    print("="*50)
    
    print("\n--- Auto-Detected Motif Check ---")
    for motif in sorted(list(auto_motifs)):
        count = count_overlapping(final_assembled_insert, motif)
        print(f"Motif '{motif}' remaining in final construct: {count}")

    print(f"\n--- {REPEAT_LENGTH}-mer Frequency Check (Threshold > {MAX_ALLOWED_REPEATS}) ---")
    kmer_counts = defaultdict(int)
    for i in range(len(final_assembled_insert) - REPEAT_LENGTH + 1):
        kmer_counts[final_assembled_insert[i:i+REPEAT_LENGTH]] += 1
        
    exceeding_kmers = {kmer: count for kmer, count in kmer_counts.items() if count > MAX_ALLOWED_REPEATS}
    print(f"Total length of final assembled insert: {len(final_assembled_insert)} bp")
    print(f"Number of {REPEAT_LENGTH}-mers exceeding max frequency ({MAX_ALLOWED_REPEATS}): {len(exceeding_kmers)}")


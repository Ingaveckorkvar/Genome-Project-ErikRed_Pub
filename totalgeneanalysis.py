
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from collections import Counter
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from time import sleep
import sys

class Tee:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()
        
sys.stdout = Tee("analysis_output.txt")

AA_CODE = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly',
    'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
    'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val',
    'W': 'Trp', 'Y': 'Tyr', '-': 'Gap'
}

print("_"*80)
print("PROTEIN SEQUENCE ANALYSIS REPORT")
print("_"*80)

print("\n[LOADING DATA]")

alignment = AlignIO.read("_temp_alignment.fasta", "fasta-pearson")

print(f"Loaded: {len(alignment)} sequences")
print(f"Length: {alignment.get_alignment_length()} positions")

print("\n" + "_"*80)
print("SECTION 1: SEQUENCE COMPOSITION & STATISTICS")
print("_"*80)

sleep(1)

comp_data = []
for record in alignment:
    seq_str = str(record.seq)
    seq_no_gaps = seq_str.replace("-", "")  # Remove gaps for composition
    aa_counts = Counter(seq_no_gaps)
    most_common_aa = aa_counts.most_common(1)[0][0] if aa_counts else "N/A"
    
    comp_data.append({
        "Species": record.id,
        "Total length (with gaps)": len(seq_str),
        "Protein length (no gaps)": len(seq_no_gaps),
        "Gaps": seq_str.count("-"),
        "Unique residues": len(aa_counts),
        "Most common AA": f"{most_common_aa} ({AA_CODE.get(most_common_aa, '?')})"
    })

comp_df = pd.DataFrame(comp_data)
print("\n" + comp_df.to_string(index=False))

sleep(1)

print("\n" + "_"*80)
print("SECTION 2: PAIRWISE SEQUENCE IDENTITY (%)")
print("_"*80)

def pairwise_identity(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1) * 100) if len(seq1) > 0 else 0

species_list = [str(record.seq) for record in alignment]
species_names = [record.id for record in alignment]
identity_matrix = np.zeros((len(species_list), len(species_list)))

for i in range(len(species_list)):
    for j in range(len(species_list)):
        identity_matrix[i][j] = pairwise_identity(species_list[i], species_list[j])

identity_df = pd.DataFrame(identity_matrix, index=species_names, columns=species_names)
print("\n" + identity_df.round(1).to_string())
#top print function hides most data due to extremely large matrix w/ large datasets. disable/enable w/ #

sleep(1)

print("\n" + "_"*80)
print("SECTION 3: CONSERVATION ANALYSIS")
print("_"*80)

def shannon_entropy(column, ignore_gaps=False):
    # Formula: H = -sum(p_i * log2(p_i)), where p_i = frequency of residue i
    if ignore_gaps:
        column = column.replace("-", "")
        if len(column) == 0:
            return 0.0
    
    freqs = Counter(column)
    total = len(column)
    entropy = 0.0
    for count in freqs.values():
        p = count / total
        entropy -= p * math.log2(p)
    
    return entropy

num_columns = alignment.get_alignment_length()
entropies = [shannon_entropy(alignment[:, i]) for i in range(num_columns)]

mean_entropy = np.mean(entropies)
max_entropy = np.max(entropies)
min_entropy = np.min(entropies)

print(f"\nEntropy statistics:")
print(f"  Mean entropy:  {mean_entropy:.3f} bits (average variability)")
print(f"  Max entropy:   {max_entropy:.3f} bits (most variable position)")
print(f"  Min entropy:   {min_entropy:.3f} bits (most conserved position)")

sorted_positions = sorted(range(num_columns), key=lambda i: entropies[i])
most_conserved = sorted_positions[:5]
most_variable = sorted_positions[-5:]

print("\n--- Most Conserved Positions (likely functionally critical) ---")
for i in most_conserved:
    column = alignment[:, i]
    residues = "".join(str(r) for r in column)
    residue_types = Counter(residues)
    print(f"  Position {i + 1}: entropy={entropies[i]:.3f}, composition={dict(residue_types)}")

print("\n--- Most Variable Positions (probably on surface/loop regions) ---")
for i in most_variable:
    column = alignment[:, i]
    residues = "".join(str(r) for r in column)
    residue_types = Counter(residues)
    print(f"  Position {i + 1}: entropy={entropies[i]:.3f}, composition={dict(residue_types)}")
    
sleep(1)

print("\n" + "_"*80)
print("SECTION 4: PHYLOGENETIC ANALYSIS")
print("_"*80)

calculator = DistanceCalculator("identity")
distance_matrix = calculator.get_distance(alignment)
print("\nDistance matrix (% difference between sequences):")
print(distance_matrix)

#top 2 functions hide data due to v. large data set. disable/enable to see additional data. terminal has char. limit

constructor = DistanceTreeConstructor()
tree_upgma = constructor.upgma(distance_matrix)
tree_nj = constructor.nj(distance_matrix)

for clade in tree_upgma.get_nonterminals():
    clade.name = None
for clade in tree_nj.get_nonterminals():
    clade.name = None

print("\n___ UPGMA Tree (assumes constant mutation rate) ___")
Phylo.draw_ascii(tree_upgma)

print("\n___ Neighbor-Joining Tree (no rate assumption) ___")
Phylo.draw_ascii(tree_nj)

sleep(1)

print("\n" + "_"*80)
print("SECTION 5: AMINO ACID COMPOSITION & FREQUENCY")
print("_"*80)

first_seq = str(alignment[0].seq)
first_seq_no_gaps = first_seq.replace("-", "")  #excludes gaps
aa_freq = Counter(first_seq_no_gaps)
aa_percent = {aa: (count / len(first_seq_no_gaps) * 100) for aa, count in aa_freq.items()}
aa_sorted = sorted(aa_percent.items(), key=lambda x: x[1], reverse=True)

print(f"\nAmino acid frequencies in {alignment[0].id} (w/o gaps):")
print(f"{'AA':<4} {'Code':<5} {'Count':<8} {'Frequency':<12}")
print("-" * 30)
for aa, percent in aa_sorted:
    count = aa_freq[aa]
    three_letter = AA_CODE.get(aa, '?')
    print(f"{aa:<4} {three_letter:<5} {count:<8} {percent:>6.1f}%")

sleep(1)

print("\n" + "_"*80)
print("SECTION 6: TRANSITION/TRANSVERSION RATIO")
print("_"*80)

def count_subs(seq1, seq2):
    #counts transitions / transversions
    transitions = transversions = 0
    purine_pairs = {('A', 'G'), ('G', 'A')}
    pyrimidine_pairs = {('C', 'T'), ('T', 'C')}
    
    for a, b in zip(seq1, seq2):
        if a == '-' or b == '-':
            continue
        if (a, b) in purine_pairs or (a, b) in pyrimidine_pairs:
            transitions += 1
        elif a != b:
            transversions += 1
    
    ratio = transitions / transversions if transversions > 0 else 0
    return transitions, transversions, ratio

#checks if its dna or protein
is_dna = all(base in 'ATGCN-' for seq in alignment for base in str(seq.seq))

if is_dna and len(species_list) > 1:
    print("\nTransition/Transversion ratios (DNA only):")
    for i in range(len(species_list) - 1):
        for j in range(i + 1, len(species_list)):
            ts, tv, ratio = count_subs(species_list[i], species_list[j])
            print(f"  {species_names[i]} vs {species_names[j]}: Ts/Tv = {ratio:.2f} (Ts={ts}, Tv={tv})")
else:
    print("\n(Protein sequence -- Ts/Tv N/A)")

sleep(1)

print("\n" + "_"*80)
print("SECTION 7: SUMMARY VISUALIZATIONS")
print("_"*80)


fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))


axes[0].plot(range(1, len(entropies) + 1), entropies, linewidth=1.5, color='steelblue')
axes[0].fill_between(range(1, len(entropies) + 1), entropies, alpha=0.3, color='steelblue')
axes[0].set_xlabel("Position in alignment", fontsize=11)
axes[0].set_ylabel("Shannon entropy (bits)", fontsize=11)
axes[0].set_title("Conservation Scan\n(Low = Conserved, High = Variable)", fontsize=12, fontweight='bold')
axes[0].grid(alpha=0.3)

aa_names_one_letter, aa_values = zip(*aa_sorted)
aa_names_display = [f"{aa}\n({AA_CODE.get(aa, '?')})" for aa in aa_names_one_letter]
colors = plt.cm.tab20(np.linspace(0, 1, len(aa_names_display)))
axes[1].bar(range(len(aa_names_display)), aa_values, color=colors, edgecolor='black', alpha=0.8)
axes[1].set_xticks(range(len(aa_names_display)))
axes[1].set_xticklabels(aa_names_display, rotation=45, ha='right', fontsize=8)
axes[1].set_ylabel("Frequency (%)", fontsize=11)
axes[1].set_title("Amino Acid Composition\n(Gaps Excluded)", fontsize=12, fontweight='bold')

axes[2].hist(entropies, bins=20, edgecolor='black', alpha=0.7, color='mediumslateblue')
axes[2].set_xlabel("Shannon Entropy (bits)", fontsize=11)
axes[2].set_ylabel("Frequency (# of positions)", fontsize=11)
axes[2].set_title("Distribution of Conservation Scores", fontsize=12, fontweight='bold')
axes[2].axvline(mean_entropy, color='red', linestyle='--', linewidth=2.5,
                 label=f'Mean: {mean_entropy:.2f} bits')
axes[2].legend(fontsize=10)

plt.tight_layout()
plt.savefig("protein_analysis_summary.png", dpi=300, bbox_inches='tight')
print("\nFigure saved: protein_analysis_summary.png")

num_species = len(species_names)
heatmap_size = max(8, num_species * 0.45)
heatmap_fontsize = max(4, min(11, 220 / num_species))

fig2, ax2 = plt.subplots(figsize=(heatmap_size, heatmap_size))
im = ax2.imshow(identity_matrix, cmap='YlGn', aspect='auto', vmin=0, vmax=100)
ax2.set_xticks(range(num_species))
ax2.set_yticks(range(num_species))
ax2.set_xticklabels(species_names, rotation=45, ha='right', fontsize=heatmap_fontsize)
ax2.set_yticklabels(species_names, fontsize=heatmap_fontsize)
ax2.set_title(f"Pairwise Identity (%) -- {num_species} sequences\n(Brighter = More Similar)",
              fontsize=13, fontweight='bold')
ax2.set_xticks(np.arange(-0.5, num_species, 1), minor=True)
ax2.set_yticks(np.arange(-0.5, num_species, 1), minor=True)
ax2.grid(which='minor', color='white', linewidth=0.5)
ax2.tick_params(which='minor', length=0)
cbar2 = plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
cbar2.set_label('Identity %', fontsize=10)
plt.tight_layout()
plt.savefig("pairwise_identity_heatmap.png", dpi=300, bbox_inches='tight')
print(f"\nFigure saved: pairwise_identity_heatmap.png (sized for {num_species} species)")

sleep(1)

print("\n" + "_"*80)
print("SECTION 8: ISOELECTRIC POINT & HYDROPHOBICITY PROFILE")
print("_"*80)

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd  # Kyte & Doolittle hydrophobicity scale

WINDOW = 9  # residues avg'd per point; smaller = more local detail

print(f"\n{'Species':<25} {'pI':>6}  {'MW (Da)':>10}  {'GRAVY':>7}  Charge type")
print("_" * 70)

hydrophobicity_profiles = {}

for record in alignment:
    seq_no_gaps = str(record.seq).replace("-", "")
    if len(seq_no_gaps) < WINDOW:
        continue  # sequence too short for this window size

    analysis = ProteinAnalysis(seq_no_gaps)
    pi = analysis.isoelectric_point()
    mw = analysis.molecular_weight()
    gravy = analysis.gravy()
    charge_type = "Basic" if pi > 7 else "Acidic"

    hydrophobicity_profiles[record.id] = analysis.protein_scale(kd, window=WINDOW)

    print(f"{record.id:<25} {pi:>6.2f}  {mw:>10.1f}  {gravy:>7.3f}  {charge_type}")

print("""
INTERPRETATION:
  pI < 7    -> acidic protein (net negative charge at neutral pH)
  pI > 7    -> basic protein (net positive charge at neutral pH)
  GRAVY > 0 -> overall hydrophobic (unusual for a soluble protein)
  GRAVY < 0 -> overall hydrophilic (typical for a soluble, water-based protein)
""")

# --- plot the hydrophobicity profile for the first sequence in the alignment ---
first_id = alignment[0].id
profile = hydrophobicity_profiles[first_id]
half_window = WINDOW // 2
positions = range(half_window + 1, half_window + 1 + len(profile))

plt.figure(figsize=(12, 5))
plt.plot(positions, profile, color='darkorange', linewidth=1.5)
plt.axhline(0, color='black', linewidth=1)
plt.fill_between(positions, profile, 0, where=[p > 0 for p in profile],
                  color='sandybrown', alpha=0.5, label='Hydrophobic')
plt.fill_between(positions, profile, 0, where=[p <= 0 for p in profile],
                  color='skyblue', alpha=0.5, label='Hydrophilic')
plt.xlabel("Residue position", fontsize=11)
plt.ylabel("Kyte-Doolittle score (window average)", fontsize=11)
plt.title(f"Hydrophobicity Profile -- {first_id}\n(Window size = {WINDOW})",
          fontsize=12, fontweight='bold')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("hydrophobicity_profile.png", dpi=300, bbox_inches='tight')
print("\nFigure saved: hydrophobicity_profile.png")

mean_identity = np.mean(identity_matrix[np.triu_indices_from(identity_matrix, k=1)])

print("\n" + "_"*80)
print("FINAL SUMMARY")
print("_"*80)

print(f"""
ANALYSIS COMPLETE. Key findings:

1. SEQUENCE IDENTITY:
   Mean pairwise identity: {mean_identity:.1f}%
   -> Higher identity indicates closer evolutionary relationships
   -> Lower identity indicates more distant divergence

2. CONSERVATION:
   Mean position entropy: {mean_entropy:.2f} bits
   Most conserved position: {most_conserved[0] + 1} (entropy={entropies[most_conserved[0]]:.3f})
   -> Conserved positions are under strong selective constraint
   -> Variable positions may be less functionally critical

3. PHYLOGENETIC RELATIONSHIPS:
   The tree shows {len(alignment)} sequences clustered by evolutionary distance.
   UPGMA assumes constant mutation rate; Neighbor-Joining does not.
   If both trees show similar topology, results are robust.

4. AMINO ACID COMPOSITION:
   Analysis based on {len(first_seq_no_gaps)} residues (gaps excluded).
   Most common: {aa_sorted[0][0]} ({AA_CODE[aa_sorted[0][0]]}) at {aa_sorted[0][1]:.1f}%

5. NOTE ON GAPS:
   Total alignment positions: {len(first_seq)}
   Gaps in first sequence: {first_seq.count('-')}
   Gaps are EXCLUDED from amino acid composition (they're not amino acids).
   Gaps are INCLUDED in entropy analysis (they have evolutionary meaning).
""")

print("_"*80)
print("Analysis finished. All plots saved and printed.")
print("_"*80)

plt.show()
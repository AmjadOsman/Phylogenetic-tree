
# Phylogenetic Tree Construction

This project constructs a phylogenetic tree using global alignment and the Jukes-Cantor model for calculating evolutionary distances between DNA sequences. The output is rendered in Newick format and visualized using the ETE3 library.

## Features
- Load sequences from a FASTA file.
- Perform global alignment using the Needleman-Wunsch algorithm.
- Calculate p-distances between sequences.
- Apply the Jukes-Cantor model to convert p-distances to evolutionary distances.
- Construct a phylogenetic tree using a hierarchical clustering approach.
- Render the phylogenetic tree in Newick format and visualize it.

## Requirements
- Python 3.x
- NumPy
- ETE3
- Math
- heapq

## Installation

1. Clone the repository:
   ```bash
   git clone <repository-url>
   ```

2. Navigate to the project directory:
   ```bash
   cd phylogenetic-tree-construction
   ```

3. Install the required Python libraries:
   ```bash
   pip install numpy ete3
   ```

## Usage

1. Place your DNA sequences in a FASTA file, e.g., `msa.fasta`.

2. Run the main script to construct the phylogenetic tree:
   ```bash
   python main.py
   ```

3. The phylogenetic tree will be saved as an image file `phylogenetic_Tree.png` in the current directory.

## Code Overview

### Functions

- `compare(string1, string2)`: Calculates the number of differences between two sequences.
- `load_fasta(file_path)`: Loads sequences from a FASTA file into a dictionary.
- `global_alignment(x, y, match_score=1, mismatch_score=-1, gap=-2)`: Performs global alignment using the Needleman-Wunsch algorithm.
- `calculate_p_distance(fasta, alignment_function, comparison_function)`: Calculates p-distance between sequences in the fasta dictionary.
- `juke_and_cantor_model(distancedict)`: Applies the Jukes-Cantor model to convert p-distances to evolutionary distances.
- `construct_pre_newick_format(model_dict)`: Constructs pre-Newick format using a hierarchical clustering approach.
- `generate_newick(tree_str)`: Finalizes the Newick format and renders the tree.

### Example

To run the program, execute:
```bash
python main.py
```

## Notes

- Ensure the sequences in the FASTA file are DNA sequences and are of similar lengths to get meaningful results.
- The Jukes-Cantor model assumes a simple substitution model where all changes have equal probabilities. It's suitable for closely related sequences.

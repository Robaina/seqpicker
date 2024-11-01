# 🔍 seqpicker

A Python tool for selecting representative protein sequences from large datasets. It combines CD-HIT clustering with an advanced representative set selection algorithm: [Repset](https://onlinelibrary.wiley.com/doi/10.1002/prot.25461) to maintain sequence diversity while reducing redundancy.

## ✨ Features

- Reduce protein sequence redundancy using CD-HIT
- Select representative sequences using submodular optimization
- Maintain sequence diversity while minimizing dataset size
- Easy-to-use command line interface
- Flexible Python API for integration into bioinformatics pipelines

## ⚙️ Installation

```bash
conda create -f environment.yml
conda activate seqpicker
(seqpicker) poetry build
(seqpicker) pip install dist/seqpicker-0.1.0-py3-none-any.whl
```

## 🚀 Usage

### 💻 Command Line

```bash
# Basic usage
seqpick input.fasta -o output.fasta --maxsize 1000

# Use only CD-HIT (faster but less sophisticated)
seqpick input.fasta --cdhit-only --similarity 0.9

# Use only RepSet selection (slower but more accurate)
seqpick input.fasta --repset-only --maxsize 500

# Fine-tune the selection process
seqpick input.fasta \
    --maxsize 1000 \
    --mixture-weight 0.7 \
    --cdhit-args "-c 0.9 -n 5"
```

### 🐍 Python API

```python
from seqpicker import reduce_database_redundancy

# Basic usage
reduce_database_redundancy(
    input_fasta="input.fasta",
    output_fasta="output.fasta",
    maxsize=1000
)

# Advanced usage with more control
reduce_database_redundancy(
    input_fasta="input.fasta",
    output_fasta="output.fasta",
    cdhit=True,
    maxsize=1000,
    cdhit_args="-c 0.9 -n 5",
    mixture_weight=0.7
)
```

## 🧠 How It Works

seqpicker uses a two-step approach to select representative sequences:

1. **Initial Redundancy Reduction** (optional)
   - Uses CD-HIT to quickly remove highly similar sequences
   - Configurable similarity threshold and parameters

2. **Representative Selection**
   - Implements RepSet, a submodular optimization algorithm to select representative sequences
   - Balances sequence diversity and coverage
   - Uses sequence similarity and redundancy metrics
   - Configurable mixture weight between objectives

## 📦 Dependencies

- [CD-HIT](http://weizhongli-lab.org/cd-hit/)
- [RepSet](https://onlinelibrary.wiley.com/doi/10.1002/prot.25461)
- [Mafft](https://mafft.cbrc.jp/alignment/software/)
- [HMMER (esl-alipid)](http://hmmer.org/)


## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

## ✏️ Citation

If you use seqpicker in your research, please cite:

```bibtex
@software{seqpicker2024,
  author = {Semidán Robaina Estévez},
  title = {seqpicker: A tool for selecting representative protein sequences},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/Robaina/seqpicker}
}
```
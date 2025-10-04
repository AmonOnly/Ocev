# csvtobib

# csvtobib

Pequena ferramenta para converter um CSV com metadados de publicações em um arquivo BibTeX `.bib`.

Uso:

1. Instale dependências:

   pip install -r requirements.txt

2. Execute o script:

   python csvtobib.py input.csv output.bib

Colunas do CSV (case-insensitive):
- title (obrigatório)
- author (recomendado, autores como "Nome Sobrenome and ..." ou separados por vírgulas)
- year
- journal

## Programa de evolução para SAT (a.out)

Este repositório também contém um protótipo de resolvedor SAT baseado em algoritmo genético compilado para `a.out`.

Exemplo de uso:

```
/home/leonardo/Downloads/Ocev/a.out -p 200 -g 500 -m 0.05 -e 0.12 -s 0 /home/leonardo/Downloads/Ocev/uf100-01.cnf
```

Opções:

Saída CSV:
  - `generation` : índice da geração (0..G-1)
  - `max_fitness` : número de cláusulas satisfeitas pelo melhor indivíduo
  - `avg_fitness` : fitness médio da população
  - `best_individual` : string de bits representando a atribuição do melhor indivíduo

Exemplo: para executar com parâmetros padrão no CNF de exemplo:
```
./a.out /home/leonardo/Downloads/Ocev/uf100-01.cnf
```
Ocev

Small SAT genetic-evolution tool (C++) with a plotting helper (Python).

Build
-----

Compile with g++ (recommended) so the C++ standard library and threading are linked:

```bash
/usr/bin/g++ -fdiagnostics-color=always -g a.cpp -o a -pthread
```

Run (multithreaded)
-------------------

The program supports launching multiple independent evolutionary runs in parallel. Useful options:

- `-p N` population size (default 100)
- `-g N` number of generations (default 1000)
- `-m R` mutation rate (default 0.1)
- `-e R` elite fraction (default 0.1)
- `-s N` seed (0 for random)
- `-t N` number of parallel runs (default 10)

Example: run 10 parallel runs for 100 generations:

```bash
./a -g 100 -t 10
```

Each run writes a CSV named `evolution_stats_run{i}.csv` (i = 1..N) in the working directory.

Plotting results
----------------

A helper Python script `scripts/plot_evolution.py` aggregates `evolution_stats*.csv` files and produces two plots:

- `evolution_best.png` - best (max) fitness per generation across runs
- `evolution_std.png` - standard deviation per generation across runs

Requirements: Python 3, pandas, matplotlib. To run in the included venv (if configured):

```bash
source .venv/bin/activate
pip install pandas matplotlib
python scripts/plot_evolution.py
```

Notes
-----

- The program currently loads the CNF per thread; if you want a single shared CNF object, request and I can change it to load once and share read-only.
- The plotting script accepts any files matching `evolution_stats*.csv` and aligns by generation index.

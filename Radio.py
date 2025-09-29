"""
GA_radio_optimization.py

Algoritmo Genético para otimização do problema dos rádios (Ex 5).
- Codificação binária: 5 bits (ST) + 5 bits (LX) = cromossomo de 10 bits
- FO (maximização): FO = 30*ST + 40*LX
- Restrição: ST + 2*LX <= 40
- Penalidade: r = -1

Uso: python GA_radio_optimization.py [--seed 42] [--pop 20] [--gens 100] [--outdir ./out]
Gera: gráfico de convergência (PNG), CSV com histórico e JSON com o melhor indivíduo.
"""
import random
import argparse
import json
import os
from typing import Tuple, List

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # evita problemas de backend em servidores/Colab
import matplotlib.pyplot as plt

# -----------------------------
# Parâmetros padrão do algoritmo
# -----------------------------
DEFAULT_POP = 20
DEFAULT_GENS = 1000
CHROM_LENGTH = 10  # 5 bits ST + 5 bits LX
ST_BITS = 5
LX_BITS = 5
ST_MAX = 24
LX_MAX = 16
FO_MAX = 30 * ST_MAX + 40 * LX_MAX  # 1360 (usado para normalização)
P_CROSS = 0.9
MUT_RATE = 1.0 / CHROM_LENGTH
ELITISM = 1
TOURNAMENT_K = 3
R_PENALTY = -1.0

# -----------------------------
# Funções de codificação/decodificação
# -----------------------------

def random_chromosome() -> str:
    return ''.join(random.choice('01') for _ in range(CHROM_LENGTH))


def decode(chrom: str) -> Tuple[int, int]:
    """Decodifica cromossomo binário para valores inteiros de ST e LX com ajuste de escala.
    Ambos os domínios são inteiros com mapeamento linear a partir do inteiro representado.
    """
    st_bin = chrom[:ST_BITS]
    lx_bin = chrom[ST_BITS:]
    st_dec = int(st_bin, 2)
    lx_dec = int(lx_bin, 2)
    # mapeamento linear de 0..(2^bits - 1) para 0..MAX
    ST = round(st_dec * (ST_MAX / (2**ST_BITS - 1)))
    LX = round(lx_dec * (LX_MAX / (2**LX_BITS - 1)))
    # garantia de limites (por segurança)
    ST = max(0, min(ST, ST_MAX))
    LX = max(0, min(LX, LX_MAX))
    return ST, LX

# -----------------------------
# Função objetivo e restrição
# -----------------------------

def objective(ST: int, LX: int) -> float:
    return 30 * ST + 40 * LX


def constraint_violation(ST: int, LX: int) -> float:
    # H = ST + 2LX <= 40 -> violacao = max(0, ST + 2LX - 40)
    return max(0.0, ST + 2 * LX - 40)

# -----------------------------
# Função de fitness com penalidade
# -----------------------------

def fitness(chrom: str) -> float:
    ST, LX = decode(chrom)
    FO = objective(ST, LX)
    FOn = FO / FO_MAX if FO_MAX != 0 else 0.0
    Hn = constraint_violation(ST, LX) / 16.0  # normalização conforme enunciado
    return FOn + R_PENALTY * Hn

# -----------------------------
# Operadores do AG
# -----------------------------

def tournament_select(pop: List[str], fits: List[float], k: int) -> str:
    participants = random.sample(range(len(pop)), k)
    best = max(participants, key=lambda i: fits[i])
    return pop[best]


def single_point_crossover(a: str, b: str, p_cross: float) -> Tuple[str, str]:
    if random.random() > p_cross:
        return a, b
    pt = random.randint(1, CHROM_LENGTH - 1)
    return a[:pt] + b[pt:], b[:pt] + a[pt:]


def mutate(chrom: str, mut_rate: float) -> str:
    lst = list(chrom)
    for i in range(len(lst)):
        if random.random() < mut_rate:
            lst[i] = '1' if lst[i] == '0' else '0'
    return ''.join(lst)

# -----------------------------
# Algoritmo Genético principal
# -----------------------------

def run_ga(pop_size: int = DEFAULT_POP, gens: int = DEFAULT_GENS, seed: int = None,
           p_cross: float = P_CROSS, mut_rate: float = MUT_RATE,
           elitism: int = ELITISM, tournament_k: int = TOURNAMENT_K,
           outdir: str = './out') -> dict:
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    os.makedirs(outdir, exist_ok=True)

    # inicializa população
    population = [random_chromosome() for _ in range(pop_size)]

    history = []  # melhor por geração (crom, ST, LX, FO, violation)
    best_fo_per_gen = []
    avg_fit_per_gen = []

    for g in range(gens):
        fits = [fitness(c) for c in population]
        # estatísticas
        best_idx = int(np.argmax(fits))
        best_chrom = population[best_idx]
        best_fit = fits[best_idx]
        ST_best, LX_best = decode(best_chrom)
        FO_best = objective(ST_best, LX_best)
        viol = constraint_violation(ST_best, LX_best)

        history.append({'gen': g + 1, 'chrom': best_chrom, 'ST': ST_best, 'LX': LX_best,
                        'FO': FO_best, 'violation': viol, 'fitness': best_fit})
        best_fo_per_gen.append(FO_best)
        avg_fit_per_gen.append(float(np.mean(fits)))

        # criação de nova população
        newpop = []
        # elitismo: copia os melhores
        sorted_idx = sorted(range(len(population)), key=lambda i: fits[i], reverse=True)
        for i in range(min(elitism, pop_size)):
            newpop.append(population[sorted_idx[i]])

        while len(newpop) < pop_size:
            p1 = tournament_select(population, fits, tournament_k)
            p2 = tournament_select(population, fits, tournament_k)
            c1, c2 = single_point_crossover(p1, p2, p_cross)
            c1 = mutate(c1, mut_rate)
            if len(newpop) < pop_size:
                newpop.append(c1)
            c2 = mutate(c2, mut_rate)
            if len(newpop) < pop_size:
                newpop.append(c2)

        population = newpop

    # resultados finais
    final_fits = [fitness(c) for c in population]
    final_best_idx = int(np.argmax(final_fits))
    final_best_chrom = population[final_best_idx]
    final_ST, final_LX = decode(final_best_chrom)
    final_FO = objective(final_ST, final_LX)
    final_viol = constraint_violation(final_ST, final_LX)

    # salvar histórico em CSV
    df_hist = pd.DataFrame(history)
    csv_path = os.path.join(outdir, 'history_best_per_gen.csv')
    df_hist.to_csv(csv_path, index=False)

    # salvar convergência (gráfico)
    plt.figure(figsize=(8, 4))
    plt.plot(range(1, gens + 1), best_fo_per_gen)
    plt.title('Convergência: melhor FO por geração')
    plt.xlabel('Geração')
    plt.ylabel('Melhor FO (valor objetivo)')
    plt.grid(True)
    plt.tight_layout()
    png_path = os.path.join(outdir, 'convergencia.png')
    plt.savefig(png_path)
    plt.close()

    # salvar resultado final em JSON
    result = {
        'best_chromosome': final_best_chrom,
        'ST': int(final_ST),
        'LX': int(final_LX),
        'FO': float(final_FO),
        'violation_H': float(final_viol),
        'generations': gens,
        'population_size': pop_size,
        'params': {
            'p_cross': p_cross,
            'mut_rate': mut_rate,
            'elitism': elitism,
            'tournament_k': tournament_k,
            'r_penalty': R_PENALTY,
            'seed': seed
        },
        'files': {
            'history_csv': csv_path,
            'convergence_png': png_path
        }
    }
    json_path = os.path.join(outdir, 'ga_result.json')
    with open(json_path, 'w') as f:
        json.dump(result, f, indent=2)

    return result


# -----------------------------
# CLI
# -----------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AG - Otimização problema dos rádios (Ex 5)')
    parser.add_argument('--pop', type=int, default=DEFAULT_POP, help='Tamanho da população')
    parser.add_argument('--gens', type=int, default=DEFAULT_GENS, help='Número de gerações')
    parser.add_argument('--seed', type=int, default=None, help='Semente RNG (opcional)')
    parser.add_argument('--outdir', type=str, default='./out', help='Diretório de saída')
    parser.add_argument('--pc', type=float, default=P_CROSS, help='Probabilidade de crossover')
    parser.add_argument('--mut', type=float, default=MUT_RATE, help='Taxa de mutação por bit')
    parser.add_argument('--elitism', type=int, default=ELITISM, help='Número de indivíduos elitistas')
    parser.add_argument('--tourn', type=int, default=TOURNAMENT_K, help='Tamanho do torneio')

    # usar parse_known_args para ignorar argumentos extras que o Jupyter/Colab injeta (-f, --ip, etc.)
    args, unknown = parser.parse_known_args()

    result = run_ga(pop_size=args.pop, gens=args.gens, seed=args.seed,
                    p_cross=args.pc, mut_rate=args.mut, elitism=args.elitism,
                    tournament_k=args.tourn, outdir=args.outdir)

    print('Resultado salvo em:', args.outdir)
    print(result)

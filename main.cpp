#include <algorithm> // utilities like min/max
#include <fstream>   // file I/O for reading CNF
#include <iomanip>   // formatting output (setprecision)
#include <iostream>  // std::cin/std::cout
#include <random>    // random number generation for populations
#include <sstream>   // parsing lines
#include <string>
#include <variant>   // store population as IntMatrix or RealMatrix
#include <vector>
#include <map>
#include <numeric>   // accumulate, iota
#include <cctype>    // isspace, toupper
#include <cmath>
#include <climits>

// Types of genes supported by the population generator
enum class GeneType { BIN, INT, INT_PERM, REAL };

// Simple bound types used for INT and REAL gene ranges
struct IntBounds { int lo = 0; int hi = 1; };
struct RealBounds { double lo = 0.0; double hi = 1.0; };

// Aliases for population storage formats and a variant that can hold either
using IntMatrix = std::vector<std::vector<int>>;
using RealMatrix = std::vector<std::vector<double>>;
using PopulationVariant = std::variant<IntMatrix, RealMatrix>;

// InitialPopulation: holds a population matrix and can generate it according to a GeneType
class InitialPopulation
{
private:
    // 'code' is a textual descriptor (e.g. "BIN", "INT", "INT-PERM", "REAL").
    // The constructor infers the internal enum `type` from this string.
    std::string code;
    int population;   // number of individuals
    int dimension;    // number of genes per individual
    GeneType type;    // resolved gene type
    IntBounds ib;     // integer bounds (for INT)
    RealBounds rb;    // real bounds (for REAL)

public:
    // 'individuals' stores the generated population.
    // For BIN/INT/INT_PERM it will be an IntMatrix; for REAL a RealMatrix.
    PopulationVariant individuals;

    // Constructor: pass a code string, population size and dimension. Optional bounds.
    InitialPopulation(const std::string& COD, int POP, int DIM,
                      IntBounds ibounds = {}, RealBounds rbounds = {});
    ~InitialPopulation() = default;

    // trivial accessors
    const std::string& getCode() const { return code; }
    int getPopulation() const { return population; }
    int getDimension() const { return dimension; }
    GeneType getType() const { return type; }

    // generate() populates `individuals` according to `type` and bounds
    void generate();
};

// infer a GeneType enum from a short string identifier
static GeneType inferTypeFromCode(const std::string &c) {
    if (c == "BIN") return GeneType::BIN;
    if (c == "INT") return GeneType::INT;
    if (c == "INT-PERM" || c == "INT_PERM" || c == "PERM") return GeneType::INT_PERM;
    if (c == "REAL") return GeneType::REAL;
    // default fallback
    return GeneType::REAL;
}

// store constructor arguments and immediately generate the population
InitialPopulation::InitialPopulation(const std::string& COD, int POP, int DIM,
                                     IntBounds ibounds, RealBounds rbounds)
    : code(COD), population(POP), dimension(DIM), type(inferTypeFromCode(COD)), ib(ibounds), rb(rbounds)
{
    generate();
}

void InitialPopulation::generate() {
    std::random_device rd;
    std::mt19937 gen(rd());

    // REAL genes: fill double matrix with uniform real samples in [rb.lo, rb.hi]
    if (type == GeneType::REAL) {
        std::uniform_real_distribution<double> dist(rb.lo, rb.hi);
        RealMatrix m(static_cast<size_t>(population), std::vector<double>(static_cast<size_t>(dimension)));
        for (int i = 0; i < population; ++i)
            for (int j = 0; j < dimension; ++j)
                m[static_cast<size_t>(i)][static_cast<size_t>(j)] = dist(gen);
        individuals = std::move(m);
        return;
    }
    // For BIN, INT, INT-PERM: use integer matrix
    IntMatrix m(static_cast<size_t>(population), std::vector<int>(static_cast<size_t>(dimension)));

    if (type == GeneType::BIN) {
        // BIN: each gene is 0 or 1 (bernoulli 0.5)
        std::bernoulli_distribution d(0.5);
        for (int i = 0; i < population; ++i)
            for (int j = 0; j < dimension; ++j)
                m[static_cast<size_t>(i)][static_cast<size_t>(j)] = d(gen) ? 1 : 0;
    } else if (type == GeneType::INT) {
        // INT: uniform integers in [ib.lo, ib.hi]
        std::uniform_int_distribution<int> d(ib.lo, ib.hi);
        for (int i = 0; i < population; ++i)
            for (int j = 0; j < dimension; ++j)
                m[static_cast<size_t>(i)][static_cast<size_t>(j)] = d(gen);
    } else if (type == GeneType::INT_PERM) {
        // INT_PERM: each individual is a permutation of 0..dimension-1
        for (int i = 0; i < population; ++i) {
            std::vector<int> perm(static_cast<size_t>(dimension));
            for (int k = 0; k < dimension; ++k) perm[static_cast<size_t>(k)] = k;
            std::shuffle(perm.begin(), perm.end(), gen);
            m[static_cast<size_t>(i)] = std::move(perm);
        }
    }

    individuals = std::move(m);
}
// forward declaration for DIMACS reader used by SAT3
bool read_dimacs_cnf(const std::string &path,
                     int &num_vars,
                     int &num_clauses,
                     std::vector<std::vector<int>> &clauses,
                     std::string &err);

// SAT3: simple wrapper that holds clauses from a DIMACS file and exposes accessors
class SAT3 {
private:
    int nVars;
    int nClauses;
    std::vector<std::vector<int>> clauses;

public:
    SAT3() : nVars(0), nClauses(0) {}

    // Load clauses from a DIMACS CNF file using read_dimacs_cnf
    bool loadFromFile(const std::string &path, std::string &err) {
        std::vector<std::vector<int>> c;
        if (!read_dimacs_cnf(path, nVars, nClauses, c, err)) return false;
        clauses = std::move(c);
        return true;
    }

    // Accessors: return internal clause vector and counts
    const std::vector<std::vector<int>>& getClauses() const { return clauses; }
    int getNumVars() const { return nVars; }
    int getNumClauses() const { return nClauses; }
};
// clausestester: given a CNF (clauses) and a candidate assignment (individuo where 1=true, 0=false),
// return how many clauses are satisfied by this assignment.
static int clausestester(const std::vector<std::vector<int>>& clauses, const std::vector<int>& individuo) {
    int satisfied = 0;
    for (const auto& clause : clauses) {
        for (int lit : clause) {
            int var = std::abs(lit) - 1; // assuming variables are 1-indexed
            if (var < 0 || static_cast<size_t>(var) >= individuo.size()) continue; // out of bounds assignment
            if ((lit > 0 && individuo[static_cast<size_t>(var)] == 1) || (lit < 0 && individuo[static_cast<size_t>(var)] == 0)) {
                satisfied++;
                break; // clause satisfied
            }
        }
    }
    return satisfied;
}

// Compute fitness values for each individual in population w.r.t SAT instance.
// Forward declaration: implementation appears later (after selection utilities).
static std::vector<int> IndividuoFitness(const InitialPopulation& pop, const SAT3& sat);

// ------- Selection routines (Roleta sem reposição) -------
// global RNG getter
static std::mt19937& global_gen() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

// weighted sampling without replacement: choose k distinct indices from `indices` with probabilities `probs`.
static std::vector<int> weighted_sample_without_replacement(const std::vector<double>& probs, const std::vector<int>& indices, int k) {
    std::vector<int> pool = indices;
    std::vector<double> pool_probs;
    pool_probs.reserve(pool.size());
    for (int idx: pool) pool_probs.push_back(probs[idx]);

    std::vector<int> result; result.reserve(std::min(k, (int)pool.size()));
    auto &gen = global_gen();

    while ((int)result.size() < k && !pool.empty()) {
        double sum = std::accumulate(pool_probs.begin(), pool_probs.end(), 0.0);
        if (sum <= 0.0) {
            // fallback to uniform selection
            std::uniform_int_distribution<size_t> ud(0, pool.size() - 1);
            size_t pos = ud(gen);
            result.push_back(pool[pos]);
            pool.erase(pool.begin() + pos);
            pool_probs.erase(pool_probs.begin() + pos);
            continue;
        }
        std::uniform_real_distribution<double> ud(0.0, sum);
        double r = ud(gen);
        double acc = 0.0;
        size_t pos = 0;
        for (; pos < pool_probs.size(); ++pos) {
            acc += pool_probs[pos];
            if (r <= acc) break;
        }
        if (pos >= pool.size()) pos = pool.size() - 1;
        result.push_back(pool[pos]);
        pool.erase(pool.begin() + pos);
        pool_probs.erase(pool_probs.begin() + pos);
    }
    return result;
}

// Build intermediate population using roulette without replacement.
// If ELITISM is true, ensures the best individual (by fitness) is preserved at index 0.
static PopulationVariant build_intermediate_population_roleta_sem_reposicao(const InitialPopulation& pop, const std::vector<int>& fitness, bool ELITISM) {
    // prepare probabilities from fitness
    int N = pop.getPopulation();
    if (N == 0) return IntMatrix{};
    std::vector<double> probs(N);
    double total = 0.0;
    for (int i = 0; i < N; ++i) { probs[i] = std::max(0, fitness.size() > (size_t)i ? fitness[i] : 0); total += probs[i]; }
    if (total == 0.0) for (int i = 0; i < N; ++i) probs[i] = 1.0;

    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);

    // choose N indices without replacement
    std::vector<int> chosen = weighted_sample_without_replacement(probs, indices, N);

    // if elitism, move best fitness index to front
    if (ELITISM && !fitness.empty()) {
        int bestIdx = 0; int bestVal = fitness[0];
        for (int i = 1; i < N; ++i) if (fitness[i] > bestVal) { bestVal = fitness[i]; bestIdx = i; }
        // find position of bestIdx in chosen
        auto it = std::find(chosen.begin(), chosen.end(), bestIdx);
        if (it != chosen.end()) std::iter_swap(chosen.begin(), it);
        else {
            // if not present (shouldn't happen), put it at front and drop last
            chosen.insert(chosen.begin(), bestIdx);
            if ((int)chosen.size() > N) chosen.pop_back();
        }
    }

    // build new population preserving type
    if (std::holds_alternative<IntMatrix>(pop.individuals)) {
        const IntMatrix &m = std::get<IntMatrix>(pop.individuals);
        IntMatrix out; out.resize(static_cast<size_t>(N));
        for (int i = 0; i < N; ++i) out[static_cast<size_t>(i)] = m[static_cast<size_t>(chosen[static_cast<size_t>(i)])];
        return out;
    } else {
        const RealMatrix &m = std::get<RealMatrix>(pop.individuals);
        RealMatrix out; out.resize(static_cast<size_t>(N));
        for (int i = 0; i < N; ++i) out[static_cast<size_t>(i)] = m[static_cast<size_t>(chosen[static_cast<size_t>(i)])];
        return out;
    }
}

// Full implementation of IndividuoFitness (was forward-declared earlier)
static std::vector<int> IndividuoFitness(const InitialPopulation& pop, const SAT3& sat) {
    std::vector<int> fitnessValues;
    if (std::holds_alternative<IntMatrix>(pop.individuals)) {
        const IntMatrix &m = std::get<IntMatrix>(pop.individuals);
        for (const auto &ind : m) {
            fitnessValues.push_back(clausestester(sat.getClauses(), ind));
        }
    } else {
        // if RealMatrix, convert each row to binary by sign (>=0 -> 1 else 0) as a heuristic
        const RealMatrix &m = std::get<RealMatrix>(pop.individuals);
        std::vector<int> tmp;
        tmp.reserve(pop.getDimension());
        for (const auto &row : m) {
            tmp.clear();
            for (double v : row) tmp.push_back(v >= 0.0 ? 1 : 0);
            fitnessValues.push_back(clausestester(sat.getClauses(), tmp));
        }
    }
    return fitnessValues;
}

// Helper to print a population (prints first N individuals)
static void printPopulationSample(const InitialPopulation& pop, int lines = 5) {
    std::cout << "Code=" << pop.getCode() << " POP=" << pop.getPopulation() << " DIM=" << pop.getDimension() << " TYPE=";
    switch (pop.getType()) {
        case GeneType::BIN: std::cout << "BIN"; break;
        case GeneType::INT: std::cout << "INT"; break;
        case GeneType::INT_PERM: std::cout << "INT-PERM"; break;
        case GeneType::REAL: std::cout << "REAL"; break;
    }
    std::cout << "\n";

    if (std::holds_alternative<IntMatrix>(pop.individuals)) {
        const IntMatrix &m = std::get<IntMatrix>(pop.individuals);
        int upto = std::min(lines, static_cast<int>(m.size()));
        for (int i = 0; i < upto; ++i) {
            for (size_t j = 0; j < m[i].size(); ++j) {
                std::cout << m[i][j];
                if (j + 1 < m[i].size()) std::cout << ' ';
            }
            std::cout << '\n';
        }
    } else {
        const RealMatrix &m = std::get<RealMatrix>(pop.individuals);
        int upto = std::min(lines, static_cast<int>(m.size()));
        for (int i = 0; i < upto; ++i) {
            for (size_t j = 0; j < m[i].size(); ++j) {
                std::cout << std::fixed << std::setprecision(4) << m[i][j];
                if (j + 1 < m[i].size()) std::cout << ' ';
            }
            std::cout << '\n';
        }
    }
}

bool read_dimacs_cnf(const std::string &path,
                     int &num_vars,
                     int &num_clauses,
                     std::vector<std::vector<int>> &clauses,
                     std::string &err)
{
    std::ifstream in(path);
    if (!in) { err = "Não foi possível abrir o ficheiro"; return false; }

    num_vars = num_clauses = 0;
    clauses.clear();
    std::string line;
    std::vector<int> cur; // acumula literais até encontrar 0
    bool header_seen = false;

    while (std::getline(in, line)) {
        // trim left
        size_t i = 0;
        while (i < line.size() && isspace((unsigned char)line[i])) ++i;
        if (i == line.size()) continue;
    if (line[i] == 'c' || line[i] == 'C') continue; // comentário

    // '%' indicates end of clauses in some DIMACS variants
    if (line[i] == '%') break;

        if (line[i] == 'p') {
            std::istringstream iss(line.substr(i));
            std::string p, fmt;
            if (!(iss >> p >> fmt >> num_vars >> num_clauses)) {
                err = "Cabeçalho 'p cnf' inválido";
                return false;
            }
            header_seen = true;
            continue;
        }

        // parse inteiros da linha (pode ser parte de cláusula; cláusula termina com 0)
        std::istringstream iss(line.substr(i));
        int lit;
        while (iss >> lit) {
            if (lit == 0) {
                if (!cur.empty()) {
                    clauses.push_back(cur);
                    cur.clear();
                } else {
                    // cláusula vazia (clause of size 0) -> cláusula impossível (opcional: tratar)
                    clauses.push_back({});
                }
            } else {
                cur.push_back(lit);
            }
        }
    }

    // if there's an unfinished clause (no trailing 0) we push it as a clause
    if (!cur.empty()) {
        clauses.push_back(cur);
        cur.clear();
    }

    // If header was present but counts differ, prefer actual read count and update num_clauses
    if (header_seen) {
        num_clauses = static_cast<int>(clauses.size());
    }

    return true;
}


static double decode_binary_to_real(const std::vector<int>& bits, double lo, double hi) {
    // interpret bits as unsigned integer (MSB at index 0)
    unsigned long long val = 0;
    for (size_t i = 0; i < bits.size(); ++i) {
        val = (val << 1) | (bits[i] ? 1ULL : 0ULL);
    }
    unsigned long long maxv;
    if (bits.size() == 0) maxv = 0;
    else if (bits.size() >= 64) maxv = ULLONG_MAX;
    else maxv = ((1ULL << bits.size()) - 1ULL);
    double ratio = (maxv == 0) ? 0.0 : (static_cast<double>(val) / static_cast<double>(maxv));
    double x = lo + ratio * (hi - lo);
    // round to 4 decimal places
    double scaled = std::round(x * 1e4) / 1e4;
    return scaled;
}

// The algebraic function from the slide (example): we'll implement f(x) = cos(20*x) - |x|/2 + fmod(x,3)/4
// (matches the earlier stray 'extras' expression). It's defined for real x; return double.
static double ex3_function(double x) {
    return std::cos(20.0 * x) - std::abs(x) / 2.0 + std::fmod(x, 3.0) / 4.0;
}

// Compute fitness for a binary-encoded population (IntMatrix) using EX3 function.
// mode="max" or "min" selects whether larger is better (maximization) or smaller is better (minimization).
static std::vector<double> EX3_fitness_from_binary(const IntMatrix &pop_bits, double lo, double hi, const std::string &mode = "max") {
    std::vector<double> out; out.reserve(pop_bits.size());
    for (const auto &ind : pop_bits) {
        double x = decode_binary_to_real(ind, lo, hi);
        double fx = ex3_function(x);
        if (mode == "min") fx = -fx; // convert minimization to maximization by negation
        out.push_back(fx);
    }
    return out;
}

// When called with EX3=LO,HI,BITS,MODE from CLI we will run a small demo and exit.
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Uso: " << argv[0] << " caminho_para_arquivo.cnf\n";
        return 1;
    }

    SAT3 sat;
    std::string err;
    if (!sat.loadFromFile(argv[1], err)) {
        std::cerr << "Erro a carregar CNF: " << err << "\n";
        return 1;
    }

    std::cout << "CNF carregado: vars=" << sat.getNumVars() << " clauses(header)=" << sat.getNumClauses() << " clauses(lidas)=" << sat.getClauses().size() << "\n";

    // ---------------------------
    // parse simple key=value args from argv[2..]
    // accepted keys (case-insensitive): COD, RUN, GEN, POP, DIM, BOUNDS, ELITISM
    // example: ./main file.cnf COD=BIN RUN=5 GEN=10 POP=50 DIM=100 BOUNDS=-5,10 ELITISM=true
    // ---------------------------
    auto up = [](std::string s){ for (char &c: s) c = static_cast<char>(toupper((unsigned char)c)); return s; };
    std::map<std::string,std::string> opts;
    for (int i = 2; i < argc; ++i) {
        std::string token = argv[i];
        auto eq = token.find('=');
        if (eq == std::string::npos) {
            // allow bare flag 'ELITISM'
            opts[up(token)] = "true";
        } else {
            std::string k = up(token.substr(0, eq));
            std::string v = token.substr(eq+1);
            opts[k] = v;
        }
    }

    std::string COD = "BIN";
    if (opts.count("COD")) COD = up(opts["COD"]);
    int RUN = opts.count("RUN") ? std::stoi(opts["RUN"]) : 1;
    int GEN = opts.count("GEN") ? std::stoi(opts["GEN"]) : 1;
    int POP = opts.count("POP") ? std::stoi(opts["POP"]) : 10;
    int DIM = opts.count("DIM") ? std::stoi(opts["DIM"]) : sat.getNumVars();
    bool ELITISM = false;
    if (opts.count("ELITISM")) {
        std::string v = up(opts["ELITISM"]);
        ELITISM = (v == "1" || v == "TRUE" || v == "YES" );
    }
    // parse bounds if provided like BOUNDS=-5,10
    IntBounds ib = {};
    RealBounds rb = {};
    if (opts.count("BOUNDS")) {
        std::string b = opts["BOUNDS"];
        auto comma = b.find(',');
        if (comma != std::string::npos) {
            std::string s1 = b.substr(0, comma);
            std::string s2 = b.substr(comma+1);
            try {
                if (COD == "REAL") {
                    rb.lo = std::stod(s1);
                    rb.hi = std::stod(s2);
                } else {
                    ib.lo = std::stoi(s1);
                    ib.hi = std::stoi(s2);
                }
            } catch(...) { /* ignore parse errors, keep defaults */ }
        }
    }

    std::cout << "Config: COD=" << COD << " RUN=" << RUN << " GEN=" << GEN << " POP=" << POP << " DIM=" << DIM << " ELITISM=" << (ELITISM?"ON":"OFF") << "\n";

    // ensure sensible defaults
    if (POP <= 0) POP = 10;
    if (DIM <= 0) DIM = sat.getNumVars();

    // Run simple random-search loop RUN x GEN with optional elitism.
    // This is a minimal example: each generation we re-generate a new population
    // and optionally preserve the best individual from the previous generation.
    std::vector<int> overallBestInd;
    int overallBestScore = -1;

    for (int run = 0; run < RUN; ++run) {
        // create population for this run
        InitialPopulation pop(COD, POP, DIM, ib, rb);
        std::vector<int> bestInd;
        int bestScore = -1;

        for (int gen = 0; gen < GEN; ++gen) {
            // evaluate current population
            auto fits = IndividuoFitness(pop, sat);
            // find best in this generation
            for (size_t i = 0; i < fits.size(); ++i) {
                if (fits[i] > bestScore) {
                    bestScore = fits[i];
                    // extract individual's integer representation (or convert real->bin)
                    if (std::holds_alternative<IntMatrix>(pop.individuals)) {
                        const IntMatrix &m = std::get<IntMatrix>(pop.individuals);
                        bestInd = m[i];
                    } else {
                        const RealMatrix &m = std::get<RealMatrix>(pop.individuals);
                        bestInd.clear(); bestInd.reserve(m[i].size());
                        for (double v : m[i]) bestInd.push_back(v >= 0.0 ? 1 : 0);
                    }
                }
            }

            // prepare next generation: fully regenerate, but if elitism preserve best in slot 0
            if (gen + 1 < GEN) {
                pop.generate();
                if (ELITISM && !bestInd.empty() && std::holds_alternative<IntMatrix>(pop.individuals)) {
                    IntMatrix &m = std::get<IntMatrix>(pop.individuals);
                    if (!m.empty()) m[0] = bestInd;
                }
            }
        }

        std::cout << "Run " << (run+1) << ": best=" << bestScore << " / " << sat.getNumClauses() << "\n";
        if (bestScore > overallBestScore) { overallBestScore = bestScore; overallBestInd = bestInd; }
    }

    std::cout << "Overall best across runs: " << overallBestScore << " / " << sat.getNumClauses() << "\n";
    if (!overallBestInd.empty()) {
        std::cout << "Best individuo (prefix 100 bits):\n";
        for (size_t i = 0; i < overallBestInd.size() && i < 100; ++i) std::cout << overallBestInd[i] << ' ';
        std::cout << '\n';
    }

    return 0;
}


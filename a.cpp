#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <thread>
#include <mutex>

struct Config {
    std::string cnfPath = "uf100-01.cnf";
    size_t population = 100;
    int generations = 1000;
    double mutationRate = 0.1;
    double eliteFraction = 0.1; // fraction of population kept as elite
    unsigned int seed = 0; // 0 means random
};
class inicial_population
{
private:
    int population;
    int dimension;

public:
    std::vector<std::vector<int>> individuals;

    inicial_population(int POP, int DIM)
        : population(POP), dimension(DIM) {}

    int getPopulation() const { return population; }
    int getDimension() const { return dimension; }

    void generate()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> d(0, 1); // Distribuição para gerar 0 ou 1

        std::vector<std::vector<int>> m(static_cast<size_t>(population), std::vector<int>(static_cast<size_t>(dimension)));

        for (int i = 0; i < population; i++)
        {
            for (int e = 0; e < dimension; e++)
            {
                m[static_cast<size_t>(i)][static_cast<size_t>(e)] = d(gen);
            }
        }

        individuals = std::move(m);
    }
};

bool read_dimacs_cnf(const std::string &path,
                     int &num_vars,
                     int &num_clauses,
                     std::vector<std::vector<int>> &clauses,
                     std::string &err)
{
    std::ifstream file(path);
    if (!file.is_open())
    {
        err = "Não foi possível abrir o arquivo.";
        return false;
    }

    auto ltrim = [](const std::string &s) {
        size_t i = 0;
        while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
        return s.substr(i);
    };

    std::string line;
    bool header_seen = false;
    while (std::getline(file, line))
    {
        line = ltrim(line);
        if (line.empty()) continue;
        if (line[0] == 'c') continue; // comment
        if (line[0] == '%') break;    // end

        if (line[0] == 'p')
        {
            std::istringstream iss(line);
            std::string p, cnf;
            if (!(iss >> p >> cnf >> num_vars >> num_clauses))
            {
                err = "Formato de cabeçalho inválido";
                return false;
            }
            clauses.clear();
            clauses.reserve(static_cast<size_t>(std::max(0, num_clauses)));
            header_seen = true;
            continue;
        }

        // cláusula
        std::istringstream iss(line);
        std::vector<int> clause;
        int lit;
        while (iss >> lit)
        {
            if (lit == 0) break;
            clause.push_back(lit);
        }
        if (!clause.empty()) clauses.push_back(std::move(clause));
    }
    file.close();

    if (header_seen && num_clauses >= 0 && static_cast<int>(clauses.size()) != num_clauses)
    {
        err = "Número de cláusulas lidas diferente do cabeçalho: header=" + std::to_string(num_clauses) +
              " parsed=" + std::to_string(clauses.size());
        return false;
    }

    return true;
}
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

    // Fitness function: returns number of satisfied clauses for an individual
    int getfitness(const std::vector<int>& individuo) const {
        return clausestester(clauses, individuo);
    }
};


static std::vector<int> IndividuoFitness(const std::vector<std::vector<int>> &pop, const SAT3 &sat)
{
    std::vector<int> fitnessValues;
    for (const auto &ind : pop)
    {
        fitnessValues.push_back(sat.getfitness(ind));
    }

    return fitnessValues;
}

static std::mt19937 &global_gen()
{
    // thread-local generator so parallel runs don't contend
    thread_local static std::mt19937 gen((std::random_device())());
    return gen;
}

static std::mutex g_print_mtx;

std::vector<int> rolleta(const std::vector<std::vector<int>> &pop, const SAT3 &sat) {
    size_t n = pop.size();
    if (n == 0) return {};
    std::vector<int> fitness = IndividuoFitness(pop, sat);
    // converter para double/long long se necessário
    long long total = std::accumulate(fitness.begin(), fitness.end(), 0LL);
    std::vector<int> selected;
    selected.reserve(n);

    if (total <= 0) {
        // fallback uniforme
        std::uniform_int_distribution<size_t> uni(0, n - 1);
        for (size_t i = 0; i < n; ++i) selected.push_back(static_cast<int>(uni(global_gen())));
        return selected;
    }

    // seleção por pares sem reposição dentro do par
    for (size_t pair = 0; pair < (n+1)/2; ++pair) {
        // recompute distribution for current weights
        std::discrete_distribution<int> dist(fitness.begin(), fitness.end());
        int p1 = dist(global_gen());
        // temporariamente zeroa p1
        int saved = fitness[static_cast<size_t>(p1)];
        fitness[static_cast<size_t>(p1)] = 0;
        // recompute distribution; if total becomes 0, fallback a uniforme (ex.: escolher aleatório)
        std::discrete_distribution<int> dist2(fitness.begin(), fitness.end());
        int p2 = dist2(global_gen());
        // restaura
        fitness[static_cast<size_t>(p1)] = saved;
        selected.push_back(p1);
        selected.push_back(p2);
        if (selected.size() >= n) break;
    }
    selected.resize(n); // garantir tamanho
    return selected;
}
std::vector<std::vector<int>> crossover(const  std::vector<std::vector<int>> &pop, const SAT3 &sat)
{   
    if (pop.empty()) return {};
    std::vector<int> selected = rolleta(pop, sat);
    std::vector<std::vector<int>> offspring;
    size_t geneCount = pop[0].size();
    if (geneCount < 2) return {};
    std::uniform_int_distribution<size_t> dist(1, geneCount - 1);

    for (size_t i = 0; i + 1 < selected.size(); i += 2)
    {
        size_t parent1Idx = static_cast<size_t>(selected[i]);
        size_t parent2Idx = static_cast<size_t>(selected[i + 1]);
        if (parent1Idx >= pop.size() || parent2Idx >= pop.size()) continue;

        size_t crossoverPoint = dist(global_gen());

        std::vector<int> child1(geneCount);
        std::vector<int> child2(geneCount);

        for (size_t j = 0; j < crossoverPoint; j++)
        {
            child1[j] = pop[parent1Idx][j];
            child2[j] = pop[parent2Idx][j];
        }
        for (size_t j = crossoverPoint; j < geneCount; j++)
        {
            child1[j] = pop[parent2Idx][j];
            child2[j] = pop[parent1Idx][j];
        }

        offspring.push_back(std::move(child1));
        offspring.push_back(std::move(child2));
    }

    return offspring;
}
std::vector<std::vector<int>> mutation(const std::vector<std::vector<int>> &pop, double mutationRate)
{
    std::vector<std::vector<int>> mutatedOffspring = pop;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (auto &individual : mutatedOffspring)
    {
        for (auto &gene : individual)
        {
            if (dist(global_gen()) < mutationRate)
            {
                gene = 1 - gene; // Inverte o gene (0 para 1 ou 1 para 0)
            }
        }
    }

    return mutatedOffspring;
}
std::vector<std::vector<int>> new_population(const std::vector<std::vector<int>> &pop, const SAT3 &sat, const Config &cfg)
{
    size_t N = pop.size();
    if (N == 0) return {};
    // Simple elitism: keep top fraction
    size_t elite = std::max<size_t>(1, static_cast<size_t>(std::ceil(cfg.eliteFraction * static_cast<double>(N))));
    std::vector<int> fitness = IndividuoFitness(pop, sat);
    std::vector<size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){ return fitness[a] > fitness[b]; });

    std::vector<std::vector<int>> newPop;
    newPop.reserve(N);
    for (size_t i = 0; i < elite && newPop.size() < N; ++i) newPop.push_back(pop[idx[i]]);

    // Generate offspring until we reach N
    while (newPop.size() < N) {
        std::vector<int> selected = rolleta(pop, sat);
        // need at least two selected indices
        if (selected.size() < 2) break;
        // use first two as parents
        std::vector<std::vector<int>> children;
        // create children from parents
        std::vector<std::vector<int>> parents = { pop[static_cast<size_t>(selected[0])], pop[static_cast<size_t>(selected[1])] };
        children = crossover(parents, sat);
        if (children.empty()) break;
    children = mutation(children, cfg.mutationRate);
        for (auto &c : children) {
            if (newPop.size() < N) newPop.push_back(std::move(c));
            else break;
        }
    }

    // If still short (rare), fill randomly
    std::uniform_int_distribution<size_t> uni(0, N - 1);
    while (newPop.size() < N) newPop.push_back(pop[uni(global_gen())]);

    return newPop;
}
static void print_usage(const char *prog) {
    std::cout << "Usage: " << prog << " [options] [cnf_path]\n";
    std::cout << "  -p N    population size (default 100)\n";
    std::cout << "  -g N    generations (default 1000)\n";
    std::cout << "  -m R    mutation rate (default 0.1)\n";
    std::cout << "  -e R    elite fraction (0..1, default 0.1)\n";
    std::cout << "  -s N    random seed (0 for random)\n";
    std::cout << "  -t N    number of parallel runs (default 10)\n";
}

int main(int argc, char **argv)
{
    Config cfg;
    int threads = 10; // default number of parallel runs
    // simple arg parsing
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-p" && i + 1 < argc) { cfg.population = static_cast<size_t>(std::stoul(argv[++i])); }
        else if (a == "-g" && i + 1 < argc) { cfg.generations = std::stoi(argv[++i]); }
        else if (a == "-m" && i + 1 < argc) { cfg.mutationRate = std::stod(argv[++i]); }
        else if (a == "-e" && i + 1 < argc) { cfg.eliteFraction = std::stod(argv[++i]); }
        else if (a == "-s" && i + 1 < argc) { cfg.seed = static_cast<unsigned int>(std::stoul(argv[++i])); }
        else if (a == "-t" && i + 1 < argc) { threads = std::max(1, std::stoi(argv[++i])); }
        else if (a == "-h" || a == "--help") { print_usage(argv[0]); return 0; }
        else if (i == argc - 1) { cfg.cnfPath = a; }
    }

    // spawn worker threads
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (int tid = 0; tid < threads; ++tid) {
        workers.emplace_back([tid, cfg]() {
            // seed per-thread generator
            if (cfg.seed != 0) {
                global_gen().seed(cfg.seed + static_cast<unsigned int>(tid));
            } else {
                std::random_device rd; global_gen().seed(rd() ^ (static_cast<unsigned int>(tid) + 0x9e3779b9));
            }

            SAT3 sat;
            std::string err;
            if (!sat.loadFromFile(cfg.cnfPath, err))
            {
                std::lock_guard<std::mutex> lk(g_print_mtx);
                std::cerr << "[run " << tid << "] Error loading CNF file ('" << cfg.cnfPath << "'): " << err << std::endl;
                return;
            }

            int dimension = sat.getNumVars();
            inicial_population pop(static_cast<int>(cfg.population), dimension);
            pop.generate();

            // prepare per-thread CSV logging
            std::ostringstream ss;
            ss << "evolution_stats_run" << tid + 1 << ".csv";
            std::string csvPath = ss.str();
            std::ofstream csv(csvPath);
            if (csv) csv << "generation,max_fitness,avg_fitness,best_individual\n";

            for (int generation = 0; generation < cfg.generations; generation++)
            {
                pop.individuals = new_population(pop.individuals, sat, cfg);
                std::vector<int> fitnessValues = IndividuoFitness(pop.individuals, sat);
                if (fitnessValues.empty()) break;
                int maxFitness = *std::max_element(fitnessValues.begin(), fitnessValues.end());
                double avgFitness = std::accumulate(fitnessValues.begin(), fitnessValues.end(), 0.0) / fitnessValues.size();

                auto it = std::max_element(fitnessValues.begin(), fitnessValues.end());
                size_t bestIdx = std::distance(fitnessValues.begin(), it);
                std::string bestStr;
                for (int bit : pop.individuals[bestIdx]) bestStr.push_back(bit ? '1' : '0');

                {
                    std::lock_guard<std::mutex> lk(g_print_mtx);
                    std::cout << "[run " << tid << "] Generation " << generation << ": Max Fitness = " << maxFitness << std::endl;
                }

                if (csv) csv << generation << "," << maxFitness << "," << avgFitness << "," << bestStr << "\n";

                if (maxFitness == sat.getNumClauses())
                {
                    std::lock_guard<std::mutex> lk(g_print_mtx);
                    std::cout << "[run " << tid << "] Solution found in generation " << generation << std::endl;
                    break;
                }
            }
            if (csv) csv.close();
        });
    }

    for (auto &t : workers) if (t.joinable()) t.join();

    return 0;
}
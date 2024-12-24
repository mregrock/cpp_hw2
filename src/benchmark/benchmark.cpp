#include "../types/defs.h"
#include "../selector/selector_bootstrap.h"
#include "../types/parse_helpers.h"
#include "../simulation/simulation_helpers.h"
#include <chrono>
#include <numeric>
#include <fstream>
#include <iomanip>

class Benchmark {
public:
    static constexpr size_t NUM_RUNS = 5;
    static constexpr size_t TICKS_PER_RUN = 200;

    struct Result {
        double total_time_ms;
        double avg_time_per_tick_ms;
        double ticks_per_second;
    };

    static void save_results(const std::vector<Result>& results, const std::string& filename) {
        double avg_total = 0, avg_per_tick = 0, avg_tps = 0;
        for (const auto& r : results) {
            avg_total += r.total_time_ms;
            avg_per_tick += r.avg_time_per_tick_ms;
            avg_tps += r.ticks_per_second;
        }
        avg_total /= results.size();
        avg_per_tick /= results.size();
        avg_tps /= results.size();

        std::ofstream out(filename);
        out << std::fixed << std::setprecision(3);
        
        for (size_t i = 0; i < results.size(); ++i) {
            out << "Run " << (i + 1) << ":\n";
            out << "  Total time: " << results[i].total_time_ms << " ms\n";
            out << "  Average time per tick: " << results[i].avg_time_per_tick_ms << " ms\n";
            out << "  Ticks per second: " << results[i].ticks_per_second << "\n\n";
        }

        out << "Average across " << results.size() << " runs:\n";
        out << "  Total time: " << avg_total << " ms\n";
        out << "  Average time per tick: " << avg_per_tick << " ms\n";
        out << "  Ticks per second: " << avg_tps << "\n";
    }
};

template<typename P, typename V, typename VF, typename Size=SizeType<dynamic_size, dynamic_size>>
void run_benchmark(const std::string& filename="field.txt") {    
    try {
        std::vector<Benchmark::Result> results;
        
        std::cerr << "Performing warm-up run...\n";
        {
            auto initial_field = read_field(filename);
            SimulationState<P, V, VF, Size::n, Size::m> state(initial_field);
                
            FluidSimulator<P, V, VF, Size::n, Size::m> simulator(state);
            simulator.run(true, Benchmark::TICKS_PER_RUN);
        }

        for (size_t run = 0; run < Benchmark::NUM_RUNS; ++run) {
            std::cerr << "Running benchmark " << (run + 1) << "/" << Benchmark::NUM_RUNS << "...\n";
            
            SimulationState<P, V, VF, Size::n, Size::m> state(initial_field);
                
            FluidSimulator<P, V, VF, Size::n, Size::m> simulator(state);
            
            auto start = std::chrono::high_resolution_clock::now();
            simulator.run(true, Benchmark::TICKS_PER_RUN);
            auto end = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            
            results.push_back({
                static_cast<double>(duration.count()),
                static_cast<double>(duration.count()) / Benchmark::TICKS_PER_RUN,
                (1000.0 * Benchmark::TICKS_PER_RUN) / duration.count()
            });
        }

        Benchmark::save_results(results, "benchmark_results.txt");
        std::cerr << "Benchmark completed. Results saved to benchmark_results.txt\n";

    } catch (const std::exception& e) {
        std::cerr << "Exception in run_simulation: " << e.what() << std::endl;
        throw;
    }
}

int main(int argc, char** argv) {
    run_benchmark<Fixed<32, 16>, Fixed<31, 17>, double, SizeType<36, 84>>();
    return 0;
}
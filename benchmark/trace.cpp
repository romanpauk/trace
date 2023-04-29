//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#include <benchmark/benchmark.h>
#include <trace/trace.h>
#include <string>
#include <iostream>

BENCHMARK_MAIN();

template< typename T > static void time_traits_get(benchmark::State& state)
{
    typename T::value_type start, end;
    for (auto _ : state) {
        start = T::get();
        end = T::get();
    }

    //std::cerr << T::diff(end, start) << std::endl;
    state.SetBytesProcessed(state.iterations()*2);
    benchmark::DoNotOptimize(start);
}

static void BM_StringCopy(benchmark::State& state)
{
    std::string x = "hellohellohellohellohellohellohellohellohello";
    for (auto _ : state) {
        std::string copy(x);
    }
}

static void BM_StringCopyTrace(benchmark::State& state)
{
    std::string x = "hellohellohellohellohellohellohellohellohello";
    for (auto _ : state) {
        TRACE("String");
        std::string copy(x);
    }
}

static void BM_MapCopy(benchmark::State& state)
{
    std::map< std::string, std::string > x = { {"A", "1"}, {"B", "2"} };
    for (auto _ : state) {
        auto copy = x;
    }
}

static void BM_MapCopyTrace(benchmark::State& state)
{
    std::map< std::string, std::string > x = { {"A", "1"}, {"B", "2"} };
    for (auto _ : state) {
        TRACE("Map");
        auto copy = x;
    }
}

#if defined(_WIN32)
BENCHMARK_TEMPLATE(time_traits_get, trace::qpc_time_traits);
#endif

#if defined(__linux__)
BENCHMARK_TEMPLATE(time_traits_get, trace::clock_gettime_time_traits< CLOCK_MONOTONIC >);
BENCHMARK_TEMPLATE(time_traits_get, trace::clock_gettime_time_traits< CLOCK_THREAD_CPUTIME_ID >);
#endif

BENCHMARK_TEMPLATE(time_traits_get, trace::default_time_traits);

BENCHMARK(BM_StringCopy);
BENCHMARK(BM_StringCopyTrace);
BENCHMARK(BM_MapCopy);
BENCHMARK(BM_MapCopyTrace);

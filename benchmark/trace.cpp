//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif

#include <benchmark/benchmark.h>
#include <trace/trace.h>
#include <string>
#include <iostream>

BENCHMARK_MAIN();

template< typename T > static void time_traits(benchmark::State& state)
{
    uint64_t start, end;
    for (auto _ : state) {
        start = T::begin();
        end = T::end();
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
BENCHMARK_TEMPLATE(time_traits, trace::time_traits< trace::qpc_counter >);
#endif

#if defined(__linux__)
BENCHMARK_TEMPLATE(time_traits, trace::time_traits< trace::clock_gettime_counter< CLOCK_MONOTONIC > >);
BENCHMARK_TEMPLATE(time_traits, trace::time_traits< trace::clock_gettime_counter< CLOCK_THREAD_CPUTIME_ID > >);
#endif

BENCHMARK_TEMPLATE(time_traits, trace::default_time_traits);

BENCHMARK(BM_StringCopy);
BENCHMARK(BM_StringCopyTrace);
BENCHMARK(BM_MapCopy);
BENCHMARK(BM_MapCopyTrace);

//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#include <benchmark/benchmark.h>
#include <trace/trace.h>
#include <string>

BENCHMARK_MAIN();

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

BENCHMARK(BM_StringCopy);
BENCHMARK(BM_StringCopyTrace);
BENCHMARK(BM_MapCopy);
BENCHMARK(BM_MapCopyTrace);
